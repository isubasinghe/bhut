//! Vulkan compute: integer histogram, separable blur, exposure, and packed RGB output.
//! Pixel bins are calculated in f64 on the CPU so bin edges match the CPU path.
use crate::{brightness, gaussian_kernel, Options, Result};
use wgpu::util::DeviceExt;

pub struct GpuDensity {
    device: wgpu::Device,
    queue: wgpu::Queue,
    bindings: wgpu::BindGroup,
    pipelines: Vec<wgpu::ComputePipeline>,
    config: wgpu::Buffer,
    counts: wgpu::Buffer,
    tiles: wgpu::Buffer,
    bins: wgpu::Buffer,
    output: wgpu::Buffer,
    readback: wgpu::Buffer,
    width: usize,
    height: usize,
    first: bool,
    packed_bins: Vec<u8>,
}

fn floats(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

impl GpuDensity {
    pub fn new(opts: &Options, particles: usize) -> Result<Self> {
        pollster::block_on(Self::create(opts, particles))
    }
    async fn create(opts: &Options, particles: usize) -> Result<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::VULKAN,
            ..Default::default()
        });
        let adapter = instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false, compatible_surface: None,
        }).await.ok_or("no Vulkan GPU found; enter nix develop, check driver/device access, or use --renderer cpu")?;
        let info = adapter.get_info();
        if info.device_type == wgpu::DeviceType::Cpu {
            return Err("Vulkan selected a software adapter; use --renderer cpu or enable GPU device access".into());
        }
        let size = opts
            .width
            .checked_mul(opts.height)
            .ok_or("GPU image too large")?;
        let limits = wgpu::Limits::default();
        if size > u32::MAX as usize / 2
            || size as u64 * 8 > limits.max_storage_buffer_binding_size as u64
            || particles as u64 * 4 > limits.max_storage_buffer_binding_size as u64
            || particles.div_ceil(256) > limits.max_compute_workgroups_per_dimension as usize
            || opts.width.div_ceil(16) > limits.max_compute_workgroups_per_dimension as usize
            || opts.height.div_ceil(16) > limits.max_compute_workgroups_per_dimension as usize
        {
            return Err("input/resolution exceeds GPU compute buffer limits; use --renderer cpu or a smaller resolution".into());
        }
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("density renderer"),
                    required_features: wgpu::Features::empty(),
                    required_limits: limits,
                },
                None,
            )
            .await?;
        eprintln!("GPU density renderer: {} ({:?})", info.name, info.backend);
        device.push_error_scope(wgpu::ErrorFilter::Validation);
        let buffer = |label, bytes, usage| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: bytes,
                usage,
                mapped_at_creation: false,
            })
        };
        let storage = wgpu::BufferUsages::STORAGE;
        let copy_dst = wgpu::BufferUsages::COPY_DST;
        let counts = buffer("histograms", size as u64 * 8, storage | copy_dst);
        let vertical = buffer("vertical blur", size as u64 * 8, storage);
        let smooth = buffer("blurred density", size as u64 * 8, storage);
        let tiles = buffer(
            "tile masks",
            (opts.width.div_ceil(32) * opts.height.div_ceil(32) * 4) as u64,
            storage | copy_dst,
        );
        let exposure = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("fixed exposure"),
            contents: &floats(&[1.0, 0.0, 0.0, 0.0]),
            usage: storage,
        });
        let output = buffer(
            "RGB",
            size as u64 * 3,
            storage | wgpu::BufferUsages::COPY_SRC,
        );
        let readback = buffer(
            "RGB readback",
            size as u64 * 3,
            wgpu::BufferUsages::MAP_READ | copy_dst,
        );
        let bins = buffer(
            "particle bins",
            (particles.max(1) * 4) as u64,
            storage | copy_dst,
        );
        let mut values = gaussian_kernel().to_vec();
        values.extend_from_slice(brightness::table());
        let lookup = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("kernel and brightness"),
            contents: &floats(&values),
            usage: storage,
        });
        let config = buffer("configuration", 32, wgpu::BufferUsages::UNIFORM | copy_dst);
        let layout_entries: Vec<_> = (0..9)
            .map(|binding| wgpu::BindGroupLayoutEntry {
                binding,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: if binding == 0 {
                        wgpu::BufferBindingType::Uniform
                    } else {
                        wgpu::BufferBindingType::Storage {
                            read_only: binding == 6 || binding == 7,
                        }
                    },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            })
            .collect();
        let bind_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: None,
            entries: &layout_entries,
        });
        let resources = [
            &config, &counts, &vertical, &smooth, &exposure, &output, &bins, &lookup, &tiles,
        ];
        let entries: Vec<_> = resources
            .iter()
            .enumerate()
            .map(|(binding, buffer)| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        let bindings = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &bind_layout,
            entries: &entries,
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[&bind_layout],
            push_constant_ranges: &[],
        });
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("density.wgsl"),
            source: wgpu::ShaderSource::Wgsl(include_str!("density.wgsl").into()),
        });
        let pipelines = [
            "histogram",
            "blur_vertical",
            "blur_horizontal",
            "measure_exposure",
            "prepare_curve",
            "colour",
        ]
        .iter()
        .map(|entry_point| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry_point),
                layout: Some(&layout),
                module: &shader,
                entry_point,
            })
        })
        .collect();
        if let Some(error) = device.pop_error_scope().await {
            return Err(format!("GPU renderer setup failed: {error}").into());
        }
        Ok(Self {
            device,
            queue,
            bindings,
            pipelines,
            config,
            counts,
            tiles,
            bins,
            output,
            readback,
            width: opts.width,
            height: opts.height,
            first: true,
            packed_bins: Vec::with_capacity(particles * 4),
        })
    }

    pub fn render(
        &mut self,
        points: &[[f64; 2]],
        groups: &[usize],
        centre: [f64; 2],
        half: f64,
        pixels: &mut [u8],
        fast: bool,
    ) -> Result<()> {
        let (w, h) = (self.width, self.height);
        let size = w * h;
        let scale = w.min(h) as f64 / (2.0 * half);
        self.packed_bins.clear();
        for (p, &group) in points.iter().zip(groups) {
            let x = (p[0] - centre[0]) * scale + w as f64 / 2.0;
            let y = (p[1] - centre[1]) * scale + h as f64 / 2.0;
            if x >= 0.0 && x <= w as f64 && y >= 0.0 && y <= h as f64 {
                let index = group * size + (y as usize).min(h - 1) * w + (x as usize).min(w - 1);
                self.packed_bins.extend((index as u32).to_le_bytes());
            }
        }
        let points = self.packed_bins.len() / 4;
        let params = [
            w as u32,
            h as u32,
            size as u32,
            points as u32,
            u32::from(fast),
            w.div_ceil(32) as u32,
            h.div_ceil(32) as u32,
            0,
        ];
        let bytes: Vec<_> = params.iter().flat_map(|v| v.to_le_bytes()).collect();
        self.queue.write_buffer(&self.config, 0, &bytes);
        if points != 0 {
            self.queue.write_buffer(&self.bins, 0, &self.packed_bins);
        }
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.clear_buffer(&self.counts, 0, None);
        encoder.clear_buffer(&self.tiles, 0, None);
        for stage in 0..6 {
            if (stage == 0 && points == 0) || ((stage == 3 || stage == 4) && !self.first) {
                continue;
            }
            // Separate passes establish storage-buffer visibility between stages.
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipelines[stage]);
            pass.set_bind_group(0, &self.bindings, &[]);
            if stage == 0 {
                pass.dispatch_workgroups(points.div_ceil(256) as u32, 1, 1);
            } else if stage == 4 {
                pass.dispatch_workgroups(1, 1, 1);
            } else if stage == 5 {
                pass.dispatch_workgroups((size / 4).div_ceil(256) as u32, 1, 1);
            } else {
                pass.dispatch_workgroups(
                    w.div_ceil(16) as u32,
                    h.div_ceil(16) as u32,
                    if stage == 1 || stage == 2 { 2 } else { 1 },
                );
            }
        }
        encoder.copy_buffer_to_buffer(&self.output, 0, &self.readback, 0, (size * 3) as u64);
        self.queue.submit(Some(encoder.finish()));
        let slice = self.readback.slice(..);
        let (tx, rx) = std::sync::mpsc::sync_channel(1);
        slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = tx.send(result);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()??;
        pixels.copy_from_slice(&slice.get_mapped_range());
        self.readback.unmap();
        self.first = false;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn shader_validates_without_gpu() {
        let source = include_str!("density.wgsl");
        let module = wgpu::naga::front::wgsl::parse_str(source)
            .unwrap_or_else(|error| panic!("{}", error.emit_to_string(source)));
        wgpu::naga::valid::Validator::new(
            wgpu::naga::valid::ValidationFlags::all(),
            wgpu::naga::valid::Capabilities::all(),
        )
        .validate(&module)
        .unwrap();
    }
}
