//! Render the simulator's native-endian binary output, streaming RGB to FFmpeg.
//! The original Matplotlib renderer is make_video.py.
mod frame_io;
use frame_io::{Backend, FrameReader};
use std::env;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread;
use std::time::Instant;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;
const COLOURS: [[f32; 3]; 2] = [[0.25, 0.70, 1.0], [1.0, 0.55, 0.20]];
const HELP: &str =
    "Make a video from simulation.bin (requires FFmpeg with drawtext and libx265/libx264).
Usage: make_video [simulation.bin] [options]
  -o, --output PATH         Output video (default: simulation.mp4; overwritten)
  --fps INTEGER             Frames per second (default: 10)
  --width INTEGER           Positive even output width (default: 3840)
  --height INTEGER          Positive even output height (default: 2160)
  --render-mode MODE        density or particles (default: density)
  --particle-size NUMBER    Marker diameter in output pixels (default: 2)
  --fraction NUMBER         Fraction of positions in view (default: 0.95)
  --elev NUMBER             Camera elevation in degrees (default: 25)
  --azim NUMBER             Camera azimuth in degrees (default: 45)
  --codec NAME              h265 or h264 (default: h265)
  --preset NAME             x265/x264 speed preset (default: fast)
  --crf NUMBER              Quality, 0–51; lower is better (default: 18)
  --io-backend NAME         auto, uring (Linux), or threaded (default: auto)
  -h, --help                Show this help";

struct Options {
    input: String,
    output: String,
    fps: u32,
    width: usize,
    height: usize,
    particles: bool,
    particle_size: f64,
    fraction: f64,
    elev: f64,
    azim: f64,
    h265: bool,
    preset: String,
    crf: u8,
    backend: Backend,
}

impl Options {
    fn parse(args: impl IntoIterator<Item = String>) -> Result<Option<Self>> {
        let mut options = Self {
            input: "simulation.bin".into(),
            output: "simulation.mp4".into(),
            fps: 10,
            width: 3840,
            height: 2160,
            particles: false,
            particle_size: 2.0,
            fraction: 0.95,
            elev: 25.0,
            azim: 45.0,
            h265: true,
            preset: "fast".into(),
            crf: 18,
            backend: Backend::Auto,
        };
        let mut args = args.into_iter();
        let mut positional = false;
        while let Some(arg) = args.next() {
            if arg == "--help" || arg == "-h" {
                return Ok(None);
            }
            if !arg.starts_with('-') {
                if positional {
                    return Err(format!("unexpected argument: {arg}").into());
                }
                options.input = arg;
                positional = true;
                continue;
            }
            let (key, inline) = arg
                .split_once('=')
                .map_or((arg.as_str(), None), |(k, v)| (k, Some(v)));
            if !matches!(
                key,
                "-o" | "--output"
                    | "--fps"
                    | "--width"
                    | "--height"
                    | "--render-mode"
                    | "--particle-size"
                    | "--fraction"
                    | "--elev"
                    | "--azim"
                    | "--codec"
                    | "--preset"
                    | "--crf"
                    | "--io-backend"
            ) {
                return Err(format!("unknown option: {key}").into());
            }
            let value = inline
                .map(str::to_owned)
                .or_else(|| args.next())
                .ok_or_else(|| format!("missing value for {key}"))?;
            match key {
                "-o" | "--output" => options.output = value,
                "--fps" => options.fps = value.parse()?,
                "--width" => options.width = value.parse()?,
                "--height" => options.height = value.parse()?,
                "--render-mode" => {
                    options.particles = match value.as_str() {
                        "particles" => true,
                        "density" => false,
                        _ => return Err("render mode must be density or particles".into()),
                    }
                }
                "--particle-size" => options.particle_size = value.parse()?,
                "--fraction" => options.fraction = value.parse()?,
                "--elev" => options.elev = value.parse()?,
                "--azim" => options.azim = value.parse()?,
                "--codec" => {
                    options.h265 = match value.as_str() {
                        "h265" => true,
                        "h264" => false,
                        _ => return Err("codec must be h265 or h264".into()),
                    }
                }
                "--preset" => options.preset = value,
                "--crf" => options.crf = value.parse()?,
                "--io-backend" => {
                    options.backend = match value.as_str() {
                        "auto" => Backend::Auto,
                        "uring" => Backend::Uring,
                        "threaded" => Backend::Threaded,
                        _ => return Err("I/O backend must be auto, uring, or threaded".into()),
                    }
                }
                _ => unreachable!(),
            }
        }
        if options.fps == 0 || !(options.fraction > 0.0 && options.fraction <= 1.0) {
            return Err("fps must be positive and fraction must be in (0, 1]".into());
        }
        if options.width == 0
            || options.height == 0
            || options.width % 2 != 0
            || options.height % 2 != 0
        {
            return Err("width and height must be positive even numbers for yuv420p output".into());
        }
        options
            .width
            .checked_mul(options.height)
            .and_then(|n| n.checked_mul(3))
            .ok_or("output dimensions are too large")?;
        if !options.particle_size.is_finite() || options.particle_size <= 0.0 {
            return Err("particle size must be finite and positive".into());
        }
        if !options.elev.is_finite() || !options.azim.is_finite() {
            return Err("camera angles must be finite".into());
        }
        if options.crf > 51 {
            return Err("CRF must be in 0–51".into());
        }
        if !matches!(
            options.preset.as_str(),
            "ultrafast"
                | "superfast"
                | "veryfast"
                | "faster"
                | "fast"
                | "medium"
                | "slow"
                | "slower"
                | "veryslow"
                | "placebo"
        ) {
            return Err("invalid x264/x265 speed preset".into());
        }
        Ok(Some(options))
    }
}

struct Simulation {
    file: File,
    n: usize,
    frames: usize,
    offset: u64,
    bytes: Vec<u8>,
    positions: Vec<[f64; 3]>,
}

impl Simulation {
    fn open(path: &str) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut header = [0; 8];
        file.read_exact(&mut header)?;
        let n = i32::from_ne_bytes(header[..4].try_into()?);
        let declared_frames = i32::from_ne_bytes(header[4..].try_into()?);
        if n < 2 || declared_frames < 1 {
            return Err("input must declare at least two particles and one frame".into());
        }
        let n = n as usize;
        let offset = 8 + 8 * n as u64;
        let frame_bytes = n.checked_mul(24).ok_or("frame is too large")?;
        // Snapshot the available complete frames, including during an active simulation.
        let available = file.metadata()?.len().saturating_sub(offset) / frame_bytes as u64;
        let frames = available.min(declared_frames as u64) as usize;
        if frames == 0 {
            return Err("input contains no complete frames".into());
        }
        Ok(Self {
            file,
            n,
            frames,
            offset,
            bytes: vec![0; frame_bytes],
            positions: vec![[0.0; 3]; n],
        })
    }

    fn rewind(&mut self) -> Result<()> {
        self.file.seek(SeekFrom::Start(self.offset))?;
        Ok(())
    }

    fn read_frame(&mut self) -> Result<&[[f64; 3]]> {
        self.file.read_exact(&mut self.bytes)?;
        for (position, bytes) in self.positions.iter_mut().zip(self.bytes.chunks_exact(24)) {
            for (coordinate, bytes) in position.iter_mut().zip(bytes.chunks_exact(8)) {
                *coordinate = f64::from_ne_bytes(bytes.try_into()?);
                if !coordinate.is_finite() {
                    return Err("input positions must be finite".into());
                }
            }
        }
        Ok(&self.positions)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

struct Projection {
    right: [f64; 3],
    up: [f64; 3],
}
impl Projection {
    fn new(elev: f64, azim: f64) -> Self {
        let (se, ce) = elev.to_radians().sin_cos();
        let (sa, ca) = azim.to_radians().sin_cos();
        Self {
            right: [-sa, ca, 0.0],
            up: [-se * ca, -se * sa, ce],
        }
    }
    fn point(&self, p: [f64; 3]) -> [f64; 2] {
        [dot(p, self.right), dot(p, self.up)]
    }
}

// Jacobi diagonalization of the symmetric 3x3 covariance matrix.
fn principal_axis(mut matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let mut axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (p, q) = [(0, 1), (0, 2), (1, 2)]
            .into_iter()
            .max_by(|&(a, b), &(c, d)| matrix[a][b].abs().total_cmp(&matrix[c][d].abs()))
            .unwrap();
        let scale = (0..3).map(|i| matrix[i][i].abs()).fold(0.0, f64::max);
        if matrix[p][q].abs() <= f64::EPSILON * scale {
            break;
        }
        let angle = 0.5 * (2.0 * matrix[p][q]).atan2(matrix[q][q] - matrix[p][p]);
        let (s, c) = angle.sin_cos();
        let (pp, qq, pq) = (matrix[p][p], matrix[q][q], matrix[p][q]);
        matrix[p][p] = c * c * pp - 2.0 * s * c * pq + s * s * qq;
        matrix[q][q] = s * s * pp + 2.0 * s * c * pq + c * c * qq;
        matrix[p][q] = 0.0;
        matrix[q][p] = 0.0;
        for k in 0..3 {
            if k != p && k != q {
                let (kp, kq) = (matrix[k][p], matrix[k][q]);
                matrix[k][p] = c * kp - s * kq;
                matrix[p][k] = matrix[k][p];
                matrix[k][q] = s * kp + c * kq;
                matrix[q][k] = matrix[k][q];
            }
            let (kp, kq) = (axes[k][p], axes[k][q]);
            axes[k][p] = c * kp - s * kq;
            axes[k][q] = s * kp + c * kq;
        }
    }
    let i = (0..3)
        .max_by(|&a, &b| matrix[a][a].total_cmp(&matrix[b][b]))
        .unwrap();
    let mut axis = [axes[0][i], axes[1][i], axes[2][i]];
    // Eigenvector signs are arbitrary; use a stable convention for group colours.
    let largest = (0..3)
        .max_by(|&a, &b| axis[a].abs().total_cmp(&axis[b].abs()))
        .unwrap();
    if axis[largest] < 0.0 {
        axis.iter_mut().for_each(|v| *v = -*v);
    }
    axis
}

fn quantile(values: &mut [f64], fraction: f64) -> f64 {
    let index = (values.len() - 1) as f64 * fraction;
    let lower = index.floor() as usize;
    let (_, value, rest) = values.select_nth_unstable_by(lower, f64::total_cmp);
    let low = *value;
    if index.fract() == 0.0 {
        return low;
    }
    let high = rest.iter().copied().min_by(f64::total_cmp).unwrap();
    low + (high - low) * index.fract()
}

fn initial_groups(pos: &[[f64; 3]]) -> Vec<usize> {
    let mut centre = [0.0; 3];
    for p in pos {
        for k in 0..3 {
            centre[k] += p[k] / pos.len() as f64;
        }
    }
    let mut covariance = [[0.0; 3]; 3];
    for p in pos {
        for i in 0..3 {
            for j in 0..3 {
                covariance[i][j] += (p[i] - centre[i]) * (p[j] - centre[j]) / pos.len() as f64;
            }
        }
    }
    let axis = principal_axis(covariance);
    let along: Vec<_> = pos
        .iter()
        .map(|p| dot(std::array::from_fn(|i| p[i] - centre[i]), axis))
        .collect();
    let median = quantile(&mut along.clone(), 0.5);
    let mut groups: Vec<_> = along.iter().map(|&v| usize::from(v > median)).collect();
    for _ in 0..50 {
        let mut counts = [0; 2];
        let mut centres = [[0.0; 3]; 2];
        for (p, &group) in pos.iter().zip(&groups) {
            counts[group] += 1;
            for k in 0..3 {
                centres[group][k] += p[k];
            }
        }
        if counts.contains(&0) {
            break;
        }
        for g in 0..2 {
            for k in 0..3 {
                centres[g][k] /= counts[g] as f64;
            }
        }
        let mut changed = false;
        for (p, group) in pos.iter().zip(&mut groups) {
            let distances: [f64; 2] =
                std::array::from_fn(|g| (0..3).map(|k| (p[k] - centres[g][k]).powi(2)).sum());
            let updated = usize::from(distances[1] < distances[0]);
            changed |= updated != *group;
            *group = updated;
        }
        if !changed {
            break;
        }
    }
    groups
}

fn fixed_view(
    sim: &mut Simulation,
    projection: &Projection,
    fraction: f64,
) -> Result<([f64; 2], f64)> {
    // Three sequential scans keep exact NumPy-style quantiles while storing only
    // one scalar per particle/frame, rather than all 3D and projected positions.
    let count = sim.n.checked_mul(sim.frames).ok_or("too many positions")?;
    let mut values = Vec::new();
    values.try_reserve_exact(count)?;
    let mut centre = [0.0; 2];
    let mut half = 0.0;
    for pass in 0..3 {
        eprintln!("Measuring fixed camera ({}/3)...", pass + 1);
        values.clear();
        sim.rewind()?;
        for _ in 0..sim.frames {
            for &p in sim.read_frame()? {
                let screen = projection.point(p);
                let value = if pass < 2 {
                    screen[pass]
                } else {
                    (screen[0] - centre[0])
                        .abs()
                        .max((screen[1] - centre[1]).abs())
                };
                if !value.is_finite() {
                    return Err("projected positions exceed the finite coordinate range".into());
                }
                values.push(value);
            }
        }
        if pass < 2 {
            centre[pass] = quantile(&mut values, 0.5);
        } else {
            half = (quantile(&mut values, fraction) * 1.05).max(1e-12);
        }
    }
    if !half.is_finite() {
        return Err("camera bounds exceed the finite coordinate range".into());
    }
    Ok((centre, half))
}

fn blur(layer: &mut [f32], scratch: &mut [f32], width: usize, height: usize) {
    let mut kernel: [f32; 9] =
        std::array::from_fn(|i| (-0.5 * ((i as f32 - 4.0) / 1.4).powi(2)).exp());
    let sum: f32 = kernel.iter().sum();
    kernel.iter_mut().for_each(|v| *v /= sum);
    // Finish each row while it is in cache instead of scanning the whole 4K
    // image nine times. Keep the accumulation order and zero padding unchanged.
    for (y, out) in scratch.chunks_exact_mut(width).enumerate() {
        out.fill(0.0);
        for (k, &weight) in kernel.iter().enumerate() {
            let source = y as isize + k as isize - 4;
            if source < 0 || source >= height as isize {
                continue;
            }
            let row = &layer[source as usize * width..(source as usize + 1) * width];
            for (out, &value) in out.iter_mut().zip(row) {
                *out += weight * value;
            }
        }
    }
    for (y, out) in layer.chunks_exact_mut(width).enumerate() {
        out.fill(0.0);
        for (k, &weight) in kernel.iter().enumerate() {
            let shift = k as isize - 4;
            let start = (-shift).max(0) as usize;
            let end = (width as isize - shift).min(width as isize).max(0) as usize;
            if start >= end {
                continue;
            }
            let source = (y * width + start) as isize + shift;
            let input = &scratch[source as usize..source as usize + end - start];
            for (out, &value) in out[start..end].iter_mut().zip(input) {
                *out += weight * value;
            }
        }
    }
}

struct Density {
    layers: [Vec<f32>; 2],
    scratch: Vec<f32>,
}
impl Density {
    fn new(size: usize) -> Self {
        Self {
            layers: [vec![0.0; size], vec![0.0; size]],
            scratch: vec![0.0; size],
        }
    }
    fn render(
        &mut self,
        points: &[[f64; 2]],
        groups: &[usize],
        centre: [f64; 2],
        half: f64,
        opts: &Options,
    ) {
        let (w, h) = (opts.width, opts.height);
        let scale = w.min(h) as f64 / (2.0 * half);
        self.layers.iter_mut().for_each(|layer| layer.fill(0.0));
        for (p, &group) in points.iter().zip(groups) {
            let x = (p[0] - centre[0]) * scale + w as f64 / 2.0;
            let y = (p[1] - centre[1]) * scale + h as f64 / 2.0;
            if x >= 0.0 && x <= w as f64 && y >= 0.0 && y <= h as f64 {
                // As in histogram2d, include the upper bound in the last bin.
                let (x, y) = ((x as usize).min(w - 1), (y as usize).min(h - 1));
                self.layers[group][y * w + x] += 1.0;
            }
        }
        for layer in &mut self.layers {
            blur(layer, &mut self.scratch, w, h);
        }
    }
    fn exposure(&self) -> f32 {
        self.layers[0]
            .iter()
            .zip(&self.layers[1])
            .map(|(a, b)| a + b)
            .fold(1.0, f32::max)
    }
    fn rgb(&self, pixels: &mut [u8], width: usize, height: usize, exposure: f32) {
        let denominator = (10.0 * exposure).ln_1p();
        for y in 0..height {
            let row = &mut pixels[(height - 1 - y) * width * 3..(height - y) * width * 3];
            for (x, pixel) in row.chunks_exact_mut(3).enumerate() {
                let a = self.layers[0][y * width + x];
                let b = self.layers[1][y * width + x];
                let total = a + b;
                if total == 0.0 {
                    pixel.fill(0);
                    continue;
                }
                let brightness = ((10.0 * total).ln_1p() / denominator).min(1.0).powf(0.6);
                for c in 0..3 {
                    pixel[c] = (255.0 * (a * COLOURS[0][c] + b * COLOURS[1][c]) / total.max(1e-12)
                        * brightness) as u8;
                }
            }
        }
    }
}

fn particles(
    pixels: &mut [u8],
    points: &[[f64; 2]],
    groups: &[usize],
    centre: [f64; 2],
    half: f64,
    opts: &Options,
) {
    pixels.fill(0);
    let (w, h) = (opts.width as f64, opts.height as f64);
    let size = (0.96 * w).min(0.88 * h);
    let (left, top) = ((w - size) / 2.0, 0.08 * h + (0.88 * h - size) / 2.0);
    let radius = opts.particle_size / 2.0;
    for (p, &group) in points.iter().zip(groups) {
        let x = left + ((p[0] - centre[0]) / half + 1.0) * size / 2.0;
        let y = top + (1.0 - (p[1] - centre[1]) / half) * size / 2.0;
        if x + radius < left
            || x - radius > left + size
            || y + radius < top
            || y - radius > top + size
        {
            continue;
        }
        let x0 = (x - radius - 1.0).max(left).max(0.0) as usize;
        let x1 = (x + radius + 1.0).min(left + size).min(w).ceil() as usize;
        let y0 = (y - radius - 1.0).max(top).max(0.0) as usize;
        let y1 = (y + radius + 1.0).min(top + size).min(h).ceil() as usize;
        for py in y0..y1 {
            for px in x0..x1 {
                let distance = (px as f64 + 0.5 - x).hypot(py as f64 + 0.5 - y);
                let alpha = (radius + 0.5 - distance).clamp(0.0, 1.0) as f32;
                let pixel = &mut pixels[(py * opts.width + px) * 3..(py * opts.width + px) * 3 + 3];
                for c in 0..3 {
                    pixel[c] =
                        (pixel[c] as f32 * (1.0 - alpha) + 255.0 * COLOURS[group][c] * alpha) as u8;
                }
            }
        }
    }
}

fn read_ahead(
    mut input: FrameReader,
    projection: &Projection,
    frames: usize,
    ready: SyncSender<Vec<[f64; 2]>>,
    recycled: Receiver<Vec<[f64; 2]>>,
) -> Result<()> {
    for _ in 0..frames {
        let Ok(mut points) = recycled.recv() else {
            break;
        };
        points.clear();
        for bytes in input.read_frame()?.chunks_exact(24) {
            let mut position = [0.0; 3];
            for (value, bytes) in position.iter_mut().zip(bytes.chunks_exact(8)) {
                *value = f64::from_ne_bytes(bytes.try_into()?);
            }
            let point = projection.point(position);
            if !position.iter().chain(point.iter()).all(|v| v.is_finite()) {
                return Err("input and projected positions must be finite".into());
            }
            points.push(point);
        }
        if ready.send(points).is_err() {
            break;
        }
    }
    Ok(())
}

fn write_frames(
    mut output: impl Write,
    ready: Receiver<Vec<u8>>,
    recycled: SyncSender<Vec<u8>>,
) -> Result<()> {
    for pixels in ready {
        output.write_all(&pixels)?;
        // The renderer may already be done. Still drain all queued frames before EOF.
        let _ = recycled.send(pixels);
    }
    output.flush()?;
    Ok(())
}

fn render_pipeline(
    input: FrameReader,
    output: impl Write + Send,
    projection: &Projection,
    groups: &[usize],
    frames: usize,
    centre: [f64; 2],
    half: f64,
    opts: &Options,
) -> Result<()> {
    // Two reusable buffers at each boundary bound memory and apply backpressure.
    // Reading/projecting and writing to FFmpeg happen independently of rendering.
    let (input_tx, input_rx) = sync_channel(1);
    let (input_free_tx, input_free_rx) = sync_channel(2);
    let (output_tx, output_rx) = sync_channel(1);
    let (output_free_tx, output_free_rx) = sync_channel(2);
    for _ in 0..2 {
        input_free_tx.send(Vec::with_capacity(groups.len()))?;
        output_free_tx.send(vec![0; opts.width * opts.height * 3])?;
    }
    thread::scope(|scope| -> Result<()> {
        let reader = scope.spawn(|| read_ahead(input, projection, frames, input_tx, input_free_rx));
        let writer = scope.spawn(|| write_frames(output, output_rx, output_free_tx));
        let rendering = (|| -> Result<()> {
            let mut density = if opts.particles {
                None
            } else {
                Some(Density::new(opts.width * opts.height))
            };
            let mut exposure = 1.0;
            for frame in 0..frames {
                let points = input_rx.recv().map_err(|_| "input worker stopped")?;
                let mut pixels = output_free_rx
                    .recv()
                    .map_err(|_| "encoder writer stopped")?;
                if let Some(density) = &mut density {
                    density.render(&points, groups, centre, half, opts);
                    if frame == 0 {
                        exposure = density.exposure();
                    }
                    density.rgb(&mut pixels, opts.width, opts.height, exposure);
                } else {
                    particles(&mut pixels, &points, groups, centre, half, opts);
                }
                let _ = input_free_tx.send(points);
                output_tx
                    .send(pixels)
                    .map_err(|_| "encoder writer stopped")?;
                eprint!("\rframe {}/{} queued", frame + 1, frames);
                io::stderr().flush()?;
            }
            eprintln!("\nFinishing encoder...");
            Ok(())
        })();
        // Disconnect both directions before joining, including on errors. A producer
        // blocked on a full queue or an empty buffer pool must be allowed to exit.
        drop(input_rx);
        drop(input_free_tx);
        drop(output_tx);
        drop(output_free_rx);
        let read_result = reader.join().map_err(|_| "input worker panicked");
        let write_result = writer.join().map_err(|_| "encoder writer panicked");
        read_result??;
        write_result??;
        rendering
    })
}

fn encode(
    sim: &mut Simulation,
    projection: &Projection,
    groups: &[usize],
    centre: [f64; 2],
    half: f64,
    opts: &Options,
) -> Result<()> {
    let caption = if opts.particles {
        "Individual particles · fixed camera"
    } else {
        "Projected particle density · logarithmic brightness · fixed camera"
    };
    let filter = format!("setsar=1,drawtext=text='Galaxy collision  |  step %{{n}}':fontcolor=white:fontsize=27:x=w*0.04:y=h*0.035,drawtext=text='Initial group 1':fontcolor=0x40b3ff:fontsize=21:x=w*0.04:y=h*0.075,drawtext=text='Initial group 2':fontcolor=0xff8c33:fontsize=21:x=w*0.24:y=h*0.075,drawtext=text='{caption}':fontcolor=0xa6a6a6:fontsize=19:x=w*0.04:y=h*0.975-th");
    let input = FrameReader::new(sim.file.try_clone()?, sim.offset, sim.n * 24, opts.backend)?;
    eprintln!(
        "Input: {}; codec: {}; preset: {}",
        input.name(),
        if opts.h265 { "H.265" } else { "H.264" },
        opts.preset
    );
    let mut command = Command::new("ffmpeg");
    command
        .args([
            "-y",
            "-loglevel",
            "error",
            "-f",
            "rawvideo",
            "-pixel_format",
            "rgb24",
            "-video_size",
        ])
        .arg(format!("{}x{}", opts.width, opts.height))
        .args([
            "-framerate",
            &opts.fps.to_string(),
            "-i",
            "pipe:0",
            "-an",
            "-c:v",
            if opts.h265 { "libx265" } else { "libx264" },
            "-preset",
            &opts.preset,
            "-crf",
            &opts.crf.to_string(),
            "-pix_fmt",
            "yuv420p",
            "-vf",
            &filter,
            "-movflags",
            "+faststart",
        ]);
    if opts.h265 {
        command.args(["-tag:v", "hvc1", "-x265-params", "log-level=error"]);
    }
    let mut child = command
        .arg(&opts.output)
        .stdin(Stdio::piped())
        .spawn()
        .map_err(|e| format!("could not start FFmpeg (enter nix develop): {e}"))?;
    let pipe = child.stdin.take().expect("piped FFmpeg stdin");
    let started = Instant::now();
    let result = render_pipeline(
        input, pipe, projection, groups, sim.frames, centre, half, opts,
    );
    let status = child.wait()?;
    result?;
    if !status.success() {
        return Err(format!("FFmpeg exited with {status}").into());
    }
    eprintln!(
        "Rendered and encoded {} frames in {:.2}s ({:.2} frames/s)",
        sim.frames,
        started.elapsed().as_secs_f64(),
        sim.frames as f64 / started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn run() -> Result<()> {
    let Some(opts) = Options::parse(env::args().skip(1))? else {
        println!("{HELP}");
        return Ok(());
    };
    if Path::new(&opts.output).exists()
        && std::fs::canonicalize(&opts.input)? == std::fs::canonicalize(&opts.output)?
    {
        return Err("input and output must be different files".into());
    }
    let mut sim = Simulation::open(&opts.input)?;
    sim.rewind()?;
    let groups = initial_groups(sim.read_frame()?);
    let projection = Projection::new(opts.elev, opts.azim);
    let (centre, half) = fixed_view(&mut sim, &projection, opts.fraction)?;
    eprintln!(
        "{} planets, {} frames, centre {:?}, half width {:.4}",
        sim.n, sim.frames, centre, half
    );
    let orange = groups.iter().filter(|&&g| g == 1).count();
    eprintln!(
        "Initial spatial groups: {} blue, {} orange",
        sim.n - orange,
        orange
    );
    encode(&mut sim, &projection, &groups, centre, half, &opts)?;
    println!("wrote {}", opts.output);
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("make_video: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolated_quantiles() {
        for (fraction, expected) in [(0.0, 1.0), (0.5, 3.5), (0.95, 9.25), (1.0, 10.0)] {
            assert!((quantile(&mut [10.0, 2.0, 1.0, 5.0], fraction) - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn projection_orientation() {
        let p = Projection::new(0.0, 0.0).point([1.0, 2.0, 3.0]);
        assert_eq!(p, [2.0, 3.0]);
    }

    #[test]
    fn clusters_follow_rotated_principal_axis() {
        let pos = [
            [-10.0, -20.0, -30.0],
            [-11.0, -19.0, -31.0],
            [10.0, 20.0, 30.0],
            [11.0, 19.0, 31.0],
        ];
        let groups = initial_groups(&pos);
        assert_eq!(groups[0], groups[1]);
        assert_eq!(groups[2], groups[3]);
        assert_ne!(groups[0], groups[2]);
        assert_eq!(initial_groups(&[[1.0; 3]; 2]), [0, 0]);
    }

    #[test]
    fn blur_preserves_interior_mass_and_handles_tiny_images() {
        let mut layer = vec![0.0; 81];
        layer[40] = 1.0;
        blur(&mut layer, &mut [0.0; 81], 9, 9);
        assert!((layer.iter().sum::<f32>() - 1.0).abs() < 1e-6);
        assert!((layer[40] - 0.0813689).abs() < 1e-6);
        assert_eq!(layer[39], layer[41]);
        let mut tiny = [1.0; 4];
        blur(&mut tiny, &mut [0.0; 4], 2, 2);
        assert!(tiny.iter().all(|v| *v > 0.0 && *v < 1.0));
    }

    #[test]
    fn invalid_options_are_rejected() {
        for args in [
            ["--width", "3"],
            ["--height", "0"],
            ["--fps", "0"],
            ["--fraction", "NaN"],
            ["--particle-size", "inf"],
            ["--elev", "NaN"],
            ["--render-mode", "other"],
            ["--codec", "other"],
            ["--io-backend", "other"],
            ["--crf", "52"],
            ["--preset", "other"],
        ] {
            assert!(Options::parse(args.map(str::to_owned)).is_err(), "{args:?}");
        }
    }

    struct TestFile(std::path::PathBuf);
    impl TestFile {
        fn new(data: &[u8]) -> io::Result<Self> {
            use std::sync::atomic::{AtomicUsize, Ordering};
            static NEXT: AtomicUsize = AtomicUsize::new(0);
            let path = env::temp_dir().join(format!(
                "bhut-pipeline-{}-{}.bin",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)?;
            file.write_all(data)?;
            Ok(Self(path))
        }
        fn reader(&self, size: usize, backend: Backend) -> io::Result<FrameReader> {
            FrameReader::new(File::open(&self.0)?, 0, size, backend)
        }
    }
    impl Drop for TestFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    #[test]
    fn frame_readers_preserve_order_and_report_short_reads() -> Result<()> {
        let fixture = TestFile::new(&[1, 2, 3, 4, 5])?;
        for backend in [Backend::Threaded, Backend::Auto, Backend::Uring] {
            let mut reader = match fixture.reader(2, backend) {
                Ok(reader) => reader,
                Err(error)
                    if backend == Backend::Uring
                        && env::var_os("BHUT_TEST_URING").is_none()
                        && matches!(
                            error.kind(),
                            io::ErrorKind::PermissionDenied | io::ErrorKind::Unsupported
                        ) =>
                {
                    eprintln!("Skipping unavailable forced io_uring: {error}");
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            assert_eq!(reader.read_frame()?, [1, 2]);
            assert_eq!(reader.read_frame()?, [3, 4]);
            assert_eq!(
                reader.read_frame().unwrap_err().kind(),
                io::ErrorKind::UnexpectedEof
            );
        }
        Ok(())
    }

    #[test]
    fn pipeline_preserves_pixels_order_and_drains_last_frame() -> Result<()> {
        let opts = Options::parse(["--width", "64", "--height", "64"].map(str::to_owned))?.unwrap();
        let projection = Projection::new(25.0, 45.0);
        let mut bytes = Vec::new();
        let mut expected = Vec::new();
        let mut density = Density::new(64 * 64);
        let mut exposure = 1.0;
        for frame in 0..7 {
            let pos = [
                [frame as f64 / 10.0, -0.5, 0.0],
                [0.0, frame as f64 / 8.0, 0.5],
            ];
            for p in pos {
                for value in p {
                    bytes.extend(value.to_ne_bytes());
                }
            }
            let points = pos.map(|p| projection.point(p));
            density.render(&points, &[0, 1], [0.0; 2], 1.0, &opts);
            if frame == 0 {
                exposure = density.exposure();
            }
            let mut pixels = vec![0; 64 * 64 * 3];
            density.rgb(&mut pixels, 64, 64, exposure);
            expected.extend(pixels);
        }
        let fixture = TestFile::new(&bytes)?;
        let mut actual = Vec::new();
        render_pipeline(
            fixture.reader(48, Backend::Auto)?,
            &mut actual,
            &projection,
            &[0, 1],
            7,
            [0.0; 2],
            1.0,
            &opts,
        )?;
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn pipeline_propagates_io_failures_and_cancels_workers() -> Result<()> {
        struct BrokenEncoder;
        impl Write for BrokenEncoder {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "encoder stopped"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let opts = Options::parse(["--width", "16", "--height", "16"].map(str::to_owned))?.unwrap();
        let projection = Projection::new(0.0, 0.0);
        let fixture = TestFile::new(&[0; 48 * 8])?;
        let error = render_pipeline(
            fixture.reader(48, Backend::Threaded)?,
            BrokenEncoder,
            &projection,
            &[0, 1],
            8,
            [0.0; 2],
            1.0,
            &opts,
        )
        .unwrap_err();
        assert!(error.to_string().contains("encoder stopped"), "{error}");
        // Asking for a ninth frame simulates a file truncated after the snapshot.
        assert!(render_pipeline(
            fixture.reader(48, Backend::Threaded)?,
            io::sink(),
            &projection,
            &[0, 1],
            9,
            [0.0; 2],
            1.0,
            &opts
        )
        .is_err());
        let bad = TestFile::new(&f64::NAN.to_ne_bytes().repeat(6))?;
        let error = render_pipeline(
            bad.reader(48, Backend::Threaded)?,
            io::sink(),
            &projection,
            &[0, 1],
            1,
            [0.0; 2],
            1.0,
            &opts,
        )
        .unwrap_err();
        assert!(error.to_string().contains("finite"), "{error}");
        Ok(())
    }

    #[test]
    fn incomplete_simulation_uses_complete_frames_and_exact_bounds() -> Result<()> {
        let path = env::temp_dir().join(format!("bhut-video-test-{}.bin", std::process::id()));
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)?;
        // Ensure the test fixture is removed even if an assertion fails.
        struct Fixture(std::path::PathBuf);
        impl Drop for Fixture {
            fn drop(&mut self) {
                let _ = std::fs::remove_file(&self.0);
            }
        }
        let _fixture = Fixture(path.clone());
        file.write_all(&2i32.to_ne_bytes())?;
        file.write_all(&5i32.to_ne_bytes())?;
        for value in [
            1.0_f64, 1.0, 0.0, -2.0, -4.0, 0.0, 2.0, 4.0, 0.0, -3.0, -6.0, 0.0, 3.0, 6.0,
        ] {
            file.write_all(&value.to_ne_bytes())?;
        }
        file.write_all(b"partial")?;
        drop(file);
        let mut sim = Simulation::open(path.to_str().unwrap())?;
        assert_eq!(sim.frames, 2);
        let (centre, half) = fixed_view(&mut sim, &Projection::new(0.0, 0.0), 0.95)?;
        assert_eq!(centre, [0.0, 0.0]);
        assert!((half - 6.3).abs() < 1e-12);
        sim.rewind()?;
        assert_eq!(sim.read_frame()?, [[0.0, -2.0, -4.0], [0.0, 2.0, 4.0]]);
        // Nonfinite input is rejected before it can corrupt a camera quantile.
        let mut file = std::fs::OpenOptions::new().write(true).open(&path)?;
        file.seek(SeekFrom::Start(24))?;
        file.write_all(&f64::NAN.to_ne_bytes())?;
        sim.rewind()?;
        assert!(sim.read_frame().is_err());
        Ok(())
    }
}
