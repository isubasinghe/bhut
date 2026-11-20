//! Render the simulator's native-endian binary output, streaming RGB to FFmpeg.
//! The original Matplotlib renderer is make_video.py.
mod brightness;
mod camera_cache;
mod frame_io;
mod gpu_density;
use camera_cache::CameraCache;
use frame_io::{Backend, FrameReader};
use rayon::prelude::*;
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
    "Make a video from simulation.bin (requires FFmpeg with drawtext and the selected encoder).
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
  --encoder NAME            software or nvenc (NVIDIA GPU; default: software)
  --preset NAME             Software: fast; NVENC: p4 (p1 fastest, p7 slowest)
  --crf NUMBER              Software quality, 0–51 (default: 18)
  --cq NUMBER               NVENC target quality, 1–51 (default: 18)
  --io-backend NAME         auto, uring (Linux), or threaded (default: auto)
  --render-threads INTEGER  Density workers (0 = automatic, default: 0)
  --no-camera-cache         Recompute bounds without reading/writing the sidecar
  --renderer NAME          cpu or gpu (Vulkan compute; default: cpu)
  --brightness MODE        fast lookup or exact math (default: fast)
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
    nvenc: bool,
    cq: u8,
    render_threads: usize,
    camera_cache: bool,
    gpu: bool,
    fast_brightness: bool,
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
            nvenc: false,
            cq: 18,
            render_threads: 0,
            camera_cache: true,
            gpu: false,
            fast_brightness: true,
        };
        let mut args = args.into_iter();
        let mut positional = false;
        let (mut preset_given, mut crf_given, mut cq_given) = (false, false, false);
        while let Some(arg) = args.next() {
            if arg == "--help" || arg == "-h" {
                return Ok(None);
            }
            if arg == "--no-camera-cache" {
                options.camera_cache = false;
                continue;
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
                    | "--encoder"
                    | "--cq"
                    | "--render-threads"
                    | "--renderer"
                    | "--brightness"
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
                "--render-threads" => options.render_threads = value.parse()?,
                "--renderer" => {
                    options.gpu = match value.as_str() {
                        "cpu" => false,
                        "gpu" => true,
                        _ => return Err("renderer must be cpu or gpu".into()),
                    }
                }
                "--brightness" => {
                    options.fast_brightness = match value.as_str() {
                        "fast" => true,
                        "exact" => false,
                        _ => return Err("brightness must be fast or exact".into()),
                    }
                }
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
                "--preset" => {
                    options.preset = value;
                    preset_given = true;
                }
                "--crf" => {
                    options.crf = value.parse()?;
                    crf_given = true;
                }
                "--cq" => {
                    options.cq = value.parse()?;
                    cq_given = true;
                }
                "--encoder" => {
                    options.nvenc = match value.as_str() {
                        "software" => false,
                        "nvenc" => true,
                        _ => return Err("encoder must be software or nvenc".into()),
                    }
                }
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
        if options.gpu && options.particles {
            return Err("GPU rendering currently requires --render-mode density".into());
        }
        if options.nvenc {
            if !preset_given {
                options.preset = "p4".into();
            }
            if crf_given {
                return Err("NVENC uses --cq, not --crf".into());
            }
            if !(1..=51).contains(&options.cq) {
                return Err("NVENC CQ must be in 1–51 (lower is better)".into());
            }
            if !matches!(
                options.preset.as_str(),
                "p1" | "p2" | "p3" | "p4" | "p5" | "p6" | "p7"
            ) {
                return Err("NVENC presets are p1 (fastest) through p7 (slowest)".into());
            }
        } else if cq_given {
            return Err("--cq requires --encoder nvenc; use --crf for software encoding".into());
        } else if !matches!(
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

    fn encoder_name(&self) -> &'static str {
        match (self.nvenc, self.h265) {
            (false, true) => "libx265",
            (false, false) => "libx264",
            (true, true) => "hevc_nvenc",
            (true, false) => "h264_nvenc",
        }
    }

    fn encoder_args(&self) -> Vec<String> {
        let mut args: Vec<String> = ["-c:v", self.encoder_name(), "-preset", &self.preset]
            .into_iter()
            .map(str::to_owned)
            .collect();
        if self.nvenc {
            args.extend([
                "-rc".into(),
                "vbr".into(),
                "-cq".into(),
                self.cq.to_string(),
                "-b:v".into(),
                "0".into(),
            ]);
        } else {
            args.extend(["-crf".into(), self.crf.to_string()]);
            if self.h265 {
                args.extend(["-x265-params".into(), "log-level=error".into()]);
            }
        }
        args
    }
}

fn check_gpu_encoder(opts: &Options) -> Result<()> {
    if !opts.nvenc {
        return Ok(());
    }
    // Exercise the hardware at the requested size before scanning a large input or
    // touching the output. Listing FFmpeg encoders alone does not test the driver.
    eprintln!("Checking {} on the GPU...", opts.encoder_name());
    let probe = Command::new("ffmpeg")
        .args(["-hide_banner", "-loglevel", "error", "-f", "lavfi", "-i"])
        .arg(format!(
            "color=size={}x{}:rate={}",
            opts.width, opts.height, opts.fps
        ))
        .args(["-frames:v", "1", "-pix_fmt", "yuv420p"])
        .args(opts.encoder_args())
        .args(["-f", "null", "-"])
        .stdin(Stdio::null())
        .output()
        .map_err(|error| format!("could not start FFmpeg for GPU check: {error}"))?;
    if !probe.status.success() {
        return Err(format!("{} GPU check failed. Check NVIDIA driver/device access and FFmpeg NVENC support, or use --encoder software.\n{}",
            opts.encoder_name(), String::from_utf8_lossy(&probe.stderr).trim()).into());
    }
    Ok(())
}

struct Simulation {
    file: File,
    n: usize,
    frames: usize,
    offset: u64,
    bytes: Vec<u8>,
    positions: Vec<[f64; 3]>,
    signature: Option<String>,
}

impl Simulation {
    fn open(path: &str) -> Result<Self> {
        let mut file = File::open(path)?;
        let signature = camera_cache::signature(&file).ok();
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
            signature,
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

fn cached_fixed_view(
    sim: &mut Simulation,
    projection: &Projection,
    opts: &Options,
) -> Result<camera_cache::Bounds> {
    let cache = if opts.camera_cache {
        sim.signature.as_ref().and_then(|stamp| {
            match CameraCache::new(
                Path::new(&opts.input),
                stamp,
                sim.n,
                sim.frames,
                [opts.elev, opts.azim, opts.fraction],
            ) {
                Ok(cache) => Some(cache),
                Err(error) => {
                    eprintln!("Camera cache disabled: {error}");
                    None
                }
            }
        })
    } else {
        None
    };
    if let Some(cache) = &cache {
        if cache.unchanged(&sim.file) {
            if let Some(bounds) = cache.load() {
                if cache.unchanged(&sim.file) {
                    eprintln!("Using cached camera bounds: {}", cache.path.display());
                    return Ok(bounds);
                }
            }
        }
    }
    let bounds = fixed_view(sim, projection, opts.fraction)?;
    if let Some(cache) = &cache {
        if cache.unchanged(&sim.file) {
            if let Err(error) = cache.store(bounds) {
                eprintln!("Could not save camera cache (rendering continues): {error}");
            }
        } else {
            eprintln!("Input changed during camera setup; bounds were not cached");
        }
    }
    Ok(bounds)
}

fn gaussian_kernel() -> [f32; 9] {
    let mut kernel: [f32; 9] =
        std::array::from_fn(|i| (-0.5 * ((i as f32 - 4.0) / 1.4).powi(2)).exp());
    let sum: f32 = kernel.iter().sum();
    kernel.iter_mut().for_each(|v| *v /= sum);
    kernel
}

const TILE: usize = 32;
fn active_ranges(mask: Option<&[u8]>, width: usize, mut visit: impl FnMut(usize, usize)) {
    let Some(mask) = mask else {
        visit(0, width);
        return;
    };
    let mut tile = 0;
    while tile < mask.len() {
        if mask[tile] == 0 {
            tile += 1;
            continue;
        }
        let start = tile;
        while tile < mask.len() && mask[tile] != 0 {
            tile += 1;
        }
        visit(start * TILE, (tile * TILE).min(width));
    }
}

#[cfg(test)]
fn blur(layer: &mut [f32], scratch: &mut [f32], width: usize, height: usize) {
    blur_tiles(layer, scratch, width, height, None);
}

fn blur_tiles(
    layer: &mut [f32],
    scratch: &mut [f32],
    width: usize,
    height: usize,
    masks: Option<(&[u8], &[u8])>,
) {
    let kernel = gaussian_kernel();
    let cols = width.div_ceil(TILE);
    // Finish each row while it is in cache instead of scanning the whole 4K
    // image nine times. Keep the accumulation order and zero padding unchanged.
    scratch
        .par_chunks_exact_mut(width)
        .enumerate()
        .with_min_len(16)
        .for_each(|(y, out)| {
            out.fill(0.0);
            let mask =
                masks.map(|(vertical, _)| &vertical[(y / TILE) * cols..(y / TILE + 1) * cols]);
            active_ranges(mask, width, |start, end| {
                for (k, &weight) in kernel.iter().enumerate() {
                    let source = y as isize + k as isize - 4;
                    if source < 0 || source >= height as isize {
                        continue;
                    }
                    let row =
                        &layer[source as usize * width + start..source as usize * width + end];
                    for (out, &value) in out[start..end].iter_mut().zip(row) {
                        *out += weight * value;
                    }
                }
            });
        });
    layer
        .par_chunks_exact_mut(width)
        .enumerate()
        .with_min_len(16)
        .for_each(|(y, out)| {
            out.fill(0.0);
            let mask =
                masks.map(|(_, final_mask)| &final_mask[(y / TILE) * cols..(y / TILE + 1) * cols]);
            active_ranges(mask, width, |x0, x1| {
                for (k, &weight) in kernel.iter().enumerate() {
                    let shift = k as isize - 4;
                    let start = ((-shift).max(0) as usize).max(x0);
                    let end =
                        ((width as isize - shift).min(width as isize).max(0) as usize).min(x1);
                    if start >= end {
                        continue;
                    }
                    let source = (y * width + start) as isize + shift;
                    let input = &scratch[source as usize..source as usize + end - start];
                    for (out, &value) in out[start..end].iter_mut().zip(input) {
                        *out += weight * value;
                    }
                }
            });
        });
}

struct Density {
    layers: [Vec<f32>; 2],
    scratch: Vec<f32>,
    vertical: [Vec<u8>; 2],
    active: [Vec<u8>; 2],
    combined: Vec<u8>,
}
impl Density {
    fn new(size: usize) -> Self {
        Self {
            layers: [vec![0.0; size], vec![0.0; size]],
            scratch: vec![0.0; size],
            vertical: [Vec::new(), Vec::new()],
            active: [Vec::new(), Vec::new()],
            combined: Vec::new(),
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
        let cols = w.div_ceil(TILE);
        let tiles = cols * h.div_ceil(TILE);
        for mask in self.vertical.iter_mut().chain(self.active.iter_mut()) {
            mask.resize(tiles, 0);
            mask.fill(0);
        }
        self.layers.iter_mut().for_each(|layer| layer.fill(0.0));
        for (p, &group) in points.iter().zip(groups) {
            let x = (p[0] - centre[0]) * scale + w as f64 / 2.0;
            let y = (p[1] - centre[1]) * scale + h as f64 / 2.0;
            if x >= 0.0 && x <= w as f64 && y >= 0.0 && y <= h as f64 {
                // As in histogram2d, include the upper bound in the last bin.
                let (x, y) = ((x as usize).min(w - 1), (y as usize).min(h - 1));
                self.layers[group][y * w + x] += 1.0;
                // Only tiles touched by the exact four-pixel kernel halo are active.
                for ty in y.saturating_sub(4) / TILE..=(y + 4).min(h - 1) / TILE {
                    self.vertical[group][ty * cols + x / TILE] = 1;
                    for tx in x.saturating_sub(4) / TILE..=(x + 4).min(w - 1) / TILE {
                        self.active[group][ty * cols + tx] = 1;
                    }
                }
            }
        }
        for group in 0..2 {
            blur_tiles(
                &mut self.layers[group],
                &mut self.scratch,
                w,
                h,
                Some((&self.vertical[group], &self.active[group])),
            );
        }
        self.combined.resize(tiles, 0);
        for (i, value) in self.combined.iter_mut().enumerate() {
            *value = self.active[0][i] | self.active[1][i];
        }
    }
    fn exposure(&self) -> f32 {
        self.layers[0]
            .iter()
            .zip(&self.layers[1])
            .map(|(a, b)| a + b)
            .fold(1.0, f32::max)
    }
    fn rgb(&self, pixels: &mut [u8], width: usize, height: usize, exposure: f32, fast: bool) {
        let curve = brightness::Brightness::new(exposure, fast);
        let cols = width.div_ceil(TILE);
        pixels
            .par_chunks_exact_mut(width * 3)
            .enumerate()
            .with_min_len(16)
            .for_each(|(screen_y, row)| {
                let y = height - 1 - screen_y;
                row.fill(0);
                active_ranges(
                    Some(&self.combined[(y / TILE) * cols..(y / TILE + 1) * cols]),
                    width,
                    |start, end| {
                        for (local_x, pixel) in
                            row[start * 3..end * 3].chunks_exact_mut(3).enumerate()
                        {
                            let x = start + local_x;
                            let a = self.layers[0][y * width + x];
                            let b = self.layers[1][y * width + x];
                            let total = a + b;
                            if total == 0.0 {
                                pixel.fill(0);
                                continue;
                            }
                            let brightness = curve.value(total);
                            for c in 0..3 {
                                pixel[c] = (255.0 * (a * COLOURS[0][c] + b * COLOURS[1][c])
                                    / total.max(1e-12)
                                    * brightness) as u8;
                            }
                        }
                    },
                );
            });
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

#[cfg(test)]
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
    let gpu = if opts.gpu {
        Some(gpu_density::GpuDensity::new(opts, groups.len())?)
    } else {
        None
    };
    render_pipeline_with_gpu(
        input, output, projection, groups, frames, centre, half, opts, gpu,
    )
}

fn render_pipeline_with_gpu(
    input: FrameReader,
    output: impl Write + Send,
    projection: &Projection,
    groups: &[usize],
    frames: usize,
    centre: [f64; 2],
    half: f64,
    opts: &Options,
    mut gpu: Option<gpu_density::GpuDensity>,
) -> Result<()> {
    // Two reusable buffers at each boundary bound memory and apply backpressure.
    // Reading/projecting and writing to FFmpeg happen independently of rendering.
    let pool = if opts.particles || opts.gpu {
        None
    } else {
        let workers = if opts.render_threads == 0 {
            thread::available_parallelism().map_or(1, |n| n.get().min(8))
        } else {
            opts.render_threads
        };
        eprintln!("Density rendering: {workers} workers");
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()?,
        )
    };
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
            let mut density = if opts.particles || opts.gpu {
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
                if let Some(gpu) = &mut gpu {
                    gpu.render(
                        &points,
                        groups,
                        centre,
                        half,
                        &mut pixels,
                        opts.fast_brightness,
                    )?;
                } else if let Some(density) = &mut density {
                    pool.as_ref().unwrap().install(|| {
                        density.render(&points, groups, centre, half, opts);
                        if frame == 0 {
                            exposure = density.exposure();
                        }
                        density.rgb(
                            &mut pixels,
                            opts.width,
                            opts.height,
                            exposure,
                            opts.fast_brightness,
                        );
                    });
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
    gpu: Option<gpu_density::GpuDensity>,
) -> Result<()> {
    let caption = if opts.particles {
        "Individual particles · fixed camera"
    } else {
        "Projected particle density · logarithmic brightness · fixed camera"
    };
    let filter = format!("setsar=1,drawtext=text='Galaxy collision  |  step %{{n}}':fontcolor=white:fontsize=27:x=w*0.04:y=h*0.035,drawtext=text='Initial group 1':fontcolor=0x40b3ff:fontsize=21:x=w*0.04:y=h*0.075,drawtext=text='Initial group 2':fontcolor=0xff8c33:fontsize=21:x=w*0.24:y=h*0.075,drawtext=text='{caption}':fontcolor=0xa6a6a6:fontsize=19:x=w*0.04:y=h*0.975-th");
    let input = FrameReader::new(sim.file.try_clone()?, sim.offset, sim.n * 24, opts.backend)?;
    eprintln!(
        "Input: {}; encoder: {}; preset: {}",
        input.name(),
        opts.encoder_name(),
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
            "-pix_fmt",
            "yuv420p",
            "-vf",
            &filter,
            "-movflags",
            "+faststart",
        ])
        .args(opts.encoder_args());
    if opts.h265 {
        command.args(["-tag:v", "hvc1"]);
    }
    let mut child = command
        .arg(&opts.output)
        .stdin(Stdio::piped())
        .spawn()
        .map_err(|e| format!("could not start FFmpeg (enter nix develop): {e}"))?;
    let pipe = child.stdin.take().expect("piped FFmpeg stdin");
    let started = Instant::now();
    let result = render_pipeline_with_gpu(
        input, pipe, projection, groups, sim.frames, centre, half, opts, gpu,
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
    check_gpu_encoder(&opts)?;
    let gpu = if opts.gpu {
        Some(gpu_density::GpuDensity::new(&opts, sim.n)?)
    } else {
        None
    };
    sim.rewind()?;
    let groups = initial_groups(sim.read_frame()?);
    let projection = Projection::new(opts.elev, opts.azim);
    let (centre, half) = cached_fixed_view(&mut sim, &projection, &opts)?;
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
    encode(&mut sim, &projection, &groups, centre, half, &opts, gpu)?;
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
    fn parallel_density_preserves_every_pixel() -> Result<()> {
        let opts =
            Options::parse(["--width", "258", "--height", "130"].map(str::to_owned))?.unwrap();
        let points: Vec<_> = (0..2000)
            .map(|i| {
                [
                    ((i * 37 % 501) as f64 - 250.0) / 200.0,
                    ((i * 73 % 401) as f64 - 200.0) / 200.0,
                ]
            })
            .collect();
        let groups: Vec<_> = (0..points.len()).map(|i| i % 2).collect();
        let render = |workers| -> Result<_> {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()?;
            Ok(pool.install(|| {
                let mut density = Density::new(opts.width * opts.height);
                density.render(&points, &groups, [0.0; 2], 1.0, &opts);
                let mut pixels = vec![0; opts.width * opts.height * 3];
                density.rgb(
                    &mut pixels,
                    opts.width,
                    opts.height,
                    density.exposure(),
                    opts.fast_brightness,
                );
                (density.layers, pixels)
            }))
        };
        let serial = render(1)?;
        for workers in [2, 4, 8] {
            assert_eq!(render(workers)?, serial);
        }
        Ok(())
    }

    #[test]
    fn sparse_tiles_match_dense_blur_at_boundaries_and_clear_old_frames() -> Result<()> {
        let opts = Options::parse(["--width", "70", "--height", "66"].map(str::to_owned))?.unwrap();
        let mut density = Density::new(70 * 66);
        let mut reference = [vec![0.0; 70 * 66], vec![0.0; 70 * 66]];
        let bins = [
            (0, 0),
            (31, 31),
            (32, 32),
            (35, 35),
            (63, 63),
            (64, 64),
            (69, 65),
        ];
        let groups: Vec<_> = (0..bins.len()).map(|i| i % 2).collect();
        let points: Vec<_> = bins
            .iter()
            .map(|&(x, y)| {
                [
                    (x as f64 + 0.5 - 35.0) / 33.0,
                    (y as f64 + 0.5 - 33.0) / 33.0,
                ]
            })
            .collect();
        for ((x, y), group) in bins.into_iter().zip(&groups) {
            reference[*group][y * 70 + x] += 1.0;
        }
        for layer in &mut reference {
            blur(layer, &mut vec![0.0; 70 * 66], 70, 66);
        }
        density.render(&points, &groups, [0.0; 2], 1.0, &opts);
        assert_eq!(density.layers, reference);
        let mut exact = vec![0; 70 * 66 * 3];
        let mut fast = exact.clone();
        density.rgb(&mut exact, 70, 66, density.exposure(), false);
        density.rgb(&mut fast, 70, 66, density.exposure(), true);
        assert!(exact.iter().zip(&fast).all(|(&a, &b)| a.abs_diff(b) <= 1));
        density.render(&[], &[], [0.0; 2], 1.0, &opts);
        density.rgb(&mut fast, 70, 66, 1.0, true);
        assert!(density.layers.iter().flatten().all(|&v| v == 0.0));
        assert!(fast.iter().all(|&v| v == 0));
        Ok(())
    }

    #[test]
    #[ignore = "requires a Vulkan GPU; run explicitly outside GPU-restricted sandboxes"]
    fn gpu_density_matches_cpu() -> Result<()> {
        let opts = Options::parse(
            ["--width", "258", "--height", "130", "--renderer", "gpu"].map(str::to_owned),
        )?
        .unwrap();
        let points: Vec<_> = (0..2000)
            .map(|i| {
                [
                    ((i * 37 % 501) as f64 - 250.0) / 200.0,
                    ((i * 73 % 401) as f64 - 200.0) / 200.0,
                ]
            })
            .collect();
        let groups: Vec<_> = (0..points.len()).map(|i| i % 2).collect();
        for fast in [false, true] {
            let mut gpu = gpu_density::GpuDensity::new(&opts, points.len())?;
            let mut cpu = Density::new(opts.width * opts.height);
            let mut rgb = vec![0; opts.width * opts.height * 3];
            let mut gpu_rgb = vec![0; opts.width * opts.height * 3];
            let mut exposure = 1.0;
            for (frame, points) in [points.as_slice(), &[], points.as_slice()]
                .into_iter()
                .enumerate()
            {
                cpu.render(points, &groups, [0.0; 2], 1.0, &opts);
                if frame == 0 {
                    exposure = cpu.exposure();
                }
                cpu.rgb(&mut rgb, opts.width, opts.height, exposure, fast);
                gpu.render(points, &groups, [0.0; 2], 1.0, &mut gpu_rgb, fast)?;
                let worst = rgb
                    .chunks_exact(3)
                    .zip(gpu_rgb.chunks_exact(3))
                    .flat_map(|(cpu, gpu)| (0..3).map(move |i| cpu[i].abs_diff(gpu[i])))
                    .max()
                    .unwrap();
                assert!(
                    worst <= 1,
                    "GPU pixel error {worst}; fast={fast}, frame={frame}"
                );
            }
        }
        Ok(())
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
            ["--encoder", "other"],
            ["--cq", "18"],
            ["--preset", "p4"],
            ["--render-threads", "-1"],
            ["--renderer", "other"],
            ["--brightness", "other"],
        ] {
            assert!(Options::parse(args.map(str::to_owned)).is_err(), "{args:?}");
        }
    }

    #[test]
    fn encoder_options_keep_software_and_nvenc_settings_separate() -> Result<()> {
        let software = Options::parse(Vec::<String>::new())?.unwrap();
        assert_eq!(software.render_threads, 0);
        assert!(software.camera_cache);
        assert!(
            !Options::parse(["--no-camera-cache".to_owned()])?
                .unwrap()
                .camera_cache
        );
        assert_eq!(software.encoder_name(), "libx265");
        assert_eq!(software.preset, "fast");
        assert!(software.encoder_args().contains(&"-crf".into()));
        assert!(!software.encoder_args().contains(&"-cq".into()));
        let gpu = Options::parse(["--encoder", "nvenc"].map(str::to_owned))?.unwrap();
        assert_eq!(gpu.encoder_name(), "hevc_nvenc");
        assert_eq!(gpu.preset, "p4");
        assert_eq!(
            gpu.encoder_args(),
            [
                "-c:v",
                "hevc_nvenc",
                "-preset",
                "p4",
                "-rc",
                "vbr",
                "-cq",
                "18",
                "-b:v",
                "0"
            ]
        );
        for codec in ["h264", "h265"] {
            let gpu = Options::parse(
                [
                    "--preset",
                    "p7",
                    "--cq",
                    "23",
                    "--codec",
                    codec,
                    "--encoder",
                    "nvenc",
                ]
                .map(str::to_owned),
            )?
            .unwrap();
            assert_eq!(
                gpu.encoder_name(),
                if codec == "h264" {
                    "h264_nvenc"
                } else {
                    "hevc_nvenc"
                }
            );
            assert_eq!(gpu.preset, "p7");
            assert_eq!(gpu.cq, 23);
            assert!(!gpu
                .encoder_args()
                .iter()
                .any(|arg| arg == "-crf" || arg == "-x265-params"));
        }
        for (option, value) in [
            ("--crf", "18"),
            ("--preset", "veryfast"),
            ("--cq", "0"),
            ("--cq", "52"),
        ] {
            assert!(
                Options::parse(["--encoder", "nvenc", option, value].map(str::to_owned)).is_err()
            );
        }
        Ok(())
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
            let mut sidecar = self.0.as_os_str().to_os_string();
            sidecar.push(".camera-cache");
            let _ = std::fs::remove_file(sidecar);
        }
    }

    fn camera_fixture() -> Result<(TestFile, Options)> {
        let mut bytes = Vec::new();
        bytes.extend(2i32.to_ne_bytes());
        bytes.extend(2i32.to_ne_bytes());
        for value in [
            1.0_f64, 1.0, 0.0, -2.0, -4.0, 0.0, 2.0, 4.0, 0.0, -3.0, -6.0, 0.0, 3.0, 6.0,
        ] {
            bytes.extend(value.to_ne_bytes());
        }
        let file = TestFile::new(&bytes)?;
        let opts = Options::parse([file.0.to_str().unwrap().to_owned()])?.unwrap();
        Ok((file, opts))
    }

    #[test]
    fn camera_cache_roundtrip_and_invalidation() -> Result<()> {
        let (fixture, opts) = camera_fixture()?;
        let mut sim = Simulation::open(&opts.input)?;
        let projection = Projection::new(opts.elev, opts.azim);
        let expected = fixed_view(&mut sim, &projection, opts.fraction)?;
        assert_eq!(cached_fixed_view(&mut sim, &projection, &opts)?, expected);
        let stamp = sim.signature.as_ref().unwrap().clone();
        let key = |angles, n, frames| CameraCache::new(&fixture.0, &stamp, n, frames, angles);
        let angles = [opts.elev, opts.azim, opts.fraction];
        let cache = key(angles, sim.n, sim.frames)?;
        assert_eq!(cache.load(), Some(expected));
        // A hit performs no reads from the simulation after opening the file.
        sim.file.seek(SeekFrom::End(0))?;
        let before = sim.file.stream_position()?;
        assert_eq!(cached_fixed_view(&mut sim, &projection, &opts)?, expected);
        assert_eq!(sim.file.stream_position()?, before);
        for i in 0..3 {
            let mut changed = angles;
            changed[i] += 0.01;
            assert!(key(changed, sim.n, sim.frames)?.load().is_none());
        }
        assert!(key(angles, sim.n + 1, sim.frames)?.load().is_none());
        assert!(key(angles, sim.n, sim.frames + 1)?.load().is_none());
        // Corruption that still represents a finite float must also be detected.
        let mut bytes = std::fs::read(&cache.path)?;
        let end = bytes.len();
        bytes[end - 32] ^= 1;
        std::fs::write(&cache.path, bytes)?;
        assert!(cache.load().is_none());
        assert_eq!(cached_fixed_view(&mut sim, &projection, &opts)?, expected);
        std::fs::write(&cache.path, b"truncated cache")?;
        assert!(cache.load().is_none());
        // Same-size changes invalidate the input signature.
        thread::sleep(std::time::Duration::from_millis(10));
        let mut edit = std::fs::OpenOptions::new().write(true).open(&fixture.0)?;
        #[cfg(unix)]
        let previous_modified = edit.metadata()?.modified()?;
        edit.seek(SeekFrom::Start(24))?;
        edit.write_all(&1.0_f64.to_ne_bytes())?;
        assert!(!cache.unchanged(&sim.file));
        #[cfg(unix)]
        {
            edit.set_times(std::fs::FileTimes::new().set_modified(previous_modified))?;
            assert!(!cache.unchanged(&sim.file));
        }
        // Replacement and growth are independently represented in the key.
        let replacement = TestFile::new(&std::fs::read(&fixture.0)?)?;
        assert_ne!(
            camera_cache::signature(&File::open(&replacement.0)?)?,
            stamp
        );
        Ok(())
    }

    #[test]
    fn camera_cache_is_optional_and_never_publishes_changed_input() -> Result<()> {
        let (fixture, mut opts) = camera_fixture()?;
        let projection = Projection::new(opts.elev, opts.azim);
        let mut sim = Simulation::open(&opts.input)?;
        let cache = CameraCache::new(
            &fixture.0,
            sim.signature.as_ref().unwrap(),
            sim.n,
            sim.frames,
            [opts.elev, opts.azim, opts.fraction],
        )?;
        opts.camera_cache = false;
        let expected = cached_fixed_view(&mut sim, &projection, &opts)?;
        assert!(!cache.path.exists());
        opts.camera_cache = true;
        // A non-writable cache location cannot prevent rendering.
        std::fs::create_dir(&cache.path)?;
        let result = cached_fixed_view(&mut sim, &projection, &opts);
        std::fs::remove_dir(&cache.path)?;
        assert_eq!(result?, expected);
        // Simulate a producer appending after the initial frame-count snapshot.
        std::fs::OpenOptions::new()
            .append(true)
            .open(&fixture.0)?
            .write_all(&[0; 48])?;
        assert_eq!(cached_fixed_view(&mut sim, &projection, &opts)?, expected);
        assert!(!cache.path.exists());
        Ok(())
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
            density.rgb(&mut pixels, 64, 64, exposure, opts.fast_brightness);
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
