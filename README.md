# A Distributed-Memory simulation of the Andromeda-Milky Way collision.

![Andromeda–Milky Way collision simulation](simulation.gif)

## What is this?
The Andromeda and Milky Way galaxies will collide in approximately 5 billion years, forming Milkdromeda or the Andromeda-Milky Way galaxy. 
This is a distributed memory simulation of this collision run on the ''Supercomputer''/HPC system at the University of Melbourne done as part of coursework for Masters. 

## How was this done?
OpenMPI was used for message passing on the supercomputer/HPC system, work was split evenly among the number of worker nodes and each worker node ran 
an implementation of the Barnes-Hut algorithm 

## I wanna read more
Thanks for your interest, I am pretty proud/happy of this implementation, my report is available [here](report.pdf)

## Render a video
Inside `nix develop`, run `just video`, or:

```sh
cargo run --release -- simulation.bin -o simulation.mp4 --fps 10
```

The Rust renderer (`make_video.rs`) streams frames directly to FFmpeg without
temporary PNGs. It defaults to H.265/HEVC (`libx265`), CRF 18, and the `fast`
encoder preset, with an `hvc1` tag for MP4 playback. `just video` uses this renderer at
30 fps (override with `just fps=10 video`). The original Python script is still
available through `just video-python` or
`python3 make_video.py simulation.bin -o simulation.mp4 --fps 10`.

Use `--codec h264` for H.264, `--preset veryfast` to favor encoding speed, or
`--preset medium` to favor compression. `--crf` sets quality (0–51, lower is
better); the same CRF does not imply identical quality between codecs. H.265
does not inherently encode faster. For example:

```sh
cargo run --release -- simulation.bin -o simulation-h265.mp4 --fps 30 --preset veryfast
```

For NVIDIA GPU encoding, run `just video-gpu` (writes `simulation.mp4`), or:

```sh
cargo run --release -- simulation.bin -o simulation-gpu.mp4 --fps 30 --encoder nvenc
```

`--encoder nvenc` selects `hevc_nvenc` for H.265, or `h264_nvenc` with
`--codec h264`. Software encoding remains the default. The GPU path defaults
to preset `p4` and CQ 18 with variable bitrate and no average bitrate target.
Use `--preset p1` for maximum encoding speed, through `--preset p7` for greater
compression effort, and `--cq 18` to set target quality (1–51; lower is better).
CQ and software CRF are different controls and their numeric values are not
directly comparable. Software preset names such as `veryfast` and `--crf` are
rejected with NVENC; `--cq` requires NVENC.

GPU encoding requires an NVENC-capable NVIDIA GPU, a compatible installed driver,
and an FFmpeg build with the selected NVENC encoder. A one-frame hardware check
at the requested resolution runs before camera scans and before opening the
output. Failure reports the driver/encoder error without silently falling back.
The density rendering, text overlays, and camera calculation still run on the
CPU; only video compression moves to the GPU. See
[NVIDIA's FFmpeg guide](https://docs.nvidia.com/video-technologies/video-codec-sdk/13.0/ffmpeg-with-nvidia-gpu/index.html).

Output defaults to 3840 × 2160 (4K UHD). Use `--width` and `--height` to change
the resolution, for example `--width 2560 --height 2560` for a square video.
Both dimensions must be positive even numbers.

The video renders a density grid matching the output resolution (3840 × 2160 by
default), with a small Gaussian blur (sigma 1.4 output pixels) and logarithmic
brightness. Each grid cell maps to one output pixel, without interpolated
enlargement. Higher resolutions require more rendering time and memory.
Blur and colour conversion run on a reusable worker pool, preserving the pixel
calculations and accumulation order. The default uses up to eight available CPU
threads. Override with `--render-threads 4`, use `--render-threads 1` for one
worker, or `--render-threads 0` for automatic selection. This setting applies to
density rendering; particle rendering is unchanged. The pool shares the density
grids rather than allocating another full image per worker.
Blue and orange track two spatial groups inferred from the initial frame; these
approximate galaxy memberships because the simulation output has no galaxy labels.
The camera and exposure stay fixed throughout the video. Exposure is set from
the first frame, so later densities above that reference saturate.

Use `--render-mode particles` to draw individual markers instead, and adjust
their diameter in pixels with `--particle-size 2`. Particles that overlap in
projection can still appear merged in dense regions in either mode.

By default, the view includes at least 95% of projected positions across all
frames. Density mode expands the view along the longer output dimension to
preserve proportions. Use `--fraction 0.99` to include more outskirts, or adjust the camera
with `--elev 25 --azim 45` (angles in degrees).

Rust retains exact camera quantiles across the complete simulation. On the first
run, three input scans use one 8-byte value per particle per complete frame.
The resulting bounds are saved beside the canonical input file, e.g.
`simulation.bin.camera-cache`. Subsequent runs with a valid cache skip all three
scans and their large temporary allocation.

The cache is invalidated by changes to the input's size, timestamps, particle or
frame count, camera elevation/azimuth, or `--fraction`. On Unix it also checks
device/inode identity and change time, detecting replacements and same-size edits
even if the modification time is restored. Resolution, codec, output path, and
worker count can change without recalculating bounds. This uses filesystem
metadata rather than rereading the entire input to hash its contents.

Cache records are versioned and checksummed, and writes use an atomic rename.
Missing, corrupt, or stale records are recomputed; a cache write failure does not
stop rendering. Bounds are not cached if the input changes during camera setup.
Use `--no-camera-cache` to bypass both reading and writing the sidecar. Cache files
are ignored by Git. A cache miss still requires the scans and their input-sized
temporary allocation before rendering can begin.

After camera setup, a bounded pipeline overlaps reading/projecting upcoming frames,
rendering, and writing completed RGB frames to FFmpeg. Linux uses the `io-uring`
crate for reads when supported; `--io-backend auto` falls back to a reader thread
if the kernel or sandbox blocks it. `--io-backend uring` requires io_uring and
reports an error if unavailable; `--io-backend threaded` selects portable reads.
Both backends run on the input worker, with one read in flight at a time and
one projected frame of lookahead. io_uring does not accelerate the blur or codec.

Rendering uses two reusable projected-frame buffers, two RGB buffers (6 bytes per
output pixel total), and density grids of about 12 bytes per pixel. Queues apply
backpressure instead of buffering the whole video. The blur processes each row
while it is in cache, preserving the original pixel calculations. The frame
counter reports frames queued; success is reported only after FFmpeg finishes.

The renderer snapshots the number of complete frames at startup. Incomplete
trailing frames and frames appended later are ignored: this is streaming file
processing, not a live tail of a running simulation. Group colours can be swapped
relative to Python because the principal axis has an arbitrary sign; particle
antialiasing and font rendering also differ slightly.

The Nix shell includes Rust, Cargo, and FFmpeg (with libx265, libx264, and drawtext), plus
NumPy and Matplotlib for the Python version. Run `cargo test` for the Rust checks,
or `cargo run --release -- --help` for all renderer options.
On a Linux host with io_uring enabled, `BHUT_TEST_URING=1 cargo test` also requires
the native reader tests to run rather than allowing an unavailable backend to skip.

## Render a GIF
After `just video`, run `just gif` inside `nix develop` to convert
`simulation.mp4` into a looping `simulation.gif`. The GIF keeps the full video
duration and aspect ratio, at 640 pixels wide and 10 frames per second.
For a smaller file, use `just --set gif_width 480 --set gif_fps 8 gif`.

Commit `simulation.gif` alongside the README to embed it with:

```markdown
![Andromeda–Milky Way collision simulation](simulation.gif)
```

## Caveats
Well I do not simulate general relativity, this is pretty simple Newtonian physics, which is mostly fine
apart from the fact there are two giant black holes at the middle of each galaxy. So this simulation is really a very rough approximation 
of what would happen, to get a more accurate solution, simulation of general relativity is required. 

## Neat facts about this implementation
I do not use Runge-Kutta! Yup, I use something called leapfrog integration which keeps a bound on the error. 
This is critical since we need to conserve momentumn (because you know Physics be like that), something Runge-Kutta cannot do.

## Grade 
20/20 




