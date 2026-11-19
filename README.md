# A Distributed-Memory simulation of the Andromeda-Milky Way collision.
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
python3 make_video.py simulation.bin -o simulation.mp4 --fps 10
```

Output defaults to 3840 × 2160 (4K UHD). Use `--width` and `--height` to change
the resolution, for example `--width 2560 --height 2560` for a square video.
Both dimensions must be positive even numbers.

The video renders a density grid matching the output resolution (3840 × 2160 by
default), with a small Gaussian blur (sigma 1.4 output pixels) and logarithmic
brightness. Each grid cell maps to one output pixel, without interpolated
enlargement. Higher resolutions require more rendering time and memory.
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
with `--elev 25 --azim 45` (angles in degrees). Rendering requires NumPy,
Matplotlib, and FFmpeg, all included in the Nix shell.

## Caveats
Well I do not simulate general relativity, this is pretty simple Newtonian physics, which is mostly fine
apart from the fact there are two giant black holes at the middle of each galaxy. So this simulation is really a very rough approximation 
of what would happen, to get a more accurate solution, simulation of general relativity is required. 

## Neat facts about this implementation
I do not use Runge-Kutta! Yup, I use something called leapfrog integration which keeps a bound on the error. 
This is critical since we need to conserve momentumn (because you know Physics be like that), something Runge-Kutta cannot do.

## Grade 
20/20 




