# Run these recipes in the nix shell: `nix develop`, then `just`.
# Or without the shell: `nix develop -c just`.

# Number of MPI processes
np := "12"
# Frames for each second of video
fps := "30"
# Smaller defaults for a GIF embedded in the README
gif_fps := "10"
gif_width := "640"

bin := "build/bhut"

# Do the whole process: compile, simulate, make the video
all: build simulate video

# Compile the simulation
build:
    mkdir -p build
    mpicxx isithasubasinghe_bhut.cpp -o {{bin}} -O3 -funroll-loops -Wall

# Run the simulation on planets.txt. Writes simulation.bin
simulate: build
    mpirun -np {{np}} {{bin}}

# Make simulation.mp4 from simulation.bin
video:
    cargo run --release -- simulation.bin -o simulation.mp4 --fps {{fps}}

# Use the original Python/Matplotlib renderer
video-python:
    python3 make_video.py simulation.bin -o simulation.mp4 --fps {{fps}}

# Make a looping simulation.gif from simulation.mp4 (run `just video` first)
gif:
    ffmpeg -hide_banner -loglevel warning -y -i simulation.mp4 \
        -filter_complex "fps={{gif_fps}},scale={{gif_width}}:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=bayer:bayer_scale=5" \
        -an -loop 0 simulation.gif

# Delete the compiled program and the outputs
clean:
    rm -rf build simulation.bin simulation.mp4 simulation.gif
