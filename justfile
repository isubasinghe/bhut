# Run these recipes in the nix shell: `nix develop`, then `just`.
# Or without the shell: `nix develop -c just`.

# Number of MPI processes
np := "12"
# Frames for each second of video
fps := "30"

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
    python3 make_video.py simulation.bin -o simulation.mp4 --fps {{fps}}

# Delete the compiled program and the outputs
clean:
    rm -rf build simulation.bin simulation.mp4
