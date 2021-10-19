#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=bhut-phys.out
#SBATCH --partition=snowy
#SBATCH --time=0:4:00
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=4G
#SBATCH --job-name=bhut-phys
module load gcc/8.3.0
module load openmpi/3.1.4

mpic++ main.cpp -o bhut -static-libstdc++ -static -O3 -fno-exceptions -march=native -fno-rtti -funroll-loops
mpirun -n 8 bhut

echo "=============================================================================================="
mpirun -n 16 bhut

echo "=============================================================================================="
mpirun -n 24 bhut


echo "=============================================================================================="
mpirun -n 32 bhut


echo "=============================================================================================="
mpirun -n 40 bhut


echo "=============================================================================================="
mpirun bhut
