#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
module load openmpi
mpirun -n 1792 ./nbody 786688 800 0.2 1