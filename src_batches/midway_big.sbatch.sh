#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
module load openmpi
mpirun -n 224 ./nbody 100352 800 0.2 1