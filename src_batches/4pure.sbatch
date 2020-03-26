#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
module load openmpi
mpirun -n 112 ./nbody_pure 100352 10 0.1 1