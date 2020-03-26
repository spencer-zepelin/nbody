#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive
module load openmpi
mpirun -n 4 ./nbody_hybrid 100352 10 0.1 28