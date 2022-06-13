#!/bin/bash
# #SBATCH -p test7d
#SBATCH -A special00005
#SBATCH -J interIsoFoam
#SBATCH -n 1
#SBATCH --mem-per-cpu=6000
#SBATCH -t 23:00:00

srun interIsoFoam
