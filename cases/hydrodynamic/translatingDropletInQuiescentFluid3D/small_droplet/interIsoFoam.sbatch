#!/bin/bash
#SBATCH -J interIsoFoam 
#SBATCH -A special00004
#SBATCH -n 1 
#SBATCH --mem-per-cpu=6000 
#SBATCH -t 23:00:00     

# -------------------------------

srun interIsoFoam 
# -no-density-equation
