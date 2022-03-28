#!/usr/bin/bash
#SBATCH --account special00005
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu 1000
#SBATCH --time 10:00:00
#SBATCH --job-name interIsoFoam-convection
srun interIsoFoam
