#!/usr/bin/bash
#SBATCH --account special00005
#SBATCH --ntasks 16
#SBATCH --mem-per-cpu 1500
#SBATCH --time 24:00:00
#SBATCH --job-name interIsoFoam-stationaryDroplet2D
#SBATCH --mail-type=NONE
srun interIsoFoam -parallel
