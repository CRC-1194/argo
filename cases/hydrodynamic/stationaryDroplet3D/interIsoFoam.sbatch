#!/usr/bin/bash
#SBATCH --account special00005
#SBATCH --ntasks 8
#SBATCH --mem-per-cpu 1500
#SBATCH --time 24:00:00
#SBATCH --job-name interIsoFoam-stationaryDroplet3D
#SBATCH --mail-type=NONE
srun interIsoFoam -parallel
