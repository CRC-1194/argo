#!/bin/bash
#SBATCH -A special00005 
#SBATCH -J locInit 
#SBATCH -n 1 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 00:20:00     

srun blockMesh && setAlphaField && decomposePar -force
