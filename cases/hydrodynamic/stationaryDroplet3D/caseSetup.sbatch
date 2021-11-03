#!/bin/bash
#SBATCH -A special00004 
#SBATCH -J locInit 
#SBATCH -n 1 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 00:20:00     

srun blockMesh && setAlphaField && decomposePar -force

# Use the line below to use the accuracy driven SMCA algorithm
# for initialization of the volume fractions
#srun blockMesh && surfaceInitVolumeFraction && decomposePar -force
