#!/bin/bash
#SBATCH -J foamPerturbMesh 
#SBATCH -n 1      
#SBATCH --mem-per-cpu=3000
#SBATCH -t 00:10:00     

# -------------------------------
srun foamPerturbMesh -alpha 0.25
