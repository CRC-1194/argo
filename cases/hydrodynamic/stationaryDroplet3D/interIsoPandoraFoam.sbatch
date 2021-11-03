#!/bin/bash
#SBATCH -A special00004 
#SBATCH -J interIsoPandoraFoam  
#SBATCH -n 8 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 01:00:00     

srun interIsoPandoraFoam -parallel
