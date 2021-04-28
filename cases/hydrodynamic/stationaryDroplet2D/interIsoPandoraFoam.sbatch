#!/bin/bash
#SBATCH -A special00005 
#SBATCH -J interIsoPandoraFoam  
#SBATCH -n 4 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 00:45:00     

srun interIsoPandoraFoam -parallel
