#!/bin/bash
#SBATCH -p test24
#SBATCH -A special00005 
#SBATCH -J interIsoPandoraFoam  
#SBATCH -n 8 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 20:00:00     

srun interIsoPandoraFoam -parallel
