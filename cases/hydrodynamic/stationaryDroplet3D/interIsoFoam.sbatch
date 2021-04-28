#!/bin/bash
#SBATCH -A special00005 
#SBATCH -J interIsoFoam  
#SBATCH -n 4 
#SBATCH --mem-per-cpu=2000 
#SBATCH -t 05:00:00     

srun interIsoFoam  -parallel 
