#!/bin/bash
#SBATCH -J CCMI 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tm38tomy
#SBATCH -n 1      # Process 
#SBATCH -c 1      # Core 
#SBATCH -C mpi    # Core 
#SBATCH --mem-per-cpu=10000   # Memory 
#SBATCH -t 4:00:00 # Hours, Minutes, Seconds: runtime 

# -------------------------------

srun voFoamTestCellCellIntersectMeshes 
