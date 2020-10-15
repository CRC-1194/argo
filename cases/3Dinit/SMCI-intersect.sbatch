#!/bin/bash
#SBATCH -J SMCI 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tm38tomy
#SBATCH -n 1      
#SBATCH -c 1      
#SBATCH -C mpi      
#SBATCH --mem-per-cpu=7000   
#SBATCH -t 0:15:00  
#SBATCH -A project01456

# -------------------------------

srun voFoamTestSurfaceCellIntersectMeshes -surfaceFile ../toolCase/sphere.stl -fixNormals 
