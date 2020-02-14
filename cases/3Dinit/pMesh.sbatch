#!/bin/bash
#SBATCH -J bMesh 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tm38tomy
#SBATCH -n 1      # Processes 
#SBATCH -C mpi      
#SBATCH --mem-per-cpu=10000   # Hauptspeicher in MByte pro MPI-Task
#SBATCH -t 00:10:00     # in Stunden, Minuten und Sekunden, oder '#SBATCH -t 10' - nur Minuten

# -------------------------------
srun  pMesh 
