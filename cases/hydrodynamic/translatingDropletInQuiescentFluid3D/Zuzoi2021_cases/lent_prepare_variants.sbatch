#!/bin/bash
##SBATCH -J interIsoFoam 
##SBATCH --mail-type=END
## Please use the full path :
##SBATCH -e /home/<username>/<project_name>/job_name.err.%j
##SBATCH -o /home/<username>/<project_name>/job_name.out.%j
##
#SBATCH --partition=test24 #partition for LB2 {test30m, test7d} are other options 
#SBATCH -n 1 # Number of cores
#SBATCH --mem-per-cpu=2000 # Main memory in MByte per MPI task
#SBATCH -t 00:35:00     # Hours, minutes and seconds, or '#SBATCH -t 10' = only minutes
#
## -------------------------------
## Afterwards you write your "payload" commands for the calculation, e.g.
lent_prepare_study_variants.py -m block Bussmann2002.parameter
