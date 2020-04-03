#!/bin/env python
import argparse
import sys
import os
from subprocess import call

parser = argparse.ArgumentParser(description='Runs the surfaceCellVofInit in each case directory that fits the pattern.')

parser.add_argument('--dir_pattern', dest="dir_pattern", type=str, required=True,
                    help='Pattern contained in the name of each initialization directory.')

parser.add_argument('--surface_file', dest="surface_file", type=str, required=True,
                    help='Surface mesh file used for the volume fraction initialization.')

parser.add_argument('--serial', dest="serial", type=bool, default=True,
                    help='Generate the mesh.')

parser.add_argument('--slurm', dest="slurm", type=bool, default=False,
                    help='Use SLURM workload manager to submit mesh generation jobs.')

args = parser.parse_args()

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and args.dir_pattern in parameter_dir]


    rel_surface_file = os.path.join(os.pardir, args.surface_file)
    print(rel_surface_file)

    vof_init_call = ["surfaceCellVofInit", 
                     "-checkVolume", 
                     "-surfaceFile", 
                     rel_surface_file] 

    if (args.serial):
        args.slurm = False
        for parameter_dir in parameter_dirs: 
            pwd = os.getcwd()
            os.chdir(parameter_dir)
            call(vof_init_call)
            os.chdir(pwd)
    elif (args.slurm):
        for parameter_dir in parameter_dirs: 
            pwd = os.getcwd()
            os.chdir(parameter_dir)
            all(["sbatch", os.path.join(os.pardir, "surfaceCellVofInit.sbatch")])
            os.chdir(pwd)
