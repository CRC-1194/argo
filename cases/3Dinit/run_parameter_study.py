#!/bin/env python
import argparse
import sys
import os
from subprocess import call

parser = argparse.ArgumentParser(description='Run an application in each case directory that fits the pattern.')

parser.add_argument('application', type=str,
                    choices=['surfaceCellVofInit', 'poFoamTestVofInit'],
                    help='Name of the application to be run, e.g. surfaceCellVofInit.')

parser.add_argument('--dir_pattern', dest="dir_pattern", type=str, required=True,
                    help='Pattern contained in the name of each initialization directory.')

parser.add_argument('--surface_file', dest="surface_file", type=str, required=True, 
                    help='Surface mesh file used for the volume fraction initialization.')

parser.add_argument('--slurm', action='store_true',
                    help='Submit each variation as a SLURM job.')

args = parser.parse_args()

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and args.dir_pattern in parameter_dir]
    parameter_dirs.sort()

    application_call = ''

    if args.application == "surfaceCellVofInit":
        application_call = [args.application, 
                            "-checkVolume", 
                            "-writeGeometry",
                            "-surfaceFile", 
                            args.surface_file] 
    elif args.application == "poFoamTestVofInit":
        application_call = [args.application, 
                            "-writeFields",
                            "-surfaceFile", 
                            args.surface_file] 

    for parameter_dir in parameter_dirs: 
        pwd = os.getcwd()
        os.chdir(parameter_dir)
        if (args.slurm):
            call(["sbatch", os.path.join(os.pardir, args.application+".sbatch")])
        else:
            call(application_call)
        os.chdir(pwd)
