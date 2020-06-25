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

parser.add_argument('--serial', dest="serial_run", action='store_true',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.add_argument('--slurm', dest="serial_run", action='store_false',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.set_defaults(serial_run=True)

args = parser.parse_args()

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and args.dir_pattern in parameter_dir]
    parameter_dirs.sort()

    vof_init_call = ["surfaceCellVofInit", 
                     "-checkVolume", 
                     "-surfaceFile", 
                     args.surface_file] 

    for parameter_dir in parameter_dirs: 
        pwd = os.getcwd()
        os.chdir(parameter_dir)
        if (args.serial_run):
            call(vof_init_call)
        else: 
            # Memory and run time can now be made dependent on mesh size. 
            base_command="srun --mem-per-cpu=1000 --time=05:00 --ntasks=1"

            variable_command=" --job-name surfaceCellVofInit " + \
                " -o %s.log " + \
                " surfaceCellVofInit -checkVolume -surfaceFile " + \
                " %s >/dev/null 2>&1 &" % (args.surface_file, args.surface_file)
            call(base_command + variable_command, shell=True)
        os.chdir(pwd)
