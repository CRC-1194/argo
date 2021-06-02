#!/usr/bin/env python3
import argparse
import sys
import os
from subprocess import call
#from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

parser = argparse.ArgumentParser(description='Run an application in each case directory that fits the pattern.')

parser.add_argument('--dir_pattern', dest="dir_pattern", type=str, required=True,
                    help='Pattern contained in the name of each initialization directory.')

parser.add_argument('--slurm', dest="slurm_run", action='store_true',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.add_argument('--solver', dest="solver", type=str, required=True,
                    help='Set the solver for the simulation, e.g. interIsoPandoraFoam.')

parser.set_defaults(serial_run=True)

args = parser.parse_args()

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and args.dir_pattern in parameter_dir]
    parameter_dirs.sort()

    for parameter_dir in parameter_dirs: 
        pwd = os.getcwd()
        os.chdir(parameter_dir)

        print(parameter_dir)

        if (args.slurm_run):
            call("sbatch", "../" + args.solver + ".sbatch")
        else: 
            call(args.solver)

        os.chdir(pwd)
