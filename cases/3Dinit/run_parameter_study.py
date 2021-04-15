#!/usr/bin/env python
import argparse
import sys
import os
from subprocess import call
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

parser = argparse.ArgumentParser(description='Run an application in each case directory that fits the pattern.')

parser.add_argument('algorithm', type=str,
                    choices=['SMCI', 'SMCA'],
                    help='Choose the VoF calculation algorithm.')

parser.add_argument('--dir_pattern', dest="dir_pattern", type=str, required=True,
                    help='Pattern contained in the name of each initialization directory.')

parser.add_argument('--surface_file', dest="surface_file", type=str, required=True, 
                    help='Surface mesh file used for the volume fraction initialization.')

parser.add_argument('--slurm', dest="slurm_run", action='store_true',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.set_defaults(serial_run=True)

args = parser.parse_args()

if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and args.dir_pattern in parameter_dir]
    parameter_dirs.sort()

    vof_init_call = ["surfaceInitVolumeFraction",
                     "-algorithm",
                     args.algorithm,
                     "-checkVolume", 
                     "-surfaceFile", 
                     args.surface_file] 

    for parameter_dir in parameter_dirs: 
        pwd = os.getcwd()
        os.chdir(parameter_dir)

        print(parameter_dir)

        # Get mesh density from blockMeshDict
        blockMeshDict = ParsedParameterFile("./system/blockMeshDict") 
        Nc = blockMeshDict["blocks"][2][0]

        if (args.slurm_run):

            if (Nc >= 256): # Use more memory for the largest resolution. 
                base_command="srun -A project01456 --mem-per-cpu=10000 --time=01:00:00 --ntasks=1"
            else:
                base_command="srun -A project01456 --mem-per-cpu=5000 --time=01:00:00 --ntasks=1"

            variable_command=" --cpu-freq=HighM1-HighM1 --job-name smciVofInit -o slurm-%s.log smciVofInit -checkVolume -surfaceFile %s >/dev/null 2>&1 &" % (args.surface_file, args.surface_file)

            print(base_command + variable_command)

            call(base_command + variable_command, shell=True)

        else: 
            call(vof_init_call)

        os.chdir(pwd)
