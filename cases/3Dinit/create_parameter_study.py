#!/bin/env python
import argparse
import sys
import os
from subprocess import call

parser = argparse.ArgumentParser(description='Generates simulation cases for parameter study using PyFoam.')

parser.add_argument('--parameter_file', dest="parameter_file", type=str, required=True,
                    help='PyFoam .parameter file')

parser.add_argument('--template_case', dest="template_case", type=str, 
                    help='OpenFOAM template case', default="templateCase")

parser.add_argument('--mesh_generator', dest="mesh_generator", type=str, default="",
                    help='Supported mesh generators: blockMesh, cartesianMesh, tetMesh, polyMesh')

parser.add_argument('--prefix', dest="prefix", type=str, default="",
                    help='String for categorizing the parameter study.')

parser.add_argument('--serial', dest="serial", type=bool, default=True,
                    help='Generate the mesh.')

parser.add_argument('--slurm', dest="slurm", type=bool, default=False,
                    help='Use SLURM workload manager to submit mesh generation jobs.')

args = parser.parse_args()

if __name__ == '__main__':
    args = parser.parse_args(sys.argv[1:])

    call_args = ["pyFoamRunParameterVariation.py", 
                 "--every-variant-one-case-execution", 
                 "--create-database", 
                 "--no-mesh-create",
                 "--no-server-process",
                 "--no-execute-solver",
                 "--no-case-setup",
                 "--cloned-case-prefix=%s" % args.prefix,
                 args.template_case, 
                 args.parameter_file]

    call(call_args)

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and "00" in parameter_dir \
                      and parameter_dir.startswith(args.prefix + args.parameter_file) \
                      and parameter_dir.endswith(args.template_case)]
    parameter_dirs.sort()

    print(parameter_dirs)

    if (args.serial):
        args.slurm = False
        for parameter_dir in parameter_dirs: 
            call([args.mesh_generator, "-case", os.path.join(os.curdir,parameter_dir)])
    elif (args.slurm):
        for parameter_dir in parameter_dirs: 
            pwd = os.path.cwd()
            os.chdir(parameter_dir)
            all(["sbatch", os.path.join(pwd, "%s.sbatch" % args.mesh_generator)])
            os.chdir(pwd)
