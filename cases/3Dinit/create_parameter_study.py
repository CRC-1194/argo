#!/bin/env python
import argparse
import sys
import os
from subprocess import Popen 
from subprocess import call

parser = argparse.ArgumentParser(description='Generates simulation cases for parameter study using PyFoam.')

parser.add_argument('--parameter_file', dest="parameter_file", type=str, required=True,
                    help='PyFoam .parameter file')

parser.add_argument('--template_case', dest="template_case", type=str, 
                    help='OpenFOAM template case', default="templateCase")

parser.add_argument('--surface', dest="surface", type=str, required=True,
                    help='Name of the surface used for initialization: available surfaces see are .geo files in the templateCase directory.')

parser.add_argument('--mesh_generator', dest="mesh_generator", type=str, default="", required=True,
                    help='Supported mesh generators: blockMesh, cartesianMesh, tetMesh, polyMesh')

parser.add_argument('--serial', dest="serial_run", action='store_true',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.add_argument('--slurm', dest="serial_run", action='store_false',
                    help='Use SLURM workload manager to submit mesh generation jobs.')

parser.set_defaults(serial_run=True)

args = parser.parse_args()

if __name__ == '__main__':
    args = parser.parse_args(sys.argv[1:])

    prefix = args.surface + "_" + args.parameter_file

    call_args = ["pyFoamRunParameterVariation.py", 
                 "--every-variant-one-case-execution", 
                 "--create-database", 
                 "--no-mesh-create",
                 "--no-server-process",
                 "--no-execute-solver",
                 "--no-case-setup",
                 "--cloned-case-prefix=%s" % prefix,
                 args.template_case, 
                 args.parameter_file]

    call(call_args)

    parameter_dirs = [parameter_dir for parameter_dir in os.listdir(os.curdir) \
                      if os.path.isdir(parameter_dir) \
                      and "00" in parameter_dir \
                      and parameter_dir.startswith(prefix) \
                      and parameter_dir.endswith(args.template_case)]
    parameter_dirs.sort()

    for parameter_dir in parameter_dirs: 
        pwd = os.getcwd()
        os.chdir(parameter_dir)
        geo_file = "%s.geo" % args.surface

        if (args.serial_run):
            call([args.mesh_generator])
            print("Using GMSH file %s for surface generation." % geo_file)
            if (os.path.exists(geo_file) and os.path.isfile(geo_file)):
                call(["gmsh", "-2", geo_file, "-o", "%s.vtk" % args.surface])
        else: 
            base_command="srun --mem-per-cpu=500 --time=00:10 --ntasks=1"

            # Submit volume mesh generation to the SLURM workload manager.
            variable_command=" --job-name %s -o %s.log %s >/dev/null 2>&1 &" % \
                    (args.mesh_generator, args.mesh_generator, args.mesh_generator)
            call(base_command + variable_command, shell=True)

            # Submit surface mesh generation to the SLURM workload manager.
            variable_command=" --job-name %s -o %s.log gmsh -2 %s.geo -o %s.vtk >/dev/null 2>&1 &" % \
                    (args.surface, args.surface, args.surface, args.surface)
            call(base_command + variable_command, shell=True)

        os.chdir(pwd)
