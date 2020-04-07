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

parser.add_argument('--surface', dest="surface", type=str, required=True,
                    help='Name of the surface used for initialization.')

parser.add_argument('--mesh_generator', dest="mesh_generator", type=str, default="", required=True,
                    help='Supported mesh generators: blockMesh, cartesianMesh, tetMesh, polyMesh')

parser.add_argument('--serial', dest="serial", type=bool, default=True,
                    help='Generate the mesh.')

parser.add_argument('--slurm', dest="slurm", type=bool, default=False,
                    help='Use SLURM workload manager to submit mesh generation jobs.')

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
        if (args.serial):
            call([args.mesh_generator])
        elif (args.slurm): 
            call(["sbatch", os.path.join(pwd, "%s.sbatch" % args.mesh_generator)])
        geo_file = "%s.geo" % args.surface
        print("Using GMSH file %s for surface generation." % geo_file)
        if (os.path.exists(geo_file) and os.path.isfile(geo_file)):
            call(["gmsh", "-2", geo_file, "-o", "%s.vtk" % args.surface])
        os.chdir(pwd)
