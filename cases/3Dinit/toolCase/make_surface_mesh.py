#!/bin/env python

import argparse
import sys
import os
from subprocess import call

parser = argparse.ArgumentParser(description='Generates surface meshes from .geo files using gmsh and stores them as .vtk and .stl files.')


parser.add_argument('--geo_file', dest="geo_file", type=str, required=True,
                    help='GMSH .geo input file for the surface mesh.')


if __name__ == '__main__':

    args = parser.parse_args(sys.argv[1:])
    base_name = "".join(args.geo_file.split('.')[:-1])
    
    print(base_name)

    # Generate STL and VTK surface meshes.
    call(["gmsh","-2", args.geo_file, "-o", "%s.vtk" % base_name])
    call(["gmsh","-2", args.geo_file, "-o", "%s.stl" % base_name])
    call(["sed","-i", r"s/Created by Gmsh/%s/" % base_name, "%s.stl" % base_name])
