#!/usr/bin/bash

gmsh -3 -algo del3d j.compfluid.2018.10.021-table3-box.geo -o box.msh
gmshToFoam box.msh
gmsh -2 j.compfluid.2018.10.021-table3-sphere.geo -o sphere.vtk  
cp system/vofInitDict-j.compfluid.2018.10.021-triSurface system/vofInitDict
