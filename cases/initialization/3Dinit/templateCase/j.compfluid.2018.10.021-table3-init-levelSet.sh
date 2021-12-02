#!/usr/bin/bash

gmsh -3 -algo del3d j.compfluid.2018.10.021-table3-box.geo -o box.msh
gmshToFoam box.msh
cp system/vofInitDict-j.compfluid.2018.10.021-levelSet system/vofInitDict
