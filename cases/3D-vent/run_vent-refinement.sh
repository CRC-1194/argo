#!/usr/bin/env bash

# PyFoam does requires a constant directory even if it is empty
mkdir constant

cd ..

pyFoamRunParameterVariation.py \
--every-variant-one-case-execution \
--create-database \
--no-mesh-create \
--no-server-process \
--no-execute-solver \
--no-case-setup \
3D-vent \
3D-vent/vent-refinement.parameter

wait

for STUDYDIR in $(find . -type d -name "vent-refinement.parameter*");
do
    cd $STUDYDIR
    cp ../3D-vent/box.fms .
    cp ../3D-vent/inside_vent.stl .
    cartesianMesh
    surfaceInitVolumeFraction -algorithm SMCA -checkVolume -surfaceFile inside_vent.stl
    cd ..
    mv $STUDYDIR 3Dinit
done
