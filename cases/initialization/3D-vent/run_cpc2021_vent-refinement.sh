#!/usr/bin/env bash

# PyFoam does requires a constant directory even if it is empty
mkdir -p constant
cp vent-refinement.parameter ..
cd ..
# Argo scripts require a default.parameter file which can be empty
touch default.parameter

argo-create-parameter-study.py vent-refinement.parameter -p CPC2021 -t 3D-vent $1
argo-initialize-parameter-study.py CPC2021-vent-refinement_0000 -m none -f prepareCase.sh
argo-run-study.py surfaceInitVolumeFraction --solver-args="-algorithm SMCA" -d CPC2021-vent-refinement_0000
argo-agglomerate-study-data.py CPC2021-vent-refinement_00000_3D-vent/vof-init-results-SMCA.csv \
    -p vent-refinement.parameter -f CPC2021-vent-refinement-SMCA

# Cleanup temporary files/copies
rm -f vent-refinement.parameter vent-refinement.parameter.database default.parameter

cd 3D-vent
