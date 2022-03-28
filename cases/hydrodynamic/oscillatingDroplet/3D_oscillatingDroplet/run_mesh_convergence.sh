#!/usr/bin/env bash
argo-create-parameter-study.py mesh_convergence.parameter -p interFoam $1
argo-create-parameter-study.py mesh_convergence.parameter -p interIsoFoam $2

argo-initilize-parameter-study.py interFoam-mesh -m blockMesh -f initFields.sh -par
argo-initilize-parameter-study.py interIsoFoam-mesh -m blockMesh -f initFields.sh -par

argo-run-study.py interFoam -d interFoam-mesh -np 8
argo-run-study.py interIsoFoam -d interIsoFoam-mesh -np 8
