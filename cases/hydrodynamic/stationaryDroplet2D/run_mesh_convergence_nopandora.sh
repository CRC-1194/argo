#!/usr/bin/env bash
argo-create-parameter-study.py mesh_convergence_nopandora.parameter -p interFoam $1
argo-create-parameter-study.py mesh_convergence_nopandora.parameter -p interIsoFoam $2

argo-initialize-parameter-study.py interFoam-mesh -m blockMesh -f initFields.sh -par
argo-initialize-parameter-study.py interIsoFoam-mesh -m blockMesh -f initFields.sh -par

argo-run-study.py interFoam -d interFoam-mesh -np 4
argo-run-study.py interIsoFoam -d interIsoFoam-mesh -np 4
