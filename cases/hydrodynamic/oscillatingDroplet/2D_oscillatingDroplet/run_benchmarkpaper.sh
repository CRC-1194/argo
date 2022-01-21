#!/usr/bin/env bash

# Create the study split by the solver to be used:
argo-create-parameter-study.py benchmarkpaper.parameter -p validation-interFoam $1
argo-create-parameter-study.py benchmarkpaper.parameter -p validation-interIsoFoam $2

# Initialize variants: create mesh and initialize fields
argo-initialize-parameter-study.py validation- -m blockMesh -f initFields.sh -par #-j

# Wait a bit to ensure that initialization jobs have finished.
echo "Wait 30 seconds for init jobs to finish."
sleep 30
wait

argo-run-study.py interFoam -d validation-interFoam -np 4 #-j
argo-run-study.py interIsoFoam -d validation-interIsoFoam -np 4 #-j
