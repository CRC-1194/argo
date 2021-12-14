#!/usr/bin/env bash

# Create the study split by the solver to be used:
argo-create-parameter-study.py study.parameter -p validation-interFoam -v 0,1,2 
argo-create-parameter-study.py study.parameter -p validation-interIsoFoam -v 3,4,5

# Initialize variants: create mesh and initialize fields
argo-initialize-parameter-study.py validation- -m blockMesh -f initFields.sh

# Wait a bit to ensure that initialization jobs have finished.
# Only relevant for runs on a cluster
echo "Wait 30 seconds for init jobs to finish."
sleep 30
wait

argo-run-study.py interFoam -d validation-interFoam -nc 2
argo-run-study.py interIsoFoam -d validation-interIsoFoam -nc 2
