#!/usr/bin/env bash

# Create the study split by the solver to be used:
argo-create-parameter-study.py study.parameter -p validation-interFoam -v 0-8 
argo-create-parameter-study.py study.parameter -p validation-interIsoFoam -v 9-17 

# Initialize variants: create mesh and initialize fields
argo-initilize-parameter-study.py validation- -m none -f initFields.sh -j -par

# Wait a bit to ensure that initialization jobs have finished.
echo "Wait 30 seconds for init jobs to finish."
sleep 30
wait

argo-run-study.py interFoam -d validation-interFoam -j -np 4
argo-run-study.py interIsoFoam -d validation-interIsoFoam -j -np 4
