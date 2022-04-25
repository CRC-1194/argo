#!/usr/bin/env bash

# Create the study split by the solver to be used:
argo-create-parameter-study.py benchmarkpaper.parameter -p benchmark-interFoam -v 0,1,2 
argo-create-parameter-study.py benchmarkpaper.parameter -p benchmark-interIsoFoam -v 3,4,5

# Initialize variants: create mesh and initialize fields
argo-initialize-parameter-study.py benchmark- -m blockMesh -f initFields.sh

# Wait a bit to ensure that initialization jobs have finished.
# Only relevant for runs on a cluster
#echo "Wait 30 seconds for init jobs to finish."
#sleep 30
#wait

nohup argo-run-study.py interFoam -d benchmark-interFoam -nc 2 > logs.interFoam &
nohup argo-run-study.py interIsoFoam -d benchmark-interIsoFoam -nc 2 > logs.interIsoFoam &
