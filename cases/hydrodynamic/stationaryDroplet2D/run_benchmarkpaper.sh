#!/usr/bin/env bash
argo-create-parameter-study.py benchmarkpaper.parameter -p interFoam $1
argo-create-parameter-study.py benchmarkpaper.parameter -p interIsoFoam $2

# Append option -j to submit initialization as SLURM jobs
argo-initialize-parameter-study.py interFoam-benchmarkpaper -m blockMesh -f initFields.sh -par -j
argo-initialize-parameter-study.py interIsoFoam-benchmarkpaper -m blockMesh -f initFields.sh -par -j

# Wait a bit to ensure that initialization jobs have finished.
echo "Wait 30 seconds for init jobs to finish."
sleep 30
wait

# Append option -j to submit simluations as SLURM jobs
argo-run-study.py interFoam -d interFoam-benchmarkpaper -np 4 -j
argo-run-study.py interIsoFoam -d interIsoFoam-benchmarkpaper -np 4 -j
