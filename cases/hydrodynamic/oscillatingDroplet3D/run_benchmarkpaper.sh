#!/usr/bin/env bash

# On a cluster: submit as jobs to the workload manager
# On a workstation: run directly
if [[ "$1" == "help" || "$1" == "-h" || -z $1 ]];
then
    echo "Please specify how to run the study:"
    echo "  slurm : submit as SLURM jobs."
    echo "  ws    : run direcly on workstation."
    echo "If you use 'ws', you can set the number of cases running in parallel"
    echo "with a second, integer argument."
    exit 0
fi

RUNNER=$1

if [[ "$RUNNER" != "slurm" && "$RUNNER" != "ws" ]];
then
    echo "Error: unknown runner. Use either 'slurm' or 'ws'."
    exit 1
fi

# Set job-flag a.k.a. SLURM flag
JOB_ARG=""
if [[ "$RUNNER" == "slurm" ]];
then
    $JOB_ARG = "-j"
fi

# Set number of cases run concurrently. Only use 1/2 as 
# interFoam and interIsoFoam are executed concurrently.
NUM_CASES=""
if [[ "$RUNNER" == "ws" && ! -z $2 ]];
then
    let HALF=$2/2
    ROUND=$( printf "%.0f" $HALF )
    NUM_CASES="--num-cases-parallel $ROUND"
fi

# Create the study split by the solver to be used:
# NOTE: change if the parameter file changes so that the variation numbers
# corresponds to interFoam / interIsoFoam
argo-create-parameter-study.py benchmarkpaper.parameter -p benchmark-interFoam -v 0-8
argo-create-parameter-study.py benchmarkpaper.parameter -p benchmark-interIsoFoam -v 9-17
argo-create-parameter-study.py benchmarkpaper.parameter -p benchmark-interFlow -v 18-35

# Initialize variants: create mesh and initialize fields
argo-initialize-parameter-study.py benchmark- -m blockMesh -f initFields.sh -par $JOB_ARG

# Wait a bit to ensure that initialization jobs have finished.
if [[ "$RUNNER" == "slurm" ]];
then
    echo "Wait 30 seconds for init jobs to finish."
    sleep 30
    wait
fi

nohup argo-run-study.py interFoam -d benchmark-interFoam -np 27 $JOB_ARG $NUM_CASES > logs.interFoam &
nohup argo-run-study.py interIsoFoam -d benchmark-interIsoFoam -np 27 $JOB_ARG $NUM_CASES > logs.interIsoFoam &
nohup argo-run-study.py interFlow -d benchmark-interFlow -np 27 $JOB_ARG $NUM_CASES > logs.interFlow &
