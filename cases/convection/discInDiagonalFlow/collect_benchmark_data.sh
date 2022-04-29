#!/usr/bin/env bash

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py benchmark-interFoam-benchmarkpaper_00000_templateCase/advectionErrors.csv -p benchmarkpaper.parameter -f disc_in_diagonal_flow_interFoam
argo-agglomerate-study-data.py benchmark-interIsoFoam-benchmarkpaper_00003_templateCase/advectionErrors.csv -p benchmarkpaper.parameter -f disc_in_diagonal_flow_interIsoFoam