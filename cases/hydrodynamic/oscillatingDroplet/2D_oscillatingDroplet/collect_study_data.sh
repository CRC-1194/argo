#!/usr/bin/env bash

# Need to collect data for interFoam and interIsoFoam separately due to different case base names.
argo-agglomerate-study-data.py validation-interFoam-study_00000_templateCase/postProcessing/interfaceHeight1/0/height.dat -p study.parameter -f oscillating_droplet_interFoam
argo-agglomerate-study-data.py validation-interIsoFoam-study_00009_templateCase/postProcessing/interfaceHeight1/0/height.dat -p study.parameter -f oscillating_droplet_interIsoFoam
