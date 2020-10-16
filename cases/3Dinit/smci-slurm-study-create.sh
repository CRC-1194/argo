#!/usr/bin/env bash

./create_parameter_study.py --slurm \
    --study_name hexmesh \
    --parameter_file blockMesh.parameter \
    --surface sphere \
    --mesh_generator blockMesh #&&
#./create_parameter_study.py --slurm \
    #--study_name hexmesh \
    #--parameter_file blockMesh.parameter \
    #--surface ellipsoid \
    #--mesh_generator blockMesh &&
#./create_parameter_study.py --slurm \
    #--study_name hexmesh-perturbed-alpha025 \
    #--parameter_file blockMesh.parameter \
    #--surface sphere \
    #--mesh_generator blockMesh &&
#./create_parameter_study.py --slurm \
    #--study_name hexmesh-perturbed-alpha025 \
    #--parameter_file blockMesh.parameter \
    #--surface ellipsoid \
    #--mesh_generator blockMesh 
