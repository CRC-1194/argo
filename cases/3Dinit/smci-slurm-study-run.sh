#!/usr/bin/env bash

./run_parameter_study.py --slurm \
    --dir_pattern hexmesh_sphere_blockMesh.parameter_000 \
    --surface_file sphere.vtk && \
./run_parameter_study.py --slurm \
    --dir_pattern hexmesh-perturbed-alpha025_sphere_blockMesh.parameter_000 \
    --surface_file sphere.vtk && \
./run_parameter_study.py --slurm \
    --dir_pattern hexmesh_ellipsoid_blockMesh.parameter_000 \
    --surface_file ellipsoid.vtk && \ 
./run_parameter_study.py --slurm \
    --dir_pattern hexmesh-perturbed-alpha025_ellipsoid_blockMesh.parameter_000 \
    --surface_file ellipsoid.vtk

