#!/usr/bin/env bash

PERTURBATION=$1

# Sphere
## Equidistant 
./run_parameter_study.py SMCA --dir_pattern refinement-equidistant --surface_file sphere.vtk 2>&1 > log.run.equidistant.sphere.refinement & 

## Perturbed 
./run_parameter_study.py SMCA --dir_pattern refinement-perturbed_alpha_"$PERTURBATION" --surface_file sphere.vtk 2>&1 > log.run.perturbed.sphere.refinement & 

wait
