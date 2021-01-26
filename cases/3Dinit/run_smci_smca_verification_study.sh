#!/usr/bin/env bash

ALG_COMMAND=$1
PERTURBATION=$2

# Sphere

## Equidistant 
./run_parameter_study.py "$ALG_COMMAND" --dir_pattern equidistant_sphere --surface_file sphere.vtk 2>&1 > log.run.equidistant.sphere."$ALG_COMMAND" & 

## Perturbed 
./run_parameter_study.py "$ALG_COMMAND" --dir_pattern perturbed_alpha_"$PERTURBATION"_sphere --surface_file sphere.vtk 2>&1 > log.run.perturbed.sphere."$ALG_COMMAND" & 


# Ellipsoid 

## Equidistant 
./run_parameter_study.py "$ALG_COMMAND" --dir_pattern equidistant_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.equidistant.ellipsoid."$ALG_COMMAND" & 

## Perturbed 
./run_parameter_study.py "$ALG_COMMAND" --dir_pattern perturbed_alpha_"$PERTURBATION"_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.perturbed.ellipsoid."$ALG_COMMAND" & 

wait
