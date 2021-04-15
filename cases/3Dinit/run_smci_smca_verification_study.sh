#!/usr/bin/env bash

ALGORITHM=$1
PERTURBATION=$2

# Sphere

## Equidistant 
./run_parameter_study.py "$ALGORITHM" --dir_pattern equidistant_sphere --surface_file sphere.vtk 2>&1 > log.run.equidistant.sphere."$ALGORITHM" & 

## Perturbed 
./run_parameter_study.py "$ALGORITHM" --dir_pattern perturbed_alpha_"$PERTURBATION"_sphere --surface_file sphere.vtk 2>&1 > log.run.perturbed.sphere."$ALGORITHM" & 


# Ellipsoid 

## Equidistant 
./run_parameter_study.py "$ALGORITHM" --dir_pattern equidistant_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.equidistant.ellipsoid."$ALGORITHM" & 

## Perturbed 
./run_parameter_study.py "$ALGORITHM" --dir_pattern perturbed_alpha_"$PERTURBATION"_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.perturbed.ellipsoid."$ALGORITHM" & 

wait
