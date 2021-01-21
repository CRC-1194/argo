#!/usr/bin/env bash

PERTURBATION=$1

# Sphere

## Equidistant 
./run_parameter_study.py smciVofInit --dir_pattern equidistant_sphere --surface_file sphere.vtk 2>&1 > log.run.equidistant.sphere & 

#
### Perturbed 
#./run_parameter_study.py smciVofInit --dir_pattern perturbed_alpha_"$PERTURBATION"_sphere --surface_file sphere.vtk 2>&1 > log.run.perturbed.sphere & 
#
## Ellipsoid 
#
### Equidistant 
#./run_parameter_study.py smciVofInit --dir_pattern equidistant_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.equidistant.ellipsoid & 
#
### Perturbed 
#./run_parameter_study.py smciVofInit --dir_pattern perturbed_alpha_"$PERTURBATION"_ellipsoid --surface_file ellipsoid.vtk 2>&1 > log.run.perturbed.ellipsoid & 
