#!/usr/bin/env bash

PERTURBATION=$1

# Sphere

## Equidistant 
./create_parameter_study.py --study_name equidistant --surface sphere --parameter_file blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.equidistant.sphere & 

#
### Perturbed 
#./create_parameter_study.py --study_name perturbed_alpha_"$PERTURBATION" --perturb_mesh "$PERTURBATION"  --surface sphere --parameter_file blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.perturbed.sphere & 
#
## Ellipsoid 
#
### Equidistant 
#./create_parameter_study.py --study_name equidistant --surface ellipsoid --parameter_file blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.equidistant.ellipsoid & 
#
### Perturbed 
#./create_parameter_study.py --study_name perturbed_alpha_"$PERTURBATION" --perturb_mesh "$PERTURBATION"  --surface ellipsoid --parameter_file blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.perturbed.ellipsoid & 
