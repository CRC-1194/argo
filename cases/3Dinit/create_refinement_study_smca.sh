#!/usr/bin/env bash

# Create a study verifying the effectiveness of mesh refinement used by the SMCA algorithm.
PERTURBATION=$1

# Equidistant
./create_parameter_study.py --study_name refinement-equidistant --surface sphere --parameter_file smca-refinement-convergence-blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.equidistant.sphere.refinement & 

# Perturbed 
./create_parameter_study.py --study_name refinement-perturbed_alpha_"$PERTURBATION" --perturb_mesh "$PERTURBATION"  --surface sphere --parameter_file smca-refinement-convergence-blockMesh.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.perturbed.sphere.refinement & 

wait
