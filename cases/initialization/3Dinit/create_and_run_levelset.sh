#!/usr/bin/env bash

# Equidistant
./create_parameter_study.py --study_name vofImplicit --surface ellipsoid --levelset --parameter_file blockMeshSmokeTestLevelset.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.levelset.ellipsoid

wait

echo "Running study..."
./run_parameter_study.py SMCA --levelset --dir_pattern vofImplicit --surface_file irrelevant 2>&1 > log.run.levelset.ellipsoid
