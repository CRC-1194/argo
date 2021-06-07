#!/usr/bin/env bash

# Equidistant, refinement level driven studies
./create_parameter_study.py --study_name levelBasedLS --surface sphere --levelset --parameter_file blockMeshSmokeTestLevelset.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.levelBasedLS.sphere
./create_parameter_study.py --study_name levelBasedLS --surface ellipsoid --levelset --parameter_file blockMeshSmokeTestLevelset.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.levelBasedLS.ellipsoid

# Equidistant, accuracy driven studies
./create_parameter_study.py --study_name accuracyDrivenLS --surface sphere --levelset --parameter_file blockMeshAccuracyDrivenLevelset.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.accuracyDrivenLS.sphere
./create_parameter_study.py --study_name accuracyDrivenLS --surface ellipsoid --levelset --parameter_file blockMeshAccuracyDrivenLevelset.parameter --mesh_generator blockMesh --template_case templateCase 2>&1 > log.create.accuracyDrivenLS.ellipsoid

wait

echo "Running study..."
./run_parameter_study.py SMCA --levelset --dir_pattern levelBasedLS_sphere --surface_file irrelevant 2>&1 > log.run.levelBasedLS.sphere
./run_parameter_study.py SMCA --levelset --dir_pattern levelBasedLS_ellipsoid --surface_file irrelevant 2>&1 > log.run.levelBasedLS.ellipsoid

./run_parameter_study.py SMCA --levelset --dir_pattern accuracyDrivenLS_sphere --surface_file irrelevant 2>&1 > log.run.accuracyDrivenLS.sphere
./run_parameter_study.py SMCA --levelset --dir_pattern accuracyDrivenLS_ellipsoid --surface_file irrelevant 2>&1 > log.run.accuracyDrivenLS.ellipsoid
