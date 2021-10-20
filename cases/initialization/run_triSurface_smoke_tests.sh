#!/usr/bin/env bash

# This script runs a subset of the CPC2021 cases (mainly low resolution variants)
# as basis for regression tests

# Source auxiliary scripts for creating and running parameter studies
source ../scripts/bashrc

# Run parameter studies
cd 3Dinit
./run_cpc2021_convergence_study_SMCIA.sh --variants="36,52"
./run_cpc2021_refinement_study_SMCA.sh --variants="1"
cd ../3D-vent
./run_cpc2021_vent-refinement.sh --variants="1"
cd ..

mkdir -p triSurface-smoke-test-results
rm -rf triSurface-smoke-test-results/*

mv 3Dinit/*.{csv,json} triSurface-smoke-test-results/
mv *.{csv,json} triSurface-smoke-test-results/
