#!/usr/bin/env bash

# This script reproduces the studies underlying the "sphere and ellipsoid"
# result section in the publication
#   ADD REFERENCE
# Computed data (figures and tables) are written to ./results

# The argument "0.25" is the parameter controlling the volume mesh perturbation
# and matches the one used in the publication
PERTURBATION=0.25

# Prepare study cases
./create_smci_smca_verification_study.sh "$PERTURBATION"
./create_refinement_study_smca.sh "$PERTURBATION"

# Ensure test are completely setup before they are run
wait

# Run studies for sphere and ellipsoid
./run_smci_smca_verification_study.sh SMCI "$PERTURBATION"
./run_smci_smca_verification_study.sh SMCA "$PERTURBATION"
./run_refinement_study_smca.sh "$PERTURBATION"

wait

# Run refinemnt study for CAD geometry
cd ../3D-vent
./run_vent-refinement.sh
cd ../3Dinit

# Ensure studies have finished before post-processing
wait

# Generate tables and figures using jupyter notebook
if [[ -z "${GEOM_VOF_INIT}" ]]; then
    mkdir -p results
    export GEOM_VOF_INIT=./results
fi

jupyter-nbconvert --execute --to=pdf smci-vof-init.ipynb
jupyter-nbconvert --execute --to=pdf smca-vof-init.ipynb
