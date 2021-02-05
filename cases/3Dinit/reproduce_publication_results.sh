#!/usr/bin/env bash

# This script reproduces the studies underlying the "sphere and ellipsoid"
# result section in the publication
#   ADD REFERENCE

# The argument "0.25" is the parameter controlling the volume mesh perturbation
# and matches the one used in the publication
PERTURBATION=0.25

# Prepare study cases
./create_smci_smca_verification_study.sh "$PERTURBATION"
./create_refinement_study_smca.sh "$PERTURBATION"

# Ensure test are completely setup before they are run
wait

# Run studies
./run_smci_smca_verification_study.sh smciVofInit "$PERTURBATION"
./run_smci_smca_verification_study.sh smcaVofInit "$PERTURBATION"
./run_refinement_study_smca.sh "$PERTURBATION"

# Ensure studies have finished before post-processing
wait

# Generate tables and figures using jupyter notebook
# TODO: set GEOM_VOF_INIT here if it is not already set? (TT)
jupyter-nbconvert --execute --to=pdf smci-vof-init.ipynb
jupyter-nbconvert --execute --to=pdf smca-vof-init.ipynb
