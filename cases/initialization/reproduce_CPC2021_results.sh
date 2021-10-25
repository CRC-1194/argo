#!/usr/bin/env bash

# This reproduces all data and plots used in publication
# https://arxiv.org/abs/2101.08511 

# Source auxiliary scripts for creating and running parameter studies
source ../scripts/bashrc

# Run parameter studies
cd 3Dinit
./run_cpc2021_convergence_study_SMCIA.sh
./run_cpc2021_refinement_study_SMCA.sh
./run_cpc2021_levelSet_comparison_j.compfluid.2018.10.021.sh
./run_cpc2021_triSurface_comparison_j.compfluid.2018.10.021.sh
cd ../3D-vent
./run_cpc2021_vent-refinement.sh
cd ..

# Move data to a single folder
if [[ -z "${SMCIA_VOF_INIT_RESULTS}" ]]; then
    mkdir -p CPC2021-results
    export SMCIA_VOF_INIT_RESULTS=.
fi
rm -rf CPC2021-results/*

mv 3Dinit/*.{csv,json} CPC2021-results/
mv *.{csv,json} CPC2021-results/
cp 3Dinit/plot_study.py CPC2021-results/
cp 3Dinit/triSurface-SMCIA-convergence.ipynb CPC2021-results/triSurface-SMCA-convergence.ipynb
cp 3Dinit/triSurface-SMCIA-convergence.ipynb CPC2021-results/triSurface-SMCI-convergence.ipynb
cp 3Dinit/triSurface-SMCA-refinement-convergence.ipynb CPC2021-results/
cp 3Dinit/j.compfluid.2018.10.021-table3.ipynb CPC2021-results/
cp 3D-vent/vent-refinement.ipynb CPC2021-results/

# Execute notebooks to generate plots
cd CPC2021-results

export VOF_INIT_ALGORITHM=SMCA
jupyter-nbconvert --execute --inplace triSurface-SMCA-convergence.ipynb
jupyter-nbconvert --execute --to=html triSurface-SMCA-convergence.ipynb

jupyter-nbconvert --execute --inplace j.compfluid.2018.10.021-table3.ipynb 
jupyter-nbconvert --execute --to=html j.compfluid.2018.10.021-table3.ipynb 

export VOF_INIT_ALGORITHM=SMCI
jupyter-nbconvert --execute --inplace triSurface-SMCI-convergence.ipynb
jupyter-nbconvert --execute --to=html triSurface-SMCI-convergence.ipynb

jupyter-nbconvert --execute --inplace triSurface-SMCA-refinement-convergence.ipynb
jupyter-nbconvert --execute --to=html triSurface-SMCA-refinement-convergence.ipynb

jupyter-nbconvert --execute --inplace vent-refinement.ipynb
jupyter-nbconvert --execute --to=html vent-refinement.ipynb
