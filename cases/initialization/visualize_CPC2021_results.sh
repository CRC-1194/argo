#!/usr/bin/env bash

# This reproduces all data and plots used in publication
# https://arxiv.org/abs/2101.08511 

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
