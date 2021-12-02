#!/usr/bin/env bash
argo-create-parameter-study.py triSurface-SMCA-refinement-convergence.parameter -p CPC2021 $1
argo-initilize-parameter-study.py CPC2021-triSurface-SMCA-refinement-convergence_00 -m blockMesh -f prepareCase.py
argo-run-study.py surfaceInitVolumeFraction --solver-args="-algorithm SMCA" -d CPC2021-triSurface-SMCA-refinement-convergence_00 -nc 2
argo-agglomerate-study-data.py CPC2021-triSurface-SMCA-refinement-convergence_00000_templateCase/vof-init-results-SMCA.csv -p triSurface-SMCA-refinement-convergence.parameter -f CPC2021-refinement-convergence-SMCA
