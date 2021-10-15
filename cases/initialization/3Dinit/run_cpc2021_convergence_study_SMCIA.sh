#!/usr/bin/env bash
argo-create-parameter-study.py triSurface-SMCIA-convergence.parameter -p CPC2021
argo-initilize-parameter-study.py CPC2021-triSurface-SMCIA-convergence_00 -m blockMesh -f prepareCase.py

# Run cases using SMCA algorithm and agglomerate data 
argo-run-study.py surfaceInitVolumeFraction --solver-args="-algorithm SMCA" \
    -d CPC2021-triSurface-SMCIA-convergence_00 -nc 2
argo-agglomerate-study-data.py CPC2021-triSurface-SMCIA-convergence_00000_templateCase/vof-init-results-SMCA.csv -p triSurface-SMCIA-convergence.parameter -f CPC2021-convergence-SMCA

# Run cases using SMCI algorithm and agglomerate data 
argo-run-study.py surfaceInitVolumeFraction --solver-args="-algorithm SMCI" \
    -d CPC2021-triSurface-SMCIA-convergence_00 -nc 2
argo-agglomerate-study-data.py CPC2021-triSurface-SMCIA-convergence_00000_templateCase/vof-init-results-SMCI.csv -p triSurface-SMCIA-convergence.parameter -f CPC2021-convergence-SMCI
