# Overview

This is a parameter study for intersecting a pseudo-2D cylinder mesh with a pseudo-2D volume mesh in OpenFOAM.

You can run the study serially, or on a cluster that is using the SLURM workload manager.

Measured parameter  | Description
--- | ---
Nt | number of cells or triangles in the tool mesh
Nb | number of cells in the base mesh
Ev | volume conservation error
Ti | initialization time, loading the mesh and fields
Te | execution time of the mesh intersection operation
Tx | intersection time 
Nx | total number of cells that are intersected 
Ax | average number of intersections per intersected cell 
Ni | number of interface cells,
Nb | number of bulk cells.

# Running 

## Cluster with the SLURM workload manager

### Prepare the input data 

    slurmCloneAndMesh cciDomain # CCI intersection 
    slurmCloneAndMesh smciDomain # SMCI intersection

CCI : cell-cell domain intersection
SMCI : surface mesh - cell domain intersection

### Submit the jobs


    ./slurmSMCIintersect

# Serial  

Run 

    serialCCintersect

or

    serialSMCintersect

# Result processing 

    python3 2D-init-error-analyze.py
