# Overview

This is a parameter study for intersecting a 3D sphere mesh with a 3D volume mesh in OpenFOAM (SMCI or CCMI). A third approach approximates volume fractions based on
the signed distance computed from the sphere mesh in combination with tetrahedral
cell refinement.

You can run the study serially, or on a cluster that is using the SLURM workload manager.

Measured parameter  | Description
--- | ---
Nt | number of cells or triangles in the tool mesh
Nb | number of cells in the base mesh
Ev | volume conservation error
Ti | initialization time, loading the mesh and fields
Te | execution time of the mesh intersection operation
Nx | total number of cell / halfspace intersections 
Ni | number of interface cells,
Nk | number of bulk cells.
Va | volume as approximated by the VoF field
Ed | volume error relative to the discrete surface volume

# Running 

**Note**:
* ensure you have build, installed and sourced the `geom-vof-init`
project properly. 
* ensure the OpenFOAM environment is set (aka `etc/bashrc` is sourced.)

Create the `meshed-sphere.vtk` and `cfmesh-sphere.stl` file in the
toolCase/ directory (requires GMsh)

    ./makeSurfaceMesh

Create the 3D tool mesh using cfMesh (only required for CCMI)

    toolCase> cartesianMesh

Clone the cases for the parameter study in serial

    ./00-serial-clone-mesh-study algorithm parameterFile

or parallel

    ./00-SLURM-clone-mesh-study algorithm parameterFile

where *algorithm* is either `smci`, `ccmi` or `povof` and *parameterFile*
is the name of the studies corresponding parameter file. Note that there are
three variants of a parameter file, one for each mesh type *hexahedral*,
*tetrahedral* and *polyhedral*.

Run the parameter variations in serial

    ./01-serial-povof parameterFile

or parallel using one of the following:

    ./01-SLURM-CCMI-intersect
    ./01-SLURM-SMCI-intersect

Note: a script using SLURM for poVof still needs to be implemented.

# Result processing 
Evaluation and visualization of the results is done using a Jupyter notebook
named *analyze-volume-fraction-initialization.ipynb*.
