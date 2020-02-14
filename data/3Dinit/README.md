# Overview

This is a parameter study for intersecting a 3D sphere mesh with a 3D volume mesh in OpenFOAM.

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
Nb | number of bulk cells.

# Running 

Create the sphere.stl in the toolCase/ directory

```
    pvpython paraview-create-sphere.py 
```

Create the 3D tool mesh using cfMesh 

```
    toolCase> cartesianMesh
```

Clone the cases for the parameter study 

```
    ./00-SLURM-clone-mesh-study hexDomain smciDomain.parameter
    ./00-SLURM-clone-mesh-study hexDomain ccmiDomain.parameter
```

Run the parameter variations 

```
    ./01-SLURM-CCMI-intersect
    ./01-SLURM-SMCI-intersect
```

# Result processing 

    python3 2D-init-error-analyze.py
