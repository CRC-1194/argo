# Overview

This directory contains the necessary files to prepare, run and evaluate parameter studies
with either a spherical or ellipsoidal interface. Such studies allow to analyze the 
convergence behaviour of the relative global volume error with respect to interface resolution for the 
SMCI/A algorithm.

You can run a study serially, or on a cluster that is using the SLURM workload manager.

# Running 

**Prerequisites**:
* ensure you have build, installed and sourced the `argo` project properly. 
* ensure the OpenFOAM environment is set (aka `etc/bashrc` is sourced.)
* ensure Gmsh, PyFoam and Jupyter Notebook are installed on your system

## Conducting parameter studies

### Reproduce publication results
To reproduce the plots and their corresponding data from the publication run
```
./reproduce_publication_results.sh
```
This will also run the SMCA refinement convergence study for the CAD geometry.
The output in form of CSV and PDF files is written to `./results`. You can specify another directory
by setting the environment variable `GEOM_VOF_INIT` to the desired path. For a concise summary of
the results in one place, open and run the Jupyter notebooks `smci-vof-init.ipynb` and 
`smca-vof-init.ipynb`.  
To remove the generated directories and files run `./clean_studies_and results.sh`.

### Create and run specific studies
Specific studies are created with
```
./create_parameter_study.py ...
```
and then run by
```
./run_parameter_study ...
```
where `...` needs to be replaced with appropriate arguments. Please call each script with the `--help` option
to get a list of available arguments. You can also have a look at the scripts
`create_smci_smca_verification_study.sh` and `run_smci_smca_verification_study.sh` to see examples on
how to use them.
An overview of the available studies and study parameters is given further below.

### Study parameters and available studies
Two interface shapes are available, a _sphere_ and an _ellipsoid_. Their parameters can be found and changed
in the `templateCase/*.geo.template` files.  
Three additional parameters can be set in the `*.parameter` file:
* `N_CELLS`: controls the number volume mesh cells. For the `blockMesh` mesh generator, the resulting number
    of cells is $n_\text{cells}^3$.
* `TRIANGLE_EDGE_LENGTH`: controls the triangle size of the surface mesh by setting Gmsh's 
    `Mesh.CharacteristicLengthMax` in `templateCase/*.geo`.
* `MAX_REFINEMENT_LEVEL` (only used by SMCA algorithm, but must be present): sets the maximum refinement level used during tetrahedral
    refinement. A value of `-1` enables the automatic computation of the maximum refinement level.

Currently, the following studies are available:
* `blockMesh.parameter`: examine algorithm convergence with surface mesh resolution for different
    volume mesh resolutions. Use with the `blockMesh` mesh generator.
* `blockMeshSmokeTest.parameter`: same idea as `blockMesh.parameter`, but with reduced number of parameter combinations for faster testing.
* `polyMesh.parameter`: examine algorithm convergence with surface mesh resolution for different
    volume mesh resolutions. Use with the `pMesh` mesh generator.
* `tetMesh.parameter`: examine algorithm convergence with surface mesh resolution for different
    volume mesh resolutions. Use with the `blockMesh` mesh generator.
* `smca-refinement-convergence-blockMesh.parameter`: examine the effectiveness of tetrahedral refinement for the SMCA algorithm by varying the maximum refinement level with fixed surface and volume mesh resolution.


### Postprocessing / evaluation
Postprocsessing is performed by two Jupyter notebooks, namely `smci-vof-init.ipynb` and `smca-vof-init.ipynb`.
Shared functionality, e.g. data aggregation and plotting, is defined in `plot_study.py`.
