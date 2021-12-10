# argo  

The "argo" project is an [OpenFOAM](https://develop.openfoam.com/Development/openfoam) module that implements unstructured Lagrangian / Eulerian Interface (LEIA) methods for multiphase flow simulations in complex geometries.

## Authors

* **Tomislav Maric** - *Development* - [MMA, TU Darmstadt](https://www.mma.tu-darmstadt.de/index/mitarbeiter_3/mitarbeiter_details_mma_43648.en.jsp)

* **Tobias Tolle** - *Development* - [MMA, TU Darmstadt](https://www.mathematik.tu-darmstadt.de/fb/personal/details/tobias_tolle.de.jsp)

* **Dirk Gründing** - *Development* - [MMA, TU Darmstadt](https://www.mma.tu-darmstadt.de/index/mitarbeiter_3/mitarbeiter_details_mma_47488.en.jsp)

## Publications 

[1] [Tolle, T., Gründing, D., Bothe, D., & Marić, T. (2021). Computing volume fractions and signed distances from arbitrary surfaces on unstructured meshes. arXiv preprint arXiv:2101.08511.](https://arxiv.org/abs/2101.08511)

[2] [Hartmann, M., Fricke, M., Weimar, L., Gründing, D., Marić, T., Bothe, D., & Hardt, S. (2021). Breakup Dynamics of Capillary Bridges on Hydrophobic Stripes. International Journal of Multiphase Flow, 103582.](https://doi.org/10.1016/j.ijmultiphaseflow.2021.103582)

## License

This project is licensed under the GPL3.0 License - see the [LICENSE.md](LICENSE.md) file for details.

## Installation

These instructions will get your copy of the project up and running on your local machine for development and testing purposes. 

`argo` is a project that builds on [OpenFOAM](https://develop.openfoam.com/Development/openfoam) so it compiles and links against OpenFOAM libraries.  

### Compilation & Installation dependencies 

* Compiler:  g++ (GCC) 10.2.0
* Build system: cmake version 3.19.3

#### Computing dependencies

Meshing 

* [gmsh](http://gmsh.info/) meshing software version 4.7.1, used for generating surface meshes
* [cfmesh](https://cfmesh.com/cfmesh/), available as OpenFOAM sub-module, used for automatic generation of unstructured volume meshes

OpenFOAM

`argo` is based on OpenFOAM, git tag OpenFOAM-v2012

To install OpenFOAM follow the [instructions on installing OpenFOAM from sources](https://develop.openfoam.com/Development/openfoam/). 

1. Check out openfoam using git. 
2. Check out the git tag 

```
    ?> git checkout OpenFOAM-v2012
```
3. Compile OpenFOAM as instructed by its documentation. Make sure you to compile OpenFOAM with
`c++17` by changing `CC = g++ -std=c++11` to `CC = g++ -std=c++17` in the
file `wmake/rules/General/Gcc/c++`.

#### Post-processing dependencies

We use [Jupyter notebooks](https://jupyter.org/) for visualization and processing of test results, and following packages (may be differently named on your Operating System) 

* python, python-pandas, python-numpy, python-jupyter

### Installing `argo`

`argo` is built using the [CMake](https://cmake.org) build system.  

Execute following command to build `argo` and install its libraries and executables in the OpenFOAM PATH structure, once you have installed all its dependencies listed above. Inside the `argo` directory, call

```
    ?> ./install.sh
```

#### Manual compilation 

`argo` can be built by directly calling `cmake`

```
?>  mkdir build && cd build 
?>  cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
?>  make && make install
```

where the flag `-DCMAKE_EXPORT_COMPILE_COMMANDS=on` is optional (it instructs CMake to create a `compile_commands.json` file).

## Tutorial: computation of volume fractions for a brain surface
This tutorial guides you through the necessary steps to compute volume fractions (and signed distances) for
a given surface mesh. Prerequisite is that you have a working installation of OpenFOAM and Argo. The result
is very similar to the initialization test case `cases/initialization/3D-brain`.

### 1. Obtain a suitable surface mesh
Here we use a surface mesh of a brain that has been obtained from a MRI brain scan which you can download from
[Thingiverse](https://www.thingiverse.com/thing:1287025) under _Download All Files_. In the archive, you
can find the file `files/andrewbrain_1.stl` which we'll use as input surface. Before this can be used as input
for the SMCIA algorithm, two conditions have to be ensured:
* The triangle normals are oriented consistently, meaning they all point either inside or outside.
* All triangles have an area greater than zero, meaning that there are no coinciding vertices.

If the conditions above are not fulfilled, SMCIA may behave strangely and/or simply crash.
Further properties of the surface mesh which are favorable for SMCIA to work correctly are:
* There are no duplicate triangles.
* The surface has no holes.

In principle SMCIA should still work correctly in the presence of duplicate triangles, but this has not been
tested yet. Even with holes in the mesh SMCI should still correctly compute inside/outside information, but
the absolute distance to the surface will be erroneous.  
In this tutorial some OpenFOAM tools are used to prepare the surface mesh, but other tools for mesh
manipulation and clean up, e.g.
[PyMesh](https://pymesh.readthedocs.io/en/latest/index.html), can be used.

#### 1.1 Check the quality of the surface mesh
**Before you start**: make sure you have sourced the OpenFOAM environment.  

OpenFOAM provides the `surfaceCheck` utility to assess the quality of a surface mesh. Run this with the
`-verbose` option for more detailled output:
```
surfaceCheck -verbose andrewbrain_1.stl
```
The most important pieces of information are:
* whether the surface has any illegal faces, e.g. `Surface has 21 illegal triangles.`
* whether the minimum triangle quality is 0, e.g. `Minimum triangle quality is 0`. This means there
    are collapsed triangles with an area of 0.

If there are no illegal faces nor collapsed triangles, the mesh quality is fine for SMCIA.
However, for our example mesh, the output shows that the mesh cannot be used as it is.

#### 1.2 Fix the surface mesh
OpenFOAM offers several tools to manipulate and clean up surface meshes. All commands related to these
surface operations start with `surface`.  
As you may have already noticed, `surfaceCheck` has written several `andrewbrain_1_*.obj` files. These are
different parts of the mesh which are not connected. In the following, `andrewbrain_1_0.obj` is used which
comprises almost all elements of the original surface and is now a single connected surface.  
Running 
```
surfaceCheck andrewbrain_1_0.obj
```
shows that there are still problems with this surface. This can be fixed with `surfaceClean` which - 
according to its description - does _Clean surface by removing baffles, sliver faces, collapsing small
edges, etc._. Besides the input and ouput file name, the required arguments are the minimal edge length
and the minimal triangle quality. Both can be set to quite low values (as long as a normal can be
computed, a triangle is fine here.):
```
surfaceClean andrewbrain_1_0.obj 1e-12 1e-6 clean_surface.obj
```
This gives a suitable surface in terms of triangle sizes. The final step of the surface preparation is
to check the consistent normal orientation.

#### 1.3 Ensure consistent normal orientation
A visual inspection of surface normals can be done with Paraview and the `Normal Glyphs` filter. In case the
normals are inconsistent, OpenFOAM's `surfaceOrient` tool can fix the orientation. In addition to the surface
itself a point located on the outside of the surface is required. `surfaceCheck` also outputs the bounding box
of the surface which can help to determine an outside point. For the surface at hand, the
point `(-100 -100 -100)` is located on the outside:
```
surfaceOrient clean_surface.obj "(-100 -100 -100)" oriented_surface.obj
```

### 2. Copy a case from the initialization cases
As SMCIA builds on OpenFOAM, a valid OpenFOAM case is required to run it. For the sake of simplicity,
you can simply copy one of the example cases, e.g. `cases/initialization/3D-pores`, as only a few changes
are required:
```
cp -r /path/to/Argo/cases/initialization/3D-pores 3D-brain
```

### 3. Configure the volume mesh using OpenFOAM's blockMesh
In Principle, SMCIA works on any mesh that is supported by OpenFOAM. The surface mesh does not even need to be
completely embedded in the volume mesh. For example, if you want to initialize a droplet sitting on a surface
that has the shape of a spherical cap, you can use a complete sphere from which only a part in located inside
the volume mesh.  
Here, again for simplicity, we use a simple block shaped domain which is meshed with
[blockMesh](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility). We re-use a blockMesh config file from an existing test case:
```
cp /path/to/Argo/cases/initialization/3Dinit/templateCase/system/blockMeshDict 3D-brain/system/
```
and adapt the vertices of the block. Based on the bounding box of `oriented_surface.obj` 
```
Bounding Box : (-71.2844 -94.3172 -78.821) (67.1883 86.5888 56.4499)
```
(e.g. given by `surfaceCheck`) we can construct a slightly larger block, e.g.
_(-80 -100 -82) (70 90 62)_ and change the entry _vertices_ in `3D-brain/system/blockMeshDict`
accordingly:
``` 
vertices        
(
    (-80  -100 -82)
    ( 70  -100 -82)
    ( 70    90 -82)
    (-80    90 -82)

    (-80  -100 62)
    ( 70  -100 62)
    ( 70    90 62)
    (-80    90 62)
);
```
The resolution is currently set to 64 in each direchtion (_block_ entry in `blockMeshDict`):
```
blocks          
(
    hex (0 1 2 3 4 5 6 7) (64 64 64) simpleGrading (1 1 1)
);
```
We change this to 100:
```
blocks          
(
    hex (0 1 2 3 4 5 6 7) (100 100 100) simpleGrading (1 1 1)
);
```
Change into the case directory `3D-brain` and run `blockMesh` to create the mesh:
```
cd 3D-brain/
blockMesh
```

### 4. Configure signed distance - / volume fraction calculation
The application for volume fraction computation is `surfaceInitVolumeFraction`, for signed distance
calculation `surfaceInitSignedDistances`. Run each with the `-help` option to show the available
parameters and their purpose. In principle, it is possible to provide all required options on the
command line. But you can also set the options in configuration files or dictionaries as they are called
by OpenFOAM. This allows you to have a persistent configuration. Note that an option can be given in the
configuration file and the command line at the same time. In this case, the command line option has
precendence.  
Volume fraction calculation is configured by the file `system/vofInitDict`. Open it in a text editor to see
the current configuration:
```
fieldName       alpha.water;
algorithm       SMCA;
writeGeometry   off;
invert          on;
writeAllFields  on;
checkVolume     on;
refinementLevel 3;

distCalc
{
    surfaceType triSurface;
    surfaceFile CAB_XD_95VC-clean.stl;
    narrowBandWidth 4.0;
}
```

Here, `fieldName` specifies the name of the volume fraction field. This field in form of a file has to be
present in the `0.org` or `0` directory of the OpenFOAM case. As this was copied from a working case, there
is no need to change it.  

With `algorithm` you can choose between the geometric intersection (SMCI) or approximation in combination
with refinement (SMCA). For this tutorial, SMCA is used.  

Leave the option `writeGeometry` turned off. If active, additional geometric information is written out, e.g.
the decomposition of each interface cell. This is intended for debugging purposes.  

`invert` changes on which side of the interface the volume fractions are set to 1, so the normals of the
surface mesh can be left unchanged. Here, leave it on as the normals of the surface mesh point to the
outside, but we want cells on the inside to have volume fractions of 1.  

If `writeAllFields` is active, all auxilliary fields are also written, meaning especially the signed
distance fields. So use this options if you want both, volume fractions and signed distances.
In this tutorial, we leave the option `on`.  

The flag `checkVolume` is intended for testing purposes and we can set it to `off`.  

`refinementLevel` is only used by the SMCA algorithm and determines how many levels of refinement are used.
Increasing the number of refinement levels increases accuracy, but also required computational time. A value
of `-1` triggers an automatic computation of the refinement level based on the size of the smallest triangle.
For this tutorial, the value `3` is fine.

`distCalc` is a so-called _sub-dictionary_ which configures the signed distance configuration. Here, only
`surfaceFile` has to be adapted to `oriented_surface.obj`, the surface mesh we prepared in step 1.3. The other
parameters can be left untouched.

If you are only interested in signed distance calculation, you can use the options described
above. The configuration dictionary is `system/signedDistanceInitDict`


### 5. Compute signed distances / volume fractions
Before the volume fractions can be computed, the surface mesh has to be placed in the case folder `3D-brain`.
So copy the file from step 1.3:
```
cp /path/to/oriented_surface.obj /path/to/case/3D-brain/
```
Furthermore, a `0` folder is required. You can simply copy the `0.org` folder. Change into the case directory
and copy it:
```
cd /path/to/case/3D-brain
cp -r 0.org 0
```
The case is ready for initilization now:
```
surfaceInitVolumeFraction
```
As the volume fraction computation has been completely configured through its dictionary (step 4), no
command line arguments need to passed to `surfaceInitVolumeFraction`.  
If you want to view the computed fields in ParaView, run
```
foamToVTK
```
This creates VTK files and writes them into the folder `VTK`.


## How to setup and run example cases
Currently, there are two groups of test cases in `argo/cases`:
* `initialization`: cases demonstrating the computation of signed distances and volume fractions from
    triangulated surfaces and level sets given by implicit functions.
* `hydrodynamic`: verification test cases for two-phase flows, e.g. _Stationary Droplet_.

The directory `scripts` contains executable scripts and Python modules to setup and execute parameter studies and collect and visualize results. Only the scripts beginning with `argo-*` are functional at the moment.

### How to use the hydrodynamic test cases
First of all, source `cases/scripts/bashrc` in the Argo folder. This makes the scripts used below available
in your shell. Use each of these scripts with `--help` to get further information on what the script does
and what options are available.

#### Short version
* Go into a test case folder, e.g. `stationaryDroplet2D`:  
```
?> cd cases/hydrodynamic/stationaryDroplet2D
```
* Run `argo-create-parameter-study.py` to create the case directories of a parameter study:
```
?> argo-create-parameter-study.py smoothedMarkerCurvature.parameter
```
* Run `argo-initilize-parameter-study.py` to create the mesh and initialize fields in each case directory:
```
?> argo-initilize-parameter-study.py smoothedMarkerCurvature_0000 -m blockMesh -f initFields.sh 
```
* Run `argo-run-study.py` to execute a solver in each case directory:
```
?> argo-run-study.py interIsoPandoraFoam -d smoothedMarkerCurvature_0000
```
* If data has been written to a `*.csv` or `*.dat` file in each case directory, e.g. by a function object, you can collect this data
with `argo-agglomerate-study-data.py`:
```
?> argo-agglomerate-study-data.py smoothedMarkerCurvature_00000_templateCase/postProcessing/minMaxU/0/fieldMinMax.dat -p smoothedMarkerCurvature.parameter
```
This assembles a multi-indexed Pandas DataFrame and saves it as `smoothedMarkerCurvature.csv`.

#### Additional information to the short version
All test cases within this folder are _templated_ so that parameter studies can be created from them.
So you find the following files and directories within a test case directory, e.g. `stationaryDroplet2D`:
* the templated case `templateCase` which follows OpenFOAM's case structure,
* a `default.parameter` file,
* one or more additional `*.parameter` files.

Using the `templateCase` directly is not possible as some parameters, e.g. the mesh resolution, just contain a
placeholder, e.g. `@!RESOLUTION!@`. These placeholders or parameters are replaced in the parameter study
creation process with values from a study parameter file, e.g. `smoothedMarkerCurvature.parameter` or from
`default.parameter`.  

**Create study case directories**  
To create a case you can actually run use the script `argo-create-parameter-study.py` which uses
`pyFoamRunParameterVariation.py` in the background. If you just what a single case rather than all directories
of a study (which can be quite a lot depending on the parameter file), use the `-v` option. This allows you
to setup only a single variant of a study, e.g. `-v 12` creates only variation number 12 of a study. To find
out which number corresponds to which parameter vector, use
```
pyFoamRunParamaterVariation.py --list-variations templateCase myStudy.parameter
```
Note that for a case directory setup in this way no meshing or preprocessing has been done yet.  

**Mesh and initialize case directories**  
The script `argo-initilize-parameter-study.py` takes care of the case initialization which comprises
three steps:
1. **Mesh creation**: create the mesh with user-prescribed meshing application via option `-m` (Required).  
One can prescribe any meshing application that is present on one's system, e.g. OpenFOAM's
`blockMesh` or `cartesianMesh` from _cfMesh_. However, please check beforehand that the required dictionaries are present.  
2. **Field initialization**: execute a script, e.g. to initialize fields like volume fraction, via option `-f` (Required).  
Each template case contains a script for further preprocessing steps, usually the initialization of the
interface for a two-phase system. A typical name for this script is `initFields.sh`.  
3. **Domain decomposition**: decompose the case for parallel simulation using OpenFOAM's `decomposePar` if `-par` option is present (Optional).

For this purpose, the script iterates over all directories that match a prescribed pattern, e.g. `myAwesomeStudy_00`, and executes each of the three steps listed above.

**Run a simulation/study**  
At this stage, the cases are ready for solver execution. So, if you are only interested in a single case,
just change into the directory and execute a suitable solver. If you want to execute a study with multiple
case directories, use `argo-run-study.py`. This tool iterates over all directories matching a given pattern,
e.g. `myAwesomeStudy_00`, and executes a prescribed solver.  
Each case can be executed in parallel with the
`--use-mpi M_MPI_PROCS` option, where `M_MPI_PROCS` is the number of MPI processes.  
Further more, `N_PROCS` cases can be
executed simultaneously with the option `--num-processes N_PROCS`. Make sure that `N_PROCS * M_MPI_PROCS`
is less or equal the number of available _physical_ cores.

**Using the scripts with SLURM**  
The scripts `argo-initialize-parameter-study.py` and `argo-run-study.py` support the submission of their
workload as jobs to the SLURM workload manager via the option `-j`. The latter expects to find a SLURM 
sbatch script named `solver_to_be_executed.sbatch` which it uses for job submission. For example, if you want
to use the solver `interIsoFoam`, there must be a `interIsoFoam.sbatch` script.


### How to use the initialization test cases
**Note**: _This section is is work in progress._  
There are two initialization applications:
* `surfaceInitVolumeFractions`
* `surfaceInitSignedDistances`

Calling them with the `-help` option gives you an overview of the options you can pass.
Have a look at the `Allrun` scripts present in most of the
cases. Furthermore, you can have a look at `system/vofInitDict` to get an idea what options are available for
the SMCI/A algorithm. 

## Examples 

### Experimental fluid interface

The case with the experimental surface from [[1](https://arxiv.org/abs/2101.08511)] and [[2](https://doi.org/10.1016/j.ijmultiphaseflow.2021.103582)] is available in 

```
argo/cases/3D-SFB1194-A02b
```

There are two scripts, `Allrunscmi` and `Allrunsmca`, that each generate the mesh and compute volume fractions, shown in the below figure 

![Volume fractions initialized from an experimental fluid interface](./etc/experimental-case-smci.png)

## Contributing

The code is maintained on [GitLab](https://gitlab.com/leia-methods/argo). Feedback,  requests to join the project, bug reports or feature requests are handled [via this email](mailto:incoming+leia-methods-argo-23940628-issue-@incoming.gitlab.com).
