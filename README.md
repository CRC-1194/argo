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
