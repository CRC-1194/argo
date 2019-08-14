# geom-vof-init

Implementation of the Surface-Cell Mesh Intersection (SMCI) and Cell-Cell Mesh Intersection (CCMI) algorithms for computing volume fractions by intersecting unstructured meshes. The SCMI algorithm computes volume fractions by intersecting a surface mesh with a volume mesh. The CCMI algorithm intersects cells of two volume meshes with each other.  

## Authors

* **Tomislav Maric** - *Development* - [MMA, TU Darmstadt](https://www.mma.tu-darmstadt.de/index/mitarbeiter_3/mitarbeiter_details_mma_43648.en.jsp)

* **Tobias Tolle** - *Development* - [MMA, TU Darmstadt](https://www.mathematik.tu-darmstadt.de/fb/personal/details/tobias_tolle.de.jsp)

## License

This project is licensed under the GPL3.0 License - see the [LICENSE.md](LICENSE.md) file for details.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

`geom-vof-init` is a project that is used together with OpenFOAM. It compiles on its own and links against OpenFOAM libraries. It consists of a library for geometrical intersections between mesh data structures in OpenFOAM: triSurfaceMesh (triangle surface mesh) and fvMesh (polyhedral volume mesh). 

### Prerequisites

List of prerequisites with tested versions in brackets:

* g++   (9.1.0)
* CMake (3.13) 

These you can install on your system using a package manager.

* OpenFOAM-plus (v1906)

To install OpenFOAM-plus, tag (version) v1906 follow the [instructions on installing OpenFOAM v1906 from sources](https://develop.openfoam.com/Development/OpenFOAM-plus).

### Installing

`geom-vof-init` is using the [CMake](https://cmake.org) build system.  

Execute following commands to build `geom-vof-init`, once you have installed OpenFOAM-plus (v1906). 

Inside the `geom-vof-init` directory


```
?>  mkdir build && cd build 
?>  cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release ..
?>  make && make install
```

Unlike OpenFOAM, CMake will save the `geom-vof-init` executable and library files in the `build` folder.

The next step is to expand the `PATH` variabe so that the executables are found and to expand the `LD_LIBRARY_PATH` so that the linker can find the mesh intersection library. 

Add something along the lines of 


```
export PATH="/path/to/your/own/geom-vof-init/build/bin":$PATH
export LD_LIBRARY_PATH="/path/to/your/own/geom-vof-init/build/lib":$LD_LIBRARY_PATH

```  

to your `.bashrc` file, then execute

```
?> source $HOME/.bashrc
```

to finish installing `geom-vof-init`. 

## Running the tests 

TODO.

## Running the applications 

TODO.

## Contributing

The code is maintained at [TU-GitLab](https://git.rwth-aachen.de/leia/geom-vof-init). Feedback in the form of contributions, bug reports or feature requests is handled there. Users external to the German TU-GitLab network can login using their github.com credentials. 
