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

The next step is to expand the `PATH` variabe so that the executables are found. Add something along the lines of 


```
export PATH="/path/to/your/own/geom-vof-init/build/bin":$PATH

```  

to your `.bashrc` file, then execute

```
?> source $HOME/.bashrc
```

and `geom-vof-init` is set to be used.  

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

