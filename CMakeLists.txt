include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(FetchContent)

cmake_minimum_required(VERSION 3.13)

project(geom-vof-init VERSION 1.0
        DESCRIPTION "Computing volume fractions by intersecting unstructured meshes"
        LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build dependencies

# - Dependency: geophase
FetchContent_Declare(geophase
    GIT_REPOSITORY          git@git.rwth-aachen.de:leia/geophase.git
    # TODO: Replace the branch with a new GIT_TAG for geophase as soon as everything runs.  
    #       geophase modifications are needed to make geometrical init work, so the geom-vof-init 
    #       branch is used. 
    #GIT_TAG                 v1.0 
    GIT_TAG                 master 
    GIT_SHALLOW             YES
    GIT_PROGRESS            YES
    USES_TERMINAL_DOWNLOAD  YES
    USES_TERMINAL_UPDATE    YES
)

FetchContent_GetProperties(geophase)
if(NOT geophase_POPULATED)
    FetchContent_Populate(geophase)
endif()

add_library(geophase INTERFACE) 
target_include_directories(geophase INTERFACE 
    ${geophase_SOURCE_DIR}/include/ 
    ${geophase_SOURCE_DIR}/include/core/
    ${geophase_SOURCE_DIR}/include/io/
    ${geophase_SOURCE_DIR}/include/io/vtk/
    ${geophase_SOURCE_DIR}/include/testing/
    ${geophase_SOURCE_DIR}/include/geometrical_algorithm/
    ${geophase_SOURCE_DIR}/include/geometrical_model/
)
# OpenMP support
find_package(OpenMP REQUIRED)

# Read build relevant OpenFOAM environment variables and cache them
set(of_wm_project_dir "$ENV{WM_PROJECT_DIR}" CACHE PATH "Path to OpenFOAM project folder.")
set(of_wm_arch "$ENV{WM_ARCH}" CACHE STRING "Architecture. Usually linux64.")
set(of_wm_arch_option "$ENV{WM_ARCH_OPTION}" CACHE STRING "Information if 32 or 64 bit operating system.")
set(of_wm_precision_option "$ENV{WM_PRECISION_OPTION}" CACHE STRING "Flag if to use single precision (SP) or double precision (DP).")
set(of_wm_label_size "$ENV{WM_LABEL_SIZE}" CACHE STRING "Size in bit to use as label type. Can be either 32 or 64.")
set(of_wm_compile_option "$ENV{WM_COMPILE_OPTION}" CACHE STRING "OpenFOAM build type: Opt, Debug, Prof.")
set(of_wm_compiler "$ENV{WM_COMPILER}" CACHE STRING "Compiler used for OpenFOAM build.")
set(of_wm_label_option "$ENV{WM_LABEL_OPTION}" CACHE STRING "Concrete Type used for label. Either Int32 or Int64.")
# WM_ARCH + WM_COMPILER + WM_PRECISION_OPTION + WM_LABEL_OPTION + WM_COMPILE_OPTION
set(of_wm_options "${of_wm_arch}${of_wm_compiler}${of_wm_precision_option}${of_wm_label_option}" CACHE STRING "Name of subfolder which contains compiled OpenFOAM libraries" FORCE)

# Determine whether specific build exists
# TODO: change build path according to build option passed to cmake command (TT)
if(IS_DIRECTORY "${of_wm_project_dir}/platforms/${of_wm_options}${of_wm_compile_option}")
    set(of_wm_options "${of_wm_options}${of_wm_compile_option}")
else()
    message(FATAL_ERROR "Path ${of_wm_project_dir}/platforms/${of_wm_options}${of_wm_compile_option} does not exist. Make sure OpenFOAM environment is set.")
endif()

set(of_lib_path "${of_wm_project_dir}/platforms/${of_wm_options}/lib" CACHE PATH "Path to compiled OpenFOAM libraries.")
set(of_src_path "${of_wm_project_dir}/src" CACHE PATH "Path to OpenFOAM/src folder")

message(STATUS "OpenFOAM lib path: ${of_lib_path}")
message(STATUS "OpenFOAM src path: ${of_src_path}")

# TODO: Check that those variables are set. Optional: check for reasonable values (TT)

# Compile definitions required for OpenFOAM
add_compile_definitions(
    WM_LABEL_SIZE=${of_wm_label_size}
    WM_${of_wm_precision_option}
    WM_ARCH_OPTION=${of_wm_arch_option}
    ${of_wm_arch}
    OPENFOAM="$ENV{WM_PROJECT_VERSION}" # Figures out OF version. TM. 
    NoRepository
)

# Required to make linking to OpenFOAM libraries work
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed")

if(MSVC)
    add_compile_options(/W4 /WX)
else(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wno-deprecated -Wno-attributes -Wno-unused-variable -Wno-unused-parameter -D_USE_MATH_DEFINES")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -Wall -Wextra -DDEBUG -pedantic -D_USE_MATH_DEFINES -DBOOST_MATH_INSTRUMENT")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb3 -O3 -Wall -Wextra -pedantic -DDEBUG -D_USE_MATH_DEFINES")
endif(MSVC)

# Required OpenFOAM libraries
find_library(OF_FINITEVOLUME finiteVolume PATHS ${of_lib_path})
find_library(OF_MESHTOOLS meshTools PATHS ${of_lib_path})
find_library(OF_OPENFOAM OpenFOAM PATHS ${of_lib_path})
find_library(OF_SURFMESH surfMesh PATHS ${of_lib_path})

# Use lib instead of lib64 for libraries.
set(CMAKE_INSTALL_LIBDIR "$ENV{FOAM_USER_LIBBIN}")
set(CMAKE_INSTALL_BINDIR "$ENV{FOAM_USER_APPBIN}")

# Build the library
# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which type is required / reasonable for this project (TT)
add_library("signedDistanceCalculator" SHARED
    src/signedDistanceCalculation/signedDistanceCalculator/signedDistanceCalculator.cpp
    src/signedDistanceCalculation/insideOutsidePropagation/insideOutsidePropagation.cpp
)

target_include_directories("signedDistanceCalculator" PUBLIC
    src/signedDistanceCalculation/signedDistanceCalculator
    src/signedDistanceCalculation/insideOutsidePropagation
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)

add_library("geometricalMeshIntersection" SHARED
    src/geomMeshIntersection/geomMeshIntersection/geomMeshIntersection.cpp
    src/geomMeshIntersection/geomSurfaceCellMeshIntersection/geomSurfaceCellMeshIntersection.cpp
    src/geomMeshIntersection/geomTriSurfaceTools/geomTriSurfaceTools.cpp
)

target_link_libraries("geometricalMeshIntersection" geophase)

target_include_directories("geometricalMeshIntersection" PUBLIC
    src/geomMeshIntersection/geomMeshIntersection
    src/geomMeshIntersection/geomReconstructError
    src/geomMeshIntersection/geomSurfaceCellMeshIntersection
    src/geomMeshIntersection/geomTriSurfaceTools
    src/include
)
# Using PUBLIC here removes the need to specify again the required include directories
# for the applications.
target_include_directories("geometricalMeshIntersection" SYSTEM PUBLIC
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)

add_library("polynomialVofInitialization" SHARED
    src/polynomialVofInitialization/polynomialVofInitialization/polynomialVofInitialization.cpp
    src/polynomialVofInitialization/tetVofCalculator/tetVofCalculator.cpp
    src/polynomialVofInitialization/orientedSurfaces/orientedPlane.cpp
    src/polynomialVofInitialization/orientedSurfaces/triSurfaceAdapter.cpp
)
target_compile_options("polynomialVofInitialization" PRIVATE "-fopenmp")
target_include_directories("polynomialVofInitialization" PUBLIC
    src/polynomialVofInitialization/polynomialVofInitialization
    src/polynomialVofInitialization/AdaptiveTetCellRefinement
    src/polynomialVofInitialization/tetVofCalculator
    src/polynomialVofInitialization/orientedSurfaces
)
target_include_directories("polynomialVofInitialization" SYSTEM PUBLIC
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
target_link_libraries("polynomialVofInitialization"
    "signedDistanceCalculator"
    ${OpenMP_CXX_LIBRARIES}
)

# There is nothing to link, it is a header library: check what CMake link does in this case. TM.
#target_link_libraries("geometricalMeshIntersection${PROJECT_BUILD_TYPE}" geophase)

# Install the libraries 
install(TARGETS "signedDistanceCalculator" LIBRARY)
install(TARGETS "geometricalMeshIntersection" LIBRARY)
install(TARGETS "polynomialVofInitialization" LIBRARY)

# Build executables

# Tests
add_executable("cellCellTestVofInit"
               applications/cellCellTestVofInit/cellCellTestVofInit.cpp) 
target_link_libraries("cellCellTestVofInit"
    "geometricalMeshIntersection"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("surfaceCellTestVofInit" 
               applications/surfaceCellTestVofInit/surfaceCellTestVofInit.cpp)
target_link_libraries("surfaceCellTestVofInit"
    "geometricalMeshIntersection"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("geophaseUnitTest"
                applications/geophaseUnitTest/geophaseUnitTest.cpp)
target_include_directories("geophaseUnitTest" PRIVATE 
    src/include
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
target_link_libraries("geophaseUnitTest" 
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    geophase
    dl
    m
)
add_test(NAME unitTests COMMAND geophaseUnitTest)
set_tests_properties(unitTests PROPERTIES LABELS "UNIT")

add_executable("polynomialVofTest${PROJECT_BUILD_TYPE}"
               applications/polynomialVofTest/polynomialVofTest.cpp)
target_link_libraries("polynomialVofTest${PROJECT_BUILD_TYPE}"
    "polynomialVofInitialization${PROJECT_BUILD_TYPE}"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

# Preprocessing
add_executable("cellCellVofInit"
               applications/cellCellVofInit/cellCellVofInit.cpp)
target_link_libraries("cellCellVofInit"
    "geometricalMeshIntersection"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("surfaceCellVofInit"
               applications/surfaceCellVofInit/surfaceCellVofInit.cpp)
target_link_libraries("surfaceCellVofInit"
    "geometricalMeshIntersection"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("smcaVofInit"
               applications/smcaVofInit/smcaVofInit.cpp)

target_link_libraries("smcaVofInit"
    "geometricalMeshIntersection"
    "polynomialVofInitialization"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("signedDistanceTestCalc"
               applications/signedDistanceTestCalc/signedDistanceTestCalc.cpp)
target_include_directories("signedDistanceTestCalc" PRIVATE
    applications/signedDistanceTestCalc
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/fileFormats/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
target_link_libraries("signedDistanceTestCalc"
    "signedDistanceCalculator"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("foamPerturbMesh"
               applications/foamPerturbMesh/foamPerturbMesh.cpp)
target_include_directories("foamPerturbMesh" PRIVATE 
    src/include
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
target_link_libraries("foamPerturbMesh"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("computeLiquidVolume"
    applications/computeLiquidVolume/computeLiquidVolume.cpp)
target_include_directories("computeLiquidVolume" PRIVATE 
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
target_link_libraries("computeLiquidVolume"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    ${OF_SURFMESH}
    dl
    m
)

# Install executables
install(TARGETS 
        "cellCellTestVofInit"
        "surfaceCellTestVofInit"
        "surfaceCellVofInit"  
        "cellCellVofInit"
        "geophaseUnitTest"
        "smcaVofInit"
        "signedDistanceTestCalc"
        "polynomialVofTest"
        "foamPerturbMesh"
        "computeLiquidVolume"
        RUNTIME)
