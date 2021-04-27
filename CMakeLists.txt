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
set(of_app_path "${of_wm_project_dir}/applications" CACHE PATH "Path to OpenFOAM/applications folder")

message(STATUS "OpenFOAM lib path: ${of_lib_path}")
message(STATUS "OpenFOAM src path: ${of_src_path}")
message(STATUS "OpenFOAM applications path: ${of_app_path}")

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

# Required OpenFOAM libraries for interIsoPandoraFoam
find_library(OF_FVOPTIONS fvOptions PATHS ${of_lib_path})
find_library(OF_SAMPLING sampling PATHS ${of_lib_path})
find_library(OF_DYNAMICFVMESH dynamicFvMesh PATHS ${of_lib_path})
find_library(OF_GEOMETRICVOF geometricVoF PATHS ${of_lib_path})
find_library(OF_INCOMPRESSIBLETRANSPORTMODELS incompressibleTransportModels PATHS ${of_lib_path})
find_library(OF_INTERFACEPROPERTIES interfaceProperties PATHS ${of_lib_path})
find_library(OF_IMMISCIBLEINCOMPRESSIBLETWOPHASEMIXTURE immiscibleIncompressibleTwoPhaseMixture PATHS ${of_lib_path})
find_library(OF_TURBULENCEMODELS turbulenceModels PATHS ${of_lib_path})
find_library(OF_INCOMPRESSIBLETURBULENCEMODELS incompressibleTransportModels PATHS ${of_lib_path})
find_library(OF_WAVEMODELS waveModels PATHS ${of_lib_path})

# Use lib instead of lib64 for libraries.
set(CMAKE_INSTALL_LIBDIR "$ENV{FOAM_USER_LIBBIN}")
set(CMAKE_INSTALL_BINDIR "$ENV{FOAM_USER_APPBIN}")

# Build the library
# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which type is required / reasonable for this project (TT)
add_library("triSurfaceImmersion" SHARED
    src/triSurfaceImmersion/geomMeshIntersection/geomMeshIntersection/geomMeshIntersection.cpp
    src/triSurfaceImmersion/geomMeshIntersection/geomTriSurfaceTools/geomTriSurfaceTools.cpp
    src/triSurfaceImmersion/implicitSurfaces/implicitSurfaces.cpp
    src/triSurfaceImmersion/polynomialVofInitialization/tetVofCalculator/tetVofCalculator.cpp
    src/triSurfaceImmersion/polynomialVofInitialization/orientedSurfaces/orientedPlane.cpp
    src/triSurfaceImmersion/polynomialVofInitialization/orientedSurfaces/triSurfaceAdapter.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/insideOutsidePropagation/insideOutsidePropagation.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/searchDistanceCalculator/searchDistanceCalculator.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/levelSetDistCalc/levelSetDistCalc.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/signedDistanceCalculator/signedDistanceCalculator.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/triSurfaceDistCalc/triSurfaceDistCalc.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/volumeFractionCalculator/volumeFractionCalculator.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellApproximation/surfaceMeshCellApproximation.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellIntersection/surfaceMeshCellIntersection.cpp
)

target_include_directories("triSurfaceImmersion" PUBLIC
    src/triSurfaceImmersion/geomMeshIntersection/geomMeshIntersection
    src/triSurfaceImmersion/geomMeshIntersection/geomReconstructError
    src/triSurfaceImmersion/geomMeshIntersection/geomTriSurfaceTools
    src/triSurfaceImmersion/implicitSurfaces
    src/triSurfaceImmersion/include
    src/triSurfaceImmersion/polynomialVofInitialization/AdaptiveTetCellRefinement
    src/triSurfaceImmersion/polynomialVofInitialization/tetVofCalculator
    src/triSurfaceImmersion/polynomialVofInitialization/orientedSurfaces
    src/triSurfaceImmersion/signedDistanceCalculation/insideOutsidePropagation
    src/triSurfaceImmersion/signedDistanceCalculation/searchDistanceCalculator
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/levelSetDistCalc
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/signedDistanceCalculator
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/triSurfaceDistCalc
    src/triSurfaceImmersion/volumeFractionCalculators/volumeFractionCalculator
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellApproximation
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellIntersection
)

target_include_directories("triSurfaceImmersion" PUBLIC SYSTEM
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)

target_link_libraries("triSurfaceImmersion"
    geophase
    ${OpenMP_CXX_LIBRARIES}
)


add_library("pandora" SHARED
    src/pandora/curvatureFunctionObject/curvatureFunctionObject.cpp
    src/pandora/pandora/pandora.cpp
    src/pandora/pandoraCurvature/pandoraCurvature/pandoraCurvature.cpp
    src/pandora/pandoraCurvature/pandoraDivGradCurvature/pandoraDivGradCurvature.cpp
    src/pandora/pandoraCurvature/pandoraDivNormalCurvature/pandoraDivNormalCurvature.cpp
    src/pandora/pandoraCurvature/pandoraSmoothedMarkerCurvature/pandoraSmoothedMarkerCurvature.cpp
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureAverageExtension/pandoraCurvatureAverageExtension.cpp
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureExtension/pandoraCurvatureExtension.cpp
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureNoExtension/pandoraCurvatureNoExtension.cpp
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureLaplaceRegularisation/pandoraCurvatureLaplaceRegularisation.cpp
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureNoRegularisation/pandoraCurvatureNoRegularisation.cpp
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureNormalRegularisation/pandoraCurvatureNormalRegularisation.cpp
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureRegularisation/pandoraCurvatureRegularisation.cpp
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureTangentialAverageRegularisation/pandoraCurvatureTangentialAverageRegularisation.cpp
)

target_include_directories("pandora" PUBLIC
    src/pandora/curvatureFunctionObject
    src/pandora/pandora
    src/pandora/pandoraCurvature/pandoraCurvature
    src/pandora/pandoraCurvature/pandoraDivGradCurvature
    src/pandora/pandoraCurvature/pandoraDivNormalCurvature
    src/pandora/pandoraCurvature/pandoraSmoothedMarkerCurvature
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureAverageExtension
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureExtension
    src/pandora/pandoraCurvatureExtension/pandoraCurvatureNoExtension
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureLaplaceRegularisation
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureNoRegularisation
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureNormalRegularisation
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureRegularisation
    src/pandora/pandoraCurvatureRegularisation/pandoraCurvatureTangentialAverageRegularisation
)

target_include_directories("pandora" PUBLIC SYSTEM
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)


# TODO: enable if OPenMP is re-enabled (TT)
#target_compile_options("triSurfaceImmersion" PRIVATE "-fopenmp")

# Install the libraries 
install(TARGETS
    "pandora"
    "triSurfaceImmersion"
    LIBRARY
)


# Build executables
add_executable("geophaseUnitTest"
                applications/geophaseUnitTest/geophaseUnitTest.cpp)
target_include_directories("geophaseUnitTest" PRIVATE 
    src/triSurfaceImmersion/include
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

# Preprocessing
add_executable("surfaceInitVolumeFraction"
               applications/surfaceInitVolumeFraction/surfaceInitVolumeFraction.cpp)
target_link_libraries("surfaceInitVolumeFraction"
    "triSurfaceImmersion"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("surfaceInitSignedDistances"
               applications/surfaceInitSignedDistances/surfaceInitSignedDistances.cpp)
target_link_libraries("surfaceInitSignedDistances"
    "triSurfaceImmersion"
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

# Pandora executables
add_executable("interIsoPandoraFoam"
    applications/interIsoPandoraFoam/interIsoPandoraFoam.cpp)
target_include_directories("interIsoPandoraFoam" PRIVATE
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/sampling/lnInclude
    ${of_src_path}/dynamicFvMesh/lnInclude
    ${of_src_path}/dynamicMesh/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
    ${of_src_path}/transportModels
    ${of_src_path}/transportModels/geometricVoF/lnInclude
    ${of_src_path}/transportModels/twoPhaseMixture/lnInclude
    ${of_src_path}/transportModels/incompressible/lnInclude
    ${of_src_path}/transportModels/interfaceProperties/lnInclude
    ${of_src_path}/TurbulenceModels/turbulenceModels/lnInclude
    ${of_src_path}/TurbulenceModels/incompressible/lnInclude
    ${of_src_path}/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/transportModels/geometricVoF/lnInclude
    ${of_app_path}/solvers/multiphase/VoF
    ${of_app_path}/solvers/multiphase/interFoam
)
target_link_libraries("interIsoPandoraFoam"
    "pandora"
    ${OF_FINITEVOLUME}
    ${OF_FVOPTIONS}
    ${OF_MESHTOOLS}
    ${OF_SAMPLING}
    ${OF_DYNAMICFVMESH}
    ${OF_GEOMETRICVOF}
    ${OF_INCOMPRESSIBLETRANSPORTMODELS}
    ${OF_INTERFACEPROPERTIES}
    ${OF_IMMISCIBLEINCOMPRESSIBLETWOPHASEMIXTURE}
    ${OF_TURBULENCEMODELS}
    ${OF_INCOMPRESSIBLETURBULENCEMODELS}
    ${OF_WAVEMODELS}
    ${OF_SURFMESH}
    dl
    m
)

add_executable("pandoraTestCurvature"
    applications/pandoraTestCurvature/pandoraTestCurvature.cpp
)
target_include_directories("pandoraTestCurvature" PRIVATE
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/transportModels/geometricVoF/lnInclude
)
target_link_libraries("pandoraTestCurvature"
    "pandora"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_SURFMESH}
    ${OF_GEOMETRICVOF}
)

# Install executables
install(TARGETS
        # Pandora
        "interIsoPandoraFoam"
        "pandoraTestCurvature"
        # TriSurfaceImmersion
        # Refactored
        "surfaceInitVolumeFraction"
        "surfaceInitSignedDistances"
        # Pre-refactor
        "geophaseUnitTest"
        "foamPerturbMesh"
        "computeLiquidVolume"
        RUNTIME)
