include(FetchContent)

cmake_minimum_required(VERSION 3.13)

project(Argo VERSION 1.0
        DESCRIPTION "The 'argo' project is an OpenFOAM module that implements unstructured Lagrangian / Eulerian Interface (LEIA) methods for multiphase flow simulations in complex geometries."
        LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

####################################################################################################
#   Fetch and build dependencies
####################################################################################################

# - Dependency: geophase
FetchContent_Declare(geophase
    GIT_REPOSITORY          https://git.rwth-aachen.de/leia/geophase.git
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


####################################################################################################
#   Read and store OpenFOAM environment and configuration
####################################################################################################

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

# Additional librariesrequired for high density ratio solvers
find_library(OF_VOFPHASETURBULENTTRANSPORTMODELS VoFphaseTurbulentTransportModels PATHS ${of_lib_path})
find_library(OF_INCOMPRESSIBLEINTERPHASETRANSPORTMODELS incompressibleInterPhaseTransportModels PATHS ${of_lib_path})

# Use lib instead of lib64 for libraries.
set(CMAKE_INSTALL_LIBDIR "$ENV{FOAM_USER_LIBBIN}")
set(CMAKE_INSTALL_BINDIR "$ENV{FOAM_USER_APPBIN}")


####################################################################################################
#   Argo libraries
####################################################################################################

# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which type is required / reasonable for this project (TT)
add_library("triSurfaceImmersion" SHARED
    src/triSurfaceImmersion/implicitSurfaces/implicitSurfaces.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/insideOutsidePropagation/insideOutsidePropagation.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/searchDistanceCalculator/searchDistanceCalculator.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/levelSetDistCalc/levelSetDistCalc.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/signedDistanceCalculator/signedDistanceCalculator.cpp
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/triSurfaceDistCalc/triSurfaceDistCalc.cpp
    src/triSurfaceImmersion/tetrahedralTools/tetVolumeFractionCalculator/tetVolumeFractionCalculator.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellApproximation/surfaceMeshCellApproximation.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellIntersection/surfaceMeshCellIntersection.cpp
    src/triSurfaceImmersion/volumeFractionCalculators/volumeFractionCalculator/volumeFractionCalculator.cpp
)
target_include_directories("triSurfaceImmersion" PUBLIC
    src/triSurfaceImmersion/implicitSurfaces
    src/triSurfaceImmersion/include
    src/triSurfaceImmersion/signedDistanceCalculation/insideOutsidePropagation
    src/triSurfaceImmersion/signedDistanceCalculation/searchDistanceCalculator
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/levelSetDistCalc
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/signedDistanceCalculator
    src/triSurfaceImmersion/signedDistanceCalculation/signedDistanceCalculator/triSurfaceDistCalc
    src/triSurfaceImmersion/tetrahedralTools/AdaptiveTetCellRefinement
    src/triSurfaceImmersion/tetrahedralTools/IntersectionCriteria
    src/triSurfaceImmersion/tetrahedralTools/tetVolumeFractionCalculator
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellApproximation
    src/triSurfaceImmersion/volumeFractionCalculators/surfaceMeshCellIntersection
    src/triSurfaceImmersion/volumeFractionCalculators/volumeFractionCalculator
)
target_include_directories("triSurfaceImmersion" SYSTEM PUBLIC
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
# TODO: enable if OPenMP is re-enabled (TT)
#target_compile_options("triSurfaceImmersion" PRIVATE "-fopenmp")


add_library("pandora" SHARED
    src/pandora/curvatureFunctionObject/curvatureFunctionObject.cpp
    src/pandora/pandora/pandora.cpp
    src/pandora/pandoraCurvature/pandoraCurvature/pandoraCurvature.cpp
    src/pandora/pandoraCurvature/pandoraDivGradCurvature/pandoraDivGradCurvature.cpp
    src/pandora/pandoraCurvature/pandoraDivNormalCurvature/pandoraDivNormalCurvature.cpp
    src/pandora/pandoraCurvature/pandoraIsosurfaceCurvature/pandoraIsosurfaceCurvature.cpp
    src/pandora/pandoraCurvature/pandoraIsosurfaceCellCurvature/pandoraIsosurfaceCellCurvature.cpp
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
    src/pandora/pandoraCurvature/pandoraIsosurfaceCurvature
    src/pandora/pandoraCurvature/pandoraIsosurfaceCellCurvature
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
target_include_directories("pandora" SYSTEM PUBLIC
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
    ${of_src_path}/sampling/lnInclude
    ${of_src_path}/surfMesh/lnInclude
)


add_library("twoPhaseValidation" SHARED
    # NOTE: commented lines give errors due to missing header files (TT)
    src/twoPhaseValidation/divFree/divFreeFunctionObject/divFreeFunctionObject.C
    src/twoPhaseValidation/divFree/fieldModels/deformation2DFieldModel/deformation2DFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/deformationFieldModel/deformationFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/divFreeFieldModel/divFreeFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/harmonicFieldModel/harmonicFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/periodicFieldModel/periodicFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/rotationFieldModel/rotationFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/shear2DFieldModel/shear2DFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/shearFieldModel/shearFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/supermanFieldModel/supermanFieldModel.C
    src/twoPhaseValidation/divFree/fieldModels/translationFieldModel/translationFieldModel.C
    #src/twoPhaseValidation/error/fieldErrorBoundedness/fieldErrorsBoundedness.C
    src/twoPhaseValidation/error/fieldErrorBoundedness/volScalarFieldErrorBoundedness.C
    src/twoPhaseValidation/error/fieldErrorComposite/fieldErrorComposites.C
    src/twoPhaseValidation/error/fieldError/fieldErrors.C
    src/twoPhaseValidation/error/fieldErrorFunctionObject/fieldErrorFunctionObject.C
    #src/twoPhaseValidation/error/fieldErrorInterfaceSharpness2/fieldErrorsInterfaceSharpness2.C
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness2/volScalarFieldErrorInterfaceSharpness2.C
    #src/twoPhaseValidation/error/fieldErrorInterfaceSharpness3/fieldErrorsInterfaceSharpness3.C
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness3/volScalarFieldErrorInterfaceSharpness3.C
    #src/twoPhaseValidation/error/fieldErrorInterfaceSharpness/fieldErrorsInterfaceSharpness.C
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness/volScalarFieldErrorInterfaceSharpness.C
    src/twoPhaseValidation/error/fieldErrorL1/fieldErrorsL1.C
    src/twoPhaseValidation/error/fieldErrorL1normalized/fieldErrorsL1normalized.C
    src/twoPhaseValidation/error/fieldErrorVolume/fieldErrorsVolume.C
    src/twoPhaseValidation/hydrodynamics/hydrodynamicFunctionObject/hydrodynamicFunctionObject.C
    src/twoPhaseValidation/hydrodynamics/oscillatingDropletFunctionObject/oscillatingDropletFunctionObject.C
    src/twoPhaseValidation/hydrodynamics/risingBubbleFunctionObject/risingBubbleFunctionObject.C
    src/twoPhaseValidation/hydrodynamics/stationaryDropletFunctionObject/stationaryDropletFunctionObject.C
)
target_include_directories("twoPhaseValidation" PUBLIC
    src/twoPhaseValidation/divFree/divFreeFunctionObject
    src/twoPhaseValidation/divFree/fieldModels/deformation2DFieldModel
    src/twoPhaseValidation/divFree/fieldModels/deformationFieldModel
    src/twoPhaseValidation/divFree/fieldModels/divFreeFieldModel
    src/twoPhaseValidation/divFree/fieldModels/harmonicFieldModel
    src/twoPhaseValidation/divFree/fieldModels/periodicFieldModel
    src/twoPhaseValidation/divFree/fieldModels/rotationFieldModel
    src/twoPhaseValidation/divFree/fieldModels/shear2DFieldModel
    src/twoPhaseValidation/divFree/fieldModels/shearFieldModel
    src/twoPhaseValidation/divFree/fieldModels/supermanFieldModel
    src/twoPhaseValidation/divFree/fieldModels/translationFieldModel
    src/twoPhaseValidation/error/fieldError
    src/twoPhaseValidation/error/fieldErrorBoundedness
    src/twoPhaseValidation/error/fieldErrorComposite
    src/twoPhaseValidation/error/fieldErrorFunctionObject
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness2
    src/twoPhaseValidation/error/fieldErrorInterfaceSharpness3
    src/twoPhaseValidation/error/fieldErrorL1
    src/twoPhaseValidation/error/fieldErrorL1normalized
    src/twoPhaseValidation/error/fieldErrorVolume
    src/twoPhaseValidation/error/testing
    src/twoPhaseValidation/hydrodynamics/hydrodynamicFunctionObject
    src/twoPhaseValidation/hydrodynamics/oscillatingDropletFunctionObject
    src/twoPhaseValidation/hydrodynamics/risingBubbleFunctionObject
    src/twoPhaseValidation/hydrodynamics/stationaryDropletFunctionObject
    src/twoPhaseValidation/meshMagnitudeTemplate
)
# Using PUBLIC here removes the need to specify again the required include directories
# for the applications.
target_include_directories("twoPhaseValidation" SYSTEM PUBLIC
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
)
#target_link_libraries("twoPhaseValidation"
#    ${OF_DEPS}
#)


# Install the libraries
install(TARGETS
    "pandora"
    "triSurfaceImmersion"
    "twoPhaseValidation"
    LIBRARY
)


####################################################################################################
#   Argo applications
####################################################################################################

# Testing only
add_executable("geophaseUnitTest"
                applications/test/geophaseUnitTest/geophaseUnitTest.cpp)
            target_include_directories("geophaseUnitTest" SYSTEM PRIVATE
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

add_executable("pandoraTestCurvature"
    applications/test/pandoraTestCurvature/pandoraTestCurvature.cpp
)
target_include_directories("pandoraTestCurvature" SYSTEM PRIVATE
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

# Utilities: preprocessing
add_executable("surfaceInitVolumeFraction"
               applications/utilities/preProcessing/surfaceInitVolumeFraction/surfaceInitVolumeFraction.cpp)
target_link_libraries("surfaceInitVolumeFraction"
    "triSurfaceImmersion"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("surfaceInitSignedDistances"
               applications/utilities/preProcessing/surfaceInitSignedDistances/surfaceInitSignedDistances.cpp)
target_link_libraries("surfaceInitSignedDistances"
    "triSurfaceImmersion"
    ${OF_FINITEVOLUME}
    ${OF_MESHTOOLS}
    ${OF_OPENFOAM}
    dl
    m
)

# Utilities: mesh manipulation
add_executable("foamPerturbMesh"
               applications/utilities/mesh/manipulation/foamPerturbMesh/foamPerturbMesh.cpp)
target_include_directories("foamPerturbMesh" SYSTEM PRIVATE
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

# Utilities: miscellaneous
add_executable("computeLiquidVolume"
    applications/utilities/miscellaneous/computeLiquidVolume/computeLiquidVolume.cpp)
target_include_directories("computeLiquidVolume" SYSTEM PRIVATE
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

# Solvers
add_executable("interIsoPandoraFoam"
    applications/solvers/interIsoPandoraFoam/interIsoPandoraFoam.cpp)
target_include_directories("interIsoPandoraFoam" SYSTEM PRIVATE
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

# OpenFOAM solvers modified for high density ratios
add_executable("interIsoRhoFoam"
    applications/solvers/interIsoRhoFoam/interIsoRhoFoam.cpp)
target_include_directories("interIsoRhoFoam" SYSTEM PRIVATE
    ${of_src_path}/dynamicFvMesh/lnInclude
    ${of_src_path}/dynamicMesh/lnInclude
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/phaseSystemModels/twoPhaseInter/incompressibleInterPhaseTransportModel/lnInclude
    ${of_src_path}/phaseSystemModels/twoPhaseInter/VoFphaseIncompressibleTurbulenceModels/lnInclude
    ${of_src_path}/sampling/lnInclude
    ${of_src_path}/surfMesh/lnInclude
    ${of_src_path}/transportModels
    ${of_src_path}/transportModels/geometricVoF/lnInclude
    ${of_src_path}/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude
    ${of_src_path}/transportModels/incompressible/lnInclude
    ${of_src_path}/transportModels/interfaceProperties/lnInclude
    ${of_src_path}/transportModels/twoPhaseMixture/lnInclude
    ${of_src_path}/TurbulenceModels/incompressible/lnInclude
    ${of_src_path}/TurbulenceModels/phaseIncompressible/lnInclude
    ${of_src_path}/TurbulenceModels/turbulenceModels/lnInclude
    ${of_app_path}/solvers/multiphase/interIsoFoam
    ${of_app_path}/solvers/multiphase/VoF
    ${of_app_path}/solvers/multiphase/interFoam
)
target_link_libraries("interIsoRhoFoam"
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
    ${OF_VOFPHASETURBULENTTRANSPORTMODELS}
    ${OF_INCOMPRESSIBLEINTERPHASETRANSPORTMODELS}
    ${OF_OPENFOAM}
    dl
    m
)

add_executable("interRhoFoam"
    applications/solvers/interRhoFoam/interRhoFoam.cpp)
target_include_directories("interRhoFoam" SYSTEM PRIVATE
    ${of_src_path}/dynamicFvMesh/lnInclude
    ${of_src_path}/finiteVolume/lnInclude
    ${of_src_path}/meshTools/lnInclude
    ${of_src_path}/OpenFOAM/lnInclude
    ${of_src_path}/OSspecific/POSIX/lnInclude
    ${of_src_path}/phaseSystemModels/twoPhaseInter/incompressibleInterPhaseTransportModel/lnInclude
    ${of_src_path}/phaseSystemModels/twoPhaseInter/VoFphaseIncompressibleTurbulenceModels/lnInclude
    ${of_src_path}/sampling/lnInclude
    ${of_src_path}/transportModels
    ${of_src_path}/transportModels/geometricVoF/lnInclude
    ${of_src_path}/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude
    ${of_src_path}/transportModels/incompressible/lnInclude
    ${of_src_path}/transportModels/interfaceProperties/lnInclude
    ${of_src_path}/transportModels/twoPhaseMixture/lnInclude
    ${of_src_path}/TurbulenceModels/incompressible/lnInclude
    ${of_src_path}/TurbulenceModels/phaseIncompressible/lnInclude
    ${of_src_path}/TurbulenceModels/turbulenceModels/lnInclude
    ${of_app_path}/solvers/multiphase/VoF
    ${of_app_path}/solvers/multiphase/interFoam
)
target_link_libraries("interRhoFoam"
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
    ${OF_VOFPHASETURBULENTTRANSPORTMODELS}
    ${OF_INCOMPRESSIBLEINTERPHASETRANSPORTMODELS}
    ${OF_OPENFOAM}
    dl
    m
)


# Install executables
install(TARGETS
        "computeLiquidVolume"
        "foamPerturbMesh"
        "geophaseUnitTest"
        "interIsoPandoraFoam"
        "interIsoRhoFoam"
        "interRhoFoam"
        "pandoraTestCurvature"
        "surfaceInitVolumeFraction"
        "surfaceInitSignedDistances"
        RUNTIME)
