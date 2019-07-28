include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(FetchContent)

cmake_minimum_required(VERSION 3.13)

project(geo-vof-initialization VERSION 0.8
                DESCRIPTION "Initialization of volume fraction fields by surface/volume mesh intersection or volume/volume mesh intersection."
                LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build dependencies
FetchContent_Declare(geometricLibrary_source
    GIT_REPOSITORY          git@git.rwth-aachen.de:leia/geophase.git
    GIT_TAG                 9360ddb7fa9951978f6236fc818ce979290c8e01
    GIT_SHALLOW             YES
    GIT_PROGRESS            YES
    USES_TERMINAL_DOWNLOAD  YES
    USES_TERMINAL_UPDATE    YES
)

FetchContent_GetProperties(geometricLibrary_source)
if(NOT geometricLibrary_source_POPULATED)
    FetchContent_Populate(geometricLibrary_source)
    add_subdirectory(${geometriclibrary_source_SOURCE_DIR})
endif()

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
    message(FATAL_ERROR "Path ${of_wm_project_dir}/platforms/${of_wm_options}${of_wm_compile_option} does not exist.")
endif()

set(of_lib_path "${of_wm_project_dir}/platforms/${of_wm_options}/lib" CACHE PATH "Path to compiled OpenFOAM libraries.")
set(of_src_path "${of_wm_project_dir}/src" CACHE PATH "Path to OpenFOAM/src folder")

message(STATUS "OpenFOAM lib path: ${of_lib_path}")
message(STATUS "OpenFOAM src path: ${of_src_path}")

# TODO:
# Check that those variables are set. Optional: check for reasonable values (TT)

# Compile definitions required for OpenFOAM
add_compile_definitions(
    WM_LABEL_SIZE=${of_wm_label_size}
    WM_${of_wm_precision_option}
    WM_ARCH_OPTION=${of_wm_arch_option}
    ${of_wm_arch}
    OPENFOAM=1806 # Change when moving to other OpenFOAM version
    NoRepository
)

# Required to make linking to OpenFOAM libraries work
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed")

if(MSVC)
    add_compile_options(/W4 /WX)
else(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -Wall -Wextra -DDEBUG -pedantic -D_USE_MATH_DEFINES -DBOOST_MATH_INSTRUMENT")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra -pedantic -Wno-deprecated -D_USE_MATH_DEFINES")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb3 -O3 -Wall -Wextra -pedantic -DDEBUG -D_USE_MATH_DEFINES")
endif(MSVC)

# Required OpenFOAM libraries
find_library(OF_FINITEVOLUME finiteVolume PATHS ${of_lib_path})
find_library(OF_MESHTOOLS meshTools PATHS ${of_lib_path})
find_library(OF_OPENFOAM OpenFOAM PATHS ${of_lib_path})

# Build the library
# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which type is required / reasonable for this project (TT)
    add_library(geometricalMeshIntersection SHARED
        src/geomMeshIntersection/geomMeshIntersection/geomMeshIntersection.cpp
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection/geomSurfaceCellMeshIntersection.cpp
        src/geomMeshIntersection/geomTriSurfaceTools/geomTriSurfaceTools.cpp
    )
    target_include_directories(geometricalMeshIntersection PUBLIC
        src/geomMeshIntersection/geomMeshIntersection
        src/geomMeshIntersection/geomReconstructError
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
        src/include
    )
    # Using PUBLIC here removes the need to specify again the required include directories
    # for the applications.
    target_include_directories(geometricalMeshIntersection SYSTEM PUBLIC
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/surfMesh/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
    )
    target_link_libraries(geometricalMeshIntersection
        geophase
    )

# Executables

    # Tests
    add_executable(voFoamTestCellCellIntersectMeshes applications/test/transport/voFoamTestCellCellIntersectMeshes/voFoamTestCellCellIntersectMeshes.cpp) 
    target_link_libraries(voFoamTestCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )

    add_executable(voFoamTestSurfaceCellIntersectMeshes applications/test/transport/voFoamTestSurfaceCellIntersectMeshes/voFoamTestSurfaceCellIntersectMeshes.cpp)
    target_link_libraries(voFoamTestSurfaceCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )

    # Preprocessing
    add_executable(voFoamCellCellIntersectMeshes applications/utilities/preProcessing/voFoamCellCellIntersectMeshes/voFoamCellCellIntersectMeshes.cpp)
    target_link_libraries(voFoamCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )

    add_executable(voFoamSurfaceCellIntersectMeshes applications/utilities/preProcessing/voFoamSurfaceCellIntersectMeshes/voFoamSurfaceCellIntersectMeshes.cpp)
    target_link_libraries(voFoamSurfaceCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )
