# For now only copied from geophase. READ Cmake book before adopting (TT)
# From the WMake build process:
# - generate a shared library (geometricalMeshIntersection)
# - build several executables:
#       * voFoamTestCellCellIntersectMeshes
#       * voFoamTestSurfaceCellIntersectMeshes
#       * voFoamCellCellIntersectMeshes
#       * voFoamSurfaceCellIntersectMeshes

cmake_minimum_required(VERSION 3.13)

project(geo-vof-initialization VERSION 0.8
                DESCRIPTION "Initialization of volume fraction fields by surface/volume mesh intersection or volume/volume mesh intersection."
                LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# TODO: read in the relevant environment variables used by the OpenFOAM build process.
# Store them by using CACHE variables.
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
# Check that those variables are set. Optional: check for reasonable values

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

# TODO: Re-enable GoogleTest once it is required (TT)
# GoogleTest integration by Craig Scott 
# https://crascit.com/2015/07/25/cmake-gtest/
#
#    # Download and unpack googletest at configure time
#    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
#    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
#        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
#    execute_process(COMMAND "${CMAKE_COMMAND}" --build .
#        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
#
#    # Prevent GoogleTest from overriding our compiler/linker options
#    # when building with Visual Studio
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#
#    # Add googletest directly to our build. This adds
#    # the following targets: gtest, gtest_main, gmock
#    # and gmock_main
#    add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
#                     "${CMAKE_BINARY_DIR}/googletest-build")
#
#    # The gtest/gmock targets carry header search path
#    # dependencies automatically when using CMake 2.8.11 or
#    # later. Otherwise we have to add them here ourselves.
#    if(CMAKE_VERSION VERSION_LESS 2.8.11)
#        include_directories("${gtest_SOURCE_DIR}/include"
#                            "${gmock_SOURCE_DIR}/include")
#    endif()

# Required OepnFOAM libraries
find_library(OF_FINITEVOLUME finiteVolume PATHS ${of_lib_path})
find_library(OF_MESHTOOLS meshTools PATHS ${of_lib_path})
find_library(OF_OPENFOAM OpenFOAM PATHS ${of_lib_path})

# Build the library
# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which type is required / reasonable for this project (TT)
    add_library(geometricalMeshIntersection SHARED
        src/geomMeshIntersection/geomMeshIntersection/geomMeshIntersection.C
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection/geomSurfaceCellMeshIntersection.C
        src/geomMeshIntersection/geomTriSurfaceTools/geomTriSurfaceTools.C
    )
    target_include_directories(geometricalMeshIntersection PUBLIC # TODO: check if PRIVATE is possible as flag
        src/geomMeshIntersection/geomMeshIntersection
        src/geomMeshIntersection/geomReconstructError
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
    )
    target_include_directories(geometricalMeshIntersection PUBLIC
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/surfMesh/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
    )

# Executables

    # Tests
    add_executable(voFoamTestCellCellIntersectMeshes applications/test/transport/voFoamTestCellCellIntersectMeshes/voFoamTestCellCellIntersectMeshes.cpp) 
    target_include_directories(voFoamTestCellCellIntersectMeshes PUBLIC #TODO: determine if PUBLIC is required or if PRIVTE is possible
        src/geomMeshIntersection/geomMeshIntersection
        src/include
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
    )
    # NOTE: the last three entries($Open_foam, dl, m) must be present
    target_link_libraries(voFoamTestCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )

    add_executable(voFoamTestSurfaceCellIntersectMeshes applications/test/transport/voFoamTestSurfaceCellIntersectMeshes/voFoamTestSurfaceCellIntersectMeshes.cpp)
    target_include_directories(voFoamTestSurfaceCellIntersectMeshes PUBLIC
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/surfMesh/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
    )
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
    target_include_directories(voFoamCellCellIntersectMeshes PUBLIC
        src/include
        src/geomMeshIntersection/geomMeshIntersection
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
    )
    target_link_libraries(voFoamCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )

    add_executable(voFoamSurfaceCellIntersectMeshes applications/utilities/preProcessing/voFoamSurfaceCellIntersectMeshes/voFoamSurfaceCellIntersectMeshes.cpp)
    target_include_directories(voFoamSurfaceCellIntersectMeshes PUBLIC
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
        ${of_src_path}/finiteVolume/lnInclude
        ${of_src_path}/meshTools/lnInclude
        ${of_src_path}/surfMesh/lnInclude
        ${of_src_path}/OpenFOAM/lnInclude
        ${of_src_path}/OSspecific/POSIX/lnInclude
    )
    target_link_libraries(voFoamSurfaceCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITEVOLUME}
        ${OF_MESHTOOLS}
        ${OF_OPENFOAM}
        dl
        m
    )
