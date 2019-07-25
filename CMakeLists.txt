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

# Path to OpenFOAM (See section 5.3 of CMake book for Cache variables)
set(ofpath "$ENV{WM_PROJECT_DIR}" CACHE PATH "Path to OpenFOAM project folder.")
# TODO: check if not empty
if(NOT IS_DIRECTORY ${ofpath})
    message(FATAL_ERROR "No valid path to OpenFOAM set. Cannot build library.")
endif()

# Define values controlled by environment variables used by OpenFOAM
set(wm_label_size "$ENV{WM_LABEL_SIZE}" CACHE STRING "Size in bit to use as label type. Can be either 32 or 64.")

# Gather macro options required for OpenFOAM source files
add_compile_definitions(WM_LABEL_SIZE=${wm_label_size})
add_compile_definitions(WM_DP)
add_compile_definitions(NoRepository)
add_compile_definitions(OPENFOAM=1806)
add_compile_definitions(linux64)
add_compile_definitions(WM_ARCH_OPTION=64)

# Experimental
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

# Require Eigen 
#find_package(Eigen3 REQUIRED)
find_library(OF_FINITE_VOLUME finiteVolume PATHS "$ENV{WM_PROJECT_DIR}/platforms/linux64GccDPInt32Opt/lib")
find_library(OF_MESH_TOOLS meshTools PATHS "$ENV{WM_PROJECT_DIR}/platforms/linux64GccDPInt32Opt/lib")
find_library(OF_OPEN_FOAM OpenFOAM PATHS "$ENV{WM_PROJECT_DIR}/platforms/linux64GccDPInt32Opt/lib")

# Build the library
# TODO: CMake book recommends to leave out the type of library (STATIC SHARED MODULE).
#       Check which time is required / reasonable for this project (TT)
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
        ${ofpath}/src/OSspecific/POSIX/lnInclude
        ${ofpath}/src/OpenFOAM/lnInclude
        ${ofpath}/src/finiteVolume/lnInclude
        ${ofpath}/src/surfMesh/lnInclude
        ${ofpath}/src/meshTools/lnInclude
    )

# Executables

    # Tests
    add_executable(voFoamTestCellCellIntersectMeshes applications/test/transport/voFoamTestCellCellIntersectMeshes/voFoamTestCellCellIntersectMeshes.cpp) 
    target_include_directories(voFoamTestCellCellIntersectMeshes PUBLIC #TODO: determine if PUBLIC is required or if PRIVTE is possible
        src/geomMeshIntersection/geomMeshIntersection
        src/include
        ${ofpath}/src/finiteVolume/lnInclude
        ${ofpath}/src/meshTools/lnInclude
        ${ofpath}/src/OpenFOAM/lnInclude
        ${ofpath}/src/OSspecific/POSIX/lnInclude
    )
    # NOTE: the last three entries($Open_foam, dl, m) must be present
    target_link_libraries(voFoamTestCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITE_VOLUME}
        ${OF_MESH_TOOLS}
        ${OF_OPEN_FOAM}
        dl
        m
    )

    add_executable(voFoamTestSurfaceCellIntersectMeshes applications/test/transport/voFoamTestSurfaceCellIntersectMeshes/voFoamTestSurfaceCellIntersectMeshes.cpp)
    target_include_directories(voFoamTestSurfaceCellIntersectMeshes PUBLIC
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
        ${ofpath}/src/finiteVolume/lnInclude
        ${ofpath}/src/meshTools/lnInclude
        ${ofpath}/src/surfMesh/lnInclude
        ${ofpath}/src/OpenFOAM/lnInclude
        ${ofpath}/src/OSspecific/POSIX/lnInclude
    )
    target_link_libraries(voFoamTestSurfaceCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITE_VOLUME}
        ${OF_MESH_TOOLS}
        ${OF_OPEN_FOAM}
        dl
        m
    )

    # Preprocessing
    add_executable(voFoamCellCellIntersectMeshes applications/utilities/preProcessing/voFoamCellCellIntersectMeshes/voFoamCellCellIntersectMeshes.cpp)
    target_include_directories(voFoamCellCellIntersectMeshes PUBLIC
        src/include
        src/geomMeshIntersection/geomMeshIntersection
        ${ofpath}/src/finiteVolume/lnInclude
        ${ofpath}/src/meshTools/lnInclude
        ${ofpath}/src/OpenFOAM/lnInclude
        ${ofpath}/src/OSspecific/POSIX/lnInclude
    )
    target_link_libraries(voFoamCellCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITE_VOLUME}
        ${OF_MESH_TOOLS}
        ${OF_OPEN_FOAM}
        dl
        m
    )

    add_executable(voFoamSurfaceCellIntersectMeshes applications/utilities/preProcessing/voFoamSurfaceCellIntersectMeshes/voFoamSurfaceCellIntersectMeshes.cpp)
    target_include_directories(voFoamSurfaceCellIntersectMeshes PUBLIC
        src/geomMeshIntersection/geomSurfaceCellMeshIntersection
        src/geomMeshIntersection/geomTriSurfaceTools
        ${ofpath}/src/finiteVolume/lnInclude
        ${ofpath}/src/meshTools/lnInclude
        ${ofpath}/src/surfMesh/lnInclude
        ${ofpath}/src/OpenFOAM/lnInclude
        ${ofpath}/src/OSspecific/POSIX/lnInclude
    )
    target_link_libraries(voFoamSurfaceCellIntersectMeshes
        geometricalMeshIntersection
        ${OF_FINITE_VOLUME}
        ${OF_MESH_TOOLS}
        ${OF_OPEN_FOAM}
        dl
        m
    )

    ## Configure the INTERFACE library 
    #    add_library(geophase INTERFACE) 
    #    target_include_directories(geophase INTERFACE 
    #        ${EIGEN3_INCLUDE_DIR}
    #        include/ 
    #        include/core/
    #        include/io/
    #        include/io/vtk/
    #        include/testing/
    #        include/geometrical_algorithm/
    #        include/geometrical_model/
    #        include/interface_positioning/
    #    )
    # 
    ## Executables and tests
    #
    #    # Unit testing
    #    add_executable(geophaseUnitTest app/test/geophaseUnitTest.cpp) 
    #    target_link_libraries(geophaseUnitTest geophase gtest)
    #    add_test(NAME unitTests COMMAND geophaseUnitTest)
    #    set_tests_properties(unitTests PROPERTIES LABELS "UNIT")
    #
    #    # Perturbation 
    #    add_executable(geophaseTestPerturb app/test/geophaseTestPerturb.cpp) 
    #    target_link_libraries(geophaseTestPerturb geophase gtest)
    #    add_test(NAME perturbationTests COMMAND geophaseTestPerturb)
    #    set_tests_properties(perturbationTests PROPERTIES LABELS "UNIT")
    #
    #    # Inside tests 
    #    add_executable(geophaseTestInside app/test/geophaseTestInside.cpp) 
    #    target_include_directories(geophaseTestInside PUBLIC /usr/include/png++)
    #    target_link_directories(geophaseTestInside PUBLIC /usr/lib)
    #    target_link_libraries(geophaseTestInside geophase gtest png)
    #    add_test(NAME insideOutsideTests COMMAND geophaseTestInside)
    #    set_tests_properties(insideOutsideTests PROPERTIES LABELS "UNIT")
    #
    #    # Area and volume
    #    add_executable(geophaseTestAreaVolume app/test/geophaseTestAreaVolume.cpp) 
    #    target_link_libraries(geophaseTestAreaVolume geophase gtest)
    #    add_test(NAME areaAndVolumeTests COMMAND geophaseTestAreaVolume)
    #    set_tests_properties(areaAndVolumeTests PROPERTIES LABELS "UNIT")
    #
    #    # Extreme points 
    #    add_executable(geophaseTestExtremePoints app/test/geophaseTestExtremePoints.cpp) 
    #    target_link_libraries(geophaseTestExtremePoints geophase gtest)
    #    add_test(NAME extremePointsTests COMMAND geophaseTestExtremePoints)
    #    set_tests_properties(extremePointsTests PROPERTIES LABELS "UNIT")
    #
    #    # Triangulation
    #    add_executable(geophaseTestTriangulation app/test/geophaseTestTriangulation.cpp) 
    #    target_link_libraries(geophaseTestTriangulation geophase gtest)
    #    add_test(NAME triangulationTests COMMAND geophaseTestTriangulation)
    #    set_tests_properties(triangulationTests PROPERTIES LABELS "UNIT")
    #
    #    # Polygon orientation 
    #    add_executable(geophaseTestOrient app/test/geophaseTestOrient.cpp) 
    #    target_link_libraries(geophaseTestOrient geophase gtest)
    #    add_test(NAME orientationTests COMMAND geophaseTestOrient)
    #    set_tests_properties(orientationTests PROPERTIES LABELS "UNIT")
    #
    #    # Intersection 
    #    add_executable(geophaseTestIntersect app/test/geophaseTestIntersect.cpp) 
    #    target_link_libraries(geophaseTestIntersect geophase gtest)
    #    add_test(NAME intersectionTests COMMAND geophaseTestIntersect)
    #    set_tests_properties(intersectionTests PROPERTIES LABELS "UNIT")
    #
    #    # Interface positioning
    #    add_executable(geophaseTestInterfacePositioning app/test/geophaseTestInterfacePositioning.cpp) 
    #    target_link_libraries(geophaseTestInterfacePositioning geophase gtest)
    #    add_test(NAME interfacePositioningUnitTests COMMAND geophaseTestInterfacePositioning --gtest_filter=UNIT*)
    #    set_tests_properties(interfacePositioningUnitTests PROPERTIES LABELS "UNIT")
    #    add_test(NAME interfacePositioningConvergenceTests COMMAND geophaseTestInterfacePositioning --gtest_filter=CONVERGENCE*)
    #    set_tests_properties(interfacePositioningConvergenceTests PROPERTIES LABELS "CONVERGENCE")
