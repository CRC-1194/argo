/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    voFoamTestGeophaseMake

Description
    Test interface to geophase.

    Unit test application for building geophase models from OpenFOAM mesh
elements.

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

// Make sure the STL assertions are enabled. They are disabled in Release mode.
#undef NDEBUG
#include <cassert>
#include <limits>

#include "Geophase.hpp"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    const double EPSILON = std::numeric_limits<double>::epsilon();

    // DO NOT CHANGE THESE POINTS, COMPARISON VALUES BELOW ARE BASED ON THEM.
    const point p0(0, 0, 0);
    const point p1(1, 0, 0);
    const point p2(0, 1, 0);
    const point p3(0, 0, 1);

    // Test equivalence of points.
    assert(equal_by_tolerance(p0, p0, EPSILON));

    // Compute an OpenFOAM triangle from OpenFOAM points.
    foamVectorPolygon foamTriangle{p0, p1, p2};

    // Compute the area normal vector of the OpenFOAM triangle.
    auto foamTriNormal = area_normal_vector(foamTriangle);
    assert(foamTriNormal == point(0, 0, 0.5));
    assert(std::abs(foamTriNormal[0]) < EPSILON);
    assert(std::abs(foamTriNormal[1]) < EPSILON);
    assert(std::abs(foamTriNormal[2] - 0.5) < EPSILON);

    // Triangle orientation test.
    auto orientedTriangle = orient(foamTriangle, point(0, 0, 1));
    auto orientedNormal = area_normal_vector(orientedTriangle);
    auto foamNormal = area_normal_vector(foamTriangle);
    assert(equal_by_tolerance(foamNormal, orientedNormal, EPSILON));

    // Write the tetrahedron out as VTK.
    auto unitTet = geophase::make_unit_tetrahedron();
    write_vtk(unitTet, "unitTet.vtk");
    std::cout << "Unit tetrahedron volume by centroid triangulation (should be "
                 "1/6.) = "
              << volume_by_vol_tri(unitTet) << std::endl;

    auto unitCentroid = centroid(unitTet);
    std::cout << "Unit tetrahedron centroid (should be (0.25,0.25,0.25)) = "
              << unitCentroid << std::endl;
    assert((std::abs(unitCentroid[0] - 0.25) < EPSILON));
    assert((std::abs(unitCentroid[1] - 0.25) < EPSILON));
    assert((std::abs(unitCentroid[2] - 0.25) < EPSILON));

    // Make a 3D tetrahedron from OpenFOAM points.
    auto foamTetrahedron =
        make_tetrahedron<foamVectorPolyhedron>(p0, p1, p2, p3);
    write_vtk(foamTetrahedron, "foamTetrahedron.vtk");

    // OpenFOAM tetrahedron volume using surface the surface integral.
    auto tetraVolumeBySurfTri = volume_by_surf_tri(foamTetrahedron);
    std::cout << "FOAM unit tetrahedron volume by surface triangulation "
                 "(should be 1/6.) = "
              << tetraVolumeBySurfTri << std::endl;
    assert((std::abs(tetraVolumeBySurfTri - 1.0 / 6.) < EPSILON));

    // Volume of the OpenFOAM tetrahedron using centroid triangulation.
    auto tetraVolumeByVolTri = volume_by_vol_tri(foamTetrahedron);
    std::cout << "FOAM unit tetrahedron volume by surface triangulation "
                 "(should be 1/6.) = "
              << tetraVolumeBySurfTri << std::endl;
    assert((std::abs(tetraVolumeByVolTri - 1.0 / 6.) < EPSILON));

    // OpenFOAM tetrahedron centroid
    auto foamTetrahedronCentroid = centroid(foamTetrahedron);
    std::cout << "Unit tetrahedron centroid = (should be (0.25,0.25,0.25))"
              << unitCentroid << std::endl;
    assert((std::abs(foamTetrahedronCentroid[0] - 0.25) < EPSILON));
    assert((std::abs(foamTetrahedronCentroid[1] - 0.25) < EPSILON));
    assert((std::abs(foamTetrahedronCentroid[2] - 0.25) < EPSILON));
    // Equal to the unit tetrahedron centroid.
    assert((std::abs(foamTetrahedronCentroid[0] - unitCentroid[0]) < EPSILON));
    assert((std::abs(foamTetrahedronCentroid[1] - unitCentroid[1]) < EPSILON));
    assert((std::abs(foamTetrahedronCentroid[2] - unitCentroid[2]) < EPSILON));

    // Make an OpenFOAM halfspace at the centroid of the unit tetrahedron.
    auto centroidHspace =
        foamHalfspace(centroid(foamTetrahedron), vector(0, 0, 1));

    // Test signed distance to halfspace.
    Info << "Centroid hspace position = " << centroidHspace.position() << nl
         << "Centroid hspace direction = " << centroidHspace.direction()
         << endl;
    assert((std::abs(signed_distance(centroidHspace, point(0, 0, 1)) - 0.75) <
        EPSILON));

    // Halve the OpenFOAM tetrahedron with the OpenFOAM halfspace.
    auto halvingHspace = foamHalfspace(point(0, 0, 0), vector(1, -1, 0));
    auto foamTetraIntersection =
        intersect_tolerance<foamPolyhedronIntersection>(
            foamTetrahedron, halvingHspace, EPSILON);
    const auto& halvedFoamTetrahedron = foamTetraIntersection.polyhedron();
    write_vtk(halvedFoamTetrahedron, "halvedFoamTetrahedron.vtk");
    double Vhalved = volume_by_surf_tri(halvedFoamTetrahedron);
    double Vfoam = volume_by_surf_tri(foamTetrahedron);
    Info << "Vhalved = " << Vhalved << nl << "Vfoam = " << Vfoam << endl;
    assert((std::abs(Vhalved - 0.5 * Vfoam) < EPSILON));

    Info << "SUCCESS" << endl;
    return 0;
}


// ************************************************************************* //
