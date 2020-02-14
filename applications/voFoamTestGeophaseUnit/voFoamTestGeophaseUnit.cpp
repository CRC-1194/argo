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

    Unit test application for building geophase models from OpenFOAM mesh elements.

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

#include "GeophaseMake.hpp"
#include "Make.hpp"
#include "Area.hpp"
#include "Volume.hpp"
#include "WriteVtkPolyData.hpp"
#include "ReadVtkPolyData.hpp"

#include "messageStream.H"
#include "vector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using geophase::vectorPolygon;
using geophase::vectorPolyhedron;

using geophase::foamVectorPolygon;
using geophase::foamVectorPolyhedron;

using geophase::make_tetrahedron;

int main(int argc, char *argv[])
{
    const double EPSILON = std::numeric_limits<double>::epsilon();

    // DO NOT CHANGE THESE POINTS, COMPARISON VALUES BELOW ARE BASED ON THEM. 
    const point p0 (0,0,0); 
    const point p1 (1,0,0);
    const point p2 (0,1,0);
    const point p3 (0,0,1);

    // Compute an OpenFOAM triangle from OpenFOAM points.
    foamVectorPolygon foamTriangle {p0,p1,p2};  

    // Compute the area normal vector of the OpenFOAM triangle.
    auto foamTriNormal = area_normal_vector(foamTriangle);  
    assert(foamTriNormal == point(0,0,0.5));

    assert(std::abs(foamTriNormal[0]) < EPSILON);
    assert(std::abs(foamTriNormal[1]) < EPSILON);
    assert(std::abs(foamTriNormal[2] - 0.5) < EPSILON);

    // TODO: Remove
    auto unitTet = geophase::make_unit_tetrahedron(); 
    write_vtk(unitTet, "unitTet.vtk");
    std::cout << volume_by_vol_tri(unitTet) << std::endl;
    // TODO: Remove

    // Make a 3D tetrahedron from OpenFOAM points.
    auto tetrahedron = make_tetrahedron<foamVectorPolyhedron>(p0,p1,p2,p3);
    write_vtk(tetrahedron, "tetrahedron.vtk");

    // Compute the volume of the OpenFOAM tetrahedron using the surface
    // integral.
    auto tetraVolumeBySurfTri = volume_by_surf_tri(tetrahedron);
    Info << tetraVolumeBySurfTri << endl;
    assert ((std::abs(tetraVolumeBySurfTri - 1.0 / 6.) < EPSILON));

    // Compute the volume of the OpenFOAM tetrahedron using the centroid
    // triangulation.
    auto tetraVolumeByVolTri = volume_by_vol_tri(tetrahedron);
    Info << tetraVolumeByVolTri << endl;
    assert ((std::abs(tetraVolumeByVolTri - 1.0 / 6.) < EPSILON));

    // Make an OpenFOAM halfspace. 

    // Intersect the OpenFOAM tetrahedron with the OpenFOAM halfspace. 

    Info<< "SUCCESS" << endl;
}


// ************************************************************************* //
