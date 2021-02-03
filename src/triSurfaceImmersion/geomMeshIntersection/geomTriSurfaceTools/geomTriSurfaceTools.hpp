/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Description
    Global functions used for testing the Surface Cell Mesh Intersection (SMCI)
    algorithm. 

    Some functions assume that the triangulated surface is star-shaped. 

    A star-shaped triangulated surface is a surface whose volume can be
    decomposed into tetrahedra with using its centroid without resulting in
    tetrahedra that have negative volumes.  

SourceFiles
    geomTriSurfaceTools.C

Author 
    Tomislav Maric, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com

\*---------------------------------------------------------------------------*/

#ifndef geomTriSurfaceTools_H
#define geomTriSurfaceTools_H

#include "triSurface.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace GeometricalTransport { 

// Orient the triangles of a star-shaped triangulated surface such that they
// are all pointing inwards. 
void orientNormalsInward(triSurface& tri);

// Compute the volume of a star-shaped triangulated surface.
// A two-dimensional triangular surface is a surface-stripe.  
scalar starSurfaceVolume
(
    const triSurface& tri, 
    const Vector<label>& solutionVector // Volume mesh solution vector.
); 

// Displace any triangulated surface without erasing the geometrical +
// topological data.
void displaceSurface(triSurface& tri, const vector& displacement); 

// Displace any triangulated surface randomly within the box 'bbox' such that
// bbox contains the axis aligned box of the surface 'tri'. Mesh solution
// vector is used to determine the inactive displacement vector components for
// a 2D surface mesh. Returns the random displacement.
vector placeSurfaceRandomlyInBox 
(
    triSurface& tri, 
    const boundBox& bbox, 
    const Vector<label> solutionVector // Volume mesh solution vector.
); 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
