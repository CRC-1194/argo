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

#include "geomTriSurfaceTools.hpp"
#include "boundBox.H"
#include <random>
//#include "Build.H"
//#include "WriteVtkPolydata.H"

namespace Foam { namespace GeometricalTransport { 

void orientNormalsInward(triSurface& tri)
{
    const pointField& triPoints = tri.localPoints(); 
    pointField& triNormals = const_cast<pointField&>(tri.faceNormals()); 

    vector surfaceMeshCentroid = sum(triPoints) / tri.nPoints();  
    const auto& triCentres = tri.faceCentres(); 

    forAll(triNormals, tI)
    {
        if (((triCentres[tI] - surfaceMeshCentroid) & triNormals[tI]) > 0) 
        {
            triNormals[tI] *= -1; 
            tri[tI].flip(); 
        }
    }
}

scalar starSurfaceVolume(const triSurface& tri, const Vector<label>& solutionVector)
{
    scalar Vs = 0; 

    const pointField& triPoints = tri.points(); 
    const pointField& triLocalPoints = tri.localPoints(); 

    const vector triCenter = sum(tri.points()) / tri.nPoints(); 

    // 3D volume: 
    // Compute the surface mesh volume using the surface mesh centroid. 
    forAll(tri, tI)
    {
        const point& p0 = triPoints[tri[tI][0]];
        const point& p1 = triPoints[tri[tI][1]];
        const point& p2 = triPoints[tri[tI][2]];

        // Assuming outward pointing normals for distance calculation!
        scalar deltaV = 1.0 / 6.0 * (p0 - triCenter) & 
            ((p1 - triCenter) ^ (p2 - triCenter)); 

        if (deltaV < 0) 
        {
            /* ADDGEOM
            write_vtk_polydata(build<arrayTriangle>(p0,p1,p2), "bad-triangle.vtk");
            Perr << "Negative mesh volume contribution!\n" 
                << "This volume computation assumes the surface mesh is star-shaped:\n" 
                << "the volume it covers can be triangulated with its centroid." 
                << "deltaV = " << deltaV << endl
                << "triCenter = " << triCenter << endl
                << "p0 = " << p0 << endl
                << "p1 = " << p1 << endl
                << "p2 = " << p2 << endl
                << Foam::abort(FatalError);
            */
        }
        else 
            Vs += deltaV;
    }

    label dimension = 0; 
    forAll(solutionVector, I) 
        dimension += solutionVector[I]; 

    if (dimension < 3)
    {
        // Pseudo 2D: update the total volume for the edges that have only a
        // single edge owner.
        
        // Find the hight of the surface mesh using the bounding box height
        // and mesh solution vector. 
        const boundBox toolBox = boundBox(tri.points());  
        const vector toolBoxMin = toolBox.min();  
        const vector toolBoxMax = toolBox.max(); 

        vector triHeightVector = 0.5 * (toolBoxMax - toolBoxMin);  
        forAll(solutionVector, sI)
            if (solutionVector[sI] == 1)
                triHeightVector[sI] = 0; 

        // For each edge. 
        const edgeList& edges = tri.edges(); 
        forAll(edges, edgeI)
        {
            if (! tri.isInternalEdge(edgeI))
            {
                const vector& ePoint0 = triLocalPoints[edges[edgeI][0]]; 
                const vector& ePoint1 = triLocalPoints[edges[edgeI][1]]; 

                // Triangulate the boundary edge using the centroid. 
                
                // Some volumes will be positive, and some negative, as this
                // loop will iterate over both "top" and "bottom" edges with
                // respect to the height vector. That is why mag is used to
                // compute the absolute value. 
                Vs +=mag
                (
                    (1.0 / 6.0) * 
                    ((ePoint0 - triCenter) ^ (ePoint1 - triCenter))
                    & triHeightVector
                ); 
            }
        }
    }

    return Vs;
}

void displaceSurface(triSurface& tri, const vector& displacement)
{
    pointField& triPoints = const_cast<pointField&>(tri.points());
    pointField& triLocalPoints = const_cast<pointField&>(tri.localPoints());

    triPoints += displacement; 
    // Displace also local points used for surface addressing.
    triLocalPoints += displacement; 
}

vector placeSurfaceRandomlyInBox
(
    triSurface& tri, 
    const boundBox& bbox, 
    const Vector<label> solutionVector
)
{

    // Use a copy of the AABB box that has to be shrunk in order to displace the triSurface
    // in such a way that the bbox contains AABB(tri). 
    boundBox  aabb(bbox); 

    const boundBox triBox = boundBox(tri.points()); 
    vector triBoxDelta = 0.5 * (triBox.max() - triBox.min()); 
    // Null the triBoxDelta component that is not used in pseudo2D. 
    forAll(triBoxDelta, I)
        if (solutionVector[I] == -1)
            triBoxDelta[I] *= 0;

    // Shrink the base box, so that the random position of the tri mesh
    // centroid is placed in a way that the base mesh bounding box contains
    // the tri mesh bounding box. 
    point& baseMin = aabb.min(); 
    point& baseMax = aabb.max(); 
    baseMin += triBoxDelta; 
    baseMax -= triBoxDelta;  

    // Initialize the random number generator and distribution. 
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Compute the random displacement with respect to the base mesh centroid.
    vector randomPosition = baseMin;  
    forAll(randomPosition, cmptI)
        randomPosition[cmptI] += dis(gen) * (baseMax[cmptI] - baseMin[cmptI]);

    // Place the tri surface such that the centroid of the tri surface box overlaps
    // with the randomly generated position within aabb. 
    vector randomDisplacement = randomPosition - (0.5 * (triBox.min() + triBox.max()));

    // Zero the displacement component in a dimension that is not used in 2D. 
    forAll(randomDisplacement, I)
        if (solutionVector[I] == -1)
            randomDisplacement[I] *= 0;

    if (! aabb.contains(randomPosition))
    {
        FatalErrorIn("voFoamTestCellCellIntersectMeshes")                                                 
            << "Random position outside of the base mesh bounding box." 
            << Foam::abort(FatalError);
    }

    // Displace the surface using the random displacement.
    displaceSurface(tri, randomDisplacement);

    return randomDisplacement;
}

}} // End namespace Foam::GeometricalTransport


// ************************************************************************* //
