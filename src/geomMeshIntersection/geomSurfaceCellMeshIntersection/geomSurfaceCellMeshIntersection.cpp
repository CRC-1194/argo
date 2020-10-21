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

\*---------------------------------------------------------------------------*/

// OpenFOAM 
#include "geomSurfaceCellMeshIntersection.hpp"
#include "calculatedFvPatchField.H"
#include "fvcAverage.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "fvc.H"

// Distance encoding
#include <bitset>
#include <cmath>

// geophase 
#include "Geophase.hpp"

// Testing
#include <cassert>
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam { namespace GeometricalTransport {

    defineTypeNameAndDebug(geomSurfaceCellMeshIntersection,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from a mesh (initialize field data)
geomSurfaceCellMeshIntersection::geomSurfaceCellMeshIntersection
(
    const fvMesh& mesh,
    const triSurface& triSurf,
    const scalar sqrDistFactor,
    const IOobject::writeOption& writeOption // Allows output for testing purposes.
)
:
    mesh_(mesh), 
    runTime_(mesh.time()),  
    cellSqrSearchDist_
    (
        IOobject
        (
            "cellSqrSearchDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            writeOption     
        ),
        fvc::average(pow(mesh.deltaCoeffs(), -2))
    ),
    cellSignedDist_
    (
        IOobject
        (
            "cellSignedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            writeOption 
        ),
        mesh,
        dimensionedScalar("cellSignedDist", dimLength,0),
        "zeroGradient"
    ),
    cellSignedDist0_("cellSignedDist0", cellSignedDist_), 
    cellNearestTriangle_(), 
    cellsToPointsInterp_(mesh),
    pMesh_(mesh_),
    pointSqrSearchDist_ 
    (
        IOobject
        (
            "pointSqrSearchDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            writeOption            
        ),
        pMesh_,
        dimensionedScalar("pointSqrSearchDist", dimLength,0),
        "zeroGradient"
    ),
    pointSignedDist_ 
    (
        IOobject
        (
            "pointSignedDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            writeOption            
        ),
        pMesh_,
        dimensionedScalar("pointSignedDist", dimLength,0),
        "zeroGradient"
    ),
    pointNearestTriangle_(),
    triSurf_(triSurf),
    triSurfSearch_(triSurf),
    intersectedCellLabels_(),
    sqrDistFactor_(max(2.0, sqrDistFactor)), 
    nTrianglesPerCell_(0.)
{
    calcSqrSearchDists(); 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void geomSurfaceCellMeshIntersection::calcSqrSearchDists()
{
    cellSqrSearchDist_ = fvc::average(pow(mesh_.deltaCoeffs(), -2));
    cellsToPointsInterp_.interpolate(cellSqrSearchDist_, pointSqrSearchDist_);
}

void geomSurfaceCellMeshIntersection::findIntersectedCells()
{
    const auto& meshCellPoints = mesh_.cellPoints();
    const auto& meshCellEdges = mesh_.cellEdges();

    intersectedCellLabels_.resize(0);

    forAll(cellSignedDist_, cellI)
    {
        const auto& cellDist = cellSignedDist_[cellI];  
        const auto& cellPoints = meshCellPoints[cellI]; 

        forAll(cellPoints, pointI)
        {
            if ((pointSignedDist_[cellPoints[pointI]] * cellDist) < 0)
            {
                intersectedCellLabels_.append(cellI); 
                break;
            }
        }
    }
    assert((intersectedCellLabels_.size() < mesh.nCells())); 
}

void geomSurfaceCellMeshIntersection::calcSignedDist()
{

    // Zero the signed distance.
    cellSignedDist_ = dimensionedScalar("cellSignedDist", dimLength, 0);

    // Use the octree and the square search distance multiplied by a distance factor 
    // to build the surface mesh / volume mesh proximity information list. 
    triSurfSearch_.findNearest(
        mesh_.C(),
        (sqrDistFactor_ * sqrDistFactor_) * cellSqrSearchDist_,
        cellNearestTriangle_
    );

    // Compute the signed distance in each cell as the distance between the triangle
    // nearest to the cell center and the cell center.
    const volVectorField& C = mesh_.C();  
    const pointField& triPoints = triSurf_.points();  
    const vectorField& triNormals = triSurf_.faceNormals(); 
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
        {
            cellSignedDist_[cellI] = 
                (C[cellI] - triPoints[triSurf_[cellHit.index()][0]]) & 
                triNormals[cellHit.index()];
        }
    }

    // Save the signed distance field given by the octree.
    cellSignedDist0_ = cellSignedDist_; 

    // Propagate the sign information into the bulk by solving a Laplace
    // equation for a single iteration for the signed distance field. 
    fvScalarMatrix distEqn
    (
        -fvm::laplacian(cellSignedDist_)
    );
    distEqn.solve(); 

    // Correct signed distances in the narrow band using 
    // octree signed distances.
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
            cellSignedDist_[cellI] = cellSignedDist0_[cellI]; 
    }

    // Once the cell-centered signed distance is computed, compute the point
    // distances only for determining intersected tetrahedra in the volume
    // calculation step. TM.
    // Use the octree and the square search distance multiplied by a distance factor 
    // to build the surface mesh / mesh points proximity information list. 
    triSurfSearch_.findNearest(
        mesh_.points(),
        (sqrDistFactor_ * sqrDistFactor_) * pointSqrSearchDist_,
        pointNearestTriangle_ 
    );

    cellsToPointsInterp_.interpolate(cellSignedDist_, pointSignedDist_);  
    const pointField& meshPoints = mesh_.points();  
    forAll(pointNearestTriangle_, pointI)
    {
        const pointIndexHit& pointHit = pointNearestTriangle_[pointI];

        if (pointHit.hit()) 
        {
            pointSignedDist_[pointI] = 
                (meshPoints[pointI] - triPoints[triSurf_[pointHit.index()][0]]) & 
                triNormals[pointHit.index()];
        }
    }

    findIntersectedCells(); 

    // Correct the volume fractions geometrically in the intersected cells. 
    forAll(intersectedCellLabels_, cellL)
    {
        auto cellK = intersectedCellLabels_[cellL];
        cellSignedDist_[cellK] = cellSignedDist0_[cellK];
    }
}

void geomSurfaceCellMeshIntersection::calcVolFraction(volScalarField& alpha)
{
    // Compute the signed distances using the octree in triSurfaceSearch. 
    calcSignedDist(); 
    

    const auto& meshCells = mesh_.cells(); 
    const auto& cellCenters = mesh_.C();  
    const auto& cellVolumes = mesh_.V(); 

    const auto& faces = mesh_.faces(); 
    const auto& faceCenters = mesh_.Cf();

    const auto& meshPoints = mesh_.points();

    const auto& octree = triSurfSearch_.tree();
    const auto& triPoints = triSurf_.points(); 
    const auto& triNormals = triSurf_.faceNormals(); 

    // Full or empty volume fractions for cells, using cell-centered signed distances.
    forAll (cellSignedDist_, cellI)
    {
        if (cellSignedDist_[cellI] > 0)
            alpha[cellI] = 1; 
        else 
            alpha[cellI] = 0;
    }

    // Encode tetrahedron signed distances into a bitset. For 4 tetrahedron
    // points: 1110 == positive, positive, positive, negative distance.
    // Barycentric triangulation of the cell is used: 
    // - First distance is the distance at the cell center. 
    // - Second, third and fourth distances are at two cell corner points. 
    std::bitset<4> dists(pow(2,5) - 1);

    // Legacy VTK output of the cut mesh geometry. 
    geophase::vtkPolyDataOStream cutMeshStream(
        geophase::vtk_file_name("cutMesh", runTime_.timeIndex())
    ); 


    // Correct the volume fractions geometrically in the intersected cells. 
    nTrianglesPerCell_ = 0;
    //forAll(intersectedCellLabels_, cellL) TODO: Remove
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
        {
            // Zero the volume fraction in the intersected cell. 
            alpha[cellI] = 0;
                
            const auto& cutCell = meshCells[cellI];  
            const auto& xC = cellCenters[cellI];
            const auto& distC = cellSignedDist_[cellI];
            // Distance encoding at the center is the sign bit of the cell dist. 
            dists[0] = !std::signbit(distC);
            // For all faces of the cut cell, 
            forAll(cutCell, faceL)
            {
                const label faceG = cutCell[faceL];
                const face& nBandFace = faces[faceG];
                const label point0 = nBandFace[0]; 
                const point& x0 = meshPoints[point0];   
                const scalar dist0 = pointSignedDist_[point0]; 
                dists[1] = !std::signbit(dist0); 

                for (label I = 1; I < (nBandFace.size() - 1); ++I)
                {
                    const label pointI = nBandFace[I]; 
                    const scalar distI = pointSignedDist_[pointI]; 
                    dists[2] = !std::signbit(distI);
                    const point& xI = meshPoints[pointI];

                    const label pointJ = nBandFace[I+1]; 
                    const scalar distJ = pointSignedDist_[pointJ]; 
                    dists[3] = !std::signbit(distJ);
                    const point& xJ = meshPoints[pointJ];

                    const vector xT = 0.25 * (xC + x0 + xI + xJ); // Tetrahedron centroid.
                    const scalar radiusT = 1.1*max(
                            max(mag(xC - xT), mag(x0 - xT)), 
                            max(mag(xI - xT), mag(xJ - xT))
                    ); // Tetrahedron radius.

                    // Intersect the tetrahedron with all the triangles from the sphere. 
                    // - Get labels of triangles that are inside the sphere that contains
                    //   the intersected tetrahedron.
                    //   // sqr radius used for search
                    auto triangleLabels = octree.findSphere(xT, radiusT*radiusT); 

                    //if ((!dists.all()) && dists.any()) // If tetrahedron is intersected
                    // If the tetrahedron sphere intersects triangles from the surface.
                    if (triangleLabels.size()) 
                    {
                        // Build the tetrahedron geometry from points. 
                        auto tetrahedron = 
                            geophase::make_tetrahedron<geophase::foamVectorPolyhedron>(xC, x0, xI, xJ);

                        // Initialize the tetrahedron intersection. 
                        geophase::foamVectorPolyhedron tetIntersection {tetrahedron};
                        nTrianglesPerCell_ += triangleLabels.size();
                        for (const auto& triangleL : triangleLabels)
                        {
                            // Intersect the tetrahedron with the triangle halfspace. 
                            tetIntersection = 
                                intersect_tolerance<geophase::foamPolyhedronIntersection>(
                                    tetIntersection, 
                                    foamHalfspace(
                                        triPoints[triSurf_[triangleL][0]], 
                                        triNormals[triangleL]
                                    )
                            ).polyhedron();
                        }
                        // Add the volume of the intersection to the phase-specific volume.  
                        alpha[cellI] += volume_by_surf_tri(tetIntersection);
                        cutMeshStream << tetIntersection;
                    }
                    else if (dists.all()) // If tetrahedron is inside surface.
                    {
                        // Add the mixed product tet volume to the alpha cell value. 
                        alpha[cellI] += (1. / 6.) * std::abs( 
                            ((x0 - xC) &  ((xI - xC) ^ (xJ - xC)))
                        ); 
                    }

                }
            }
            alpha[cellI] /= cellVolumes[cellI];
        }
    }
    nTrianglesPerCell_ /= intersectedCellLabels_.size();

    // Correct volume fractions. 
    const auto& cellPointLists = mesh_.cellPoints();
    forAll(alpha, cellI)
    {
        if (alpha[cellI] > 1.) 
            alpha[cellI] = 1.;

        const auto& cellPoints = cellPointLists[cellI]; 
        bool cutCell = false; 
        forAll(cellPoints, pointI)
        {
            if (pointSignedDist_[cellPoints[pointI]] > 0)
            {
                cutCell = true; 
                break;
            }
        }

        // if the cell is inside phase 1 
        if (!cutCell && (cellSignedDist_[cellI] < 0))
            alpha[cellI] = 1.;
    }
}

void geomSurfaceCellMeshIntersection::writeFields() const
{
    cellSqrSearchDist_.write();
    cellSignedDist_.write(); 
    cellSignedDist0_.write(); 

    pointSqrSearchDist_.write(); 
    pointSignedDist_.write(); 
}

}} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
