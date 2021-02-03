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
#include "surfaceInterpolate.H"
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
    const bool writeGeo, 
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
    faceSignedDist_(
        IOobject
        (
            "faceSignedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            writeOption 
        ),
        mesh,
        dimensionedScalar("faceSignedDist", dimLength,0)
    ),
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
    writeGeo_(writeGeo),
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
    assert((intersectedCellLabels_.size() < mesh_.nCells())); 
}

void geomSurfaceCellMeshIntersection::calcSignedDist()
{

    // Zero the signed distance.
    cellSignedDist_ = dimensionedScalar("cellSignedDist", dimLength, 0);

    // Find nearest triangles to cells using octree subdivision.
    triSurfSearch_.findNearest(
        mesh_.C(),
        (sqrDistFactor_ * sqrDistFactor_) * cellSqrSearchDist_,
        cellNearestTriangle_
    );

    // Compute distances from cell centers to nearest triangles. 
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

    // Propagate inside-outside information by approximately solving 
    // the Laplace equation for the signed distance.
    fvScalarMatrix distEqn
    (
        -fvm::laplacian(cellSignedDist_)
    );
    distEqn.solve(); 

    // Reset the signed distance in narrow band to geometrical distance.
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
            cellSignedDist_[cellI] = cellSignedDist0_[cellI]; 
    }

    // Calculate the face centered signed distances by linear 
    // interpolation: we only need the sign of the distance.
    faceSignedDist_ = fvc::interpolate(cellSignedDist_, "linear");


    // Interpolate the cell-centered signed distances to cell-corner points
    cellsToPointsInterp_.interpolate(cellSignedDist_, pointSignedDist_);  

    // Find triangles nearest to cell corner-points.
    triSurfSearch_.findNearest(
        mesh_.points(),
        (sqrDistFactor_ * sqrDistFactor_) * pointSqrSearchDist_,
        pointNearestTriangle_ 
    );

    // Correct the cell-corner signed distances with geometrical distances.
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

    // Correct the volume fractions geometrically in intersected cells. 
    nTrianglesPerCell_ = 0;

    findIntersectedCells(); 
    forAll(intersectedCellLabels_, cellJ) 
    {
        const auto cellI = intersectedCellLabels_[cellJ];
        alpha[cellI] = 0;
        const auto& cutCell = meshCells[cellI];  

        // Cell center and its distance as tetrahedron input.
        const auto& xC = cellCenters[cellI];
        const auto& distC = cellSignedDist_[cellI];
        dists[0] = !std::signbit(distC);

        // For all faces of the cut cell, 
        forAll(cutCell, faceL)
        {
            const label faceG = cutCell[faceL];

            // Face center & distance for tetrahedron input.
            const auto& xF = faceCenters[faceG]; 
            const scalar distF = faceSignedDist_[faceG]; 
            dists[1] = !std::signbit(distF); 

            const face& nBandFace = faces[faceG];
            for (label I = 0; I < nBandFace.size(); ++I)
            {
                // Current face-point & dist for tetrahedron input 
                const label pointI0 = nBandFace[I]; 
                const scalar distI0 = pointSignedDist_[pointI0]; 
                dists[2] = !std::signbit(distI0);
                const point& xI0 = meshPoints[pointI0];

                // Next face-point & distance for tetrahedron input 
                const label pointI1 = nBandFace[(I+1) % nBandFace.size()]; 
                const scalar distI1 = pointSignedDist_[pointI1]; 
                dists[3] = !std::signbit(distI1);
                const point& xI1 = meshPoints[pointI1];

                // Tetrahedron centroid.
                const vector xT = 0.25 * (xC + xF + xI0 + xI1); 

                // Tetrahedron sphere radius.
                const scalar radiusT = max(
                        max(mag(xC - xT), mag(xF - xT)), 
                        max(mag(xI0 - xT), mag(xI1 - xT))
                ); 

                // Triangles interesecting tetrahedron sphere. 
                auto triangleLabels = octree.findSphere(xT, radiusT*radiusT); 
                if (dists.all()) // If tetrahedron is inside surface.
                {
                    // Add the mixed product tet volume to the alpha cell value. 
                    alpha[cellI] += (1. / 6.) * std::abs( 
                        ((xF - xC) &  ((xI0 - xC) ^ (xI1 - xC)))
                    ); 
                }
                else if (triangleLabels.size() > 0) 
                {
                    // Initialize the tetrahedron intersection. 
                    geophase::foamVectorPolyhedron tetIntersection {
                        geophase::make_tetrahedron<geophase::foamVectorPolyhedron>(xC, xF, xI0, xI1)
                    };
                    nTrianglesPerCell_ += triangleLabels.size();
                    for (const auto& triangleL : triangleLabels)
                    {
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
                    if (writeGeo_)
                        cutMeshStream << tetIntersection;
                }

            }
        }
        alpha[cellI] /= cellVolumes[cellI];
    }
    nTrianglesPerCell_ /= intersectedCellLabels_.size();
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
