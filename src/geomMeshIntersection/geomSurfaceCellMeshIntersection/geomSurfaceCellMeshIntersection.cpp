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
#include "Polyhedron.hpp"
#include "Make.hpp"
#include "GeophaseMake.hpp"

#include <cassert>

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
    sqrDistFactor_(max(2.0, sqrDistFactor)) 
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

    Info << intersectedCellLabels_.size() << endl;
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
}

void geomSurfaceCellMeshIntersection::calcVolFraction(volScalarField& alpha)
{
    // Compute the signed distances using the octree in triSurfaceSearch. 
    calcSignedDist(); 

    // ALGORITHM 
    // Compute the volume fraction based on the cell centered signed distance
    // field in the bulk, and perform geometrical intersections of the cells
    // cut by the surface. 

    // For every cell in the narrow band. 
    
        // Use a centroid triangulation for the cell to build tets.  
        // Every tet that has a sign change in the signed distance of its
        // two points if it is intersected.

        // Those that have only positive distances are inside and therefore
        // full, otherwise they are empty.
        
        // For each tet that is intersected, build a bounding sphere using
        // the maximal centroid-vertex length as radius.  
        
        // For the tet sphere, get triangles from the STL surface that
        // intersect the sphere. 
        
        // For each triangle intersecting the tet box, intersect the tetrahedron. 
       
        // Add the volume of the intersection to the cell. 
    
    // FIXME: Use CMAKE Debug release flag.
    // In the testing mode, prepare the mesh / surface mesh intersection output.
//#ifdef 
    //vtk_polydata_stream cutCellStream(prependVtkFileName("cutCells", runTime_.timeIndex())); 
//#endif

    const auto& meshCells = mesh_.cells(); 
    const auto& cellCenters = mesh_.C();  
    const auto& cellVolumes = mesh_.V(); 

    const auto& faces = mesh_.faces(); 
    const auto& faceCenters = mesh_.Cf();

    const auto& meshPoints = mesh_.points();

    const auto& octree = triSurfSearch_.tree();
    const auto& triPoints = triSurf_.points(); 
    const auto& triNormals = triSurf_.faceNormals(); 

    
    // Set volume fraction values based on the cell signed distance. 
    forAll (cellSignedDist_, cellI)
    {
        if (cellSignedDist_[cellI] < 0)
            alpha[cellI] = 1; 
        else 
            alpha[cellI] = 0;
    }

    // Encode tetrahedron signed distances into a bitset.
    std::bitset<4> dists(pow(2,5) - 1);

    // Correct the volume fractions geometrically in the intersected cells. 
    forAll(intersectedCellLabels_, cellL)
    {
        // Global cell label.
        const label cellI = intersectedCellLabels_[cellL];
        Info << "Cell " << cellI << endl;

        // Use the geometrical signed distance in the intersected cells.
        cellSignedDist_[cellI] = cellSignedDist0_[cellI];
        alpha[cellI] = 0;

        // TODO: Check this part. TM.
        // Correct boundary oscillation using the octree distance field.
        //if (cellSignedDist0_[cellI] < 0)
            //alpha[cellI] = 1; 
        //else if (cellSignedDist0_[cellI] > 0)
            //alpha[cellI] = 0; 
            
        const auto& cutCell = meshCells[cellI];  
        const auto& xCell = cellCenters[cellI];
        const auto& distCell = cellSignedDist_[cellI];
        dists[0] = std::signbit(distCell);
        forAll(cutCell, faceI)
        {
            const label faceG = cutCell[faceI];
            const face& nBandFace = faces[faceG];
            const label point0 = nBandFace[0]; 
            const point& x0 = meshPoints[point0];   
            const scalar dist0 = pointSignedDist_[point0]; 
            dists[1] = std::signbit(dist0); 
            for (label pointI = 1; pointI < nBandFace.size() - 1; ++pointI)
            {
                const label pointM = nBandFace[pointI]; 
                const scalar distM = pointSignedDist_[pointM]; 
                // Set the inside / outside flag for pointM
                dists[2] = std::signbit(distM);
                const point& xM = meshPoints[pointM];

                const label pointN = nBandFace[pointI + 1];
                const scalar distN = pointSignedDist_[pointN]; 
                // Set the inside / outside flag for pointN 
                dists[3] = std::signbit(distN);
                const point& xN = meshPoints[pointN];


                // If the tet is intersected
                if ((!dists.all()) && dists.any())
                {
                    const vector xTet = 0.25 * (xCell + x0 + xM + xN);
                    const scalar radTet = max(
                            max(mag(xCell - xTet), mag(x0 - xTet)), 
                            max(mag(xM - xTet), mag(xN - xTet))
                    );
                    
                    // Fetch labels of triangles that intersect the tet sphere. 
                    auto triangleLabels = octree.findSphere(xTet, radTet);
                    
                    // TODO: Make the tetrahedron build work in the test app. 
                    // Build the tetrahedron from foam vectors. 
                    //auto tetrahedron = 
                        //geophase::make_tetrahedron3D<geophase::foamVectorPolyhedron>(
                        //xCell, 
                        //point0, 
                        //pointM,
                        //pointN
                    //);

                    // For all triangles in the sphere.

                        // Intersect the tet with the triangle halfspace. 
                        
                    // Add the intersected volume to the cell volume. 
                }
                // Else if the tet is inside
                else if (dists.all())
                {
                    // Add the mixed product tet volume to the alpha cell value. 
                    scalar tetVolume = (1. / 6.) * std::abs( 
                        ((x0 - xCell) &  ((xM - xCell) ^ (xN - xCell)))); 
                    std::cout << dists << ":" << tetVolume << std::endl;
                    alpha[cellI] += tetVolume; 
                }
            }
        }
        std::cout << "Volume = " << alpha[cellI] << std::endl;
        alpha[cellI] /= cellVolumes[cellI];
        std::cout << "Alpha = " << alpha[cellI] << std::endl;
    }

    //// FIXME: old code, remove this. TM.
            //////// Obtain cell points.
            //////const pointField cellPoints = cells[cellI].points(
                    //////mesh_.faces(), 
                    //////mesh_.points()
            //////);  

            //////// Find all the surface mesh triangles in the bounding box of a
            //////// cell in a narrow band. 
            //////labelList cellTriangles = octree.findBox(treeBoundBox(cellPoints)); 

            //////// FIXME: Add geophase geometry. 
            //////// If there are triangles in the cell.
            //////if (!cellTriangles.empty())
            //////{
                //////// Triangulate the cell. 
                
                ////////auto surfCellIntersection 
                    ////////= make<triangulationIntersection>(cellI, mesh_); 

                ////////triangulationIntersection cellIntersection
                ////////(
                    ////////barycentric_triangulate<tetrahedronVector>
                    ////////(
                        ////////make<pointVectorVector>(cellI, mesh_)
                    ////////)
                ////////);

                //////// For all triangles in a cell
                //////for (const auto triLabel : cellTriangles)
                //////{
                    //////// Intersect the cell triangulation with the triangle halfspace. 
                    //////surfCellIntersection = intersect_tolerance<triangulationIntersection>(
                        //////triangulationIntersection, 
                        ////////halfspace(triPoints[triangles[triLabel][0]], triNormals[triLabel])
                        //////halfspace(triPoints[tri[triLabel][0]], -triNormals[triLabel])
                    //////);  
                    //////Nx_++; 
                //////}

                //////// Set the volume fraction of the cell.  
                //////alpha[cellI] = min(1, volume(cellIntersection) / V[cellI]); 

//////// Visualize all the cut cells as well as each cell, a piece of the surface that intersects
//////// the cell, as well as the intersection result per cell.
//////#ifdef TESTING 
                //////cutCellStream << cellIntersection;
                //////// Uncomment for debugging. 
                //////// For every cell that is intersected: 
                //////// - write the part of the triSurface that intersects it into an stl file
                //////List<labelledTri> cellTriangleGeo(cellTriangles.size()); 
                //////forAll(cellTriangleGeo, I)
                    ////////cellTriangleGeo[I] = triangles[cellTriangles[I]];
                    //////cellTriangleGeo[I] = tri[cellTriangles[I]];
                //////triSurface cellSurface(cellTriangleGeo, triPoints); 
                //////cellSurface.write(appendSuffix("cellSurface", cellI) + ".stl");
                //////// - write the intersection into a .vtk file. 
                //////write_vtk_polydata(cellIntersection, appendSuffix("cellIntersection", cellI) + ".vtk");
                //////write_vtk_polydata
                //////(
                    //////build<pointVectorVector>(cellI, mesh_), 
                    //////appendSuffix("cell", cellI) + ".vtk"
                //////);
//////#endif
            ////}
        ////}
    ////}
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
