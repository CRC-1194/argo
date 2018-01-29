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

#include "geomSurfaceCellMeshIntersection.H"
#include "calculatedFvPatchField.H"
#include "fvcAverage.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "Geometry.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam { namespace GeometricalTransport {

    defineTypeNameAndDebug(geomSurfaceCellMeshIntersection,0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geomSurfaceCellMeshIntersection::geomSurfaceCellMeshIntersection
(
    const fvMesh& mesh,
    const scalar sqrDistFactor,
    const IOobject::writeOption& wo  // Allows output for testing purposes.
)
:
    mesh_(mesh), 
    runTime_(mesh.time()),  
    sqrSearchDist_
    (
        IOobject
        (
            "sqrSearchDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        fvc::average(pow(mesh.deltaCoeffs(), -2))
    ),
    signedDist_
    (
        IOobject
        (
            "signedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        mesh,
        dimensionedScalar("signedDist", dimLength,0),
        "zeroGradient"
    ),
    signedDist0_("signedDist0", signedDist_),  
    lambda_
    (
        IOobject
        (
            "lambda", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        mesh,
        // Solve a Laplace equation to propagate the sign.
        dimensionedScalar ("lambda", sqr(dimLength) * pow(dimTime,-1), 1)
    ), 
    triSurfPtr_(), // WARNING: Uninitialized pointer, watch out for this when inheriting. 
    cellNearestTriangle_(), 
    sqrDistFactor_(max(3.0, sqrDistFactor)) // Narrow band minimal width = 3 cells. 
{
    cellNearestTriangle_.reserve(signedDist_.size());
}

geomSurfaceCellMeshIntersection::geomSurfaceCellMeshIntersection
(
    const fvMesh& mesh,
    const fileName& triName,
    const scalar sqrDistFactor,
    const IOobject::writeOption& wo // Allows output for testing purposes.
)
:
    geomSurfaceCellMeshIntersection(mesh, sqrDistFactor, wo) 
{
    // Construct the surface from name.
    triSurfPtr_ = new triSurface(triName); 
} 

geomSurfaceCellMeshIntersection::geomSurfaceCellMeshIntersection
(
    const fvMesh& mesh,
    const triSurface& tri,
    const scalar sqrDistFactor,
    const IOobject::writeOption& wo // Allows output for testing purposes.
)
:
    geomSurfaceCellMeshIntersection(mesh, sqrDistFactor, wo) 
{
    // Copy-construct the surface. 
    triSurfPtr_ = new triSurface(tri);
} 

geomSurfaceCellMeshIntersection::geomSurfaceCellMeshIntersection
(
    const geomSurfaceCellMeshIntersection& copy 
)
:
    mesh_(copy.mesh_),
    runTime_(copy.runTime_), 
    sqrSearchDist_(copy.sqrSearchDist_), 
    signedDist_(copy.signedDist_), 
    signedDist0_(copy.signedDist0_), 
    lambda_(copy.lambda_),
    triSurfPtr_(new triSurface(copy.triSurfPtr_())), // Deep copy. 
    cellNearestTriangle_(copy.cellNearestTriangle_), 
    sqrDistFactor_(copy.sqrDistFactor_)
{} 

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void geomSurfaceCellMeshIntersection::calcSqrSearchDist()
{
    sqrSearchDist_ = fvc::average(pow(mesh_.deltaCoeffs(), -2));
}

void geomSurfaceCellMeshIntersection::calcSignedDist(
    const triSurface& tri, 
    const triSurfaceSearch& triSearch
)
{
    // Zero the signed distance.
    signedDist_ = dimensionedScalar("signedDist", dimLength, 0);

    // Use the octree and the square search distance multiplied by a distance factor 
    // to build the surface mesh / volume mesh proximity information list. 
    triSearch.findNearest(
        mesh_.C(),
        sqrDistFactor_ * sqrDistFactor_ * sqrSearchDist_,
        cellNearestTriangle_
    );

    // Compute the signed distance in each cell as the distance between the triangle
    // nearest to the cell center and the cell center.
    const volVectorField& C = mesh_.C();  
    const pointField& triPoints = tri.points();  
    const vectorField& triNormals = tri.faceNormals(); 
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
        {
            signedDist_[cellI] = 
                (C[cellI] - triPoints[tri[cellHit.index()][0]]) & 
                triNormals[cellHit.index()];
        }
    }

    // Save the signed distance field given by the octree.
    signedDist0_ = signedDist_; 

    // Propagate the sign information into the bulk by solving a Laplace
    // equation for a single iteration for the signed distance field. 
    fvScalarMatrix distEqn
    (
        -fvm::laplacian(lambda_, signedDist_)
    );
    distEqn.solve(); 
}

void geomSurfaceCellMeshIntersection::calcSignedDist()
{
    if (triSurfPtr_.empty())
    {
         FatalErrorInFunction
            << "Empty triSurfPtr_."
            << abort(FatalError);
    }
    else
    {
        triSurfaceSearch triSearch (triSurfPtr_()); 
        calcSignedDist(triSurfPtr_(), triSearch); 
    }
}

void geomSurfaceCellMeshIntersection::calcVolFraction(
    volScalarField& alpha, 
    const triSurface& tri
)
{
    // Build the octree around the triSurface. 
    triSurfaceSearch triSearch(tri);

    // Compute the signed distance using the octree in triSurfaceSearch. 
    calcSignedDist(tri, triSearch); 

    // Compute the volume fraction based on the signed distance field in 
    // the bulk, and perform geometrical intersections of the cells cut 
    // by the surface. 
    
    // In the test mode, prepare the mesh / surface mesh intersection output.
#ifdef TESTING
    vtk_polydata_stream cutCellStream(prependVtkFileName("cutCells", runTime_.timeIndex())); 
#endif

    const cellList& cells = mesh_.cells(); 
    const auto& octree = triSearch.tree();
    const auto& triPoints = tri.points(); 
    const auto& triNormals = tri.faceNormals(); 
    const auto& V = mesh_.V(); 
    forAll(cellNearestTriangle_, cellI)
    {
        if (signedDist_[cellI] > 0)
            alpha[cellI] = 1; 

        // Correct boundary oscillation using the octree distance field.
        if (signedDist0_[cellI] < 0)
            alpha[cellI] = 0; 

        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        // If a cell is in a narrow band.
        if (cellHit.hit()) 
        {
            // Get the cell points.
            const pointField cellPoints = cells[cellI].points(mesh_.faces(), mesh_.points());  
            // Find all the triangles in the bounding box of a cell in a narrow band. 
            labelList cellTriangles = octree.findBox(treeBoundBox(cellPoints)); 

            // If there are triangles in the cell.
            if (!cellTriangles.empty())
            {
                // Initialize the cell intersection as the barycentri triangulation of
                // a cell: necessary for non-convex polyhedral cells.
                triangulationIntersection cellIntersection
                (
                    barycentric_triangulate<tetrahedronVector>
                    (
                        build<pointVectorVector>(cellI, mesh_)
                    )
                );

                // For all triangles in a cell
                for (const auto triLabel : cellTriangles)
                {
                    // Intersect the cell triangulation with the triangle halfspace. 
                    cellIntersection = intersect<triangulationIntersection>(
                        cellIntersection, 
                        //halfspace(triPoints[triangles[triLabel][0]], triNormals[triLabel])
                        halfspace(triPoints[tri[triLabel][0]], triNormals[triLabel])
                    );  
                }

                // Set the volume fraction of the cell.  
                alpha[cellI] = min(1, volume(cellIntersection) / V[cellI]); 

// Visualize all the cut cells as well as each cell, a piece of the surface that intersects
// the cell, as well as the intersection result per cell.
#ifdef TESTING 
                cutCellStream << cellIntersection;
                // Uncomment for debugging. 
                // For every cell that is intersected: 
                // - write the part of the triSurface that intersects it into an stl file
                List<labelledTri> cellTriangleGeo(cellTriangles.size()); 
                forAll(cellTriangleGeo, I)
                    //cellTriangleGeo[I] = triangles[cellTriangles[I]];
                    cellTriangleGeo[I] = tri[cellTriangles[I]];
                triSurface cellSurface(cellTriangleGeo, triPoints); 
                cellSurface.write(appendSuffix("cellSurface", cellI) + ".stl");
                // - write the intersection into a .vtk file. 
                write_vtk_polydata(cellIntersection, appendSuffix("cellIntersection", cellI) + ".vtk");
                write_vtk_polydata
                (
                    build<pointVectorVector>(cellI,mesh), 
                    appendSuffix("cell", cellI) + ".vtk"
                );
#endif
            }
        }
    }
}

void geomSurfaceCellMeshIntersection::calcVolFraction(volScalarField& alpha)
{
    calcSignedDist(); 
    calcVolFraction(alpha, triSurfPtr_()); 
}

void geomSurfaceCellMeshIntersection::writeFields() const
{
    sqrSearchDist_.write();
    signedDist_.write(); 
    signedDist0_.write(); 
}

// * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * * //


void geomSurfaceCellMeshIntersection::operator=(const geomSurfaceCellMeshIntersection& rhs)
{
    if (&rhs == this)
    {
         FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
    sqrSearchDist_ = rhs.sqrSearchDist_; 
    signedDist_ = rhs.signedDist_; 
    signedDist0_ = rhs.signedDist0_; 
    lambda_ = rhs.lambda_; 
    // Deep copy on assignment.
    triSurfPtr_ = new triSurface(rhs.triSurfPtr_()); 
    cellNearestTriangle_ = rhs.cellNearestTriangle_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
