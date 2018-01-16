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
    voFoamSetSurfaceFraction 

Description
    Intersects a surface mesh with the volume mesh by intersecting cell
    triangulations with a set of nearest surface triangles that intersect the cell. 

    The volume fraction is then given as the ratio of the volume of the intersection 
    and the volume of the cell.  

    Test description: 
    
    Use bounding boxes of both meshes to randomly place the tool mesh within
    the base mesh such that the tool mesh bounding box is contained within the
    base mesh bounding box.  

    Intersect two meshes to obtain the volume fraction of the tool mesh in 
    the base mesh.

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Geometry.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "OFstream.H"
#include <set>
#include <omp.h>

#ifndef _OPENMP
    #define omp_get_wtime() 0
#endif

#include "Random.H"

using namespace GeometricalTransport;

#include "surfaceMeshIntersection.H"

typedef TriangulationIntersection<tetrahedronVector,pointVectorVector> triangulationIntersection;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "createOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label distanceFactor = args.optionLookupOrDefault<scalar>("distanceFactor", 1.5); 
    const label radiusFactor = args.optionLookupOrDefault<scalar>("radiusFactor", 1.71); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "surfaceCellMeshIntersection.dat"); 

    fileName triFile = args.path() + "/surface.stl";

    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    // Triangle surface
    triSurface tri(triFile);

    #include "createFields.H"

    // Initialize mesh geometrical data. This data uses lazy evaluation:
    // it is stored with raw pointers set to null until initialized. If
    // every thread calls these functions, a race condition leads to 
    // undefined behavior. TM. 
    // TODO: Check for large cases if this needs speeding up with 
    // tasks. TM.
    #pragma omp parallel 
    #pragma omp single nowait
    {
        tri.faceFaces();
        tri.faceCentres(); 
        tri.faceNormals(); 
        mesh.deltaCoeffs();
        mesh.cells();
        mesh.faces();
        mesh.Sf();
        mesh.C();
        mesh.V(); 
        mesh.faceOwner();
    }

    const pointField& triPoints = tri.points();

    auto& triNormals = const_cast<pointField&>(tri.faceNormals());
    auto& triangles = static_cast<List<labelledTri>& >(tri);

    const Switch fixNormals = args.optionFound("fixNormals");

    // Correct normal orientation based on the centroid position. 
    // Assumes a star-shaped tri-surface that can be oriented using the centroid.
    if(fixNormals)
    {
        vector surfaceMeshCentroid = sum(triPoints) / tri.nPoints();  
        const auto& triCentres = tri.faceCentres(); 

        #pragma omp parallel for 
        forAll(triNormals, tI)
        {
            if (((triCentres[tI] - surfaceMeshCentroid) & triNormals[tI]) > 0) 
            {
                triNormals[tI] *= -1; 
                triangles[tI].flip(); 
            }
        }
    }

    // RANDOM tool mesh positioning. 
    // Position the tool mesh centroid at the base mesh centroid. 
    // NOTE: triSurface is not decomposed!
    label nBaseCells = mesh.nCells();  
    Pstream::gather(nBaseCells, sumOp<label>()); 
    const vector baseCenter = gSum(mesh.C()) / nBaseCells;  

    vector toolCenter = sum(tri.points()) / tri.size(); 
    // This invalidates surface normals!
    tri.movePoints(tri.points() + baseCenter - toolCenter); 

    boundBox baseBox = mesh.bounds(); 
    const boundBox toolBox = boundBox(tri.points()); 
    const Vector<label> solutionVector = mesh.solutionD(); 
    vector toolBoxDelta = toolBox.max() - toolBox.min(); 
    // Null the toolBoxDelta in a dimension that is not used in pseudo2D. 
    forAll(toolBoxDelta, I)
    {
        if (solutionVector[I] == -1)
            toolBoxDelta[I] *= 0;
    }
    point& baseMin = baseBox.min(); 
    point& baseMax = baseBox.max(); 
    // Shrink the base mesh box, so that the random position of the tool mesh
    // centroid is placed in a way that the base mesh bounding box contains  
    // the tool mesh bounding box. 
    baseMin += toolBoxDelta; 
    baseMax -= toolBoxDelta;  

    // Generate a random tool mesh centroid position within the shrunk base // mesh bounding box. 
    Random r(1e04); 

    // Open the error file for measurement output..
    OFstream errorFile(dataFileName); 
    // Nt : number of cells in the tool mesh 
    // Nb : number of cells in the tool mesh 
    // Ev : volume conservation error: 
    // |tool mesh volume from volume fraction - tool mesh volume | / tool mesh volume. 
    // Te : execution time of the CCI mesh intersection operation.
    errorFile << "Nt,Nb,Ev,Te\n"; 

    while(runTime.run())
    {
        // Zero the volume fraction and signed distance fields.
        alpha = dimensionedScalar("alpha", dimless, 0);  
        signedDist = dimensionedScalar("signedDist", dimLength, 0);

        // Compute the random displacement with respect to the base mesh centroid.
        const vector randomPosition = r.position(baseMin, baseMax); 
        vector randomDisplacement = randomPosition - baseCenter; 
        // Null the displacement in a dimension that is not used in pseudo2D. 
        forAll(randomDisplacement, I)
        {
            if (solutionVector[I] == -1)
                randomDisplacement[I] *= 0;
        }

        if (! baseBox.contains(randomPosition))
        {
            // Remove
            //Info << baseMin << ":" << baseMax << ":" << randomPosition << nl;
            FatalErrorIn("voFoamTestCellCellIntersectMeshes")                                                 
                << "Random position outside of the base mesh bounding box." 
                << Foam::abort(FatalError);
        }
        tri.movePoints(triPoints + randomDisplacement); 

        tri.write(appendSuffix("triMesh", runTime.timeIndex()) + ".stl");

        // Moving a triSurface clears all geometrical data, so it needs to be
        // re-calculated outside of a parallel region, because the calculation
        // loop is not threaded in OpenFOAM.
        tri.faceNormals(); 

        const scalar t0 = omp_get_wtime(); 

        // Build the octree around the triSurface. 
        triSurfaceMesh triMesh(
            IOobject(
                "triMesh",
                "constant",
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tri 
        );

        DynamicList<pointIndexHit> cellNearestTriangle;
        cellNearestTriangle.reserve(mesh.nCells()); 

        calcSearchFields(sqrSearchDist,sqrCellRadius,distanceFactor,radiusFactor); 

        triMesh.findNearest(
            mesh.C(),
            sqrSearchDist,
            cellNearestTriangle 
        );

        // Use octree search to find a nearest triangle for each cell point.
        DynamicList<pointIndexHit> pointNearestTriangle;
        pointNearestTriangle.reserve(mesh.nPoints()); 

        // Compute the marker field. 
        const auto& C = mesh.C();
        const auto& V = mesh.V(); 

#ifdef TESTING
        vtk_polydata_stream cutCellStream(prependVtkFileName("cutCells", runTime.timeIndex())); 
#endif
        const auto& octree = triMesh.tree();

        // Distance calculation and setting of alpha1 must precede geometrical 
        // intersections because there may exist cells with positive signed 
        // distance that are intersected and whose alpha1 value must then 
        // be subsequently corrected. TM.
        const auto& triConstNormals = tri.faceNormals(); 
        
        #pragma omp parallel for schedule (dynamic) 
        forAll(cellNearestTriangle, cellI)
        {
            const pointIndexHit& cellHit = cellNearestTriangle[cellI];

            if (cellHit.hit()) 
            {
                signedDist[cellI] = 
                    (C[cellI] - triPoints[tri[cellHit.index()][0]]) & 
                    triConstNormals[cellHit.index()];
            }
        }

#ifdef TESTING
        volScalarField signedDist0("signedDist0", signedDist); 
        signedDist0.write(); 
#endif

        // Solve a laplace equation to propagate the sign.
        dimensionedScalar lambda ("lambda", sqr(dimLength) * pow(dimTime,-1), 1);
        fvScalarMatrix distEqn
        (
            -fvm::laplacian(lambda, signedDist)
        );
        distEqn.solve(); 

        // Set the fill level of all cells with the positive
        // signed distance to 1. 
        #pragma omp parallel for schedule (dynamic) 
        forAll(signedDist, cellI)
        {
            if (signedDist[cellI] > 0)
                alpha[cellI] = 1; 
        }

        //#pragma omp parallel for schedule (dynamic) 
        forAll(cellNearestTriangle, cellI)
        {
            const pointIndexHit& cellHit = cellNearestTriangle[cellI];

            if (cellHit.hit()) 
            {
                labelList cellTriangles = octree.findSphere(C[cellI],sqrCellRadius[cellI]);


                if (!cellTriangles.empty())
                {
                    // Uncomment for debugging. This will write out a tri-surface chunk 
                    // for every cell that is intersected with it.  
                    /*
                    // Write the triSurface chunk for cellI.
                    List<labelledTri> cellTriangleGeo(cellTriangles.size()); 
                    forAll(cellTriangleGeo, I)
                        cellTriangleGeo[I] = triangles[cellTriangles[I]];
                    triSurface cellSurface(cellTriangleGeo, triPoints); 
                    cellSurface.write(appendSuffix("cellSurface", cellI) + ".stl");
                    // Write the polyhedron for cellI. 
                    build<pointVectorVector>(cellI, mesh)
                    TODO: 
                    */

                    triangulationIntersection cellIntersection
                    (
                        barycentric_triangulate<tetrahedronVector>
                        (
                            build<pointVectorVector>(cellI, mesh)
                        )
                    );

                    // Intersect the cell triangulation with all the triangle halfspaces. 
                    for (const auto triLabel : cellTriangles)
                    {
                        cellIntersection = intersect<triangulationIntersection>(
                            cellIntersection, 
                            halfspace(triPoints[triangles[triLabel][0]], triConstNormals[triLabel])
                        );  
                    }
                    // Calculate volume fraction value that is bounded from above by 1. 
                    alpha[cellI] = min(1, volume(cellIntersection) / V[cellI]); 
#ifdef TESTING
                    #pragma omp critical
                    cutCellStream << cellIntersection;
#endif
                }
            }
        }

        const scalar t1 = omp_get_wtime(); 

        // Write the volume fraction field.
        alpha.write(); 
        signedDist.write(); 
        
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        const scalar Te = t1 - t0; 

        toolCenter = sum(tri.points()) / tri.size(); 
        scalar Vs = 0; 

        // 3D and pseudo 2D.
        // Compute the surface mesh volume using the surface mesh centroid. 
        // Requires a star-shaped surface mesh: one that can be decomposed into tets using 
        // the centroid.
        forAll(tri, tI)
        {
            const point& p0 = triPoints[triangles[tI][0]];
            const point& p1 = triPoints[triangles[tI][1]];
            const point& p2 = triPoints[triangles[tI][2]];

            // Assuming inward pointing normals for distance calculation!
            // Volumes will be negative!
            scalar deltaV = -1.0 / 6.0 * (p0 - toolCenter) & 
                ((p1 - toolCenter) ^ (p2 - toolCenter)); 

            if (deltaV < 0) 
            {
                write_vtk_polydata(build<arrayTriangle>(p0,p1,p2), "bad-triangle.vtk");
                Perr << "ERRROR:" 
                    << "Negative mesh volume contribution!"  
                    << "deltaV = " << deltaV << endl
                    << "toolCenter = " << toolCenter << endl
                    << "p0 = " << p0 << endl
                    << "p1 = " << p1 << endl
                    << "p2 = " << p2 << endl
                    << Foam::abort(FatalError);
            }
            else 
                Vs += deltaV;
        }

        if (mesh.nSolutionD() == 2)
        {
            Info << "This is a 2D case!" << endl;
            // Pseudo 2D: update the total volume for the edges that have only a single
            // edge owner.
            
            // Find the hight of the surface mesh using the bounding box height
            // and mesh solution vector. 
            const boundBox toolBox = boundBox(tri.points());  
            const vector toolMin = toolBox.min();  
            const vector toolMax = toolBox.max(); 

            vector toolHeightVector = 0.5 * (toolMax - toolMin);  
            forAll(solutionVector, sI)
                if (solutionVector[sI] == 1)
                    toolHeightVector[sI] = 0; 

            // For each edge. 
            const edgeList& edges = tri.edges(); 
            const pointField& localPoints = tri.localPoints();
            forAll(edges, edgeI)
            {
                if (! tri.isInternalEdge(edgeI))
                {
                    const vector& ePoint0 = localPoints[edges[edgeI][0]]; 
                    const vector& ePoint1 = localPoints[edges[edgeI][1]]; 

                    Vs +=mag
                    (
                        (1.0 / 6.0) * 
                        ((ePoint0 - toolCenter) ^ (ePoint1 - toolCenter))
                        & toolHeightVector
                    ); 
                }
            }
        }

        const scalar Valpha = sum(alpha * mesh.V()).value();
        const scalar Ev = mag(Vs - Valpha) / Vs; 

        Info<< "Volume from surface mesh = " << Vs << endl;
        Info<< "Volume from volume fraction = " << Valpha << endl;
        Info<< "Volume error = " << Ev << endl; 
        Info<< "Execution time = " << Te << endl << endl;

        // Move the tri surface back to the original position. 
        tri.movePoints(tri.points() - randomDisplacement); 

        errorFile << triMesh.size() << "," 
            << mesh.nCells() << ","
            << Ev << ","  
            << Te << "\n";  

        runTime++; 
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
