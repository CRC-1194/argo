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
//#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"
#include <set>

#ifndef _OPENMP
    #define omp_get_wtime() 0
#endif

// Time measurement
#include <chrono>
// Random tri surface mesh positioning
#include <random>

using namespace GeometricalTransport;
using namespace std::chrono;

typedef TriangulationIntersection<tetrahedronVector,pointVectorVector> triangulationIntersection;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point p0 = high_resolution_clock::now();

    #include "createOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label sqrDistFactor = args.optionLookupOrDefault<scalar>("sqrDistFactor", 3); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "surfaceCellMeshIntersection.csv"); 

    fileName triFile = args.path() + "/surface.stl";

    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    // Triangle surface
    triSurface tri(triFile);

    #include "createFields.H"

    high_resolution_clock::time_point p1 = high_resolution_clock::now();

    // Initialize mesh geometrical data. This data uses lazy evaluation:
    // it is stored with raw pointers set to null until initialized. If
    // every thread calls these functions, a race condition leads to 
    // undefined behavior. TM. 
    // TODO: Check for large cases if this needs speeding up with 
    // tasks. TM.
    tri.faceNormals(); 
    mesh.deltaCoeffs();
    mesh.cells();
    mesh.faces();
    mesh.Sf();
    mesh.C();
    mesh.V(); 
    mesh.faceOwner();

    pointField& triPoints = const_cast<pointField&>(tri.points());
    pointField& triLocalPoints = const_cast<pointField&>(tri.localPoints());

    auto& triNormals = const_cast<pointField&>(tri.faceNormals());
    auto& triangles = static_cast<List<labelledTri>& >(tri);

    const Switch fixNormals = args.optionFound("fixNormals");

    // Correct normal orientation based on the centroid position. 
    // Assumes a star-shaped tri-surface that can be oriented using the centroid.
    if(fixNormals)
    {
        vector surfaceMeshCentroid = sum(triPoints) / tri.nPoints();  
        const auto& triCentres = tri.faceCentres(); 

        //#pragma omp for 
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
    label nBasePoints = mesh.nPoints();  
    Pstream::gather(nBasePoints, sumOp<label>()); 
    const vector baseCenter = gSum(mesh.points()) / nBasePoints;  
    vector toolCenter = sum(triPoints) / tri.nPoints(); 

#ifdef TESTING
    Info << "Base mesh centroid = " << baseCenter << nl 
        << "Initial tool mesh centroid = "  << toolCenter << nl 
        << "Displacement = " << baseCenter - toolCenter << endl;

    Info << triPoints.size()  << " : " << triLocalPoints.size() << endl;
#endif 

    // Move the triSurface 
    triPoints += (baseCenter - toolCenter);
    // Displace also local points used for addressing.
    triLocalPoints += (baseCenter - toolCenter);

#ifdef TESTING // Write the initial tri surface and test the centroid position.
    tri.write(appendSuffix("trisurface-initial", runTime.timeIndex()) + ".stl");

    // New tool center after motion.
    toolCenter = sum(triPoints) / tri.nPoints(); 
    const scalar centroidDiff = mag(baseCenter - toolCenter);  

    // Tolerance increased because of the roundoff error.
    if (centroidDiff > 1e03 * SMALL)
        FatalErrorIn("voFoamTestCellCellIntersectMeshes")                                                 
            << "Failed initial tool mesh positioning at the base mesh centroid." << nl 
            << "Final tool mesh centroid = " << toolCenter << nl 
            << "Final centroid diff = " << baseCenter - toolCenter << nl 
            << "Final centroid diff magnitude = " << centroidDiff << nl 
            << Foam::abort(FatalError);
#endif

    boundBox baseBox = mesh.bounds(); 
    const boundBox toolBox = boundBox(tri.points()); 
    const Vector<label> solutionVector = mesh.solutionD(); 
    vector toolBoxDelta = 0.5 * (toolBox.max() - toolBox.min()); 
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

    // Random number initialization.
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);

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
        // Zero the volume fraction, signed distance and radius.
        alpha = dimensionedScalar("alpha", dimless, 0);  

        // Compute the random displacement with respect to the base mesh centroid.
        vector randomPosition = baseMin;  
        forAll(randomPosition, cmptI)
            randomPosition[cmptI] += dis(gen) * (baseMax[cmptI] - baseMin[cmptI]);

        vector randomDisplacement = randomPosition - baseCenter; 
        // TODO: Remove, debugging
        Info << "Random displacement = " << randomDisplacement << endl;
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
        triPoints += randomDisplacement;
        triLocalPoints += randomDisplacement;

#ifdef TESTING
        tri.write(appendSuffix("trisurface", runTime.timeIndex()) + ".stl");
#endif
        //Moving a triSurface clears all geometrical data, so it needs to be
        //re-calculated outside of a parallel region, because the calculation
        //loop is not threaded in OpenFOAM.

        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        // Build the octree around the triSurface. 
        triSurfaceSearch triMesh(tri);

        DynamicList<pointIndexHit> cellNearestTriangle;
        cellNearestTriangle.reserve(mesh.nCells()); 

        sqrSearchDist = fvc::average(Foam::pow(mesh.deltaCoeffs(), -2));

        high_resolution_clock::time_point q0 = high_resolution_clock::now();

        triMesh.findNearest(
            mesh.C(),
            sqrDistFactor * sqrDistFactor * sqrSearchDist,
            cellNearestTriangle 
        );

        high_resolution_clock::time_point q1 = high_resolution_clock::now();

        // Use a single-cell wide search distance for intersections. 
        // Widen the squared search distance for sign propagation to ensure 
        // numerical stability when solving the Laplace equation.

        // Use octree search to find a nearest triangle for each cell point.
        //DynamicList<pointIndexHit> pointNearestTriangle;
        //pointNearestTriangle.reserve(mesh.nPoints()); 

        // Compute the marker field. 
        const auto& C = mesh.C();
        const auto& V = mesh.V(); 


        // Distance calculation and setting of alpha1 must precede geometrical 
        // intersections because there may exist cells with positive signed 
        // distance that are intersected and whose alpha1 value must then 
        // be subsequently corrected. TM.
        //const auto& triConstNormals = tri.faceNormals(); 

        forAll(cellNearestTriangle, cellI)
        {
            const pointIndexHit& cellHit = cellNearestTriangle[cellI];

            if (cellHit.hit()) 
            {
                signedDist[cellI] = 
                    (C[cellI] - triPoints[tri[cellHit.index()][0]]) & 
                    triNormals[cellHit.index()];
            }
        }

        // Initial signed distance field given by the octree.
        volScalarField signedDist0("signedDist0", signedDist); 

        high_resolution_clock::time_point r0 = high_resolution_clock::now();
        fvScalarMatrix distEqn
        (
            -fvm::laplacian(lambda, signedDist)
        );
        distEqn.solve(); 
        high_resolution_clock::time_point r1 = high_resolution_clock::now();

        high_resolution_clock::time_point s0 = high_resolution_clock::now();
        const cellList& cells = mesh.cells(); 
        const auto& octree = triMesh.tree();
        forAll(cellNearestTriangle, cellI)
        {
            if (signedDist[cellI] > 0)
                alpha[cellI] = 1; 
            // Correct boundary oscillation using the initial field.
            // NOTE: not connected with above test!
            if (signedDist0[cellI] < 0)
                alpha[cellI] = 0; 

            const pointIndexHit& cellHit = cellNearestTriangle[cellI];

            if (cellHit.hit()) 
            {
                const pointField cellPoints = cells[cellI].points(mesh.faces(), mesh.points());  
                labelList cellTriangles = octree.findBox(treeBoundBox(cellPoints)); 

                if (!cellTriangles.empty())
                {

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
                            halfspace(triPoints[triangles[triLabel][0]], triNormals[triLabel])
                        );  
                    }
                    // Calculate volume fraction value that is bounded from above by 1. 
                    alpha[cellI] = min(1, volume(cellIntersection) / V[cellI]); 
#ifdef TESTING
                    cutCellStream << cellIntersection;
                    // Uncomment for debugging. 
                    // For every cell that is intersected: 
                    // - write the part of the triSurface that intersects it into an stl file
                    List<labelledTri> cellTriangleGeo(cellTriangles.size()); 
                    forAll(cellTriangleGeo, I)
                        cellTriangleGeo[I] = triangles[cellTriangles[I]];
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
        high_resolution_clock::time_point s1 = high_resolution_clock::now();

        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        // alpha is not written in the test application. 
#ifdef TESTING
        alpha.write(); 
        signedDist.write();
        signedDist0.write();
#endif
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        toolCenter = sum(tri.points()) / tri.nPoints(); 
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
                Perr << "Negative mesh volume contribution!"  
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
            // Pseudo 2D:fields update the total volume for the edges that have only a single
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
            forAll(edges, edgeI)
            {
                if (! tri.isInternalEdge(edgeI))
                {
                    const vector& ePoint0 = triLocalPoints[edges[edgeI][0]]; 
                    const vector& ePoint1 = triLocalPoints[edges[edgeI][1]]; 

                    Vs +=mag
                    (
                        (1.0 / 6.0) * 
                        ((ePoint0 - toolCenter) ^ (ePoint1 - toolCenter))
                        & toolHeightVector
                    ); 
                }
            }
        }

        const scalar Te = duration_cast<nanoseconds>(t1 - t0).count() / 1e09;
        const scalar Ti = duration_cast<nanoseconds>(p1 - p0).count() / 1e09;
        const scalar To = duration_cast<nanoseconds>(q1 - q0).count() / 1e09;
        const scalar Tl = duration_cast<nanoseconds>(r1 - r0).count() / 1e09;
        const scalar Tx = duration_cast<nanoseconds>(s1 - s0).count() / 1e09;
        const scalar Valpha = sum(alpha * mesh.V()).value();
        const scalar Ev = mag(Vs - Valpha) / Vs; 

        Info<< "Volume from surface mesh = " << Vs << endl;
        Info<< "Volume from volume fraction = " << Valpha << endl;
        Info<< "Volume error = " << Ev << endl; 
        Info<< "Initialization time = " << Ti << endl;
        Info<< "Calculation time = " << Te << endl;
        Info<< "Octree search time = " << To << endl;
        Info<< "Laplace solution time = " << Tl << endl;
        Info<< "Intersection time = " << Tx << endl;


        errorFile << tri.size() << "," 
            << mesh.nCells() << ","
            << Ev << ","  
            << Te << "\n";  
            
        // Move the tri surface back to the original position. 
        triPoints -= randomDisplacement;
        triLocalPoints -= randomDisplacement;

        runTime++; 
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
