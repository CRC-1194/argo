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
    triangulations with a set of nearest triangles that intersect the cell. The
    volume of the intersection divided by the volume of the cell is the volume
    fraction of the surface in the cell.  

Author
    Tomislav Maric maric@mma.tu-darmstadt.de, tomislav@sourceflux.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Geometry.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "OFstream.H"
#include <set>
#include <omp.h>

using namespace GeometricalTransport;

#include "surfaceMeshIntersection.H"

typedef TriangulationIntersection<tetrahedronVector,pointVectorVector> triangulationIntersection;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#ifdef TESTING
    double t0 = omp_get_wtime(); 
#endif

    #include "createOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label distanceFactor = args.optionLookupOrDefault<scalar>("distanceFactor", 1.5); 
    const label radiusFactor = args.optionLookupOrDefault<scalar>("radiusFactor", 1.71); 

    fileName triFile = args.path() + "/surface.stl";

    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    // Triangle surface
    triSurface tri(triFile);

    // Triangle surface + octree
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
        mesh.deltaCoeffs();
        mesh.cells();
        mesh.faces();
        mesh.Sf();
        mesh.C();
        mesh.V(); 
        mesh.faceOwner();
    }

#ifdef TESTING
    double t1 = omp_get_wtime(); 
#endif 

    const auto& triPoints = tri.points();
    auto& triNormals = const_cast<pointField&>(tri.faceNormals());
    auto& triangles = static_cast<List<labelledTri>& >(tri);

    const Switch fixNormals = args.optionFound("fixNormals");
    const Switch checkVolume  = args.optionFound("checkVolume");
    const Switch writeCutCells = args.optionFound("writeCutCells");

    vector surfaceMeshCentroid (GREAT,GREAT,GREAT); 
           
    if (fixNormals || checkVolume)
    {
        surfaceMeshCentroid = Foam::sum(triPoints); 
        surfaceMeshCentroid /= tri.nPoints(); 
    }

    if(fixNormals)
    {
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

    scalar vol = 0; 
    // Compute the surfaceMeshCentroid volume from the surface mesh. 
    if (checkVolume)
    {
        forAll(triangles, tI)
            vol +=  -1.0 / 6.0 * (triPoints[triangles[tI][0]] - surfaceMeshCentroid) & 
             ((triPoints[triangles[tI][1]] - surfaceMeshCentroid) ^ 
              (triPoints[triangles[tI][2]] - surfaceMeshCentroid)); 
    } 

#ifdef TESTING
    const scalar t2 = omp_get_wtime(); 
#endif

    // Use octree search to find a nearest triangle for each cell.
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

    
#ifdef TESTING
    const scalar t3 = omp_get_wtime(); 
#endif

    // Compute the marker field. 
    const auto& C = mesh.C();
    const auto& V = mesh.V(); 

    vtk_polydata_stream cutCellStream(runTime.path() + "/cutCells.vtk");

    const auto& octree = triMesh.tree();

    // Distance calculation and setting of alpha1 must precede geometrical 
    // intersections because there may exist cells with positive signed 
    // distance that are intersected and whose alpha1 value must then 
    // be subsequently corrected. TM.
    #pragma omp parallel for schedule (dynamic) 
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

#ifdef TESTING
    volScalarField signedDist0 ("signedDist0", signedDist); 
    signedDist0.write();
#endif

    // Solve a laplace equation to propagate the sign.
    dimensionedScalar lambda ("lambda", sqr(dimLength) * pow(dimTime,-1), 1);
    fvScalarMatrix distEqn
    (
        fvm::ddt(signedDist) - fvm::laplacian(lambda, signedDist)
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

#ifdef TESTING
    const scalar t4 = omp_get_wtime(); 
#endif

    #pragma omp parallel for schedule (dynamic) 
    forAll(cellNearestTriangle, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle[cellI];

        if (cellHit.hit()) 
        {
            labelList cellTriangles = octree.findSphere(C[cellI],sqrCellRadius[cellI]);

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

                // Calculate volume fraction value that is upper bounded by 1. 
                alpha[cellI] = min(1, volume(cellIntersection) / V[cellI]); 

                // If the value is an interface value set it, otherwise use 1 that was
                // set by the distance field calculation.

                if (writeCutCells)
                {
                    #pragma omp critical
                    cutCellStream << cellIntersection;
                }
            }
        }
    }

#ifdef TESTING
    const scalar t5 = omp_get_wtime(); 
#endif

#ifdef TESTING
    signedDist.write();
#endif

    // Write the volume fraction field.
    alpha.write(); 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //OFstream(args.path() + "/voFoamSetSurfaceFraction.dat");

#ifdef TESTING
    const scalar init = t1 - t0; 
    const scalar volNormals = t2 - t1;  
    const scalar octreeNearest = t3 - t2; 
    const scalar distanceSolution = t4 - t3; 
    const scalar intersections = t5 - t4; 

    Info<< "\nData initialization time = " << init << endl;
    Info<< "Volume and normal checks = " << volNormals << endl;
    Info<< "Octree nearest distance calculation time = " << octreeNearest << endl;
    Info<< "Intersection time = " << intersections << endl;
    Info<< "Distance field solution time = " << distanceSolution << endl;
    Info<< "Total calculation time = " << t5 - t0 << endl << endl;

#endif

    if (checkVolume)
    {
        Info<< "Volume from surface mesh = " << vol << endl;
        const scalar alphaVol = sum(alpha * mesh.V()).value();
        Info<< "Volume from volume fraction = " << alphaVol << endl;
        Info<< "Volume error = " << (vol - alphaVol) / vol << endl;
    }

    Info<< "End" << endl;

#ifdef TESTING
    // TODO: Output the volume error and timing data into a file. 
#endif

    return 0;
}


// ************************************************************************* //
