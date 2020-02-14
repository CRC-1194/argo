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
#include "OFstream.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"

#include "geomSurfaceCellMeshIntersection.hpp"
#include "geomTriSurfaceTools.hpp"

// Time measurement
#include <chrono>
// Random surface mesh positioning
#include <random>

using namespace GeometricalTransport;
using namespace std::chrono;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    /*
    // Timing point 
    high_resolution_clock::time_point p0 = high_resolution_clock::now();

    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label sqrDistFactor = args.optionLookupOrDefault<scalar>("sqrDistFactor", 3); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "smci.csv"); 

    fileName triFile = args.path() + "/surface.stl";
    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    // Triangle surface
    triSurface tri(triFile);

    #include "createFields.hpp"

    // Timing point 
    high_resolution_clock::time_point p1 = high_resolution_clock::now();

    geomSurfaceCellMeshIntersection meshIntersection(mesh, sqrDistFactor); 

    const Switch fixNormals = args.optionFound("fixNormals");

    if (fixNormals)
        orientNormalsInward(tri); 

    // Position the tool mesh centroid at the base mesh centroid. 
    label nBasePoints = mesh.nPoints();  
    Pstream::gather(nBasePoints, sumOp<label>()); 
    const vector baseCenter = gSum(mesh.points()) / nBasePoints;  
    vector toolCenter = sum(tri.points()) / tri.nPoints(); 
    displaceSurface(tri, baseCenter - toolCenter); 

#ifdef TESTING // Write the initial surface and test the centroid position.
    tri.write(appendSuffix("trisurface-initial", runTime.timeIndex()) + ".stl");
#endif

    // Open the error file for measurement output..
    OFstream errorFile(dataFileName); 
    // Nt : number of cells in the tool mesh 
    // Nb : number of cells in the base mesh 
    // Ev : volume conservation error: 
    // |tool mesh volume from volume fraction - tool mesh volume | / tool mesh volume. 
    // Ti : initialization time, loading the mesh and fields. 
    // Te : execution time.
    // Nx : total number of halfspace / cell intersections. 
    // Ni : number of interface cells, 
    // Nk : number of bulk cells.
    errorFile << "Nt,Nb,Ev,Ti,Te,Nx,Ni,Nk\n"; 

    // Compute the search distances once: the mesh is not moving, nor is it
    // topologically changed. 
    meshIntersection.calcSqrSearchDist(); 

    while(runTime.run())
    {
        // Zero the volume fraction and signed distance for each test.
        alpha = dimensionedScalar("alpha", dimless, 0); 
        signedDist = dimensionedScalar("signedDist", dimLength,0);

        // Position the surface mesh randomly within the bounding box of 
        // the volume mesh.
        vector randomDisplacement = 
            placeSurfaceRandomlyInBox(tri, mesh.bounds(), mesh.solutionD()); 

#ifdef TESTING
        tri.write(appendSuffix("trisurface", runTime.timeIndex()) + ".stl");
#endif
        // Measurement point
        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        // Compute the volume fraction field.
        meshIntersection.calcVolFraction(alpha, tri); 

        // Measurement point
        high_resolution_clock::time_point t1 = high_resolution_clock::now();


#ifdef TESTING
        meshIntersection.writeFields(); 
        alpha.write(); 
#endif
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        const scalar Ti = duration_cast<nanoseconds>(p1 - p0).count() / 1e09;
        const scalar Te = duration_cast<nanoseconds>(t1 - t0).count() / 1e09;

        const scalar Vs = starSurfaceVolume(tri, mesh.solutionD()); 
        const scalar Valpha = sum(alpha * mesh.V()).value();
        const scalar Ev = mag(Vs - Valpha) / Vs; 

        Info<< "Volume from surface mesh = " << Vs << endl;
        Info<< "Volume from volume fraction = " << Valpha << endl;
        Info<< "Volume error = " << Ev << endl; 
        Info<< "Initialization time = " << Ti << endl;
        Info<< "Calculation time = " << Te << endl;

        label Nb = 0; 
        label Ni = 0; 

        forAll(alpha, cellI)
        {
            if ((alpha[cellI] > 0) && (alpha[cellI] < 1))
                ++Ni; 
            else
                ++Nb;
        }

        errorFile << tri.size() << "," 
            << mesh.nCells() << ","
            << Ev << ","  
            << Ti << "," 
            << Te << "," 
            << meshIntersection.Nx() << ","
            << Ni << ","
            << Nb << endl;
            
        // Move the surface back to its original position. 
        displaceSurface(tri, -randomDisplacement);

        runTime++; 
    }

    Info<< "End" << endl;

    */
    return 0;
}


// ************************************************************************* //
