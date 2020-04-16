/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
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
    poFoamTestVofInit

Description
    Test the volume fraction initialization using the polynomial approximation from
    signed distances as described in
        "From level set to volume of fluid and back again at second‐order accuracy"
        Detrixhe and Aslam, 2016.
    Adaptive tetrahedral refinement is used to increase the accuracy of the approach.

Author
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

#include "geomIOhelperFunctions.hpp"
#include "geomTriSurfaceTools.hpp"
#include "polynomialVofInitialization.hpp"

// Time measurement
#include <chrono>
// Random surface mesh positioning
#include <random>

using namespace Foam::PolynomialVof;
using namespace Foam::GeometricalTransport;
using namespace std::chrono;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
scalar volume_surface_integral(const triSurface& triSurf)
{
    scalar Vsurf = 0; 

    const auto& triSurfPoints = triSurf.points(); 
    forAll(triSurf, triangleI) 
    {
        const auto& Sf = triSurf.Sf()[triangleI]; 
        const auto& triangle = triSurf[triangleI]; 
        Vsurf += dot(
            -Sf,  // Surface normals are oriented into the phase. 
            (triSurfPoints[triangle[0]] + 
            triSurfPoints[triangle[1]] + 
            triSurfPoints[triangle[2]])
        );
    }

    return Vsurf/9.; 
}

int main(int argc, char *argv[])
{
    // Timing point 
    high_resolution_clock::time_point p0 = high_resolution_clock::now();

    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // * * * * Configuration * * * *
    // Precedence: commandline option > dictionary value > default
    
    // Read from dictionary if present
    IOdictionary initDict
                 (
                    IOobject(
                       "vofInitDict",
                       "system",
                       mesh.time(),
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                    )
                 );

    word fieldName = initDict.getOrDefault<word>("fieldName", "alpha.water");
    label refinementLevel = initDict.getOrDefault<label>("refinementLevel", -1);
    fileName surfaceFile = initDict.getOrDefault<fileName>("surfaceFile", args.path() + "/surface.stl");
    word dataFileName = initDict.getOrDefault<word>("dataFile", "polynomialVofInitResults.csv");
    Switch writeFields = initDict.getOrDefault<Switch>("writeFields", false);
    Switch randomPlacement = initDict.getOrDefault<Switch>("randomPlacement", false);

    // Comand line args
    args.readIfPresent<word>("fieldName", fieldName);
    args.readIfPresent<label>("refinementLevel", refinementLevel);
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    args.readIfPresent<word>("dataFile", dataFileName);
    writeFields = args.found("writeFields");
    randomPlacement = args.found("random");

    // Print configuration
    Info << "Test configuration:"
         << "\n\tfieldName: " << fieldName
         << "\n\trefinementLevel: " << refinementLevel
         << "\n\tsurfaceFile: " << surfaceFile
         << "\n\tdataFile: " << dataFileName
         << "\n\twriteFields: " << writeFields
         << "\n\trandomPlacement: " << randomPlacement
         << endl;

    
    // Initialization
    triSurface surface{surfaceFile};
    #include "createFields.hpp"

    // Timing point 
    high_resolution_clock::time_point p1 = high_resolution_clock::now();

    // Position the tool mesh centroid at the base mesh centroid. 
    if (randomPlacement)
    {
        label nBasePoints = mesh.nPoints();  
        const vector baseCenter = sum(mesh.points()) / nBasePoints;  
        vector toolCenter = sum(surface.points()) / surface.nPoints(); 
        displaceSurface(surface, baseCenter - toolCenter); 
    }

    if (writeFields)
    {
        surface.write(appendSuffix("trisurface-initial", runTime.timeIndex()) + ".stl");
    }

    // Open the error file for measurement output..
    std::ofstream errorFile;
    errorFile.open(dataFileName); 
    errorFile << "N_CELLS,N_TRIANGLES,N_INTERFACE_CELLS," 
              << "VOLUME_FROM_VOLUME_FRACTION,"
              << "VOLUME_FROM_SURFACE_INTEGRAL,"
              << "VOLUME_ERROR_FROM_SURFACE_INTEGRAL,"
              << "CPU_TIME_MICROSECONDS," 
              << "MAX_REFINEMENT_LEVEL" << "\n";

    while(runTime.run())
    {
        // Computation
        
        // Zero the volume fraction for each test.
        alpha = dimensionedScalar("alpha", dimless, 0); 

        // Position the surface mesh randomly within the bounding box of 
        // the volume mesh.
        vector randomDisplacement{0,0,0};

        if (randomPlacement)
        {
            randomDisplacement = 
                placeSurfaceRandomlyInBox(surface, mesh.bounds(), mesh.solutionD()); 
        }

        if (writeFields)
        {
            surface.write(appendSuffix("trisurface", runTime.timeIndex()) + ".stl");
        }
        
        // Measurement point
        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        // Compute the volume fraction field.
        polynomialVofInitialization vofInit{mesh, surface, 3.0, IOobject::AUTO_WRITE, refinementLevel};
        vofInit.calcVolFraction(alpha);

        // Measurement point
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        // Evaluation
        if (writeFields)
        {
            vofInit.writeFields(); 
            alpha.write(); 
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        const scalar Ti = duration_cast<nanoseconds>(p1 - p0).count() / 1e03;
        const scalar Te = duration_cast<nanoseconds>(t1 - t0).count() / 1e03;

        scalar V_surf_int = volume_surface_integral(surface);
        const scalar V_alpha = sum(alpha * mesh.V()).value();
        const scalar Ev = mag(V_alpha - V_surf_int) / V_surf_int; 

        Info<< "Volume from surface integral = " << V_surf_int << endl;
        Info<< "Volume from volume fraction = "  << V_alpha << endl;
        Info<< "Volume error = "                 << Ev << endl; 
        Info<< "Initialization time = " << Ti << " microseconds"<< endl;
        Info<< "Calculation time = "    << Te << " microseconds" << endl;

        label N_interface_cells = 0; 

        forAll(alpha, cellI)
        {
            if ((alpha[cellI] > 0) && (alpha[cellI] < 1))
                ++N_interface_cells; 
        }

        errorFile << mesh.nCells() << ","
            << surface.size() << "," 
            << N_interface_cells << ","
            << std::setprecision(20)
            << V_alpha << ","
            << V_surf_int << ","
            << Ev << ","  
            << Ti+Te << "," 
            << vofInit.maxRefinementLevel()
            << std::endl;
            
        // Move the surface back to its original position. 
        displaceSurface(surface, -randomDisplacement);

        runTime++; 
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
