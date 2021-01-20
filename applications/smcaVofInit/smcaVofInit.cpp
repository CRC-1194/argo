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
    poFoamVofInit

Description
    Initialize a volume fraction field from a triangulated surface
    using the polynomial approximation from signed distances as described in
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
#include "pointFields.H"
#include "triSurface.H"
#include "OFstream.H"

#include "polynomialVofInitialization.hpp"
#include "geomTriSurfaceTools.hpp"

#include <chrono>
#include <iomanip>

using namespace Foam::PolynomialVof;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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

    // Common options
    word fieldName = initDict.getOrDefault<word>("fieldName", "alpha.water");
    label refinementLevel = initDict.getOrDefault<label>("refinementLevel", -1);
    scalar narrowBandWidth = initDict.getOrDefault<scalar>("narrowBandWidth", 4.0);
    fileName surfaceFile = initDict.getOrDefault<fileName>("surfaceFile", args.path() + "/surface.stl");
    Switch invertVolumeFraction = initDict.getOrDefault<Switch>("invert", false);
    // Testing / debugging related options
    Switch writeFields = initDict.getOrDefault<Switch>("writeFields", false);
    Switch writeTets = initDict.getOrDefault<Switch>("writeTets", false);
    Switch checkVolume = initDict.getOrDefault<Switch>("checkVolume", false);

    // Comand line args
    args.readIfPresent<word>("fieldName", fieldName);
    args.readIfPresent<label>("refinementLevel", refinementLevel);
    args.readIfPresent<scalar>("narrowBandWidth", narrowBandWidth);
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    invertVolumeFraction = args.found("invert");
    writeFields = args.found("writeFields");
    writeTets = args.found("writeTets");
    checkVolume = args.found("checkVolume");

    // Print configuration
    Info << "Configuration:"
         << "\n\tfieldName: " << fieldName
         << "\n\trefinementLevel: " << refinementLevel
         << "\n\tsurfaceFile: " << surfaceFile
         << "\n\tinvert volume fraction: " << invertVolumeFraction
         << "\n\twriteFields: " << writeFields
         << "\n\twriteTets: " << writeTets
         << "\n\tcheckVolume: " << checkVolume
         << endl;

    // Initialization
    #include "createFields.hpp"
    triSurface surface{surfaceFile};
    alpha = dimensionedScalar("alpha", dimless, 0); 

    // Compute the volume fraction field.
    auto ctime0 = std::chrono::steady_clock::now();
    polynomialVofInitialization polyVofInit{mesh, surface, narrowBandWidth, IOobject::AUTO_WRITE, refinementLevel}; 
    polyVofInit.calcVolFraction(alpha, writeTets);
    auto ctime1 = std::chrono::steady_clock::now();
    auto calcTime = 
        std::chrono::duration_cast<std::chrono::microseconds>(ctime1-ctime0).count();

    if (invertVolumeFraction)
    {
        alpha = 1.0 - alpha;
    }

    alpha.write(); 

    // Begin testing and debugging
    if (writeFields)
    {
        polyVofInit.writeFields();
    }

    if (checkVolume)
    {
        // This is taken from surfaceCellVofInit (TT)
        if (Pstream::myProcNo() == 0) // Only compute on master rank in parallel. Thanks to TM.
        {
            scalar Vsurf = 0; 

            const auto& surfacePoints = surface.points(); 
            forAll(surface, triangleI) 
            {
                const auto& Sf = surface.Sf()[triangleI]; 
                const auto& triangle = surface[triangleI]; 
                Vsurf += dot(
                    -Sf,  // Surface normals are oriented into the phase. 
                    (surfacePoints[triangle[0]] + 
                    surfacePoints[triangle[1]] + 
                    surfacePoints[triangle[2]])
                );
            }

            Vsurf = 1./9.*mag(Vsurf); 

            scalar Valpha = gSum((mesh.V() * alpha)()); 
            scalar Evsurf = std::abs(Valpha - Vsurf) / Vsurf; 

            std::cout << std::setprecision(20) 
                << "Volume by volume fraction = " << Valpha << nl
                << "Volume of the surface mesh by divergence theorem (only closed surfaces!) = " << Vsurf << nl 
                << "Volume error from surface interval = " << Evsurf << nl;

            std::ofstream errorFile; 
            errorFile.open("smcaVofInit.csv"); 
            errorFile << "N_CELLS,"
                << "N_TRIANGLES,"
                << "VOLUME_FROM_VOLUME_FRACTION,"
                << "VOLUME_FROM_SURFACE_INTEGRAL,"
                << "VOLUME_ERROR_FROM_SURFACE_INTEGRAL,"
                << "CPU_TIME_MICROSECONDS," 
                << "MAX_REFINEMENT_LEVEL" << "\n"
                << mesh.nCells() << "," 
                << surface.size() << "," 
                << std::setprecision(20) 
                << Valpha << "," 
                << Vsurf << "," 
                << Evsurf << ","
                << calcTime << ","
                << polyVofInit.maxRefinementLevel() << "\n";
        }
    }

    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
