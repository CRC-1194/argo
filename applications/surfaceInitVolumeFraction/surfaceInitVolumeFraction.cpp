/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    surfaceInitVolumeFraction

Description
    Initialize a volume fraction field from a triangulated surface.
    
    Computation of volume fractions in interface cells is either performed by
    the SMCI or SMCA algorithm.

    This application implements the SMCI/A algorithm described in

    Reference
    \verbatim
        Tolle, T., Gründing, D., Bothe, D., & Marić, T. (2021).
        Computing volume fractions and signed distances from arbitrary surfaces
        on unstructured meshes.
        arXiv preprint arXiv:2101.08511.
    \endverbatim

\*---------------------------------------------------------------------------*/

// STD headers
#include <chrono>
#include <iomanip>

// OpenFOAM headers
#include "fvCFD.H"
#include "triSurface.H"

// Argo headers
#include "volumeFractionCalculator.hpp"

// TODO: fix usage of name space (TT)
using namespace Foam; //::TriSurfaceImmersion;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * Configuration * * * *
    // Precedence: commandline option > dictionary value > default

    // rmlater
    // Idea: read dictionary, add default values for parameters which have not been set
    //      and overwrite parameterswhich are set via command line arguments
    
    // Read from dictionary if present
    IOdictionary initDict
                 {
                    IOobject{
                       "vofInitDict",
                       "system",
                       mesh.time(),
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                    }
                 };

    // Common options
    auto fieldName = initDict.getOrDefault<word>("fieldName", "alpha.water");
    auto refinementLevel = initDict.getOrDefault<label>("refinementLevel", -1);
    auto narrowBandWidth = initDict.getOrDefault<scalar>("narrowBandWidth", 4.0);
    auto surfaceFile = initDict.getOrDefault<fileName>("surfaceFile", args.path() + "/surface.stl");
    auto invertVolumeFraction = initDict.getOrDefault<Switch>("invert", false);
    auto vofCalcType = initDict.getOrDefault<word>("type", "SMCI");
    // Testing / debugging related options
    auto writeAllFields = initDict.getOrDefault<Switch>("writeAllFields", false);
    auto writeGeometry = initDict.getOrDefault<Switch>("writeGeometry", false);
    auto checkVolume = initDict.getOrDefault<Switch>("checkVolume", false);

    // Comand line args
    args.readIfPresent<word>("fieldName", fieldName);
    args.readIfPresent<label>("refinementLevel", refinementLevel);
    args.readIfPresent<scalar>("narrowBandWidth", narrowBandWidth);
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    args.readIfPresent<word>("type", vofCalcType);
    invertVolumeFraction = args.found("invert");
    writeAllFields = args.found("writeAllFields");
    writeGeometry = args.found("writeGeometry");
    checkVolume = args.found("checkVolume");

    // Print configuration
    Info << "Configuration:"
         << "\n\tfieldName: " << fieldName
         << "\n\tsurfaceFile: " << surfaceFile
         << "\n\tvof calculation: " << vofCalcType
         << "\n\trefinementLevel: " << refinementLevel
         << "\n\tnarrowbandWidth: " << narrowBandWidth
         << "\n\tinvert volume fraction: " << invertVolumeFraction
         << "\n\twriteAllFields: " << writeAllFields
         << "\n\twriteGeometry: " << writeGeometry
         << "\n\tcheckVolume: " << checkVolume
         << endl;
    
    // Write config to dictionary passed to RTS selection
    initDict.add("fieldName", fieldName);
    initDict.add("refinementLevel", refinementLevel);
    initDict.add("narrowBandWidth", narrowBandWidth);
    initDict.add("surfaceFile", surfaceFile);
    initDict.add("invert", invertVolumeFraction);
    initDict.add("type", vofCalcType);
    initDict.add("writeAllFields", writeAllFields);
    initDict.add("writeGeometry", writeGeometry);
    initDict.add("checkVolume", checkVolume);

    // Initialization
    #include "createFields.hpp"
    triSurface surface{surfaceFile};
    alpha = dimensionedScalar("alpha", dimless, 0); 

    // Compute the volume fraction field.
    auto ctime0 = std::chrono::steady_clock::now();
    auto vofCalcPtr =
        TriSurfaceImmersion::volumeFractionCalculator::New(initDict, mesh, surface);
    vofCalcPtr->printTypeName();
    vofCalcPtr->calcVolumeFraction(alpha);
    auto ctime1 = std::chrono::steady_clock::now();
    auto calcTime = 
        std::chrono::duration_cast<std::chrono::microseconds>(ctime1-ctime0).count();

    if (invertVolumeFraction)
    {
        alpha = 1.0 - alpha;
    }

    alpha.write(); 

    // Begin testing and debugging
    if (writeAllFields)
    {
        //polyVofInit.writeFields();
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
                << calcTime << "," << "\n";
                //<< polyVofInit.maxRefinementLevel() << "\n";
        }
    }

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

//} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
