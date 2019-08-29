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

    word fieldName = initDict.getOrDefault<word>("fieldName", "alpha.water");
    label refinementLevel = initDict.getOrDefault<label>("refinementLevel", -1);
    fileName surfaceFile = initDict.getOrDefault<fileName>("surfaceFile", args.path() + "/surface.stl");
    Switch writeFields = initDict.getOrDefault<Switch>("writeFields", false);
    Switch invertVolumeFraction = initDict.getOrDefault<Switch>("invert", false);

    // Comand line args
    args.readIfPresent<word>("fieldName", fieldName);
    args.readIfPresent<label>("refinementLevel", refinementLevel);
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    args.readIfPresent<Switch>("writeFields", writeFields);
    args.readIfPresent<Switch>("invert", invertVolumeFraction);

    // Print configuration
    Info << "Configuration:"
         << "\n\tfieldName: " << fieldName
         << "\n\trefinementLevel: " << refinementLevel
         << "\n\tsurfaceFile: " << surfaceFile
         << "\n\twriteFields: " << writeFields
         << "\n\tinvert volume fraction: " << invertVolumeFraction
         << endl;

    // Initialization
    #include "createFields.hpp"
    triSurface surface{surfaceFile};
    alpha = dimensionedScalar("alpha", dimless, 0); 

    // Compute the volume fraction field.
    polynomialVofInitialization polyVofInit{mesh, surface, 3.0, IOobject::AUTO_WRITE, refinementLevel}; 
    polyVofInit.calcVolFraction(alpha);

    if (invertVolumeFraction)
    {
        alpha = 1.0 - alpha;
    }

    if (writeFields)
    {
        polyVofInit.writeFields();
    }

    scalar Valpha = sum(alpha * mesh.V()).value();

    Info << "Overall volume given by " << fieldName
         << ": " << Valpha << endl;

    alpha.write(); 
        
    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
