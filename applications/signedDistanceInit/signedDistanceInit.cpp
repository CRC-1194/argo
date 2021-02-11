/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 AUTHOR,AFFILIATION
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
    signedDistanceCalc

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurface.H"

#include "insideOutsidePropagation.hpp"
#include "signedDistanceCalculator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
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
                       "signedDistanceInitDict",
                       "system",
                       mesh.time(),
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                    )
                 );

    // Read options from dictionary
    word fieldName = initDict.getOrDefault<word>("fieldName", "signedDistance");
    fileName surfaceFile = initDict.getOrDefault<fileName>("surfaceFile", args.path() + "/surface.stl");
    scalar searchDistanceFactor = initDict.getOrDefault<scalar>("searchDistanceFactor", 3.0);
    Switch propagateInsideOutside = initDict.getOrDefault<Switch>("propagateInsideOutside", true);
    Switch invertInsideOutside = initDict.getOrDefault<Switch>("invert", false);

    // Comand line args
    args.readIfPresent<word>("fieldName", fieldName);
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    args.readIfPresent<scalar>("searchDistanceFactor", searchDistanceFactor);
    propagateInsideOutside = args.found("propagateInsideOutside");
    invertInsideOutside = args.found("invert");

    // Print configuration
    Info << "Configuration:"
         << "\n\tfieldName: " << fieldName
         << "\n\tsurfaceFile: " << surfaceFile
         << "\n\tsearchDistanceFactor: " << searchDistanceFactor
         << "\n\tpropagateInsideOutside: " << propagateInsideOutside
         << "\n\tinvert inside outside: " << invertInsideOutside
         << endl;

    // Initialization
    #include "createFields.hpp"
    triSurface surface{surfaceFile};
    SigDistCalc::signedDistanceCalculator sig_dist_calc{surface};

    if (searchDistanceFactor <= 0.0)
    {
        searchDistanceFactor = mesh.bounds().mag();
    }
    signedDistance.primitiveFieldRef() = sig_dist_calc.signed_distance(mesh.C(), sqrSearchDist*searchDistanceFactor*searchDistanceFactor, 0.0);

    if (propagateInsideOutside)
    {
        signedDistance = SigDistCalc::insideOutsidePropagation{}.propagate_inside_outside(signedDistance);
    }
    
    if (invertInsideOutside)
    {
        signedDistance *= -1.0;
    }

    signedDistance.write();

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //