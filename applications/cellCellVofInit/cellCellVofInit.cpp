/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Tomislav Maric and TU Darmstadt 
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

Description
    Compute the volume fraction field in the base mesh by intersecting the base
    mesh with the tool mesh. 

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"

#include "geomMeshIntersection.hpp"

using namespace Foam;
using namespace GeometricalTransport;

#include <omp.h>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field." 
    ); 

    argList::addBoolOption
    (
        "writeGeometry", 
        "Write the cell/cell intersections." 
    ); 

    argList::addOption
    (
        "toolCase", 
        "pathName",
        "Path to the tool case. Default '../toolCase'" 
    ); 

    argList::addOption
    (
        "dataFile", 
        "geomIntersectMeshes.dat",
        "Name of file used to save the reported errors." 
    ); 

    argList args(argc, argv);

    #include "createMeshes.hpp"

    const word fieldName = args.getOrDefault<word>("fieldName", "alpha.water"); 
    const word dataFile = args.getOrDefault<word>("dataFile", "geomIntersectMeshes.dat"); 

    Switch writeGeometry = true;  
    if (! args.found("writeGeometry"))
        writeGeometry = false; 

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTimeBase.timeName(),
            baseMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        baseMesh, 
        dimensionedScalar ("zero", dimless, 0)
    );

    geomMeshIntersection intersect(baseMesh, toolMesh, writeGeometry);
    intersect.setVolFraction(alpha1);

    alpha1.write(); 
    intersect.report(Info, alpha1); 

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
