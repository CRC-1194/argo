/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Similar to createFields.H, this file initializes two meshes for use with
   the meshIntersection class. 

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

// TODO: Generalize the options more, so that all combinations of options can
// be applied. TM.

// Initialize the base time and mesh.

fileName baseCase = args.getOrDefault<fileName>("case", "./"); 
baseCase.clean();
const fileName baseRoot = baseCase.path();
baseCase = baseCase.name();

Info << "baseRoot = " << baseRoot << endl
    << "baseCase = " << baseCase << endl;

Time runTimeBase
(
    Time::controlDictName,
    baseRoot,
    baseCase 
);

Info << "Reading base mesh for time = " << runTimeBase.timeName() << nl;

Info<< "Creating base mesh\n" << endl;

fvMesh baseMesh
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTimeBase.timeName(),
        runTimeBase, 
        IOobject::MUST_READ, 
        IOobject::AUTO_WRITE
    )
);

// Initialize the tool time and mesh.

fileName toolCase = args.getOrDefault<fileName>("toolCase", "../toolCase"); 
toolCase.clean(); 
const fileName toolRoot = toolCase.path();
toolCase = toolCase.name(); 

Time runTimeTool
(
    Time::controlDictName,
    toolRoot,
    toolCase
);

Info << "toolRoot = " << toolRoot << endl
    << "toolCase = " << toolCase << endl;

Info << "Reading tool mesh for time = " << runTimeTool.timeName() << nl;

Info<< "Creating tool mesh\n" << endl;

fvMesh toolMesh
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTimeTool.timeName(),
        runTimeTool, 
        IOobject::MUST_READ, 
        IOobject::AUTO_WRITE
    )
);

// ************************************************************************* //
