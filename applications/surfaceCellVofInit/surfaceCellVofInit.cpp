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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurface.H"
#include "OFstream.H"

#include "geomSurfaceCellMeshIntersection.hpp"
#include "geomTriSurfaceTools.hpp"

using namespace GeometricalTransport;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label sqrDistFactor = args.optionLookupOrDefault<scalar>("sqrDistFactor", 2); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "surfaceCellMeshIntersection.csv"); 

    fileName triFile = args.path() + "/meshed-surface.vtk";
    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    triSurface triSurf(triFile);

    geomSurfaceCellMeshIntersection meshIntersection(mesh, triSurf, sqrDistFactor); 

    // Volume fraction field
    volScalarField alpha
    (
        IOobject
        (
            fieldName, 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh, 
        dimensionedScalar(fieldName, dimless, 0)
    );

    // Compute the volume fraction field.
    meshIntersection.calcVolFraction(alpha); 

    alpha.write(); 

    // FIXME: Make a runtime option -testing out of this. 
    meshIntersection.writeFields(); 

    // TODO: Make a runtime option out of this. 
    //if (checkVolume)
    //{
        //Info << "Cannot check volume 
        //const scalar Vs = starSurfaceVolume(meshIntersection.surface(), mesh.solutionD()); 
        //const scalar Valpha = sum(alpha * mesh.V()).value();
        //const scalar Ev = mag(Vs - Valpha) / Vs; 

        //Info<< "Volume from surface mesh = " << Vs << endl;
        //Info<< "Volume from volume fraction = " << Valpha << endl;
        //Info<< "Volume error = " << Ev << endl; 
    //}
        
    Info<< "End" << endl;
}


// ************************************************************************* //
