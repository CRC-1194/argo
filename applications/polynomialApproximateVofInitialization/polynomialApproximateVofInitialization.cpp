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

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label sqrDistFactor = args.optionLookupOrDefault<scalar>("sqrDistFactor", 3); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "surfaceCellMeshIntersection.csv"); 

    fileName triFile = args.path() + "/surface.stl";
    if (args.optionFound("surfaceFile"))
    {
        triFile = args.optionRead<fileName>("surfaceFile");
    }

    #include "createFields.hpp"

    triSurface surface{triFile};

    polynomialVofInitialization polyVofInit(mesh, surface, sqrDistFactor, IOobject::AUTO_WRITE, -1); 

    // TODO: Remove fix normals, this needs to be done externally, by picking an outside point 
    // and using surfaceOrient or another program. The STL must have consistent normals. TM. 
    const Switch fixNormals = args.optionFound("fixNormals");

    if (fixNormals)
    {
        Foam::GeometricalTransport::orientNormalsInward(surface);
    }

    // Compute the search distances once: the mesh is not moving, nor is it
    // topologically changed. 
    polyVofInit.calcSqrSearchDist(); 

    // Zero the volume fraction and signed distance for each test.
    alpha = dimensionedScalar("alpha", dimless, 0); 
    signedDist = dimensionedScalar("signedDist", dimLength,0);

    // Compute the signed distance field based on the surface octree.
    polyVofInit.calcSignedDist(); 

    // Compute the volume fraction field.
    polyVofInit.calcVolFraction(alpha); 

    alpha.write(); 

    // FIXME: Make a runtime option -testing out of this. 
//#ifdef TESTING
    polyVofInit.writeFields(); 
//#endif

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

    return 0;
}
// ************************************************************************* //
