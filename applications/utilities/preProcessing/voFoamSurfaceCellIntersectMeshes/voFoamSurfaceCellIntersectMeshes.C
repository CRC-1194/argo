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
#include "geomSurfaceCellMeshIntersection.H"
#include "geomTriSurfaceTools.H"
#include "OFstream.H"

using namespace GeometricalTransport;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "createOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const label sqrDistFactor = args.optionLookupOrDefault<scalar>("sqrDistFactor", 3); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFileName", "surfaceCellMeshIntersection.csv"); 

    fileName triFile = args.path() + "/surface.stl";
    if (args.optionFound("surfaceFile"))
        triFile = args.optionRead<fileName>("surfaceFile");

    #include "createFields.H"

    geomSurfaceCellMeshIntersection meshIntersection(mesh, triFile, sqrDistFactor); 

    const Switch fixNormals = args.optionFound("fixNormals");

    if (fixNormals)
        orientNormalsInward(meshIntersection.surfaceRef()); 

    // Open the error file for measurement output..
    OFstream errorFile(dataFileName); 
    // Nt : number of cells in the tool mesh 
    // Nb : number of cells in the tool mesh 
    // Ev : volume conservation error: 
    // |tool mesh volume from volume fraction - tool mesh volume | / tool mesh volume. 
    // Te : execution time of the CCI mesh intersection operation.
    errorFile << "Nt,Nb,Ev,Ti,Te,Tx\n"; 

    // Compute the search distances once: the mesh is not moving, nor is it
    // topologically changed. 
    meshIntersection.calcSqrSearchDist(); 

    // Zero the volume fraction and signed distance for each test.
    alpha = dimensionedScalar("alpha", dimless, 0); 
    signedDist = dimensionedScalar("signedDist", dimLength,0);

    // Compute the signed distance field based on the surface octree.
    meshIntersection.calcSignedDist(); 

    // Compute the volume fraction field.
    meshIntersection.calcVolFraction(alpha); 

    alpha.write(); 

#ifdef TESTING
    meshIntersection.writeFields(); 
#endif

    const scalar Vs = starSurfaceVolume(meshIntersection.surface(), mesh.solutionD()); 
    const scalar Valpha = sum(alpha * mesh.V()).value();
    const scalar Ev = mag(Vs - Valpha) / Vs; 

    Info<< "Volume from surface mesh = " << Vs << endl;
    Info<< "Volume from volume fraction = " << Valpha << endl;
    Info<< "Volume error = " << Ev << endl; 

    errorFile << meshIntersection.surface().size() << "," 
        << mesh.nCells() << ","
        << Ev << nl; 
        
    Info<< "End" << endl;
}


// ************************************************************************* //
