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
#include <iomanip>
#include <chrono>

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

    const word fieldName = 
        args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 

    const label sqrDistFactor = 
        args.optionLookupOrDefault<scalar>("sqrDistFactor", 2); 

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
    auto ctime0 = std::chrono::steady_clock::now();
    meshIntersection.calcVolFraction(alpha); 
    auto ctime1 = std::chrono::steady_clock::now();
    auto calcTime = 
        std::chrono::duration_cast<std::chrono::microseconds>(ctime1-ctime0).count();
    // Write the volume fraction field.
    alpha.write(); 

    if (args.optionFound("writeAllFields"))
        meshIntersection.writeFields(); 


    if (args.optionFound("checkVolume"))
    {
        if (Pstream::myProcNo() == 0) // Only compute on master rank in parallel. TM.
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

            Vsurf *= 1./9.; 

            scalar Valpha = gSum((mesh.V() * alpha)()); 
            scalar Evsurf = std::abs(Valpha - Vsurf) / Vsurf; 

            std::cout << std::setprecision(20) 
                << "Volume by volume fraction = " << Valpha << nl
                << "Volume of the surface mesh by divergence theorem (only closed surfaces!) = " << Vsurf << nl 
                << "Volume error from surface interval = " << Evsurf << nl;

            std::ofstream errorFile; 
            errorFile.open("surfaceCellVofInit.csv"); 
            errorFile << "N_CELLS,"
                << "N_TRIANGLES,"
                << "VOLUME_FROM_VOLUME_FRACTION,"
                << "VOLUME_FROM_SURFACE_INTEGRAL,"
                << "VOLUME_ERROR_FROM_SURFACE_INTEGRAL,"
                << "CPU_TIME_MICROSECONDSi," 
                << "N_TRIANGLES_PER_CELL" << "\n"
                << mesh.nCells() << "," 
                << triSurf.size() << "," 
                << std::setprecision(20) 
                << Valpha << "," 
                << Vsurf << "," 
                << Evsurf << ","
                << calcTime << ","
                << meshIntersection.nTrianglesPerCell() << "\n";
        }
    }

    Info<< "End" << endl;
}


// ************************************************************************* //
