/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Tomislav Maric, TU Darmstadt
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
    foamTestIsoSurface 

Description

    Iso-surfaces 
    
        1. isoSurfaceTopo
        
    Test  

        1. FaceNormals
        2. PointNormals

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "error.H"
#include "fvCFD.H"

#include "isoSurfaceTopo.H"
#include "vtkSurfaceWriter.H"
#include "volPointInterpolation.H"
#include <IOobject.H>
#include <cstdlib>
#include <iomanip>

using namespace Foam::surfaceWriters;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "field",
        "string",
        "Name of the volume (cell-centered) field whose curvature is approximated."
    );
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word fieldName = args.get<word>("field");
    #include "createFields.H"

    volPointInterpolation vpInterp(mesh); 

    tmp<pointScalarField> psiPointTmp = vpInterp.interpolate(psi);
    pointScalarField& psiPoint = psiPointTmp.ref();
    psiPoint.rename(fieldName + ".point");
    psiPoint.write();

    isoSurfaceParams isoParams
    (
        isoSurfaceParams::algorithmType::ALGO_DEFAULT,
        isoSurfaceParams::filterType::CELL
    );
    isoParams.snap(false);
    isoParams.mergeTol(1e-03);
    isoSurfaceTopo isoTopo(
        mesh, 
        psi, 
        psiPoint, 
        0.5, 
        isoParams
    );
    Info << "isoSurface algorithm = " << int(isoTopo.algorithm()) << endl;
    Info << "isoSurface filter = " << int(isoTopo.filter()) << endl;
    Info << "isoSurface snap = " << Switch(isoTopo.snap()) << endl;
    Info << "isoSurface merge tolerance = " << isoTopo.mergeTol() << endl;
    std::stringstream ss; 
    ss << "isoTopo-" << std::setw(8) << std::setfill('0') 
        << Pstream::myProcNo() << ".vtk";
    isoTopo.write(ss.str());

    const auto& points = isoTopo.points();
    const auto& pointNormals = isoTopo.pointNormals();
    scalarField pointNormalErrors (pointNormals.size(),0);

    const auto& faceCenters = isoTopo.Cf(); 
    const auto& faceNormals = isoTopo.faceNormals(); 
    scalarField faceNormalErrors (faceNormals.size(),1);

    const vector origin = setAlphaFields.get<vector>("origin");
    
    forAll(points, pointI)
    {
        vector exactNormal = origin - points[pointI]; 
        exactNormal /= mag(exactNormal);
        pointNormalErrors[pointI] = mag(pointNormals[pointI] - exactNormal);
    }

    forAll(isoTopo, faceI)
    {
        vector exactNormal = origin - faceCenters[faceI];
        exactNormal /= mag(exactNormal);
        faceNormalErrors[faceI] = mag(faceNormals[faceI] - exactNormal);
    }
    
    // Write fields
    vtkWriter writer(isoTopo.points(), isoTopo, "isoSurfaceTopo");
    writer.nFields(4);

    // 1: point normals and errors 
    writer.isPointData(true);
    writer.write("pointNormal", pointNormals);
    writer.write("pointNormalError", pointNormalErrors);

    // 2: face normals and errors
    writer.isPointData(false);
    writer.write("faceNormal", faceNormals);
    writer.write("faceNormalError", faceNormalErrors);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //