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
    foamIsoSurfaceNearestDist 

Description
    Compute the nearest signed-distance from the Finite Volume mesh and and 
    iso-surface (triangulated surface mesh). 

Authors
    Tomislav Maric, MMA, TU Darmstadt, maric@mma.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "error.H"
#include "fvCFD.H"

#include "isoSurfaceTopo.H"
#include "fvcAverage.H"
#include "triSurfaceMesh.H"
#include "volPointInterpolation.H"
#include <cstdlib>
#include <iomanip>

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

    // Compute an isoSurface from the field fieldName

    volPointInterpolation vpInterp(mesh); 
    tmp<pointScalarField> psiPointTmp = vpInterp.interpolate(psi);
    pointScalarField& psiPoint = psiPointTmp.ref();
    psiPoint.rename(fieldName + ".point");
    psiPoint.write();

    isoSurfaceTopo isoTopo(
        mesh, 
        psi, 
        psiPoint, 
        0.5
    );

    std::stringstream ss; 
    ss << "isoTopo-" << std::setw(8) << std::setfill('0') 
        << Pstream::myProcNo() << ".vtk";
    isoTopo.write(ss.str());

    // Compute the nearest distance to the isoSurface
    // - Compute the squared searchRadius 
    sqrSearchRadius = Foam::sqr(
        2*fvc::average(
            Foam::pow(mesh.deltaCoeffs(), -1)
        ) // search radius = 2 * \Delta x 
    );
    
    // - The map Cell -> Nearest Triangle, Triangle in iso-surface.  
    DynamicList<pointIndexHit> cellTriangleNearest;
    const auto& cellCenters = mesh.C();
    
    // TODO: WARNING, from this point on, we have triangulated isoTopo and have 
    // multiple triangles per cell, meaning multiple triangle normals per cell. 
    // The transfer of isoTopo normals to interface cells must happen before. TM. 
    isoTopo.triangulate();
    
    // Construct triSurface for octree queries. 
    // triSurface (const List< labelledTri > &triangles, const pointField &pts)
    List<labelledTri> isoTopoTriangles(isoTopo.size(), labelledTri(-1,-1,-1,0));
    forAll(isoTopoTriangles, triI)
    {
        // After isoTopo.triangulate() isoTopo faces are triangles. 
        // triSurface needs labelledTri as a triangle data type.
        for(char i = 0; i < 3; ++i)
            isoTopoTriangles[triI][i] = isoTopo[triI][i];
    }
    
    // Construct octree
    
    triSurface isoTriSurface(isoTopoTriangles, isoTopo.points());
    isoTriSurface.write("isoTopoTriSurface.vtk");
    
    // triSurfaceMesh is constructed for octree queries. 
    // TODO: figure out how to use indexedOctree for the search query using
    // isoSurfaceTopo directly, if this is a computational bottlenechk. Profile. TM.
    triSurfaceMesh triMesh(
        IOobject(
            "triSurfaceMesh",
            "triSurfaceMesh",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        isoTriSurface 
    );

    triMesh.findNearest(
        cellCenters, // Find distance between iso-surface and cell centers.
        sqrSearchRadius, // squared search radius
        cellTriangleNearest 
    );

    // TODO: check if unit normals are available already. @tmaric
    const auto& faceNormals = isoTopo.faceNormals();
    
    // Calculate the distance to the cell center as the distance to the 
    // nearest front triangle. 
    forAll(cellTriangleNearest, cellI)
    {
        const pointIndexHit& h = cellTriangleNearest[cellI];

        if (h.hit())
        {
            const auto& normalVector = faceNormals[h.index()]; 
            auto distanceVector = cellCenters[cellI] - h.hitPoint(); 
            auto distSign = sign(normalVector & distanceVector); 
            sigDist[cellI] = distSign * mag(distanceVector);  
        }
    }
    
    sigDist.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
