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

#include "isoSurfaceCell.H"
#include "isoSurfacePoint.H"
#include "isoSurfaceTopo.H"
#include "fvcAverage.H"
#include "triSurfaceMesh.H"
#include "volPointInterpolation.H"
#include <cstdlib>
#include <iomanip>

#include "surfaceIteratorIso.H"
#include "zoneDistribute.H"
#include "tensor2D.H"

#define iter 2
#define rr 1

#include "functions.H"
#include "inverseDistanceInterpolate.H"
#include "interpolateFunctions.H"

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

    // - iso-surface field: can be signed-distance, or another field
    volScalarField psi
    (
        IOobject
        (
            fieldName, 
            runTime.timeName(), 
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField marker
    (
        IOobject
        (
            "marker", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("marker", dimless, -1)
    );

    volVectorField normExac
    (
        IOobject
        (
            "normExac", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("normExac", dimless, Zero)
    );

    volVectorField normCalc
    (
        IOobject
        (
            "normCalc", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("normCalc", dimless, Zero)
    );

    volVectorField interfaceNormals
    (
        IOobject
        (
            "interfaceNormals", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("interfaceNormals", dimless, Zero)
    );

    volVectorField interfaceCentres
    (
        IOobject
        (
            "interfaceCentres", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("interfaceCentres", dimless, Zero)
    );

    vector sphereCentre(0.005, 0.005, 0.005); // Sphere centre

    // Compute an isoSurface from the field fieldName
    volPointInterpolation vpInterp(mesh); 

    tmp<pointScalarField> psiPointTmp = vpInterp.interpolate(psi);
    pointScalarField& psiPoint = psiPointTmp.ref();
    psiPoint.rename(fieldName + ".point");
    psiPoint.write();

    isoSurfaceParams isoParams
    (
        isoSurfaceParams::algorithmType::ALGO_DEFAULT,
        isoSurfaceParams::filterType::DIAGCELL
    );
    //isoParams.snap(false);

    isoSurfaceTopo isoTopo
    (
        mesh, 
        psi, 
        psiPoint, 
        0.5,
        isoParams
    );

    const auto& triToCell = isoTopo.meshCells();
    //const auto& triNormals = isoTopo.Sf();
    //const auto& triCentres = isoTopo.Cf();

    forAll (triToCell, i)
    {
        label cellI = triToCell[i];
        marker[cellI] = 0;
        interfaceNormals[cellI] = sphereCentre - mesh.C()[cellI];
    }

    /*
    forAll (psi, i)
    {
        if (psi[i] > 1e-8 && psi[i] < 1 - 1e-8)
        //if (psi[i] > SMALL && psi[i] < 1 - SMALL)
        {
            marker[i] = 0;
            interfaceNormals[i] = sphereCentre - mesh.C()[i];
        }
    }
    */

    marker.correctBoundaryConditions();
    interfaceNormals.correctBoundaryConditions();

    volVectorField averagedNormals_ = interfaceNormals /
    (
        mag(interfaceNormals) + 
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );
    averagedNormals_.correctBoundaryConditions();

    const auto& own = mesh.owner();
    const auto& nei = mesh.neighbour();

    zoneDistribute& distribute = zoneDistribute::New(mesh);

    for (label i = 0; i < 2; ++i)
    {
        forAll (nei, fid)
        {
            if (marker[own[fid]] == i && marker[nei[fid]] == -1)
            {
                marker[nei[fid]] = i + 1;
            }
            if (marker[nei[fid]] == i && marker[own[fid]] == -1)
            {
                marker[own[fid]] = i + 1;
            }
        }
        marker.correctBoundaryConditions();

        auto& meshBoundary = marker.boundaryFieldRef();
        for (auto& mb : meshBoundary)
        {
            if (isA<processorFvPatch>(mb.patch()))
            {
                mb.initEvaluate();
                mb.evaluate();

                const auto& faceToCell = mb.patch().faceCells();
                const auto& neibrValue = mb.patchNeighbourField().cref();
                
                forAll(mb, j)
                {
                    if (neibrValue[j] == i && marker[faceToCell[j]] == -1)
                    {
                        marker[faceToCell[j]] = i + 1;
                    }
                }
            }
        }
        marker.correctBoundaryConditions();
    }

    boolList updateZone(mesh.nCells(), false);
    forAll (updateZone, uzi)
    {
        if (marker[uzi] == -1) continue;
        if (marker[uzi] ==  0) continue;
        updateZone[uzi] = true;
    }
    distribute.updateStencil(updateZone);

    for (label i = 0; i < 2; ++i)
    {
        boolList zone(mesh.nCells(), false);
        forAll (zone, zi)
        {
            if (marker[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }

        distribute.setUpCommforZone(zone, false);
        Map<vector> mapMC = 
            distribute.getDatafromOtherProc(zone, mesh.C());
        Map<vector> mapNormals = 
            distribute.getDatafromOtherProc(zone, averagedNormals_);
        Map<scalar> mapMarkers = 
            distribute.getDatafromOtherProc(zone, marker);

        const labelListList& stencil = distribute.getStencil();

        volVectorField avgNormTmp = averagedNormals_;

        forAll (marker, cellI)
        {
            if (marker[cellI] != i + 1) continue;

            avgNormTmp[cellI] = Zero;

            point p = mesh.C()[cellI];

            DynamicField<vector> centres;
            DynamicField<scalar> valuesX;
            DynamicField<scalar> valuesY;
            DynamicField<scalar> valuesZ;

            for (const label gblIdx : stencil[cellI])
            {
                scalar mk = distribute.getValue(marker, mapMarkers, gblIdx);
                if (mk == -1 || mk == i + 1) continue;

                vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);
                if (mag(n) != 0)
                {
                    n /= mag(n);

                    vector centre = distribute.getValue(mesh.C(), mapMC, gblIdx);

                    vector dist = centre - p;
                    vector distToSurf = dist & n / mag(n) * n;
                    vector verticalDist = dist - distToSurf;

                    vector cc = p - verticalDist;

                    centres.append(cc);
                    valuesX.append(n.x());
                    valuesY.append(n.y());
                    valuesZ.append(n.z());
                }
            }

            avgNormTmp[cellI][0] = interpolateSecondOrder(p, centres, valuesX, rr);
            avgNormTmp[cellI][1] = interpolateSecondOrder(p, centres, valuesY, rr);
            avgNormTmp[cellI][2] = interpolateSecondOrder(p, centres, valuesZ, rr);
        }
        avgNormTmp.correctBoundaryConditions();

        averagedNormals_ = avgNormTmp / 
        (
            mag(avgNormTmp) + 
            dimensionedScalar(
                "SMALL", avgNormTmp.dimensions(), SMALL
            )
        );
        averagedNormals_.correctBoundaryConditions();
    }

    scalar absError = 0;
    scalar l1 = 0;
    scalar l2 = 0;
    scalar lInf = 0;
    label count = 0;
    forAll(marker, cellI)
    {
        if (marker[cellI] == -1) continue;
        if (marker[cellI] == 0) continue;

        vector n1 = sphereCentre - mesh.C()[cellI];
        n1 /= mag(n1);
        normExac[cellI] = n1;

        vector n2 = averagedNormals_[cellI];
        n2 /= mag(n2);
        normCalc[cellI] = n2;

        if (marker[cellI] != 1) continue;

        absError = 1 - (n1 & n2);
        l1 += absError / mag(n1);
        l2 += Foam::sqr(absError) / magSqr(n1);
        count++;
        if (absError > lInf)
            lInf = absError;
    }
    reduce(l1, sumOp<scalar>());
    reduce(l2, sumOp<scalar>());
    reduce(count, sumOp<label>());
    reduce(lInf, maxOp<scalar>());

    Info<<"count = "<<count<<nl;
    if (count > 0)
    {
        Info<<"l1 = "<<l1/count*100<<nl;
        Info<<"l2 = "<<Foam::sqrt(l2/count)*100<<nl;
        Info<<"lInf = "<<lInf<<nl;
    }

    //marker.write();
    //normExac.write();
    //normCalc.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
