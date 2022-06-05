/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "pandoraDivNormalCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "error.H"

#include "processorFvPatchField.H"
#include "zoneDistribute.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraDivNormalCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraDivNormalCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraDivNormalCurvature::pandoraDivNormalCurvature
(
    const fvMesh& mesh, 
    const dictionary& dict
)
    :
        pandoraCurvature(mesh, dict), 
        fieldName_(curvatureDict_.get<word>("normalField")),
        nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 3)), 
        nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 3)), 
        averagedNormals_ 
        (
            IOobject
            (
                "averagedNormals", 
                mesh.time().timeName(), 
                mesh,
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ),
            mesh, 
            dimensionedVector("averagedNormals", dimless, vector(0,0,0))
        )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraDivNormalCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (meshDb.found(fieldName_))
    {
        const volVectorField& interfaceNormals = 
            mesh().lookupObject<volVectorField>(fieldName_);
        const volVectorField& interfaceCentres = 
            mesh().lookupObject<volVectorField>("interfaceCentre.dispersed");
        const volScalarField& isInterfaceCell = 
            mesh().lookupObject<volScalarField>("isInterfaceCell");
        const volScalarField& RDF = 
            mesh().lookupObject<volScalarField>("RDF");
        const volScalarField& psi = 
            mesh().lookupObject<volScalarField>("alpha.dispersed");

vector sphereCentre(0.005, 0.005, 0.005); // Sphere centre
scalar radius = 0.002; // Sphere radius

        averagedNormals_ == interfaceNormals /
        (
            mag(interfaceNormals) + 
            dimensionedScalar(
                "SMALL", interfaceNormals.dimensions(), SMALL
            )
        );
        volScalarField rdf = RDF;

        volScalarField markers = psi;
        forAll (markers, cellI)
        {
            markers[cellI] = -1;

            if (mag(averagedNormals_[cellI]) != 0)
            {
                markers[cellI] = 0;
            }
        }

        const auto& owner = mesh().owner();
        const auto& neibr = mesh().neighbour();

        zoneDistribute& distribute = zoneDistribute::New(mesh());

        boolList surfCells(mesh().nCells(), false);
        forAll (surfCells, zi)
        {
            if (markers[zi] == 0)
            {
                surfCells[zi] = true;
            }
        }
        distribute.setUpCommforZone(surfCells, true);

        for (label i = 0; i < nPropagate_; ++i)
        {
            forAll (neibr, fid)
            {
                if (markers[owner[fid]] == i && markers[neibr[fid]] == -1)
                {
                    markers[neibr[fid]] = i + 1;
                }
                if (markers[neibr[fid]] == i && markers[owner[fid]] == -1)
                {
                    markers[owner[fid]] = i + 1;
                }
            }

            auto& meshBoundary = markers.boundaryFieldRef();
            for (auto& mb : meshBoundary)
            {
                if (isA<processorFvPatch>(mb.patch()))
                {
                    const auto& faceToCell = mb.patch().faceCells();
                    const auto& neibrValue = mb.patchNeighbourField().cref();
                    
                    forAll(mb, j)
                    {
                        if (neibrValue[j] == i && markers[faceToCell[j]] == -1)
                        {
                            mb[j] = i + 1;
                        }
                    }
                }
            }

            boolList zone(mesh().nCells(), false);
            forAll (zone, zi)
            {
                if (markers[zi] == i + 1)
                {
                    zone[zi] = true;
                }
            }
            distribute.setUpCommforZone(zone, true);
            Map<vector> mapCentres = 
                distribute.getDatafromOtherProc(zone, mesh().C());
            Map<vector> mapNormals = 
                distribute.getDatafromOtherProc(zone, averagedNormals_);
            Map<scalar> mapMarkers = 
                distribute.getDatafromOtherProc(zone, markers);

            const labelListList& stencil = distribute.getStencil();

            volVectorField avgNormTmp = averagedNormals_;

            forAll (markers, cellI)
            {
                if (markers[cellI] != i + 1) continue;

                avgNormTmp[cellI] = Zero;
                const point& p = mesh().C()[cellI];

                label count = 0;
                for (const label gblIdx : stencil[cellI])
                {
                    scalar marker = distribute.getValue(markers, mapMarkers, gblIdx);
                    if (marker == -1 || marker == i + 1) continue;

                    vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);

                        vector centre = distribute.getValue(mesh().C(), mapCentres, gblIdx);
                        vector dist = centre - p;
                        vector distToSurf = dist & n / mag(n) * n;
                        vector verticalDist = dist - distToSurf;
                        avgNormTmp[cellI] += n / max(mag(verticalDist), SMALL);
                        count++;
                    }
                }
            }

            averagedNormals_ = avgNormTmp / 
            (
                mag(avgNormTmp) + 
                dimensionedScalar(
                    "SMALL", avgNormTmp.dimensions(), SMALL
                )
            );

            averagedNormals_.correctBoundaryConditions();
        }

        for (label i = 0; i < nAverage_; i++)
        {
            averagedNormals_ = fvc::average(averagedNormals_);
            averagedNormals_ /=
            (
                mag(averagedNormals_) + 
                dimensionedScalar(
                    "SMALL", averagedNormals_.dimensions(), SMALL
                )
            );

            averagedNormals_.correctBoundaryConditions();
        }

markers.rename("cellMarkers");
if (mesh().time().writeTime())
    markers.write();

scalar absError = 0;
scalar l1 = 0;
scalar l2 = 0;
scalar lInf = 0;
label count = 0;
volVectorField exactNormals = averagedNormals_;
volVectorField calcNormals = averagedNormals_;
scalar delta_x = max(pow(mesh().deltaCoeffs(), -1)).value();
scalar minArea = 0.001 * delta_x * delta_x;
forAll(averagedNormals_, i)
{
    if (mag(averagedNormals_[i]) != 0)
    {
        vector n1 = sphereCentre - mesh().C()[i];
        n1 /= mag(n1);
        exactNormals[i] = n1;

        vector n2 = averagedNormals_[i];
        n2 /= mag(n2);
        calcNormals[i] = n2;

        //if (markers[i] != 0) continue;

        absError = mag(n1 - n2);
        l1 += absError / mag(n1);
        l2 += sqr(absError) / magSqr(n1);
        count++;
        if (absError > lInf)
            lInf = absError;
    }
}
reduce(l1, sumOp<scalar>());
reduce(l2, sumOp<scalar>());
reduce(count, sumOp<label>());
reduce(lInf, maxOp<scalar>());
Info<<"l1 = "<<l1/count<<nl;
Info<<"l2 = "<<sqrt(l2/count)<<nl;
Info<<"lInf = "<<lInf<<nl;
exactNormals.rename("normExct");
if (cellCurvature_.time().writeTime())
    exactNormals.write();
calcNormals.rename("normCalc");
if (cellCurvature_.time().writeTime())
    calcNormals.write();

vector sumNorm = Zero;
forAll(averagedNormals_, i)
{
    if 
    (
        averagedNormals_[i][0] > 0
     && averagedNormals_[i][1] > 0
     && averagedNormals_[i][2] > 0
    )
    {
        sumNorm += averagedNormals_[i];
    }
}
reduce(sumNorm, sumOp<vector>());
Info<<"sumNorm = "<<sumNorm<<nl;

        cellCurvature_ = -fvc::div(averagedNormals_);

volScalarField cellCurv = cellCurvature_;
cellCurv.rename("cellCurv");
if(mesh().time().writeTime())
    cellCurv.write();

scalar sumCurv1 = 0;
forAll(cellCurvature_, i)
    sumCurv1 += cellCurvature_[i];
reduce(sumCurv1, sumOp<scalar>());
Info<<"sumCurv1 = "<<sumCurv1<<nl;

        forAll (cellCurvature_, cellI)
        {
            if (markers[cellI] == -1 || markers[cellI] == 2)
            {
                cellCurvature_[cellI] = 0;
            }
            else
            {
                cellCurvature_[cellI] = 2.0 / (2.0 / (cellCurvature_[cellI] + SMALL) + rdf[cellI] * 1.0);
            }
        }
        cellCurvature_.correctBoundaryConditions();

scalar sumCurv2 = 0;
forAll(cellCurvature_, i)
    sumCurv2 += cellCurvature_[i];
reduce(sumCurv2, sumOp<scalar>());
Info<<"sumCurv2 = "<<sumCurv2<<nl;
    }
    else
    {
        FatalErrorInFunction
            << "pandoraDivNormalCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
	    << "Available registered fields are : \n" 
	    << mesh().names() 
            << abort(FatalError);
    }


    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
