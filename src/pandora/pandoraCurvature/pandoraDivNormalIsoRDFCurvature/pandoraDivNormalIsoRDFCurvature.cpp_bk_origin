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

#include "pandoraDivNormalIsoRDFCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "error.H"

#include "processorFvPatchField.H"
#include "reconstructionSchemes.H"
#include "zoneDistribute.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"
#include "triSurfaceMesh.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraDivNormalIsoRDFCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraDivNormalIsoRDFCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraDivNormalIsoRDFCurvature::pandoraDivNormalIsoRDFCurvature
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
{
    Info<<"Selecting divNormal curvature"<<nl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraDivNormalIsoRDFCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (meshDb.found(fieldName_))
    {
        //const volVectorField& interfaceNormals = 
        //    mesh().lookupObject<volVectorField>(fieldName_);

        reconstructionSchemes* surf = 
            mesh().getObjectPtr<reconstructionSchemes>("reconstructionScheme");

        const volVectorField& interfaceNormals = surf->normal();

volVectorField in = interfaceNormals;
in.rename("interfaceNormal");
if(mesh().time().writeTime())
    in.write();

        // Look up alpha field. 
        const volScalarField& vf = mesh().lookupObject<volScalarField>(fieldName_);

        // Construct isoSurfaceTopo. 
        volPointInterpolation vpInterp(mesh());
        tmp<pointScalarField> pfTmp = vpInterp.interpolate(vf);
        pointScalarField& pf = pfTmp.ref();

        isoSurfaceTopo isoTopo(
            mesh(),
            vf,
            pf,
            0.5
        );
        const auto& triToCell = isoTopo.meshCells();

        averagedNormals_ = vector::zero;
        forAll (triToCell, i)
        {
            label cellI = triToCell[i];
            averagedNormals_[cellI] = interfaceNormals[cellI];
        }

        averagedNormals_ /=
        (
            mag(averagedNormals_) + 
            dimensionedScalar(
                "SMALL", averagedNormals_.dimensions(), SMALL
            )
        );
        averagedNormals_.correctBoundaryConditions();

        volScalarField markers
        (
            IOobject
            (
                "cellMarkers",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("cellMarkers", dimless, -1)
        );

        forAll (markers, cellI)
        {
            if (mag(averagedNormals_[cellI]) != 0)
            {
                markers[cellI] = 0;
            }
        }
        markers.correctBoundaryConditions();

        const auto& own = mesh().owner();
        const auto& nei = mesh().neighbour();

        zoneDistribute& distribute = zoneDistribute::New(mesh());

        for (label i = 0; i < nPropagate_; ++i)
        {
            forAll (nei, fid)
            {
                if (markers[own[fid]] == i && markers[nei[fid]] == -1)
                {
                    markers[nei[fid]] = i + 1;
                }
                if (markers[nei[fid]] == i && markers[own[fid]] == -1)
                {
                    markers[own[fid]] = i + 1;
                }
            }
            markers.correctBoundaryConditions();

            auto& meshBoundary = markers.boundaryFieldRef();
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
                        if (neibrValue[j] == i && markers[faceToCell[j]] == -1)
                        {
                            markers[faceToCell[j]] = i + 1;
                        }
                    }
                }
            }
            markers.correctBoundaryConditions();
        }

        boolList updateZone(mesh().nCells(), false);
        forAll (updateZone, uzi)
        {
            if (markers[uzi] == -1) continue;
            if (markers[uzi] ==  0) continue;
            updateZone[uzi] = true;
        }
        distribute.updateStencil(updateZone);

        for (label i = 0; i < nPropagate_; ++i)
        {
            boolList zone(mesh().nCells(), false);
            forAll (zone, zi)
            {
                if (markers[zi] == i + 1)
                {
                    zone[zi] = true;
                }
            }

            distribute.setUpCommforZone(zone, false);
            Map<vector> mapMC = 
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

                for (const label gblIdx : stencil[cellI])
                {
                    scalar marker = distribute.getValue(markers, mapMarkers, gblIdx);
                    if (marker == -1 || marker == i + 1) continue;

                    vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);

                        vector centre = distribute.getValue(mesh().C(), mapMC, gblIdx);
                        vector dist = centre - p;
                        vector distToSurf = dist & n / mag(n) * n;
                        vector verticalDist = dist - distToSurf;
                        avgNormTmp[cellI] += n / max(mag(verticalDist), SMALL);
                    }
                }
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

        //#include "error.hpp"

        cellCurvature_ == -fvc::div(averagedNormals_);

        const volScalarField& rdf = 
            mesh().lookupObject<volScalarField>("RDF");

        forAll (cellCurvature_, cellI)
        {
            if (markers[cellI] == -1 || markers[cellI] == nPropagate_)
                cellCurvature_[cellI] = 0;
            else
                cellCurvature_[cellI] = 2.0 / (2.0 / (cellCurvature_[cellI] + SMALL) + rdf[cellI] * 1.0);
        }

        cellCurvature_.correctBoundaryConditions();
    }
    else
    {
        FatalErrorInFunction
            << "pandoraDivNormalIsoRDFCurvature::cellCurvature \n"
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
