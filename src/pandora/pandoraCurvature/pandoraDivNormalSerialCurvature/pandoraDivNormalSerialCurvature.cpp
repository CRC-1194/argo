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

#include "pandoraDivNormalSerialCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "error.H"

#include "processorFvPatchField.H"
#include "reconstructionSchemes.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraDivNormalSerialCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraDivNormalSerialCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraDivNormalSerialCurvature::pandoraDivNormalSerialCurvature
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

volScalarField& pandoraDivNormalSerialCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (meshDb.found(fieldName_))
    {
        //const volVectorField& interfaceNormals = 
        //    mesh().lookupObject<volVectorField>(fieldName_);
        reconstructionSchemes* surf = 
            mesh().getObjectPtr<reconstructionSchemes>("reconstructionScheme");

        const volVectorField& interfaceNormals = surf->normal();

        averagedNormals_ = interfaceNormals /
        (
            mag(interfaceNormals) + 
            dimensionedScalar(
                "SMALL", interfaceNormals.dimensions(), SMALL
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
            if (surf->interfaceCell()[cellI])
            //if (mag(averagedNormals_[cellI]) != 0)
            {
                markers[cellI] = 0;
            }
        }
        markers.correctBoundaryConditions();

        const auto& owner = mesh().owner();
        const auto& neibr = mesh().neighbour();

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
            markers.correctBoundaryConditions();
        }

        for (label i = 0; i < nPropagate_; ++i)
        {
            volVectorField avgNormTmp = averagedNormals_;
            forAll (markers, cellI)
            {
                if (markers[cellI] != i + 1) continue;

                avgNormTmp[cellI] = Zero;
                const point& p = mesh().C()[cellI];

                labelHashSet hashSet;

                const labelList& cps = mesh().cellPoints()[cellI];
                forAll (cps, pI)
                {
                    const labelList& pcs = mesh().pointCells()[cps[pI]];
                    forAll (pcs, cI)
                    {
                        label cellJ = pcs[cI];

                        if (markers[cellJ] == -1) continue;
                        if (markers[cellJ] == i + 1) continue;

                        if (hashSet.found(cellJ)) continue;
                        hashSet.insert(cellJ);

                        vector n = averagedNormals_[cellJ];
                        if (mag(n) != 0)
                        {
                            n /= mag(n);

                            point centre{Zero};
                            centre = mesh().C()[cellJ];

                            vector dist = centre - p;

                            vector distToSurf = dist & n / mag(n) * n;
                            vector verticalDist = dist - distToSurf;
                            avgNormTmp[cellI] += n / max(mag(verticalDist), SMALL);
                        }
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

        #include "error.hpp"

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
            << "pandoraDivNormalSerialCurvature::cellCurvature \n"
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
