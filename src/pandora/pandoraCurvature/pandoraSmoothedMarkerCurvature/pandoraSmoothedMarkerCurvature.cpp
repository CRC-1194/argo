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

#include "pandoraSmoothedMarkerCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "error.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraSmoothedMarkerCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraSmoothedMarkerCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraSmoothedMarkerCurvature::pandoraSmoothedMarkerCurvature
(
    const fvMesh& mesh, 
    const dictionary& dict
)
    :
        pandoraCurvature(mesh, dict), 
        markerFieldName_(curvatureDict_.get<word>("markerField")),
        markerTolerance_(
            curvatureDict_.getOrDefault<scalar>("markerTolerance", 1e-03)
        ),
        nMarkerAverages_(
            curvatureDict_.getOrDefault<label>("nMarkerAverages", 3)
        ), 
        nNormalAverages_(
            curvatureDict_.getOrDefault<label>("nNormalAverages", 3)
        )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraSmoothedMarkerCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (meshDb.found(markerFieldName_))
    {
        const volScalarField& marker = 
            mesh().lookupObject<volScalarField>(markerFieldName_);

        volScalarField averageMarker(
            marker.name() + "-average", 
            marker 
        );

        for (label i = 0; i < nMarkerAverages_; ++i)
        {
            averageMarker = fvc::average(averageMarker);
        }

        volVectorField normals (
            "normals", 
            fvc::grad(averageMarker)
        );

        normals /= Foam::mag(normals) + 
            dimensionedScalar(
                "SMALL", normals.dimensions(), SMALL
            );

        volVectorField averagedNormals ("averagedNormals", normals);

        for (label i = 0; i < nNormalAverages_; ++i)
        {
            averagedNormals = fvc::average(averagedNormals);

            // Normalize averaged normal vectors.
            averagedNormals /= mag(averagedNormals) + 
                dimensionedScalar(
                    "SMALL", averagedNormals.dimensions(), SMALL
                ); 

            // Replace smoothed normal vectors in interface cells 
            // with large interface areas by normal vectors given by 
            // div (grad(averagedMarker) / magGrad(averagedMarker)). TM.
            forAll(marker, cellI)
            {
                // FIXME: Set 0.1 to attribute. 
                if ((marker[cellI] > 1e-01) && 
                    (marker[cellI] < (1. - 1e-01)))
                {
                    averagedNormals[cellI] = normals[cellI]; 
                }
            }
        }

        cellCurvature_ = -fvc::div(averagedNormals);
    }
    else
    {
        FatalErrorInFunction
            << "pandoraSmoothedMarkerCurvature::cellCurvature \n"
            << "Field " << markerFieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
