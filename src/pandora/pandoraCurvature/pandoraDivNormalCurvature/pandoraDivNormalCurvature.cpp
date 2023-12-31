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
#include "error.H"

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
        nAverages_(curvatureDict_.getOrDefault<label>("nAverages", 3)), 
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

        // Use averaging to propagate the interface normal vectors  
        // into the bulk.
        averagedNormals_ == interfaceNormals /
        (
            mag(interfaceNormals) + 
            dimensionedScalar(
                "SMALL", interfaceNormals.dimensions(), SMALL
            )
        );
        for (label i = 0; i < nAverages_; ++i)
        {
            averagedNormals_ = fvc::average(averagedNormals_);

            // Re-set the smoothed normal vectors in interface cells. 
            forAll(interfaceNormals, cellI)
            {
                if ((interfaceNormals[cellI][0] != 0) ||
                    (interfaceNormals[cellI][1] != 0) ||
                    (interfaceNormals[cellI][2] != 0)) 
                {
                    averagedNormals_[cellI] = interfaceNormals[cellI] / 
                        mag(interfaceNormals[cellI]);
                }
            }
        }

        // Normalize normal vectors.
        averagedNormals_ /= mag(averagedNormals_) + 
            dimensionedScalar("SMALL", averagedNormals_.dimensions(), SMALL); 

        cellCurvature_ = -fvc::div(averagedNormals_);
    }
    else
    {
        FatalErrorInFunction
            << "pandoraDivNormalCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
