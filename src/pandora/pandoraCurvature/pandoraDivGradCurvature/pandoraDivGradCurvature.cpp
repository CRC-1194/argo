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

#include "pandoraDivGradCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "error.H"
#include "fvcGrad.H"
#include "fvcDiv.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraDivGradCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraDivGradCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraDivGradCurvature::pandoraDivGradCurvature
(
    const fvMesh& mesh, 
    const dictionary& dict
)
    :
        pandoraCurvature(mesh, dict), 
        fieldName_(this->curvatureDict().get<word>("fieldName"))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraDivGradCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (meshDb.found(fieldName_))
    {
        const auto& vf = 
            mesh().lookupObject<volScalarField>(fieldName_);

        volVectorField vfGrad{fvc::grad(vf)}; 
        volScalarField vfMagGrad{mag(vfGrad)}; 
        vfMagGrad += dimensionedScalar(
            "vfMagGrad", 
            vfMagGrad.dimensions(),
            SMALL
        );
        cellCurvature_ = -fvc::div(vfGrad / vfMagGrad);
    }
    else
    {
        FatalErrorInFunction
            << "pandoraDivGradCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
