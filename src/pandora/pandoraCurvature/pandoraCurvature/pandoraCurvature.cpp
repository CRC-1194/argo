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

#include "pandoraCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraCurvature, false);
defineRunTimeSelectionTable(pandoraCurvature, Dictionary)

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

autoPtr<pandoraCurvature> pandoraCurvature::New 
(
    const fvMesh& mesh, 
    const dictionary& dict 
)
{
    const word name = dict.subDict("curvature").get<word>("type"); 

    auto* ctorPtr = DictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "pandoraCurvature",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<pandoraCurvature>(ctorPtr(mesh, dict));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraCurvature::pandoraCurvature
(
    const fvMesh& mesh, 
    const dictionary& dict
)
    :
        curvatureDict_(dict.subDict("curvature")),
        cellCurvature_
        (
            IOobject
            (
                "cellCurvature", 
                mesh.time().timeName(), 
                mesh, 
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ), 
            mesh, 
            dimensionedScalar("cellCurvature", dimless/dimLength, Zero)
        )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraCurvature::cellCurvature() 
{
    return cellCurvature_;
}

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
