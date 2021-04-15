/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

#include "volumeFractionCalculator.hpp"

#include "DynamicList.H"
#include "IOobject.H"
#include "addToRunTimeSelectionTable.H" 
#include "dictionary.H"
#include "fvc.H"
#include "pointIndexHit.H"
#include "pointMesh.H"

#include "insideOutsidePropagation.hpp"
#include "signedDistanceCalculator.hpp"
#include "volFieldsFwd.H"

namespace Foam::TriSurfaceImmersion {

    defineTypeNameAndDebug(volumeFractionCalculator, 0);
    defineRunTimeSelectionTable(volumeFractionCalculator, Dictionary)


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
volumeFractionCalculator::volumeFractionCalculator
(
    const dictionary& configDict,
    const fvMesh& mesh
)
:
    mesh_{mesh},
    runTime_{mesh.time()},
    pMesh_{mesh},
    writeGeometry_(configDict.get<Switch>("writeGeometry"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
autoPtr<volumeFractionCalculator>
volumeFractionCalculator::New
(
    const dictionary& configDict,
    const fvMesh& mesh
)
{
    const word name = configDict.get<word>("algorithm");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "volumeFractionCalculator::New(const word& name)"
        )   << "Unknown volumeFractionCalculator type "
            << name << nl << nl
            << "Valid volumeFractionCalculators are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<volumeFractionCalculator>(cstrIter()(configDict, mesh));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void volumeFractionCalculator::bulkVolumeFraction(volScalarField& alpha) const
{
    alpha = pos(this->sigDistCalc().cellSignedDist());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //