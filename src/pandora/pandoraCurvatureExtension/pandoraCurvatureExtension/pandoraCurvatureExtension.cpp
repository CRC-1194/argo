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

#include "pandoraCurvatureExtension.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraCurvatureExtension, false);
defineRunTimeSelectionTable(pandoraCurvatureExtension, Dictionary);

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

autoPtr<pandoraCurvatureExtension> pandoraCurvatureExtension::New 
(
    const dictionary& dict 
)
{
    const dictionary extensionDict = dict.subDict("curvatureExtension");
    const word name = extensionDict.get<word>("type"); 

    // Find the constructor pointer for the model in the constructor table.
    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "Foam::pandoraCurvatureExtension::New(const word&, const Dictionary&)" 
        )   << "Unknown curvature model type "
            << name << nl << nl
            << "Valid curvature models are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<pandoraCurvatureExtension> (cstrIter()(extensionDict));
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraCurvatureExtension::pandoraCurvatureExtension(const dictionary& dict)
:
    extensionDict_{dict}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
