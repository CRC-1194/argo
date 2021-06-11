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

#include "pandoraCurvatureNormalRegularisation.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraCurvatureNormalRegularisation, false);
addToRunTimeSelectionTable(pandoraCurvatureRegularisation, pandoraCurvatureNormalRegularisation, Dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraCurvatureNormalRegularisation::pandoraCurvatureNormalRegularisation(const dictionary& dict)
    :
        pandoraCurvatureRegularisation{dict}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void pandoraCurvatureNormalRegularisation::regularise
(
    volScalarField& curvature,
    const boolList&// Unused argument
)
{
    // Average curvature with fvc::average ("Laplace smoothing").
    const label nAverages = this->regularisationDict().get<label>("nAverages");
    for (int i = 0; i < nAverages; ++i)
        curvature = fvc::average(curvature);
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
