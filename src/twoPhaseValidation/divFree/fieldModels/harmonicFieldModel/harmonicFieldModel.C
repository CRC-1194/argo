/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2011 Tomislav Maric 
     \\/     M anipulation  |
-------------------------------------------------------------------------------

License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Multiplies the velocity with cos(pi t / T) as introduced by Rider and Kothe.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#include "harmonicFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(harmonicFieldModel, 0);
addToRunTimeSelectionTable(divFreeFieldModel, harmonicFieldModel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vector harmonicFieldModel::velocity(point X, scalar t) const
{
    using namespace Foam::constant::mathematical;

    return getDecoratedMotion().velocity(X,t) * Foam::cos((pi * t)  / getPeriod());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

harmonicFieldModel::harmonicFieldModel
(
    const Time& time,
    const dictionary& dict
)
:
    periodicFieldModel(time, dict)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
