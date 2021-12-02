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
    A simple rotation field model.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "rotationFieldModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rotationFieldModel, 0);
addToRunTimeSelectionTable(divFreeFieldModel, rotationFieldModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotationFieldModel::rotationFieldModel
(
    const Time& time,
    const dictionary& dict
)
:
    translationFieldModel(time, dict),
    axisPoint_(dict.lookup("axisPoint"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector rotationFieldModel::velocity(point X, scalar) const
{
    // Compute the rotation velocity at a point X using the angular velocity V
    // and the axis point Ap.

    const vector ApX = X - axisPoint_;

    const vector& v = getVector();

    const vector vHat = v / (mag(v) + SMALL);

    const vector r = ApX  - (ApX & vHat) * vHat;

    return v ^ r;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
