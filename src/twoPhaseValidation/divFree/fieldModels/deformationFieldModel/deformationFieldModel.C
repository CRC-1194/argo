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
    3D deformation divergence free field model defined in [1, example 11.2] 

    [1] R. LeVeque, High-resolution conservative algorithms for advection in
    incompressible flow, SIAM J. Numer. Anal. 33, 627 (1996).

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "deformationFieldModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(deformationFieldModel, 0);
addToRunTimeSelectionTable(divFreeFieldModel, deformationFieldModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

deformationFieldModel::deformationFieldModel
(
    const Time& time,
    const dictionary& dict
)
:
    divFreeFieldModel(time, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector deformationFieldModel::velocity(point X, scalar) const
{
    using namespace Foam::constant::mathematical;

    scalar Ux =   2 *sqr(sin(pi * X.x())) * sin(2 * pi * X.y()) * sin(2 * pi * X.z());

    scalar Uy = -1 * sin(2 * pi * X.x()) * sqr(sin(pi * X.y())) * sin(2 * pi * X.z());

    scalar Uz = -1 * sin(2 * pi * X.x()) * sin(2 * pi * X.y()) * sqr(sin(pi * X.z()));

    return vector(Ux, Uy, Uz);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
