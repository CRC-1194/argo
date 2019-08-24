/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
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


namespace Foam {
namespace PolynomialVof {

#include "orientedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void orientedSurface::updateDistanceToOrigin()
{
    distanceOrigin_ = unitNormal_ & refPoint_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
orientedSurface::orientedSurface(const point& refPoint, const vector& normal, scalar refLength)
:
    refPoint_{refPoint},
    unitNormal_{unitNormal},
    refLength_{refLength}
{
    updateDistanceToOrigin();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar orientedSurface::signedDistance(const point& trialPoint) const
{
    return (trialPoint & unitNormal_) - distanceOrigin_;
}

scalar orientedSurface::referenceLength() const
{
    return refLength_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
