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


#include "orientedPlane.hpp"

namespace Foam::TriSurfaceImmersion {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void orientedPlane::updateDistanceToOrigin()
{
    distance_origin_ = unit_normal_ & ref_point_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
orientedPlane::orientedPlane(const point& refPoint, const vector& normal, scalar refLength)
:
    ref_point_{refPoint},
    unit_normal_{normal},
    ref_length_{refLength},
    distance_origin_{0.0}
{
    updateDistanceToOrigin();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar orientedPlane::signedDistance(const point& trialPoint) const
{
    return (trialPoint & unit_normal_) - distance_origin_;
}

scalar orientedPlane::referenceLength() const
{
    return ref_length_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
