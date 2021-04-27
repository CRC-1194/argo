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

#include "triSurfaceAdapter.hpp"

namespace Foam::TriSurfaceImmersion {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar triSurfaceAdapter::compute_ref_length()
{
    // Use minimum edge length as reference length
    const auto& edges = surface_.edges();
    const auto& vertices = surface_.localPoints();

    scalar min_length = 1.0e15;

    for (const auto& anEdge : edges)
    {
        min_length = min(min_length, mag(vertices[anEdge.start()] - vertices[anEdge.end()]));
    }

    return min_length;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
triSurfaceAdapter::triSurfaceAdapter
(
    const triSurface& surface,
    const triSurfaceSearch& search,
    vector span
)
:
    surface_{surface},
    search_{search},
    span_{span}
{
    ref_length_ = compute_ref_length();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar triSurfaceAdapter::signedDistance(const point& trialPoint) const
{
    const auto& normals = surface_.faceNormals();
    auto hitInfo = search_.nearest(trialPoint, span_);

    if (!hitInfo.hit())
    {
        FatalErrorInFunction
            << "No closest point on surface found: "
            << "trial point = " << trialPoint
            << ", span vector = " << span_
            << exit(FatalError);
    }

    vector delta_v{trialPoint - hitInfo.hitPoint()};

    return mag(delta_v)*sign(delta_v&normals[hitInfo.index()]);
}

scalar triSurfaceAdapter::referenceLength() const
{
    return ref_length_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
