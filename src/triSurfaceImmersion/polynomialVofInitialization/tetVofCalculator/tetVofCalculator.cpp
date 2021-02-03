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

#include "tetVofCalculator.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace Foam {
namespace PolynomialVof {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
label tetVofCalculator::count_negative_entries() const
{
    label n_negative_entries = 0;

    for (const auto entry : signed_distance_buffer_)
    {
        if (entry < 0.0) ++n_negative_entries;
    }

    return n_negative_entries;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar tetVofCalculator::volume(const indexedTet& t, const std::vector<point>& p) const
{
    return mag((p[t[1]] - p[t[0]]) & ((p[t[2]] - p[t[0]])^(p[t[3]] - p[t[0]])))/6.0;
}

scalar tetVofCalculator::vof
       (
            const indexedTet& tet,
            const std::vector<scalar>& signed_distance
       ) const
{
    // This function implements the actual model of Detrixhe and Aslam (TT)
    scalar volFraction = 0.0;

    for (int idx = 0; idx !=4; ++idx)
    {
        signed_distance_buffer_[idx] = signed_distance[tet[idx]];
    }

    std::sort(signed_distance_buffer_.begin(), signed_distance_buffer_.end());
    const auto& d = signed_distance_buffer_;

    auto negative_distances = count_negative_entries();

    if (negative_distances == 4)
    {
        volFraction = 0.0;
    }
    else if (negative_distances == 3)
    {
        volFraction = std::pow(d[3],3) /
               ((d[3] - d[0]) * (d[3] - d[1]) * (d[3] - d[2]));
    }
    else if (negative_distances == 2)
    {
        volFraction = ( d[0]*d[1] * (d[2]*d[2] + d[2]*d[3] + d[3]*d[3])
                + d[2]*d[3] * (d[2]*d[3] - (d[0]+d[1])*(d[2]+d[3])) )
                / ((d[0]-d[2]) * (d[1]-d[2]) * (d[0]-d[3]) * (d[1]-d[3]));
    }
    else if (negative_distances == 1)
    {
        volFraction = 1.0 + std::pow(d[0],3) / 
                ( (d[1]-d[0]) * (d[2]-d[0]) * (d[3]-d[0]) );
    }
    else
    {
        volFraction = 1.0;
    }

    assert (volFraction >= 0.0 && volFraction <= 1.0);
    return volFraction;
}

scalar tetVofCalculator::omega_plus_volume
       (
            const indexedTet& tet,
            const std::vector<scalar>& signed_distance,
            const std::vector<point>& points
       ) const
{
    auto volume_fraction = vof(tet, signed_distance);

    return volume_fraction*volume(tet, points);
}

std::vector<scalar> tetVofCalculator::vof
                    (
                        const std::vector<indexedTet>& tets,
                        const std::vector<scalar>& signed_distance
                    ) const
{
    std::vector<scalar> volume_fractions(tets.size());

    for (uint idx = 0; idx != tets.size(); ++idx)
    {
        volume_fractions[idx] = vof(tets[idx], signed_distance);
    }

    return volume_fractions;
}

std::vector<scalar> tetVofCalculator::omega_plus_volume
                    (
                        const std::vector<indexedTet>& tets,
                        const std::vector<scalar>& signed_distance,
                        const std::vector<point>& points
                    ) const
{
    std::vector<scalar> plus_volumes(tets.size());

    for (uint idx = 0; idx != tets.size(); ++idx)
    {
        plus_volumes[idx] = omega_plus_volume(tets[idx], signed_distance, points);
    }

    return plus_volumes;
}

scalar tetVofCalculator::accumulated_omega_plus_volume
       (
            const std::vector<indexedTet>& tets,
            const std::vector<scalar>& signed_distance,
            const std::vector<point>& points
       ) const
{
    auto v = omega_plus_volume(tets, signed_distance, points);
    return std::accumulate(v.begin(), v.end(), 0.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
