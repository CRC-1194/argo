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

#include <algorithm>

namespace Foam {
namespace PolynomialVof {

#include "adaptiveTetCellRefinement.hpp"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void adaptiveTetCellRefinement::compute_decomposition()
{
}

label adaptiveTetCellRefinement::compute_max_refinement_level()
{
    return 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
adaptiveTetCellRefinement
(
    const orientedSurface& surface,
    const std::vector<point> points,
    const std::vector<scalar> signed_distance,
    const std::vector<indexedTet> tets
)
    :
    surface_{surface},
    points_{points},
    signed_distance_{signed_distance},
    tets_{tets},
    refinement_required_{tets_.size(), false},
    edge_to_point_id_{},
    level_to_pointid_range_{},
    level_to_tetid_range_{},
    max_refinement_level_{0}
{
    max_refinement_level_ = compute_max_refinement_level();

    level_to_pointid_range_.resize(max_refinement_level_);
    level_to_pointid_range_[0] = indexTuple{0, points_.size()};

    level_to_tetid_range_.resize(max_refinement_level_);
    level_to_tetid_range_[0] = indexTuple{0, tets_.size()};
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
const std::vector<point>& adaptiveTetCellRefinement::points()
{
    return points_;
}

const std::vector<scalar>& signed_distance()
{
    return signed_distance_;
}

std::vector<indexedTet> resulting_tets()
{
    std::vector<indexedTet> final_tets{};

    n_tets = std::count(refinement_required_.begin(), refinement_required_.end(), false);

    final_tets.resize(n_tets);

    int count = 0;
    for (int idx = 0; idx != refinement_required_.size(); ++idx)
    {
         if (refinement_required_[idx] == false)
         {
             final_tets[count] = tets_[idx];
         }
    }

    return final_tets;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
