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

#include "adaptiveTetCellRefinement.hpp"

#include <algorithm>

namespace Foam {
namespace PolynomialVof {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
label adaptiveTetCellRefinement::compute_max_refinement_level()
{
    return 1;
}

void adaptiveTetCellRefinement::compute_decomposition()
{
    if (decomposition_performed_)
    {
        return;
    }

    for (int level = 0; level != max_refinement_level_; ++level)
    {
        auto n_refined_tets = flag_tets_for_refinement(level);
        update_tet_container_sizes(level, n_refined_tets*n_tets_from_decomposition);
        update_edge_to_point_map(level);
        create_refined_tets(level);
        compute_signed_distances(level);
    }
}

label adaptiveTetCellRefinement::flag_tets_for_refinement(int level) 
{
    auto [start, end] = level_to_tetid_range_[level];

    auto n_tets_to_refine = 0;

    for (auto idx = start; idx != end; ++idx)
    {
        refinement_required_[idx] = has_to_be_refined(tets_[idx]);

        if (refinement_required_[idx])
        {
            ++n_tets_to_refine;
        }
    }

    return n_tets_to_refine; 
}

bool adaptiveTetCellRefinement::has_to_be_refined(const indexedTet& tet) const
{
    auto [max_dist_sqr, max_p_id] = maxiumum_distance_sqr_and_pointid(tet);

    // Bounding sphere criterion (TT)
    for (const auto p_id : tet)
    {
        if (distance_squared(points_[p_id], points_[max_p_id]) >= max_dist_sqr)
        {
            return true;
        }
    }
    
    return false;
}

std::tuple<scalar, label> adaptiveTetCellRefinement::maxiumum_distance_sqr_and_pointid(const indexedTet& tet) const
{
    scalar max_dist_sqr{0.0};
    label max_p_id{tet[0]};

    for (const auto p_id : tet)
    {
        if ((signed_distance_[p_id]*signed_distance_[p_id]) > max_dist_sqr)
        {
            max_dist_sqr = signed_distance_[p_id]*signed_distance_[p_id];
            max_p_id = p_id;
        }
    }

    return std::make_tuple(max_dist_sqr, max_p_id);
}

scalar adaptiveTetCellRefinement::distance_squared(const point& p_a, const point& p_b) const
{
    return ((p_a - p_b)&(p_a - p_b));
}

void adaptiveTetCellRefinement::update_tet_container_sizes(int level, int n_new_tets)
{
    // Update vector sizes related to the number of tets
    level_to_tetid_range_[level + 1] = indexTuple{tets_.size(), tets_.size() + n_new_tets};
    tets_.resize(tets_.size() + n_new_tets); 
    refinement_required_.resize(refinement_required_.size() + n_new_tets);

    // NOTE: only tets whose refine flag is false are returned by the public
    // member functions. Thus, default initialize the newly added fields to
    // false.
    auto [next_start, next_end] = level_to_tetid_range_[level + 1];
    for (auto idx = next_start; idx != next_end; ++idx)
    {
        refinement_required_[idx] = false;
    }
}

void adaptiveTetCellRefinement::add_to_map(std::array<edge, 6> tet_edges)
{
    for (const auto tet_edge : tet_edges)
    {
        edge_to_point_id_[tet_edge] = 0;
    }
}

void adaptiveTetCellRefinement::update_edge_to_point_map(int level)
{
    edge_to_point_id_.clear();

    auto [start, end] = level_to_tetid_range_[level];

    for (auto idx = start; idx != end; ++idx)
    {
        if (refinement_required_[idx])
        {
            add_to_map(edges(tets_[idx]));
        }
    }

    // Update container sizes
    auto n_new_points = edge_to_point_id_.size();
    level_to_pointid_range_[level + 1] = indexTuple{points_.size(), points_.size() + n_new_points};
    points_.resize(points_.size() + n_new_points);
    signed_distance_.resize(signed_distance_.size() + n_new_points);

    // Add global point ids to the mapping and compute the new points
    label global_point_id = std::get<0>(level_to_pointid_range_[level + 1]);

    for (auto& edge_to_point : edge_to_point_id_)
    {
        edge_to_point.second = global_point_id;
        auto [p1_id, p2_id] = edge_to_point.first;
        points_[global_point_id] = 0.5*(points_[p1_id] + points_[p2_id]);

        ++global_point_id;
    }
}

std::array<adaptiveTetCellRefinement::edge, 6> adaptiveTetCellRefinement::edges(const indexedTet& tet) const
{
    return std::array<edge, 6>{
        edge{tet[0], tet[1]}, edge{tet[0], tet[2]}, edge{tet[0], tet[3]},
        edge{tet[1], tet[2]}, edge{tet[1], tet[3]},
        edge{tet[2], tet[3]}
    };
}

void adaptiveTetCellRefinement::create_refined_tets(int level)
{
    auto [start, end] = level_to_tetid_range_[level];
    auto refined_tet_id = end;

    for (auto idx = start; idx != end; ++idx)
    {
        if (refinement_required_[idx])
        {
            decompose_and_add_new_tets(tets_[idx], refined_tet_id); 
            refined_tet_id += n_tets_from_decomposition;
        }
    }
}

void adaptiveTetCellRefinement::decompose_and_add_new_tets(const indexedTet& tet, label start_id)
{
    auto tet_edges = edges(tet);

    // Translate edges into point ids
    std::array<label, 6> pids{};

    for (int idx = 0; idx != 6; ++idx)
    {
        pids[idx] = edge_to_point_id_[tet_edges[idx]];
    }

    // Define new tets
    // Tets formed from one existing point and three points from refinement
    tets_[start_id]     = indexedTet{tet[0], pids[0], pids[1], pids[2]};
    tets_[start_id + 1] = indexedTet{tet[1], pids[0], pids[3], pids[4]};
    tets_[start_id + 2] = indexedTet{tet[2], pids[1], pids[3], pids[5]};
    tets_[start_id + 3] = indexedTet{tet[3], pids[2], pids[4], pids[5]};

    // Decomposition of the octaeder has to be done carefully to avoid
    // overly warped tets (TT)
    // TODO: Not clear if taking the minimal distance results in more regular
    // tets. Needs testing. (TT)
    // NOTE: first tests look promosing (TT)
    auto d0 = distance_squared(points_[pids[0]], points_[pids[5]]);
    auto d1 = distance_squared(points_[pids[1]], points_[pids[4]]);
    auto d2 = distance_squared(points_[pids[2]], points_[pids[3]]);

    // Tets solely constituted by by points from refinement
    if (d0 < d1 && d0 < d2)
    {
        tets_[start_id + 4] = indexedTet{pids[0], pids[1], pids[2], pids[5]};
        tets_[start_id + 5] = indexedTet{pids[0], pids[1], pids[3], pids[5]};
        tets_[start_id + 6] = indexedTet{pids[0], pids[2], pids[4], pids[5]};
        tets_[start_id + 7] = indexedTet{pids[0], pids[3], pids[4], pids[5]};
    }
    else if (d1 < d0 && d1 < d2)
    {
        tets_[start_id + 4] = indexedTet{pids[1], pids[0], pids[2], pids[4]};
        tets_[start_id + 5] = indexedTet{pids[1], pids[2], pids[5], pids[4]};
        tets_[start_id + 6] = indexedTet{pids[1], pids[5], pids[3], pids[4]};
        tets_[start_id + 7] = indexedTet{pids[1], pids[3], pids[0], pids[4]};
    }
    else
    {
        tets_[start_id + 4] = indexedTet{pids[2], pids[0], pids[1], pids[3]};
        tets_[start_id + 5] = indexedTet{pids[2], pids[1], pids[5], pids[3]};
        tets_[start_id + 6] = indexedTet{pids[2], pids[5], pids[4], pids[3]};
        tets_[start_id + 7] = indexedTet{pids[2], pids[4], pids[0], pids[3]};
    }
}

void adaptiveTetCellRefinement::compute_signed_distances(int level)
{
    auto [start, end] = level_to_pointid_range_[level + 1];

    for (auto idx = start; idx != end; ++idx)
    {
        signed_distance_[idx] = surface_.signedDistance(points_[idx]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
adaptiveTetCellRefinement::adaptiveTetCellRefinement
(
    const orientedPlane& surface,
    const std::vector<point> points,
    const std::vector<scalar> signed_distance,
    const std::vector<indexedTet> tets,
    const label max_refine_level = -1
)
    :
    surface_{surface},
    points_{points},
    signed_distance_{signed_distance},
    tets_{tets},
    refinement_required_(tets_.size(), false),
    edge_to_point_id_{},
    level_to_pointid_range_{},
    level_to_tetid_range_{},
    max_refinement_level_{max_refine_level}
{
    if (max_refinement_level_ < 0)
    {
        max_refinement_level_ = compute_max_refinement_level();
    }

    level_to_pointid_range_.resize(max_refinement_level_ + 1);
    level_to_pointid_range_[0] = indexTuple{0, points_.size()};

    level_to_tetid_range_.resize(max_refinement_level_ + 1);
    level_to_tetid_range_[0] = indexTuple{0, tets_.size()};
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
const std::vector<point>& adaptiveTetCellRefinement::points()
{
    compute_decomposition();

    return points_;
}

const std::vector<scalar>& adaptiveTetCellRefinement::signed_distance()
{
    compute_decomposition();

    return signed_distance_;
}

std::vector<indexedTet> adaptiveTetCellRefinement::resulting_tets()
{
    compute_decomposition();

    std::vector<indexedTet> final_tets{};

    auto n_tets = std::count(refinement_required_.begin(), refinement_required_.end(), false);

    final_tets.resize(n_tets);

    int count = 0;
    for (int idx(refinement_required_.size()-1); idx != -1; --idx)
    {
         if (refinement_required_[idx] == false)
         {
             final_tets[count] = tets_[idx];
             ++count;
         }
    }

    return final_tets;
}

void adaptiveTetCellRefinement::print_level_infos() const
{
    Info << "Number of tets and points per level\n";
    for (auto level = 0; level <= max_refinement_level_; ++level)
    {
        auto [tstart, tend] = level_to_tetid_range_[level];
        auto [pstart, pend] = level_to_pointid_range_[level];
        Info << "Level = " << level << ":\n"
             << "\tn_tets = " << (tend - tstart) << "\n"
             << "\tn_points = " << (pend - pstart) << "\n";

    }
}

void adaptiveTetCellRefinement::print_tets() const
{
    Info << "Indexed tets (n = " << tets_.size() << " in total):\n";

    for (uint idx = 0; idx != tets_.size(); ++idx)
    {
        const auto& tet = tets_[idx];
        Info << "tet_id = " << idx << ": "
             << tet[0] << "\t" << tet[1] << "\t"
             << tet[2] << "\t" << tet[3] << "\n";
    }
}

void adaptiveTetCellRefinement::print_points() const
{
    Info << "Points:\n";

    for (uint idx = 0; idx != points_.size(); ++idx)
    {
        const auto& p = points_[idx];

        Info << "p_id = " << idx << ": " << p << "\n";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
