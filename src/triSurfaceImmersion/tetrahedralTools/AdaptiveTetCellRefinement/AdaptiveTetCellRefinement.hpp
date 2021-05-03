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

Class
    Foam::adaptiveTetCellRefinement

Description

SourceFiles
    adaptiveTetCellRefinementI.H
    adaptiveTetCellRefinement.C
    adaptiveTetCellRefinementIO.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveTetCellRefinement_H
#define adaptiveTetCellRefinement_H

#include "fvCFD.H"
#include "tmp.H"

#include <algorithm>
#include <array>
#include <filesystem>
#include <map>
#include <utility>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {

/*---------------------------------------------------------------------------*\
                         Class adaptiveTetCellRefinement Declaration
\*---------------------------------------------------------------------------*/

using indexedTet = std::array<label, 4>;

template<class LevelSet, template<class LevelSet2> class RefinementCriterion>
class adaptiveTetCellRefinement
{
private:

    // Typedefs
    using indexTuple = std::tuple<label, label>;
    using edge = indexTuple;

    // Private data
    const LevelSet& levelSet_;
    RefinementCriterion<LevelSet> criterion_;

    std::vector<point> points_;
    std::vector<scalar> levelSetValues_;
    std::vector<indexedTet> tets_;
    std::vector<bool> refinement_required_;
    std::map<edge, label> edge_to_point_id_;

    // NOTE: the second entry in the tuple denotes the index of the first element
    // which is NOT part of the level, aka "[begin, end)" (TT)
    std::vector<indexTuple> level_to_pointid_range_;
    std::vector<indexTuple> level_to_tetid_range_;

    label max_refinement_level_;

    bool decomposition_performed_ = false;
    const label n_tets_from_decomposition = 8;


    // Private Member Functions
    std::array<scalar, 6> edge_lengths(const indexedTet& tet) const;
    label compute_max_refinement_level();
    void compute_decomposition();
    label flag_tets_for_refinement(int level);
    void update_tet_container_sizes(int level, int n_new_tets);
    void update_edge_to_point_map(int level);
    void add_to_map(std::array<edge, 6> tet_edges);
    std::array<edge, 6> edges(const indexedTet& tet) const;
    void create_refined_tets(int level);
    void decompose_and_add_new_tets(const indexedTet& tet, label start_id);
    void computeLevelSetValues(int level);
    void save_decomposition_as_vtk
    (
        const std::vector<indexedTet>& tets,
        const std::vector<point>& points,
        const std::vector<scalar>& signed_distance,
        const std::vector<label>& refinement_levels,
        std::string file_name
    ) const;


public:

    // Constructors
    adaptiveTetCellRefinement
    (
        const LevelSet& levelSet,
        std::vector<point> points,
        std::vector<scalar> levelSetValues,
        std::vector<indexedTet> tets,
        label max_refine_level = -1
    );
    

    // Member Functions
    const std::vector<point>& points();
    const std::vector<scalar>& signedDistance();
    std::vector<indexedTet> resultingTets();

    label refinementLevel() const;
    std::vector<label> refinementLevels(label n_tets);
    void writeTets(label cellID);

    void printLevelInfos() const;
    void printTets() const;
    void printPoints() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AdaptiveTetCellRefinementI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
