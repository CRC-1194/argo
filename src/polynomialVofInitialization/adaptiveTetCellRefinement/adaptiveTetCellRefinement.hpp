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

#include "orientedSurface.H"

#include "fvCFD.H"

#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace PolynomialVof {

using indexedTet = std::array<label, 4>;
using indexTuple = std::pair<label, label>;
using edge = indexTuple;

/*---------------------------------------------------------------------------*\
                         Class adaptiveTetCellRefinement Declaration
\*---------------------------------------------------------------------------*/

class adaptiveTetCellRefinement
{
    // Private data
    orientedSurface surface_;

    std::vector<point> points_;
    std::vector<scalar> signed_distance_;

    std::vector<indexedTet> tets_;
    std::vector<bool> refinement_required_;

    std::unordered_map<edge, label> edge_to_point_id_;

    // NOTE: the second entry in the tuple denotes the index of the first element
    // which is NOT part of the level, aka "[begin, end)" (TT)
    std::vector<indexTuple> level_to_pointid_range_;
    std::vector<indexTuple> level_to_tetid_range_;

    label max_refinement_level_;


    // Private Member Functions
    void compute_decomposition();
    label compute_max_refinement_level();

public:

    // Constructors
    adaptiveTetCellRefinement
    (
        const orientedSurface& surface,
        const std::vector<point> points,
        const std::vector<scalar> signed_distance,
        const std::vector<indexedTet> tets
    );
    

    // Member Functions
    const std::vector<point>& points();
    const std::vector<scalar>& signed_distance();
    std::vector<indexedTet> resulting_tets();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
