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
    Foam::tetVofCalculator

Description

SourceFiles
    tetVofCalculatorI.H
    tetVofCalculator.C
    tetVofCalculatorIO.C

\*---------------------------------------------------------------------------*/

#ifndef tetVofCalculator_H
#define tetVofCalculator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AdaptiveTetCellRefinement.hpp"

#include <array>

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                         Class tetVofCalculator Declaration
\*---------------------------------------------------------------------------*/

class tetVofCalculator
{
    // Private data
    mutable std::array<scalar, 4> signed_distance_buffer_;

    // Private Member Functions
    label count_negative_entries() const;


public:

    // Constructors
    
    // Member Functions
    scalar volume(const indexedTet& t, const std::vector<point>& p) const;

    scalar vof
           (
                const indexedTet& tet,
                const std::vector<scalar>& signed_distance
           ) const;
    scalar omega_plus_volume
           (
                const indexedTet& tet,
                const std::vector<scalar>& signed_distance,
                const std::vector<point>& points
           ) const;
    std::vector<scalar> vof
                        (
                            const std::vector<indexedTet>& tets,
                            const std::vector<scalar>& signed_distance
                        ) const;
    std::vector<scalar> omega_plus_volume
                        (
                            const std::vector<indexedTet>& tets,
                            const std::vector<scalar>& signed_distance,
                            const std::vector<point>& points
                        ) const;
    scalar accumulated_omega_plus_volume
           (
                const std::vector<indexedTet>& tets,
                const std::vector<scalar>& signed_distance,
                const std::vector<point>& points
           ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
