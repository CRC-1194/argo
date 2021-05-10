/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::TriSurfaceImmersion::boundingBallCriterion
    Foam::TriSurfaceImmersion::signCriterion

Description
    Criteria to decide whether a cell is intersected by an interface based
    on the level set values or the signed distances defined at the cell
    corner points.

SourceFiles
    refinementCriteriaI.hpp

\*---------------------------------------------------------------------------*/

#ifndef refinementCriteria_H
#define refinementCriteria_H

#include "fvCFD.H"

#include "AdaptiveTetCellRefinement.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

template<class LevelSet>
class boundingBallCriterion
{
public:

    // Generic interface member functions
    //- Determine wheter the given tetrahedron needs refinement.
    //  This function evaluates to true if the there is no bounding ball
    //  centred around a tetrahedron vertex that is smaller than the ball
    //  defined by the vertex and its signed distance as radius.
    static bool needsRefinement(const indexedTet& tet,
        const std::vector<point>& points,
        const std::vector<scalar>& signedDistances);

    //- Return signed distance of p to the interface.
    static scalar levelSetValue(const LevelSet& ls, const point& p);

    // Criterion specific member functions
    //- Return a tuple with the maximum distance squared and vertex label.
    static std::tuple<scalar, label> maxDistSqrAndPointID(
        const indexedTet& tet, const std::vector<scalar>& signedDistances);
};


template<class LevelSet>
class signCriterion
{
public:

    // Generic interface member functions
    //- Determine wheter the given tetrahedron needs refinement.
    //  This function evaluates to true if the level set values do not have
    //  all the same sign.
    static bool needsRefinemeent(const indexedTet& tet,
        const std::vector<point>& points,
        const std::vector<scalar>& levelSetValues);

    //- Return level set value at point p. Not necessarily a signed distance.
    static scalar levelSetValue(const LevelSet& ls, const point& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RefinementCriteriaI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
