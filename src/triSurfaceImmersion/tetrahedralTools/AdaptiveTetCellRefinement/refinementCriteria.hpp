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

namespace Foam::TriSurfaceImmersion {

template<class LevelSet>
class boundingBallCriterion
{
public:
    // The following member function defines the informal interface
    // expected by the AdaptiveRefinement class
    static bool needsRefinement
    (
        const indexedTet& tet,
        const std::vector<point>& points,
        const std::vector<scalar>& signedDistances
    );
    // The criterion determines whether an arbitrary level set field or
    // a signed distance field is required. Thus, the corresponding function
    // is called through the criterion class
    static scalar levelSetValue(const LevelSet& ls, const point& p);
    
    // Criterion specific functions
    static std::tuple<scalar, label> maxDistSqrAndPointID
    (
        const indexedTet& tet,
        const std::vector<scalar>& signedDistances
    );
};

template<class LevelSet>
class signCriterion
{
public:
    static bool needsRefinemeent
    (
        const indexedTet& tet,
        const std::vector<point>& points,
        const std::vector<scalar>& levelSetValues
    );
    static scalar levelSetValue(const LevelSet& ls, const point& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refinementCriteriaI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //