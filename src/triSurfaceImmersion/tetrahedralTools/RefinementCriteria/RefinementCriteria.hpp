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
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

/*
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
*/

//--- Refactored intersection criteria functions
struct boundingBallCriterion{};
struct signCriterion{};

// --- boundingBallCriterion functions
template<class IndexedPolyhedron,
    template<class PointType> class PointContainer,
    template<class ValueType> class LevelSetValueContainer>
bool considerIntersected(const point& refPoint,
    scalar maxDistSqr,
    const IndexedPolyhedron& poly,
    const PointContainer<point>& points,
    const LevelSetValueContainer<scalar>& values,
    const boundingBallCriterion& criterion)
{
    return std::any_of(poly.begin(), poly.end(),
        [&](label pID){return Foam::magSqr(points[pID] - refPoint) >= maxDistSqr;}
    );
}


template<class IndexSet, template<class ValueType> class ValueContainer>
std::tuple<scalar,label> maxDistSqrAndPointID(const IndexSet& ids,
    const ValueContainer<scalar>& values)
{
    label maxID = *(std::max_element(ids.begin(), ids.end(),
        [&](label a, label b){
            return (Foam::magSqr(values[a]) < Foam::magSqr(values[b]));}
    ));

    return std::make_tuple(Foam::magSqr(values[maxID]), maxID);
}


template<class IndexedPolyhedron,
    template<class PointType> class PointContainer,
    template<class ValueType> class LevelSetValueContainer>
bool considerIntersected(const IndexedPolyhedron& poly,
    const PointContainer<point>& points,
    const LevelSetValueContainer<scalar>& values,
    const boundingBallCriterion& criterion)
{
    const auto [maxDistSqr, pID] = maxDistSqrAndPointID(poly, values);

    return considerIntersected(points[pID],
        maxDistSqr, poly, points, values, criterion);
}


template<class LevelSet>
scalar levelSetValue(const LevelSet& ls,
    const point& p,
    const boundingBallCriterion& criterion)
{
    return ls.signedDistance(p);
}


// --- signCriterion functions
template<template<class IndexType> class IndexedPolyhedron,
    template<class PointType> class PointContainer,
    template<class ValueType> class LevelSetValueContainer>
bool considerIntersected(const point& refPoint,
    scalar refValue,
    const IndexedPolyhedron<label>& poly,
    const PointContainer<point>& points,
    const LevelSetValueContainer<scalar>& values,
    const signCriterion& criterion)
{
    return std::any_of(poly.begin(), poly.end(),
        [&](label pID){return values[pID]*refValue < 0.0;}
    );
}


template<template<class IndexType> class IndexedPolyhedron,
    template<class PointType> class PointContainer,
    template<class ValueType> class LevelSetValueContainer>
bool considerIntersected(const IndexedPolyhedron<label>& poly,
    const PointContainer<point>& points,
    const LevelSetValueContainer<scalar>& values,
    const signCriterion& criterion)
{
    return considerIntersected(points[poly[0]],
        values[poly[0]], poly, points, values, criterion);
}


template<class LevelSet>
scalar levelSetValue(const LevelSet& ls,
    const point& p,
    const signCriterion& criterion)
{
    return ls.value(p);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "RefinementCriteriaI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
