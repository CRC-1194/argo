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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class LevelSet>
bool Foam::TriSurfaceImmersion::boundingBallCriterion<LevelSet>::needsRefinement
(
    const indexedTet& tet,
    const std::vector<point>& points,
    const std::vector<scalar>& signedDistances
)
{
    const auto [max_dist_sqr, max_p_id] = maxDistSqrAndPointID(tet, signedDistances);

    // Bounding ball criterion (TT)
    for (const auto p_id : tet)
    {
        if (Foam::magSqr(points[p_id] - points[max_p_id]) >= max_dist_sqr)
        {
            return true;
        }
    }
    
    return false;
}


template<class LevelSet>
scalar Foam::TriSurfaceImmersion::boundingBallCriterion<LevelSet>::levelSetValue(const LevelSet& ls, const point& p)
{
    return ls.signedDistance(p);
}

template<class LevelSet>
std::tuple<scalar, label> Foam::TriSurfaceImmersion::boundingBallCriterion<LevelSet>::maxDistSqrAndPointID
(
    const indexedTet& tet,
    const std::vector<scalar>& signedDistances
)
{
    scalar max_dist_sqr{0.0};
    label max_p_id{tet[0]};

    for (const auto p_id : tet)
    {
        if ((signedDistances[p_id]*signedDistances[p_id]) > max_dist_sqr)
        {
            max_dist_sqr = signedDistances[p_id]*signedDistances[p_id];
            max_p_id = p_id;
        }
    }

    return std::make_tuple(max_dist_sqr, max_p_id);
}


template<class LevelSet>
bool Foam::TriSurfaceImmersion::signCriterion<LevelSet>::needsRefinemeent
(
    const indexedTet &tet,
    const std::vector<point>& points,
    const std::vector<scalar>& levelSetValues
)
{
    scalar lsSign = sign(levelSetValues[tet[0]]);

    for (const auto pI : tet)
    {
        if (sign(levelSetValues[pI]) != lsSign)
        {
            return true;
        }
    }

    return false;
}


template<class LevelSet>
scalar Foam::TriSurfaceImmersion::signCriterion<LevelSet>::levelSetValue(const LevelSet& ls, const point& p)
{
    return ls.value(p);
}

// ************************************************************************* //