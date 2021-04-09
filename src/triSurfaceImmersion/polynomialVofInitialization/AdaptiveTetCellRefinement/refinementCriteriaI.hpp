template<class Surface>
bool Foam::TriSurfaceImmersion::boundingBallCriterion<Surface>::needsRefinement
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


template<class Surface>
scalar Foam::TriSurfaceImmersion::boundingBallCriterion<Surface>::levelSetValue(const Surface& surface, const point& x)
{
    return surface.signedDistance(x);
}

template<class Surface>
std::tuple<scalar, label> Foam::TriSurfaceImmersion::boundingBallCriterion<Surface>::maxDistSqrAndPointID
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


template<class Surface>
bool Foam::TriSurfaceImmersion::signCriterion<Surface>::needsRefinemeent
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


template<class Surface>
scalar Foam::TriSurfaceImmersion::signCriterion<Surface>::levelSetValue(const Surface& surface, const point& x)
{
    return surface.value(x);
}