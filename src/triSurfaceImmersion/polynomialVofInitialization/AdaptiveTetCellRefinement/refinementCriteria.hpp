#ifndef refinementCriteria_H
#define refinementCriteria_H

#include "fvCFD.H"

#include "AdaptiveTetCellRefinement.hpp"

namespace Foam::TriSurfaceImmersion {

template<class Surface>
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
    static scalar levelSetValue(const Surface&, const point&);
    
    // Criterion specific functions
    static std::tuple<scalar, label> maxDistSqrAndPointID
    (
        const indexedTet& tet,
        const std::vector<scalar>& signedDistances
    );
};

template<class Surface>
class signCriterion
{
public:
    static bool needsRefinemeent
    (
        const indexedTet& tet,
        const std::vector<point>& points,
        const std::vector<scalar>& levelSetValues
    );
    static scalar levelSetValue(const Surface&, const point&);
};

} // End namespace Foam::TriSurfaceImmersion

#include "refinementCriteriaI.hpp"

#endif