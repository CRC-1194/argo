/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tobias Tolle, TU Darmstadt
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
    Foam::signedDistanceCalculator

Author
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Thermo-Fluids and Interfaces
    TU Darmstadt
    Germany

Description
    Given an oriented triangular surface, compute the signed distance of a given
    point or point cloud to the surface.
    The absolute distance is computed as the Euclidean distance between a given
    point and its closest point on the surface.
    The sign (inside/outside) is determined using normal information
    interpolated from the triangle vertices by a linear shape function.

SourceFiles
    signedDistanceCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef signedDistanceCalculator_H
#define signedDistanceCalculator_H

#include "triSurface.H"
#include "triSurfaceSearch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace SigDistCalc{


/*---------------------------------------------------------------------------*\
                         Class signedDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class signedDistanceCalculator
{
    // Private Data
    const triSurface& surface_;
    triSurfaceSearch surface_search_;
    vectorField vertex_normals_;


    // Private Member Functions
    void compute_vertex_normals();

public:

    // Constructors
    signedDistanceCalculator(const triSurface& surface);


    //- Destructor
    ~signedDistanceCalculator() = default;


    // Member Functions

    // Make point to nearest triangle information available to the caller of this
    // member function
    scalarField signed_distance
                (
                    DynamicList<pointIndexHit>& point_to_nearest_triangle,
                    const pointField& pf,
                    const scalarField& search_dist_sqr,
                    const scalar out_of_search_domain=0.0
                ) const;
    scalarField signed_distance
                (
                    const pointField& pf,
                    const scalarField& search_dist_sqr,
                    const scalar out_of_search_domain=0.0
                ) const;
    std::tuple<pointIndexHit, scalar> signed_distance(const point& p, const scalar search_dist_sqr) const;
    scalar signed_distance(const point& p) const;
    vector normal_at_surface(const pointIndexHit& hit_info) const;

    // Access
    const vectorField& vertex_normals() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
