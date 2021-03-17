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

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                         Class signedDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class signedDistanceCalculator
{
    // Private Data
    const triSurface& surface_;
    triSurfaceSearch surfaceSearch_;
    vectorField vertexNormals_;


    // Private Member Functions
    void computeVertexNormals();

public:

    // Constructors
    explicit signedDistanceCalculator(const triSurface& surface);


    // Member Functions

    // Make point to nearest triangle information available to the caller of this
    // member function
    scalarField 
    signedDistance
    (
        DynamicList<pointIndexHit>& pointToNearestTriangle,
        const pointField& pf,
        const scalarField& searchDistSqr,
        scalar outOfSearchDomain=0.0
    ) const;

    scalarField
    signedDistance
    (
        const pointField& pf,
        const scalarField& searchDistSqr,
        scalar outOfSearchDomain=0.0
    ) const;

    std::tuple<pointIndexHit, scalar> signedDistance(const point& p, scalar searchDistSqr) const;
    scalar signedDistance(const point& p) const;
    vector normalAtSurface(const pointIndexHit& hitInfo) const;

    // Access
    const vectorField& vertexNormals() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
