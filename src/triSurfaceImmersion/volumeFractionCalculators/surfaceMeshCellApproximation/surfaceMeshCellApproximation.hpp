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
    Foam::TriSurfaceImmersion::surfaceMeshCellApproximation

Description
    Implements the surface-mesh-cell-approximation algorithm (SMCA) using 
    tetrahedral refinement and tet volume fraction approximation
    to compute volume fractions of interface cells.

SourceFiles
    surfaceMeshCellApproximation.cpp
    surfaceMeshCellApproximationI.hpp

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshCellApproximation_H
#define surfaceMeshCellApproximation_H

#include <vector>

#include "signedDistanceCalculator.hpp"
#include "volFieldsFwd.H"

#include "AdaptiveTetCellRefinement.hpp"
#include "volumeFractionCalculator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\ 
               Class surfaceMeshCellApproximation Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshCellApproximation
:
    public volumeFractionCalculator
{
private:

    // Private types
    using cellDecompositionTuple =
        std::tuple<std::vector<indexedTet>, std::vector<point>, std::vector<scalar>>;
    struct searchSphere
    {
        vector centre;
        scalar radiusSquared;
    };

    // Private data
    autoPtr<signedDistanceCalculator> sigDistCalcPtr_;
    std::vector<label> interfaceCellIDs_;
    label maxAllowedRefinementLevel_;
    label maxUsedRefinementLevel_ = 0;
    
    // Member functions
    bool intersectionPossible(label cellID) const;
    cellDecompositionTuple decomposeCell(label cellID) const;
    label nTets(label cellID) const;
    //searchSphere cellInterfaceSearchSphere(label cellID) const;
    //triSurface surfaceSubset(vectorField& vertexNormals, label cellID) const;


public:

    // Static Data Members
    TypeName ("SMCA")


    // Constructors
    surfaceMeshCellApproximation
    (
        const dictionary& configDict,
        const fvMesh& mesh
    );


    // Member Functions

    //- Access

    inline double nTrianglesPerCell() const override;

    inline label nIntersectedCells() const override;

    inline label maxRefinementLevel() const override;

    const signedDistanceCalculator& sigDistCalc() const override;

    //- Computation

    void calcVolumeFraction(volScalarField& alpha) override; 

    void findIntersectedCells() override;

    //- Write

    void writeFields() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMeshCellApproximationI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
