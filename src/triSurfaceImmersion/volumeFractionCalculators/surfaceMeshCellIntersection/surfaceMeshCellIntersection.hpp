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
    Foam::TriSurfaceImmersion::surfaceMeshCellIntersection

Description
    Implements the surface-mesh-cell-intersection algorithm (SMCI)
    using geometrical intersections between surface triangles and cells
    to compute volume fractions of interface cells.

SourceFiles
    surfaceMeshCellIntersection.cpp
    surfaceMeshCellIntersectionI.hpp

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshCellIntersection_H
#define surfaceMeshCellIntersection_H

#include "triSurfaceDistCalc.hpp"
#include "volumeFractionCalculator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{
/*---------------------------------------------------------------------------*\
                 Class surfaceMeshCellIntersection Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshCellIntersection : public volumeFractionCalculator
{
private:

    // Private Data
    //- Reference to triSurface signed distance calculator
    triSurfaceDistCalc sigDistCalc_;

    using dynamicLabelList = DynamicList<label>;
    //- Stores IDs of all cells intersected by the interface
    dynamicLabelList intersectedCellLabels_;

    //- Average number of triangles per interface cell
    label nTrianglesPerCell_;

    // Private Member Functions
    //- Compute volume fractions for intersected cells
    void interfaceCellVolumeFraction(volScalarField& alpha);

public:

    // Static Data Members
    TypeName("SMCI");


    // Constructors
    surfaceMeshCellIntersection(
        const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    //- Return average number of triangles per cell
    inline scalar nTrianglesPerCell() const override;

    //- Return number of cells intersected by the interface
    inline label nIntersectedCells() const override;

    //- Return the maximum refinement level used for tetrahedral refinement
    //  Note: the SMCI algorithm does not use tetrahedral refinement yet.
    inline label maxRefinementLevel() const override;

    //- Return reference to the signed distance calculator.
    //  For the SMCI algorithm, this is guaranteed to be a
    //  triSurface signed distance calculator.
    const signedDistanceCalculator& sigDistCalc() const override;

    //- Compute the volume fraction field alpha
    void calcVolumeFraction(volScalarField& alpha) override;

    //- Find all cells intersected by the interface
    void findIntersectedCells() override;

    //- Write additional fields used during the computation of volume fractions
    void writeFields() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMeshCellIntersectionI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
