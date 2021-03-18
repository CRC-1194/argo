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
    surfaceMeshCellApproximation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshCellApproximation_H
#define surfaceMeshCellApproximation_H

#include "volFieldsFwd.H"
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

    // Private Data
    label nIntersectedCells_;
    label nTrianglesPerCell_;
    label maxRefinementLevel_;

public:

    // Static Data Members
    TypeName ("SMCA")

    // Generated Methods

//        //- Default construct
//        surfaceMeshCellApproximation() = default;
//
//        //- Copy construct
//        surfaceMeshCellApproximation(const surfaceMeshCellApproximation&) = default;
//
//        //- Copy assignment
//        surfaceMeshCellApproximation& operator=(const surfaceMeshCellApproximation&) = default;


    // Constructors
    surfaceMeshCellApproximation
    (
        const dictionary& configDict,
        const fvMesh& mesh,
        const triSurface& surface
    );


    //- Destructor
    //~surfaceMeshCellApproximation() override = default;


    // Member Functions

    //- Access

    inline const double nTrianglesPerCell() const override;

    inline const label nIntersectedCells() const override;

    inline const label maxRefinementLevel() const override;

    //- Computation

    void calcVolumeFraction(volScalarField& alpha) override; 

    void findIntersectedCells() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceMeshCellApproximationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
