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
    surfaceMeshCellIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshCellIntersection_H
#define surfaceMeshCellIntersection_H

#include "volFieldsFwd.H"
#include "volumeFractionCalculator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                 Class surfaceMeshCellIntersection Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshCellIntersection
:
    public volumeFractionCalculator
{
private:

    // Private Data


public:

    // Static Data Members
    TypeName ("SMCI")

    // Generated Methods

//        //- Default construct
//        surfaceMeshCellIntersection() = default;
//
//        //- Copy construct
//        surfaceMeshCellIntersection(const surfaceMeshCellIntersection&) = default;
//
//        //- Copy assignment
//        surfaceMeshCellIntersection& operator=(const surfaceMeshCellIntersection&) = default;


    // Constructors
    explicit surfaceMeshCellIntersection(const dictionary& configDict);


    //- Destructor
    //~surfaceMeshCellIntersection() override = default;


    // Member Functions
    void calcVolumeFraction(volScalarField& alpha) override; 

    // Access

    // Check

    // Edit

    // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "surfaceMeshCellIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
