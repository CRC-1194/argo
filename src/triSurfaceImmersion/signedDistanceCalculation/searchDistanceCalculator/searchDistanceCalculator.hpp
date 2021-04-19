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
    Foam::searchDistanceCalculator

Description
    Compute search distances at cell centres and cell corner points of a
    mesh. These distances are used to define search spheres for octree based
    search operations.

SourceFiles
    searchDistanceCalculatorI.H
    searchDistanceCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef searchDistanceCalculator_H
#define searchDistanceCalculator_H

#include "pointFields.H"
#include "pointFieldsFwd.H"
#include "volFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                         Class searchDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class searchDistanceCalculator
{
private:
    
    //- Reference to the volume mesh
    const fvMesh& mesh_;

    //- Set width of the narrow band in number of cells
    const scalar searchDistFactor_;

    //- Squared search distance field in cell centers. 
    volScalarField cellSqrSearchDist_;

    //- Point mesh constituted by cell corner points
    pointMesh pMesh_;

    //- Inverse Distance Interpolation : cell centers to cell corners. 
    volPointInterpolation cellsToPointsInterp_;

    //- Squared search distance at cell corner points. 
    pointScalarField pointSqrSearchDist_;  

public:

    // Constructors
    explicit searchDistanceCalculator(const fvMesh& mesh, scalar searchDistFactor=4.0);

    // Member Functions

    //- Access
    inline const fvMesh& mesh() const;
    inline scalar searchDistFactor() const;
    inline const volScalarField& cellSqrSearchDist() const;
    inline const pointScalarField& pointSqrSearchDist() const;

    //- Write
    void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "searchDistanceCalculatorI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
