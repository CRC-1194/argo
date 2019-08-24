/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
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
    Foam::orientedSurface

Description

SourceFiles
    orientedSurfaceI.H
    orientedSurface.C
    orientedSurfaceIO.C

\*---------------------------------------------------------------------------*/

#ifndef orientedSurface_H
#define orientedSurface_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace PolynomialVof {


/*---------------------------------------------------------------------------*\
                         Class orientedSurface Declaration
\*---------------------------------------------------------------------------*/

class orientedSurface
{
    // Private data
    point refPoint_;
    vector unitNormal_;
    scalar ref_length_; // This is required for the automatic determination of
                        // the refinement level (TT)

    scalar distanceOrigin_;


    // Private Member Functions
    void updateDistanceToOrigin();


public:

    // Constructors
    orientedSurface(const point& refPoint, const vector& unitNormal, scalar refLength);


    // Member Functions
    scalar signedDistance(const point& trialPoint) const;
    scalar referenceLength() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PolynomialVof

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
