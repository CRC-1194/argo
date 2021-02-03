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
    Foam::triSurfaceAdapter

Description

SourceFiles
    triSurfaceAdapterI.H
    triSurfaceAdapter.C
    triSurfaceAdapterIO.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceAdapter_H
#define triSurfaceAdapter_H

#include "triSurface.H"
#include "triSurfaceSearch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace PolynomialVof {


/*---------------------------------------------------------------------------*\
                         Class triSurfaceAdapter Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceAdapter
{
    // Private data
    const triSurface& surface_;
    const triSurfaceSearch& search_;
    const vector span_;
    scalar ref_length_;


    // Private Member Functions
    scalar compute_ref_length();

public:


    // Constructors
    triSurfaceAdapter(const triSurface&, const triSurfaceSearch&, vector); 

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
