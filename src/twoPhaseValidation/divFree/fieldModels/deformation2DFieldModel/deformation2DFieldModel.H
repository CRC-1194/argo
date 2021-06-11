/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2011 Tomislav Maric 
     \\/     M anipulation  |
-------------------------------------------------------------------------------

License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

    2D deformation divergence free velocity field described with a stream
    function in [1, equation 17] 

    [1] Rider, W. J., & Kothe, D. B. (1998). Reconstructing Volume Tracking.
    Journal of Computational Physics, 141(2), 112–152.
    http://doi.org/10.1006/jcph.1998.5906

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 

    Mathematical Modeling and Analysis
    Center of Smart Interfaces
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef deformation2DFieldModel_H
#define deformation2DFieldModel_H

#include "divFreeFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Pseudo 2D shear velocity function object
\*---------------------------------------------------------------------------*/

class deformation2DFieldModel
    :
        public divFreeFieldModel
{

public:

    //- Runtime type information
    TypeName("deformation2D");

    // Constructors
    deformation2DFieldModel
    (
        const Time& time,
        const dictionary& dict
    );

    // Member Functions

    // Compute velocity at point X in time t.
    vector velocity(point X, scalar t) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
