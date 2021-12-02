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
    3D deformation divergence free field model defined in [1, example 11.2] 

    [1] R. LeVeque, High-resolution conservative algorithms for advection in
    incompressible flow, SIAM J. Numer. Anal. 33, 627 (1996).

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef deformationFieldModel_H
#define deformationFieldModel_H

#include "divFreeFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Pseudo 2D shear velocity function object
\*---------------------------------------------------------------------------*/

class deformationFieldModel
    :
        public divFreeFieldModel
{

public:

    //- Runtime type information
    TypeName("deformation");

    // Constructors

    //- Construct for given objectRegistry and dictionary.
    deformationFieldModel
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
