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

    A 3D divergence-free shear field model introduced in [1] 

    [1] Liovic, P., Rudman, M., Liow, J. L., Lakehal, D., & Kothe, D. (2006). A
    3D unsplit-advection volume tracking algorithm with planarity-preserving
    interface reconstruction. Computers and Fluids, 35(10), 1011–1032.
    http://doi.org/10.1016/j.compfluid.2005.09.003

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef shearFieldModel_H
#define shearFieldModel_H

#include "shear2DFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   3D shear velocity function object
\*---------------------------------------------------------------------------*/

class shearFieldModel
    :
        public shear2DFieldModel
{

    scalar R_;
    scalar x0_;
    scalar y0_;

public:

    //- Runtime type information
    TypeName("shear");

    // Constructors

    shearFieldModel
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
