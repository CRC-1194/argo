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
    A class template for L1 VoF advection error calculation. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef fieldErrorL1_H
#define fieldErrorL1_H

#include "fieldError.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldErrorL1 Declaration
\*---------------------------------------------------------------------------*/

template<typename GeomField>
class fieldErrorL1
:
    public fieldError<GeomField>
{

    // FIXME: Wrong if GeomField == Surface field.  
    // Use OF standard template parameters for this class. TM
    tmp<volScalarField> EgTmp_; 

public:

    TypeName("E_g");

    // Constructors

        //- Construct from components
        explicit fieldErrorL1(scalar errorTolerance=SMALL);

        //- Construct as copy
        fieldErrorL1(const fieldErrorL1&);

    //- Destructor

        virtual ~fieldErrorL1();

    // Member Functions

        virtual void computeError(
            const GeomField& initialField,
            const GeomField& currentField
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldErrorL1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
