/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2015 Daniel Deising 
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
    Sharpness volume field error.  

Author
    Daniel Deising deising@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef volScalarFieldErrorInterfaceSharpness3_H
#define volScalarFieldErrorInterfaceSharpness3_H

#include "fieldError.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class volScalarFieldErrorInterfaceSharpness3 Declaration
\*---------------------------------------------------------------------------*/

class volScalarFieldErrorInterfaceSharpness3
:
    public fieldError<volScalarField>
{
public:

    TypeName("fieldErrorInterfaceSharpness3");

    // Constructors

        //- Construct from components
        explicit volScalarFieldErrorInterfaceSharpness3(scalar errorTolerance=SMALL);

        //- Construct as copy
        volScalarFieldErrorInterfaceSharpness3(const volScalarFieldErrorInterfaceSharpness3&);

    //- Destructor

        virtual ~volScalarFieldErrorInterfaceSharpness3();

    // Member Functions

        virtual void computeError(
            const volScalarField& initialField,
            const volScalarField& currentField
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "fieldErrorInterfaceSharpness3.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
