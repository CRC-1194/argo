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


#include "meshMagnitude.H"
#include "volScalarFieldErrorInterfaceSharpness.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldErrorsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
    defineTypeNameAndDebug (volScalarFieldErrorInterfaceSharpness, 0);
    addToRunTimeSelectionTable(volScalarFieldError, volScalarFieldErrorInterfaceSharpness, Scalar);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorInterfaceSharpness::volScalarFieldErrorInterfaceSharpness
(
    scalar errorTolerance
)
:
    fieldError(errorTolerance)
{}

Foam::volScalarFieldErrorInterfaceSharpness::volScalarFieldErrorInterfaceSharpness
(
    const volScalarFieldErrorInterfaceSharpness& copy
)
:
    fieldError(copy)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorInterfaceSharpness::~volScalarFieldErrorInterfaceSharpness()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volScalarFieldErrorInterfaceSharpness::computeError(
    const volScalarField& initialField,
    const volScalarField& currentField
)
{
    const scalarField& M = getMeshMagnitude(initialField);

    //this->setErrorValue(
    //        sum(M * mag(initialField - currentField)) /
    //        sum(M * mag(initialField))
    this->setErrorValue(
            scalar(1.) - (sum(M * mag(currentField - scalar(0.5))) /
            sum(scalar(0.5) * M))
    );
}

// ************************************************************************* //
