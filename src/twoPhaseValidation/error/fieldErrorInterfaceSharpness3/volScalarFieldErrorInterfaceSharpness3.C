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
#include "volScalarFieldErrorInterfaceSharpness3.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldErrorsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
    defineTypeNameAndDebug (volScalarFieldErrorInterfaceSharpness3, 0);
    addToRunTimeSelectionTable(volScalarFieldError, volScalarFieldErrorInterfaceSharpness3, Scalar);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorInterfaceSharpness3::volScalarFieldErrorInterfaceSharpness3
(
    scalar errorTolerance
)
:
    fieldError(errorTolerance)
{}

Foam::volScalarFieldErrorInterfaceSharpness3::volScalarFieldErrorInterfaceSharpness3
(
    const volScalarFieldErrorInterfaceSharpness3& copy
)
:
    fieldError(copy)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorInterfaceSharpness3::~volScalarFieldErrorInterfaceSharpness3()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volScalarFieldErrorInterfaceSharpness3::computeError(
    const volScalarField& initialField,
    const volScalarField& currentField
)
{
    const scalarField& M = getMeshMagnitude(initialField);

    scalar tol = 0.001; //ToDo: make it an optional parameter at function call!
    volScalarField interfaceMarkerInit{pos(initialField - tol) * pos(scalar(1.) - tol - initialField)};
    volScalarField interfaceMarkerCurr{pos(currentField - tol) * pos(scalar(1.) - tol - currentField)};
    
    this->setErrorValue(
            sum(M * (interfaceMarkerInit - interfaceMarkerCurr)) /
            sum(interfaceMarkerInit * M)
            //mag(sum(M * (interfaceMarkerInit - interfaceMarkerCurr))) /
            //sum(interfaceMarkerInit * M)
    );
}

// ************************************************************************* //
