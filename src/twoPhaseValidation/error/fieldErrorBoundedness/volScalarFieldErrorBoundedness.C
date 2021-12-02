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
    Class template instantiation for the boundedness field error.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "meshMagnitude.H"
#include "volScalarFieldErrorBoundedness.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldErrorsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
    defineTypeNameAndDebug (volScalarFieldErrorBoundedness, 0);
    addToRunTimeSelectionTable(volScalarFieldError, volScalarFieldErrorBoundedness, Scalar);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorBoundedness::volScalarFieldErrorBoundedness(
    scalar errorTolerance,
    scalar lowerBound,
    scalar upperBound
)
:
    fieldError(errorTolerance),
    lowerBound_(lowerBound),
    upperBound_(upperBound)
{}

Foam::volScalarFieldErrorBoundedness::volScalarFieldErrorBoundedness(const volScalarFieldErrorBoundedness& copy)
:
    fieldError(copy)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::volScalarFieldErrorBoundedness::~volScalarFieldErrorBoundedness()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volScalarFieldErrorBoundedness::computeError(
    const volScalarField&,
    const volScalarField& currentField
)
{
    auto Eb = max(max(0,gMax((currentField - upperBound_)())), 
                  max(0,gMax((lowerBound_ - currentField)())));

    auto alphaMax = gMax(currentField); 
    auto alphaMin = gMin(currentField); 

    Info << "Max [" << currentField.name() << "]: " << alphaMax << endl;
    Info << "Min [" << currentField.name() << "]: " << alphaMin << endl;

#ifdef TESTING
    errorField.write(); 
#endif

    this->setErrorValue(Eb);
}

// ************************************************************************* //
