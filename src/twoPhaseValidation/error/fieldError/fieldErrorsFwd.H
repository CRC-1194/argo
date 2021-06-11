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
    Class template for field error calculation.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef fieldErrorsFwd_H
#define fieldErrorsFwd_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fieldError.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    typedef fieldError<volScalarField> volScalarFieldError;
    typedef fieldError<volVectorField> volVectorFieldError;
    typedef fieldError<volTensorField> volTensorFieldError;
    typedef fieldError<volSymmTensorField> volSymmTensorFieldError;

    typedef fieldError<surfaceScalarField> surfaceScalarFieldError;
    typedef fieldError<surfaceVectorField> surfaceVectorFieldError;
    typedef fieldError<surfaceTensorField> surfaceTensorFieldError;
    typedef fieldError<surfaceSymmTensorField> surfaceSymmTensorFieldError;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
