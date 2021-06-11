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
    Class template for the normalized L1 advection error calculation. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/



#ifndef fieldErrorsVolumeFwd_H
#define fieldErrorsVolumeFwd_H

#include "fieldErrorsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    template<typename GeomField> class fieldErrorL1normalized;

    typedef fieldErrorL1normalized<volScalarField> volScalarFieldErrorL1normalized;
    typedef fieldErrorL1normalized<volVectorField> volVectorFieldErrorL1normalized;
    typedef fieldErrorL1normalized<volTensorField> volTensorFieldErrorL1normalized;
    typedef fieldErrorL1normalized<volSymmTensorField> volSymmTensorFieldErrorL1normalized;

    typedef fieldErrorL1normalized<surfaceScalarField> surfaceScalarFieldErrorL1normalized;
    typedef fieldErrorL1normalized<surfaceVectorField> surfaceVectorFieldErrorL1normalized;
    typedef fieldErrorL1normalized<surfaceTensorField> surfaceTensorFieldErrorL1normalized;
    typedef fieldErrorL1normalized<surfaceSymmTensorField> surfaceSymmTensorFieldErrorL1normalized;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
