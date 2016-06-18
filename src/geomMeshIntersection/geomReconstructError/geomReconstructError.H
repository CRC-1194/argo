/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::geomReconstructError

Description
    A class for computing the standard errors of the geometrical two phase
    flow methods (VoF, MoF) with piecewise planar interface description.

    Volume of symmetric difference. 

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

SourceFiles
    geomReconstructError.C

\*---------------------------------------------------------------------------*/

#ifndef geomReconstructError_H
#define geomReconstructError_H

#include "geomMeshIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace GeometricalTransport {

    //- Volume of symmetric difference calculator.
    template<typename MeshIntersection, typename GeomInterface>
    tmp<volScalarField> volSymmDiff(
        MeshIntersection const& intersection, 
        GeomInterface const& interface, 
        const volScalarField& alpha1
    ); 

} // End namespace GeometricalTransport 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "geomReconstructError.C"

#endif

// ************************************************************************* //
