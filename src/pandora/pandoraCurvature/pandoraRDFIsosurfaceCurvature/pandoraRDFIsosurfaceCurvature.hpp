/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::pandoraRDFIsosurfaceCurvature

Description

SourceFiles
    pandoraRDFIsosurfaceCurvatureI.H
    pandoraRDFIsosurfaceCurvature.C
    pandoraRDFIsosurfaceCurvatureIO.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraRDFIsosurfaceCurvature_H
#define pandoraRDFIsosurfaceCurvature_H

#include "isoSurfacePoint.H"
#include "pandoraCurvature.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class pandoraRDFIsosurfaceCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraRDFIsosurfaceCurvature
:
    public pandoraCurvature 
{
    // Private Data
    const word fieldName_;
    const label nPropagate_;
    const label nAverage_;


    // Private Member Functions


public:

    // Static Data Members
    TypeName ("RDFisosurface");

    // Constructors

        //- Construct null
        pandoraRDFIsosurfaceCurvature();

        //- Construct from components
        pandoraRDFIsosurfaceCurvature(const fvMesh& mesh, const dictionary& dict);

        //- Construct as copy
        pandoraRDFIsosurfaceCurvature(const pandoraRDFIsosurfaceCurvature&) = default;


    //- Destructor
    virtual ~pandoraRDFIsosurfaceCurvature() = default;


    // Member Functions
    virtual volScalarField& cellCurvature(); 


    // Member Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "pandoraRDFIsosurfaceCurvatureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
