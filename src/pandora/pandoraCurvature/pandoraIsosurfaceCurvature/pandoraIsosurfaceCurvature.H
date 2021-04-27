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
    Foam::pandoraIsosurfaceCurvature

Description

SourceFiles
    pandoraIsosurfaceCurvatureI.H
    pandoraIsosurfaceCurvature.C
    pandoraIsosurfaceCurvatureIO.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraIsosurfaceCurvature_H
#define pandoraIsosurfaceCurvature_H

#include "isoSurface.H"
#include "pandoraCurvature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class pandoraIsosurfaceCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraIsosurfaceCurvature
:
    public pandoraCurvature 
{
    // Private Data
    const word fieldName_;


    // Private Member Functions


public:

    // Static Data Members
    TypeName ("isosurface");

    // Constructors

        //- Construct null
        pandoraIsosurfaceCurvature();

        //- Construct from components
        pandoraIsosurfaceCurvature(const fvMesh& mesh, const dictionary& dict);

        //- Construct as copy
        pandoraIsosurfaceCurvature(const pandoraIsosurfaceCurvature&) = default;


    //- Destructor
    virtual ~pandoraIsosurfaceCurvature() = default;


    // Member Functions
    virtual volScalarField& cellCurvature(); 


    // Member Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "pandoraIsosurfaceCurvatureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
