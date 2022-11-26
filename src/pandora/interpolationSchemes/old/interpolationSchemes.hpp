/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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
    Foam::interpolationSchemes

Description

SourceFiles
    interpolationSchemesI.H
    interpolationSchemes.C
    interpolationSchemesIO.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationSchemes_H
#define interpolationSchemes_H

#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class interpolationSchemes Declaration
\*---------------------------------------------------------------------------*/

class interpolationSchemes
{
    // Private Data

        //- Reference to the fvMesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        interpolationSchemes(const interpolationSchemes&) = delete;

        //- No copy assignment
        void operator=(const interpolationSchemes&) = delete;

        //- Calculate 3D least square gradient
        vector gradLeastSquare
        (
            const vector& p, 
            const scalar& alphap,
            const List<vector>& c,
            const scalarList& alphac
        );

        //- 3D inverse distance weighting interpolation
        scalar inverseDistanceInterpolate
        (
            const vector& p, 
            const List<vector>& c,
            const scalarList& psi,
            scalar r
        );

        //- 3D high order inverse distance weighting interpolation
        scalar interpolateSecondOrder
        (
            const vector& p, 
            const List<vector>& c,
            const scalarList& psi,
            scalar r
        );

public:

    // Constructors

        //- Construct from components
        interpolationSchemes(const fvMesh&);


    //- Destructor
    ~interpolationSchemes();


    // Member Functions

        //- IDW interpolation function
        scalar IDWinterp
        (
            const vector& p,
            const List<vector>& c,
            const scalarList& psi,
            scalar r = 1.0
        );

        //- IDeC interpolation function
        scalar IDeCinterp
        (
            const vector& p,
            const List<vector>& c,
            const scalarList& psi,
            scalar r = 1.0
        );
        
        //- Least square interpolation function
        scalar LSfitting
        (
            const vector& p,
            const List<vector>& c,
            const scalarList& psi,
            scalar r = 1.0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
