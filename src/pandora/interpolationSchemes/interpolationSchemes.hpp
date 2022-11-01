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

        //- Get the dimension for 2D case
        Vector2D<label> getDimension();

        //- Convert a 3D vector to 2D
        void convert
        (
            vector2D & v2d, 
            const vector & v3d,  
            const Vector2D<label> & index
        );

        //- Convert a 3D vector field to 2D
        void convert
        (
            Field<vector2D> & v2d, 
            const vectorField & v3d,  
            const Vector2D<label> & index
        );

        //- Convert a 2D vector to 3D
        void convert
        (
            vector & v3d, 
            const vector2D & v2d,  
            const Vector2D<label> & index
        );

        //- Convert a 2D vector field to 3D
        void convert
        (
            vectorField & v3d,
            const Field<vector2D> & v2d, 
            const Vector2D<label> & index
        );

        //- Calculate 3D least square gradient
        vector gradLeastSquare
        (
            const vector & p, 
            const scalar & alphap,
            const Field<vector> & c,
            const scalarField & alphac
        );

        //- Calculate 2D least square gradient
        vector2D gradLeastSquare
        (
            const vector2D & p, 
            const scalar & alphap,
            const Field<vector2D> & c,
            const scalarField & alphac
        );

        //- 3D inverse distance weighting interpolation
        scalar inverseDistanceInterpolate
        (
            const vector & p, 
            const vectorField & c,
            const scalarField & psi,
            const scalar & r 
        );

        scalar inverseDistanceInterpolate
        (
            const vector & p, 
            const vector & n,
            const vectorField & c,
            const scalarField & psi,
            const scalar & r 
        );

        //- 2D inverse distance weighting interpolation
        scalar inverseDistanceInterpolate
        (
            const vector2D & p, 
            const Field<vector2D> & c,
            const scalarField & psi,
            const scalar & r 
        );

        //- 3D high order inverse distance weighting interpolation
        scalar interpolateSecondOrder
        (
            const vector & p, 
            const vectorField & c,
            const scalarField & psi,
            const scalar & r 
        );

        scalar interpolateSecondOrder
        (
            const vector & p, 
            const vector & n,
            const vectorField & c,
            const scalarField & psi,
            const scalar & r 
        );

        //- 2D high order inverse distance weighting interpolation
        scalar interpolateSecondOrder
        (
            const vector2D & p, 
            const Field<vector2D> & c,
            const scalarField & psi,
            const scalar & r 
        );


public:

    // Constructors

        //- Construct from components
        interpolationSchemes
        (
            const fvMesh& mesh
        );


    //- Destructor
    ~interpolationSchemes();


    // Member Functions

        //- IDW interpolation function
        scalar IDWinterpolate
        (
            const vector& p,
            const vectorField& c,
            const scalarField& psi,
            const label& r
        );

        //- IDeC interpolation function
        scalar IDeCinterpolate
        (
            const vector& p,
            const vectorField& c,
            const scalarField& psi,
            const label& r
        );
        
        scalar IDeCinterpolate
        (
            const vector& p,
            const vector& n,
            const vectorField& c,
            const scalarField& psi,
            const label& r
        );
        
        //- Least square interpolation function
        scalar LSinterpolate
        (
            const vector& p,
            const vectorField& c,
            const scalarField& psi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
