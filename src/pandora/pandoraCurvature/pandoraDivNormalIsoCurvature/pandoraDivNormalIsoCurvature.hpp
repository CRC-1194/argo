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
    Foam::pandoraDivNormalIsoCurvature

Description

SourceFiles
    pandoraDivNormalIsoCurvatureI.H
    pandoraDivNormalIsoCurvature.C
    pandoraDivNormalIsoCurvatureIO.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraDivNormalIsoCurvature_H
#define pandoraDivNormalIsoCurvature_H

#include "isoSurfacePoint.H"
#include "isoSurfaceTopo.H"
#include "pandoraCurvature.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class pandoraDivNormalIsoCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraDivNormalIsoCurvature
:
    public pandoraCurvature 
{
    // Private Data
    const word fieldName_;
    const label nPropagate_;
    const label nAverage_;
    const scalar range_;
    const label nLayer_;
    volVectorField normals_;


    // Private Member Functions
    void calcKappa 
    (
        const isoSurfaceTopo& isoTopo,
        const scalar& iso,
        const scalar& delta_x,
        vectorField& normal,
        scalarField& isoList,
        scalarField& kappa,
        scalarField& scaledKappa,
        scalarField& area,
        labelList& marker
    );

    vector largestAreaNormal
    (
        const vectorField& normalList,
        const scalarList& isoList,
        const scalarList& kappaList,
        const scalarList& scaledKappaList,
        const scalarList& areaList,
        scalar& maxAreaIso,
        scalar& maxAreaKappa,
        scalar& maxAreaScaledKappa
    );


public:

    // Static Data Members
    TypeName ("divNormalIso");

    // Constructors

        //- Construct null
        pandoraDivNormalIsoCurvature();

        //- Construct from components
        pandoraDivNormalIsoCurvature(const fvMesh& mesh, const dictionary& dict);

        //- Construct as copy
        pandoraDivNormalIsoCurvature(const pandoraDivNormalIsoCurvature&) = default;


    //- Destructor
    virtual ~pandoraDivNormalIsoCurvature() = default;


    // Member Functions
    virtual volScalarField& cellCurvature(); 


    // Member Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "pandoraDivNormalIsoCurvatureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
