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
    Foam::pandoraIntegralCurvature

Description

SourceFiles
    pandoraIntegralCurvatureI.H
    pandoraIntegralCurvature.C
    pandoraIntegralCurvatureIO.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraIntegralCurvature_H
#define pandoraIntegralCurvature_H

#include "isoSurfacePoint.H"
#include "isoSurfaceTopo.H"
#include "pandoraCurvature.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class pandoraIntegralCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraIntegralCurvature
:
    public pandoraCurvature 
{
    // Private Data
    const word fieldName_;
    const label nPropagate_;
    const label nAverage_;
    const scalar range_;
    const label nLayer_;


    // Private Member Functions
    void calcKappa 
    (
        const isoSurfaceTopo& isoTopo,
        const volScalarField& rdf,
        const scalar& iso,
        const scalar& delta_x,
        scalarField& kappa,
        scalarField& area,
        labelList& marker
    );

    scalar largestAreaCurv
    (
        const scalarList& kappaList,
        const scalarList& areaList
    );


public:

    // Static Data Members
    TypeName ("integral");

    // Constructors

        //- Construct null
        pandoraIntegralCurvature();

        //- Construct from components
        pandoraIntegralCurvature(const fvMesh& mesh, const dictionary& dict);

        //- Construct as copy
        pandoraIntegralCurvature(const pandoraIntegralCurvature&) = default;


    //- Destructor
    virtual ~pandoraIntegralCurvature() = default;


    // Member Functions
    virtual volScalarField& cellCurvature(); 


    // Member Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "pandoraIntegralCurvatureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
