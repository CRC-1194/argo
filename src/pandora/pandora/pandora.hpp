/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tomislav Maric, Tobias Tolle, TU Darmstadt
                       Anja Lippert, BOSCH CR 
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
    Foam::pandora

Description
    Curvature and surface tension model library.  

Usage
     Currently linked with the interIsoPandoraFoam solver and only available
     with this solver. 

     In system/fvSolution, use 

     \verbatim
        pandora
        {
            curvature
            {
                type divNormal;
                nAverages 10; 
                normalField "recon::normal"; 
            }
            curvatureExtension
            {
                type noExtension;
            }
            curvatureRegularisation
            {
                type noRegularisation;
            }
        }
     \endverbatim

SourceFiles
    pandoraI.H
    pandora.C
    pandoraIO.C

\*---------------------------------------------------------------------------*/

#ifndef pandora_H
#define pandora_H

#include "fvMesh.H"
#include "surfaceFields.H"
#include "typeInfo.H"

#include "pandoraCurvature.hpp"
#include "pandoraCurvatureExtension.hpp"
#include "pandoraCurvatureRegularisation.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandora Declaration
\*---------------------------------------------------------------------------*/

class pandora
{
    // Private Data

    const dictionary& pandoraDict_; 

    autoPtr<pandoraCurvature> curvPtr_; 
    autoPtr<pandoraCurvatureRegularisation> curvRegularisationPtr_;
    autoPtr<pandoraCurvatureExtension> curvExtensionPtr_; 

    dimensionedScalar sigma_; 

    surfaceScalarField fSigma_;

public:

    // Static Data Members

    TypeName ("pandora");

    // Constructors

        //- Construct from components
        pandora(const fvMesh& mesh);

        //- Construct as copy
        pandora(const pandora&) = default;

    //- Destructor
    virtual ~pandora() = default;

    // Member Functions

    const surfaceScalarField& surfaceTensionForce
    (
        const volScalarField& indicator
    ); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
