/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "pandora.hpp"

#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandora, false);

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pandora::pandora(const fvMesh& mesh)
:
    pandoraDict_(mesh.solutionDict().subDict("pandora")),
    curvPtr_(pandoraCurvature::New(mesh, pandoraDict_)),
    curvRegularisationPtr_(
        pandoraCurvatureRegularisation::New(pandoraDict_)
    ),
    curvExtensionPtr_(
        pandoraCurvatureExtension::New(pandoraDict_)
    ), 
    sigma_
    (
        dimensionedScalar(
            "sigma",
            dimForce / dimLength,
            IOdictionary 
            (
                IOobject
                (
                    "transportProperties",
                    "constant", 
                    mesh,
                    IOobject::MUST_READ, 
                    IOobject::NO_WRITE
                )
            ).get<scalar>("sigma")
        )
    ),
    fSigma_
    (
        IOobject
        (
            "fSigma", 
            mesh.time().timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh, 
        dimensionedScalar("fSigma", dimForce / pow(dimLength,3), 0) 
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const surfaceScalarField& pandora::surfaceTensionForce
(
    const volScalarField& indicator
)
{
    // TODO (TT): move this to an appropriate member function if it works
    volScalarField isInterfaceCell{indicator};
    forAll(isInterfaceCell, I)
    {
        if ((0.1 < isInterfaceCell[I]) && (isInterfaceCell[I] < 0.9))
        {
            isInterfaceCell[I] = 1.0;
        }
        else
        {
            isInterfaceCell[I] = 0.0;
        }
    }
    isInterfaceCell.correctBoundaryConditions();

    volScalarField& cellCurvature = curvPtr_->cellCurvature();

    curvRegularisationPtr_->regularise(cellCurvature, isInterfaceCell); 

    curvExtensionPtr_->extend(cellCurvature, isInterfaceCell); 

    fSigma_ = 
        sigma_ * fvc::interpolate(cellCurvature) * fvc::snGrad(indicator);

    return fSigma_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
