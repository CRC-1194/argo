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
#include "processorFvPatch.H"

#include "isoSurfaceTopo.H"
#include "volPointInterpolation.H"

#include "reconstructionSchemes.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandora, false);

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void pandora::updateInterfaceCells(const volScalarField& indicator)
{
    // isInterfaceCells already up-to-date
    if (lastUpdate_ == indicator.mesh().time().timeIndex())
    {
        return;
    }

    /*
     *  Central idea for determining whether a cell is an interface cell
     *  (a.k.a. a cell intersected by the fluid interface) is the following:
     *  assuming the interface is located at points for which the indicator value
     *  is 0.5, we need to find those cells which contain such points. Here we
     *  use the cell connection across faces to find such points. Although
     *  this does not guarantee to detect all interface cells, this approach
     *  relies on the available parallelization provided by OpenFOAM and
     *  detection may only fail for underresolved interfaces (TT).
     
    // Detect via face based indicator values
    auto indicatorFaceTmp = fvc::interpolate(indicator);
    const auto& alphaf= indicatorFaceTmp.cref();

    const auto& owner = indicator.mesh().owner();
    const auto& neighbour = indicator.mesh().neighbour();

    forAll(alphaf, fid)
    {
        if (mag(alphaf[fid] - 0.5) < SMALL)
        {
            isInterfaceCell_[owner[fid]] = 1.0;
            isInterfaceCell_[neighbour[fid]] = 1.0;
            continue;
        }

        if ((alphaf[fid] < 0.5 && 0.5 < indicator[owner[fid]]) || 
            (alphaf[fid] > 0.5 && 0.5 > indicator[owner[fid]]))
        {
            isInterfaceCell_[owner[fid]] = 1.0;
        }

        if ((alphaf[fid] < 0.5 && 0.5 < indicator[neighbour[fid]]) || 
            (alphaf[fid] > 0.5 && 0.5 > indicator[neighbour[fid]]))
        {
            isInterfaceCell_[neighbour[fid]] = 1.0;
        }
    }

    const auto& meshBoundary = alphaf.boundaryField();

    for (const auto& alphab : meshBoundary)
    {
        if (isA<processorFvPatch>(alphab.patch()))
        {
            const auto& faceToCell = alphab.patch().faceCells();

            forAll(alphab, I)
            {
                if ((alphab[I] < 0.5 && 0.5 < indicator[faceToCell[I]]) || 
                    (alphab[I] > 0.5 && 0.5 > indicator[faceToCell[I]]))
                {
                    isInterfaceCell_[faceToCell[I]] = 1.0;
                }
            }
        }
    }

    // Detect via indicator threshold
    forAll(isInterfaceCell_, cid)
    {
        const scalar eps = 0.01;
        if (eps < indicator[cid] && indicator[cid] < 1.0-eps)
        {
            isInterfaceCell_[cid] = 1.0;
        }
    }
    isInterfaceCell_.correctBoundaryConditions();

    // Close potential gaps betwwen interface cells. Some of the 
    // curvature regularisation classes need a proper connectivity, meaning
    // that interface cells share a face.
    auto iCellFaceTmp = fvc::interpolate(isInterfaceCell_);
    const auto& iCellFace = iCellFaceTmp.cref();
    volScalarField faceCount{"faceCount", isInterfaceCell_};
    faceCount = 0.0;

    forAll(iCellFace, fid)
    {
        if (iCellFace[fid] > 0.0)
        {
            faceCount[owner[fid]] += 1.0;
            faceCount[neighbour[fid]] += 1.0;
        }
    }

    const auto& iCellBoundaries = iCellFace.boundaryField();

    for (const auto& iCellb: iCellBoundaries)
    {
        if (isA<processorFvPatch>(iCellb.patch()))
        {
            const auto& faceToCell = iCellb.patch().faceCells();

            forAll(iCellb, I)
            {
                if (iCellb[I] > 0.0)
                {
                    faceCount[faceToCell[I]] += 1.0;
                }
            }
        }
    }

    forAll(isInterfaceCell_, cid)
    {
        const scalar smallEps = 1.0e-8;
        if ((faceCount[cid] >= 2.0) && (smallEps < indicator[cid]) && (indicator[cid] < (1.0-smallEps)))
        {
            isInterfaceCell_[cid] = 1.0;
        }
    }
    isInterfaceCell_.correctBoundaryConditions();
     */

    lastUpdate_ = indicator.mesh().time().timeIndex();

    reconstructionSchemes* surf = 
        indicator.mesh().getObjectPtr<reconstructionSchemes>("reconstructionScheme");

    const boolList& ics = surf->interfaceCell();

    forAll (ics, i)
    {
        if (ics[i])
            isInterfaceCell_[i] = 1.0;
        else
            isInterfaceCell_[i] = 0.0;
    }

    isInterfaceCell_.correctBoundaryConditions();
    /*
    */
}


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
    ),
    isInterfaceCell_{
        IOobject
        (
            "isInterfaceCell",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar{"zero", dimless, 0}
    },
    lastUpdate_{-1}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const surfaceScalarField& pandora::surfaceTensionForce
(
    const volScalarField& indicator
)
{
    volScalarField& cellCurvature = curvPtr_->cellCurvature();

    updateInterfaceCells(indicator);

    curvRegularisationPtr_->regularise(cellCurvature, isInterfaceCell_); 

    curvExtensionPtr_->extend(cellCurvature, isInterfaceCell_); 

    fSigma_ = 
        sigma_ * fvc::interpolate(cellCurvature) * fvc::snGrad(indicator);

    return fSigma_;
}


const volScalarField& pandora::isInterfaceCell
(
    const volScalarField& indicator
)
{
    updateInterfaceCells(indicator);

    return isInterfaceCell_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
