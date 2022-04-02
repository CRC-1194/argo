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

\*---------------------------------------------------------------------------*/

#include "pandoraIsosurfaceCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "error.H"
#include "dimensionedScalarFwd.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "volFieldsFwd.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraIsosurfaceCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraIsosurfaceCurvature, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
pandoraIsosurfaceCurvature::pandoraIsosurfaceCurvature
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pandoraCurvature{mesh, dict},
    fieldName_{curvatureDict_.get<word>("fieldName")},
    nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 3)),
    nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 3))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField& pandoraIsosurfaceCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (!meshDb.found(fieldName_))
    {
        FatalErrorInFunction
            << "pandoraIsosurfaceCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    // Constuct isoSurface. 
    const volScalarField& vf = mesh().lookupObject<volScalarField>(fieldName_);

    volPointInterpolation vpInterp(mesh());

    tmp<pointScalarField> pfTmp = vpInterp.interpolate(vf);
    pointScalarField& pf = pfTmp.ref();

    isoSurfaceParams isoParams
    (
        isoSurfaceParams::algorithmType::ALGO_DEFAULT,
        isoSurfaceParams::filterType::DIAGCELL
    );
    //isoParams.snap(false);
    //isoParams.mergeTol(1e-22);

    isoSurfaceTopo alphaIso(
        mesh(),
        vf,
        pf,
        0.5,
        isoParams
    );
    alphaIso.write("alphaIso.vtk");

    // Get meshCells and normals from isoSurface. 
    const auto& triToCell = alphaIso.meshCells();
    const auto& triNormals = alphaIso.Sf();

    volVectorField interfaceNormals
    (
        IOobject
        (
            "interfaceNormals",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(
            "interfaceNormals",
            dimless,
            vector::zero
        )
    );

    forAll(triToCell, i)
    {
        label cellI = triToCell[i];
        interfaceNormals[cellI] = triNormals[i];
    }

    // Use averaging to propagate the normals into the bulk. 
    volVectorField averagedNormals = interfaceNormals /
    (
        mag(interfaceNormals) +
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );

    // Propagate the normals into the bulk. 
    for (label i = 0; i < nPropagate_; i++)
    {
        averagedNormals = fvc::average(averagedNormals);

        averagedNormals /= mag(averagedNormals) + 
            dimensionedScalar("SMALL", averagedNormals.dimensions(), SMALL);

        // Re-set the smoothed normal vectors in interface cells. 
        forAll(interfaceNormals, cellI)
        {
            if ((interfaceNormals[cellI][0] != 0) ||
                (interfaceNormals[cellI][1] != 0) ||
                (interfaceNormals[cellI][2] != 0))
            {
                averagedNormals[cellI] = interfaceNormals[cellI] / 
                    mag(interfaceNormals[cellI]);
            }
        }
    }

    // Smooth the normals
    for (label i = 0; i < nAverage_; i++)
    {
        averagedNormals = fvc::average(averagedNormals);

        // Normalize the normal vectors. 
        averagedNormals /= mag(averagedNormals) +
            dimensionedScalar("SMALL", averagedNormals.dimensions(), SMALL);
    }

    cellCurvature_ = -fvc::div(averagedNormals);

    return cellCurvature_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
