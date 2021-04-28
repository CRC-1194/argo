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

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "dimensionedScalarFwd.H"
#include "isoSurface.H"
#include "pandoraIsosurfaceCurvature.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "volFieldsFwd.H"

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
    fieldName_{dict.get<word>("fieldName")}
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

    const volScalarField& vf = mesh().lookupObject<volScalarField>(fieldName_);
    
    pointMesh pmesh{mesh()}; 

    pointScalarField pf 
    (
        IOobject
        (
            "pf", 
            mesh().time().timeName(), 
            mesh(), 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pmesh,
        dimensionedScalar{"zero", dimless, 0.0}
    );

    isoSurface alphaIso{vf, pf, 0.5, false};

    const auto& triToCell = alphaIso.meshCells();
    const auto& triNormals = alphaIso.Sf();
        

    return cellCurvature_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
