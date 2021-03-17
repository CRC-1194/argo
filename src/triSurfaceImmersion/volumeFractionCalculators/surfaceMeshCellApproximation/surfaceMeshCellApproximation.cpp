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

#include "surfaceMeshCellApproximation.hpp"
#include "volumeFractionCalculator.hpp"

#include "addToRunTimeSelectionTable.H"

namespace Foam::TriSurfaceImmersion {

    defineTypeNameAndDebug(surfaceMeshCellApproximation, 0);
    addToRunTimeSelectionTable(volumeFractionCalculator, surfaceMeshCellApproximation, Dictionary);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
surfaceMeshCellApproximation::surfaceMeshCellApproximation
(
    const dictionary& configDict,
    const fvMesh& mesh,
    const triSurface& surface
)
:
    volumeFractionCalculator{configDict, mesh, surface}
    {}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void surfaceMeshCellApproximation::calcVolumeFraction(volScalarField& alpha)
{
}

void surfaceMeshCellApproximation::findIntersectedCells()
{
}

}  // namespace Foam::TriSurfaceImmersion

// ************************************************************************* //