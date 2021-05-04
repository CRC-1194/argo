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

#ifndef surfaceMeshCellApproximationI_H
#define surfaceMeshCellApproximationI_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double surfaceMeshCellApproximation::nTrianglesPerCell() const
{
    return sigDistCalcPtr_->nSurfaceElements()/interfaceCellIDs_.size();
}

label surfaceMeshCellApproximation::nIntersectedCells() const
{
    return interfaceCellIDs_.size();
}

label surfaceMeshCellApproximation::maxRefinementLevel() const
{
    return maxUsedRefinementLevel_;
}

const signedDistanceCalculator& surfaceMeshCellApproximation::sigDistCalc() const
{
    return *sigDistCalcPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //