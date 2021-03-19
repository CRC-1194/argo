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

#ifndef volumeFractionCalculatorI_H
#define volumeFractionCalculatorI_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Time& volumeFractionCalculator::time() const
{
    return runTime_;
}

const fvMesh& volumeFractionCalculator::mesh() const
{
    return mesh_;
}

const triSurface& volumeFractionCalculator::surface() const
{
    return surface_;
}

const volScalarField& volumeFractionCalculator::cellSignedDist() const
{
    return cellSignedDist_;
}

const volScalarField& volumeFractionCalculator::cellSignedDist0() const
{
    return cellSignedDist0_;
}

const DynamicList<pointIndexHit>& volumeFractionCalculator::cellNearestTriangle() const
{
    return cellNearestTriangle_;
}

const Foam::pointScalarField& volumeFractionCalculator::pointSignedDist() const
{
    return pointSignedDist_;
}

const DynamicList<pointIndexHit>& volumeFractionCalculator::pointNearestTriangle() const
{
    return pointNearestTriangle_;
}

const signedDistanceCalculator& volumeFractionCalculator::signedDistCalc() const
{
    return sigDistCalc_;
}

const searchDistanceCalculator& volumeFractionCalculator::searchDistCalc() const
{
    return searchDistCalc_;
}

bool volumeFractionCalculator::writeGeometry() const
{
    return writeGeometry_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
