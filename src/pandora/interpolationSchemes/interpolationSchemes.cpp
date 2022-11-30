/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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

#include "interpolationSchemes.hpp"

#include "leastSquareGrad.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::vector Foam::interpolationSchemes::gradLeastSquare
(
    const vector& p, 
    const scalar& alphap,
    const List<vector>& c,
    const scalarList& alphac
)
{
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());

    DynamicList<vector> centres(100);
    DynamicList<scalar> values(100);

    forAll (c, i)
    {
        centres.append(c[i] - p);
        values.append(alphac[i] - alphap);
    }

    return lsGrad.grad(centres, values);
}

Foam::scalar Foam::interpolationSchemes::inverseDistanceInterpolate
(
    const vector& p, 
    const List<vector>& c,
    const scalarList& psi,
    scalar r
)
{
    scalar w = 0; 
    scalar sum = 0;
    forAll(psi, i)
    {
        scalar a = Foam::pow(mag(c[i] - p), r);
        scalar d = 1.0 / max(a, SMALL);
        sum += psi[i]*d; 
        w += d; 
    }

    return sum/w;
}

Foam::scalar Foam::interpolationSchemes::interpolateSecondOrder 
(
    const vector& p, 
    const List<vector>& c, 
    const scalarList& psi, 
    scalar r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, c, psi, r); 

    if (mag(psi_p) < 1e-8)
        return psi_p;

    scalar psi_p0 = psi_p; 
    scalar psi_p2 = psi_p + 1; 

    label loop = 0;
    scalar delta = fabs((psi_p - psi_p2) / (psi_p2 + SMALL));
    while(delta > 1e-5)
    {
        psi_p2 = psi_p;

        vector grad = gradLeastSquare(p, psi_p2, c, psi); 

        scalarList psi_bar(psi.size(), 0); 
        forAll(psi, i)
        {
            psi_bar[i] = psi[i] - (grad & (c[i] - p));   
        }                

        psi_p = inverseDistanceInterpolate(p, c, psi_bar, r);

        delta = fabs((psi_p - psi_p2) / (psi_p2 + SMALL));

        if ((++loop) > 20) break; 
    }

    scalar minValue = min(psi);
    scalar maxValue = max(psi);
    if (psi_p < minValue || psi_p > maxValue)
        return psi_p0;
    else
        return psi_p;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolationSchemes::interpolationSchemes
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolationSchemes::~interpolationSchemes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::interpolationSchemes::IDWinterp
(
    const vector& p,
    const List<vector>& c,
    const scalarList& psi,
    scalar r
)
{
    if (c.size() == 1)
        return psi[0];
    else
        return inverseDistanceInterpolate(p, c, psi, r);
}

Foam::scalar Foam::interpolationSchemes::IDeCinterp
(
    const vector& p,
    const List<vector>& c,
    const scalarList& psi,
    scalar r
)
{
    if (c.size() == 1)
        return psi[0];

    else if (
                (mesh_.nSolutionD() == 3 && c.size() <= 4) ||
                (mesh_.nSolutionD() == 2 && c.size() <= 2)
            )
        return inverseDistanceInterpolate(p, c, psi, r);

    else
        return interpolateSecondOrder(p, c, psi, r);
}

Foam::scalar Foam::interpolationSchemes::LSfitting
(
    const vector& p,
    const List<vector>& c,
    const scalarList& psi,
    scalar r
)
{
    if (c.size() == 1)
        return psi[0];

    else if (
                (mesh_.nSolutionD() == 3 && c.size() <= 4) ||
                (mesh_.nSolutionD() == 2 && c.size() <= 2)
            )
        return inverseDistanceInterpolate(p, c, psi, r);

    multiDimPolyFitter<scalar> polyFitter("polyDegree1", mesh_.geometricD());

    DynamicList<vector> centres;
    DynamicList<scalar> values;

    forAll (c, i)
    {
        centres.append(c[i] - p);
        values.append(psi[i]);
    }

    List<scalar> fitData = polyFitter.fitData
    (
        centres,
        values
    );

    return fitData[0];
}


// ************************************************************************* //
