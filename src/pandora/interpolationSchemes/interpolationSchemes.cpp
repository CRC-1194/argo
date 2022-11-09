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
#include "tensor2D.H"

#include "leastSquareGrad.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::Vector2D<Foam::label> Foam::interpolationSchemes::getDimension()
{
    const Vector<label> & solution = mesh_.solutionD(); 
    label i = 0; 
    Vector2D<label> index; 
    forAll(solution, j)
    {
        if(solution[j] == -1) continue; 
        index[i] = j;
        i++;
    }
    return index;
}

void Foam::interpolationSchemes::convert
(
    vector2D & v2d, 
    const vector & v3d,  
    const Vector2D<label> & index
)
{
    forAll(index, i)
    {
        v2d[i] = v3d[index[i]];
    }
}

void Foam::interpolationSchemes::convert
(
    Field<vector2D> & v2d, 
    const vectorField & v3d,  
    const Vector2D<label> & index
)
{
    forAll(v3d, i)
    {
        convert(v2d[i], v3d[i], index);
    }
}

void Foam::interpolationSchemes::convert
(
    vector & v3d, 
    const vector2D & v2d,  
    const Vector2D<label> & index
)
{
    forAll(index, i)
    {
        v3d[index[i]] = v2d[i];
    }
}

void Foam::interpolationSchemes::convert
(
    vectorField & v3d,
    const Field<vector2D> & v2d, 
    const Vector2D<label> & index
)
{
    forAll(v3d, i)
    {
        convert(v3d[i], v2d[i], index);
    }
}

Foam::vector Foam::interpolationSchemes::gradLeastSquare
(
    const vector & p, 
    const scalar & alphap,
    const Field<vector> & c,
    const scalarField & alphac
)
{
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());

    DynamicField<vector> centres(100);
    DynamicField<scalar> values(100);

    forAll (c, i)
    {
        centres.append(c[i]);
        values.append(alphac[i]);
    }

    centres -= p;
    values -= alphap;

    return lsGrad.grad(centres, values);
}

/*
Foam::vector Foam::interpolationSchemes::gradLeastSquare
(
    const vector & p, 
    const scalar & alphap,
    const Field<vector> & c,
    const scalarField & alphac
)
{
    scalarField deltax(c.size()); 
    scalarField deltay(c.size()); 
    scalarField deltaz(c.size()); 
    scalarField deltaAlpha(c.size());
    scalarField w(c.size());

    forAll(c, i)
    {
        deltax[i] = c[i][0] - p[0]; 
        deltay[i] = c[i][1] - p[1];
        deltaz[i] = c[i][2] - p[2];
        deltaAlpha[i] = alphac[i] - alphap;
        w[i] = 1.0/max(mag(c[i] - p), SMALL);
    }

    tensor t; 
    t.xx() = sum(w*w*deltax*deltax); 
    t.xy() = sum(w*w*deltax*deltay); 
    t.xz() = sum(w*w*deltax*deltaz); 

    t.yx() = sum(w*w*deltay*deltax); 
    t.yy() = sum(w*w*deltay*deltay);
    t.yz() = sum(w*w*deltay*deltaz);

    t.zx() = sum(w*w*deltaz*deltax); 
    t.zy() = sum(w*w*deltaz*deltay);
    t.zz() = sum(w*w*deltaz*deltaz);

    vector b; 
    b.x() = sum(w*w*deltax*deltaAlpha); 
    b.y() = sum(w*w*deltay*deltaAlpha); 
    b.z() = sum(w*w*deltaz*deltaAlpha); 

    return (inv(t) & b); 
}
*/

Foam::vector2D Foam::interpolationSchemes::gradLeastSquare
(
    const vector2D & p, 
    const scalar & alphap,
    const Field<vector2D> & c,
    const scalarField & alphac
)
{
    scalarField deltax(c.size()); 
    scalarField deltay(c.size()); 
    scalarField deltaAlpha(c.size());
    scalarField w(c.size());

    forAll(c, i)
    {
        deltax[i] = c[i][0] - p[0]; 
        deltay[i] = c[i][1] - p[1];
        deltaAlpha[i] = alphac[i] - alphap;
        w[i] = 1.0/mag(c[i] - p);
    }

    tensor2D t; 
    t.xx() = sum(w*w*deltax*deltax); 
    t.xy() = sum(w*w*deltax*deltay); 
    t.yx() = sum(w*w*deltay*deltax); 
    t.yy() = sum(w*w*deltay*deltay);

    vector2D b; 
    b.x() = sum(w*w*deltax*deltaAlpha); 
    b.y() = sum(w*w*deltay*deltaAlpha); 

    tensor2D invt = inv(t);

    return (invt & b); 
}

Foam::scalar Foam::interpolationSchemes::inverseDistanceInterpolate
(
    const vector & p, 
    const vectorField & c,
    const scalarField & psi,
    const scalar & r 
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

Foam::scalar Foam::interpolationSchemes::inverseDistanceInterpolate
(
    const vector & p, 
    const vector & n,
    const vectorField & c,
    const scalarField & psi,
    const scalar & r 
)
{
    scalar w = 0; 
    scalar sum = 0;
    forAll(psi, i)
    {
        scalar a = Foam::pow(mag((c[i] - p) & n), r);
        scalar d = 1.0 / max(a, SMALL);
        sum += psi[i]*d; 
        w += d; 
    }

    return sum/w;
}

Foam::scalar Foam::interpolationSchemes::inverseDistanceInterpolate
(
    const vector2D & p, 
    const Field<vector2D> & c,
    const scalarField & psi,
    const scalar & r 
)
{
    scalar w = 0; 
    scalar sum = 0;
    forAll(psi, i)
    {
        scalar d = 1.0 / Foam::pow(mag(c[i] - p), r);
        sum += psi[i]*d;
        w += d; 
    }

    return sum/w;
}

Foam::scalar Foam::interpolationSchemes::interpolateSecondOrder 
(
    const vector & p, 
    const vectorField & c, 
    const scalarField & psi, 
    const scalar & r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, c, psi, r); 

    if (mag(psi_p) < 1e-8)
        return psi_p;

    scalar psi_p2 = psi_p + 1; 

    scalar delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    while(delta > 1e-5)
    {
        psi_p2 = psi_p;

        vector grad;
        if(c.size() == 1)
        {
            grad.x() = (psi[0] - psi_p2) / (c[0].x() - p.x());
            grad.y() = (psi[0] - psi_p2) / (c[0].y() - p.y());
            grad.z() = (psi[0] - psi_p2) / (c[0].z() - p.z());
        }
        else
        {
            grad = gradLeastSquare(p, psi_p2, c, psi); 
        }

        scalarField psi_bar(psi.size()); 
        forAll(psi, i)
        {
            psi_bar[i] = psi[i]  - (grad & (c[i] - p));   
        }                

        psi_p = inverseDistanceInterpolate(p, c, psi_bar, r);

        delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    }
    return psi_p;
}

Foam::scalar Foam::interpolationSchemes::interpolateSecondOrder 
(
    const vector & p, 
    const vector & n, 
    const vectorField & c, 
    const scalarField & psi, 
    const scalar & r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, n, c, psi, r); 

    if (mag(psi_p) < 1e-8)
        return psi_p;

    scalar psi_p2 = psi_p + 1; 

    scalar delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    while(delta > 1e-5)
    {
        psi_p2 = psi_p;

        vector grad;
        if(c.size() == 1)
        {
            grad.x() = (psi[0] - psi_p2) / (c[0].x() - p.x());
            grad.y() = (psi[0] - psi_p2) / (c[0].y() - p.y());
            grad.z() = (psi[0] - psi_p2) / (c[0].z() - p.z());
        }
        else
        {
            grad = gradLeastSquare(p, psi_p2, c, psi); 
        }

        scalarField psi_bar(psi.size()); 
        forAll(psi, i)
        {
            //vector pt = (c[i] - p);
            vector pt = (c[i] - p) & n * n;
            psi_bar[i] = psi[i] - (grad & pt);   
        }                

        psi_p = inverseDistanceInterpolate(p, n, c, psi_bar, r);

        delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    }
    return psi_p;
}

Foam::scalar Foam::interpolationSchemes::interpolateSecondOrder 
(
    const vector2D & p, 
    const Field<vector2D> & c, 
    const scalarField & psi, 
    const scalar & r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, c, psi, r); 

    if (mag(psi_p) < 1e-8)
        return psi_p;

    scalar psi_p2 = psi_p + 1; 

    scalar delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    while(delta > 1e-5)
    {
        psi_p2 = psi_p;

        vector2D grad;

        if(c.size() == 1)
        {
            grad.x() = (psi[0] - psi_p2) / (c[0].x() - p.x());
            grad.y() = (psi[0] - psi_p2) / (c[0].y() - p.y());
        }
        else
        {
            grad = gradLeastSquare(p, psi_p2, c, psi); 
        }

        scalarField psi_bar(psi.size()); 

        forAll(psi, i)
        {
            psi_bar[i] = psi[i]  - (grad & (c[i] - p));   
        }                

        psi_p = inverseDistanceInterpolate(p, c, psi_bar, r);

        delta = fabs((psi_p-psi_p2) / (psi_p2+SMALL));
    }

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

Foam::scalar Foam::interpolationSchemes::IDWinterpolate
(
    const vector& p,
    const vectorField& c,
    const scalarField& psi,
    const label& r
)
{
    if(mesh_.nSolutionD() == 2)
    {
        Vector2D<label> index = getDimension();
        Field<vector2D> c2d(c.size());
        convert(c2d, c, index);
        vector2D p2d;
        convert(p2d, p, index);
        return inverseDistanceInterpolate(p2d, c2d, psi, r);
    }
    else
    {
        return inverseDistanceInterpolate(p, c, psi, r);
    }
}

Foam::scalar Foam::interpolationSchemes::IDeCinterpolate
(
    const vector& p,
    const vectorField& c,
    const scalarField& psi,
    const label& r
)
{
/*
    if(mesh_.nSolutionD() == 2)
    {
        Vector2D<label> index = getDimension();
        Field<vector2D> c2d(c.size());
        convert(c2d, c, index);
        vector2D p2d;
        convert(p2d, p, index);
        return interpolateSecondOrder(p2d, c2d, psi, r);
    }
    else
*/
    {
        return interpolateSecondOrder(p, c, psi, r);
    }
}

Foam::scalar Foam::interpolationSchemes::IDeCinterpolate
(
    const vector& p,
    const vector& n,
    const vectorField& c,
    const scalarField& psi,
    const label& r
)
{
    if(mesh_.nSolutionD() == 2)
    {
        Vector2D<label> index = getDimension();
        Field<vector2D> c2d(c.size());
        convert(c2d, c, index);
        vector2D p2d;
        convert(p2d, p, index);
        return interpolateSecondOrder(p2d, c2d, psi, r);
    }
    else
    {
        return interpolateSecondOrder(p, n, c, psi, r);
    }
}

Foam::scalar Foam::interpolationSchemes::LSinterpolate
(
    const vector& p,
    const vectorField& c,
    const scalarField& psi
)
{
    multiDimPolyFitter<scalar> polyFitter("polyDegree1", mesh_.geometricD());

    DynamicField<vector> centres;
    DynamicField<scalar> values;

    forAll (c, i)
    {
        centres.append(c[i]);
        values.append(psi[i]);
    }

    centres -= p;

    List<scalar> fitData = polyFitter.fitData
    (
        centres,
        values
    );

    return fitData[0];
}


// ************************************************************************* //
