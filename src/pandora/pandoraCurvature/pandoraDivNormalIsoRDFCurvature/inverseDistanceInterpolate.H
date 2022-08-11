
scalar inverseDistanceInterpolate
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

scalar interpolateSecondOrder 
(
    const vector & p, 
    const vectorField & c, 
    const scalarField & psi, 
    const scalar & r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, c, psi, r); 
    scalar psi_p2 = psi_p + 1; 

    while(fabs(psi_p - psi_p2) > 1e-5)
    //for(int j=0; j<iter; j++)
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
    }
    return psi_p;
}


scalar inverseDistanceInterpolate
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

scalar interpolateSecondOrder 
(
    const vector2D & p, 
    const Field<vector2D> & c, 
    const scalarField & psi, 
    const scalar & r
)
{
    scalar psi_p = inverseDistanceInterpolate(p, c, psi, r); 
    scalar psi_p2 = psi_p + 1; 
    scalarField pt(psi.size());

    while(fabs(psi_p - psi_p2) > 1e-5)
    //for(int j=0; j<iter; j++)
    {
        psi_p2 = psi_p;

        vector2D grad;

        if(c.size() == 1)
        {
            grad.x() = (psi[0] - psi_p2) / (c[0].x() - p.x());
            grad.y() = (psi[0] - psi_p2) / (c[0].y() - p.y());
        }

        else
            grad = gradLeastSquare(p, psi_p2, c, psi); 

        scalarField psi_bar(psi.size()); 

        forAll(psi, i)
        {
            psi_bar[i] = psi[i]  - (grad & (c[i] - p));   
        }                
        pt = psi_bar;
        psi_p = inverseDistanceInterpolate(p, c, psi_bar, r);
    }

    return psi_p;
}

