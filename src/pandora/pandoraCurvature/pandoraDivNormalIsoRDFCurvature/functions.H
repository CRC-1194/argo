
vector2D gradLeastSquare
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
    t.xx() = sum(w*deltax*deltax); 
    t.xy() = sum(w*deltax*deltay); 
    t.yx() = sum(w*deltay*deltax); 
    t.yy() = sum(w*deltay*deltay);

    vector2D b; 
    b.x() = sum(w*deltax*deltaAlpha); 
    b.y() = sum(w*deltay*deltaAlpha); 

    tensor2D invt = inv(t);

    return (invt & b); 
}

vector gradLeastSquare
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
    t.xx() = sum(w*deltax*deltax); 
    t.xy() = sum(w*deltax*deltay); 
    t.xz() = sum(w*deltax*deltaz); 

    t.yx() = sum(w*deltay*deltax); 
    t.yy() = sum(w*deltay*deltay);
    t.yz() = sum(w*deltay*deltaz);

    t.zx() = sum(w*deltaz*deltax); 
    t.zy() = sum(w*deltaz*deltay);
    t.zz() = sum(w*deltaz*deltaz);

    vector b; 
    b.x() = sum(w*deltax*deltaAlpha); 
    b.y() = sum(w*deltay*deltaAlpha); 
    b.z() = sum(w*deltaz*deltaAlpha); 

    return (inv(t) & b); 
}

void convert
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

void convert
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

void convert
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

void convert
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

Vector2D<label> getDimension
(
    const fvMesh & mesh 
)
{
    const Vector<label> & solution = mesh.solutionD(); 
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

