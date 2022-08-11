
scalar interpolateFunction
(
    const fvMesh & mesh,
    const vector & p,
    const vectorField & c,
    const scalarField & psi,
    const scalar & r
)
{
    if(mesh.nSolutionD() == 2)
    {
        Vector2D<label> index = getDimension(mesh);
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

scalar secondOrderInterpolate
(
    const fvMesh & mesh,
    const vector & p,
    const vectorField & c,
    const scalarField & psi,
    const scalar & r
)
{
    if(mesh.nSolutionD() == 2)
    {
        Vector2D<label> index = getDimension(mesh);
        Field<vector2D> c2d(c.size());
        convert(c2d, c, index);
        vector2D p2d;
        convert(p2d, p, index);
        return interpolateSecondOrder(p2d, c2d, psi, r);
    }
    else
    {
        return interpolateSecondOrder(p, c, psi, r);
    }
}
