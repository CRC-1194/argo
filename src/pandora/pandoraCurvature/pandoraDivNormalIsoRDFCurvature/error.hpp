
vector sphereCentre(0.005, 0.005, 0.005); // Sphere centre
scalar radius = 0.002; // Sphere radius
scalar absError = 0;
scalar l1 = 0;
scalar l2 = 0;
scalar lInf = 0;
label count = 0;
volVectorField exactNormals = averagedNormals_;
volVectorField calcNormals = averagedNormals_;
forAll(averagedNormals_, i)
{
    if (mag(averagedNormals_[i]) != 0)
    {
        if (markers[i] == 0) continue;

        vector n1 = sphereCentre - mesh().C()[i];
        n1 /= mag(n1);
        exactNormals[i] = n1;

        vector n2 = averagedNormals_[i];
        n2 /= mag(n2);
        calcNormals[i] = n2;

        absError = 1 - (n1 & n2);

        l1 += absError / mag(n1);
        l2 += sqr(absError) / magSqr(n1);
        count++;
        if (absError > lInf)
            lInf = absError;
    }
}
reduce(l1, sumOp<scalar>());
reduce(l2, sumOp<scalar>());
reduce(count, sumOp<label>());
reduce(lInf, maxOp<scalar>());

Info<<"l1 = "<<l1/count*100<<nl;
Info<<"l2 = "<<sqrt(l2/count)*100<<nl;
Info<<"lInf = "<<lInf*100<<nl;

exactNormals.rename("normExct");
if (cellCurvature_.time().writeTime())
    exactNormals.write();
calcNormals.rename("normCalc");
if (cellCurvature_.time().writeTime())
    calcNormals.write();

markers.rename("cellMarker");
if (mesh().time().writeTime())
    markers.write();
