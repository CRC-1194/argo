
scalar absError = 0;
scalar l1 = 0;
scalar l2 = 0;
scalar lInf = 0;
label count = 0;
volVectorField exactNormals = averagedNormals_;
volVectorField calcNormals = averagedNormals_;
forAll(markers, i)
{
    if (markers[i] == -1) continue;
    //if (markers[i] != 0) continue;

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
    if ((absError/mag(n1)) > lInf)
        lInf = absError/mag(n1);
}
reduce(l1, sumOp<scalar>());
reduce(l2, sumOp<scalar>());
reduce(count, sumOp<label>());
reduce(lInf, maxOp<scalar>());

Info<<"count = "<<count<<nl;
if (count > 0)
{
    Info<<"l1 = "<<l1/count*100<<nl;
    Info<<"l2 = "<<sqrt(l2/count)*100<<nl;
    Info<<"lInf = "<<lInf*100<<nl;
}

scalar absErrorRdf = 0;
scalar l1Rdf = 0;
scalar l2Rdf = 0;
scalar lInfRdf = 0;
label countRdf = 0;
forAll (rdf, i)
{
    if (markers[i] == -1 || markers[i] == nPropagate_) continue;

    scalar exactRdf = sphereRadius - mag(sphereCentre - mesh().C()[i]);
    absErrorRdf = fabs(rdf[i] - exactRdf);
    l1Rdf += absErrorRdf / fabs(exactRdf);
    l2Rdf += sqr(absErrorRdf / exactRdf);
    countRdf++;
    scalar inf = absErrorRdf / fabs(exactRdf);
    if (inf > lInfRdf)
        lInfRdf = inf;
}

Info<<"countRdf = "<<countRdf<<nl;
if (countRdf > 0)
{
    Info<<"l1Rdf = "<<l1Rdf/count*100<<nl;
    Info<<"l2Rdf = "<<sqrt(l2Rdf/count)*100<<nl;
    Info<<"lInfRdf = "<<lInfRdf*100<<nl;
}

/*
exactNormals.rename("normExct");
if (cellCurvature_.time().writeTime())
    exactNormals.write();
calcNormals.rename("normCalc");
if (cellCurvature_.time().writeTime())
    calcNormals.write();
*/

markers.rename("cellMarker");
if (mesh().time().writeTime())
    markers.write();
