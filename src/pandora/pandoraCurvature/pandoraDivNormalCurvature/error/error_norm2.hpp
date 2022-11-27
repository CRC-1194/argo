{

//vector sphereCentre(1.0, 0.0, 1.0);

OFstream outFile("postProcessing/reconstructionError/error_norm2.dat");
outFile<<"# Time, Length, LNormal1, LNormal2, LNormalInf"<<nl;

scalar absError = 0;
scalar l1 = 0;
scalar l2 = 0;
scalar lInf = 0;
label count = 0;
volVectorField exactNormals = averagedNormals_;
volVectorField calcNormals = averagedNormals_;
forAll(averagedNormals_, cellI)
{
    if (cellDistLevel_[cellI] != 2) continue;

    vector n1 = sphereCentre - mesh.C()[cellI];
    n1 /= mag(n1);
    exactNormals[cellI] = n1;

    vector n2 = averagedNormals_[cellI];
    n2 /= mag(n2);
    calcNormals[cellI] = n2;

//averagedNormals_[cellI] = n1;

    absError = 1 - (n1 & n2);
    l1 += absError / mag(n1);
    l2 += sqr(absError) / magSqr(n1);
    count++;
    if ((absError/mag(n1)) > lInf)
    {
        lInf = absError/mag(n1);
    }
}
reduce(l1, sumOp<scalar>());
reduce(l2, sumOp<scalar>());
reduce(count, sumOp<label>());
reduce(lInf, maxOp<scalar>());

//Info<<"count = "<<count<<nl;
if (count > 0)
{
    Info<<"norm2"<<nl;
    Info<<"l1 = "<<l1/count<<nl;
    Info<<"l2 = "<<sqrt(l2/count)<<nl;
    Info<<"lInf = "<<lInf<<nl;
}

outFile<<mesh.time().value()
       <<" "<<average(mag(mesh.delta())()).value()
       <<" "<<l1/count
       <<" "<<sqrt(l2/count)
       <<" "<<lInf
       <<endl;

}
