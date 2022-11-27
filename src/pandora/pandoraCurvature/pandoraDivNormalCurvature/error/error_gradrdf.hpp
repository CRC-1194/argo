{

//vector sphereCentre(1.0, 0.0, 1.0);

OFstream outFile("postProcessing/error_gradrdf.dat");
outFile<<"# Time, Length, LNormal1, LNormal2, LNormalInf"<<nl;

scalar absError = 0;
scalar l1 = 0;
scalar l2 = 0;
scalar lInf = 0;
label count = 0;
volVectorField exactNormals = gradRDF;
volVectorField calcNormals = gradRDF;
forAll(gradRDF, i)
{
    //if (!nextToInter[i]) continue;
    //if (markers[i] != 0) continue;
    if (cellDistLevel_[i] != 0) continue;

    vector n1 = sphereCentre - mesh.C()[i];
    n1 /= mag(n1);
    exactNormals[i] = n1;

    vector n2 = gradRDF[i];
    n2 /= mag(n2);
    calcNormals[i] = n2;

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
    Info<<"gradrdf"<<nl;
    Info<<"l1gradRDF = "<<l1/count<<nl;
    Info<<"l2gradRDF = "<<sqrt(l2/count)<<nl;
    Info<<"lIgradRDF = "<<lInf<<nl;
}

outFile<<mesh.time().value()
       <<" "<<average(mag(mesh.delta())()).value()
       <<" "<<l1/count
       <<" "<<sqrt(l2/count)
       <<" "<<lInf
       <<endl;

}
