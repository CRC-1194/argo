{

OFstream outFile("postProcessing/error_curv0.dat");
outFile<<"# Time, Length, LCurv1, LCurv2, LCurvInf"<<nl;

scalar absErrorCurv = 0;
scalar l1Curv = 0;
scalar l2Curv = 0;
scalar lInfCurv = 0;
label countCurv = 0;
forAll (cellCurvature_, i)
{
    if (cellDistLevel_[i] == -1) continue;
    if (cellDistLevel_[i] ==  2) continue;
    if (cellDistLevel_[i] ==  3) continue;

    //3D
    scalar exactCurv = 2.0 / mag(sphereCentre - mesh.C()[i]);
    //2D
    //scalar exactCurv = 1.0 / mag(sphereCentre - mesh.C()[i]);
    absErrorCurv = fabs(cellCurvature_[i] - exactCurv);
    l1Curv += absErrorCurv / exactCurv;
    l2Curv += sqr(absErrorCurv / exactCurv);
    countCurv++;
    scalar inf = absErrorCurv;
    if (inf > lInfCurv)
    {
        lInfCurv = inf;
    }
}
reduce(l1Curv, sumOp<scalar>());
reduce(l2Curv, sumOp<scalar>());
reduce(countCurv, sumOp<label>());
reduce(lInfCurv, maxOp<scalar>());

if (countCurv > 0)
{
    Info<<"curv0"<<nl;
    Info<<"l1Curv = "<<l1Curv/countCurv<<nl;
    Info<<"l2Curv = "<<sqrt(l2Curv/countCurv)<<nl;
    Info<<"lICurv = "<<lInfCurv<<nl;
}

outFile<<mesh.time().value()
       <<" "<<average(mag(mesh.delta())()).value()
       <<" "<<l1Curv/countCurv
       <<" "<<sqrt(l2Curv/countCurv)
       <<" "<<lInfCurv
       <<endl;

}
