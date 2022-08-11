
    scalar absErrorCurv = 0;
    scalar l1Curv = 0;
    scalar l2Curv = 0;
    scalar lInfCurv = 0;
    label countCurv = 0;
    forAll(markers, cellI)
    {
        if (markers[cellI] == -1) continue;
        //if (markers[cellI] == 0) continue;

//Info<<cellCurvature_[cellI]<<nl;
        absErrorCurv = fabs(1000 - cellCurvature_[cellI]);
        l1Curv += absErrorCurv / 1000;
        l2Curv += Foam::sqr(absErrorCurv/1000);
        countCurv++;
        if (absErrorCurv > lInfCurv)
        {
            lInfCurv = absErrorCurv;
        }
    }

    Info<<"countCurv = "<<countCurv<<nl;
    if (countCurv > 0)
    {
        Info<<"l1Curv = "<<l1Curv/countCurv<<nl;
        Info<<"l2Curv = "<<Foam::sqrt(l2Curv/countCurv)<<nl;
        Info<<"lInfCurv = "<<lInfCurv<<nl;
    }

