//vector sphereCentre(1.0, 0.0, 1.0);

{

OFstream outFile("postProcessing/error_rdf0.dat");
outFile<<"# Time, Length, LRDF1, LRDF2, LRDFInf"<<nl;

//scalar sphereRadius = 0.5;
scalar absErrorRdf = 0;
scalar l1Rdf = 0;
scalar l2Rdf = 0;
scalar lInfRdf = 0;
label countRdf = 0;
forAll (RDF, i)
{
    if (!nextToInter_[i]) continue;
    //if (markers[i] != 0) continue;
    if (cellDistLevel_[i] != 0) continue;

    scalar exactRdf = sphereRadius - mag(sphereCentre - mesh.C()[i]);
    absErrorRdf = fabs(RDF[i] - exactRdf);
    l1Rdf += absErrorRdf / sphereRadius;
    l2Rdf += sqr(absErrorRdf / sphereRadius);
    countRdf++;
    scalar inf = absErrorRdf / sphereRadius;
    if (inf > lInfRdf)
    {
        lInfRdf = inf;
    }
}
reduce(l1Rdf, sumOp<scalar>());
reduce(l2Rdf, sumOp<scalar>());
reduce(countRdf, sumOp<label>());
reduce(lInfRdf, maxOp<scalar>());

//Info<<"countRdf = "<<countRdf<<nl;
if (countRdf > 0)
{
    Info<<"rdf0"<<nl;
    Info<<"l1Rdf = "<<l1Rdf/countRdf<<nl;
    Info<<"l2Rdf = "<<sqrt(l2Rdf/countRdf)<<nl;
    Info<<"lInfRdf = "<<lInfRdf<<nl;
}

outFile<<mesh.time().value()
       <<" "<<average(mag(mesh.delta())()).value()
       <<" "<<l1Rdf/countRdf
       <<" "<<sqrt(l2Rdf/countRdf)
       <<" "<<lInfRdf
       <<endl;
}
