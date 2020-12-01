    volScalarField signedDistance
    (
        IOobject
        (
            fieldName, 
            runTime.timeName(), 
            mesh, 
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField sqrSearchDist
    (
        IOobject
        (
            "sqrSearchDist", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::average(pow(mesh.deltaCoeffs(), -2))
    );
