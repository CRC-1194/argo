    // Volume fraction field
    volScalarField alpha
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

    // Squared search distance field for sign propagation.
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
        mesh,
        dimensionedScalar("sqrSearchDist", Foam::pow(dimLength,2),0)
    );

    // Squared search distance for intersection. 
    volScalarField sqrInterDist
    (
        IOobject
        (
            "sqrInterDist", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("sqrInterDist", Foam::pow(dimLength,2),0)
    );


    // The signed distance field 
    volScalarField signedDist 
    (
        IOobject
        (
            "signedDist", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("signedDist", dimLength,0),
        "zeroGradient"
    );
