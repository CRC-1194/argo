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

    // Initial signed distance field given by the octree.
    volScalarField signedDist0("signedDist0", signedDist); 

    // The diffusion coefficient field used for sign propagation. 
    surfaceScalarField lambda 
    (
        IOobject
        (
            "lambda", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        // Solve a Laplace equation to propagate the sign.
        dimensionedScalar ("lambda", sqr(dimLength) * pow(dimTime,-1), 1)
    );
