    volScalarField psi
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
    
    volVectorField fSigma
    (
        IOobject
        (
            "fSigma", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fSigma", dimForce / pow(dimLength,3), vector(0,0,0)) 
    );

    surfaceVectorField fSigmaFace 
    (
        IOobject
        (
            "fSigmaFace", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fSigmaFace", dimForce / pow(dimLength,3), vector(0,0,0)) 
    );

    volScalarField kappa 
    (
        IOobject
        (
            "kappa", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", pow(dimLength, -1), 0.)
    );

    // Young-Laplace equation
    
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh", 
            runTime.timeName(), 
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    