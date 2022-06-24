    // - iso-surface field: can be signed-distance, or another field
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

    // - signed-distance field 
    volScalarField sigDist 
    (
        IOobject
        (
            "sigDist", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sigDist", dimLength, 0.)
    );
    
    // - Squared octree search radius
    volScalarField sqrSearchRadius 
    (
        IOobject
        (
            "sqrSearchRadius", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sqrSearchRadius", Foam::sqr(dimLength), 0)
    );
