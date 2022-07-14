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
    
    IOdictionary setAlphaFields
    (
        IOobject
        (
            "setAlphaFieldDict",
            "system",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::NO_WRITE
        )
    );
