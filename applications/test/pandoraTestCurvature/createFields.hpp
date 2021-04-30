volScalarField alpha
(
    IOobject
    (
        fieldname,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField curvatureRequired
(
    IOobject
    (
        "cell_curvature_required",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar{"noname", dimless, 0.0}
);
volScalarField curvatureErrorField
(
    IOobject
    (
        "curvature_errors_interface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar{"noname", dimless, 0.0}
);
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "createPhi.H"
isoAdvection advector(alpha, phi, U);
