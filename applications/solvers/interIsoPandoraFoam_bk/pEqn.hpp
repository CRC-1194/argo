{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
    );
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

scalar sumPhiHbyA = 0;
forAll (phiHbyA, i)
    if (phiHbyA[i] > 0)
        sumPhiHbyA += phiHbyA[i];
reduce(sumPhiHbyA, sumOp<scalar>());
Info<<"sumPhiHbyA = "<<sumPhiHbyA<<nl;

scalar sumFSigma = 0;
surfaceScalarField stf = mixture.surfaceTensionForce();
forAll (stf, i)
    if (stf[i] > 0)
        sumFSigma += stf[i];
//forAll (fSigma, i)
//    if (fSigma[i] > 0)
//        sumFSigma += fSigma[i];
reduce(sumFSigma, sumOp<scalar>());
Info<<"sumFSigma = "<<sumFSigma<<nl;

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
            //pandoraModel.surfaceTensionForce(alpha1)
            //fSigma
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

scalar sumPhig = 0;
forAll (phig, i)
    if (phig[i] > 0)
        sumPhig += phig[i];
reduce(sumPhig, sumOp<scalar>());
Info<<"sumPhig = "<<sumPhig<<nl;

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            //U = fvc::reconstruct(phi);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}