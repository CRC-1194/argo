fvScalarMatrix YoungLaplaceEqn
(
    fvm::laplacian(p_rgh) == 
        fvc::div(fSigma*mesh.magSf())
);
YoungLaplaceEqn.solve();
p_rgh.write();
