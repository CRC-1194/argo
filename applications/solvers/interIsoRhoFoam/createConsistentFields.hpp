const dimensionedScalar nu1 = mixture.nuModel1().viscosityProperties().get<dimensionedScalar>("nu");
const dimensionedScalar nu2 = mixture.nuModel2().viscosityProperties().get<dimensionedScalar>("nu");
surfaceScalarField muf(mixture.muf());
surfaceScalarField rhof("rhof", fvc::interpolate(rho));
surfaceScalarField alphaface("alphaface", fvc::interpolate(alpha1));
cutFaceAdvect cutfaceInfo(mesh, alpha1);
dimensionedScalar areaDim("areaDim",dimensionSet(0,2,0,0,0,0,0),1.0);
