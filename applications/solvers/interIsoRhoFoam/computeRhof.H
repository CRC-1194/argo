// rhof == fvc::interpolate(rho);
rhof == upwind<scalar>(phi.mesh(), phi).interpolate(rho);

// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
const auto& faceNeighbour = mesh.neighbour();

alphaf == upwind<scalar>(phi.mesh(), phi).interpolate(alpha1);

//Jun: reference to reconstructionSchemes

// const reconstructionSchemes& surf = advector.surf();-->
reconstructionSchemes& surf = advector.surf(); 
surf.reconstruct(); 

const boolList& isASurfaceCell = surf.interfaceCell();
const surfaceScalarField& magSf = mesh.magSf();

forAll(alphaf, faceI) // alphaf geometric calculation
{ 
    const scalar& sfArea = magSf[faceI];
    
    if (phi[faceI] > 0)
    {
        if(isASurfaceCell[faceOwner[faceI]])
        {
            // cell interface centre
            const point x0 = surf.centre()[faceOwner[faceI]]; 
            // cell interface unit normal
            vector n0 = surf.normal()[faceOwner[faceI]]; 
            n0 /= (mag(n0)+SMALL);
            cutFace.calcSubFace(faceI, n0, x0);
            alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea; 
        } 
    } 
    else 
    {
        if(isASurfaceCell[faceNeighbour[faceI]])
        {
            // cell interface centre
            const point x0 = surf.centre()[faceNeighbour[faceI]]; 
            // cell interface unit normal
            vector n0 = surf.normal()[faceNeighbour[faceI]]; 
            n0 /= (mag(n0)+SMALL);
            cutFace.calcSubFace(faceI, n0, x0);
            alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea; 
        }        
    }
}

rhof == alphaf*rho1 + (1 - alphaf)*rho2;
muf == alphaf*rho1*nu1 + (1 - alphaf)*rho2*nu2;

// TODO: alphaf geometric upwinding in parallel. 
//const surfaceScalarField::Boundary& phib = phi.boundaryField();
//// const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
//surfaceScalarField::Boundary& alphafb = alphaf.boundaryFieldRef();
//const surfaceScalarField::Boundary& magSfb = magSf.boundaryField();
//surfaceScalarField::Boundary& rhofb = rhof.boundaryFieldRef();
//surfaceScalarField::Boundary& mufb = muf.boundaryFieldRef();

//forAll(mesh.boundaryMesh(), patchID) 
//{
    //const UList<label>& boundaryCells = mesh.boundaryMesh()[patchID].faceCells();
    //forAll (mesh.boundaryMesh()[patchID],faceI)
    //{
        //// const label& faceI = mesh.boundaryMesh()[patchID].start() + face;

         //alphafb[patchID][faceI] =  alpha1[boundaryCells[faceI]]; //alphafUpwind.boundaryField()[patchID][faceI];

        //if(isASurfaceCell[boundaryCells[faceI]] && (mag(phib[patchID][faceI]) > 10*SMALL))
        //{
           //const point x0 = surf.centre()[boundaryCells[faceI]]; // cell interface centre
           //vector n0 = surf.normal()[boundaryCells[faceI]]; // cell interface unit normal
           //n0 /= (mag(n0)+SMALL);
           //cutFace.calcSubFace(faceI, n0, x0);
           //alphafb[patchID][faceI] = mag(cutFace.subFaceArea()) / magSfb[patchID][faceI]; 
        //}

        //// Info << "boundary patch name: " << mesh.boundaryMesh()[patchID].name() <<"; boundary cell alpha: "<< alpha1[boundaryCells[faceI]] << "; " << endl;
        //// Info << "is a surface cell? : " << isASurfaceCell[boundaryCells[faceI]] <<"; alphaf: " << alphafb[patchID][faceI] << endl;

        //rhofb[patchID][faceI] = alphafb[patchID][faceI]*rho1.value() + (1 - alphafb[patchID][faceI])*rho2.value();
        //mufb[patchID][faceI] = alphafb[patchID][faceI]*rho1.value()*nu1.value() +(1 - alphafb[patchID][faceI])*rho2.value()*nu2.value();     
    //}
//}
