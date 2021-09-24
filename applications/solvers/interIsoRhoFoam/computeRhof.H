// FIXME: This shouldn't be boilerplate, refactor. TM. 
// Compute rhoPhi using the old approach since it is valid for all bulk cells
rhof.oldTime() == rhof;
rhof == fvc::interpolate(rho);
// --> alphaf == fvc::interpolate(alpha1);

// #include "SortableList.H"
// #include <assert.h>

// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
const auto& faceNeighbour = mesh.neighbour();
// //const auto& cellList = mesh.cells();
// //const auto& phiIn = phi.primitiveField();

const surfaceScalarField& alphafUpwind = upwind<scalar>(phi.mesh(), phi).interpolate(alpha1);

// // const auto& meshVertices = mesh.points();
// // const auto& meshFaces = mesh.faces();

 const reconstructionSchemes& surf = advector.surf(); //Jun: reference to reconstructionSchemes
// // const DynamicField<label>& interfaceLabels = surf.interfaceLabels(); // Jun: list of interface cells
 const boolList& isASurfaceCell = surf.interfaceCell();

forAll(rhoPhi, faceI) 
{ 
    const scalar& sfArea = mesh.magSf()[faceI];
    alphaf[faceI] = alphafUpwind[faceI];
    
    if (phi[faceI] > 10*SMALL)
    {
        if(isASurfaceCell[faceOwner[faceI]])
        {
           const point x0 = surf.centre()[faceOwner[faceI]]; // cell interface centre
           vector n0 = -surf.normal()[faceOwner[faceI]]; // cell interface unit normal
           n0 /= (mag(n0)+SMALL);
           cutFace.calcSubFace(faceI, n0, x0);
           alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea; 
        } 
    } 
    else if (phi[faceI] < -10*SMALL)
    {
        if(isASurfaceCell[faceNeighbour[faceI]])
        {
           const point x0 = surf.centre()[faceNeighbour[faceI]]; // cell interface centre
           vector n0 = -surf.normal()[faceNeighbour[faceI]]; // cell interface unit normal
           n0 /= (mag(n0)+SMALL);
           cutFace.calcSubFace(faceI, n0, x0);
           alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea; 
        }        
    }
    else
    {
        if (isASurfaceCell[faceOwner[faceI]] && !(isASurfaceCell[faceNeighbour[faceI]]))
        {
            alphaf[faceI] = alpha1[faceNeighbour[faceI]];
        }
        else if(isASurfaceCell[faceNeighbour[faceI]] && !(isASurfaceCell[faceOwner[faceI]]))
        {
            alphaf[faceI] = alpha1[faceOwner[faceI]];
        }
    } 

    rhof[faceI] = alphaf[faceI]*rho1.value() + (1 - alphaf[faceI])*rho2.value();
    muf[faceI] = alphaf[faceI]*rho1.value()*nu1.value() +(1 - alphaf[faceI])*rho2.value()*nu2.value(); 
}    
