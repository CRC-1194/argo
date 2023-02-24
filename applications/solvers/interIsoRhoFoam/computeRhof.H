// rhof == fvc::interpolate(rho);
// rhof == upwind<scalar>(phi.mesh(), phi).interpolate(rho);
{
// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
// const auto& faceNeighbour = mesh.neighbour();
const cellList& cellFaces = mesh.cells();

alphaf == upwind<scalar>(phi.mesh(), phi).interpolate(alpha1);

//Jun: reference to reconstructionSchemes

// const reconstructionSchemes& surf = advector.surf();-->
// stop reconstructing and see the different alphaPhi/phi and alphaf

// reconstructionSchemes& surf = advector.surf(); 
//alpha1_cr == 0.5*(alpha1 + alpha1.oldTime());
alpha1_cr == alpha1;
advector_cr.surf().reconstruct();

/*----------compute Rhof---------------*/
const boolList& isASurfaceCell = advector_cr.surf().interfaceCell();

//-----> compute Rhof in internal field

const scalarField& phiIn = phi.primitiveField();
const scalarField& magSfIn = mesh.magSf().primitiveField();
DynamicList<label> bsFaces(label(0.2*mesh.nBoundaryFaces()));
DynamicList<vector> bsx0(bsFaces.size());
DynamicList<vector> bsn0(bsFaces.size());

const DynamicField<label>& interfaceLabels = advector_cr.surf().interfaceLabels(); 

forAll(interfaceLabels, i)
{
    const label celli = interfaceLabels[i];
    if (mag(advector_cr.surf().normal()[celli]) != 0)
    {
        const point x0 = advector_cr.surf().centre()[celli];
        vector n0 = advector_cr.surf().normal()[celli];
        n0 /= (mag(n0));
        const cell& celliFaces = cellFaces[celli];
        forAll(celliFaces, fi)
        {
            const label facei = celliFaces[fi];

            if (mesh.isInternalFace(facei))
            {
                bool isDownwindFace = false;

                if (celli == faceOwner[facei])
                {
                    if (phiIn[facei] >= 0)
                    { isDownwindFace = true; }
                }
                else
                {
                    if (phiIn[facei] < 0)
                    { isDownwindFace = true; }
                }

                if (isDownwindFace)
                {
                    cutFace.calcSubFace(facei, n0, x0);
                    alphaf[facei] = mag(cutFace.subFaceArea()) / magSfIn[facei];
                }
            }
            else
            {
                bsFaces.append(facei);
                bsx0.append(x0);
                bsn0.append(n0);
            }
        }
    }
}

//-----> compute Rhof in boundary field
const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
const surfaceScalarField::Boundary& phib = phi.boundaryField();
const surfaceScalarField::Boundary& magSfb = mesh.magSf().boundaryField();
surfaceScalarField::Boundary& alphafb = alphaf.boundaryFieldRef();
const label nInternalFaces = mesh.nInternalFaces();
DynamicList<label> procPatchLabels(mesh.boundary().size());
forAll(boundaryMesh, patchi)
{
   if
   (
       isA<processorPolyPatch>(boundaryMesh[patchi])
    && boundaryMesh[patchi].size() > 0
   )
   {
       procPatchLabels.append(patchi);
   }
}
List<DynamicList<label>> surfaceCellFacesOnProcPatches(0);
surfaceCellFacesOnProcPatches.resize(boundaryMesh.size());

forAll(bsFaces, i)
{
    const label facei = bsFaces[i];
    const label patchi = boundaryMesh.patchID()[facei - nInternalFaces];
    const polyPatch& pp = boundaryMesh[patchi];
    const label start = boundaryMesh[patchi].start();
    const label patchFacei = facei - start;

    if (phib[patchi].size())
    {
        const scalar phiP = phib[patchi][patchFacei];

        if(phiP >= 0)
        {
            const scalar magSf = magSfb[patchi][patchFacei];

            cutFace.calcSubFace(facei, bsn0[i], bsx0[i]);
            alphafb[patchi][patchFacei] = mag(cutFace.subFaceArea()) / magSf;

            // Handling upwind cyclic boundary patches 
            const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);
            if (cpp)
            {
                label neiPatchID = cpp->neighbPatchID();
                alphafb[neiPatchID][patchFacei] = alphafb[patchi][patchFacei];
            }
                
            // Check if the face is on processor patch and append it to
            // the list if necessary
            if(isA<processorPolyPatch>(pp) && pp.size())
            {
                surfaceCellFacesOnProcPatches[patchi].append(patchFacei);
            }
        }
    }
}
//#include"processorBC.H"
syncProcPatches(mesh, procPatchLabels, surfaceCellFacesOnProcPatches, alphaf, phi);
bsFaces.clear();
bsx0.clear();
bsn0.clear();
}
/*
forAll(alphaf, faceI) // alphaf geometric calculation
{ 
    const scalar& sfArea = magSf[faceI];
    
    if (phi[faceI] > 0)
    {
        if(isASurfaceCell[faceOwner[faceI]])
        {
            // cell interface centre
            const point x0 = advector_cr.surf().centre()[faceOwner[faceI]]; 
            // cell interface unit normal
            vector n0 = advector_cr.surf().normal()[faceOwner[faceI]];
	    if (mag(n0) != 0)
	    { 
                n0 /= mag(n0);
                cutFace.calcSubFace(faceI, n0, x0);
                alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea;
	    }*/
/*            const DynamicList<point>& subFacePoints = cutFace.subFacePoints();
            const label nPoints = subFacePoints.size();
            if (nPoints > 0)
            {
                isoFacePts.append(subFacePoints);
            } */
/*        } 
    } 
    else 
    {
        if(isASurfaceCell[faceNeighbour[faceI]])
        {
            // cell interface centre
            const point x0 = advector_cr.surf().centre()[faceNeighbour[faceI]]; 
            // cell interface unit normal
            vector n0 = advector_cr.surf().normal()[faceNeighbour[faceI]];
	    if (mag(n0) != 0)
	    { 
                n0 /= mag(n0);
                cutFace.calcSubFace(faceI, n0, x0);
                alphaf[faceI] = mag(cutFace.subFaceArea()) / sfArea;
             } 

/*            const label nPoints = subFacePoints.size();
            if (nPoints > 0)
            {
                isoFacePts.append(subFacePoints);
            }*/
 
/*        }        
    }
}*/

// writeIsoFaces(mesh_,isoFacePts,writeIsoFacesToFile);

Info<< "### display alphaf: "
                << "  Min(" << alphaf.name() << ") = " << min(alphaf).value()
                << "  Max(" << alphaf.name() << ") = " << max(alphaf).value()
                << endl;

rhof == alphaf*rho1 + (1 - alphaf)*rho2;
//muf == alphaf*rho1*nu1 + (1 - alphaf)*rho2*nu2;
rhoPhi == rhof*phi;


/*----------------write isofaces--------*/
/* JUN: comment out this part because we can use functionObject surface 
   to do the same work.
DynamicList<List<point>> isoFacePts;
const fvMesh& mesh_ = alpha1.mesh();
const dictionary dict_= mesh_.solverDict(alpha1.name());
// bool writeIsoFacesToFile = dict_.getOrDefault("writeIsoFaces", false); 
const DynamicField<label>& interfaceLabels = surf.interfaceLabels();
forAll(interfaceLabels, i)
{
    const label celli = interfaceLabels[i];
    vector normali = -surf.normal()[celli];
    if (mag(normali) == 0)
        continue;
    }
    sIterPLIC.vofCutCell(celli, alpha1[celli], 1e-8, 100, normali);
    isoFacePts.append(sIterPLIC.facePoints());
}
writeIsoFaces(mesh_,isoFacePts,writeIsoFacesToFile);
*/


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
