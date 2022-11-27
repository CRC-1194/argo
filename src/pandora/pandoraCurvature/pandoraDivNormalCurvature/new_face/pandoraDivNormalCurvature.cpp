/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pandoraDivNormalCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSmooth.H"
#include "surfaceInterpolate.H"
#include "error.H"

#include "processorFvPatchField.H"
#include "reconstructionSchemes.H"
#include "reconstructedDistanceFunction.H"
#include "tensor2D.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"
#include "triSurfaceMesh.H"

#include "wedgePolyPatch.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraDivNormalCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraDivNormalCurvature, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraDivNormalCurvature::pandoraDivNormalCurvature
(
    const fvMesh& mesh, 
    const dictionary& dict
)
    :
        pandoraCurvature(mesh, dict), 
        fieldName_(curvatureDict_.get<word>("normalField")),
        nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 6)), 
        nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 0)), 
        averagedNormals_ 
        (
            IOobject
            (
                "averagedNormals", 
                mesh.time().timeName(), 
                mesh,
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ),
            mesh, 
            dimensionedVector("averagedNormals", dimless, vector(0,0,0))
        ),
        curvFromTr_(dict.lookupOrDefault("curvFromTr",true)),
        markers_
        (
            IOobject
            (
                "cellMarkers",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("cellMarkers", dimless, -1),
            "calculated"
        ),
        interCells_(mesh.nCells(), false),
        cellDistLevel_
        (
            IOobject
            (
                "cellDistLevel_",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("cellDistLevel_", dimless, -1),
            "calculated"
        ),
        nextToInter_(mesh.nCells(), false),
        counts_(mesh.nCells(), 0),
        distribute_(zoneDistribute::New(mesh)),
        interp_(mesh),
        index_(mesh.nSolutionD(), -1),
        maxCount_(nPropagate_, 0.0)
{
    Info<<"Selecting divNormal curvature"<<nl;

    // curvature from trace does not work with wedges
    bool wedge = false;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    for (const polyPatch& pp : patches)
    {
        if(isA<wedgePolyPatch>(pp))
        {
            wedge = true;
        }
    }

    if(wedge)
    {
        curvFromTr_ = false;
    }

    // Get the dimension of the case. 
    const labelVector& solutionD = mesh.solutionD();
    label k = 0;
    forAll (solutionD, j)
    {
        if (solutionD[j] == 1)
        {
            index_[k++] = j;
        }
    }
}

// * * * * * * * * * * Protect Member Functions  * * * * * * * * * * * * * * //

void pandoraDivNormalCurvature::normalise(vectorField& vec)
{
    forAll (vec, i)
    {
        if (mag(vec[i]) < ROOTVSMALL)
            vec[i] = vector::zero;
        else
            vec[i] /= mag(vec[i]);
    }
}

void Foam::pandoraDivNormalCurvature::updateMarkersAndCounts()
{
    const fvMesh& mesh = cellCurvature_.mesh();

    forAll (markers_, cellI)
    {
        if (cellDistLevel_[cellI] == 0)
        {
            markers_[cellI] = 0;
        }
        else
        {
            markers_[cellI] = -1;
        }
    }
    markers_.correctBoundaryConditions();

    const auto& own = mesh.owner();
    const auto& nei = mesh.neighbour();

    for (label i = 0; i < nPropagate_; ++i)
    {
        labelField counts(mesh.nCells(), 0);
        volScalarField markersTmp = markers_;

        forAll (nei, fid)
        {
            if (
                    markers_[own[fid]] == i 
                 && markers_[nei[fid]] == -1 
                 && nextToInter_[nei[fid]]
               )
            {
                markersTmp[nei[fid]] = i + 1;
                counts[nei[fid]]++;
            }
            if (
                    markers_[nei[fid]] == i 
                 && markers_[own[fid]] == -1 
                 && nextToInter_[own[fid]]
               )
            {
                markersTmp[own[fid]] = i + 1;
                counts[own[fid]]++;
            }
        }
        markersTmp.correctBoundaryConditions();

        auto& meshBoundary = markers_.boundaryFieldRef();
        for (auto& mb : meshBoundary)
        {
            if (isA<processorFvPatch>(mb.patch()))
            {
                mb.initEvaluate();
                mb.evaluate();

                const auto& faceToCell = mb.patch().faceCells();
                const auto& neibrValue = mb.patchNeighbourField().cref();
                
                forAll(mb, j)
                {
                    if (
                            neibrValue[j] == i 
                         && markers_[faceToCell[j]] == -1 
                         && nextToInter_[faceToCell[j]]
                       )
                    {
                        markersTmp[faceToCell[j]] = i + 1;
                        counts[faceToCell[j]]++;
                    }
                }
            }
        }
        markersTmp.correctBoundaryConditions();
        markers_ = markersTmp;
        markers_.correctBoundaryConditions();

        maxCount_[i] = max(counts);
        reduce(maxCount_[i], maxOp<label>());
        counts_ += counts;
    }
}

void Foam::pandoraDivNormalCurvature::normalPropagate
(
    const bool& needUpdate, 
    volVectorField& interCellNormals
)
{
    const fvMesh& mesh = cellCurvature_.mesh();

    for (label i = 0; i < nPropagate_; ++i)
    {
        boolList zone(mesh.nCells(), false);
        forAll (zone, zi)
        {
            if (markers_[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }
        distribute_.setUpCommforZone(zone, false);
        Map<vector> mapMC = 
            distribute_.getDatafromOtherProc(zone, mesh.C());
        const labelListList& stencil = distribute_.getStencil();

        // Perform the interpolation in order from largest to smallest
        // neighbour cell numbers. Guoliang
        for (label j = maxCount_[i]; j > 0; j--)
        {
            Map<vector> mapNormals = 
                distribute_.getDatafromOtherProc(zone, interCellNormals);

            volVectorField avgNormTmp = interCellNormals;

            forAll (markers_, cellI)
            {
                if (!zone[cellI]) continue;
                if (counts_[cellI] != j) continue;

                avgNormTmp[cellI] = vector::zero;

                point p = mesh.C()[cellI];

                DynamicList<vector> centres;
                List<DynamicList<scalar>> values(index_.size());

                for (const label gblIdx : stencil[cellI])
                {
                    vector n = distribute_.getValue(interCellNormals, mapNormals, gblIdx);

                    if (mag(n) != 0)
                    {
                        n /= mag(n);

                        vector centre = distribute_.getValue(mesh.C(), mapMC, gblIdx);

                        vector dist = centre - p;
                        vector distToSurf = dist & n / mag(n) * n;
                        vector verticalDist = dist - distToSurf;

                        if (mag(verticalDist) < 1e-8)
                        {
                            avgNormTmp[cellI] = n;
                            centres.clearStorage();
                            break;
                        }

                        vector cc = p - verticalDist;
                        centres.append(cc);

                        forAll (index_, k)
                        {
                            values[k].append(n[index_[k]]);
                        }
                    }
                }

                if (centres.capacity() == 0) continue;

                centres.shrink();
                forAll (values, vi)
                    values[vi].shrink();

                // The choice of these two depends, but the results are close.  
                // On coarse mesh, IDeC is better. On fine mesh, LS is better. 
                // Need further investigation. 
                forAll (index_, k)
                {
                    avgNormTmp[cellI][index_[k]] = interp_.IDeCinterp(p, centres, values[k]);
                    //avgNormTmp[cellI][index_[k]] = interp_.LSfitting(p, centres, values[k]);
                }
            }
            normalise(avgNormTmp);
            avgNormTmp.correctBoundaryConditions();

            interCellNormals = avgNormTmp;
            interCellNormals.correctBoundaryConditions();
        }
    }
}

void Foam::pandoraDivNormalCurvature::curvInterpolate
(
    const volVectorField& interfaceCentres,
    const volScalarField& RDF
)
{
    const fvMesh& mesh = cellCurvature_.mesh();
    volScalarField curvature("curvature" ,cellCurvature_);
    curvature.correctBoundaryConditions();

    boolList zone(mesh.nCells(), false);
    forAll (zone, zi)
    {
        if (markers_[zi] == 0)
        {
            zone[zi] = true;
        }
    }

    distribute_.setUpCommforZone(zone, false);
    Map<vector> mapMC = 
        distribute_.getDatafromOtherProc(zone, mesh.C());
    Map<scalar> mapCurv = 
        distribute_.getDatafromOtherProc(zone, curvature);

    const labelListList& stencil = distribute_.getStencil();
    const volScalarField& alpha = mesh.lookupObject<volScalarField>(fieldName_);

    forAll(cellCurvature_, cellI)
    {
        if (!zone[cellI])
        {
            cellCurvature_[cellI] = 0;
            continue;
        }

        if (alpha[cellI] < 0.99 && alpha[cellI] > 0.01)
        {
             cellCurvature_[cellI] = 2.0 / 
                 (2.0 / (cellCurvature_[cellI] + ROOTVSMALL) + RDF[cellI]);
        }
        else
        {
            vector p = interfaceCentres[cellI];

            DynamicList<vector> points;
            DynamicList<scalar> values;

            for (const label gblIdx : stencil[cellI])
            {
               scalar curv = distribute_.getValue(curvature, mapCurv, gblIdx);
               vector cent = distribute_.getValue(mesh.C(), mapMC, gblIdx);

               points.append(cent);
               values.append(curv);
            }

            points.shrink();
            values.shrink();

            // LS gives better results than IDeC. 
            cellCurvature_[cellI] = interp_.IDeCinterp(p, points, values);
            //cellCurvature_[cellI] = interp_.LSfitting(p, points, values);
        }
    }
    cellCurvature_.correctBoundaryConditions();
}

void Foam::pandoraDivNormalCurvature::curvAverage()
{
    labelField count{cellCurvature_.size(), 0};
    scalarField curvatureSum{cellCurvature_.size(), 0.0};

    const fvMesh& mesh = cellCurvature_.mesh();
    const auto& own = mesh.owner();
    const auto& nei = mesh.neighbour();

    for (label iter = 0; iter != nAverage_; ++iter)
    {
        auto faceCurvatureTmp = fvc::interpolate(cellCurvature_);
        const auto& faceCurvature = faceCurvatureTmp.cref();

        forAll (faceCurvature, fid)
        {
            if ((markers_[own[fid]] == 0) && (markers_[nei[fid]] == 0))
            {
                // Curvature okay: face shared by interface cells
                count[own[fid]] += 1;
                curvatureSum[own[fid]] += faceCurvature[fid];

                count[nei[fid]] += 1;
                curvatureSum[nei[fid]] += faceCurvature[fid];
            }
        }

        // Iterate processor boundaries
        const auto& isInterfaceCellBoundary = markers_.boundaryField();
        const auto& faceCurvatureBoundary = faceCurvature.boundaryField();

        forAll(isInterfaceCellBoundary, patchID)
        {
            const auto& isInterfaceCellPatch = isInterfaceCellBoundary[patchID];
            const auto& faceCurvaturePatch = faceCurvatureBoundary[patchID];

            if (isType<processorFvPatch>(isInterfaceCellPatch.patch()))
            {
                // Values of isInterfaceCell do not change here. Ensure at computation of
                // isInterFaceCell that processor neighbour fields are up-to-date (TT)
                const auto& faceToCell = isInterfaceCellPatch.patch().faceCells();
                const auto& nei = isInterfaceCellPatch.patchNeighbourField().cref();

                forAll(isInterfaceCellPatch, I)
                {
                    if ((markers_[faceToCell[I]] == 0) && (nei[I] == 0))
                    {
                        // Curvature okay: face shared by interface cells
                        count[faceToCell[I]] += 1;
                        curvatureSum[faceToCell[I]] += faceCurvaturePatch[I];
                    }
                }
            }
        }

        forAll(count, cid)
        {
            if (count[cid] > 0)
            {
                cellCurvature_[cid] = curvatureSum[cid]/count[cid];
            }
        }

        count = 0;
        curvatureSum = 0.0;

        // Update processor neighbour values of curvature field
        cellCurvature_.correctBoundaryConditions();
    }
}

void Foam::pandoraDivNormalCurvature::curvExtend
(
    const volVectorField& interfaceCentres,
    const volVectorField& interfaceNormals
)
{
    const fvMesh& mesh = cellCurvature_.mesh();
    volScalarField curvature = cellCurvature_;

    boolList zone(mesh.nCells(), false);
    forAll (zone, zi)
    {
        if (markers_[zi] == 1)
        {
            zone[zi] = true;
        }
    }

    distribute_.setUpCommforZone(zone, false);
    Map<vector> mapCentres = 
        distribute_.getDatafromOtherProc(zone, interfaceCentres);
    Map<vector> mapNormals = 
        distribute_.getDatafromOtherProc(zone, interfaceNormals);
    Map<scalar> mapCurvs = 
        distribute_.getDatafromOtherProc(zone, cellCurvature_);

    const labelListList& stencil = distribute_.getStencil();

    forAll (markers_, cellI)
    {
        if (!zone[cellI]) continue;

        point p = mesh.C()[cellI];

        DynamicList<vector> points;
        DynamicList<scalar> values;

        for (const label gblIdx : stencil[cellI])
        {
            vector n = distribute_.getValue(interfaceNormals, mapNormals, gblIdx);

            if (mag(n) != 0)
            {
                n /= mag(n);

                vector centre = 
                    distribute_.getValue(interfaceCentres, mapCentres, gblIdx);
                scalar curv = 
                    distribute_.getValue(cellCurvature_, mapCurvs, gblIdx);

                vector dist = centre - p;
                vector distToSurf = dist & n / mag(n) * n;
                vector verticalDist = dist - distToSurf;

                if (mag(verticalDist) < 1e-8)
                {
                    curvature[cellI] = curv;
                    points.clearStorage();
                    break;
                }

                vector cc = p - verticalDist;

                points.append(cc);
                values.append(curv);
            }
        }

        if (points.capacity() == 0) continue;

        points.shrink();
        values.shrink();

        // The choice of these two depends, but the results are close.  
        // On coarse mesh, LS is better. On fine mesh, IDeC is better. 
        // Update: problem with LS, suggest using IDeC. 
        //curvature[cellI] = interp_.IDWinterp(p, points, values);
        curvature[cellI] = interp_.IDeCinterp(p, points, values);
        //curvature[cellI] = interp_.LSfitting(p, points, values);
    }
    curvature.correctBoundaryConditions();

    cellCurvature_ = curvature;
    cellCurvature_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraDivNormalCurvature::cellCurvature()
{
    const fvMesh& mesh = cellCurvature_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (!meshDb.found(fieldName_))
    {
        FatalErrorInFunction
            << "pandoraDivNormalCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
	    << "Available registered fields are : \n" 
	    << mesh.names() 
            << abort(FatalError);
    }

    reconstructionSchemes& surf =
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct(false);

    const volVectorField& interfaceNormals = surf.normal();
    const volVectorField& interfaceCentres = surf.centre();
    const boolList& interfaceCells = surf.interfaceCell();

    reconstructedDistanceFunction& RDF = 
        mesh.lookupObjectRef<reconstructedDistanceFunction>("RDF");

    // Check if the interface cells are updated
    bool needUpdate = false;
    if (interCells_ != interfaceCells) 
    {
        needUpdate = true;
    }
    reduce(needUpdate, maxOp<bool>());

    if (needUpdate)
    {
        counts_ = 0.0;
        interCells_ = interfaceCells;
        RDF.markCellsNearSurf(interfaceCells, 2);
        cellDistLevel_ == RDF.cellDistLevel();
        nextToInter_ = RDF.nextToInterface();
        distribute_.updateStencil(nextToInter_);
        updateMarkersAndCounts();
    }

    volVectorField gradRDF(fvc::grad(RDF));
    normalise(gradRDF);
    gradRDF.correctBoundaryConditions();

vector sphereCentre(0.005, 0.005, 0.005);
scalar sphereRadius = 0.002; // Sphere radius

#include "error/error_rdf0.hpp"
#include "error/error_rdf1.hpp"
#include "error/error_gradrdf.hpp"

    averagedNormals_ = gradRDF;  
    forAll (averagedNormals_, i)
    {
        if (cellDistLevel_[i] != 0)
        {
            averagedNormals_[i] = vector::zero;
        }
    }
    //normalise(averagedNormals_);
    averagedNormals_.correctBoundaryConditions();

    // Propagate the interface normals to the narrow band
    normalPropagate(needUpdate, averagedNormals_);

#include "error/error_norm1.hpp"
#include "error/error_norm2.hpp"

//#include "error.hpp"

    surfaceVectorField interfaceVec("interfaceVec",fvc::interpolate((averagedNormals_)));

    surfaceVectorField normalVec("normalVec" ,interfaceVec);
    normalise(normalVec);

    // correct contact angle
    //correctContactAngle(normalVec.boundaryFieldRef(), interfaceVec.boundaryFieldRef());

    if (curvFromTr_)
    {
        const fvBoundaryMesh& boundary = mesh.boundary();

        forAll(boundary, patchi)
        {
            fvPatchVectorField& nHatp = averagedNormals_.boundaryFieldRef()[patchi];
            nHatp = normalVec.boundaryFieldRef()[patchi];
        }
    }

    // Face unit interface normal flux
    surfaceScalarField nHatf_ = normalVec & Sf;

    // Simple expression for curvature
    if (curvFromTr_)
    {
        cellCurvature_ = -tr(fvc::grad(averagedNormals_));
    }
    else
    {
        cellCurvature_ = -fvc::div(nHatf_);
    }
    cellCurvature_.correctBoundaryConditions();

#include "error/error_curv0.hpp"

    // Interpolate curvature from cell centres to PLIC centres
    curvInterpolate(interfaceCentres, RDF);

    // Laplace averaging of the curvature
    if (nAverage_ > 0)
    {
        curvAverage();
    }

#include "error/error_curv1.hpp"

    // Extend the interface curvature to the first layer
    curvExtend(interfaceCentres, interfaceNormals);

#include "error/error_curv2.hpp"

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
