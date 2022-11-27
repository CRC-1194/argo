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
#include "zoneDistribute.H"
#include "tensor2D.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"
#include "triSurfaceMesh.H"

#include "interpolationSchemes.hpp"

namespace Foam {

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
        )
{
    Info<<"Selecting divNormal curvature"<<nl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

volScalarField& pandoraDivNormalCurvature::cellCurvature()
{
    const fvMesh& mesh = cellCurvature_.mesh();
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

    interpolationSchemes interp(mesh);

    reconstructionSchemes& surf = 
        mesh.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");
    surf.reconstruct(false);

    reconstructedDistanceFunction& RDF = 
        mesh.lookupObjectRef<reconstructedDistanceFunction>("RDF");

    const volVectorField& interfaceNormals = surf.normal();
    const volVectorField& interfaceCentres = surf.centre();

    RDF.markCellsNearSurf(surf.interfaceCell(), 2);
    zoneDistribute& distribute = zoneDistribute::New(mesh);
    const boolList& nextToInter = RDF.nextToInterface();
    distribute.updateStencil(nextToInter);

    const boolList& nextToInter_ = RDF.nextToInterface();
    const volScalarField& cellDistLevel_ = RDF.cellDistLevel();

    /*
    volScalarField rdf = RDF;
    forAll (rdf, i)
        rdf[i] = 0.0;
    rdf.correctBoundaryConditions();
    if (RDF.internalField() == rdf.internalField())
    {
        RDF.constructRDF
        (
            RDF.nextToInterface(),
            interfaceCentres,
            interfaceNormals,
            distribute,
            false
        );
    }
    RDF.correctBoundaryConditions();
    */

volScalarField rdf = RDF;
rdf.rename("rdf");
if (mesh.time().writeTime())
    rdf.write();

    volVectorField gradRDF(fvc::grad(RDF));
    normalise(gradRDF);
    gradRDF.correctBoundaryConditions();

vector sphereCentre(0.005, 0.005, 0.005);
scalar sphereRadius = 0.002; // Sphere radius

//#include "error/error_rdf0.hpp"
//#include "error/error_rdf1.hpp"
//#include "error/error_gradrdf.hpp"

    // Mark interface markers
    volScalarField markers
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
        dimensionedScalar("cellMarkers", dimless, -1)
    );

    // Count interpolated cell's neighbour cell numbers
    volScalarField counts
    (
        IOobject
        (
            "cellCounts",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellCounts", dimless, 0)
    );

    forAll (markers, cellI)
    {
        if (mag(interfaceNormals[cellI]) != 0)
        {
            markers[cellI] = 0;
        }
    }
    markers.correctBoundaryConditions();

    const auto& own = mesh.owner();
    const auto& nei = mesh.neighbour();

    labelList maxCount(nPropagate_);
    for (label i = 0; i < nPropagate_; ++i)
    {
        volScalarField markersTmp = markers;

        forAll (nei, fid)
        {
            if (markers[own[fid]] == i && markers[nei[fid]] == -1 && nextToInter[nei[fid]])
            {
                //markers[nei[fid]] = i + 1;
                markersTmp[nei[fid]] = i + 1;
                counts[nei[fid]]++;
            }
            if (markers[nei[fid]] == i && markers[own[fid]] == -1 && nextToInter[own[fid]])
            {
                //markers[own[fid]] = i + 1;
                markersTmp[own[fid]] = i + 1;
                counts[own[fid]]++;
            }
        }
        //markers.correctBoundaryConditions();
        markersTmp.correctBoundaryConditions();

        auto& meshBoundary = markers.boundaryFieldRef();
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
                    if (neibrValue[j] == i && markers[faceToCell[j]] == -1 && nextToInter[faceToCell[j]])
                    {
                        //markers[faceToCell[j]] = i + 1;
                        markersTmp[faceToCell[j]] = i + 1;
                        counts[faceToCell[j]]++;
                    }
                }
            }
        }
        markersTmp.correctBoundaryConditions();
        counts.correctBoundaryConditions();

        markers = markersTmp;
        markers.correctBoundaryConditions();

        maxCount[i] = max(counts).value();
    }

    volScalarField markers_ = markers;

    averagedNormals_ = gradRDF; 
    forAll (averagedNormals_, i)
    {
        if (markers[i] != 0)
        {
            averagedNormals_[i] = vector::zero;
        }
    }
    averagedNormals_.correctBoundaryConditions();


    // Interface normals propagate. 
    for (label i = 0; i < nPropagate_; ++i)
    {
        boolList zone(mesh.nCells(), false);
        forAll (zone, zi)
        {
            if (!nextToInter[zi]) continue;

            if (markers[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }

        distribute.setUpCommforZone(zone, false);
        Map<vector> mapMC = 
            distribute.getDatafromOtherProc(zone, mesh.C());

        const labelListList& stencil = distribute.getStencil();

        // Perform the interpolation in order from largest to smallest
        // neighbour cell numbers. Guoliang
        for (label j = maxCount[i]; j > 0; --j)
        {
            Map<vector> mapNormals = 
                distribute.getDatafromOtherProc(zone, averagedNormals_);
            volVectorField avgNormTmp = averagedNormals_;

            forAll (markers, cellI)
            {
                if (markers[cellI] != i + 1) continue;
                if (!nextToInter[cellI]) continue;
                if (counts[cellI] != j) continue;

                avgNormTmp[cellI] = vector::zero;

                point p = mesh.C()[cellI];

                DynamicField<vector> centres;
                DynamicField<scalar> valuesX;
                DynamicField<scalar> valuesY;
                DynamicField<scalar> valuesZ;

                for (const label gblIdx : stencil[cellI])
                {
                    vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);

                    if (mag(n) != 0)
                    {
                        n /= mag(n);

                        vector centre = distribute.getValue(mesh.C(), mapMC, gblIdx);

                        vector dist = centre - p;
                        vector distToSurf = dist & n / mag(n) * n;
                        vector verticalDist = dist - distToSurf;

                        vector cc = p - verticalDist;

                        centres.append(cc);
                        valuesX.append(n.x());
                        valuesY.append(n.y());
                        valuesZ.append(n.z());
                    }
                }

                centres.shrink();
                valuesX.shrink();
                valuesY.shrink();
                valuesZ.shrink();

                if (centres.capacity() == 0)
                {
                    avgNormTmp[cellI] = vector::zero;
                }

                else if (centres.capacity() == 1)
                {
                    avgNormTmp[cellI][0] = valuesX[0];
                    avgNormTmp[cellI][1] = valuesY[0];
                    avgNormTmp[cellI][2] = valuesZ[0];
                }

                else if (centres.capacity() >= 2 && centres.capacity() <= 4)
                {
                    avgNormTmp[cellI][0] = interp.IDWinterp(p, centres, valuesX, 1);
                    avgNormTmp[cellI][1] = interp.IDWinterp(p, centres, valuesY, 1);
                    avgNormTmp[cellI][2] = interp.IDWinterp(p, centres, valuesZ, 1);
                }

                else
                {
                    avgNormTmp[cellI][0] = interp.IDeCinterp(p, centres, valuesX, 1);
                    avgNormTmp[cellI][1] = interp.IDeCinterp(p, centres, valuesY, 1);
                    avgNormTmp[cellI][2] = interp.IDeCinterp(p, centres, valuesZ, 1);
                }
            }
            normalise(avgNormTmp);
            avgNormTmp.correctBoundaryConditions();

            averagedNormals_ = avgNormTmp;
            averagedNormals_.correctBoundaryConditions();
        }
    }

//#include "error/error_norm1.hpp"
//#include "error/error_norm2.hpp"

    cellCurvature_ = -tr(fvc::grad(averagedNormals_));
    //cellCurvature_ = -fvc::div(averagedNormals_);
    cellCurvature_.correctBoundaryConditions();

//#include "error/error_curv0.hpp"

volScalarField curv1 = cellCurvature_;
curv1.rename("curv1");
if (mesh.time().writeTime())
    curv1.write();


    forAll (cellCurvature_, cellI)
    {
        if (markers[cellI] == 0)
        {
            cellCurvature_[cellI] = 2.0 / 
                (2.0 / (cellCurvature_[cellI] + ROOTVSMALL) + RDF[cellI]);
        }
        else
        {
            cellCurvature_[cellI] = 0;
        }
    }
    cellCurvature_.correctBoundaryConditions();
/*
*/

/*
volScalarField curv2 = cellCurvature_;
curv2.rename("curv2");
if (mesh.time().writeTime())
    curv2.write();








{
    // Interpolate curvature from cell centres to PLIC centres
    volScalarField curvature("curvature" ,cellCurvature_);
    curvature.correctBoundaryConditions();

    boolList zone(mesh.nCells(), false);
    forAll (zone, zi)
    {
        if (markers[zi] == 0)
        {
            zone[zi] = true;
        }
    }

    distribute.setUpCommforZone(zone, false);
    Map<vector> mapMC = 
        distribute.getDatafromOtherProc(zone, mesh.C());
    Map<scalar> mapCurv = 
        distribute.getDatafromOtherProc(zone, curvature);

    const labelListList& stencil = distribute.getStencil();

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

           DynamicField<vector> points;
           DynamicField<scalar> values;

           for (const label gblIdx : stencil[cellI])
           {
               vector cent = distribute.getValue(mesh.C(), mapMC, gblIdx);
               scalar curv = distribute.getValue(curvature, mapCurv, gblIdx);

               points.append(cent);
               values.append(curv);
           }

           points.shrink();
           values.shrink();

           if (points.capacity() == 0)
           {
               cellCurvature_[cellI] = 0.0;
           }

           else if (points.capacity() == 1)
           {
               cellCurvature_[cellI] = values[0];
           }

           else if (points.capacity() >= 2 && points.capacity() <= 4)
           {
               cellCurvature_[cellI] = interp.IDWinterp(p, points, values, 1);
           }

           else
           {
               cellCurvature_[cellI] = interp.IDeCinterp(p, points, values, 1);
           }
       }
    }
    cellCurvature_.correctBoundaryConditions();
}
*/

{
    labelField count{cellCurvature_.size(), 0};
    scalarField curvatureSum{cellCurvature_.size(), 0.0};

    for (label iter = 0; iter != nAverage_; ++iter)
    {
        auto faceCurvatureTmp = fvc::interpolate(cellCurvature_);
        const auto& faceCurvature = faceCurvatureTmp.cref();

        forAll (faceCurvature, fid)
        {
            if ((markers[own[fid]] == 0) && (markers[nei[fid]] == 0))
            {
                // Curvature okay: face shared by interface cells
                count[own[fid]] += 1;
                curvatureSum[own[fid]] += faceCurvature[fid];

                count[nei[fid]] += 1;
                curvatureSum[nei[fid]] += faceCurvature[fid];
            }
        }

        // Iterate processor boundaries
        const auto& isInterfaceCellBoundary = markers.boundaryField();
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
                    if ((markers[faceToCell[I]] == 0) && (nei[I] == 0))
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

//#include "error/error_curv1.hpp"
/*

volScalarField curv3 = cellCurvature_;
curv3.rename("curv3");
if (mesh.time().writeTime())
    curv3.write();
*/

{
    for (label i = 0; i < 1; ++i)
    {
        volScalarField curvature = cellCurvature_;

        boolList zone(mesh.nCells(), false);
        forAll (zone, zi)
        {
            if (markers[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }

        distribute.setUpCommforZone(zone, false);
        Map<vector> mapCentres = 
            distribute.getDatafromOtherProc(zone, interfaceCentres);
        Map<vector> mapNormals = 
            distribute.getDatafromOtherProc(zone, interfaceNormals);
        Map<scalar> mapCurvs = 
            distribute.getDatafromOtherProc(zone, cellCurvature_);

        const labelListList& stencil = distribute.getStencil();

        forAll (markers, cellI)
        {
            if (!zone[cellI]) continue;

            point p = mesh.C()[cellI];

            DynamicField<vector> points;
            DynamicField<scalar> values;

            for (const label gblIdx : stencil[cellI])
            {
                vector n = distribute.getValue(interfaceNormals, mapNormals, gblIdx);

                if (mag(n) != 0)
                {
                    n /= mag(n);

                    vector centre = 
                        distribute.getValue(interfaceCentres, mapCentres, gblIdx);
                    scalar curv = 
                        distribute.getValue(cellCurvature_, mapCurvs, gblIdx);

                    vector dist = centre - p;
                    vector distToSurf = dist & n / mag(n) * n;
                    vector verticalDist = dist - distToSurf;

                    vector cc = p - verticalDist;

                    points.append(cc);
                    values.append(curv);
                }
            }

            points.shrink();
            values.shrink();

            if (points.capacity() == 0)
            {
                curvature[cellI] = 0.0;
            }

            else if (points.capacity() == 1)
            {
                curvature[cellI] = values[0];
            }

            else if (points.capacity() >= 2 && points.capacity() <= 4)
            {
                curvature[cellI] = interp.IDWinterp(p, points, values, 1);
            }

            else
            {
                //curvature[cellI] = interp.IDWinterp(p, points, values, 1);
                curvature[cellI] = interp.IDeCinterp(p, points, values, 1);
            }
        }
        curvature.correctBoundaryConditions();

        cellCurvature_ = curvature;
        cellCurvature_.correctBoundaryConditions();
    }
}

//#include "error/error_curv2.hpp"
/*

volScalarField curv4 = cellCurvature_;
curv4.rename("curv4");
if (mesh.time().writeTime())
    curv4.write();
*/








markers.rename("cellMarker");
if (mesh.time().writeTime())
    markers.write();

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
