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
#include "zoneDistribute.H"
#include "tensor2D.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"
#include "triSurfaceMesh.H"

#include "interpolationSchemes.hpp"

#define iter 5
#define rr 1

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
        nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 3)), 
        nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 3)), 
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
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (!meshDb.found(fieldName_))
    {
        FatalErrorInFunction
            << "pandoraDivNormalCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
	    << "Available registered fields are : \n" 
	    << mesh().names() 
            << abort(FatalError);
    }

    interpolationSchemes interpolator(mesh());

    //const volVectorField& interfaceNormals = 
    //    mesh().lookupObject<volVectorField>(fieldName_);

    reconstructionSchemes* surf = 
        mesh().getObjectPtr<reconstructionSchemes>("reconstructionScheme");

vector sphereCentre(0.005, 0.005, 0.005);
scalar sphereRadius = 0.002; // Sphere radius

    const volVectorField& interfaceNormals = surf->normal();
    const volVectorField& interfaceCentres = surf->centre();

    volVectorField ifn = interfaceNormals / 
    (
        mag(interfaceNormals) + 
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );
    ifn.correctBoundaryConditions();

    // Mark interface markers
    volScalarField markers
    (
        IOobject
        (
            "cellMarkers",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("cellMarkers", dimless, -1)
    );

    forAll (markers, cellI)
    {
        if (mag(interfaceNormals[cellI]) != 0)
        {
            markers[cellI] = 0;
        }
    }
    markers.correctBoundaryConditions();

    const auto& own = mesh().owner();
    const auto& nei = mesh().neighbour();

    for (label i = 0; i < nPropagate_; ++i)
    {
        forAll (nei, fid)
        {
            if (markers[own[fid]] == i && markers[nei[fid]] == -1)
            {
                markers[nei[fid]] = i + 1;
            }
            if (markers[nei[fid]] == i && markers[own[fid]] == -1)
            {
                markers[own[fid]] = i + 1;
            }
        }
        markers.correctBoundaryConditions();

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
                    if (neibrValue[j] == i && markers[faceToCell[j]] == -1)
                    {
                        markers[faceToCell[j]] = i + 1;
                    }
                }
            }
        }
        markers.correctBoundaryConditions();
    }

    zoneDistribute& distribute = zoneDistribute::New(mesh());

    boolList updateZone(mesh().nCells(), false);
    forAll (updateZone, uzi)
    {
        if (markers[uzi] == -1) continue;
        updateZone[uzi] = true;
    }
    distribute.updateStencil(updateZone);

    averagedNormals_ = interfaceNormals /  
    (
        mag(interfaceNormals) + 
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );
    averagedNormals_.correctBoundaryConditions();

/*
*/
// PLIC normals refinement.
{
    boolList zone(mesh().nCells(), false);
    forAll (zone, zi)
    {
        if (markers[zi] == 0)
        {
            zone[zi] = true;
        }
    }
    distribute.setUpCommforZone(zone, false);
    const labelListList& stencil = distribute.getStencil();

    Map<vector> mapCentres = 
        distribute.getDatafromOtherProc(zone, interfaceCentres);

    Map<vector> mapIfn = 
        distribute.getDatafromOtherProc(zone, ifn);
    Map<vector> mapMC = 
        distribute.getDatafromOtherProc(zone, mesh().C());

    forAll (markers, cellI)
    {
        averagedNormals_[cellI] = Zero;

        if (markers[cellI] != 0) continue;

        point p = mesh().C()[cellI];

        DynamicField<vector> centres;
        DynamicField<scalar> valuesX;
        DynamicField<scalar> valuesY;
        DynamicField<scalar> valuesZ;

        for (const label gblIdx : stencil[cellI])
        {
            vector n = distribute.getValue(ifn, mapIfn, gblIdx);

            if (mag(n) != 0)
            {
                n /= mag(n);

                //vector centre = distribute.getValue(mesh().C(), mapMC, gblIdx);
                vector centre = distribute.getValue(interfaceCentres, mapCentres, gblIdx);

                vector dist = centre - p;
                vector distToSurf = dist & n / mag(n) * n;
                vector verticalDist = dist - distToSurf;
                scalar weight = 1 / max(mag(verticalDist), SMALL);

                vector cc = p - verticalDist;

                centres.append(cc);
                valuesX.append(n.x());
                valuesY.append(n.y());
                valuesZ.append(n.z());
            }
        }

        //averagedNormals_[cellI][0] = interpolator.IDWinterpolate(p, centres, valuesX, rr);
        //averagedNormals_[cellI][1] = interpolator.IDWinterpolate(p, centres, valuesY, rr);
        //averagedNormals_[cellI][2] = interpolator.IDWinterpolate(p, centres, valuesZ, rr);

        //averagedNormals_[cellI][0] = interpolator.IDeCinterpolate(p, centres, valuesX, rr);
        //averagedNormals_[cellI][1] = interpolator.IDeCinterpolate(p, centres, valuesY, rr);
        //averagedNormals_[cellI][2] = interpolator.IDeCinterpolate(p, centres, valuesZ, rr);

        averagedNormals_[cellI][0] = interpolator.LSinterpolate(p, centres, valuesX);
        averagedNormals_[cellI][1] = interpolator.LSinterpolate(p, centres, valuesY);
        averagedNormals_[cellI][2] = interpolator.LSinterpolate(p, centres, valuesZ);
    }

    averagedNormals_ /=  
    (
        mag(averagedNormals_) + 
        dimensionedScalar(
            "SMALL", averagedNormals_.dimensions(), SMALL
        )
    );
    averagedNormals_.correctBoundaryConditions();
}
    
    // Interface normals propagate. 
    for (label i = 0; i < nPropagate_; ++i)
    {
        boolList zone(mesh().nCells(), false);
        forAll (zone, zi)
        {
            if (markers[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }

        distribute.setUpCommforZone(zone, false);
        Map<vector> mapMC = 
            distribute.getDatafromOtherProc(zone, mesh().C());
        Map<vector> mapNormals = 
            distribute.getDatafromOtherProc(zone, averagedNormals_);

        const labelListList& stencil = distribute.getStencil();

        volVectorField avgNormTmp = averagedNormals_;
        forAll (markers, cellI)
        {
            if (markers[cellI] != i + 1) continue;

            avgNormTmp[cellI] = Zero;

            point p = mesh().C()[cellI];

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

                    vector centre = distribute.getValue(mesh().C(), mapMC, gblIdx);

                    vector dist = centre - p;
                    vector distToSurf = dist & n / mag(n) * n;
                    vector verticalDist = dist - distToSurf;
                    scalar weight = 1 / max(mag(verticalDist), SMALL);

                    vector cc = p - verticalDist;

                    centres.append(cc);
                    valuesX.append(n.x());
                    valuesY.append(n.y());
                    valuesZ.append(n.z());
                }
            }

            /*
            avgNormTmp[cellI][0] = interpolator.IDWinterpolate(p, centres, valuesX, rr);
            avgNormTmp[cellI][1] = interpolator.IDWinterpolate(p, centres, valuesY, rr);
            avgNormTmp[cellI][2] = interpolator.IDWinterpolate(p, centres, valuesZ, rr);
            */

            /*
            avgNormTmp[cellI][0] = interpolator.IDeCinterpolate(p, centres, valuesX, rr);
            avgNormTmp[cellI][1] = interpolator.IDeCinterpolate(p, centres, valuesY, rr);
            avgNormTmp[cellI][2] = interpolator.IDeCinterpolate(p, centres, valuesZ, rr);
            */

            avgNormTmp[cellI][0] = interpolator.LSinterpolate(p, centres, valuesX);
            avgNormTmp[cellI][1] = interpolator.LSinterpolate(p, centres, valuesY);
            avgNormTmp[cellI][2] = interpolator.LSinterpolate(p, centres, valuesZ);
            /*
            */
        }
        avgNormTmp.correctBoundaryConditions();

        averagedNormals_ = avgNormTmp / 
        (
            mag(avgNormTmp) + 
            dimensionedScalar(
                "SMALL", avgNormTmp.dimensions(), SMALL
            )
        );
        averagedNormals_.correctBoundaryConditions();
    }

    for (label i = 0; i < nAverage_; i++)
    {
        averagedNormals_ = fvc::average(averagedNormals_);

        averagedNormals_ /=
        (
            mag(averagedNormals_) + 
            dimensionedScalar(
                "SMALL", averagedNormals_.dimensions(), SMALL
            )
        );

        averagedNormals_.correctBoundaryConditions();
    }

    cellCurvature_ == -fvc::div(averagedNormals_);

    const volScalarField& RDF = mesh().lookupObject<volScalarField>("RDF");
    volScalarField rdf = RDF;

/*
*/
    forAll (cellCurvature_, cellI)
    {
        if (markers[cellI] == -1 || markers[cellI] == nPropagate_)
        {
            cellCurvature_[cellI] = 0;
        }
        else
        {
            cellCurvature_[cellI] = 2.0 / 
                (2.0 / (cellCurvature_[cellI] + SMALL) + rdf[cellI] * 1.0);
        }
    }

    cellCurvature_.correctBoundaryConditions();

    //#include "error.hpp"

    return cellCurvature_;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
