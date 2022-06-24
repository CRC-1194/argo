/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

#include "pandoraRDFIsosurfaceCurvature.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "error.H"
#include "dimensionedScalarFwd.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "volFieldsFwd.H"

#include "volPointInterpolation.H"
#include "isoSurfaceTopo.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraRDFIsosurfaceCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraRDFIsosurfaceCurvature, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
pandoraRDFIsosurfaceCurvature::pandoraRDFIsosurfaceCurvature
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pandoraCurvature{mesh, dict},
    fieldName_{curvatureDict_.get<word>("fieldName")},
    nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 3)),
    nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 3))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField& pandoraRDFIsosurfaceCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (!meshDb.found(fieldName_))
    {
        FatalErrorInFunction
            << "pandoraRDFIsosurfaceCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    // Constuct isoSurface. 
    volScalarField vof = mesh().lookupObject<volScalarField>("alpha.dispersed");
    volScalarField rdf = mesh().lookupObject<volScalarField>(fieldName_);
    scalar maxRDF = Foam::max(rdf).value();
    scalar minRDF = Foam::min(rdf).value();

    volScalarField rdfTmp = rdf;
    forAll (rdf, cellI)
    {
        if (mag(rdf[cellI]) > SMALL 
              && vof[cellI] > SMALL 
              && vof[cellI] < 1 - SMALL) continue;

        if (vof[cellI] < SMALL)
            rdfTmp[cellI] = minRDF;
        else if (vof[cellI] > 1 - SMALL)
            rdfTmp[cellI] = maxRDF;
    }
    rdf = rdfTmp;

    rdf.rename("rdf1");
    if(mesh().time().writeTime())
        rdf.write();

    for (label i = 0; i < 1; i++)
        rdf = fvc::average(rdf);

    rdf.rename("rdf2");
    if(mesh().time().writeTime())
        rdf.write();

    volPointInterpolation vpInterp(mesh());

    tmp<pointScalarField> rdfpTmp = vpInterp.interpolate(rdf);
    pointScalarField& rdfp = rdfpTmp.ref();

    scalar delta_x = max(pow(mesh().deltaCoeffs(), -1)).value();

    isoSurfaceTopo RDFIso000(
        mesh(),
        rdf,
        rdfp,
        0//,
        //isoParams
    );
    if(mesh().time().writeTime())
        RDFIso000.write("RDFIso000.vtk");

    isoSurfaceTopo RDFIso001(
        mesh(),
        rdf,
        rdfp,
        delta_x * 0.5//,
        //isoParams
    );
    if(mesh().time().writeTime())
        RDFIso001.write("RDFIso001.vtk");

    isoSurfaceTopo RDFIso002(
        mesh(),
        rdf,
        rdfp,
        - delta_x * 0.5//,
        //isoParams
    );
    if(mesh().time().writeTime())
        RDFIso002.write("RDFIso002.vtk");

    // Get meshCells and normals from isoSurface. 
    const auto& triToCell000 = RDFIso000.meshCells();
    const auto& triNormals000 = RDFIso000.Sf();

    const auto& triToCell001 = RDFIso001.meshCells();
    const auto& triNormals001 = RDFIso001.Sf();

    const auto& triToCell002 = RDFIso002.meshCells();
    const auto& triNormals002 = RDFIso002.Sf();

    volVectorField interfaceNormals
    (
        IOobject
        (
            "interfaceNormals",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(
            "interfaceNormals",
            dimless,
            vector::zero
        )
    );

    forAll(cellCurvature_, cellI)
        cellCurvature_[cellI] = 0.0;

    forAll(triToCell001, i)
    {
        label cellI = triToCell001[i];
        interfaceNormals[cellI] = triNormals001[i];
        //cellCurvature_[cellI] = 2000.0;
    }
    interfaceNormals.rename("interfaceNormals1");
    if(mesh().time().writeTime())
        interfaceNormals.write();

    forAll(triToCell002, i)
    {
        label cellI = triToCell002[i];
        interfaceNormals[cellI] = triNormals002[i];
        //cellCurvature_[cellI] = 2000.0;
    }
    interfaceNormals.rename("interfaceNormals2");
    if(mesh().time().writeTime())
        interfaceNormals.write();

    /*
    forAll(triToCell000, i)
    {
        label cellI = triToCell000[i];
        interfaceNormals[cellI] = triNormals000[i];
        //cellCurvature_[cellI] = 2000.0;
    }
    interfaceNormals.rename("interfaceNormals3");
    if(mesh().time().writeTime())
        interfaceNormals.write();

    vector c(0.005, 0.005, 0.005);
    forAll(interfaceNormals, i)
    {
        vector mc = mesh().C()[i];
        interfaceNormals[i] = c - mc;
    }
    interfaceNormals.rename("interfaceNormals4");
    if(mesh().time().writeTime())
        interfaceNormals.write();

    for (label i = 0; i < nPropagate_; i++)
    {
        auto cc = cellCurvature_;
        forAll(cc, cellI)
        {
            if (cc[cellI] > SMALL)
            {
                const auto& neibrs = mesh().cellCells()[cellI];
                forAll(neibrs, i)
                {
                    label cellJ = neibrs[i];
                    cellCurvature_[cellJ] = 2000.0;
                }
            }
        }
    }
    */

    // Normalize the normal vectors. 
    volVectorField averagedNormals = interfaceNormals /
    (
        mag(interfaceNormals) +
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );
    averagedNormals.rename("normalizedNormals");
    if(mesh().time().writeTime())
        averagedNormals.write();

    // Propagate the normals into the bulk. 
    interfaceNormals == averagedNormals;
    for (label i = 0; i < nPropagate_; i++)
    {
        averagedNormals == fvc::average(averagedNormals);

        averagedNormals /= mag(averagedNormals) + 
            dimensionedScalar("SMALL", averagedNormals.dimensions(), SMALL);

        // Re-set the smoothed normal vectors in interface cells. 
        forAll(interfaceNormals, cellI)
        {
            if ((interfaceNormals[cellI][0] != 0) ||
                (interfaceNormals[cellI][1] != 0) ||
                (interfaceNormals[cellI][2] != 0))
            {
                averagedNormals[cellI] = interfaceNormals[cellI] / 
                    mag(interfaceNormals[cellI]);
            }
        }
    }
    averagedNormals.rename("averagedNormals");
    if(mesh().time().writeTime())
        averagedNormals.write();

    // Smooth the normals
    for (label i = 0; i < nAverage_; i++)
    {
        averagedNormals == fvc::average(averagedNormals);

        // Normalize the normal vectors. 
        averagedNormals /= mag(averagedNormals) +
            dimensionedScalar("SMALL", averagedNormals.dimensions(), SMALL);
    }
    averagedNormals.rename("smoothedNormals");
    if(mesh().time().writeTime())
        averagedNormals.write();

    cellCurvature_ = - fvc::div(averagedNormals);

    /*
    const volScalarField& marker = 
        mesh().lookupObject<volScalarField>("isInterfaceCell");
    forAll(marker, cellI)
        if (marker[cellI] < SMALL)
            cellCurvature_[cellI] = 0.0;

    forAll(cellCurvature_, cellI)
        cellCurvature_[cellI] = 2000.0;
    */

    return cellCurvature_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
