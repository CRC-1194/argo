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

#include "pandoraIntegralCurvature.hpp"

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

defineTypeNameAndDebug(pandoraIntegralCurvature, false);
addToRunTimeSelectionTable(pandoraCurvature, pandoraIntegralCurvature, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void pandoraIntegralCurvature::calcKappa
(
    const isoSurfaceTopo& isoTopo,
    const volScalarField& rdf,
    const scalar& iso,
    const scalar& delta_x,
    scalarField& kappa,
    scalarField& area,
    labelList& marker
)
{
    const auto& faceAreas = isoTopo.faceAreas();
    const auto& pointNormals = isoTopo.pointNormals();
    const auto& points = isoTopo.points(); 
    const auto& meshCells = isoTopo.meshCells();
    const auto& magFaceAreas = isoTopo.magFaceAreas();
    const auto& V = mesh().V();
    forAll(isoTopo, faceI)
    {
        const auto& isoFace = isoTopo[faceI]; 
        vector kappaSum (0, 0, 0);
        for (decltype(isoFace.size()) pI0 = 0; pI0 < isoFace.size(); ++pI0)
        {
            auto pI1 = (pI0 + 1) % isoFace.size(); 
            const auto& pG0 = isoFace[pI0]; // Global p0 label
            const auto& pG1 = isoFace[pI1]; // Global p1 label

            const auto& p0 = points[pG0];
            const auto& p1 = points[pG1];

            vector edgeNormal = 0.5 * (pointNormals[pG0] + pointNormals[pG1]); 
            edgeNormal /= mag(edgeNormal);
            vector edgeVector = p1 - p0;
            kappaSum += (edgeVector ^ edgeNormal); 
        }

        label cellId = meshCells[faceI];
        
        area[cellId] = Foam::sqr(magFaceAreas[faceI]);

        marker[cellId] = 1;

        kappa[cellId] = (kappaSum & faceAreas[faceI]) / area[cellId];
        if (mag(kappa[cellId]) > SMALL)
        {
            scalar radius = 2.0 / kappa[cellId] + iso * delta_x;
            kappa[cellId] = 2.0 / radius;
        }
    }
}

scalar pandoraIntegralCurvature::largestAreaCurv
(
    const scalarList& kappaList,
    const scalarList& areaList
)
{
    label which = 0;
    scalar max = -1.0 / SMALL;
    forAll (areaList, i)
    {
        if (areaList[i] > max)
        {
            which = i;
            max = areaList[i];
        }
    }
    return kappaList[which];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
pandoraIntegralCurvature::pandoraIntegralCurvature
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pandoraCurvature{mesh, dict},
    fieldName_{curvatureDict_.get<word>("fieldName")},
    nPropagate_(curvatureDict_.getOrDefault<label>("nPropagate", 3)),
    nAverage_(curvatureDict_.getOrDefault<label>("nAverage", 3)),
    range_(curvatureDict_.getOrDefault<scalar>("range", 1)),
    nLayer_(curvatureDict_.getOrDefault<label>("nLayer", 3))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField& pandoraIntegralCurvature::cellCurvature()
{
    const auto& meshDb = cellCurvature_.mesh().thisDb();
    if (!meshDb.found(fieldName_))
    {
        FatalErrorInFunction
            << "pandoraIntegralCurvature::cellCurvature \n"
            << "Field " << fieldName_ << " not in mesh registry." 
            << abort(FatalError);
    }

    // Constuct isoSurface. 
    const volScalarField& vf = mesh().lookupObject<volScalarField>(fieldName_);

    const volScalarField& rdf = mesh().lookupObject<volScalarField>("RDF");
    const volScalarField& iic = mesh().lookupObject<volScalarField>("isInterfaceCell");

    scalar maxRdf = Foam::max(rdf).value();
    scalar minRdf = Foam::min(rdf).value();
    volScalarField rdfExt = vf;
    forAll (rdf, cellI)
    {
        rdfExt[cellI] = rdf[cellI];

        /*
        //vector c(0.005, 0.005, 0.005);
        //rdfExt[cellI] = 0.002 - mag(mesh().C()[cellI] - c);
        */

        if (iic[cellI] > 0.5)
        {
            continue;
        }

        if (mag(rdf[cellI]) > SMALL)
        {
            continue;
        }

        if (vf[cellI] < 0.5)
        {
            rdfExt[cellI] = minRdf;
        }
        else if (vf[cellI] > 0.5)
        {
            rdfExt[cellI] = maxRdf;
        }
    }

    scalar delta_x = max(pow(mesh().deltaCoeffs(), -1)).value();
    rdfExt /= delta_x;
    rdfExt.rename("rdf");
    if (mesh().time().writeTime())
        rdfExt.write();
    
    volPointInterpolation vpInterp(mesh());

    tmp<pointScalarField> pfTmp = vpInterp.interpolate(rdfExt);
    pointScalarField& pf = pfTmp.ref();

    labelList marker(mesh().nCells(), 0);

    DynamicList<scalarField> kappaList;
    DynamicList<scalarField> areaList;
    for (label i = 0; i < nLayer_; i++)
    {
        //scalar iso = -range_;
        scalar iso = -range_ + i * 2.0 * range_ / (nLayer_ - 1);
        isoSurfaceTopo isoTopo
        (
            mesh(),
            rdfExt,
            pf,
            iso
        );
        isoTopo.write("isoTopo"+Foam::name(i)+".vtk");
        scalarField kappa(mesh().nCells(), 0.0);
        scalarField area(mesh().nCells(), 0.0);
        calcKappa(isoTopo, rdfExt, iso, delta_x, kappa, area, marker);
        kappaList.append(kappa);
        areaList.append(area);
    }

    forAll (cellCurvature_, cellI)
    {
        if (marker[cellI] == 0)
        {
            cellCurvature_[cellI] = 0.0;
            continue;
        }
        else
        {
            scalarList kl(kappaList.size());
            scalarList al(areaList.size());
            forAll (kl, i)
            {
                kl[i] = kappaList[i][cellI];
                al[i] = areaList[i][cellI];
            }
            cellCurvature_[cellI] = largestAreaCurv(kl, al);
            //cellCurvature_[cellI] = 1000.0;
if (vf[cellI] > 1e-8 && vf[cellI] < 1 - 1e-8)
    Info<<cellCurvature_[cellI]<<nl;
        }
    }

    return cellCurvature_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //