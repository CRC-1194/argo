/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 AUTHOR,AFFILIATION
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

Application
    pandoraTestCurvature

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoAdvection.H"

#include "pandora.hpp"
#include "pandoraCurvature.hpp"
#include "pandoraCurvatureExtension.hpp"
#include "pandoraCurvatureRegularisation.hpp"

struct errorNorms
{
    scalar linf;
    scalar l2;
    scalar l1;

    errorNorms() = default;
};

scalar linfNorm(const volScalarField& field)
{
    return max(mag(field)).value();
}

scalar l2Norm(const volScalarField& field)
{
    return Foam::sqrt(sum(pow(field,2))/field.size()).value();
}

scalar l1Norm(const volScalarField& field)
{
    return (sum(field)/field.size()).value();
}

errorNorms computeErrorNorms(const volScalarField& field)
{
    errorNorms errors{};

    errors.linf = linfNorm(field);
    errors.l2 = l2Norm(field);
    errors.l1 = l1Norm(field);

    return errors;
}

void filterCurvatureErrors(volScalarField& errors, const volScalarField& isInterfaceCell)
{
    forAll(errors, cid)
    {
        if (isInterfaceCell[cid] == 0.0)
        {
            errors[cid] = 0.0;
        }
    }
}

void writeErrors(OFstream& stream, const errorNorms& errors, const word& header, const label nCells)
{
    stream << header << nCells << ","
           << errors.linf << ","
           << errors.l2 << ","
           << errors.l1 << "\n";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (!args.found("exactCurvature")) 
        FatalErrorInFunction
            << "Cannot compute curvature errors without exact value." 
            << abort(FatalError);

    word fieldname = "alpha.water";
    args.readIfPresent<word>("fieldName", fieldname);
    dimensionedScalar exactCurvature
    {
        "exactCurvature",
        dimless / dimLength,
        args.get<scalar>("exactCurvature")
    };

    #include "createFields.hpp"
    
    // Mark cells whose curvature influences surface tension
    auto surfaceTensionFacesTmp = fvc::snGrad(alpha);
    const auto& surfaceTensionFaces = surfaceTensionFacesTmp.ref();
    const auto& neighbour = mesh.neighbour();
    const auto& owner = mesh.owner();

    forAll(surfaceTensionFaces, fid)
    {
        if (mag(surfaceTensionFaces[fid]) > 0.0)
        {
            curvatureRequired[owner[fid]] = 1.0;
            curvatureRequired[neighbour[fid]] = 1.0;
        }
    }
    curvatureRequired.write();

    dictionary pandoraDict{mesh.solutionDict().subDict("pandora")};
    pandora pandoraObj{mesh};
    const auto& isInterfaceCell = pandoraObj.isInterfaceCell(alpha);
    autoPtr<pandoraCurvature> curvPtr(pandoraCurvature::New(mesh, pandoraDict));
    autoPtr<pandoraCurvatureRegularisation> regularisationPtr_{
        pandoraCurvatureRegularisation::New(pandoraDict)
    };
    autoPtr<pandoraCurvatureExtension> extensionPtr_{
        pandoraCurvatureExtension::New(pandoraDict)
    }; 

    // Save error norms
    errorNorms curvatureModelErrors{};
    errorNorms curvatureRegularisedErrors{};
    errorNorms curvatureExtensionErrors{};

    // Curvature as given by the curvature model
    volScalarField& cellCurvature = curvPtr->cellCurvature();
    cellCurvature.write();
    curvatureErrorField = (cellCurvature - exactCurvature)/exactCurvature;
    filterCurvatureErrors(curvatureErrorField, isInterfaceCell);
    curvatureModelErrors = computeErrorNorms(curvatureErrorField);
    curvatureErrorField.write();

    // Regularisation
    cellCurvature.rename("curvature_regularised");
    regularisationPtr_->regularise(cellCurvature, isInterfaceCell);
    cellCurvature.write();
    curvatureErrorField = (cellCurvature - exactCurvature)/exactCurvature;
    filterCurvatureErrors(curvatureErrorField, isInterfaceCell);
    curvatureRegularisedErrors = computeErrorNorms(curvatureErrorField);
    curvatureErrorField.rename("curvature_errors_regularised");
    curvatureErrorField.write();

    // Extension
    cellCurvature.rename("curvature_extended");
    extensionPtr_->extend(cellCurvature, isInterfaceCell);
    cellCurvature.write();
    curvatureErrorField = (cellCurvature - exactCurvature)/exactCurvature*curvatureRequired;
    curvatureExtensionErrors = computeErrorNorms(curvatureErrorField);
    curvatureErrorField.rename("curvature_errors_extension");
    curvatureErrorField.write();

    // Save norms
    OFstream curvatureInterfaceData("curvatureInterfaceErrors.csv");
    OFstream curvatureInterfaceRegularisedData("curvatureInterfaceRegularisedErrors.csv");
    OFstream curvatureNarrowbandData("curvatureNarrowbandErrors.csv");
    word header = "NCELLS,LINF,L2,L1\n";
    
    writeErrors(curvatureInterfaceData, curvatureModelErrors, header, mesh.nCells());
    writeErrors(curvatureInterfaceRegularisedData, curvatureRegularisedErrors, header, mesh.nCells());
    writeErrors(curvatureNarrowbandData, curvatureExtensionErrors, header, mesh.nCells());

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
