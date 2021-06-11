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

#include "curvatureFunctionObject.hpp"

#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvc.H"
#include "messageStream.H"
#include "volFieldsFwd.H"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(curvatureFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, curvatureFunctionObject, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
const volScalarField& curvatureFunctionObject::getField(const word& fieldName) const
{
    // Get the reference to the object registry.
    const objectRegistry& reg
    {
        time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
    };

    if (! reg.foundObject<volScalarField>(fieldName))
    {
        FatalErrorIn("OscillatingDropletFunctionObject::start()")
            << "Field " << fieldName << " is not registered." << endl;
    }

    // Get the reference to the field.
    const volScalarField& field =
        reg.lookupObject<volScalarField>(fieldName);

    return field;
}

fileName curvatureFunctionObject::assembleFileName() const
{
    fileName dataFileName = time_.rootPath() + "/" + time_.globalCaseName() + "/"
        + "curvatureData.csv";

    return dataFileName;
}

std::vector<scalar> curvatureFunctionObject::filterCurvature() const
{
    std::vector<scalar> curvatureSquashed{};
    const auto& volFraction = getField(vofFieldName_);
    const auto& curvature = getField(curvatureFieldName_);
    std::vector<scalar> curvatureFiltered(volFraction.size(), 0.0);

    const auto& mesh = volFraction.mesh();
    const auto& owner = mesh.owner();
    const auto& neighbor = mesh.neighbour();

    auto gradAlphaFTmp = fvc::snGrad(volFraction);
    const auto& gradAlphaF = gradAlphaFTmp.ref();

    forAll(gradAlphaF, fid)
    {
        if (mag(gradAlphaF[fid] != 0.0))
        {
            curvatureFiltered[owner[fid]] = curvature[owner[fid]];
            curvatureFiltered[neighbor[fid]] = curvature[neighbor[fid]];
        }
    }

    curvatureSquashed.reserve(curvatureFiltered.size()/10);
    for (auto val : curvatureFiltered)
    {
        if (val != 0.0)
        {
            curvatureSquashed.push_back(val);
        }
    }

    return curvatureSquashed;
}

void curvatureFunctionObject::clearDataBuffer()
{
    timeSteps_.clear();
    meanCurvature_.clear();
    standardDeviation_.clear();
    minCurvature_.clear();
    maxCurvature_.clear();
    errorL1_.clear();
    errorL2_.clear();
    errorLinf_.clear();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
curvatureFunctionObject::curvatureFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject{name},
    time_{time},
    vofFieldName_{dict.get<word>("vofFieldName")},
    curvatureFieldName_{dict.get<word>("curvatureFieldName")},
    exactCurvature_{dict.get<scalar>("exactCurvature")}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool curvatureFunctionObject::execute()
{
    auto K = filterCurvature();
    auto meanCurvature = std::accumulate(K.begin(), K.end(), 0.0)/K.size();
    scalar standardDeviation = 0.0;
    for (const auto val : K)
    {
        standardDeviation += std::pow(val - meanCurvature,2.0);
    }
    standardDeviation = std::sqrt(standardDeviation/K.size());

    timeSteps_.push_back(time_.timeOutputValue());
    meanCurvature_.push_back(meanCurvature);
    standardDeviation_.push_back(standardDeviation);
    minCurvature_.push_back(*std::min_element(K.begin(), K.end()));
    maxCurvature_.push_back(*std::max_element(K.begin(), K.end()));

    std::vector<scalar> relCurvatureErrors(K.size());
    scalar l1 = 0.0;
    scalar l2 = 0.0;

    for (unsigned int idx=0; idx != K.size(); ++idx)
    {
        relCurvatureErrors[idx] = mag((K[idx] - exactCurvature_)/exactCurvature_);
        l1 += relCurvatureErrors[idx];
        l2 += relCurvatureErrors[idx]*relCurvatureErrors[idx];
    }

    errorL1_.push_back(l1/relCurvatureErrors.size());
    errorL2_.push_back(std::sqrt(l2/relCurvatureErrors.size()));
    errorLinf_.push_back(*std::max_element(relCurvatureErrors.begin(), relCurvatureErrors.end()));

    return true;
}

bool curvatureFunctionObject::write()
{
    std::ofstream outFile;
    std::string sep = ",";

    if (firstWrite_)
    {
        outFile.open(assembleFileName(), ios_base::out);
        outFile << "time" << sep
                << "mean_curvature" << sep
                << "standard_deviation" << sep
                << "minimum_curvature" << sep
                << "maximum_curvature" << sep
                << "error_l1" << sep
                << "error_l2" << sep
                << "error_linf" << "\n";
        firstWrite_ = false;
    }
    else
    {
        outFile.open(assembleFileName(), ios_base::app);
    }

    for (unsigned int idx=0; idx != timeSteps_.size(); ++idx)
    {
        outFile << timeSteps_[idx] << sep
                << meanCurvature_[idx] << sep
                << standardDeviation_[idx] << sep
                << minCurvature_[idx] << sep
                << maxCurvature_[idx] << sep
                << errorL1_[idx] << sep
                << errorL2_[idx] << sep
                << errorLinf_[idx] << "\n";
    }

    outFile.close();

    clearDataBuffer();

    return true;
}

bool curvatureFunctionObject::read(const dictionary&)
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
