/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2017 Tobias Tolle 
     \\/     M anipulation  |
-------------------------------------------------------------------------------

License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::hydrodynamicFunctionObject

SourceFiles
    hydrodynamicFunctionsObjectI.H
    hydrodynamicFunctionObject.C

Description
    This is the base class for several function objects intended for
    evaluation of hydrodynamic test cases.

Author
    Tobias Tolle tolle@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#include "hydrodynamicFunctionObject.H"

#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(hydrodynamicFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, hydrodynamicFunctionObject, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void hydrodynamicFunctionObject::setElapsedTime()
{
    scalarMetrics_[elapsedTime_][time_.timeIndex()] = time_.timeOutputValue();
}

fileName hydrodynamicFunctionObject::assembleFileName() const
{
    fileName dataFileName = time_.rootPath() + "/" + time_.globalCaseName() + "/"
        + type() + "Results.csv";

    return dataFileName;
}

void hydrodynamicFunctionObject::computeCentreOfMass()
{
    auto markerField = getCurrentField<scalar, fvPatchField, volMesh>(markerFieldName_);
    const auto& C = markerField.mesh().C();

    auto centreOfMass = dispersedPhaseVolumeAverage(C);

    vectorMetrics_[centreOfMass_][timeIndex()] = centreOfMass.value();
}

word hydrodynamicFunctionObject::scalarMetricHeader() const
{
    word header{""};

    for(auto const& metricField: scalarMetrics_)
    {
        header += metricField.first + ",";
    }

    return header;
}

word hydrodynamicFunctionObject::vectorMetricHeader() const
{
    word header{""};

    for(auto const& metricField: vectorMetrics_)
    {
        header += metricField.first + "x" + ",";
        header += metricField.first + "y" + ",";
        header += metricField.first + "z" + ",";
    }

    return header;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

label hydrodynamicFunctionObject::timeIndex() const
{
    return time_.timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrodynamicFunctionObject::hydrodynamicFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject{name},
    time_{time},
    elapsedTime_{"time"},
    centreOfMass_{"centre-of-mass-"},
    markerFieldName_{dict.get<word>("phaseIndicatorFieldName")},
    alphaValueDispersedPhase_{dict.get<scalar>("alphaValueDispersedPhase")}
{
    // Add the time to scalar metrics as it is always required 
    addMetric<scalar>(elapsedTime_);

    // Add centre of mass as it is relevant for multiple test cases
    addMetric<vector>(centreOfMass_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Time& hydrodynamicFunctionObject::time() const
{
    return time_;
}

scalar hydrodynamicFunctionObject::physicalTime(label timeStep) 
{
    scalar pointInTime = scalarMetrics_[elapsedTime_][timeStep];

    return pointInTime;
}

bool hydrodynamicFunctionObject::execute()
{
    adjustFieldSizes(scalarMetrics_);
    adjustFieldSizes(vectorMetrics_);

    setElapsedTime();
    computeCentreOfMass();

    evaluateMetrics();

    return true;
}

bool hydrodynamicFunctionObject::write()
{
    char separator = ',';

    OFstream dataFile(assembleFileName());

    dataFile << scalarMetricHeader()
             << vectorMetricHeader() << endl;

    label timeSteps = scalarMetrics_[elapsedTime_].size();

    // Skip first entry since it is not obtained by evaluation but by
    // default initialization
    for (int index = 1; index < timeSteps; index++)
    {
        for (auto const& metricField: scalarMetrics_)
        {
            dataFile << metricField.second[index] << separator;
        } 

        for (auto const& metricField: vectorMetrics_)
        {
            // Add a column for each component of vector metrics.
            // This simplifies postprocessing (TT).
            const auto& vectorData = metricField.second[index];
            
            for (direction I = 0; I != 3; ++I)
            {
                dataFile << vectorData[I] << separator;
            }
        } 
        
        dataFile << endl;
    }

    return true;
}

bool hydrodynamicFunctionObject::end()
{
    write();

    execute();

    return true;
}


bool hydrodynamicFunctionObject::read(const dictionary& dict)
{
    markerFieldName_ = word{dict.get<word>("phaseIndictorFieldName")};
    alphaValueDispersedPhase_ = dict.get<scalar>("alphaValueDispersedPhase");

    // TODO: purpose of returned value for this function is not known (TT)
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
