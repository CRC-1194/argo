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
    Foam::risingBubbleFunctionObject

SourceFiles
    risingBubbleFunctionObject.C

Description
    Computes the centre of mass of a rising bubble and its volume averaged
    velocity

Author
    Tobias Tolle tolle@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#include "risingBubbleFunctionObject.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(risingBubbleFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, risingBubbleFunctionObject, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void risingBubbleFunctionObject::computeAverageVelocity()
{
    auto U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);

    auto averageU = dispersedPhaseVolumeAverage(U);

    vectorMetrics_[averageVelocityDirection_][timeIndex()] = (averageU
                                             / (mag(averageU) + SMALL)).value();
    scalarMetrics_[averageVelocityMagnitude_][timeIndex()] = (mag(averageU)).value();
}

void risingBubbleFunctionObject::evaluateMetrics()
{
    computeAverageVelocity();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
risingBubbleFunctionObject::risingBubbleFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    hydrodynamicFunctionObject{name, time, dict},
    velocityFieldName_{dict.get<word>("velocityFieldName")},
    averageVelocityDirection_{"average velocity direction"},
    averageVelocityMagnitude_{"average velocity magnitude"}
{
    addMetric<scalar>(averageVelocityMagnitude_);

    addMetric<vector>(averageVelocityDirection_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool risingBubbleFunctionObject::read(const dictionary& dict)
{
    hydrodynamicFunctionObject::read(dict);
    velocityFieldName_ = word{dict.get<word>("velocityFieldName")};

    return execute();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
