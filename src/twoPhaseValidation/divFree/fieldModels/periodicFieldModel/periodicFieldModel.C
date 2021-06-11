/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2011 Tomislav Maric 
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

Description
    Periodical switch of the velocity sign. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "periodicFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(periodicFieldModel, 0);
addToRunTimeSelectionTable(divFreeFieldModel, periodicFieldModel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void periodicFieldModel::initFieldModelPointer
(
    const Time& time,
    const dictionary& dict
)
{
    const word name = dict.get<word>("baseType"); 
    motionFieldModelPtr_ = divFreeFieldModel::New(name, time, dict); 
}

int periodicFieldModel::velocitySign(scalar time) const
{
    // Remove the phase shift from the time.
    time += phaseShift_;

    // Get the multiplicator.
    auto N = static_cast<int>(std::round(time / period_));

    // Get the local time.
    scalar timeLocal = time - N * period_;

    // If the local time value lies in a new half-period.
    if (timeLocal <= (period_ / 2))
    {
        return 1;
    }

    return -1;
}

vector periodicFieldModel::velocity(point X, scalar time) const
{
    using namespace Foam::constant::mathematical;

    return motionFieldModelPtr_->velocity(X, time) * velocitySign(time);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

periodicFieldModel::periodicFieldModel
(
    const Time& time,
    const dictionary& dict
)
:
    divFreeFieldModel(time, dict),
    motionFieldModelPtr_(),
    period_(dict.lookupOrDefault<scalar>("period", 1.0)),
    phaseShift_(dict.lookupOrDefault<scalar>("phaseShift", 0))
{
    initFieldModelPointer(time, dict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
