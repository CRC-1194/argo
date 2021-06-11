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
    Foam::oscillatingDropletFunctionObject

SourceFiles
    oscillatingDropletFunctionObject.C

Description
    Evaluates the kinetic energy, the surface energy, the total energy and the
    bounding box (aka semi axes) of an oscillating droplet for each time step.
    From this data the period of a single oscillation and the decay of the amplitude
    can be computed and compared to analytical solutions.

Author
    Tobias Tolle tolle@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "oscillatingDropletFunctionObject.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(oscillatingDropletFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, oscillatingDropletFunctionObject, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
const triSurface& oscillatingDropletFunctionObject::getFront() const
{
    auto U = hydrodynamicFunctionObject::getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& front = U.mesh().lookupObject<triSurface>("front");

    return front;
}

const dictionary& oscillatingDropletFunctionObject::getTransportProperties() const
{
    auto U = hydrodynamicFunctionObject::getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& transportProperties =
        U.db().lookupObject<dictionary>("transportProperties");

    return transportProperties;
}

void oscillatingDropletFunctionObject::computeKineticEnergy()
{
    dimensionedScalar kineticEnergyDroplet{};
    dimensionedScalar kineticEnergyTotal{};

    // NOTE: only computes the energy of the droplet phase
    auto U = hydrodynamicFunctionObject::getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    auto markerField = hydrodynamicFunctionObject::getCurrentField<scalar, fvPatchField, volMesh>(markerFieldName_);
    const auto& cellVolume = U.mesh().V();

    // Assumption: the droplet phase has a higher density than the ambient phase
    const auto& transportProperties = getTransportProperties();
    dimensionedScalar rho0 = transportProperties.subDict("air").get<dimensionedScalar>("rho");
    dimensionedScalar rho1 = transportProperties.subDict("water").get<dimensionedScalar>("rho");

    auto rhoDroplet = (rho0 > rho1) ? rho0 : rho1;

    if (alphaValueDispersedPhase_ == 1.0)
    {
        kineticEnergyDroplet = 0.5*rhoDroplet.value()*sum(markerField*magSqr(U)*cellVolume);
    }
    else if (alphaValueDispersedPhase_ == 0.0)
    {
        kineticEnergyDroplet = 0.5*rhoDroplet.value()*sum((1.0 - markerField)*magSqr(U)*cellVolume);
    }
    else
    {
        FatalErrorIn("oscillatingDropletFunctionObject::computeKineticEnergy()")
            << "Value " << alphaValueDispersedPhase_ << " is invalid. "
            << "Choose either 1.0 or 0.0" << endl;
    }

    const auto& rhoField = hydrodynamicFunctionObject::getCurrentField<scalar, fvPatchField, volMesh>("rho");
    kineticEnergyTotal = sum(0.5*rhoField*magSqr(U)*cellVolume);

    scalarMetrics_[kineticEnergyDroplet_][timeIndex()] = kineticEnergyDroplet.value(); 
    scalarMetrics_[kineticEnergyTotal_][timeIndex()] = kineticEnergyTotal.value(); 
}

void oscillatingDropletFunctionObject::computeSurfaceEnergy()
{
    const auto& transportProperties = getTransportProperties();
    const dimensionedScalar sigma = transportProperties.get<dimensionedScalar>("sigma");

    const auto& front = getFront();
    const auto& faceAreas = front.magSf();

    scalar surfaceArea = 0.0;

    forAll(faceAreas, I)
    {
        surfaceArea += faceAreas[I];
    }

    scalarMetrics_[surfaceEnergy_][timeIndex()] = sigma.value()*surfaceArea;
}

void oscillatingDropletFunctionObject::computeTotalEnergy()
{
    scalarMetrics_[totalEnergy_][timeIndex()] = scalarMetrics_[kineticEnergyDroplet_][timeIndex()]
                                            + scalarMetrics_[surfaceEnergy_][timeIndex()];
}

void oscillatingDropletFunctionObject::computeSemiAxes()
{
    // Look up front...
    auto U = hydrodynamicFunctionObject::getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& front = U.mesh().lookupObject<triSurface>("front");

    // This should be a member function of triSurface or triSurfaceFront
    // returning the bounding box
    point minimum{};
    point maximum{};

    auto vertices = front.points();

    for (direction I = 0; I < 3; I++)
    {
        minimum[I] = min(vertices.component(I));
        maximum[I] = max(vertices.component(I));
    }

    vectorMetrics_[semiAxes_][timeIndex()] = 0.5*(maximum - minimum);
}

void oscillatingDropletFunctionObject::evaluateMetrics()
{
    computeKineticEnergy();
    computeSurfaceEnergy();
    computeTotalEnergy();
    computeSemiAxes();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
oscillatingDropletFunctionObject::oscillatingDropletFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    hydrodynamicFunctionObject{name, time, dict},
    velocityFieldName_{dict.get<word>("velocityFieldName")},
    kineticEnergyDroplet_{"kinetic_energy_droplet"},
    kineticEnergyTotal_{"kinetic_energy_total"},
    surfaceEnergy_{"surface_energy"},
    totalEnergy_{"total_energy"},
    semiAxes_{"semi-axes-"}
{
    addMetric<scalar>(kineticEnergyDroplet_);
    addMetric<scalar>(kineticEnergyTotal_);
    addMetric<scalar>(surfaceEnergy_);
    addMetric<scalar>(totalEnergy_);
    addMetric<vector>(semiAxes_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool oscillatingDropletFunctionObject::read(const dictionary& dict)
{
    hydrodynamicFunctionObject::read(dict);
    velocityFieldName_ = word{dict.get<word>("velocityFieldName")};

    return execute();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
