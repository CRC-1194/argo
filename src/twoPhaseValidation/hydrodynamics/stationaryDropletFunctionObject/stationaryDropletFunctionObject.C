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
    Foam::stationaryDropletFunctionObject

SourceFiles
    stationaryDropletFunctionObject.C

Description
    Evaluates pressure and velocity errors for a stationary droplet test case.
    The evaluation is based on the following publication:
        "A balanced-force algorithm for continuous and sharp interfacial
         surface tension models within a volume tracking framework,
         Francois et al. (2006)"
Author
    Tobias Tolle tolle@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "stationaryDropletFunctionObject.H"

#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stationaryDropletFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, stationaryDropletFunctionObject, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void stationaryDropletFunctionObject::updateCentre(const Time& runTime)
{
    centre_ += runTime.deltaT()*backgroundVelocity_;
}

scalar stationaryDropletFunctionObject::magUMax() const
{
    const volVectorField& U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    return max(mag(U - backgroundVelocity_)).value();
}

void stationaryDropletFunctionObject::saveConstants()
{
    // Look-up of material properties is messy...
    const auto& mesh = getCurrentField<scalar, fvPatchField, volMesh>(pressureFieldName_).mesh();
    //const auto& transportProperties = 
    //    time().lookupObject<dictionary>("transportProperties");
    const auto& U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& transportProperties =
        U.db().lookupObject<dictionary>("transportProperties");
    const dimensionedScalar rhoAir = transportProperties.subDict("air").get<dimensionedScalar>("rho");
    const dimensionedScalar rhoWater = transportProperties.subDict("water").get<dimensionedScalar>("rho");

    scalar densityRatio = rhoWater.value() / rhoAir.value();
    if (densityRatio < 1.0)
    {
        densityRatio = 1 / densityRatio;
    }

    scalar h = max(mag(mesh.delta())).value();

    scalarMetrics_[meshSpacing_][timeIndex()] = h;
    scalarMetrics_[densityRatio_][timeIndex()] = densityRatio;
}

void stationaryDropletFunctionObject::computeVelocityErrors()
{
    const volVectorField& U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);

    scalarMetrics_[L1VelocityError_][timeIndex()] = sum(mag(U - backgroundVelocity_)).value() / U.size();
    scalarMetrics_[L2VelocityError_][timeIndex()] = std::sqrt(sum(magSqr(U - backgroundVelocity_)).value()) / std::sqrt(U.size());
    scalarMetrics_[LInfVelocityError_][timeIndex()] = magUMax();

    // Find cell with maximum relative velocity
    label maxUCell = -1;
    scalar magU = 0.0;

    forAll(U, cellID)
    {
        if (mag(U[cellID] - backgroundVelocity_.value()) > magU)
        {
            maxUCell = cellID;
            magU = mag(U[cellID] - backgroundVelocity_.value());
        }
    }

    scalarMetrics_[maxUCellID_][timeIndex()] = maxUCell;
}

void stationaryDropletFunctionObject::computePressureErrors()
{
    const volScalarField& p = getCurrentField<scalar, fvPatchField, volMesh>(pressureFieldName_);
    const auto& mesh = p.mesh();
    const auto& C = mesh.C();
    //const auto& transportProperties = 
    //    time().lookupObject<dictionary>("transportProperties");
    const auto& U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& transportProperties =
        U.db().lookupObject<dictionary>("transportProperties");
    const scalar sigma{dimensionedScalar{"sigma", transportProperties}.value()};

    updateCentre(mesh.time());

    scalar exactDeltaP = sigma / radius_.value() + SMALL;
    if (is3Dcase(mesh))
    {
        exactDeltaP *= 2;
    }

    dimensionedScalar averagePDroplet = sum(p * pos(radius_ - mag(C - centre_)))
                                / (sum(pos(radius_ - mag(C - centre_))) + SMALL);
    dimensionedScalar averagePAmbient = sum(p * pos(mag(C - centre_) - radius_))
                                / (sum(pos(mag(C - centre_) - radius_)) + SMALL);
    scalarMetrics_[totalPressureError_][timeIndex()] = mag(mag(averagePDroplet - averagePAmbient).value() - exactDeltaP) / exactDeltaP;

    averagePDroplet = sum(p * pos(0.5*radius_ - mag(C - centre_)))
                                / (sum(pos(0.5*radius_ - mag(C - centre_))) + SMALL);
    averagePAmbient = sum(p * pos(mag(C - centre_) - 1.5*radius_))
                                / (sum(pos(mag(C - centre_) - 1.5*radius_)) + SMALL);
    scalarMetrics_[partialPressureError_][timeIndex()] = mag(mag(averagePDroplet - averagePAmbient).value() - exactDeltaP) / exactDeltaP;

    scalar maxDeltaP = mag(max(p) - min(p)).value();
    scalarMetrics_[maxDeltaPressureError_][timeIndex()] = mag(maxDeltaP - exactDeltaP) / exactDeltaP;
}

void stationaryDropletFunctionObject::computeCapillaryNumber()
{
    //const auto& transportProperties = 
    //    time().lookupObject<dictionary>("transportProperties");
    const auto& U = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    const auto& transportProperties =
        U.db().lookupObject<dictionary>("transportProperties");
    const auto& dropletProperties = transportProperties.subDict(dropletPhaseName_);

    const dimensionedScalar rho = dropletProperties.get<dimensionedScalar>("rho");
    const dimensionedScalar nu = dropletProperties.get<dimensionedScalar>("nu");
    const dimensionedScalar sigma{"sigma", transportProperties};
    const auto& markerField = getCurrentField<scalar, fvPatchField, volMesh>("alpha."+dropletPhaseName_);
    const auto& velocity = getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);

    dimensionedScalar UMaxDroplet{max(mag(markerField * (velocity - backgroundVelocity_)))};

    scalarMetrics_[LInfDropletPhaseVelocityError_][timeIndex()]
        = mag(UMaxDroplet.value());

    // FIXME: this calculation is only valid if both phases have the
    // same material properties nu and rho. 
    // The correct way would be to look up the cell index of the maximum
    // velocity and use this cell's material properties to compute
    // the capillary number (TT)
    if (sigma.value() > 0.0)
    {
        scalarMetrics_[capillaryNumber_][timeIndex()]
            = magUMax()*(rho*nu/sigma).value();

        scalarMetrics_[capillaryNumberDroplet_][timeIndex()]
            = (UMaxDroplet * rho * nu / sigma).value();
    }
    else
    {
        scalarMetrics_[capillaryNumber_][timeIndex()]
            = 0.0;

        scalarMetrics_[capillaryNumberDroplet_][timeIndex()]
            = 0.0;
    }
}

void stationaryDropletFunctionObject::computeFrontVolume()
{
    const volScalarField& p = getCurrentField<scalar, fvPatchField, volMesh>(pressureFieldName_);
    const auto& mesh = p.mesh();

    // Enable front volume calculation if Lent is used.
    if (mesh.foundObject<triSurface>("front"))
    {
        const auto& front = mesh.lookupObject<triSurface>("front");
        scalar volume = 0.0;

        // Compute geometric centre
        const auto& vertices = front.localPoints();

        point geometricCentre{0,0,0};

        for (const auto& vertex : vertices)
        {
            geometricCentre += vertex;
        }

        geometricCentre /= vertices.size();

        // Compute volume using tet decomposition
        const auto& trias = front.localFaces();

        for (const auto& tria : trias)
        {
            volume += mag(((vertices[tria[1]] - vertices[tria[0]]) ^ (vertices[tria[2]] - vertices[tria[0]])) & (geometricCentre - vertices[tria[0]]));
        }
        
        // Compute relative error in terms of the front volume of the first time step
        if (timeIndex() == 1)
        {
            scalarMetrics_[frontVolumeError_][timeIndex()] = 0.0;
            initialFrontVolume_ = volume;
        }
        else
        {
            scalarMetrics_[frontVolumeError_][timeIndex()] = mag(volume - initialFrontVolume_)/initialFrontVolume_;
        }
    }
    else
    {
        scalarMetrics_[frontVolumeError_][timeIndex()] = -1.0;
    }
}

void stationaryDropletFunctionObject::computeMassAndKineticEnergy()
{
    dimensionedScalar sumAlpha{0};
    dimensionedScalar sumKinE{0};

    // NOTE: only computes the energy of the droplet phase
    auto U = hydrodynamicFunctionObject::getCurrentField<vector, fvPatchField, volMesh>(velocityFieldName_);
    auto markerField = hydrodynamicFunctionObject::getCurrentField<scalar, fvPatchField, volMesh>(markerFieldName_);
    const auto& cellVolume = U.mesh().V();

    // Assumption: the droplet phase has a higher density than the ambient phase
/*    const auto& transportProperties = getTransportProperties();
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
*/
    const auto& rhoField = hydrodynamicFunctionObject::getCurrentField<scalar, fvPatchField, volMesh>("rho");
    const auto& KinE = 0.5*rhoField*cellVolume*magSqr(U);
   /* 
    forAll(markerField, cellID)
    {
        if(markerField[cellID]>1e-6)
        {
            sumAlpha += markerField[cellID]*cellVolume[cellID];
            sumKinE += 0.5*rhoField[cellID]*cellVolume[cellID]*magSqr(U[cellID]);
        }
    }
    */
    sumAlpha = gSum(markerField);
    sumKinE = gSum(KinE);

    scalarMetrics_[sumAlpha_][timeIndex()] = sumAlpha.value();
    scalarMetrics_[sumKinE_][timeIndex()] = sumKinE.value();
}

void stationaryDropletFunctionObject::evaluateMetrics()
{
    saveConstants();
    computeVelocityErrors();
    computePressureErrors();
    computeCapillaryNumber();
    computeFrontVolume();
    computeMassAndKineticEnergy();
}

bool stationaryDropletFunctionObject::is3Dcase(const fvMesh& mesh)
{
    bool is3D = true;

    forAll(mesh.boundary(), patchI)
    {
        const auto& fvp = mesh.boundary()[patchI];

        if (fvp.type() == "empty")
        {
            is3D = false;
        }
    }

    return is3D;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stationaryDropletFunctionObject::stationaryDropletFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    hydrodynamicFunctionObject{name, time, dict},
    velocityFieldName_{dict.get<word>("velocityFieldName")},
    pressureFieldName_{dict.get<word>("pressureFieldName")},
    dropletPhaseName_{dict.get<word>("dropletPhaseName")},
    radius_{"Radius", dimLength, dict.get<scalar>("radius")},
    centre_{"Centre", dimLength, dict.get<vector>("centre")},
    backgroundVelocity_{"uniformU", dimLength/dimTime,
        vector(dict.lookupOrDefault("backgroundVelocity", vector{0,0,0}))},
    meshSpacing_{"h"},
    densityRatio_{"density ratio"},
    L1VelocityError_{"L1 velocity error"},
    L2VelocityError_{"L2 velocity error"},
    LInfVelocityError_{"Linf velocity error"},
    LInfDropletPhaseVelocityError_{"Linf droplet velocity error"},
    totalPressureError_{"total pressure error"},
    partialPressureError_{"partial pressure error"},
    maxDeltaPressureError_{"max delta P error"},
    capillaryNumber_{"capillary number (whole domain)"},
    capillaryNumberDroplet_{"capillary number (droplet)"},
    frontVolumeError_{"front volume"},
    maxUCellID_{"Max mag(U) cell ID"},
    sumAlpha_{"sumAlpha"},
    sumKinE_{"sumKinE"},
    initialFrontVolume_{0.0}
{
    addMetric<scalar>(meshSpacing_);
    addMetric<scalar>(densityRatio_);
    addMetric<scalar>(L1VelocityError_);
    addMetric<scalar>(L2VelocityError_);
    addMetric<scalar>(LInfVelocityError_);
    addMetric<scalar>(LInfDropletPhaseVelocityError_);
    addMetric<scalar>(totalPressureError_);
    addMetric<scalar>(partialPressureError_);
    addMetric<scalar>(maxDeltaPressureError_);
    addMetric<scalar>(capillaryNumber_);
    addMetric<scalar>(capillaryNumberDroplet_);
    addMetric<scalar>(frontVolumeError_);
    addMetric<scalar>(maxUCellID_);
    addMetric<scalar>(sumAlpha_);
    addMetric<scalar>(sumKinE_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool stationaryDropletFunctionObject::read(const dictionary& dict)
{
    velocityFieldName_ = dict.get<word>("velocityFieldName");
    pressureFieldName_ = dict.get<word>("pressureFieldName");
    radius_ = dimensionedScalar("Radius", dimLength, dict.get<scalar>("radius"));
    centre_ = dimensionedVector("Centre", dimLength, dict.get<vector>("centre"));
    backgroundVelocity_ = dimensionedVector(
        "uniformU", dimLength/dimTime,
        dict.lookupOrDefault<vector>("backgroundVelocity", vector{0,0,0})
    );

    return execute();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
