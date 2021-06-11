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
    Foam::hudrodynamicFunctionObject

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


namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<typename Type>
void hydrodynamicFunctionObject::adjustFieldSizes(metricDataMap<Type>& metricMap)
{
    if (!metricMap.empty())
    {
        for(auto& metricField: metricMap)
        {
            // .second accesses the actual field while .first accesses the key
            // of the unordered_map
            metricField.second.resize(time_.timeIndex() + 1);
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<typename Type, template<class> class PatchField, typename GeoMesh>
const GeometricField<Type, PatchField, GeoMesh>& hydrodynamicFunctionObject::getCurrentField(word fieldName) const
{
    // Get the reference to the object registry.
    const objectRegistry& reg
    (
        time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
    );

    if (! reg.foundObject<GeometricField<Type, PatchField, GeoMesh>>(fieldName))
    {
        FatalErrorIn("OscillatingDropletFunctionObject::start()")
            << "Field " << fieldName << " is not registered." << endl;
    }

    // Get the reference to the field.
    const GeometricField<Type, PatchField, GeoMesh>& field =
        reg.lookupObject<GeometricField<Type, PatchField, GeoMesh>>(fieldName);

    return field;
}

template<typename Type>
void hydrodynamicFunctionObject::addMetric(std::string metricName)
{
    // NOTE: since funciton object evaluation is not performed for initial
    // data, the metrics for initial data have to be initialized to some 
    // safe value. Here, simply zero ist used
    if (typeid(Type) == typeid(scalar))
    {
        scalarMetrics_[metricName].resize(1);
        scalarMetrics_[metricName][0] = 0.0;
    }
    else if (typeid(Type) == typeid(vector))
    {
        vectorMetrics_[metricName].resize(1);
        vectorMetrics_[metricName][0] = vector{0,0,0};
    }
}

// * * * * * * * * * * * *  Public Member Functions    * * * * * * * * * * * //
template<typename Type>
dimensioned<Type> hydrodynamicFunctionObject::dispersedPhaseVolumeAverage(const GeometricField<Type, fvPatchField, volMesh>& quantityField) const
{
    auto markerField = getCurrentField<scalar, fvPatchField, volMesh>(markerFieldName_);
    const auto& mesh = markerField.mesh();
    const auto& V = mesh.V();

    dimensioned<Type> average{"zero", dimless, Type::zero};

    if (alphaValueDispersedPhase_ == 1.0)
    {
        average = sum(markerField * V * quantityField) / sum(markerField * V);
    }
    else if (alphaValueDispersedPhase_ == 0.0)
    {
        average = sum((1.0 - markerField) * V * quantityField) / sum((1.0 - markerField) * V);
    }
    else
    {
        FatalErrorIn("hydrodynamicFunctionObject::dispersedPhaseVolumeAverage(...)")
            << "Value " << alphaValueDispersedPhase_ << " is invalid. "
            << "Choose either 1.0 or 0.0" << endl;
    }

    return average;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
