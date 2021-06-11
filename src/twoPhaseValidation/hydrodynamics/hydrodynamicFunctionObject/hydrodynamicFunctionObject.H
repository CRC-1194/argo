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

#ifndef hydrodynamicFunctionObject_H
#define hydrodynamicFunctionObject_H

#include <string>
#include <unordered_map>

#include "dictionary.H"
#include "DynamicField.H"
#include "functionObject.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hydrodynamicFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class hydrodynamicFunctionObject
:
    public functionObject
{
    template<typename Type>
    using metricDataMap = std::unordered_map<std::string, DynamicField<Type>>;
    using dynamicScalarField = DynamicField<scalar>;
    using dynamicVectorField = DynamicField<vector>;
    // TODO: if necessary, add tensorial data (TT)

    // Private data
    const Time& time_;

    // Private Member Functions
    template<typename Type>
    void adjustFieldSizes(metricDataMap<Type>&);
    void setElapsedTime();

    // Place all evaluation in this function, will be called in execute
    virtual void evaluateMetrics(){};
    
    // Aux functions for writing the results
    fileName assembleFileName() const;
    word scalarMetricHeader() const;
    word vectorMetricHeader() const;
    //template<typename Type>
    //word metricHeader(const metricDataMap<Type>&) const;

protected:

    // This is intended to to store the various metrics computed by
    // the specialized function objects
    // TODO: make the data private and use functions for access? (TT)
    std::string elapsedTime_;
    std::string centreOfMass_;
    word markerFieldName_;
    scalar alphaValueDispersedPhase_;
    metricDataMap<scalar> scalarMetrics_;
    metricDataMap<vector> vectorMetrics_;

    template<typename Type, template<class> class PatchField, typename GeoMesh>
    const GeometricField<Type, PatchField, GeoMesh>& getCurrentField(word) const;

    label timeIndex() const;
    virtual void computeCentreOfMass();

    template<typename Type>
    void addMetric(std::string metricName);

public:

    // Runtime type information
    TypeName("hydrodynamicFunctionObject")


    // Constructors
    hydrodynamicFunctionObject
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );


    //- Destructor
    ~hydrodynamicFunctionObject(){};


    // Member Functions
        const Time& time() const;
        // TODO: why does giving this function the const suffix throw a
        // compiler error?
        scalar physicalTime(label);
    
        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual executeControl behaviour and
        //  forces execution (used in post-processing mode)
        virtual bool execute();

        //- Called at each ++ or += of the time-loop.
        //  postProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        virtual bool read(const dictionary&);


        template<typename Type>
        dimensioned<Type> dispersedPhaseVolumeAverage(const GeometricField<Type, fvPatchField, volMesh>& quantityField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hydrodynamicFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
