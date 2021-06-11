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

#ifndef stationaryDropletFunctionObject_H
#define stationaryDropletFunctionObject_H

#include <string>

#include "dictionary.H"
#include "DynamicField.H"
#include "functionObject.H"
#include "point.H"
#include "volFields.H"

#include "hydrodynamicFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stationaryDropletFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stationaryDropletFunctionObject
:
    public hydrodynamicFunctionObject
{

    // Private data
    word velocityFieldName_;
    word pressureFieldName_;
    word dropletPhaseName_;
    dimensionedScalar radius_;
    dimensionedVector centre_;
    dimensionedVector backgroundVelocity_;

    // Keys to store/access data in map
    std::string meshSpacing_;
    std::string densityRatio_;

    std::string L1VelocityError_;
    std::string L2VelocityError_;
    std::string LInfVelocityError_;
    std::string LInfDropletPhaseVelocityError_;

    std::string totalPressureError_;
    std::string partialPressureError_;
    std::string maxDeltaPressureError_;

    std::string capillaryNumber_;
    std::string capillaryNumberDroplet_;

    std::string frontVolumeError_;
    std::string maxUCellID_;

    scalar initialFrontVolume_;

    // Private Member Functions
    virtual void updateCentre(const Time& runTime);

    virtual scalar magUMax() const;

    virtual void saveConstants();
    virtual void computeVelocityErrors();
    virtual void computePressureErrors();
    virtual void computeCapillaryNumber();
    virtual void computeFrontVolume();
    virtual void evaluateMetrics();

    bool is3Dcase(const fvMesh& mesh);


public:

    // Runtime type information
    TypeName("stationaryDroplet")

    // Constructors
    stationaryDropletFunctionObject
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );


    //- Destructor
    ~stationaryDropletFunctionObject(){};

    // Member functions
    virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
