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

#ifndef oscillatingDropletFunctionObject_H
#define oscillatingDropletFunctionObject_H

#include <string>

#include "dictionary.H"
#include "volFields.H"
#include "triSurface.H"

#include "hydrodynamicFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class oscillatingDropletFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class oscillatingDropletFunctionObject
:
    //public functionObject
    public hydrodynamicFunctionObject
{

    // Private data
    word velocityFieldName_;

    // String keys used to access the data stored in a map
    std::string kineticEnergyDroplet_;
    std::string kineticEnergyTotal_;
    std::string surfaceEnergy_;
    std::string totalEnergy_;
    std::string semiAxes_;


    // Private Member Functions
    virtual const triSurface& getFront() const;
    virtual const dictionary& getTransportProperties() const;

    virtual void computeKineticEnergy();
    virtual void computeSurfaceEnergy();
    virtual void computeTotalEnergy();
    virtual void computeSemiAxes();
    virtual void evaluateMetrics();


public:
    
    // Runtime type information
    TypeName("oscillatingDroplet")

    // Constructors
    oscillatingDropletFunctionObject
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );


    //- Destructor
    ~oscillatingDropletFunctionObject(){};

    virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
