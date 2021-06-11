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

#ifndef risingBubbleFunctionObject_H
#define risingBubbleFunctionObject_H

#include <string>

#include "dictionary.H"
#include "volFields.H"

#include "hydrodynamicFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class risingBubbleFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class risingBubbleFunctionObject
:
    //public functionObject
    public hydrodynamicFunctionObject
{

    // Private data
    word velocityFieldName_;

    // String keys used to access the data stored in a map
    std::string averageVelocityDirection_;
    std::string averageVelocityMagnitude_;


    // Private Member Functions
    virtual void computeAverageVelocity();
    virtual void evaluateMetrics();


public:
    
    // Runtime type information
    TypeName("risingBubble")

    // Constructors
    risingBubbleFunctionObject
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );


    //- Destructor
    ~risingBubbleFunctionObject(){};

    virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
