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
    Function object for computing divergence-free velocity and flux fields.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef divFreeFunctionObject_H
#define divFreeFunctionObject_H

#include "stringList.H"
#include "pointField.H"
#include "functionObject.H"
#include "divFreeFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               A functionObject adapter for the divFreeModel 
\*---------------------------------------------------------------------------*/

class divFreeFunctionObject
    :
        public functionObject
{
    autoPtr<divFreeFieldModel> divFreeModelPtr_; 

public:

    //- Runtime type information
    TypeName("divFree");

    // Constructors
    divFreeFunctionObject (
        const word& name,
        const Time& time,
        const dictionary& dict
    );

    // Member Functions

    //- Read attribute data from dictionary
    bool read(const dictionary&);

    //- Called at each ++ or += of the time-loop.
    //  postProcess overrides the usual writeControl behaviour and
    //  forces writing always (used in post-processing mode)
    inline virtual bool write() 
    { 
        return true;
    };  

    //- Execute the "executeCalls" at each time-step
    bool execute();

    //- Called at the start of the time-loop
    bool start();

    //- Execute the "endCalls" at the final time-loop
    bool end();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
