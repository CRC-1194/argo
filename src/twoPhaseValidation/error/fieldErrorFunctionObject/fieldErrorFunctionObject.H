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
    Evaluates standard advection errors for interface capturing methods.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef fieldErrorFunctionObject_H
#define fieldErrorFunctionObject_H

#include "functionObject.H"
#include "volFields.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class fieldErrorFunctionObject
    :
        public functionObject
{
    typedef DynamicField<scalar> dynamicScalarField;

    // Private data
        bool started_;
        scalar alphaDropletPhase_;

        //- Pointer to the initial volume fraction filed.
        autoPtr<volScalarField> initialFieldPtr_;

        //- Time
        const Time& time_;

        //- Volume fraction field name.
        word fieldName_;

        //- Volume conservation errror.
        dynamicScalarField errorVolume_;

        //- L1 error.
        dynamicScalarField errorL1_;

        //- L1 error normalized.
        dynamicScalarField errorL1norm_;

        //- Elapsed time.
        dynamicScalarField elapsedTime_;

        //- Initial sum of the volume fraction.
        scalar initialFieldSum_;

   // Private Member Functions

        // Error calculation
        void calcErrorVolume();

        void calcErrorL1();

        void calcErrorL1norm();

        const volScalarField& getCurrentField() const;

        void adjustFieldSizes();

        void setElapsedTime();

public:

    //- Runtime type information
    TypeName("fieldErrors");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fieldErrorFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

    //- Destructor
     ~fieldErrorFunctionObject();

    // Member Functions

        //- Access time
        const Time& getTime() const
        {
            return time_;
        }

        //- Access errors
        const dynamicScalarField& getErrorVolume() const
        {
            return errorVolume_;
        }

        const dynamicScalarField& getErrorL1norm() const
        {
            return errorL1norm_;
        }

        const dynamicScalarField& getErrorL1() const
        {
            return errorL1_;
        }

        //- Access elapsed time
        const dynamicScalarField& getElapsedTime() const
        {
            return elapsedTime_;
        }

        //- Read the system calls
         bool read(const dictionary&);

        //- Execute the "executeCalls" at each time-step
         bool execute();

        //- Called at the start of the time-loop
         bool start();

        //- Execute the "endCalls" at the final time-loop
         bool end();

        //- Write, execute the "writeCalls"
         bool write();

        //- Update for changes of mesh
         void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
         void movePoints(const polyMesh&);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
