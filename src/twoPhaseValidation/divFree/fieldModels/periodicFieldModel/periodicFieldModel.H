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

#ifndef periodicFieldModel_H
#define periodicFieldModel_H

#include "divFreeFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class periodicFieldModel
    :
        public divFreeFieldModel
{

    // Decorated interface motion function object.
    autoPtr<divFreeFieldModel> motionFieldModelPtr_;

    // Oscillation period.
    scalar period_;

    // Phase shift.
    scalar phaseShift_;

protected:

        // Intialize the decorated motion.
        void initFieldModelPointer
        (
            const Time& time,
            const dictionary& dict
        );

public:

    //- Runtime type information
    TypeName("periodic");

    // Constructors
    
    periodicFieldModel
    (
        const Time& time,
        const dictionary& dict
    );

   // Member functions

    // Compute the velocity sign for the time value.
    int velocitySign(scalar t) const;

    // Compute velocity at point X in time t.
    vector velocity(point X, scalar t) const;

    //- Access the oscillation period.
    scalar getPeriod() const
    {
        return period_;
    }

    //- Access the phase shift.
    scalar getPhaseShift() const
    {
        return phaseShift_;
    }

    //- Access the decorated motion.
    const divFreeFieldModel& getDecoratedMotion() const
    {
        return motionFieldModelPtr_();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
