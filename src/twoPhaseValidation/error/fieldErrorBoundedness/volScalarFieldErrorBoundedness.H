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
    Class template instantiation for the boundedness field error.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef volScalarFieldErrorBoundedness_H
#define volScalarFieldErrorBoundedness_H

#include "fieldError.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class volScalarFieldErrorBoundedness Declaration
\*---------------------------------------------------------------------------*/

class volScalarFieldErrorBoundedness 
:
    public fieldError<volScalarField>
{

    scalar lowerBound_;
    scalar upperBound_;

public:

    TypeName("E_b");

    // Constructors

        //- Construct from components
        explicit volScalarFieldErrorBoundedness(
            scalar errorTolerance=SMALL,
            scalar lowerBound=0,
            scalar upperBound=1
        );

        //- Construct as copy
        volScalarFieldErrorBoundedness(const volScalarFieldErrorBoundedness&);

    //- Destructor

        virtual ~volScalarFieldErrorBoundedness();

    // Member Functions

        virtual void computeError(
            const volScalarField& initialField,
            const volScalarField& currentField
        );

        void setUpperBound(scalar upperBound)
        {
            upperBound_ = upperBound;
        }
        void setLowerBound(scalar lowerBound)
        {
            lowerBound_ = lowerBound;
        }

        scalar getUpperBound() const
        {
            return upperBound_;
        }

        scalar getLowerBound() const
        {
            return lowerBound_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
