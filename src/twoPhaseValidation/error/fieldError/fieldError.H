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
    Class template for field error calculation.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef fieldError_H
#define fieldError_H

#include "scalar.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldError Declaration
\*---------------------------------------------------------------------------*/

template<typename GeomField>
class fieldError
{
    // Private Data

    scalar errorTolerance_;
    scalar errorValue_;

protected:

    void setErrorValue(scalar val)
    {
        errorValue_ = val;
    }

public:

    TypeName("fieldError");

    declareRunTimeSelectionTable
    ( 	  	
        autoPtr,
        fieldError,
        Scalar,
        (
            scalar errorTolerance
        ),
        (errorTolerance)	
	) 		

    // Constructors
        explicit fieldError(scalar errorTolerance=SMALL);

        fieldError(const fieldError&);

    // Selectors
        static autoPtr<fieldError> New(
            const word& name,
            scalar errorTolerance
        );

    //- Destructor

        virtual ~fieldError();

    // Member Functions

        virtual scalar errorTolerance() const
        {
            return errorTolerance_;
        }

        virtual void setErrorTolerance(scalar tol)
        {
            errorTolerance_ = tol;
        }

        virtual scalar errorValue() const
        {
            return errorValue_;
        }

        virtual bool isErrorAcceptable() const
        {
            if (errorValue_ <= errorTolerance_)
            {
                return true;
            }

            return false;
        }

        virtual void computeError
        (
            const GeomField& initialField,
            const GeomField& currentField
        ) = 0;

    // Member Operators

        void operator=(const fieldError&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFieldErrorBase(GeomField)                                 \
defineNamedTemplateTypeNameAndDebug(GeomField##Error, 0);             \
defineTemplateRunTimeSelectionTable(GeomField##Error, Scalar)        \

#define makeFieldError(GeomField, Name)                               \
                                                                      \
defineNamedTemplateTypeNameAndDebug(GeomField##Error##Name, 0);       \
addToRunTimeSelectionTable                                            \
(                                                                     \
    GeomField##Error,                                                 \
    GeomField##Error##Name,                                           \
    Scalar                                                            \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldError.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
