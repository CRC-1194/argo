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
    Class template for composing the output of the volume field error calculation.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#ifndef fieldErrorComposite_H
#define fieldErrorComposite_H

#include "fieldError.H"
#include "PtrList.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldErrorComposite Declaration
\*---------------------------------------------------------------------------*/

template<typename GeomField>
class fieldErrorComposite
:
    public PtrList<fieldError<GeomField> >
{

public:

    TypeName("fieldErrorComposite");

    typedef fieldError<GeomField> fieldErrorType;
    typedef autoPtr<fieldError<GeomField> > fieldErrorPtrType;
    typedef PtrList<fieldError<GeomField> > containerType;
    typedef fieldErrorComposite<GeomField> errorCompositeType;

    static const word errorTypeName();
    static word getErrorType(const dictionary& dict);

    static const word errorToleranceName();
    static scalar getErrorTolerance(const dictionary& dict);

    // Selectors

    template<typename Tolerances>
    static autoPtr<fieldErrorComposite<GeomField> > New(
        const wordList& errorNames,
        const Tolerances& errorTolerances
    );

    static autoPtr<fieldErrorComposite<GeomField> > New(
        const dictionary& errorsDict
    );

    //- Destructor
    virtual ~fieldErrorComposite();

    // Member Functions

    void appendError(const dictionary& errorDict);

    void computeErrors(
        const GeomField& initialField,
        const GeomField& currentField
    );

    scalarField errorValues() const;

    scalarField errorTolerances() const;

    bool areAllErrorsAcceptable() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFieldErrorComposite(GeomField, Name)                     \
                                                                     \
defineNamedTemplateTypeNameAndDebug(GeomField##Error##Name, 0);      \

template<typename GeomField>
Ostream& operator<<(Ostream& s, const fieldErrorComposite<GeomField>& errors)
{
   forAll (errors, I)
   {
       s << errors[I].errorValue() << " ";
   }

   return s;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldErrorComposite.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
