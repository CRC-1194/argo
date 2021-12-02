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


#include "meshMagnitude.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<typename GeomField>
template<typename Tolerances>
Foam::autoPtr<Foam::fieldErrorComposite<GeomField> >
Foam::fieldErrorComposite<GeomField>::New(
    const wordList& errorTypes,
    const Tolerances& errorTolerances
)
{
    autoPtr<errorCompositeType> errorCompositePtr (
         new errorCompositeType()
    );

    errorCompositeType& errorComposite = errorCompositePtr();

    errorComposite.resize(errorTypes.size());

    // Adjust the errorTolerances number to the number of errors.
    // All tolerances that are not defined have SMALL as the default value.
    Tolerances errorTolerancesResized (errorTolerances);
    errorTolerancesResized.resize(errorTypes.size());
    errorTolerancesResized = SMALL;

    forAll(errorTolerances, I)
    {
        errorTolerancesResized[I] = errorTolerances[I];
    }

    // Use errorTypes and errorTolerances to initialize errors.
    forAll (errorTypes, I)
    {
        fieldErrorPtrType errorPtr =
            fieldErrorType::New(errorTypes[I], errorTolerancesResized[I]);

        errorComposite.set(I, errorPtr.ptr());
    }

    return errorCompositePtr;
}

template<typename GeomField>
Foam::autoPtr<Foam::fieldErrorComposite<GeomField> >
Foam::fieldErrorComposite<GeomField>::New(
    const dictionary& errorsDict
)
{
    autoPtr<errorCompositeType> errorCompositePtr (
         new errorCompositeType()
    );

    DynamicList<word> errorTypes;
    DynamicList<scalar> errorTolerances;

    // Loop over error sub dictionaries, and fill errorTypes and
    // errorTolerances.
    forAllConstIter(dictionary, errorsDict, iter)
    {
        const entry& item = *iter;

        if (item.isDict())
        {
            const dictionary& errorDict = errorsDict.subDict(item.keyword());
            errorTypes.append(getErrorType(errorDict));
            errorTolerances.append(getErrorTolerance(errorDict));
         }
    }

    return fieldErrorComposite<GeomField>::New(errorTypes, errorTolerances);
}

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

template<typename GeomField>
const Foam::word Foam::fieldErrorComposite<GeomField>::errorTypeName()
{
    static word errorTypeName = "errorType";
    return errorTypeName;
}

template<typename GeomField>
Foam::word
Foam::fieldErrorComposite<GeomField>::getErrorType(const Foam::dictionary& errorDict)
{
    return errorDict.get<word>(errorCompositeType::errorTypeName());
}

template<typename GeomField>
const Foam::word Foam::fieldErrorComposite<GeomField>::errorToleranceName()
{
    static word errorToleranceName = "errorTolerance";
    return errorToleranceName;
}

template<typename GeomField>
Foam::scalar
Foam::fieldErrorComposite<GeomField>::getErrorTolerance(const dictionary& errorDict)
{
    return readScalar(
        errorDict.lookup(
            errorCompositeType::errorToleranceName()
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::fieldErrorComposite<GeomField>::~fieldErrorComposite()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename GeomField>
void Foam::fieldErrorComposite<GeomField>::computeErrors(
    const GeomField& initialField,
    const GeomField& currentField
)
{
    containerType& errors = *this;

    forAll (errors, I)
    {
        errors[I].computeError(initialField, currentField);
    }
}

template<typename GeomField>
void Foam::fieldErrorComposite<GeomField>::appendError(
        const Foam::dictionary& errorDict
)
{
    word errorType(getErrorType(errorDict));
    scalar errorTolerance(getErrorTolerance(errorDict));

    fieldErrorPtrType fieldErrorPtr(
        fieldErrorType::New(errorType, errorTolerance)
    );

    this->resize(this->size() + 1);
    this->set(this->size() - 1, fieldErrorPtr.ptr());
}

template<typename GeomField>
Foam::scalarField Foam::fieldErrorComposite<GeomField>::errorValues() const
{
    scalarField errorValues(this->size());

    containerType& errors = *this;

    forAll (errors, I)
    {
        errorValues[I] = errors[I].errorValue();
    }

    return errorValues;
}

template<typename GeomField>
Foam::scalarField Foam::fieldErrorComposite<GeomField>::errorTolerances() const
{
    scalarField errorTolerances(this->size());

    containerType& errors = *this;

    forAll (errors, I)
    {
        errorTolerances[I] = errors[I].errorTolerance();
    }

    return errorTolerances;
}

template<typename GeomField>
bool Foam::fieldErrorComposite<GeomField>::areAllErrorsAcceptable() const
{
    containerType& errors = *this;

    forAll (errors, I)
    {
        if (! errors[I].isErrorAcceptable())
        {
            return false;
        }
    }

    return true;
}

// ************************************************************************* //
