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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::fieldError<GeomField>::fieldError(scalar errorTolerance)
:
    errorTolerance_(errorTolerance),
    errorValue_(0)
{}

template<typename GeomField>
Foam::fieldError<GeomField>::fieldError(
   const fieldError<GeomField>& copy
)
:
    errorTolerance_(copy.errorTolerance_),
    errorValue_(copy.errorValue_)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::autoPtr<Foam::fieldError<GeomField> >
Foam::fieldError<GeomField>::New(
    const word& name,
    scalar errorTolerance
)
{
    auto* ctorPtr = ScalarConstructorTable(name);

    if (!ctorPtr)
    {
        FatalErrorIn (
            "fieldError::New(const word&, const scalar)"
        )   << "Unknown fieldError type "
            << name << nl << nl
            << "Valid fieldErrors are : " << endl
            << ScalarConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return Foam::autoPtr<Foam::fieldError<GeomField>>(ctorPtr(errorTolerance));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::fieldError<GeomField>::~fieldError()
{}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<typename GeomField>
void Foam::fieldError<GeomField>::operator=(
        const fieldError<GeomField>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs) {
        FatalErrorIn(
            "Foam::fieldError::operator=(const Foam::fieldError&)"
        ) << "Attempted assignment to self" << abort(FatalError);
    }

    errorTolerance_ = rhs.errorTolerance_;
    errorValue_ = rhs.errorValue_;
}

// ************************************************************************* //
