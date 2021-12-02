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
    Test functions for the fieldError class.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef fieldErrorTests_H
#define fieldErrorTests_H

namespace Foam {

    template<typename GeomField, template <typename> class FieldError>
    void testFieldError(const GeomField& f, FieldError<GeomField> & e)
    {
        e.computeError(f, f);

        Info << e.type() << " error value: " << e.errorValue() << endl;

        Info << e.type() << " error tolerance: " << e.errorTolerance() << endl;

        if (e.isErrorAcceptable())
        {
            Info << e.type() << ": PASS" << endl;
        }
        else
        {
            Info << e.type() << ": FAIL" << endl;
        }
    }

    template<typename GeomField>
    void selectAndTestFieldError(const GeomField& f, word errorName, scalar errorTolerance)
    {
        autoPtr<fieldError<GeomField> > ePtr =
            fieldError<GeomField>::New(errorName, errorTolerance);

        fieldError<GeomField>& e = ePtr();

        testFieldError(f, e);
    }

    template<typename GeomField>
    void selectAndTestFieldError(const GeomField& f, const Time& time)
    {
        autoPtr<fieldError<GeomField> > ePtr =
            fieldError<GeomField>::New(time);

        fieldError<GeomField>& e = ePtr();

        testFieldError(f, e);
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
