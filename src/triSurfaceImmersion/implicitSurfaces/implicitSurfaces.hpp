/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::TriSurfaceImmersion::

Description
    Implicit surfaces described by a level set.

\*---------------------------------------------------------------------------*/

#ifndef implicitSurfaces_H
#define implicitSurfaces_H

#include "ITstream.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "vector.H"

namespace Foam::TriSurfaceImmersion
{
/*---------------------------------------------------------------------------*\
                    Class implicitSurface Declaration
\*---------------------------------------------------------------------------*/

class implicitSurface
{
public:
    // Static Data Members
    TypeName("implicitSurface");

    declareRunTimeSelectionTable(
        autoPtr, implicitSurface, ITstream, (ITstream & is), (is));

    declareRunTimeSelectionTable(autoPtr,
        implicitSurface,
        Dictionary,
        (const dictionary& configDict),
        (configDict));


    // Constructors
    //- Default constructor
    implicitSurface() = default;

    //- Construct from ITstream
    explicit implicitSurface(ITstream&){};

    //- Construct from dictionary
    explicit implicitSurface(const dictionary&){};


    // Selectors
    static autoPtr<implicitSurface> New(const word& name, ITstream& is);
    static autoPtr<implicitSurface> New(const dictionary& configDict);


    // Member functions
    //- Returns scalar value at point x
    virtual scalar value(const vector&) const = 0;

    //- Returns scalar value at point x
    virtual scalar operator()(const vector&) const = 0;

    //- Return gradient at point x
    virtual vector grad(const vector&) const = 0;

    //- Return the volume enclosed by the surface.
    //  If the surface is not closed, a negative value is returned.
    virtual scalar volume() const;
};


class plane : public implicitSurface
{
    //- Reference point determining the postion of the plane
    vector position_;

    //- Normal orientation of the plane
    vector normal_;

public:

    // Static Data Members
    TypeName("plane");


    // Constructors
    //- Construct from reference point and normal vector
    plane(vector position, vector normal);

    //- Contruct from ITstream
    explicit plane(ITstream& is);

    //- Construct from dictionary
    explicit plane(const dictionary& configDict);

    
    // Member functions
    scalar value(const vector& x) const override;

    scalar operator()(const vector& x) const override;

    vector grad(const vector& x) const override;

    //- Return plane's reference point
    vector position() const;

    //- Return plane unit normal
    vector normal() const;
};


class sphere : public implicitSurface
{
    //- Sphere centre
    vector center_;

    //- Sphere radius
    scalar radius_;

public:

    // Static Data Members
    TypeName("sphere");


    // Constructors
    //- Construct from centre and radius
    sphere(vector center, scalar radius);

    //- Construct from ITstream
    explicit sphere(ITstream& is);

    //- Construct from dictionary
    explicit sphere(const dictionary& configDict);

    
    // Member functions
    scalar value(const vector& x) const override;

    scalar operator()(const vector& x) const override;

    vector grad(const vector& x) const override;

    scalar volume() const override;

    //- Return sphere's centre
    vector center() const;

    //- Return sphere's radius
    scalar radius() const;
};


class ellipsoid : public implicitSurface
{
    //- Ellipsoid centre
    vector center_;

    //- Ellipsoid semi axes in x, y and z direction
    vector axes_;

    //- Ellipsoid semi axes squared. Used for Computation.
    vector axesSqr_;

    //- Compute the squared semi axes
    void setAxesSqr(const vector& axes);

public:

    // Static Data members
    TypeName("ellipsoid");


    // Constructors
    //- Construct from centre and vector of semi axes
    ellipsoid(vector center, vector axes);

    //- Construct from ITstream
    explicit ellipsoid(ITstream& is);

    //- Construct from dictionary
    explicit ellipsoid(const dictionary& configDict);


    // Member functions
    scalar value(const vector& x) const override;

    scalar operator()(const vector& x) const override;

    vector grad(const vector& x) const override;

    scalar volume() const override;

    //- Return ellipsoid's centre
    vector center() const;

    //- Return ellipsoid's semi axes
    vector axes() const;
};


class sinc : public implicitSurface
{
    //- Sinc function's origin
    vector origin_;

    //- Amplitude
    scalar amplitude_;

    //- Frequency
    scalar omega_;

public:

    // Static Data members
    TypeName("sinc");


    // Constructors
    //- Construct from orign, amplitude and frequency
    sinc(vector origin, scalar amplitude, scalar omega);

    //- Construct from ITstream
    explicit sinc(ITstream& is);

    //- Construct from dictionary
    explicit sinc(const dictionary& configDict);


    //Member functions
    scalar value(const vector& x) const override;

    scalar operator()(const vector& x) const override;

    vector grad(const vector& x) const override;

    //- Return function's origin
    vector origin() const;

    //- Return amplitude
    scalar amplitude() const;

    //- Return frequency
    scalar omega() const;
};


/* Disbaled for now (TT)
class sincScaled : public implicitSurface // TODO (TM): Scale the amplitude
{
    vector origin_;
    scalar amplitude_;
    scalar omega_;

public:
    TypeName("sincScaled");

    sincScaled(vector origin, scalar amplitude, scalar omega);

    explicit sincScaled(ITstream& is);

    explicit sincScaled(const dictionary& configDict);

    scalar value(const vector& x) const override;

    scalar operator()(const vector& x) const override;

    vector grad(const vector& x) const override;

    vector origin() const;

    scalar amplitude() const;

    scalar omega() const;
};
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
