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

Description
    Implicit surfaces used for NN input data generation. 

\*---------------------------------------------------------------------------*/

#include "implicitSurfaces.hpp"
#include "addToRunTimeSelectionTable.H"

namespace Foam::TriSurfaceImmersion {

// * * * * * * * * * * * * Class implicitSurface  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(implicitSurface, false);
defineRunTimeSelectionTable(implicitSurface, ITstream);
defineRunTimeSelectionTable(implicitSurface, Dictionary);

autoPtr<implicitSurface> implicitSurface::New 
(
    const word& name, 
    ITstream is
)
{
    // Find the constructor pointer for the model in the constructor table.
    ITstreamConstructorTable::iterator cstrIter =
        ITstreamConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == ITstreamConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "AI::implicitSurface::New(const word&, ITstream&&)"
        )   << "Unknown implicitSurface type "
            << name << nl << nl
            << "Valid implicitSurfaces are : " << endl
            << ITstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<implicitSurface>
        (cstrIter()(is));
}


autoPtr<implicitSurface> implicitSurface::New 
(
    const dictionary& configDict
)
{
    const auto name = configDict.get<word>("type");
    // Find the constructor pointer for the model in the constructor table.
    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "Foam::TriSurfaceImmersion::implicitSurface::New(const word&, const dictionary&)"
        )   << "Unknown implicitSurface type "
            << name << nl << nl
            << "Valid implicitSurfaces are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object. 
    return autoPtr<implicitSurface>
        (cstrIter()(configDict));
}


scalar implicitSurface::signedDistance(const vector& x) const
{
    return sign(this->value(x))*mag(x - this->closestPoint(x));
}

// * * * * * * * * * * * * Class plane  * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(plane, false);
addToRunTimeSelectionTable(implicitSurface, plane, ITstream);
addToRunTimeSelectionTable(implicitSurface, plane, Dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plane::plane(vector position, vector normal)
: 
    position_(position), 
    normal_(normal)
{
    normal_ /= Foam::mag(normal_);
}

plane::plane(ITstream is)
{
    is >> position_; 
    is >> normal_; 
}

plane::plane(const dictionary& configDict)
:
    position_{configDict.get<vector>("position")},
    normal_{configDict.get<vector>("normal")}
{
    normal_ /= Foam::mag(normal_);
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar plane::value(const vector& x) const
{
    return Foam::dot(x - position_, normal_);
}

scalar plane::operator()(const vector& x) const
{
    return value(x); 
}

vector plane::closestPoint(const vector& x) const
{
    return x - value(x)*normal_;
}

 vector plane::grad(const vector& x) const
{
    return normal_; 
}

vector plane::position() const
{
    return position_; 
}

vector plane::normal() const
{
    return normal_; 
}

// * * * * * * * * * * * * Class sphere * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sphere, false);
addToRunTimeSelectionTable(implicitSurface, sphere, ITstream);
addToRunTimeSelectionTable(implicitSurface, sphere, Dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphere::sphere(vector center, scalar radius)
    : 
        center_(center), 
        radius_(radius)
{}

sphere::sphere(ITstream is)
{
    is >> center_; 
    is >> radius_; 
}

sphere::sphere(const dictionary& configDict)
:
    center_{configDict.get<vector>("center")},
    radius_{configDict.get<scalar>("radius")}
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar sphere::value(const vector& x) const
{
    return Foam::mag(x - center_) - radius_; 
}

scalar sphere::operator()(const vector& x) const
{
    return value(x); 
}

vector sphere::closestPoint(const vector& x) const
{
    auto delta = x - center_;
    auto dist = mag(delta);

    if (dist > SMALL)
    {
        return center_ + radius_*delta/mag(delta);
    }
    else
    {
        // x coincides with centre: every point on the sphere can be 
        // considered closest point
        return center_ + radius_*vector{1,0,0};
    }
}

vector sphere::grad(const vector& x) const
{
    scalar x0c0 = x[0] - center_[0];
    scalar x1c1 = x[1] - center_[1];
    scalar x2c2 = x[2] - center_[2];

    return vector(x0c0, x1c1, x2c2) / 
        sqrt(x0c0*x0c0 + x1c1*x1c1 + x2c2*x2c2);
}

vector sphere::center() const
{
    return center_; 
}

scalar sphere::radius() const
{
    return radius_; 
}

// * * * * * * * * * * * * Class ellipsoid * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ellipsoid, false);
addToRunTimeSelectionTable(implicitSurface, ellipsoid, ITstream);
addToRunTimeSelectionTable(implicitSurface, ellipsoid, Dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ellipsoid::ellipsoid(vector center, vector axes)
    : 
        center_(center), 
        axes_(axes)
{
    setAxesSqr(axes);
}

ellipsoid::ellipsoid(ITstream is)
{
    is >> center_; 
    is >> axes_; 
    setAxesSqr(axes_);
}

ellipsoid::ellipsoid(const dictionary& configDict)
:
    center_{configDict.get<vector>("center")},
    axes_{configDict.get<vector>("axes")}
{
    setAxesSqr(axes_);
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void ellipsoid::setAxesSqr(const vector& axes)
{
    axesSqr_[0] = Foam::sqr(axes[0]);
    axesSqr_[1] = Foam::sqr(axes[1]);
    axesSqr_[2] = Foam::sqr(axes[2]);
}

scalar ellipsoid::value(const vector& x) const
{
    return Foam::sqr(x[0] - center_[0]) / axesSqr_[0] + 
           Foam::sqr(x[1] - center_[1]) / axesSqr_[1] + 
           Foam::sqr(x[2] - center_[2]) / axesSqr_[2] - 1;
}

scalar ellipsoid::operator()(const vector& x) const
{
    return value(x); 
}

vector ellipsoid::closestPoint(const vector& x) const
{
    return vector{1,0,0};
}

vector ellipsoid::grad(const vector& x) const
{
    return 2*vector
    (
        (x[0] - center_[0])/axesSqr_[0], 
        (x[1] - center_[1])/axesSqr_[1], 
        (x[2] - center_[2])/axesSqr_[2]
    );
}

vector ellipsoid::center() const
{
    return center_; 
}

vector ellipsoid::axes() const
{
    return axes_; 
}

// * * * * * * * * * * * * Class sinc * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sinc, false);
addToRunTimeSelectionTable(implicitSurface, sinc, ITstream);
addToRunTimeSelectionTable(implicitSurface, sinc, Dictionary);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sinc::sinc(vector origin, scalar amplitude, scalar omega)
    : 
        origin_(origin), 
        amplitude_(amplitude), 
        omega_(omega)
{}

sinc::sinc(ITstream is)
{
    is >> origin_; 
    is >> amplitude_; 
    is >> omega_;
}

sinc::sinc(const dictionary& configDict)
:
    origin_{configDict.get<vector>("origin")},
    amplitude_{configDict.get<scalar>("amplitude")},
    omega_{configDict.get<scalar>("omega")}
{}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

scalar sinc::value(const vector& x) const
{
    double r = Foam::sqrt
    (
        Foam::sqr(x[0] - origin_[0]) + 
        Foam::sqr(x[1] - origin_[1]) 
    );

    if (r < std::numeric_limits<double>::min())
        return amplitude_; 
    else 
    {
        return x[2] - origin_[2] - amplitude_ * sin(omega_ * r) / (omega_ * r);
    }
}

scalar sinc::operator()(const vector& x) const
{
    return value(x); 
}

vector sinc::closestPoint(const vector& x) const
{
    return vector{1,0,0};
}

vector sinc::grad(const vector& x) const
{
    const scalar& A = amplitude_; 
    const scalar& O0 = origin_[0];
    const scalar& O1 = origin_[1];

    const scalar& x0 = x[0];
    const scalar& x1 = x[1];

    return vector // Expression calculated in sympy.
    (
        A*(O0 - x0)*(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)),
            A*(O1 - x1)*(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 3.0/2.0)*cos(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))) - (pow(O0 - x0, 2) + pow(O1 - x1, 2))*sin(omega_*sqrt(pow(O0 - x0, 2) + pow(O1 - x1, 2))))/(omega_*pow(pow(O0 - x0, 2) + pow(O1 - x1, 2), 5.0/2.0)),
            1
        );
}

vector sinc::origin() const
{
    return origin_; 
}

scalar sinc::amplitude() const
{
    return amplitude_; 
}

scalar sinc::omega() const
{
    return omega_; 
}

// * * * * * * * * * * * * Class sincScaled * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sincScaled, false);
addToRunTimeSelectionTable(implicitSurface, sincScaled, ITstream);

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sincScaled::sincScaled(vector origin, scalar amplitude, scalar omega)
    : 
        origin_(origin), 
        amplitude_(amplitude), 
        omega_(omega)
{}

sincScaled::sincScaled(ITstream is)
{
    is >> origin_; 
    is >> amplitude_; 
    is >> omega_;
}

sincScaled::sincScaled(const dictionary& configDict)
:
    origin_{configDict.get<vector>("origin")},
    amplitude_{configDict.get<scalar>("amplitude")},
    omega_{configDict.get<scalar>("omega")}
{}
    

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //


scalar sincScaled::value(const vector& x) const // TODO scale the amplitude
{
    double r = Foam::sqrt
    (
        Foam::sqr(x[0] - origin_[0]) + 
        Foam::sqr(x[1] - origin_[1]) 
    );

    if (r < std::numeric_limits<double>::min())
        return amplitude_; 
    else 
    {
        //decltype(auto) factor = [](double A, double z)
        //{
            //if ((z < -A) || (z > A))
                //return 0.; 
            //else
                //return A - 2*pow(z,2)/A + pow(z,4)/pow(A,3);
        //};

        //scalar z = x[2] - origin_[2];
        //return z - amplitude_ * factor(amplitude_,z)*sin(omega_ * r) / (omega_ * r);
        return x[2] - origin_[2] - amplitude_ * sin(omega_ * r) / (omega_ * r);
    }
}

scalar sincScaled::operator()(const vector& x) const
{
    return value(x); 
}

vector sincScaled::closestPoint(const vector& x) const
{
    return vector{1,0,0};
}

vector sincScaled::grad(const vector& x) const
{
    //const scalar& A = amplitude_; 
    ////const scalar& omega = omega_;
    //const scalar& O0 = origin_[0];
    //const scalar& O1 = origin_[1];

    //const scalar& x0 = x[0];
    //const scalar& x1 = x[1];

    // FIXME: insert sympy expression. 
    return vector
    (
        GREAT, GREAT, GREAT
    );
}

vector sincScaled::origin() const
{
    return origin_; 
}

scalar sincScaled::amplitude() const
{
    return amplitude_; 
}

scalar sincScaled::omega() const
{
    return omega_; 
}

} // End namespace Foam::TriSurfaceImmersion

