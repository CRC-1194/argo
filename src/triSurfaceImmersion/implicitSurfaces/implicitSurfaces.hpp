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

#ifndef implicitSurfaces_H
#define implicitSurfaces_H

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "ITstream.H"
#include "vector.H"

namespace Foam::TriSurfaceImmersion {

    class implicitSurface
    {
        public: 

            TypeName("implicitSurface");

            declareRunTimeSelectionTable
            (
                autoPtr,
                implicitSurface, 
                ITstream, 
                (
                    ITstream& is
                ), 
                (is)
            )

            declareRunTimeSelectionTable
            (
                autoPtr,
                implicitSurface, 
                Dictionary, 
                (
                    const dictionary& configDict
                ), 
                (configDict)
            )

            static autoPtr<implicitSurface> New(
                const word& name, 
                ITstream& is
            );

            static autoPtr<implicitSurface> New(
                const dictionary& configDict
            );

            implicitSurface() = default;

            explicit implicitSurface(ITstream&) {};
            explicit implicitSurface(const dictionary&) {};

            virtual scalar value(const vector&) const = 0;
            virtual scalar operator()(const vector&) const = 0;
            virtual vector grad(const vector&) const = 0;
            virtual scalar volume() const;
    };

    class plane : public implicitSurface
    {
        vector position_; 
        vector normal_; 

        public:

            TypeName ("plane");

            plane(vector position, vector normal);

            explicit plane(ITstream& is);

            explicit plane(const dictionary& configDict);

            scalar value(const vector& x) const override;

            scalar operator()(const vector& x) const override;
            
            vector grad(const vector& x) const override;

            vector position() const;

            vector normal() const;
    };

    class sphere : public implicitSurface
    {
        vector center_; 
        scalar radius_; 

        public:

            TypeName ("sphere");

            sphere(vector center, scalar radius);

            explicit sphere(ITstream& is);

            explicit sphere(const dictionary& configDict);

            scalar value(const vector& x) const override;

            scalar operator()(const vector& x) const override;

            vector grad(const vector& x) const override;

            scalar volume() const override;

            vector center() const;

            scalar radius() const;
    };

    class ellipsoid : public implicitSurface
    {
        vector center_; 
        vector axes_; 
        vector axesSqr_;

        void setAxesSqr(const vector& axes);

        public:

            TypeName ("ellipsoid");

            ellipsoid(vector center, vector axes);

            explicit ellipsoid(ITstream& is);

            explicit ellipsoid(const dictionary& configDict);

            scalar value(const vector& x) const override;

            scalar operator()(const vector& x) const override;
            
            vector grad(const vector& x) const override;

            scalar volume() const override;

            vector center() const;

            vector axes() const;
    };

    class sinc : public implicitSurface
    {
        vector origin_; 
        scalar amplitude_; 
        scalar omega_;

        public:

            TypeName ("sinc");

            sinc(vector origin, scalar amplitude, scalar omega);

            explicit sinc(ITstream& is);

            explicit sinc(const dictionary& configDict);

            scalar value(const vector& x) const override;

            scalar operator()(const vector& x) const override;
            
            vector grad(const vector& x) const override;

            vector origin() const;

            scalar amplitude() const;

            scalar omega() const;
    };

    class sincScaled : public implicitSurface // TODO (TM): Scale the amplitude 
    {
        vector origin_; 
        scalar amplitude_; 
        scalar omega_;

        public:

            TypeName ("sincScaled");

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

} // End namespace Foam::TriSurfaceImmersion

#endif