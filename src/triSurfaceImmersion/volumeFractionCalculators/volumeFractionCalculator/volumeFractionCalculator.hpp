/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::volumeFractionCalculator

Description
    Interface class enabling runtime type selection for volume fraction
    calculators, namely surface-mesh-cell-intersection and
    surface-mesh-cell-approximation.

SourceFiles
    volumeFractionCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef volumeFractionCalculator_H
#define volumeFractionCalculator_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {

/*---------------------------------------------------------------------------*\
                 Class volumeFractionCalculator Declaration
\*---------------------------------------------------------------------------*/

class volumeFractionCalculator
{
private:
    // Private Data

protected:


public:

    TypeName("volumeFractionCalculatorInterface");

    declareRunTimeSelectionTable (
        autoPtr,
        volumeFractionCalculator,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )

    // Static Data Members


    // Generated Methods

//        //- Default construct
//        volumeFractionCalculator() = default;
//
//        //- Copy construct
//        volumeFractionCalculator(const volumeFractionCalculator&) = default;
//
//        //- Copy assignment
//        volumeFractionCalculator& operator=(const volumeFractionCalculator&) = default;


    // Constructors
    explicit volumeFractionCalculator(const dictionary& configDict) {};


    // Selectors
    static autoPtr<volumeFractionCalculator> New(const dictionary& configDict);


    //- Destructor
    virtual ~volumeFractionCalculator() = default;


    // Member Functions
    void printTypeName() const;

    virtual void calcVolumeFraction(volScalarField& ) = 0;

    // Access

    // Check

    // Edit

    // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "volumeFractionCalculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
