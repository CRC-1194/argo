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
    Foam::curvatureFunctionObject

Description

SourceFiles
    curvatureFunctionObjectI.H
    curvatureFunctionObject.C
    curvatureFunctionObjectIO.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureFunctionObject_H
#define curvatureFunctionObject_H

#include "functionObject.H"
#include "messageStream.H"
#include "typeInfo.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {


/*---------------------------------------------------------------------------*\
                         Class curvatureFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class curvatureFunctionObject
:
    public functionObject
{
    // Private Data
    const Time& time_;
    word vofFieldName_;
    word curvatureFieldName_;
    scalar exactCurvature_;
    bool firstWrite_ = true;

    // Curvature metrics
    std::vector<scalar> timeSteps_;
    std::vector<scalar> meanCurvature_;
    std::vector<scalar> standardDeviation_;
    std::vector<scalar> minCurvature_;
    std::vector<scalar> maxCurvature_;
    std::vector<scalar> errorL1_;
    std::vector<scalar> errorL2_;
    std::vector<scalar> errorLinf_;

    // Private Member Functions
    const volScalarField& getField(const word& fieldName) const;
    fileName assembleFileName() const;
    std::vector<scalar> filterCurvature() const;
    void clearDataBuffer();

public:

    TypeName ("curvatureEvaluation");

    // Constructors
    curvatureFunctionObject(const word& name, const Time& time, const dictionary& dict);


    //- Destructor
    virtual ~curvatureFunctionObject() = default;


    // Member Functions
    virtual bool execute() final;
    virtual bool write() final;
    virtual bool read(const dictionary& dict) final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "curvatureFunctionObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
