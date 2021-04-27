/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tomislav Maric, Tobias Tolle, TU Darmstadt
                       Anja Lippert, BOSCH CR 
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
    Foam::pandoraCurvatureExtension

Description
    Curvature models. 

SourceFiles
    pandoraCurvatureExtension.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraCurvatureExtension_H
#define pandoraCurvatureExtension_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraCurvatureExtension Declaration
\*---------------------------------------------------------------------------*/

class pandoraCurvatureExtension
{
    // Private Data
    const dictionary& extensionDict_;
    
    // Private Member Functions

public:

    // Static Data Members

    TypeName ("pandoraCurvatureExtension");

    declareRunTimeSelectionTable
    (
        autoPtr,
        pandoraCurvatureExtension, 
        Dictionary, 
        (
            const dictionary& dict 
        ), 
        (dict)
    );

    // Constructors

    //- Construct from components
    pandoraCurvatureExtension(const dictionary& dict);

    //- Construct as copy
    pandoraCurvatureExtension(const pandoraCurvatureExtension&) = default;

    // Selectors
    static autoPtr<pandoraCurvatureExtension> New(
        const dictionary& dict 
    );

    //- Destructor
    virtual ~pandoraCurvatureExtension() = default;

    // Member Functions
    virtual void extend(volScalarField& curvature, const boolList& isInterfaceCell) = 0;
    const dictionary& extensionDict() const
    {
        return extensionDict_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
