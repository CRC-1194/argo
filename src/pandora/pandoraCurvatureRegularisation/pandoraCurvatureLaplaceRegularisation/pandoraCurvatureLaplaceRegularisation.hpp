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
    Foam::pandoraCurvatureLaplaceRegularisation

Description
    Laplace smoothing for the curvature. 

SourceFiles
    pandoraCurvatureLaplaceRegularisation.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraCurvatureLaplaceRegularisation_H
#define pandoraCurvatureLaplaceRegularisation_H

#include "pandoraCurvatureRegularisation.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraCurvatureLaplaceRegularisation Declaration
\*---------------------------------------------------------------------------*/

class pandoraCurvatureLaplaceRegularisation
    :
        public pandoraCurvatureRegularisation
{

public:

    // Static Data Members

    TypeName ("Laplace");

    // Constructors

    //- Construct from components
    pandoraCurvatureLaplaceRegularisation(const dictionary& dict);

    //- Construct as copy
    pandoraCurvatureLaplaceRegularisation(const pandoraCurvatureLaplaceRegularisation&) = default;

    //- Destructor
    virtual ~pandoraCurvatureLaplaceRegularisation() = default;

    // Member Functions

    virtual void regularise(volScalarField& curvature, const volScalarField& isInterfaceCell); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
