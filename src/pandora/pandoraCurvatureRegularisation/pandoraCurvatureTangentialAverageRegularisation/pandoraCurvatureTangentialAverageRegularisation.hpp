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
    Foam::pandoraCurvatureTangentialAverageRegularisation

Description
    Laplace smoothing for the curvature. 

SourceFiles
    pandoraCurvatureTangentialAverageRegularisation.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraCurvatureTangentialAverageRegularisation_H
#define pandoraCurvatureTangentialAverageRegularisation_H

#include "pandoraCurvatureRegularisation.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraCurvatureTangentialAverageRegularisation Declaration
\*---------------------------------------------------------------------------*/

class pandoraCurvatureTangentialAverageRegularisation
    :
        public pandoraCurvatureRegularisation
{
    label nAveragingIterations_;

public:

    // Static Data Members

    TypeName ("tangentialAverage");

    // Constructors

    //- Construct from components
    pandoraCurvatureTangentialAverageRegularisation(const dictionary& dict);

    //- Construct as copy
    pandoraCurvatureTangentialAverageRegularisation(const pandoraCurvatureTangentialAverageRegularisation&) = default;

    //- Destructor
    virtual ~pandoraCurvatureTangentialAverageRegularisation() = default;

    // Member Functions
    virtual void regularise(volScalarField& curvature, const boolList& isInterfaceCell); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
