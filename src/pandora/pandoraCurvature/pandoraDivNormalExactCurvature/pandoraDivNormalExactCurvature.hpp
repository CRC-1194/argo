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
    Foam::pandoraDivNormalExactCurvature

Description
    Extend the interface normal vector field by diffusion and calculate the 
    curvature as div(normal). Re-sets the normals in interface cells during 
    smoothing. 

SourceFiles
    pandoraDivNormalExactCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraDivNormalExactCurvature_H
#define pandoraDivNormalExactCurvature_H

#include "pandoraCurvature.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraDivNormalExactCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraDivNormalExactCurvature
    :
        public pandoraCurvature
{

protected:

    const word fieldName_; 
    const label nPropagate_;
    const label nAverage_; 
    volVectorField averagedNormals_; 

public:

    // Static Data Members

    TypeName ("divNormalExact");

    // Constructors

    //- Construct from components
    pandoraDivNormalExactCurvature(const fvMesh& mesh, const dictionary& dict);

    //- Construct as copy
    pandoraDivNormalExactCurvature(const pandoraDivNormalExactCurvature&) = default;

    //- Destructor
    virtual ~pandoraDivNormalExactCurvature() = default;

    // Member Functions
    
    virtual volScalarField& cellCurvature(); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
