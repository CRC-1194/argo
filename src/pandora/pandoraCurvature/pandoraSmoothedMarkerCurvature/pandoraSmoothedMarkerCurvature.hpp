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
    Foam::pandoraSmoothedMarkerCurvature

Description
    Extend the interface normal vector field by diffusion and calculate the 
    curvature as div(normal). Re-sets the normals in interface cells during 
    smoothing. 

SourceFiles
    pandoraSmoothedMarkerCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraSmoothedMarkerCurvature_H
#define pandoraSmoothedMarkerCurvature_H

#include "pandoraCurvature.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraSmoothedMarkerCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraSmoothedMarkerCurvature
    :
        public pandoraCurvature
{

protected:

    const word   markerFieldName_; 
    const scalar markerTolerance_; 
    const label  nMarkerAverages_;
    const label  nNormalAverages_;

public:

    // Static Data Members

    TypeName ("smoothedMarker");

    // Constructors

    //- Construct from components
    pandoraSmoothedMarkerCurvature(const fvMesh& mesh, const dictionary& dict);

    //- Construct as copy
    pandoraSmoothedMarkerCurvature(const pandoraSmoothedMarkerCurvature&) = default;

    //- Destructor
    virtual ~pandoraSmoothedMarkerCurvature() = default;

    // Member Functions
    
    virtual volScalarField& cellCurvature(); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
