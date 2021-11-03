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
    Foam::pandoraCurvature

Description
    Curvature models. 

SourceFiles
    pandoraCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef pandoraCurvature_H
#define pandoraCurvature_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pandoraCurvature Declaration
\*---------------------------------------------------------------------------*/

class pandoraCurvature
{

protected:

    const dictionary& curvatureDict_;
    volScalarField cellCurvature_;

public:

    // Static Data Members

    TypeName ("pandoraCurvature");

    declareRunTimeSelectionTable
    (
        autoPtr,
        pandoraCurvature, 
        Dictionary, 
        (
            const fvMesh& mesh,
            const dictionary& dict 
        ), 
        (mesh, dict)
    )

    // Constructors

    //- Construct from components
    pandoraCurvature(const fvMesh& mesh, const dictionary& dict);

    //- Construct as copy
    pandoraCurvature(const pandoraCurvature&) = default;

    // Selectors
    static autoPtr<pandoraCurvature> New(
        const fvMesh& mesh,
        const dictionary& dict 
    );

    //- Destructor
    virtual ~pandoraCurvature() = default;

    // Member Functions
    virtual volScalarField& cellCurvature(); 

    const fvMesh& mesh() const
    {
        return cellCurvature_.mesh();
    }

    const dictionary& curvatureDict() const
    {
        return curvatureDict_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
