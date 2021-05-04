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
    Foam::TriSurfaceImmersion::volumeFractionCalculator

Description
    Interface class enabling runtime type selection for volume fraction
    calculators, namely surface-mesh-cell-intersection and
    surface-mesh-cell-approximation.

SourceFiles
    volumeFractionCalculator.cpp
    volumeFractionCalculatorI.hpp

\*---------------------------------------------------------------------------*/

#ifndef volumeFractionCalculator_H
#define volumeFractionCalculator_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "pointIndexHit.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"
#include "Time.H"
#include "triSurface.H"
#include "typeInfo.H"
#include "volFields.H"

#include "signedDistanceCalculator.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {

/*---------------------------------------------------------------------------*\
                 Class volumeFractionCalculator Declaration
\*---------------------------------------------------------------------------*/

class volumeFractionCalculator
{
private:

    // Private Data
    //- A reference to the mesh.  
    const fvMesh& mesh_; 

    //- A reference to time.
    const Time& runTime_;  

    //- Point mesh required to construct point fields
    pointMesh pMesh_;

    //- Write additional geometric objects used in the computation of
    //  volume fraction of interface cells
    const bool writeGeometry_; 


public:

    // Static Data Members
    TypeName("volumeFractionCalculatorInterface");

    declareRunTimeSelectionTable (
        autoPtr,
        volumeFractionCalculator,
        Dictionary,
        (
            const dictionary& configDict,
            const fvMesh& mesh
        ),
        (configDict, mesh)
    )


    // Constructors
    explicit volumeFractionCalculator
    (
        const dictionary& configDict,
        const fvMesh& mesh
    );


    // Selectors
    static autoPtr<volumeFractionCalculator>
    New
    (
        const dictionary& configDict,
        const fvMesh& mesh
    );


    // Member Functions
    //- Access

    inline const Time& time() const;

    inline const fvMesh& mesh() const;

    inline bool writeGeometry() const;

    virtual double nTrianglesPerCell() const = 0;

    virtual label nIntersectedCells() const = 0;

    virtual label maxRefinementLevel() const = 0;

    virtual const signedDistanceCalculator& sigDistCalc() const = 0;


    //- Computation

    virtual void findIntersectedCells() = 0;

    static void bulkVolumeFraction(volScalarField& alpha, const volScalarField& inOut);

    virtual void calcVolumeFraction(volScalarField& alpha) = 0;


    //- Write 

    virtual void writeFields() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volumeFractionCalculatorI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //