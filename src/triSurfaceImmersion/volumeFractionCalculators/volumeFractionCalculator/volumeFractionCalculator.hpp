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
#include "fvMesh.H"
#include "pointFields.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"
#include "Time.H"
#include "triSurface.H"
#include "typeInfo.H"
#include "volFields.H"

#include "searchDistanceCalculator.hpp"
#include "signedDistanceCalculator.hpp"
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
    //- A reference to the mesh.  
    const fvMesh& mesh_; 

    //- A reference to time.
    const Time& runTime_;  

    //- A reference to the surface mesh
    const triSurface& surface_;

    //- Point mesh required to construct point fields
    pointMesh pMesh_;

    // Signed distances
    searchDistanceCalculator searchDistCalc_;
    signedDistanceCalculator sigDistCalc_;
    
    //- Signed distance at cell centers. 
    volScalarField cellSignedDist_; 
    //- Initial signed distance field given by the octree, used to correct the 
    //  signed distance propagated by the solution of the Laplace equation. 
    volScalarField cellSignedDist0_;  

    //- Signed distance at cell corner points. 
    pointScalarField pointSignedDist_;

    //- Write additional geometric objects used in the computation of
    //  volume fraction of interface cells
    const bool writeGeometry_; 


protected:


public:

    // Static Data Members
    TypeName("volumeFractionCalculatorInterface");

    declareRunTimeSelectionTable (
        autoPtr,
        volumeFractionCalculator,
        Dictionary,
        (
            const dictionary& configDict,
            const fvMesh& mesh,
            const triSurface& surface
        ),
        (configDict, mesh, surface)
    )


    // Constructors
    explicit volumeFractionCalculator
    (
        const dictionary& configDict,
        const fvMesh& mesh,
        const triSurface& surface
    );


    // Selectors
    static autoPtr<volumeFractionCalculator>
    New
    (
        const dictionary& configDict,
        const fvMesh& mesh,
        const triSurface& surface
    );


    // Member Functions
    //- Access

    inline const Time& time() const;

    inline const fvMesh& mesh() const;

    inline const triSurface& surface() const;

    inline const volScalarField& cellSignedDist() const; 

    inline const volScalarField& cellSignedDist0() const; 

    inline const pointScalarField& pointSignedDist() const;

    inline const signedDistanceCalculator& signedDistCalc() const;

    inline bool writeGeometry() const;

    virtual const double nTrianglesPerCell() const = 0;

    virtual const label nIntersectedCells() const = 0;

    virtual const label maxRefinementLevel() const = 0;


    //- Computation

    void calcSignedDist();  

    virtual void findIntersectedCells() = 0;

    void bulkVolumeFraction(volScalarField& alpha);

    virtual void calcVolumeFraction(volScalarField& alpha) = 0;


    //- Write 

    virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volumeFractionCalculatorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
