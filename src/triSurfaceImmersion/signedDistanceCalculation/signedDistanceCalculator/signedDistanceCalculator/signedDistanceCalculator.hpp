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
    Foam::TriSurfaceImmersion::signedDistanceCalculator

Description
    Unified interface class for signed distance calculation based on
    triangulated surfaces and implicit level set functions.

SourceFiles
    signedDistanceCalculatorI.H
    signedDistanceCalculator.C
    signedDistanceCalculatorIO.C

\*---------------------------------------------------------------------------*/

#ifndef signedDistanceCalculator_H
#define signedDistanceCalculator_H

#include "DynamicList.H"
#include "pointFields.H"
#include "pointIndexHit.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                         Class signedDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class signedDistanceCalculator
{
private:
    // Private Data
    const dictionary& dict_;
    const fvMesh& mesh_;
    pointMesh pMesh_;
    scalar narrowBandWidth_;
    scalar outOfNarrowBandValue_;


    // Private Member Functions

protected:

    // Create shared fields here?
    //- Signed distance at cell centers. 
    volScalarField cellSignedDist_; 
    //- Initial signed distance field given by the octree, used to correct the 
    //  signed distance propagated by the solution of the Laplace equation. 
    volScalarField cellSignedDist0_;  
    //- Nearest surface triangle to a cell centre in the narrow band
    DynamicList<pointIndexHit> cellNearestTriangle_;

    //- Signed distance at cell corner points. 
    pointScalarField pointSignedDist_;
    //- Nearest surface triangle to a cell centre in the narrow band
    DynamicList<pointIndexHit> pointNearestTriangle_;


public:

    // Static Data Members
    TypeName("signedDistanceCalculatorInterface");

    declareRunTimeSelectionTable (
        autoPtr,
        signedDistanceCalculator,
        Dictionary,
        (
            const dictionary& configDict,
            const fvMesh& mesh
        ),
        (configDict, mesh)
    )


    // Constructors

        //- Construct from components
        explicit signedDistanceCalculator(const dictionary& configDict, const fvMesh& mesh);


    // Selectors

        //- Select default constructed
        static autoPtr<signedDistanceCalculator> New(const dictionary& configDict, const fvMesh& mesh);


    // Member Functions

    // Access
    inline const dictionary& configDict() const;
    inline const fvMesh& mesh() const;
    inline scalar narrowBandWidth() const;
    inline scalar outOfNarrowBandValue() const;
    inline const DynamicList<pointIndexHit>& cellClosestPoint() const;
    inline const DynamicList<pointIndexHit>& pointClosestPoint() const;
    inline const volScalarField& cellSignedDist() const;
    inline const volScalarField& cellSignedDist0() const;
    inline const pointScalarField& pointSignedDist() const;

    // Check
    virtual scalar signedDistance(const point& x) const = 0;
    virtual scalar referenceLength() const = 0;
    virtual label nSurfaceElements() const = 0;
    virtual scalar surfaceEnclosedVolume() const = 0;

    // Edit
    void outOfNarrowBandValue(scalar value);
    void narrowBandWidth(scalar width);

    // Write
    virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "signedDistanceCalculatorI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
