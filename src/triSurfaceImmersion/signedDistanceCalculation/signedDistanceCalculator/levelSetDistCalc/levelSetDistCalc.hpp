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
    Foam::TriSurfaceImmersion::levelSetDistCalc

Description
    Compute the signed distance to a surface implicitly given by a level
    set.

SourceFiles
    levelSetDistCalcI.hpp
    levelSetDistCalc.cpp

\*---------------------------------------------------------------------------*/

#ifndef levelSetDistCalc_H
#define levelSetDistCalc_H

#include "pointFieldsFwd.H"
#include "signedDistanceCalculator.hpp"

#include "implicitSurfaces.hpp"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion {


/*---------------------------------------------------------------------------*\
                         Class levelSetDistCalc Declaration
\*---------------------------------------------------------------------------*/

class levelSetDistCalc
:
    public signedDistanceCalculator
{
    // Private Data
    autoPtr<implicitSurface> surfacePtr_;

    label maxIt_;
    scalar epsilon_;

    volScalarField cellLevelSetValues_;
    pointScalarField pointLevelSetValues_;
    std::vector<label> narrowBandCells_;


    // Private Member Functions
    point surfacePoint(const point& p) const;
    point closestPoint(const point& p) const;

    void computeMeshLevelSetValues();
    void identifyNarrowBandCells();
    void computeSignedDistances();
    void setInsideOutside();

public:

    TypeName("levelSet");


    // Constructors

    levelSetDistCalc(const dictionary& configDict, const fvMesh& mesh);


    // Member Functions
    scalar signedDistance(const point& x) const override;
    scalar referenceLength() const override;
    label nSurfaceElements() const override;
    scalar surfaceEnclosedVolume() const override;

    // Access
    inline const implicitSurface& surface() const;
    inline label maxIter() const;
    inline scalar epsilon() const;
    inline const volScalarField& cellLevelSetValues() const;
    inline const pointScalarField& pointLevelSetValues() const;

    // Write
    void writeFields() const override;


    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "levelSetDistCalcI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //