/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::geomSurfaceCellMeshIntersection

Description


\*---------------------------------------------------------------------------*/

#ifndef geomSurfaceCellMeshIntersectionI_H
#define geomSurfaceCellMeshIntersectionI_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Time& geomSurfaceCellMeshIntersection::time() const
{
    return runTime_; 
}

const volScalarField& geomSurfaceCellMeshIntersection::sqrSearchDist() const
{
    return sqrSearchDist_;  
}

const volScalarField& geomSurfaceCellMeshIntersection::signedDist() const
{
    return signedDist_;  
}

const volScalarField& geomSurfaceCellMeshIntersection::signedDist0() const
{
    return signedDist0_;  
}

const surfaceScalarField& geomSurfaceCellMeshIntersection::lambda() const
{
    return lambda_;  
}

const triSurface& geomSurfaceCellMeshIntersection::surface() const
{
    return triSurfPtr_(); 
}

triSurface& geomSurfaceCellMeshIntersection::surfaceRef() 
{
    return triSurfPtr_(); 
}

label geomSurfaceCellMeshIntersection::Nx() const
{
    return Nx_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
