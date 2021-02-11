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

const volScalarField& geomSurfaceCellMeshIntersection::cellSqrSearchDist() const
{
    return cellSqrSearchDist_;  
}

const volScalarField& geomSurfaceCellMeshIntersection::cellSignedDist() const
{
    return cellSignedDist_;  
}

const volScalarField& geomSurfaceCellMeshIntersection::cellSignedDist0() const
{
    return cellSignedDist0_;  
}

const triSurface& geomSurfaceCellMeshIntersection::surface() const
{
    return triSurf_; 
}

const double geomSurfaceCellMeshIntersection::nTrianglesPerCell() const
{
    return nTrianglesPerCell_; 
}

const label geomSurfaceCellMeshIntersection::nIntersectedCells() const
{
    return intersectedCellLabels_.size();  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //