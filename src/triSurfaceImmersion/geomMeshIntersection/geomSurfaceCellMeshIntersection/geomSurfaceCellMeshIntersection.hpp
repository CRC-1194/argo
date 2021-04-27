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

SourceFiles
    geomSurfaceCellMeshIntersectionI.H
    geomSurfaceCellMeshIntersection.C
    geomSurfaceCellMeshIntersectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef geomSurfaceCellMeshIntersection_H
#define geomSurfaceCellMeshIntersection_H

// OpenFOAM 
#include "fvMesh.H"
#include "volMesh.H"
#include "Time.H"
#include "surfaceFields.H"
#include "typeInfo.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "autoPtr.H" 
#include "pointFields.H"
#include <memory>
#include "volPointInterpolation.H"
#include "DynamicList.H"

#include "Geophase.hpp"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace GeometricalTransport {

/*---------------------------------------------------------------------------*\
             Class geomSurfaceCellMeshIntersection Declaration
\*---------------------------------------------------------------------------*/

class geomSurfaceCellMeshIntersection
{
    //- A reference to the mesh.  
    const fvMesh& mesh_; 

    //- A reference to time.
    const Time& runTime_;  

    // Cell-center distances
    
    //- Squared search distance field in cell centers. 
    volScalarField cellSqrSearchDist_;  
    //- Signed distance at cell centers. 
    volScalarField cellSignedDist_; 
    //- Initial signed distance field given by the octree, used to correct the 
    // signed distance propagated by the solution of the Laplace equation. 
    volScalarField cellSignedDist0_;  
    //- Information used to store the surface proximity information for each cell. 
    DynamicList<pointIndexHit> cellNearestTriangle_;

    // Face-center distances 
    //- Signed distance at face centers. 
    surfaceScalarField faceSignedDist_; 

    // Point signed distances 
    
    //- Inverse Distance Interpolation : cell centers to cell corners. 
    volPointInterpolation cellsToPointsInterp_;

    pointMesh pMesh_;
    //- Squared search distance at cell corner points. 
    pointScalarField pointSqrSearchDist_;  
    //- Signed distance at cell corner points. 
    pointScalarField pointSignedDist_;
    //- Nearest triangle to a cell corner point. 
    DynamicList<pointIndexHit> pointNearestTriangle_;

    // Triangulated surface.
    
    // Surface mesh
    const triSurface& triSurf_;
    // Octree of the triangulated surface.
    triSurfaceSearch triSurfSearch_;
    
    // Indexes of the intersected cells.
    using dynamicLabelList = DynamicList<label>; 
    dynamicLabelList intersectedCellLabels_;

    // Factor used to extend the narrow band by N cells. 
    // If sqrDistanceFactor = 2, the narrow band is extended by 2 cells. 
    const scalar sqrDistFactor_; 

    const bool writeGeo_; 

    // Average number of triangles per cell, used for complexity analysis.
    scalar nTrianglesPerCell_;  

public:

    TypeName ("surfaceCellMeshIntersection"); 

    // Constructors
    
    //- Construct from a mesh (initialize field data)
    geomSurfaceCellMeshIntersection
    (
        const fvMesh& mesh, 
        const triSurface& triSurf, 
        const scalar sqrDistFactor = 3, // Narrow band minimal width = 3 cells. 
        const bool writeGeo = false,
        const IOobject::writeOption& wo = IOobject::NO_WRITE
    ); 

    //- Destructor
    virtual ~geomSurfaceCellMeshIntersection() = default;

    // Member Functions

    //- Access

    const Time& time() const;

    const volScalarField& cellSqrSearchDist() const; 

    const pointScalarField& pointSqrSearchDist() const; 

    const volScalarField& cellSignedDist() const; 

    const volScalarField& cellSignedDist0() const; 

    const triSurface& surface() const;

    const double nTrianglesPerCell() const;

    const label nIntersectedCells() const;

    //- Computation

    // Squared search distance calculation. 
    virtual void calcSqrSearchDists();  

    // Signed distance computation using triSurfaceSearch and an external triSurface.
    virtual void calcSignedDist();  

    virtual void findIntersectedCells(); 

    // Volume fraction calculation.
    virtual void calcVolFraction(volScalarField& alpha);

    //- Write 
    virtual void writeFields() const;

};

#include "geomSurfaceCellMeshIntersectionI.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
