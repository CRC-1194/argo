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

#include "fvMesh.H"
#include "volMesh.H"
#include "Time.H"
#include "surfaceFields.H"
#include "typeInfo.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "autoPtr.H" 
//#include "Geometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace GeometricalTransport {

/*---------------------------------------------------------------------------*\
             Class geomSurfaceCellMeshIntersection Declaration
\*---------------------------------------------------------------------------*/

class geomSurfaceCellMeshIntersection
{
    // Private data
    
        //- A reference to the mesh.  
        const fvMesh& mesh_; 

        //- A reference to time.
        const Time& runTime_;  

        //- Squared search distance field used for the sign propagation.
        volScalarField sqrSearchDist_;  

        //- The signed distance field used for bulk cell identification. 
        volScalarField signedDist_; 

        // Initial signed distance field given by the octree, used for 
        // marker field calculation near boundaries where the Laplace 
        // equation solution may become unstable.
        volScalarField signedDist0_;  

        // Laplace sign equation coefficient, defaults to 1.  
        //surfaceScalarField lambda_; 

        // Triangular surface mesh with the octree used for spatial queries. 
        // AutoPtr used because of constructor delegation. 
        // Copy construction and assignment perform deep copies.
        // tmp<triSurface> would be better, but triSurface does not inherit from refCount.
        autoPtr<triSurface> triSurfPtr_;

        // Information used to store the surface proximity information for each cell. 
        DynamicList<pointIndexHit> cellNearestTriangle_;

        // Factor used to extend the narrow band by N cells. 
        // If sqrDistanceFactor = 2, the narrow band is extended by 2 cells. 
        const scalar sqrDistFactor_; 

        // Number of intersected cells. 
        label Nx_;
        // Average number of triangles (intersections) per intersected cell. 

public:
    // ADDGEOM
    //typedef TriangulationIntersection<tetrahedronVector,pointVectorVector> triangulationIntersection;

    TypeName ("surfaceCellMeshIntersection"); 

    // Constructors
    
        //- Construct from a mesh (initialize field data)
        geomSurfaceCellMeshIntersection
        (
            const fvMesh& mesh, 
            // Per default, don't write fields down. Writing is used for testing only.
            const scalar sqrDistFactor = 3, // Narrow band minimal width = 3 cells. 
            const IOobject::writeOption& wo = IOobject::NO_WRITE
        ); 

        //- Construct from a mesh and the file name for the surface file. 
        geomSurfaceCellMeshIntersection
        (
            const fvMesh& mesh, 
            const fileName& triName, 
            const scalar sqrDistFactor = 3, // Narrow band minimal width = 3 cells.
            const IOobject::writeOption& wo = IOobject::NO_WRITE 
        );

        //- Construct from a mesh and a triSurface 
        geomSurfaceCellMeshIntersection(
            const fvMesh& mesh, 
            const triSurface& triSurf,
            const scalar sqrDistFactor = 3, // Narrow band minimal width = 3 cells.
            const IOobject::writeOption& wo = IOobject::NO_WRITE
        );

        //- Construct copy
        geomSurfaceCellMeshIntersection(const geomSurfaceCellMeshIntersection& copy); 


    //- Destructor
        virtual ~geomSurfaceCellMeshIntersection() = default;

    // Member Functions

        //- Access

        const Time& time() const;

        const volScalarField& sqrSearchDist() const; 

        const volScalarField& signedDist() const; 

        const volScalarField& signedDist0() const; 

        //const surfaceScalarField& lambda() const; 

        const triSurface& surface() const;

        triSurface& surfaceRef(); 

        label Nx() const; 

        //- Computation

        // Squared search distance calculation function. 
        virtual void calcSqrSearchDist();  

        // Signed distance computation using triSurfaceSearch and an external triSurface.
        virtual void calcSignedDist(const triSurface& tri, const triSurfaceSearch& triSearch);  

        // Signed distance computation using triSurfaceSearch and the internal triSurface.
        virtual void calcSignedDist();  

        // Volume fraction calculation.
        virtual void calcVolFraction(volScalarField& alpha, const triSurface& tri);

        virtual void calcVolFraction(volScalarField& alpha);

        //- Write 
        virtual void writeFields() const;

    // Operators 
        void operator=(const geomSurfaceCellMeshIntersection& rhs);
};

#include "geomSurfaceCellMeshIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
