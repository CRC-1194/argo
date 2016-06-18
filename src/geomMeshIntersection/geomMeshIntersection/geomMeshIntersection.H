/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::geomMeshIntersection

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

Description
    Class for a geometrical intersection of two meshes. 
    The intersection test is O(nm) where n and m are the mesh sizes. 
    Better implementation requires the use of a AABB tree for intersection
    testing.

    Used for error estimation and pre-processing of the voFoam solver.

SourceFiles
    geomMeshIntersection.C

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "boundBox.H"
#include "volFields.H"
#include "Ostream.H"
#include "Geometry.H"
#include "geomTransportControl.H"

#include <list>
#include <memory>

#ifndef geomMeshIntersection_H
#define geomMeshIntersection_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GeometricalTransport
{

/*---------------------------------------------------------------------------*\
                         Class geomMeshIntersection Declaration
\*---------------------------------------------------------------------------*/

class geomMeshIntersection
{
    public:

        typedef std::vector<boundBox> boundBoxSeq;
        typedef std::vector<label> labelSeq;
        typedef std::vector<labelSeq> labelSeqSeq; 

        typedef pointVectorVector polyhedron; 
        typedef PolygonSequenceIntersection<polyhedron> polyhedronIntersection; 
        typedef std::vector<polyhedron> polyhedronSeq;
        typedef std::vector<polyhedronSeq> polyhedronSeqSeq;

    private:

    // Private data 
        //- References to intersected meshes 
        const fvMesh& baseMesh_;
        const fvMesh& toolMesh_;

        // Base mesh Axis Aligned Bound Boxes
        boundBoxSeq baseAABBs_;

        // Input mesh Axis Aligned Bound Boxes
        boundBoxSeq toolAABBs_;

        //- AABB of baseMesh_ intersects AABB of toolMesh_
        labelSeqSeq AABBintersects_;

        //- List of polyhedron sets 
        polyhedronSeqSeq cellPolyhedra_;

        //- Vol fraction field
        volScalarField volFraction_;

    // Private Member Functions 
    
        //- Compute the bounding boxes for both meshes.
        void computeBoundBoxes();

        //- Calculate intersections between bounding boxes O(m * n), 
        //  m = number of cells in baseMesh_, and n -//- in toolMesh_ 
        void intersectAABB();

        // Geometrically intersect polyhedra and store them.
        void intersectPolyhedra();

        // Set the volume fraction based on the intersected polyhedra.
        void setIntersectedVolFraction();

public:

    // Constructors

        //- Construct from components
        geomMeshIntersection(const fvMesh& baseMesh, const fvMesh& toolMesh);

    //- Destructor
        ~geomMeshIntersection() = default;

    // Member Functions

        // Access
        const volScalarField& volFraction(); 

        const polyhedronSeqSeq& cellPolyhedra() const
        {
            return cellPolyhedra_;
        };

        const Time& baseTime() const
        {
            return baseMesh_.time();
        }

        const fvMesh& baseMesh() const
        {
            return baseMesh_;
        }

        // Write
        Ostream& report (Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GeometricalTransport 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global utility functions for mesh intersection.

void getRootCase(fileName& casePath);

// Set the application arguments for mesh intersection.
void setMeshIntersectionArgs(int argc, char* argv[]);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
