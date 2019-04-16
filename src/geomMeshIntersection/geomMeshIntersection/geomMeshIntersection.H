/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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

        // FIXME: Too long running times for large cases. Split into AABBintersects_
        // and AABBinsides_. Intersect the intersects_ and sum the volume from the 
        // insides. TM. 
        //- AABB of baseMesh_ intersects AABB of toolMesh_
        labelSeqSeq AABBintersects_;

        //- List of polyhedron sets 
        polyhedronSeqSeq cellPolyhedra_;

        //- Enable/disable VTK IO for the intersection result. 
        Switch writeGeometry_; 

        // Number of intersected cells.
        label Nx_; 

        // Average number of intersections / intersected cell. 
public:

    // Constructors

        //- Construct from components
        geomMeshIntersection
        (
            const fvMesh& baseMesh, 
            const fvMesh& toolMesh, 
            const Switch& writeGeo=false
        );

    //- Destructor
        ~geomMeshIntersection() = default;

    // Member Functions

        void setVolFraction(volScalarField& volFraction); 

        // Access

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

        label Nx() const 
        {
            return Nx_; 
        }

        // Write
        Ostream& report (Ostream& os, const volScalarField& volFraction) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::GeometricalTransport 

#endif

// ************************************************************************* //
