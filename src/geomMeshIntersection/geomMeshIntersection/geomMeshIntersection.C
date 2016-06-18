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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "geomMeshIntersection.H"
#include "argList.H"
#include <omp.h>

namespace Foam
{

namespace GeometricalTransport
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geomMeshIntersection::geomMeshIntersection(const fvMesh& baseMesh, const fvMesh& toolMesh)
:
    baseMesh_(baseMesh),
    toolMesh_(toolMesh), 
    baseAABBs_(baseMesh.nCells()), 
    toolAABBs_(toolMesh.nCells()), 
    AABBintersects_(baseMesh.nCells()),
    cellPolyhedra_(baseMesh.nCells()), 
    volFraction_(
        IOobject
        (
            "alphaMesh", 
            baseMesh_.time().timeName(), 
            baseMesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ), 
        baseMesh_, 
        dimensionedScalar
        (
            "zero", 
            dimless, 
            pTraits<scalar>::zero
        )
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& geomMeshIntersection::volFraction() 
{

    baseMesh_.time().cpuTimeIncrement();

    const cellList& bCells = baseMesh_.cells(); 
    const pointField& bPoints = baseMesh_.points(); 
    const faceList& bFaces = baseMesh_.faces();
    const auto& Vb = baseMesh_.V();

    const cellList& tCells = toolMesh_.cells(); 
    const pointField& tPoints = toolMesh_.points(); 
    const faceList& tFaces = toolMesh_.faces();

    Info << "Setting the volume fraction by mesh intersection..." << endl;
    double time = omp_get_wtime(); 
    #pragma omp parallel 
    {
        // Computing bounding boxes for the base mesh.
        #pragma omp for schedule(dynamic) nowait
        for(label i = 0; i < baseMesh_.nCells(); ++i) 
            baseAABBs_[i] = boundBox(bPoints, bCells[i].labels(bFaces));

        // Computing bounding boxes for the tool mesh.
        #pragma omp for schedule(dynamic) 
        for(label i = 0; i < toolMesh_.nCells(); ++i) 
            toolAABBs_[i] = boundBox(tPoints, tCells[i].labels(tFaces));

        // Computing bounding boxes intersections. 
        #pragma omp for schedule(dynamic)
        for (decltype (baseAABBs_.size()) i = 0; i < baseAABBs_.size(); ++i)
            for (decltype(toolAABBs_.size()) j = 0; j < toolAABBs_.size(); ++j)
                if (baseAABBs_[i].overlaps(toolAABBs_[j]))
                    AABBintersects_[i].push_back(j);

        // Intersecting cells whose AABBs intersect
        #pragma omp for schedule(dynamic)
        for(decltype(AABBintersects_.size()) i = 0; i < AABBintersects_.size(); ++i)
        {
            auto baseCellHspaces = build<halfspaceVector>(i, baseMesh_);
            polyhedronSeq results;  
            for(const auto j : AABBintersects_[i])
            {
                auto inputCellPolyhedron = build<polyhedron>(j, toolMesh_);
                auto result = 
                    intersect<polyhedronIntersection>(baseCellHspaces, inputCellPolyhedron);

                if (result.polyhedron_size() > 3)
                    results.push_back(result.polyhedron()); 
            }
            cellPolyhedra_[i] = results; 
        }

        // Setting the volume fraction : volume(intersected polyhedron) / cellVolume 
        // Done for the base mesh. 
        #pragma omp for schedule(dynamic)
        for(decltype(cellPolyhedra_.size()) i = 0;  i < cellPolyhedra_.size(); ++i)
        {
            if (cellPolyhedra_[i].size() > 0)
            {
                for (const auto& poly : cellPolyhedra_[i])
                    volFraction_[i] += volume(poly) / Vb[i];
            }
        }
    }
    Info << "Done in " << omp_get_wtime() - time << " seconds." << endl;

    return volFraction_;  
}


Ostream& geomMeshIntersection::report (Ostream& os) const
{
    // Report volume fraction field errors.
    os << "Tool mesh volume: " << sum(toolMesh_.V()).value() << endl; 

    // Compute the tool mesh volume from the volume fraction provided by 
    // the intersection of two meshes..
    scalar toolMeshVolume = sum(baseMesh_.V() * volFraction_).value();

    os << "Tool mesh volume by volume fraction: " 
        << toolMeshVolume << endl;

    // Total volume of the tool mesh (cell magnitudes).
    scalar toolVolSum = sum(toolMesh_.V()).value();  

    // Total volume of the tool mesh (from the vol. fraction field on the
    // base mesh).
    scalar toolVolFromVolFracSum = sum(baseMesh_.V() * volFraction_).value(); 

    // Relative difference of the tool mesh volumes (cells, vol. fraction).
    scalar volFracRelError = mag(toolVolSum - toolVolFromVolFracSum) / toolVolSum; 

    os << "Relative error in volume fraction: " << volFracRelError << endl;

    // Report numerical boundedness error for the volume fraction.
    //auto epsilonB = max(max(mag(min(0.0,volFractionPtr_()))));
    auto boundednessError = max(
        max(mag(min(0.0,volFraction_))),
        max(mag(min(0.0,1.0-volFraction_)))
    ).value();
    os << "Maximal boundedness error: " << boundednessError << "\nDone.\n";

    // Uncomment to debug.
    //Info << "Writing the intersection polyhedra to VTK. " << endl;
    //const auto& runTime = baseMesh_.time(); 
    //vtk_polydata_stream polyhedronStream(runTime.path() + "/cutPolyhedra.vtk");

    //for (const auto& polys : cellPolyhedra_)
        //for (const auto& poly: polys)
            //polyhedronStream << poly;
    //Info << "Done." << endl;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GeometricalTransport 

// Global utility function for mesh intersection.
void setMeshIntersectionArgs(int argc, char* argv[])
{
    // Help information for applications that use geomIntersectMeshes.
    argList::addNote
    (
        "Intersect the base mesh with the tool mesh and store the volume fraction given by the intersection on the base mesh."
    );

    // No MPI parallel implementation is available at this point. 
    argList::noParallel();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

