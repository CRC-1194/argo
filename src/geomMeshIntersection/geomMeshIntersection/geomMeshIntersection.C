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

namespace Foam
{

namespace GeometricalTransport
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geomMeshIntersection::geomMeshIntersection
(
    const fvMesh& baseMesh, 
    const fvMesh& toolMesh, 
    const Switch& writeGeo 
)
:
    baseMesh_(baseMesh),
    toolMesh_(toolMesh), 
    baseAABBs_(baseMesh.nCells()), 
    toolAABBs_(toolMesh.nCells()), 
    AABBintersects_(baseMesh.nCells()),
    cellPolyhedra_(baseMesh.nCells()), 
    writeGeometry_(writeGeo),
    Nx_(0)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void geomMeshIntersection::setVolFraction(volScalarField& volFraction) 
{
    // Zero the measurement values. 
    Nx_ = 0; 

    baseMesh_.time().cpuTimeIncrement();

    const cellList& bCells = baseMesh_.cells(); 
    const pointField& bPoints = baseMesh_.points(); 
    const faceList& bFaces = baseMesh_.faces();
    const auto& Vb = baseMesh_.V();

    const cellList& tCells = toolMesh_.cells(); 
    const pointField& tPoints = toolMesh_.points(); 
    const faceList& tFaces = toolMesh_.faces();

    volFraction = dimensionedScalar(
            volFraction.name(), 
            volFraction.dimensions(), 
            scalar(0)
    );

    Info << "Setting the volume fraction by mesh intersection." << endl; 

    // Computing bounding boxes for the base mesh.
    for(label i = 0; i < baseMesh_.nCells(); ++i) 
        baseAABBs_[i] = boundBox(bPoints, bCells[i].labels(bFaces));

    // Computing bounding boxes for the tool mesh.
    for(label i = 0; i < toolMesh_.nCells(); ++i) 
        toolAABBs_[i] = boundBox(tPoints, tCells[i].labels(tFaces));

    // Computing bounding boxes intersections. 
    for (decltype (baseAABBs_.size()) i = 0; i < baseAABBs_.size(); ++i)
    {
        for (decltype(toolAABBs_.size()) j = 0; j < toolAABBs_.size(); ++j)
        {
            if (baseAABBs_[i].overlaps(toolAABBs_[j]))
                AABBintersects_[i].push_back(j);
        }
    }


    // Intersecting cells whose AABBs intersect
    // FIXME: Set volume directly, do not store the cell intersects, better efficiency. TM.
    for(decltype(AABBintersects_.size()) i = 0; i < AABBintersects_.size(); ++i)
    {
        auto baseCellHspaces = build<halfspaceVector>(i, baseMesh_);

        // TODO: Remove, debugging: 
        //auto baseCellPoly = build<polyhedron>(i, baseMesh_); 

        //for (unsigned int i = 0; i < baseCellPoly.size(); ++i)
        //{
            //if ((normal_area_vec(baseCellPoly[i]) & baseCellHspaces[i].direction()) > 0)
                //Info << "ERROR: hspace outward pointing " << baseMesh_.time().timeIndex() << " " << i << nl;
        //}

        polyhedronSeq results;  

        if (AABBintersects_[i].size())
        {
            autoPtr<vtk_polydata_stream> cellio; 
            autoPtr<vtk_polydata_stream> cutio; 
            autoPtr<vtk_polydata_stream> toolio; 

            if (writeGeometry_)
            {
                cellio = new vtk_polydata_stream(prependVtkFileName("cell", baseMesh_.time().timeIndex(), i));
                toolio = new vtk_polydata_stream(prependVtkFileName("toolCell", baseMesh_.time().timeIndex(), i));
                cutio = new vtk_polydata_stream(prependVtkFileName("toolCellCut", baseMesh_.time().timeIndex(), i));
                cellio->streamGeometry(build<pointVectorVector>(i, baseMesh_)); 
            }

            for(const auto j : AABBintersects_[i])
            {
                auto inputCellPolyhedron = build<polyhedron>(j, toolMesh_);
                auto result = 
                    intersect<polyhedronIntersection>(baseCellHspaces, inputCellPolyhedron);



                ++Nx_; 

                if (writeGeometry_) 
                {
                    toolio->streamGeometry(inputCellPolyhedron);  
                    cutio->streamGeometry(result.polyhedron());
                }

                if (result.polyhedron_size() > 3)
                    results.push_back(result.polyhedron()); 
            }
            cellPolyhedra_[i] = results; 
        }
    }

    // Setting the volume fraction : volume(intersected polyhedron) / cellVolume 
    // Done for the base mesh. 
    for(decltype(cellPolyhedra_.size()) i = 0;  i < cellPolyhedra_.size(); ++i)
    {
        if (cellPolyhedra_[i].size() > 0)
        {
            for (const auto& poly : cellPolyhedra_[i])
                volFraction[i] += volume(poly) / Vb[i];
        }
    }

    forAll(volFraction, cellI)
    {
        if (volFraction[cellI] > 1)
            volFraction[cellI] = 1.0; 
        else if (volFraction[cellI] > (1 - 5e-09)) 
            volFraction[cellI] = 1; 
    }

    // Correct BC field values after processing cell values.  
    volFraction.correctBoundaryConditions(); 
}


Ostream& geomMeshIntersection::report (Ostream& os, const volScalarField& volFraction) const
{
    // Report volume fraction field errors.
    os << "Tool mesh volume: " << sum(toolMesh_.V()).value() << endl; 

    // Compute the tool mesh volume from the volume fraction provided by 
    // the intersection of two meshes..
    scalar toolMeshVolume = sum(baseMesh_.V() * volFraction).value();

    os << "Tool mesh volume by volume fraction: " 
        << toolMeshVolume << endl;

    // Total volume of the tool mesh (cell magnitudes).
    scalar toolVolSum = sum(toolMesh_.V()).value();  

    // Total volume of the tool mesh (from the vol. fraction field on the
    // base mesh).
    scalar toolVolFromVolFracSum = sum(baseMesh_.V() * volFraction).value(); 

    // Relative difference of the tool mesh volumes (cells, vol. fraction).
    scalar volFracRelError = mag(toolVolSum - toolVolFromVolFracSum) / toolVolSum; 

    os << "Relative error in volume fraction: " << volFracRelError << endl;

    // Report numerical boundedness error for the volume fraction.
    auto boundednessError = max(
        max(mag(min(0.0,volFraction))),
        max(mag(min(0.0,1.0-volFraction)))
    ).value();
    os << "Maximal boundedness error: " << boundednessError << "\nDone.\n";

    if (writeGeometry_)
    {
        Info << "Writing the intersection polyhedra to VTK. " << endl;
        const auto& runTime = baseMesh_.time(); 
        vtk_polydata_stream polyhedronStream(
            prependVtkFileName(runTime.path() + "/cutPolyhedra", runTime.timeIndex())
        );

        for (const auto& polys : cellPolyhedra_)
        {
            if (polys.size())
            {
                for (const auto& poly: polys)
                    polyhedronStream << poly;
            }
        }
        Info << "Done." << endl;
    }

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam:GeometricalTransport 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

