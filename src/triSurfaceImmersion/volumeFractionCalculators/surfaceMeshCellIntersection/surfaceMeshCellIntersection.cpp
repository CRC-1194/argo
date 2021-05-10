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

\*---------------------------------------------------------------------------*/

#include "surfaceMeshCellIntersection.hpp"

#include <bitset>
#include <cassert>
#include <cmath>
#include <limits>

#include "addToRunTimeSelectionTable.H"

#include "Geophase.hpp"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(surfaceMeshCellIntersection, 0);
addToRunTimeSelectionTable(
    volumeFractionCalculator, surfaceMeshCellIntersection, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void surfaceMeshCellIntersection::interfaceCellVolumeFraction(
    volScalarField& alpha)
{
    const auto& mesh = this->mesh();
    const auto& meshCells = mesh.cells();
    const auto& cellCenters = mesh.C();
    const auto& cellVolumes = mesh.V();

    const auto& faces = mesh.faces();
    const auto& faceCenters = mesh.Cf();

    const auto& meshPoints = mesh.points();

    const auto& cellSignedDist = sigDistCalc_.cellSignedDist();
    const auto& pointSignedDist = sigDistCalc_.pointSignedDist();

    const auto& octree = sigDistCalc_.surfaceSearch().tree();
    const auto& triSurf = sigDistCalc_.surface();
    const auto& triPoints = triSurf.points();
    const auto& triNormals = triSurf.faceNormals();

    // Encode tetrahedron signed distances into a bitset. For 4 tetrahedron
    // points: 1110 == positive, positive, positive, negative distance.
    // Barycentric triangulation of the cell is used:
    // - First distance is the distance at the cell center.
    // - Second, third and fourth distances are at two cell corner points.
    std::bitset<4> dists(pow(2, 5) - 1);

    // Legacy VTK output of the cut mesh geometry.
    geophase::vtkPolyDataOStream cutMeshStream(
        geophase::vtk_file_name("cutMesh", this->time().timeIndex()));

    // Correct the volume fractions geometrically in intersected cells.
    nTrianglesPerCell_ = 0;

    forAll(intersectedCellLabels_, cellJ)
    {
        const auto cellI = intersectedCellLabels_[cellJ];
        alpha[cellI] = 0;
        const auto& cutCell = meshCells[cellI];

        // Cell center and its distance as tetrahedron input.
        const auto& xC = cellCenters[cellI];
        const auto& distC = cellSignedDist[cellI];
        dists[0] = !std::signbit(distC);

        // For all faces of the cut cell,
        forAll(cutCell, faceL)
        {
            const label faceG = cutCell[faceL];

            // Face center & distance for tetrahedron input.
            const auto& xF = faceCenters[faceG];
            const scalar distF = sigDistCalc_.signedDistance(
                faceCenters[faceG]); // faceSignedDist_[faceG];
            dists[1] = !std::signbit(distF);

            const face& nBandFace = faces[faceG];
            for (label I = 0; I < nBandFace.size(); ++I)
            {
                // Current face-point & dist for tetrahedron input
                const label pointI0 = nBandFace[I];
                const scalar distI0 = pointSignedDist[pointI0];
                dists[2] = !std::signbit(distI0);
                const point& xI0 = meshPoints[pointI0];

                // Next face-point & distance for tetrahedron input
                const label pointI1 = nBandFace[(I + 1) % nBandFace.size()];
                const scalar distI1 = pointSignedDist[pointI1];
                dists[3] = !std::signbit(distI1);
                const point& xI1 = meshPoints[pointI1];

                // Tetrahedron centroid.
                const vector xT = 0.25 * (xC + xF + xI0 + xI1);

                // Tetrahedron sphere radius.
                const scalar radiusT = max(max(mag(xC - xT), mag(xF - xT)),
                    max(mag(xI0 - xT), mag(xI1 - xT)));

                // Triangles interesecting tetrahedron sphere.
                auto triangleLabels = octree.findSphere(xT, radiusT * radiusT);
                if (dists.all()) // If tetrahedron is inside surface.
                {
                    // Add the mixed product tet volume to the alpha cell value.
                    alpha[cellI] += (1. / 6.) *
                        std::abs(((xF - xC) & ((xI0 - xC) ^ (xI1 - xC))));
                }
                else if (!triangleLabels.empty())
                {
                    // Initialize the tetrahedron intersection.
                    geophase::foamVectorPolyhedron tetIntersection{
                        geophase::make_tetrahedron<
                            geophase::foamVectorPolyhedron>(xC, xF, xI0, xI1)};
                    nTrianglesPerCell_ += triangleLabels.size();
                    for (const auto& triangleL : triangleLabels)
                    {
                        tetIntersection = intersect_tolerance<
                            geophase::foamPolyhedronIntersection>(
                            tetIntersection,
                            foamHalfspace(triPoints[triSurf[triangleL][0]],
                                triNormals[triangleL]))
                                              .polyhedron();
                    }
                    // Add the volume of the intersection to the phase-specific
                    // volume.
                    alpha[cellI] += volume_by_surf_tri(tetIntersection);
                    if (this->writeGeometry())
                    {
                        cutMeshStream << tetIntersection;
                    }
                }
            }
        }
        alpha[cellI] /= cellVolumes[cellI];
    }
    nTrianglesPerCell_ /= intersectedCellLabels_.size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
surfaceMeshCellIntersection::surfaceMeshCellIntersection(
    const dictionary& configDict, const fvMesh& mesh)
    : volumeFractionCalculator{configDict, mesh},
      sigDistCalc_{configDict.subDict("distCalc"), mesh},
      intersectedCellLabels_(0), nTrianglesPerCell_{0}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void surfaceMeshCellIntersection::calcVolumeFraction(volScalarField& alpha)
{
    bulkVolumeFraction(alpha);
    findIntersectedCells();
    interfaceCellVolumeFraction(alpha);
}


void surfaceMeshCellIntersection::findIntersectedCells()
{
    const auto& mesh = this->mesh();
    const auto& meshCellPoints = mesh.cellPoints();
    const auto& cellClosestTriangle = sigDistCalc_.cellClosestPoint();
    const auto& cellSignedDist = sigDistCalc_.cellSignedDist();
    const auto& pointSignedDist = sigDistCalc_.pointSignedDist();

    intersectedCellLabels_.resize(0);

    forAll(cellSignedDist, cellI)
    {
        // No hit-> cell not within narrow band
        if (!cellClosestTriangle[cellI].hit())
        {
            continue;
        }

        // Todo (TT): this is essentially the sign refinement criterion.
        // Reuse it here.
        const auto& cellDist = cellSignedDist[cellI];
        const auto& cellPoints = meshCellPoints[cellI];

        forAll(cellPoints, pointI)
        {
            if ((pointSignedDist[cellPoints[pointI]] * cellDist) < 0)
            {
                intersectedCellLabels_.append(cellI);
                break;
            }
        }
    }
    assert((intersectedCellLabels_.size() < mesh.nCells()));
    assert(intersectedCellLabels_.size() > 0);
}


void surfaceMeshCellIntersection::writeFields() const
{
    sigDistCalc_.writeFields();

    // Write identified interface cells as field
    volScalarField interfaceCells{
        "interface_cells", sigDistCalc_.cellSignedDist()};
    interfaceCells = dimensionedScalar{"interface_cells", dimLength, 0};

    for (const auto cellI : intersectedCellLabels_)
    {
        interfaceCells[cellI] = 1.0;
    }

    interfaceCells.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
