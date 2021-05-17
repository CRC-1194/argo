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

#include "surfaceMeshCellApproximation.hpp"

#include <algorithm>
#include <cassert>

#include "addToRunTimeSelectionTable.H"

#include "IntersectionCriteria.hpp"
#include "signedDistanceCalculator.hpp"
#include "tetVolumeFractionCalculator.hpp"
#include "triSurfaceDistCalc.hpp"

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(surfaceMeshCellApproximation, 0);
addToRunTimeSelectionTable(
    volumeFractionCalculator, surfaceMeshCellApproximation, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
surfaceMeshCellApproximation::cellDecompositionTuple
surfaceMeshCellApproximation::decomposeCell(const label cellID) const
{
    const auto& mesh = this->mesh();
    const auto& thisCell = mesh.cells()[cellID];
    const auto& cellVertexIDs = mesh.cellPoints()[cellID];
    const auto& vertices = mesh.points();
    const auto& pointSignedDist = sigDistCalcPtr_->pointSignedDist();
    const auto& cellSignedDist = sigDistCalcPtr_->cellSignedDist0();

    std::vector<indexedTet> tets(nTets(cellID));

    // Using a barycentric decomposition, the number of unique points
    // is the sum of n_cell_vertices + n_cell_faces + 1 (the cell centre) (TT).
    std::vector<point> points(cellVertexIDs.size() + thisCell.size() + 1);
    std::vector<scalar> sd(points.size());
    std::map<label, label> globalToLocal{};

    // Add vertices to points and their signed distance
    forAll(cellVertexIDs, I)
    {
        points[I] = vertices[cellVertexIDs[I]];
        sd[I] = pointSignedDist[cellVertexIDs[I]];
        globalToLocal[cellVertexIDs[I]] = I;
    }

    // Add the cell centre
    label centre_id = cellVertexIDs.size();
    points[centre_id] = mesh.C()[cellID];
    sd[centre_id] = cellSignedDist[cellID];

    // Add face centres and build the indexed tets
    const auto& faces = mesh.faces();
    label face_centre_id = centre_id + 1;
    label idx_tet = 0;
    for (const auto face_id : thisCell)
    {
        points[face_centre_id] = mesh.Cf()[face_id];
        sd[face_centre_id] =
            sigDistCalcPtr_->signedDistance(mesh.Cf()[face_id]);

        for (const auto& anEdge : faces[face_id].edges())
        {
            tets[idx_tet] = indexedTet{centre_id,
                face_centre_id,
                globalToLocal[anEdge[0]],
                globalToLocal[anEdge[1]]};
            ++idx_tet;
        }
        ++face_centre_id;
    }

    // Signed distance plausibility check
    for (uint idx = 0; idx != points.size(); ++idx)
    {
        assert(mag(sd[idx]) < mesh.bounds().mag());
    }

    return std::make_tuple(tets, points, sd);
}


label surfaceMeshCellApproximation::nTets(const label cellID) const
{
    label nTet = 0;

    const auto& thisCell = this->mesh().cells()[cellID];
    const auto& faces = this->mesh().faces();

    for (const auto faceID : thisCell)
    {
        nTet += faces[faceID].nEdges();
    }

    return nTet;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
surfaceMeshCellApproximation::surfaceMeshCellApproximation(
    const dictionary& configDict, const fvMesh& mesh)
    : volumeFractionCalculator{configDict, mesh},
      sigDistCalcPtr_{
          signedDistanceCalculator::New(configDict.subDict("distCalc"), mesh)},
      interfaceCellIDs_{}, maxAllowedRefinementLevel_{
                               configDict.get<label>("refinementLevel")}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void surfaceMeshCellApproximation::calcVolumeFraction(volScalarField& alpha)
{
    bulkVolumeFraction(alpha);
    findIntersectedCells();

    Info << "Computing volume fraction for interface cells..." << endl;
    Info << "Number of cells flagged as interface cells: "
         << interfaceCellIDs_.size() << endl;

    const auto& V = this->mesh().V();
    label max_refine = 0;

    // TODO (TT): OpenMP disabled for now. Loop does not execute correct with
    // more than two threads. See issue on GitLab.
    //#pragma omp parallel for reduction(max:max_refine)
    for (const auto cellID : interfaceCellIDs_)
    {
        auto [tets, points, signed_dist] = decomposeCell(cellID);

        adaptiveTetCellRefinement<signedDistanceCalculator,
            boundingBallCriterion>
            refiner{this->sigDistCalc(),
                points,
                signed_dist,
                tets,
                maxAllowedRefinementLevel_};
        tetVolumeFractionCalculator vofCalc{};
        alpha[cellID] =
            vofCalc.accumulatedOmegaPlusVolume(refiner.resultingTets(),
                refiner.signedDistance(),
                refiner.points()) /
            V[cellID];

        // Bound volume fraction field
        alpha[cellID] = max(min(alpha[cellID], 1.0), 0.0);

        max_refine = std::max(refiner.refinementLevel(), max_refine);

        if (this->writeGeometry())
        {
            refiner.writeTets(cellID);
        }
    }

    maxUsedRefinementLevel_ = max_refine;

    Info << "Finished volume fraction calculation" << endl;
}


void surfaceMeshCellApproximation::findIntersectedCells()
{
    const auto& cellClosestPoint = sigDistCalcPtr_->cellClosestPoint();
    const auto& cellSignedDist = sigDistCalcPtr_->cellSignedDist0();
    const auto& centres = this->mesh().C();
    const auto& points = this->mesh().points();
    const auto& meshCellPoints = this->mesh().cellPoints();

    forAll(cellClosestPoint, cellI)
    {
        auto distSqr = pow(cellSignedDist[cellI], 2.0);

        if
        (
            cellClosestPoint[cellI].hit()
            &&
            considerIntersected(centres[cellI], distSqr, meshCellPoints[cellI],
                points, std::vector<scalar>{}, boundingBallCriterion{})
        )
        {
            interfaceCellIDs_.push_back(cellI);
        }
    }
}


void surfaceMeshCellApproximation::writeFields() const
{
    sigDistCalcPtr_->writeFields();

    // Write identified interface cells as field
    volScalarField interfaceCells{
        "interfaceCells", sigDistCalcPtr_->cellSignedDist()};
    interfaceCells = dimensionedScalar{"interfaceCells", dimLength, 0};

    for (const auto idx : interfaceCellIDs_)
    {
        interfaceCells[idx] = 1.0;
    }

    interfaceCells.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
