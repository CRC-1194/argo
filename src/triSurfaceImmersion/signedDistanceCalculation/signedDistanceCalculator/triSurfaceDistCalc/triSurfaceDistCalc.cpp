/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tobias Tolle, TU Darmstadt
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

#include "triSurfaceDistCalc.hpp"

#include "addToRunTimeSelectionTable.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "triSurfaceTools.H"

#include "IntersectionCriteria.hpp"
#include "insideOutsidePropagation.hpp"

#include <algorithm>

namespace Foam::TriSurfaceImmersion
{

defineTypeNameAndDebug(triSurfaceDistCalc, 0);
addToRunTimeSelectionTable(
    signedDistanceCalculator, triSurfaceDistCalc, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void triSurfaceDistCalc::computeVertexNormals()
{
    /* The vertex normals are computed as a weighted sum of normals
     * of the adjacent triangles. The weights are the triangle angles at the
     * vertex at hand.
     * See
     *      "Computing Vertex Normals from Polygonal Facets",
     *      G. Thürmer & C. Wüthrich (2012)
     *      https://doi.org/10.1080/10867651.1998.10487487
     * for details. There is a proof that this gives correct inside/outside
     * information by J. Baerentzen and H. Aanaes, Technical University of
     * Denmark
     */
    const auto& tri_normals = surface_.faceNormals();
    const auto& vertex_to_faces = surface_.pointFaces();
    const auto& vertices = surface_.localPoints();

    forAll(vertices, v_id)
    {
        label vid_a{v_id};
        label vid_b{0};
        label vid_c{0};

        for (const auto fid : vertex_to_faces[v_id])
        {
            triSurfaceTools::otherVertices(surface_, fid, vid_a, vid_b, vid_c);
            vector v1{vertices[vid_b] - vertices[vid_a]};
            vector v2{vertices[vid_c] - vertices[vid_a]};
            scalar alpha{Foam::acos(
                std::clamp((v1 & v2) / (mag(v1) * mag(v2)), -1.0, 1.0))};

            vertexNormals_[v_id] += alpha * tri_normals[fid];
        }

        vertexNormals_[v_id] /= mag(vertexNormals_[v_id]) + SMALL;
    }
}


void triSurfaceDistCalc::computeSignedDistances()
{
    cellSignedDist0_.primitiveFieldRef() = signedDistance(cellNearestTriangle_,
        this->mesh().C(),
        searchDistCalc_.cellSqrSearchDist(),
        this->outOfNarrowBandValue());

    findIntersectedCells();

    cellSignedDist_ =
        insideOutsidePropagation::propagateInsideOutside(cellSignedDist0_);

    //findSpuriousSignSwitches();
    fixSpuriousSigns();

    pointSignedDist_.primitiveFieldRef() = signedDistance(pointNearestTriangle_,
        this->mesh().points(),
        searchDistCalc_.pointSqrSearchDist(),
        this->outOfNarrowBandValue());
}

// TODO: copied here from SMCA volume fraction calculation for testing. If it works,
// refactor! (TT)
void triSurfaceDistCalc::findIntersectedCells()
{
    const auto& centres = this->mesh().C();
    const auto& points = this->mesh().points();
    const auto& meshCellPoints = this->mesh().cellPoints();

    forAll(cellNearestTriangle_, cellI)
    {
        auto distSqr = pow(cellSignedDist0_[cellI], 2.0);

        if
        (
            cellNearestTriangle_[cellI].hit()
            &&
            considerIntersected(centres[cellI], distSqr, meshCellPoints[cellI],
                points, std::vector<scalar>{}, boundingBallCriterion{})
        )
        {
            isInterfaceCell_[cellI] = true;
        }
        else
        {
            isInterfaceCell_[cellI] = false;
        }
    }
}

std::vector<label> triSurfaceDistCalc::findSpuriousSignSwitches()
{
    std::vector<label> signSwitchFaces{};

    const auto& o = this->mesh().owner();
    const auto& n = this->mesh().neighbour();

    forAll(n, fid)
    {
        auto flag = !isInterfaceCell_[o[fid]] && !isInterfaceCell_[n[fid]];
        if (flag && (cellSignedDist_[o[fid]]*cellSignedDist_[n[fid]] < 0.0))
        {
            signSwitchFaces.push_back(fid);
        }
    }

    Info<< "There are " << signSwitchFaces.size() << " spurious cells with sign switches."
        << endl;  

    return signSwitchFaces;
}


void triSurfaceDistCalc::fixSpuriousSigns()
{
    auto switchFaces = findSpuriousSignSwitches();
    const label maxFixIterations = 100;
    label count = 0;

    const auto& owner = this->mesh().owner();
    const auto& neighbour = this->mesh().neighbour();
    const auto& cellToVertices = this->mesh().cellPoints();
    const auto& vertexToCells = this->mesh().pointCells();
    const auto& V = this->mesh().V();

    while ((!switchFaces.empty()) && (count != maxFixIterations))
    {
        for (auto fid : switchFaces)
        {
            auto oid = owner[fid];
            auto nid = neighbour[fid];

            std::vector<label> ownerCells{};
            std::vector<label> neighbourCells{};

            for (auto vid : cellToVertices[oid])
            {
                for (auto cid : vertexToCells[vid])
                {
                    ownerCells.push_back(cid);
                }
            }

            for (auto vid : cellToVertices[nid])
            {
                for (auto cid : vertexToCells[vid])
                {
                    neighbourCells.push_back(cid);
                }
            }

            // Remove duplicates from cell ID lists
            std::sort(ownerCells.begin(), ownerCells.end());
            std::sort(neighbourCells.begin(), neighbourCells.end());

            auto last = std::unique(ownerCells.begin(), ownerCells.end());
            ownerCells.erase(last, ownerCells.end());
            last = std::unique(neighbourCells.begin(), neighbourCells.end());
            neighbourCells.erase(last, neighbourCells.end());

            // Create volume weighted vote on the sign to choose
            scalar ownerSignWeight = 0.0;
            scalar neighbourSignWeight = 0.0;

            for (auto cid : ownerCells)
            {
                if (isInterfaceCell_[cid])
                {
                    continue;
                }
                ownerSignWeight += V[cid]*sign(cellSignedDist_[cid])*pos(cellSignedDist_[oid]
                                        *cellSignedDist_[cid]);
            }

            for (auto cid : neighbourCells)
            {
                if (isInterfaceCell_[cid])
                {
                    continue;
                }
                neighbourSignWeight += V[cid]*sign(cellSignedDist_[cid])*pos(cellSignedDist_[nid]
                                        *cellSignedDist_[cid]);
            }


            // Case 1: switch sign for the face neighbour
            if (mag(ownerSignWeight) > mag(neighbourSignWeight))
            {
                cellSignedDist_[nid] *= -1.0;
                cellSignedDist0_[nid] *= -1.0; 
            }
            // Case 2: switch sign for the face owner
            else if (mag(ownerSignWeight) < mag(neighbourSignWeight))
            {
                cellSignedDist_[oid] *= -1.0; 
                cellSignedDist0_[oid] *= -1.0;
            }
        }

        switchFaces = findSpuriousSignSwitches();
        ++count;
    }

    // Debugging
    switchFaces = findSpuriousSignSwitches();
    spuriousCells_ = 0.0;

    for (auto fid : switchFaces)
    {
        spuriousCells_[owner[fid]] = 1.0;
        spuriousCells_[neighbour[fid]] = 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceDistCalc::triSurfaceDistCalc(
    const dictionary& configDict, const fvMesh& mesh)
    : signedDistanceCalculator{configDict, mesh}, searchDistCalc_{mesh,
                                                      this->narrowBandWidth()},
      surface_{mesh.time().path() + "/" + configDict.get<fileName>("surfaceFile")},
      surfaceSearch_{surface_}, vertexNormals_{
                                    surface_.nPoints(), vector{0, 0, 0}},
      isInterfaceCell_(mesh.nCells()),
      spuriousCells_{IOobject("spuriousCells",
                          mesh.time().timeName(),
                          mesh,
                          IOobject::NO_READ,
                          IOobject::NO_WRITE),
          mesh,
          dimensionedScalar("spuriousCells", dimless, 0),
          "zeroGradient"}
{
    computeVertexNormals();
    computeSignedDistances();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalarField triSurfaceDistCalc::signedDistance(
    DynamicList<pointIndexHit>& pointToNearestTriangle,
    const pointField& pf,
    const scalarField& searchDistSqr,
    const scalar outOfSearchDomain) const
{
    scalarField distances(pf.size(), outOfSearchDomain);

    // Use the octree and the square search distance to build the
    // surface mesh / point field proximity information list.
    pointToNearestTriangle.reserve(pf.size());
    surfaceSearch_.findNearest(pf, searchDistSqr, pointToNearestTriangle);

    forAll(pointToNearestTriangle, p_id)
    {
        const pointIndexHit& hitInfo = pointToNearestTriangle[p_id];

        if (hitInfo.hit())
        {
            vector delta_v{pf[p_id] - hitInfo.hitPoint()};
            auto snormal = normalAtSurface(hitInfo);
            distances[p_id] = mag(delta_v) * sign(snormal & delta_v);
        }
    }

    return distances;
}


scalarField triSurfaceDistCalc::signedDistance(const pointField& pf,
    const scalarField& searchDistSqr,
    const scalar outOfSearchDomain) const
{
    DynamicList<pointIndexHit> pointToNearestTriangle{};
    pointToNearestTriangle.reserve(pf.size());

    return signedDistance(
        pointToNearestTriangle, pf, searchDistSqr, outOfSearchDomain);
}


std::tuple<pointIndexHit, scalar> triSurfaceDistCalc::signedDistance(
    const point& p, const scalar searchDistSqr) const
{
    scalar distance{1e15};

    const auto hitInfo =
        surfaceSearch_.nearest(p, vector{Foam::sqrt(searchDistSqr), 0, 0});

    if (hitInfo.hit())
    {
        vector delta_v{p - hitInfo.hitPoint()};
        auto snormal = normalAtSurface(hitInfo);
        distance = mag(delta_v) * sign(snormal & delta_v);
    }

    return std::make_tuple(hitInfo, distance);
}


scalar triSurfaceDistCalc::signedDistance(const point& p) const
{
    return std::get<1>(signedDistance(p, 1e15));
}


vector triSurfaceDistCalc::normalAtSurface(const pointIndexHit& hitInfo) const
{
    vector normal{0, 0, 0};

    const auto& fnormals = surface_.faceNormals();
    const auto& v = surface_.localPoints();

    // Transformation to local coordinate system:
    // - origin: point a (first point of triangle)
    // - x-axis: point a to point b (second point of triangle)
    // - y-axis: cross product of triangle normal and ex
    const auto tri_hit = surface_.localFaces()[hitInfo.index()];
    const auto a_to_b = v[tri_hit[1]] - v[tri_hit[0]];
    const auto a_to_c = v[tri_hit[2]] - v[tri_hit[0]];
    const auto a_to_hit = hitInfo.hitPoint() - v[tri_hit[0]];
    const auto ex = a_to_b / mag(a_to_b);
    const auto ey = fnormals[hitInfo.index()] ^ ex;
    const vector b{a_to_b & ex, 0, 0};
    const vector c{a_to_c & ex, a_to_c & ey, 0};
    const vector h_l{a_to_hit & ex, a_to_hit & ey, 0};

    // This is a linear shape function for a triangle where the vertices are:
    // - 1) the origin (0,0)
    // - 2) a point on the x-axis (t,0)
    // - 3) an arbitrary point (u,v)
    normal = vertexNormals_[tri_hit[0]] *
            (1.0 - h_l.x() / b.x() + h_l.y() / c.y() * (c.x() / b.x() - 1.0)) +
        vertexNormals_[tri_hit[1]] *
            (h_l.x() / b.x() - h_l.y() * c.x() / (b.x() * c.y())) +
        vertexNormals_[tri_hit[2]] * (h_l.y() / c.y());

    return normal;
}


scalar triSurfaceDistCalc::referenceLength() const
{
    // Use minimum edge length as reference length
    const auto& edges = surface_.edges();
    const auto& vertices = surface_.localPoints();

    scalar min_length = 1.0e15;

    for (const auto& anEdge : edges)
    {
        min_length = min(
            min_length, mag(vertices[anEdge.start()] - vertices[anEdge.end()]));
    }

    return min_length;
}


label triSurfaceDistCalc::nSurfaceElements() const
{
    return surface_.size();
}


scalar triSurfaceDistCalc::surfaceEnclosedVolume() const
{
    scalar Vsurf = 0;

    const auto& surfacePoints = surface_.points();
    forAll(surface_, triangleI)
    {
        const auto& Sf = surface_.Sf()[triangleI];
        const auto& triangle = surface_[triangleI];
        Vsurf += dot(-Sf, // Surface normals are oriented into the phase.
            (surfacePoints[triangle[0]] + surfacePoints[triangle[1]] +
                surfacePoints[triangle[2]]));
    }

    Vsurf = 1. / 9. * mag(Vsurf);

    return Vsurf;
}


void triSurfaceDistCalc::writeFields() const
{
    this->signedDistanceCalculator::writeFields();
    searchDistCalc_.writeFields();
    spuriousCells_.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
