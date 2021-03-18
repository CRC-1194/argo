#include "polynomialVofInitialization.hpp"

#include "insideOutsidePropagation.hpp"
#include "orientedPlane.hpp"
#include "tetVofCalculator.hpp"
#include "triSurfaceAdapter.hpp"

#include "fvc.H"
#include "fvm.H"

#include <algorithm>
#include <cassert>
#include <map>
#include <omp.h>

namespace Foam::TriSurfaceImmersion {

// Private member functions
void polynomialVofInitialization::setBulkFractions(volScalarField& alpha) const
{
    forAll(alpha, idx)
    {
        alpha[idx] = pos(signedDistance_[idx]);
    }
}

void polynomialVofInitialization::identifyInterfaceCells()
{
    const auto& cell_vertices = mesh_.points();
    const auto& map_cell_to_vertices = mesh_.cellPoints();
    const auto& cell_centres = mesh_.C();

    forAll(cellNearestTriangle_, idx)
    {
        if (
            cellNearestTriangle_[idx].hit()
            && 
            intersectionPossible(cell_centres[idx],
             signedDistance0_[idx]*signedDistance0_[idx],
             map_cell_to_vertices[idx],
             cell_vertices
            )
        )
        {
            interfaceCells_.push_back(idx);
        }
    }
}

bool polynomialVofInitialization::intersectionPossible
     (
        const point& centre,
        const scalar distSqr,
        const labelList& point_IDs,
        const pointList& points
     ) const
{
    // Idea: if all vertices of a cell have a distance less than the signed distance
    // of the centre to the interface, then there is no intersection possible (TT)
    for(const auto p_id : point_IDs)
    {
        auto v = points[p_id] - centre;
        if ((v & v) >= distSqr)
        {
            return true;
        }
    }

    return false;
}

void polynomialVofInitialization::calcVertexSignedDistance()
{
    // Compute the signed distance at mesh vertices in the narrow band
    const auto& points = mesh_.points();
    const auto& map_cell_to_vertex = mesh_.cellPoints();

    forAll (signedDistance0_, cell_id)
    {
        // Skip cells outside of narrow band
        if (signedDistance0_[cell_id] == 0.0)
        {
            continue;
        }

        const auto& vertex_ids = map_cell_to_vertex[cell_id];

        for (const auto v_id : vertex_ids)
        {
            // Avoid duplicate computation (TT)
            if (vertexSignedDistance_[v_id] == 1.0e15)
            {
                // Vertices are guaranteed to be close to the interface.
                // So the search distance can be large (TT)
                auto [hit_info, distance] = sig_dist_calc_.signedDistance(points[v_id], 1e15);
                vertexNearestTriangle_[v_id] = hit_info;
                vertexSignedDistance_[v_id] = distance;
            }
        }
    }
}

void polynomialVofInitialization::calcFaceSignedDistance()
{
    const auto& points = mesh_.Cf();
    const auto& cells = mesh_.cells();
    const pointField& triPoints = surface_.points();  
    const vectorField& triNormals = surface_.faceNormals(); 

    for (const auto cell_id : interfaceCells_)
    {
        for (const auto face_id : cells[cell_id])
        {
            if (faceSignedDistance_[face_id] == 1.0e15)
            {
                faceSignedDistance_[face_id] = sig_dist_calc_.signedDistance(points[face_id]);
            }
        }
    }
}

polynomialVofInitialization::cellDecompositionTuple polynomialVofInitialization::decomposeCell(const label cell_id) const
{
    const auto& thisCell = mesh_.cells()[cell_id];
    const auto& cellVertexIDs = mesh_.cellPoints()[cell_id];
    const auto& vertices = mesh_.points();

    std::vector<indexedTet> tets(n_tets(cell_id));
    // Using a barycentric decomposition, the number of unique points
    // is the sum of n_cell_vertices + n_cell_faces + 1 (the cell centre) (TT).
    std::vector<point> points(cellVertexIDs.size() + thisCell.size() + 1);
    std::vector<scalar> sd(points.size());
    std::map<label, label> globalToLocal{};
    
    // Add vertices to points and their signed distance
    forAll(cellVertexIDs, idx)
    {
        points[idx] = vertices[cellVertexIDs[idx]];
        sd[idx] = vertexSignedDistance_[cellVertexIDs[idx]];
        globalToLocal[cellVertexIDs[idx]] = idx;
    }

    // Add the cell centre
    label centre_id = cellVertexIDs.size();
    points[centre_id] = mesh_.C()[cell_id];
    sd[centre_id] = signedDistance0_[cell_id];

    // Add face centres and build the indexed tets
    const auto& faces = mesh_.faces();
    label face_centre_id = centre_id + 1;
    label idx_tet = 0;
    for (const auto face_id : thisCell)
    {
        points[face_centre_id] = mesh_.Cf()[face_id];
        sd[face_centre_id] = faceSignedDistance_[face_id];

        for (const auto& anEdge : faces[face_id].edges())
        {
            tets[idx_tet] = indexedTet{centre_id, face_centre_id,
                             globalToLocal[anEdge[0]], globalToLocal[anEdge[1]]};
            ++idx_tet;
        }
        ++face_centre_id;
    }

    // Signed distance plausibility check
    for (uint idx = 0; idx != points.size(); ++idx)
    {
        assert(mag(sd[idx]) < mesh_.bounds().mag());
    }

    return std::make_tuple(tets, points, sd);
}

label polynomialVofInitialization::n_tets(const label cell_id) const
{
    label n_tet = 0;

    const auto& thisCell = mesh_.cells()[cell_id];
    const auto& faces = mesh_.faces();

    for (const auto face_id : thisCell)
    {
        n_tet += faces[face_id].nEdges();
    }

    return n_tet;
}

polynomialVofInitialization::searchSphere polynomialVofInitialization::cellInterfaceSearchSphere(const label cell_id) const
{
    const auto& cellToVertex = mesh_.cellPoints()[cell_id];
    std::vector<vector> closestPoints(cellToVertex.size());

    for (auto idx = 0; idx != cellToVertex.size(); ++idx)
    {
        closestPoints[idx] = vertexNearestTriangle_[cellToVertex[idx]].hitPoint();
    }

    point centre = std::accumulate(closestPoints.begin(), closestPoints.end(), vector{0,0,0})/closestPoints.size();

    scalar radiusSquared = 0.0;

    for (const auto v : closestPoints)
    {
        radiusSquared = std::max(radiusSquared, (v - centre)&(v - centre)); 
    }

    // Special case: all cell vertices have the same closest point.
    //      Thus, the bounding sphere of this point set has a radius of
    //      zero. In this case, relate the search radius to the edge length
    //      of the closest triangle. (TT)
    if (radiusSquared < SMALL)
    {
        const auto& v = surface_.points();
        const auto& t = surface_[cellNearestTriangle_[cell_id].index()];
        radiusSquared = ((v[t[0]] - centre)&(v[t[0]] - centre)) +
                        ((v[t[1]] - centre)&(v[t[1]] - centre)) +
                        ((v[t[2]] - centre)&(v[t[2]] - centre));
    }

    assert(radiusSquared > 1.0e-15 && "Radius of search square is zero.");

    return searchSphere{centre, radiusSquared};
}

triSurface polynomialVofInitialization::surfaceSubset(const label cell_id) const
{
    auto boundingSphere = cellInterfaceSearchSphere(cell_id);
    auto trisInSphere = triSearch_.tree().findSphere(boundingSphere.centre, boundingSphere.radiusSquared);
    boolList includeTri(surface_.size(), false);
    labelList pointMap{};
    labelList faceMap{};
    for (const auto idx : trisInSphere)
    {
        includeTri[idx] = true;
    }

    assert(trisInSphere.size() > 0 && "Surface subset is empty set.");

    return surface_.subsetMesh(includeTri, pointMap, faceMap);
}

// Constructors
polynomialVofInitialization::polynomialVofInitialization
(
    const fvMesh& mesh, 
    const triSurface& surface,
    const scalar sqrDistFactor,
    const IOobject::writeOption& wo,  // Allows output for testing purposes.
    const label max_refine = -1
)
:
    mesh_{mesh},
    runTime_{mesh_.time()},
    surface_{surface},
    triSearch_{surface},
    sig_dist_calc_{surface},
    sqrDistFactor_{max(3.0, sqrDistFactor)},
    cellNearestTriangle_{},
    vertexNearestTriangle_{},
    sqrSearchDist_
    (
        IOobject
        (
            "sqrSearchDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        fvc::average(pow(mesh.deltaCoeffs(), -2))
    ),
    signedDistance_
    (
        IOobject
        (
            "cellSignedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        mesh,
        dimensionedScalar("signedDist", dimLength, 0.0),
        "zeroGradient"
    ),
    signedDistance0_("cellSignedDistance0", signedDistance_), 
    faceSignedDistance_{mesh_.faces().size(), 1.0e15},
    vertexSignedDistance_{mesh_.points().size(), 1.0e15},
    interfaceCells_{},
    max_refinement_level_{max_refine},
    distances_initialized_{false}
{
    cellNearestTriangle_.reserve(signedDistance_.size());
    vertexNearestTriangle_.reserve(vertexSignedDistance_.size());
}

// Public member functions
//- Access
const fvMesh& polynomialVofInitialization::mesh() const
{
    return mesh_;
}

const Time& polynomialVofInitialization::time() const
{
    return runTime_;
}

const triSurface& polynomialVofInitialization::surface() const
{
    return surface_;
}

const volScalarField& polynomialVofInitialization::sqrSearchDist() const
{
    return sqrSearchDist_;
}

const volScalarField& polynomialVofInitialization::signedDistance() const
{
    return signedDistance_;
}

const volScalarField& polynomialVofInitialization::signedDistance0() const
{
    return signedDistance0_;
}

label polynomialVofInitialization::maxRefinementLevel() const
{
    return max_used_refinement_level_;
}

//- Computation
void polynomialVofInitialization::calcSqrSearchDist()
{
    sqrSearchDist_ = fvc::average(pow(mesh_.deltaCoeffs(), -2));
}

void polynomialVofInitialization::calcSignedDist()
{
    signedDistance0_.primitiveFieldRef() = sig_dist_calc_.signedDistance(cellNearestTriangle_, mesh_.C(), sqrSearchDist_*sqrDistFactor_*sqrDistFactor_, 0.0);

    signedDistance_ = insideOutsidePropagation::propagateInsideOutside(signedDistance0_);
}

void polynomialVofInitialization::initializeDistances()
{
    if (distances_initialized_) return;

    Info << "Computing signed distance in narrow band..." << endl;

    calcSqrSearchDist();
    calcSignedDist();
    calcVertexSignedDistance();

    Info << "Identifying interface cells..." << endl;

    identifyInterfaceCells();
    calcFaceSignedDistance();

    distances_initialized_ = true;
}

void polynomialVofInitialization::calcVolFraction(volScalarField& alpha, const bool writeTets)
{
    initializeDistances();

    setBulkFractions(alpha);

    Info << "Computing volume fraction for interface cells..." << endl;
    Info << "Number of cells flagged as interface cells: "
         << interfaceCells_.size() << endl;

    const auto& V = mesh_.V();
    label max_refine = 0;

    // TODO: OpenMP disbaled for now. Loop does not execute correct with
    // more than two threads. See issue on GitLab (TT)
    //#pragma omp parallel for reduction(max:max_refine)
    for (const auto cell_id : interfaceCells_)
    {
        auto subsetSurface= surfaceSubset(cell_id);
        triSurfaceSearch subsetSearch{subsetSurface};
        scalar s = 2.0*sqrDistFactor_*Foam::sqrt(sqrSearchDist_[cell_id]);
        triSurfaceAdapter adapter{subsetSurface, subsetSearch, vector{s, s, s}};
        
        auto [tets, points, signed_dist] = decomposeCell(cell_id);

        adaptiveTetCellRefinement<triSurfaceAdapter> refiner
                                                     {
                                                         adapter,
                                                         points,
                                                         signed_dist,
                                                         tets,
                                                         max_refinement_level_,
                                                         writeTets,
                                                         cell_id
                                                     };
        tetVofCalculator vofCalc{};
        alpha[cell_id] = vofCalc.accumulated_omega_plus_volume(refiner.resulting_tets(), refiner.signed_distance(), refiner.points()) / V[cell_id]; 

        // Limit volume fraction field
        alpha[cell_id] = max(min(alpha[cell_id], 1.0), 0.0);

        max_refine = std::max(refiner.refinement_level(), max_refine);
    }

    max_used_refinement_level_ = max_refine;
    Info << "Finished volume fraction calculation" << endl;
}

//- Write
void polynomialVofInitialization::writeFields() const
{
    sqrSearchDist_.write();
    signedDistance_.write();
    signedDistance0_.write();

    // Write identified interface cells as field
    volScalarField interfaceCells{"interfaceCells", signedDistance_};
    interfaceCells = dimensionedScalar{"interfaceCells", dimLength, 0};

    for(const auto idx : interfaceCells_)
    {
        interfaceCells[idx] = 1.0;
    }

    interfaceCells.write();
}

} // End namespace Foam::TriSurfaceImmersion
