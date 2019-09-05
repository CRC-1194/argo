#include "polynomialVofInitialization.hpp"

#include "orientedPlane.hpp"
#include "tetVofCalculator.hpp"
#include "triSurfaceAdapter.hpp"

#include "fvc.H"
#include "fvm.H"

#include <algorithm>
#include <cassert>
#include <map>
#include <omp.h>

namespace Foam {
namespace PolynomialVof {

// Private member functions
void polynomialVofInitialization::setBulkFractions(volScalarField& alpha) const
{
    // NOTE: the diffused signed distance field should only be used
    // where the original signed distance has not been computed.
    // Diffusion may cause cells in interface proximity to change sign (TT).
    forAll(alpha, idx)
    {
        if (signedDistance0_[idx] != 0.0)
        {
            alpha[idx] = 0.5*(sign(signedDistance0_[idx]) + 1.0);
        }
        else
        {
            alpha[idx] = 0.5*(sign(signedDistance_[idx]) + 1.0);
        }
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
                &&
                // The conditions below are only valid for a piecewise planar
                // interface (TT)
                (
                    signSwitches(idx)
                    ||
                    cellContainsVertex(idx)
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

bool polynomialVofInitialization::signSwitches(const label cell_id) const
{
    const auto& vertex_ids = mesh_.cellPoints()[cell_id];
    scalar first_dist = vertexSignedDistance_[vertex_ids[0]];

    for (const auto v_id : vertex_ids)
    {
        if (first_dist*vertexSignedDistance_[v_id] < 0.0)
        {
            return true;
        }
    }

    return false;
}

bool polynomialVofInitialization::cellContainsVertex(const label cell_id) const
{
    const auto& tris = surface_.localFaces();
    const auto& vertices = surface_.localPoints();

    auto hitInfo = cellNearestTriangle_[cell_id];

    if (hitInfo.hit())
    {
        const auto aTri = tris[hitInfo.index()];

        for (const auto v_id : aTri)
        {
            if (mesh_.pointInCell(vertices[v_id], cell_id))
            {
                return true;
            }
        }
    }

    return false;
}

void polynomialVofInitialization::calcVertexSignedDistance()
{
    // Compute the signed distance at mesh vertices in the narrow band
    const auto& points = mesh_.points();
    const auto& map_cell_to_vertex = mesh_.cellPoints();
    const auto& triPoints = surface_.points();  
    const auto& triNormals = surface_.faceNormals(); 

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
                auto hit_info = triSearch_.nearest(points[v_id], vector{1.0e8, 1.0e8, 1.0e8});
                vertexSignedDistance_[v_id] = 
                    (points[v_id] - triPoints[surface_[hit_info.index()][0]]) & 
                    triNormals[hit_info.index()];
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
                auto hit_info = triSearch_.nearest(points[face_id], vector{1.0e8, 1.0e8, 1.0e8});
                faceSignedDistance_[face_id] = 
                    (points[face_id] - triPoints[surface_[hit_info.index()][0]]) & 
                    triNormals[hit_info.index()];
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

triSurface polynomialVofInitialization::surfaceSubset(const label cell_id) const
{
    const auto& C = mesh_.C();

    // TODO: rethink search: bbox vs sphere, required box length / radius (TT)
    auto trisInSphere = triSearch_.tree().findSphere(C[cell_id], 2.25*sqrSearchDist_[cell_id]);
    boolList includeTri(surface_.size(), false);
    labelList pointMap{};
    labelList faceMap{};
    for (const auto idx : trisInSphere)
    {
        includeTri[idx] = true;
    }

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
    sqrDistFactor_{max(3.0, sqrDistFactor)},
    cellNearestTriangle_{},
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
            "signedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            wo 
        ),
        mesh,
        dimensionedScalar("signedDist", dimLength, 0.0),
        "zeroGradient"
    ),
    signedDistance0_("signedDistance0", signedDistance_), 
    faceSignedDistance_{mesh_.faces().size(), 1.0e15},
    vertexSignedDistance_{mesh_.points().size(), 1.0e15},
    interfaceCells_{},
    max_refinement_level_{max_refine},
    distances_initialized_{false}
{
    cellNearestTriangle_.reserve(signedDistance_.size());
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

//- Computation
void polynomialVofInitialization::calcSqrSearchDist()
{
    sqrSearchDist_ = fvc::average(pow(mesh_.deltaCoeffs(), -2));
}

void polynomialVofInitialization::calcSignedDist()
{
    // Zero the signed distance.
    signedDistance_ = dimensionedScalar{"signedDist", dimLength, 0};

    // Use the octree and the square search distance multiplied by a distance factor 
    // to build the surface mesh / volume mesh proximity information list. 
    triSearch_.findNearest(
        mesh_.C(),
        sqrDistFactor_ * sqrDistFactor_ * sqrSearchDist_,
        cellNearestTriangle_
    );

    // Compute the signed distance in each cell as the distance between the triangle
    // nearest to the cell center and the cell center.
    const volVectorField& C = mesh_.C();  
    const pointField& triPoints = surface_.points();  
    const vectorField& triNormals = surface_.faceNormals(); 
    forAll(cellNearestTriangle_, cellI)
    {
        const pointIndexHit& cellHit = cellNearestTriangle_[cellI];

        if (cellHit.hit()) 
        {
            signedDistance_[cellI] = 
                (C[cellI] - triPoints[surface_[cellHit.index()][0]]) & 
                triNormals[cellHit.index()];
        }
    }

    // Save the signed distance field given by the octree.
    signedDistance0_ = signedDistance_; 

    // Propagate the sign information into the bulk by solving a Laplace
    // equation for a single iteration for the signed distance field. 
    fvScalarMatrix distEqn
    (
        -fvm::laplacian(signedDistance_)
    );
    distEqn.solve(); 
}

void polynomialVofInitialization::initializeDistances()
{
    if (distances_initialized_) return;

    Info << "Computing signed distance in narrow band..." << endl;

    calcSqrSearchDist();
    calcSignedDist();
    calcVertexSignedDistance();
    identifyInterfaceCells();
    calcFaceSignedDistance();

    distances_initialized_ = true;
}

void polynomialVofInitialization::calcVolFraction(volScalarField& alpha)
{
    initializeDistances();

    setBulkFractions(alpha);

    Info << "Computing volume fraction for interface cells..." << endl;
    Info << "Number of cells flagged as interface cells: "
         << interfaceCells_.size() << endl;

    const auto& V = mesh_.V();

    #pragma omp parallel for
    for (const auto cell_id : interfaceCells_)
    {
        auto subsetSurface= surfaceSubset(cell_id);
        triSurfaceSearch subsetSearch{subsetSurface};
        scalar s = 3.0*Foam::sqrt(sqrSearchDist_[cell_id]);
        triSurfaceAdapter adapter{subsetSurface, subsetSearch, vector{s, s, s}};

        auto [tets, points, signed_dist] = decomposeCell(cell_id);

        adaptiveTetCellRefinement<triSurfaceAdapter> refiner{adapter, points, signed_dist, tets, max_refinement_level_};
        tetVofCalculator vofCalc{};
        alpha[cell_id] = vofCalc.accumulated_omega_plus_volume(refiner.resulting_tets(), refiner.signed_distance(), refiner.points()) / V[cell_id]; 

        // Limit volume fraction field
        alpha[cell_id] = max(min(alpha[cell_id], 1.0), 0.0);
    }

    Info << "Finished volume fraction calculation" << endl;
}

//- Write
void polynomialVofInitialization::writeFields() const
{
    sqrSearchDist_.write();
    signedDistance_.write();
    signedDistance0_.write();

    // Write identified interface cells as field
    volScalarField interfaceCells{"interface_cells", signedDistance_};
    interfaceCells = dimensionedScalar{"interface_cells", dimLength, 0};

    for(const auto idx : interfaceCells_)
    {
        interfaceCells[idx] = 1.0;
    }

    interfaceCells.write();
}

// End namespace PolynomialVof
}
// End namespace Foam
}
