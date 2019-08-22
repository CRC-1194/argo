#include "polynomialVofInitialization.hpp"

#include "fvc.H"
#include "fvm.H"
#include "pointFields.H"

namespace Foam {
namespace PolynomialVof {

// Private member functions
void polynomialVofInitialization::setBulkFractions(volScalarField& alpha) const
{
    forAll(alpha, idx)
    {
        alpha[idx] = 0.5*(sign(signedDistance_[idx]) + 1.0);
    }
}

void polynomialVofInitialization::identifyInterfaceCells()
{
    const auto& cell_vertices = mesh_.points();
    const auto& map_cell_to_vertices = mesh_.cellPoints();
    const auto& cell_centres = mesh_.C();

    forAll(cellNearestTriangle_, idx)
    {
        if (cellNearestTriangle_[idx].hit())
        {
            if (intersectionPossible(cell_centres[idx],
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
    // Compute the signed distance at mesh vertices. Only vertices of cells,
    // which previously have been identified as interface cells, are considered.
    const auto& points = pMesh_.mesh().points();
    const auto& map_cell_to_vertex = mesh_.cellPoints();
    const pointField& triPoints = surface_.points();  
    const vectorField& triNormals = surface_.faceNormals(); 

    for (const auto cell_id : interfaceCells_)
    {
        const auto& vertex_ids = map_cell_to_vertex[cell_id];

        for (const auto v_id : vertex_ids)
        {
            // Avoid duplicate computation (TT)
            if (vertexSignedDistance_[v_id] == 0.0)
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

// Constructors
polynomialVofInitialization::polynomialVofInitialization
(
    const fvMesh& mesh, 
    const triSurface& surface,
    const scalar sqrDistFactor,
    const IOobject::writeOption& wo  // Allows output for testing purposes.
)
:
    mesh_{mesh},
    runTime_{mesh_.time()},
    surface_{surface},
    pMesh_{mesh},
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
    vertexSignedDistance_   
    (
        IOobject
        (
            "vertexSignedDistance",
            runTime_.timeName(),
            mesh,
            IOobject::NO_READ,
            wo
        ),
        pMesh_,
        dimensionedScalar{"signedDistance", dimLength, 0.0},
        "zeroGradient"
    ),
    interfaceCells_{}
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

void polynomialVofInitialization::calcVolFraction(volScalarField& alpha)
{
    Info << "Computing signed distance in narrow band..." << endl;
    calcSqrSearchDist();
    calcSignedDist();
    setBulkFractions(alpha);

    Info << "Identifiying interface cells..." << endl;
    identifyInterfaceCells();

    Info << "Computing vertex signed distance for interface cells..." << endl;
    calcVertexSignedDistance();

    Info << "Computing volume fraction for interface cells..." << endl;
    const auto& C = mesh_.C();
    const auto& map_cell_to_vertex = mesh_.cellPoints();

    for (const auto cell_id : interfaceCells_)
    {
        // TODO: continue here
    }
}

//- Write
void polynomialVofInitialization::writeFields() const
{
    sqrSearchDist_.write();
    signedDistance_.write();
    signedDistance0_.write();
    vertexSignedDistance_.write();

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
