#include "polynomialVofInitialization.hpp"

#include "fvc.H"
#include "fvm.H"

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
        dimensionedScalar("signedDist", dimLength,0),
        "zeroGradient"
    ),
    signedDistance0_("signedDistance0", signedDistance_), 
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

    triSurfaceSearch triSearch{surface_};

    // Use the octree and the square search distance multiplied by a distance factor 
    // to build the surface mesh / volume mesh proximity information list. 
    triSearch.findNearest(
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
    calcSqrSearchDist();
    calcSignedDist();
    setBulkFractions(alpha);
}

//- Write
void polynomialVofInitialization::writeFields() const
{
    sqrSearchDist_.write();
    signedDistance_.write();
    signedDistance0_.write();

    // Write identified interface cells as field
    /*
    volScalarField interfaceCells{"interface_cells", signedDistance_};
    interfaceCells = dimensionedScalar{"interface_cells", dimLength, 0};

    for(const auto idx : interfaceCells_)
    {
        interfaceCells[idx] = 1.0;
    }

    interfaceCells.write();
    */
}

// End namespace PolynomialVof
}
// End namespace Foam
}
