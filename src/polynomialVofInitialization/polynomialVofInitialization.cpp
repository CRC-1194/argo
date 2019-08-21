#include "polynomialVofInitialization.hpp"

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
    const volScalarField& propagatedSignedDistance,
    const volScalarField& originalSignedDistance
)
:
    mesh_{mesh},
    surface_{surface},
    signedDistance_{propagatedSignedDistance},
    signedDistance0_{originalSignedDistance},
    interfaceCells_{}
{
}

// Public member functions
void polynomialVofInitialization::computeVolumeFractions(volScalarField& alpha)
{
    setBulkFractions(alpha);
    
    // Identify interface cells
    
}
    
// End namespace PolynomialVof
}
// End namespace Foam
}
