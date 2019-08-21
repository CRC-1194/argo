#ifndef polynomialVofInitialization_hpp
#define polynomialVofInitialization_hpp


// OpenFOAM includes
#include "fvMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "volFields.H"

// std
#include <vector>

namespace Foam {
namespace PolynomialVof {

class polynomialVofInitialization
{
    private:
        
        // Data
        const fvMesh& mesh_;
        const Time& runTime_;
        const triSurface& surface_;
        
        // Factor used to extend the narrow band by N cells. 
        // If sqrDistanceFactor = 2, the narrow band is extended by 2 cells. 
        const scalar sqrDistFactor_; 

        DynamicList<pointIndexHit> cellNearestTriangle_;

        volScalarField sqrSearchDist_;  
        volScalarField signedDistance_;
        volScalarField signedDistance0_;

        std::vector<label> interfaceCells_;

        // Member functions
        void setBulkFractions(volScalarField&) const;

    public:

        // Constructors
        polynomialVofInitialization
        (
            const fvMesh&,
            const triSurface&,
            const volScalarField&,
            const volScalarField&
        );
        
        // Member functions
        void computeVolumeFractions(volScalarField&);
};

// End namespace PolynomialVof
}
// End namespace Foam
}

#endif
