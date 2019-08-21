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
            const scalar,
            const IOobject::writeOption&
        );
        
        // Member functions

        //- Access
        const fvMesh& mesh() const;
        const Time& time() const;
        const triSurface& surface() const;
        const volScalarField& sqrSearchDist() const; 
        const volScalarField& signedDistance() const; 
        const volScalarField& signedDistance0() const; 

        //- Computation

        // Squared search distance calculation function. 
        void calcSqrSearchDist();  

        // Signed distance computation using triSurfaceSearch and an external triSurface.
        void calcSignedDist(const triSurface& tri, const triSurfaceSearch& triSearch);  

        // Signed distance computation using triSurfaceSearch and the internal triSurface.
        void calcSignedDist();  

        // Volume fraction calculation.
        virtual void calcVolFraction(volScalarField& alpha);

        //- Write 
        virtual void writeFields() const;
};

// End namespace PolynomialVof
}
// End namespace Foam
}

#endif
