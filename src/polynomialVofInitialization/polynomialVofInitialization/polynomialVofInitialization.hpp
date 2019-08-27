#ifndef polynomialVofInitialization_hpp
#define polynomialVofInitialization_hpp

#include "AdaptiveTetCellRefinement.hpp"
#include "tetVofCalculator.hpp"

// OpenFOAM includes
#include "fvMesh.H"
#include "pointList.H"
#include "pointMesh.H"
#include "surfaceFields.H"
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
        using cellDecompositionTuple = std::tuple<std::vector<indexedTet>,
                                                  std::vector<point>,
                                                  std::vector<scalar>>;
        
        // Data
        const fvMesh& mesh_;
        const Time& runTime_;
        const triSurface& surface_;
        pointMesh pMesh_;
        triSurfaceSearch triSearch_;
        tetVofCalculator vofCalc_;
        
        // Factor used to extend the narrow band by N cells. 
        // If sqrDistanceFactor = 2, the narrow band is extended by 2 cells. 
        const scalar sqrDistFactor_; 

        DynamicList<pointIndexHit> cellNearestTriangle_;

        volScalarField sqrSearchDist_;  
        volScalarField signedDistance_;
        volScalarField signedDistance0_;
        surfaceScalarField faceSignedDistance_;
        pointScalarField vertexSignedDistance_;

        std::vector<label> interfaceCells_;

        // Member functions
        void setBulkFractions(volScalarField&) const;
        void identifyInterfaceCells();
        bool intersectionPossible
             (
                const point& centre,
                const scalar distSqr,
                const labelList& point_IDs,
                const pointList& points
             ) const;
        void calcVertexSignedDistance();
        void calcFaceSignedDistance();
        cellDecompositionTuple decomposeCell(const label cell_id) const;
        label n_tets(const label cell_id) const;


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

        // Signed distance computation using triSurfaceSearch and the internal triSurface.
        void calcSignedDist();  

        // Volume fraction calculation.
        void calcVolFraction(volScalarField& alpha);

        //- Write 
        void writeFields() const;
};

// End namespace PolynomialVof
}
// End namespace Foam
}

#endif
