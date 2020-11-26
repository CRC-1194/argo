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

#include "insideOutsidePropagation.hpp"

#include "gaussLaplacianScheme.H"
#include "linear.H"
#include "skewCorrectedSnGrad.H"

namespace Foam {
namespace SigDistCalc {

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<volScalarField> insideOutsidePropagation::propagate_inside_outside(const volScalarField& signed_distance) const
{
    // Hardcode Laplacian scheme and solver settings
    tmp<volScalarField> tmp_in_out_field{new volScalarField{signed_distance}};
    auto& in_out_field = tmp_in_out_field.ref();
    in_out_field.rename("inside_outside");

    surfaceScalarField Gamma
    (
        IOobject
        (
            "Diffusion coefficient",
            signed_distance.time().constant(),
            signed_distance.mesh(),
            IOobject::NO_READ
        ),
        signed_distance.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    fvScalarMatrix propagationEqn
    (
        fv::gaussLaplacianScheme<scalar,scalar>
        (
            in_out_field.mesh(),
            tmp<surfaceInterpolationScheme<scalar>>{new linear<scalar>{in_out_field.mesh()}}, 
            tmp<fv::snGradScheme<scalar>>{new fv::skewCorrectedSnGrad<scalar>(in_out_field.mesh())}
        ).fvmLaplacian(Gamma, in_out_field)
    );

    IStringStream solverInput{
        "solver PCG; preconditioner DIC; tolerance 0.5; relTol 0; minIter 1; maxIter 1;"
    };
    dictionary solverControl{solverInput};

    for (auto iteration = 0; iteration != 3; ++iteration)
    {
        propagationEqn.solve(solverControl);

        // Reset signed distance in the narrow band
        forAll(signed_distance, idx)
        {
            if (signed_distance[idx] != 0.0)
            {
                in_out_field[idx] = signed_distance[idx];
            }
        }
    }

    return tmp_in_out_field;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SigDistCalc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
