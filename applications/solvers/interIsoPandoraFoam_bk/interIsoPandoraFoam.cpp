/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
    Copyright (C) 2019 DLR
    Copyright (C) 2020 TU Darmstadt
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

Application
    interIsoPandoraFoam

Group
    grpMultiphaseSolvers

Description
    Solver derived from interIsoFoam for two incompressible, isothermal immiscible
    fluids 

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    isoAdvector code supplied by Johan Roenby, STROMNING (2018)

    Surface-tension force regularization by 
    Tomislav Maric, Tobias Tolle, and Anja Lippert (2020)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

#include "pandora.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    argList::addBoolOption
    (
	"YoungLaplace",
	"Solve Young-Laplace Equation in the first time step.",
	false 
    );

    // TODO (TT): including this file raises compilation errors I do not know how to fix 
    // yet. Left out for now as it is not required to make the solver functional
    //#include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.hpp"

    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        // In the first time step
        bool YoungLaplace = args.found("YoungLaplace");
        if (YoungLaplace)
        {
            // Solve the Young-Laplace equation
            // for a better estimate fo the 
            // initial pressure.
            Info << "Solving Young-Laplace equation to recover "
            << "the pressure jump in the first time step." << endl;
            #include "YoungLaplaceEqn.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                if (isA<dynamicRefineFvMesh>(mesh))
                {
                    advector.surf().reconstruct();
                }

                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (isA<dynamicRefineFvMesh>(mesh))
                    {
                        advector.surf().mapAlphaField();
                        alpha2 = 1.0 - alpha1;
                        alpha2.correctBoundaryConditions();
                        rho == alpha1*rho1 + alpha2*rho2;
                        rho.correctBoundaryConditions();
                        rho.oldTime() = rho;
                        alpha2.oldTime() = alpha2;
                    }

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            runTime.cpuTimeIncrement();
            #include "alphaControls.hpp"
            #include "alphaEqnSubCycle.hpp"

            // Compute the surface tension force after solving for
            // volume fractions. 
            // Reconstruct the interface in the new time step.
	    // FIXME: next time-step the reconstruction will be repeated, this 
	    //        will cause communication overhead and needs fixing. TM.
            advector.surf().reconstruct();

            fSigma = pandoraModel.surfaceTensionForce(alpha1);
            fSigmaCell = fvc::reconstruct(fSigma * mesh.magSf());

scalar sum0 = Zero;
scalar sum1 = Zero;
scalar sum2 = Zero;
forAll (fSigmaCell, i)
{
    if (fSigmaCell[i][0] > 0)
        sum0 += fSigmaCell[i][0];
    if (fSigmaCell[i][1] > 0)
        sum1 += fSigmaCell[i][1];
    if (fSigmaCell[i][2] > 0)
        sum2 += fSigmaCell[i][2];
}
reduce(sum0, sumOp<scalar>());
reduce(sum1, sumOp<scalar>());
reduce(sum2, sumOp<scalar>());
vector sum(sum0, sum1, sum2);
Info<<"sumfSigmaCell = "<<sum<<nl;

if (runTime.writeTime())
{
    volScalarField magfSigmaCell = mag(fSigmaCell);
    magfSigmaCell.rename("magfSigmaCell");
    magfSigmaCell.write();
}

            //Info << "alphaEqn solution time : " 
            //    << runTime.cpuTimeIncrement() << endl;

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            //runTime.cpuTimeIncrement(); 
            #include "UEqn.hpp"
            //Info << "Ueqn time : " 
            //    << runTime.cpuTimeIncrement() << endl;

            //runTime.cpuTimeIncrement(); 

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.hpp"
            }

            //Info << "pEqn time : " 
            //    << runTime.cpuTimeIncrement() << endl;

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        //runTime.cpuTimeIncrement(); 
        runTime.write();
        //Info << "Write time : " 
        //    << runTime.cpuTimeIncrement() << endl;

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //