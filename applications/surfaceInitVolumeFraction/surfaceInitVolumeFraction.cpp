/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    surfaceInitVolumeFraction

Description
    Initialize a volume fraction field from a triangulated surface.
    
    Computation of volume fractions in interface cells is either performed by
    the SMCI or SMCA algorithm.

    This application implements the SMCI/A algorithm described in

    Reference
    \verbatim
        Tolle, T., Gründing, D., Bothe, D., & Marić, T. (2021).
        Computing volume fractions and signed distances from arbitrary surfaces
        on unstructured meshes.
        arXiv preprint arXiv:2101.08511.
    \endverbatim

\*---------------------------------------------------------------------------*/

// STD headers
#include <chrono>
#include <iomanip>

// OpenFOAM headers
#include "fvCFD.H"
#include "fvcAverage.H"
#include "pointFieldsFwd.H"
#include "pointMesh.H"
#include "surfaceInterpolate.H"
#include "triSurface.H"

// Argo headers
#include "volumeFractionCalculator.hpp"

// Headers for prototype
#include "AdaptiveTetCellRefinement.hpp"
#include "implicitSurfaces.hpp"
#include "refinementCriteria.hpp"

using namespace Foam::TriSurfaceImmersion;

template<class T>
T setOptionByPrecedence(dictionary& dict, const argList& args, const word keyword, T def)
{
    def = dict.getOrDefault<T>(keyword, def);
    args.readIfPresent<T>(keyword, def);
    dict.add(keyword, def, true);

    return def;
}

template<>
Switch setOptionByPrecedence(dictionary& dict, const argList& args, const word keyword, Switch def)
{
    def = dict.getOrDefault<Switch>(keyword, def);
    if (args.found(keyword))
    {
        def = true;
    }
    dict.add(keyword, def, true);

    return def;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read from dictionary if present
    IOdictionary initDict
                 {
                    IOobject{
                       "vofInitDict",
                       "system",
                       mesh.time(),
                       IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE
                    }
                 };

    // * * * * Configuration * * * *
    // Precedence: commandline option > dictionary value > default

    // Configure signed distance calculator
    auto& distDict = initDict.subDictOrAdd("distCalc");
    setOptionByPrecedence<word>(distDict, args, "surfaceType", "triSurface");
    setOptionByPrecedence<scalar>(distDict, args, "narrowBandWidth", 4.0);
    setOptionByPrecedence<scalar>(distDict, args, "bulkValue", 0.0);
    setOptionByPrecedence<fileName>(distDict, args, "surfaceFile", args.path() + "/surface.stl");
    
    // Configure volume fraction calculator
    auto fieldName = 
        setOptionByPrecedence<word>(initDict, args, "fieldName", "alpha.water");
    auto algName = 
        setOptionByPrecedence<word>(initDict, args, "algorithm", "SMCI");
    setOptionByPrecedence<label>(initDict, args, "refinementLevel", -1);
    setOptionByPrecedence<Switch>(initDict, args, "writeGeometry", false);
    // TODO (TT): move invert option to signed distance calculator. This 
    // automatically causes an inversion of volume fractions.
    auto invertVolumeFraction = 
        setOptionByPrecedence<Switch>(initDict, args, "invert", false);
    auto writeAllFields = 
        setOptionByPrecedence<Switch>(initDict, args, "writeAllFields", false);
    auto checkVolume =
        setOptionByPrecedence<Switch>(initDict, args, "checkVolume", false);

    Info<< "<------------------------------------------>"
        << "\nConfiguration:" << initDict
        << "<------------------------------------------>"
        << endl;
    
    // Initialization
    #include "createFields.hpp"

    // Compute the volume fraction field.
    auto ctime0 = std::chrono::steady_clock::now();
    auto vofCalcPtr = volumeFractionCalculator::New(initDict, mesh);
    vofCalcPtr->calcVolumeFraction(alpha);
    auto ctime1 = std::chrono::steady_clock::now();
    auto calcTime = 
        std::chrono::duration_cast<std::chrono::microseconds>(ctime1-ctime0).count();

    if (invertVolumeFraction)
    {
        alpha = 1.0 - alpha;
    }

    alpha.write(); 

    // Begin testing and debugging
    if (writeAllFields)
    {
        vofCalcPtr->writeFields();
    }

    if (checkVolume)
    {
        if (Pstream::myProcNo() == 0) // Only compute on master rank in parallel. Thanks to TM.
        {
            scalar Vsurf = vofCalcPtr->sigDistCalc().surfaceEnclosedVolume(); 

            scalar Valpha = gSum((mesh.V() * alpha)()); 
            scalar Evsurf = std::abs(Valpha - Vsurf) / Vsurf; 

            std::cout << std::setprecision(20) 
                << "Volume by volume fraction = " << Valpha << nl
                << "Volume of the surface mesh by divergence theorem (only closed surfaces!) = " << Vsurf << nl 
                << "Volume error from surface integral = " << Evsurf << nl;

            std::ofstream errorFile; 
            // TODO (TT): make file name algorithm dependent
            errorFile.open("vof-init-results-" + algName + ".csv"); 
            errorFile << "N_CELLS,"
                << "N_TRIANGLES,"
                << "VOLUME_FROM_VOLUME_FRACTION,"
                << "VOLUME_FROM_SURFACE_INTEGRAL,"
                << "VOLUME_ERROR_FROM_SURFACE_INTEGRAL,"
                << "CPU_TIME_MICROSECONDS," 
                << "MAX_REFINEMENT_LEVEL" << "\n"
                << mesh.nCells() << "," 
                << vofCalcPtr->sigDistCalc().nSurfaceElements() << "," 
                << std::setprecision(20) 
                << Valpha << "," 
                << Vsurf << "," 
                << Evsurf << ","
                << calcTime << ","
                << vofCalcPtr->maxRefinementLevel() << "\n";
        }
    }

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

//} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //
