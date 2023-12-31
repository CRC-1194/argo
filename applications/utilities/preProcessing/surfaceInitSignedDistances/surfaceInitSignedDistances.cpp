/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 2020 AUTHOR,AFFILIATION
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
    surfaceInitSignedDistances

Description
    Compute signed distances from an oriented triangulated surface or
    a level set given by an implicit function.

    Relies on the implementation of the SMCI/A algorithm described in

    Reference
    \verbatim
        Tolle, T., Gründing, D., Bothe, D., & Marić, T. (2021).
        Computing volume fractions and signed distances from arbitrary surfaces
        on unstructured meshes.
        arXiv preprint arXiv:2101.08511.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "signedDistanceCalculator.hpp"

using namespace Foam::TriSurfaceImmersion;

template<class T>
T setOptionByPrecedence(
    dictionary& dict, const argList& args, const word keyword, T def)
{
    def = dict.getOrDefault<T>(keyword, def);
    args.readIfPresent<T>(keyword, def);
    dict.add(keyword, def, true);

    return def;
}

template<>
Switch setOptionByPrecedence(
    dictionary& dict, const argList& args, const word keyword, Switch def)
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

int main(int argc, char* argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // * * * * Configuration * * * *
    // Precedence: commandline option > dictionary value > default

    // Read from dictionary if present
    IOdictionary initDict(IOobject("signedDistanceInitDict",
        "system",
        mesh.time(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE));

    // Read options from dictionary
    auto fieldName = setOptionByPrecedence<word>(
        initDict, args, "fieldName", "signedDistance");
    setOptionByPrecedence<word>(initDict, args, "surfaceType", "triSurface");
    setOptionByPrecedence<fileName>(
        initDict, args, "surfaceFile", "/surface.stl");
    setOptionByPrecedence<scalar>(initDict, args, "narrowBandWidth", -1.0);
    setOptionByPrecedence(initDict, args, "bulkValue", 0.0);
    auto propagateInsideOutside = setOptionByPrecedence<Switch>(
        initDict, args, "propagateInsideOutside", false);
    auto invertInsideOutside =
        setOptionByPrecedence<Switch>(initDict, args, "invert", false);
    auto writeAllFields =
        setOptionByPrecedence<Switch>(initDict, args, "writeAllFields", false);

    // Print configuration
    Info<< "<------------------------------------------>"
        << "\nConfiguration:" << initDict
        << "<------------------------------------------>" << endl;

    // Initialization
    #include "createFields.hpp"

    autoPtr<signedDistanceCalculator> sigDistCalcPtr{
        signedDistanceCalculator::New(initDict, mesh)};

    if (propagateInsideOutside)
    {
        signedDistance = sigDistCalcPtr->cellSignedDist();
    }
    else
    {
        signedDistance = sigDistCalcPtr->cellSignedDist0();
    }

    if (invertInsideOutside)
    {
        signedDistance *= -1.0;
    }

    signedDistance.write();

    // TT: TODO: include testing of inside/outside computation here, analogue to
    // '-checkVolume' option of vof-init application
    if (writeAllFields)
    {
        sigDistCalcPtr->writeFields();
    }

    Info << nl;
    runTime.printExecutionTime(Info);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
