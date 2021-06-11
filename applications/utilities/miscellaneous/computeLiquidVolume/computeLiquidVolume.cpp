/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Tomislav Maric and TU Darmstadt
     \\/     M anipulation  |
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

Description
    Obtain the volume enclosed between an stl surface and the x-y plane.
    The computation assumes a convex surface and assumes a graph representation
    of the interface.

Author
    Dirk Gründing
    gruending@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "fvCFD.H"

#include "triSurface.H"

#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
scalar xyPrism(const pointField& tri)
{
    pointField opts = tri;

    for (auto& p : opts)
    {
        p.z() = 0;
    }

    Foam::triFace face(0, 1, 2);
    return face.sweptVol(opts, tri);
}

// Main program:
int main(int argc, char* argv[])
{
    argList::addOption("stl",
        "path/to/stl/file.stl"
        "name of the stl file");

    argList args(argc, argv);

#include "createTime.H"

    // file name of stl is input of the utility
    word fileName("");
    if (args.found("stl"))
    {
        fileName = args.get<word>("stl");
    }
    else
    {
        FatalError << "No stl-file provided.";
    }

    std::cout << "Reading: " << fileName << std::endl;
    triSurface surfaceMesh(fileName);
    // Info << surfaceMesh;
    // Information on the triSurface cout format:
    // triSurface is first writing the points to cout and then a tuple for each
    // triangle is the stl. the first entry of the tuple is a triple that
    // contains the point ids of the triangle. The second entry is an integer
    // that gives the region number of the triangle in the OpenFOAM
    // PrimitivePatch as can be found in the labelledTri: "Triangle with
    // additional region number."

    scalar volume(0);
    for (auto& triangle : surfaceMesh)
    {
        const pointField triPoints = triangle.points(surfaceMesh.points());
        volume += xyPrism(triPoints);
    }
    std::cout << "Volume with respect to the x-y-plane: " << volume
              << std::endl;

    Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
