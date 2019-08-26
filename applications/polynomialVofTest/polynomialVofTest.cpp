/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Application
    polynomialVofTest

Description
    Test the classes and functions involved in the polynomial Vof
    initialization.

Author
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

#include "adaptiveTetCellRefinement.hpp"
#include "orientedPlane.hpp"

#include <vector>

using namespace Foam::PolynomialVof;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"

    const point ref_point = args.optionLookupOrDefault<point>("refpoint", point{0.5, 0.5, 0.5});
    const vector normal = args.optionLookupOrDefault<vector>("normal", vector{0, 0, 1});
    const label refinement_level = args.optionLookupOrDefault<label>("reflevel", 1);

    // Define plane as interface
    orientedPlane plane{ref_point, normal, 1.0};

    // Define unit tet in positive quadrant
    std::vector<point> points{point{0,0,0}, point{0,0,1}, point{0,1,0}, point{1,0,0}};
    std::vector<indexedTet> tets{indexedTet{0, 1, 2, 3}};
    std::vector<scalar> signed_distance{};

    for (const auto& vertex : points)
    {
        signed_distance.push_back(plane.signedDistance(vertex));
    }

    adaptiveTetCellRefinement tet_refiner{plane, points, signed_distance, tets, refinement_level};

    auto refined_tets = tet_refiner.resulting_tets();

    Info << "Number of tets after refinement: " << refined_tets.size() << endl;

    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
