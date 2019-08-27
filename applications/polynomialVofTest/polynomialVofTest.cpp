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

#include "AdaptiveTetCellRefinement.hpp"
#include "orientedPlane.hpp"
#include "tetVofCalculator.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace Foam::PolynomialVof;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void save_to_vtk
(
    const std::vector<indexedTet>& tets,
    const std::vector<point>& points,
    const std::vector<scalar>& signed_distance,
    const std::vector<scalar>& vof,
    std::string file_name = "decomposed_tet.vtk"
)
{
    // Use VTK legacy format for now for the sak of simplicity (TT)
    std::ofstream out_file;
    out_file.open(file_name);

    // Header
    out_file << "# vtk DataFile Version 3.0\n";

    out_file << "Decomposition of a single tetragedron.\n";
    out_file << "ASCII\n";
    out_file << "DATASET UNSTRUCTURED_GRID\n";

    // Write points
    out_file << "POINTS " << std::to_string(points.size()) << " double\n";
    for (const auto& p : points)
    {
        out_file << std::to_string(p[0]) << " "
                 << std::to_string(p[1]) << " "
                 << std::to_string(p[2]) << "\n";
    }

    // Write tets
    out_file << "CELLS " << std::to_string(tets.size()) << " "
             << std::to_string(5*tets.size()) << "\n";
    for (const auto& tet : tets)
    {
        out_file << "4 "
                 << std::to_string(tet[0]) << " "
                 << std::to_string(tet[1]) << " "
                 << std::to_string(tet[2]) << " "
                 << std::to_string(tet[3]) << "\n";
    }
    out_file << "CELL_TYPES " << std::to_string(tets.size()) << "\n";
    for (uint idx = 0; idx != tets.size(); ++idx)
    {
        out_file << "10\n";
    }

    // Write signed distance
    out_file << "POINT_DATA " << std::to_string(points.size()) << "\n";
    out_file << "Scalars signed_distance double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (const auto phi : signed_distance)
    {
        out_file << std::to_string(phi) << "\n";
    }

    // Write volume fraction as cell data
    out_file << "CELL_DATA " << std::to_string(tets.size()) << "\n";
    out_file << "Scalars vof double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (const auto alpha : vof)
    {
        out_file << std::to_string(alpha) << "\n";
    }

    out_file.close();
}

scalar tet_volume(const indexedTet& t, const std::vector<point>& p)
{
    return mag((p[t[1]] - p[t[0]]) & ((p[t[2]] - p[t[0]])^(p[t[3]] - p[t[0]])))/6.0;
}

// Main program:
int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"

    const point ref_point = args.getOrDefault<point>("refpoint", point{0.5, 0.5, 0.3});
    vector normal = args.getOrDefault<vector>("normal", vector{0, 0, 1});
    const label refinement_level = args.getOrDefault<label>("reflevel", 1);

    // Define plane as interface
    if (mag(normal) > SMALL)
    {
        normal = normal/(mag(normal));
    }
    else
    {
        exit(1);
    }
    orientedPlane plane{ref_point, normal, 1.0};

    // Define unit tet in positive quadrant
    std::vector<point> points{point{0,0,0}, point{0,0,1}, point{0,1,0}, point{1,0,0}};
    std::vector<indexedTet> tets{indexedTet{0, 1, 2, 3}};
    std::vector<scalar> signed_distance{};

    for (const auto& vertex : points)
    {
        signed_distance.push_back(plane.signedDistance(vertex));
    }

    adaptiveTetCellRefinement<orientedPlane> tet_refiner{plane, points, signed_distance, tets, refinement_level};

    auto refined_tets = tet_refiner.resulting_tets();

    tet_refiner.print_level_infos();


    Info << "Volume plausibility check" << endl;
    auto vol_exact = tet_volume(tets[0], points);
    scalar vol_added = 0;

    for (const auto& tet : refined_tets)
    {
        vol_added += tet_volume(tet, tet_refiner.points());
    }

    Info << "Volume exact: " << vol_exact << "; volume added up = "
         << vol_added << "\n"
         << "Relative difference: " << mag(vol_added - vol_exact)/vol_exact << endl;

    // VoF calculation
    tetVofCalculator vof_calc{};

    auto volume_fractions = vof_calc.vof(refined_tets, tet_refiner.signed_distance());

    auto plus_volume = vof_calc.accumulated_omega_plus_volume(refined_tets, tet_refiner.signed_distance(), tet_refiner.points());

    Info << "Plus volume fraction = " << plus_volume/vol_exact << endl;

    save_to_vtk(refined_tets, tet_refiner.points(), tet_refiner.signed_distance(), volume_fractions);

    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
