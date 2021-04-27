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
    signedDistanceTestCalc

Description
    Given oriented, triangulated surface, compute the signed distances to a
    Cartesian point cloud in a bounding box of the surface.

Author
    Tobias Tolle
    tolle@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"

#include "legacyVtkOutput.H"
#include "triSurfaceDistCalc.hpp"

#include <chrono>
#include <iomanip>

using namespace Foam::TriSurfaceImmersion;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
pointField generate_points(const boundBox& bb, const scalar bb_scale, const int n)
{
    const point bb_centre{bb.centre()};
    const vector span{bb.span()};
    const point pmin{bb_centre - 0.5*bb_scale*span};
    const point pmax{bb_centre + 0.5*bb_scale*span};

    pointField test_points(n*n*n);
    auto d_by_n = (pmax - pmin)/(n-1);

    for (int i = 0; i != n; ++i)
    {
        for (int k = 0; k != n; ++k)
        {
            for (int l = 0; l != n; ++l)
            {
                test_points[i*n*n + k*n + l] = pmin +
                    vector{i*d_by_n.x(), k*d_by_n.y(), l*d_by_n.z()};
            }
        }
    }

    return test_points;
}

// Main program:
int main(int argc, char *argv[])
{
    #include "createOptions.hpp"
    #include "setRootCase.H"

    // Default options
    fileName surfaceFile = "surface.vtk";
    label npoints_per_dim = 10;
    scalar r_search = -1.0;
    scalar bb_scale_factor = 3.0;

    // Comand line args
    args.readIfPresent<fileName>("surfaceFile", surfaceFile);
    args.readIfPresent<label>("npoints", npoints_per_dim);
    args.readIfPresent<scalar>("rsearch", r_search);
    args.readIfPresent<scalar>("bbfactor", bb_scale_factor);

    // Initialize surface
    triSurface surface{surfaceFile};
    triSurfaceSearch search{surface};
    const auto& surface_bb = search.tree().bb();

    if (r_search < 0)
    {
        r_search = 0.5*bb_scale_factor*surface_bb.mag();
    };

    // New calculation using the class
    triSurfaceDistCalc sig_dist_calc{surface};
    const auto test_points = generate_points(surface_bb, bb_scale_factor, npoints_per_dim);
    const scalarField search_dist_sqr(test_points.size(), r_search*r_search);

    auto signed_distance = sig_dist_calc.signedDistance(test_points, search_dist_sqr, 1e15);

    // Write point sets and associated fields
    std::string case_path = args.rootPath() + "/" + args.globalCaseName();
    auto surface_name = surfaceFile.substr(0, surfaceFile.find("."));

    auto point_cloud_vtk_stream = vtk_stream<std::ofstream>(
        case_path + "/signed_distance_" + surface_name,
        test_points
    );
    write_to_vtk_stream(point_cloud_vtk_stream, signed_distance, "signed_distance");

    Info<< "End" << endl;

    return 0;
}
// ************************************************************************* //
