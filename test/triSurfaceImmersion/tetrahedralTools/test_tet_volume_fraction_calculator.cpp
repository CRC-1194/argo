#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <quaternion.H>

#include <tetVolumeFractionCalculator.hpp>

using Points = std::vector<point>;
using Labels = std::array<label, 4>;

template<typename ctype = double>
struct Tetrahedron
{
    Points points;
    Labels labels;
    ctype volume;
    ctype height;
};

template<typename ctype>
Tetrahedron<ctype> make_tetrahedron(ctype scale)
{
    return {
        {{
            {0.0,       0.0,       0.0},
            {scale,     0.0,       0.0},
            {0.0,       scale,     0.0},
            {scale*0.5, scale*0.5, scale},
        }},                              // points
        {0, 1, 2, 3},                    // labels
        1.0/3.0*(0.5*scale*scale)*scale, // volume
        scale                            // height
    };
}

template<typename ctype, typename Points>
void rotate(Points& points, const point& axis, const ctype angle)
{
    Foam::quaternion rotation{axis, angle};

    std::for_each(points.begin(),
                  points.end(),
                  [&] (auto& p) { rotation.transform(p); });
}

template<typename ctype>
std::vector<ctype> get_signed_dists_full(const Tetrahedron<ctype>& tet)
{
    std::vector<ctype> signed_distances(tet.points.size(), tet.height);
    signed_distances[tet.points.size() - 1] = 0.0;
    return signed_distances;
}

template<typename ctype>
std::vector<ctype> get_signed_dists_empty(const Tetrahedron<ctype>& tet)
{
    std::vector<ctype> signed_distances(tet.points.size(), 0.0);
    signed_distances[tet.points.size() - 1] = -1.0*tet.height;
    return signed_distances;
}

template<typename ctype>
std::vector<ctype> get_signed_dists_halved(const Tetrahedron<ctype>& tet)
{
    const auto halfHeight = tet.height/2.0;
    std::vector<ctype> signed_distances(tet.points.size(), halfHeight);
    signed_distances[tet.points.size() - 1] = -1.0*halfHeight;
    return signed_distances;
}

template<typename Points>
void displace(Points& points, const point& origin)
{
    std::for_each(points.begin(),
                  points.end(),
                  [&] (auto& p) { p += origin; });
}

template<typename ctype>
bool compare_equal(ctype a, ctype b, ctype eps)
{
    using std::abs;
    return abs(a - b) < eps;
}

template<typename ctype>
void check_equality(ctype computed,
                    ctype expected,
                    ctype eps,
                    const std::string& quantityName)
{
    if (!compare_equal(computed, expected, eps))
    {
        std::cout << "Computed vs expected " << quantityName << ": "
                  << computed << " / " << expected << std::endl;
        throw std::runtime_error(quantityName + " mismatch detected");
    }
}

template<typename ctype>
void check_volume(const Tetrahedron<ctype>& tet)
{
    using namespace Foam::TriSurfaceImmersion;
    const auto v = tetVolumeFractionCalculator::volume(tet.labels, tet.points);
    const auto expected = tet.volume;
    check_equality(v, expected, expected*1e-8, "volume");
}

template<typename ctype>
void check_volume_fraction(const Tetrahedron<ctype>& tet,
                           const std::vector<ctype>& signed_distances,
                           const ctype expected,
                           const ctype eps)
{
    using namespace Foam::TriSurfaceImmersion;
    const tetVolumeFractionCalculator calculator{};
    const auto f = calculator.volumeFraction(tet.labels, signed_distances);
    check_equality(f, expected, eps, "volume fraction");

    // test function taking multiple tets
    std::vector<Labels> tet_list{tet.labels, tet.labels, tet.labels};
    const auto fractions = calculator.volumeFractions(tet_list, signed_distances);
    std::for_each(
        fractions.begin(),
        fractions.end(),
        [&] (const auto& f) { return check_equality(f, expected, eps, "volume fractions"); }
    );
}

template<typename ctype>
void check_omega_plus_volume(const Tetrahedron<ctype>& tet,
                             const std::vector<ctype>& signed_distances,
                             const ctype expected,
                             const ctype eps)
{
    using namespace Foam::TriSurfaceImmersion;
    const tetVolumeFractionCalculator calculator{};
    const auto v = calculator.omegaPlusVolume(tet.labels, signed_distances, tet.points);
    check_equality(v, expected, eps, "omegaPlusVolume");

    // test function taking multiple tets
    std::vector<Labels> tet_list{tet.labels, tet.labels, tet.labels};
    const auto volumes = calculator.omegaPlusVolumes(tet_list, signed_distances, tet.points);
    std::for_each(
        volumes.begin(),
        volumes.end(),
        [&] (const auto& v) { return check_equality(v, expected, eps, "omegaPlusVolumes"); }
    );
}

template<typename ctype>
void do_checks(const Tetrahedron<ctype>& tet)
{
    check_volume(tet);

    check_volume_fraction(tet, get_signed_dists_empty(tet),  0.0,   1e-8);
    check_volume_fraction(tet, get_signed_dists_full(tet),   1.0,   1e-8);
    check_volume_fraction(tet, get_signed_dists_halved(tet), 0.875, 1e-8);

    const auto v = tet.volume;
    check_omega_plus_volume(tet, get_signed_dists_empty(tet),  0.0,     1e-8);
    check_omega_plus_volume(tet, get_signed_dists_full(tet),   v,       v*1e-8);
    check_omega_plus_volume(tet, get_signed_dists_halved(tet), 0.875*v, v*1e-8);
}

int main()
{
    using ctype = double;

    for (const auto sf : std::vector<ctype>({1e-5, 1e-2, 1, 1e2, 1e5}))
    {
        std::cout << "Testing with scale factor: " << sf << std::endl;

        const point rot_axis{sf, sf, sf};
        const point origin{sf*10.0, sf*10.0, sf*10.0};

        std::cout << " -- testing raw tet" << std::endl;
        const auto tet = make_tetrahedron(sf);
        do_checks(tet);

        std::cout << " -- testing displaced tet" << std::endl;
        auto tet_displaced = tet;
        displace(tet_displaced.points, origin);
        do_checks(tet_displaced);

        std::cout << " -- testing rotated tet" << std::endl;
        auto tet_rotated = tet;
        rotate(tet_rotated.points, rot_axis, M_PI/4.0);
        do_checks(tet_rotated);

        std::cout << " -- testing rotated/displaced tet" << std::endl;
        auto tet_rot_displaced = tet_rotated;
        displace(tet_rot_displaced.points, origin);
        do_checks(tet_rot_displaced);

    }

    return 0;
}
