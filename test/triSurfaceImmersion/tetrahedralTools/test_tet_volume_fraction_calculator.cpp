#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <quaternion.H>

#include <tetVolumeFractionCalculator.hpp>

template<typename ctype = double>
struct Tetrahedron
{
    std::vector<point> points;
    std::array<label, 4> labels;
    ctype volume;
};

template<typename ctype = double>
Tetrahedron<ctype> make_tetrahedron(ctype scale)
{
    return {
        {{
            {0.0,       0.0,       0.0},
            {scale,     0.0,       0.0},
            {0.0,       scale,     0.0},
            {scale*0.5, scale*0.5, scale},
        }},                             // points
        {0, 1, 2, 3},                   // labels
        1.0/3.0*(0.5*scale*scale)*scale // volume
    };
}

template<typename ctype = double>
bool compare_equal(ctype a, ctype b, ctype eps)
{
    using std::abs;
    return abs(a - b) < eps;
}

template<typename ctype>
void check_volume(const Tetrahedron<ctype>& tet)
{
    using namespace Foam::TriSurfaceImmersion;
    const auto v = tetVolumeFractionCalculator::volume(tet.labels, tet.points);
    const auto expected = tet.volume;

    if (!compare_equal(v, expected, expected*1e-8))
    {
        std::cout << "Computed vs expected volume: "
                  << v << " / " << expected << std::endl;
        throw std::runtime_error("volume mismatch detected");
    }
}

template<typename ctype, typename Points>
void rotate(Points& points, const point& axis, const ctype angle)
{
    Foam::quaternion rotation{axis, angle};

    std::for_each(points.begin(),
                  points.end(),
                  [&] (auto& p) { rotation.transform(p); });
}

template<typename Points>
void displace(Points& points, const point& origin)
{
    std::for_each(points.begin(),
                  points.end(),
                  [&] (auto& p) { p += origin; });
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
        check_volume(tet);

        std::cout << " -- testing displaced tet" << std::endl;
        auto tet_displaced = tet;
        displace(tet_displaced.points, origin);
        check_volume(tet_displaced);

        std::cout << " -- testing rotated tet" << std::endl;
        auto tet_rotated = tet;
        rotate(tet_rotated.points, rot_axis, M_PI/4.0);
        check_volume(tet_rotated);

        std::cout << " -- testing rotated/displaced tet" << std::endl;
        auto tet_rot_displaced = tet_rotated;
        displace(tet_rot_displaced.points, origin);
        check_volume(tet_rot_displaced);
    }

    return 0;
}
