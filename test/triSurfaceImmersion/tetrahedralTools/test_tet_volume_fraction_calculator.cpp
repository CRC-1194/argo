#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <tetVolumeFractionCalculator.hpp>

template<typename ctype = double>
struct Tetrahedron
{
    std::vector<point> points;
    std::array<label, 4> labels;
    ctype volume;
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
        }},                             // points
        {0, 1, 2, 3},                   // labels
        1.0/3.0*(0.5*scale*scale)*scale // volume
    };
}

template<typename ctype>
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
        std::string msg = "volume mismatch detected (computed/expected): ";
        msg += std::to_string(v);
        msg += " / ";
        msg += std::to_string(expected);
        throw std::runtime_error(msg);
    }
}

int main()
{
    using ctype = double;

    for (const auto sf : std::vector<ctype>({1e-5, 1e-2, 1, 1e2, 1e5}))
    {
        std::cout << "Testing with scale factor: " << sf << std::endl;
        check_volume(make_tetrahedron(sf));
    }

    return 0;
}
