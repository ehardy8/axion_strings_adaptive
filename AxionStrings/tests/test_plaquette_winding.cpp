#include "doctest.h"

#include "../PlaquetteWinding.hpp"

#include <numbers>

namespace
{
constexpr double pi = std::numbers::pi;
}

TEST_CASE("plaquette_winding detects a +1 winding")
{
    CHECK(plaquette_winding(0.0, pi / 2, pi, 3 * pi / 2) == 1);
}

TEST_CASE("plaquette_winding detects a -1 winding")
{
    CHECK(plaquette_winding(0.0, -pi / 2, -pi, -3 * pi / 2) == -1);
}

TEST_CASE("plaquette_winding is zero when the phase barely varies")
{
    CHECK(plaquette_winding(0.1, 0.2, 0.15, 0.05) == 0);
    CHECK(plaquette_winding(0.0, 0.0, 0.0, 0.0) == 0);
}

TEST_CASE("plaquette_winding is invariant under an overall phase shift")
{
    const double shift = 1.234;
    CHECK(plaquette_winding(shift, shift + pi / 2, shift + pi,
                            shift + 3 * pi / 2) == 1);
}

TEST_CASE("plaquette_winding is zero for phases oscillating across the "
          "+-pi branch cut with no net winding")
{
    CHECK(plaquette_winding(3.0, -3.0, 3.0, -3.0) == 0);
}

TEST_CASE("plaquette_winding detects winding 2 when actually present")
{
    CHECK(plaquette_winding(0.0, pi, 2 * pi, 3 * pi) == 2);
}

TEST_CASE("plaquette_winding does not depend on which corner is listed first"
          "(cyclic invariance)")
{
    const int w = plaquette_winding(0.0, pi / 2, pi, 3 * pi / 2);
    const int w_rotated =
        plaquette_winding(pi / 2, pi, 3 * pi / 2, 2 * pi + 0.0);
    CHECK(w == w_rotated);
}
