#include "doctest.h"

#include "../TangentField.hpp"

#include <cmath>

TEST_CASE("tangent_from_gradients: orthogonal unit gradients give the "
         "textbook cross product, unit-normalised")
{
    // grad(psi1) = x-hat, grad(psi2) = y-hat -> x-hat cross y-hat = z-hat.
    const TangentVector t =
        tangent_from_gradients(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    REQUIRE(t.valid);
    CHECK(t.x == doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(t.y == doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(t.z == doctest::Approx(1.0).epsilon(1.0e-12));
}

TEST_CASE("tangent_from_gradients: result is always unit-normalised and "
         "perpendicular to both input gradients")
{
    // An arbitrary, non-axis-aligned, non-orthogonal pair of gradients --
    // a genuine check of the formula, not a special case that happens to
    // simplify.
    const double g1x = 2.0, g1y = -1.0, g1z = 0.5;
    const double g2x = 0.3, g2y = 1.7, g2z = -2.1;
    const TangentVector t =
        tangent_from_gradients(g1x, g1y, g1z, g2x, g2y, g2z);
    REQUIRE(t.valid);

    const double mag = std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
    CHECK(mag == doctest::Approx(1.0).epsilon(1.0e-10));

    const double dot1 = t.x * g1x + t.y * g1y + t.z * g1z;
    const double dot2 = t.x * g2x + t.y * g2y + t.z * g2z;
    CHECK(dot1 == doctest::Approx(0.0).epsilon(1.0e-10));
    CHECK(dot2 == doctest::Approx(0.0).epsilon(1.0e-10));
}

TEST_CASE("tangent_from_gradients: parallel gradients give no defined "
         "direction (valid = false), not a spurious zero vector treated "
         "as real")
{
    // grad(psi2) = 3 * grad(psi1) -- the phase gradient direction is
    // degenerate here (no winding), so no tangent is defined.
    const TangentVector t =
        tangent_from_gradients(1.0, 2.0, -1.0, 3.0, 6.0, -3.0);
    CHECK_FALSE(t.valid);
}

TEST_CASE("tangent_from_gradients: a vanishing gradient gives valid = "
         "false rather than dividing by zero")
{
    const TangentVector t =
        tangent_from_gradients(0.0, 0.0, 0.0, 1.0, 2.0, 3.0);
    CHECK_FALSE(t.valid);
}
