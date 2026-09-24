#include "doctest.h"

#include "../PlaquetteCrossing.hpp"

#include <cmath>

TEST_CASE("find_bilinear_zero: exact linear fields recover their exact "
         "interior zero")
{
    // f(u,v) = u - u0, g(u,v) = v - v0 -- already exactly bilinear (no uv
    // term), so the corner values are just u-u0/v-v0 evaluated at the 4
    // unit-square corners, and the solver should recover (u0,v0) to
    // floating-point precision, not just approximately.
    const double u0 = 0.37, v0 = 0.62;
    const double f00 = 0.0 - u0, f10 = 1.0 - u0, f11 = 1.0 - u0,
                f01 = 0.0 - u0;
    const double g00 = 0.0 - v0, g10 = 0.0 - v0, g11 = 1.0 - v0,
                g01 = 1.0 - v0;

    const BilinearZero z =
        find_bilinear_zero(f00, f10, f11, f01, g00, g10, g11, g01);
    REQUIRE(z.valid);
    CHECK(z.u == doctest::Approx(u0).epsilon(1.0e-10));
    CHECK(z.v == doctest::Approx(v0).epsilon(1.0e-10));
}

TEST_CASE("find_bilinear_zero: a genuine single-winding vortex plaquette "
         "gives one interior zero")
{
    // psi = rho * e^{i theta}, theta winding once around the square, zero
    // placed off-centre (not at u=v=0.5, so this is a real check of the
    // interpolation, not a symmetric special case). rho > 0 everywhere on
    // the corners themselves (the corners are never exactly at the core).
    const double u0 = 0.3, v0 = 0.7;
    // A genuine bilinear (uv) twist on top of the linear zero, so this
    // isn't just the same case as the pure-linear test above.
    auto f = [&](double u, double v)
    { return (u - u0) + 0.15 * (u - u0) * (v - v0); };
    auto g = [&](double u, double v)
    { return (v - v0) - 0.1 * (u - u0) * (v - v0); };

    const double f00 = f(0, 0), f10 = f(1, 0), f11 = f(1, 1), f01 = f(0, 1);
    const double g00 = g(0, 0), g10 = g(1, 0), g11 = g(1, 1), g01 = g(0, 1);

    const BilinearZero z =
        find_bilinear_zero(f00, f10, f11, f01, g00, g10, g11, g01);
    REQUIRE(z.valid);
    CHECK(z.u == doctest::Approx(u0).epsilon(1.0e-8));
    CHECK(z.v == doctest::Approx(v0).epsilon(1.0e-8));
}

TEST_CASE("find_bilinear_zero: an unpierced plaquette (no interior zero) "
         "reports invalid rather than an out-of-range or fabricated "
         "answer")
{
    // f,g both strictly positive at all 4 corners -- no zero crossing
    // possible anywhere, let alone inside the unit square.
    const BilinearZero z =
        find_bilinear_zero(1.0, 1.2, 1.1, 0.9, 2.0, 2.1, 1.9, 2.2);
    CHECK_FALSE(z.valid);
}

TEST_CASE("find_bilinear_zero: zero sitting exactly on a corner is still "
         "recovered at the boundary of the unit square")
{
    const double f00 = 0.0, f10 = 1.0, f11 = 1.0, f01 = 0.0; // f = u
    const double g00 = 0.0, g10 = 0.0, g11 = 1.0, g01 = 1.0; // g = v
    const BilinearZero z =
        find_bilinear_zero(f00, f10, f11, f01, g00, g10, g11, g01);
    REQUIRE(z.valid);
    CHECK(z.u == doctest::Approx(0.0).epsilon(1.0e-10));
    CHECK(z.v == doctest::Approx(0.0).epsilon(1.0e-10));
}
