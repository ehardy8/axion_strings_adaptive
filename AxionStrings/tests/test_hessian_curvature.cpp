#include "doctest.h"

#include "../HessianCurvature.hpp"

#include <cmath>

TEST_CASE("hessian_curvature: two linear level sets (a straight line) give "
         "kappa = 0")
{
    // psi1 = x, psi2 = y -- the curve {x=0, y=0} is the z-axis, exactly
    // straight. Both Hessians are identically zero (linear fields).
    const HessianCurvatureResult r = hessian_curvature(
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // grad psi1, grad psi2
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // Hess psi1
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0  // Hess psi2
    );
    REQUIRE(r.valid);
    CHECK(r.kappa == doctest::Approx(0.0).epsilon(1.0e-12));
}

TEST_CASE("hessian_curvature: a circle of radius R written as an "
         "intersection of two level sets gives kappa = 1/R exactly")
{
    // psi1 = y (pins the y=0 plane); psi2 = x^2 + (z-R)^2 - R^2 (a circle
    // of radius R centred at (0,0,R) in that plane, passing through the
    // origin, tangent to the z-axis there). Evaluated at the origin,
    // analytically:
    //   grad(psi1) = (0,1,0), grad(psi2) = (2x, 0, 2(z-R)) = (0,0,-2R)
    //   Hess(psi1) = 0
    //   Hess(psi2): d2/dx2 = 2, d2/dz2 = 2, everything else 0
    // Worked by hand before writing this test (see HessianCurvature.hpp's
    // header comment) -- gives kappa = 1/R independent of the arbitrary
    // choice of R, a genuine check of the formula, not a coincidence of
    // a particular radius.
    for (const double R : {1.0, 3.7, 0.25})
    {
        const double g1x = 0.0, g1y = 1.0, g1z = 0.0;
        const double g2x = 0.0, g2y = 0.0, g2z = -2.0 * R;

        const double h1xx = 0.0, h1yy = 0.0, h1zz = 0.0, h1xy = 0.0,
                    h1yz = 0.0, h1xz = 0.0;
        const double h2xx = 2.0, h2yy = 0.0, h2zz = 2.0, h2xy = 0.0,
                    h2yz = 0.0, h2xz = 0.0;

        const HessianCurvatureResult r =
            hessian_curvature(g1x, g1y, g1z, g2x, g2y, g2z, h1xx, h1yy,
                              h1zz, h1xy, h1yz, h1xz, h2xx, h2yy, h2zz,
                              h2xy, h2yz, h2xz);
        REQUIRE(r.valid);
        CHECK(r.kappa == doctest::Approx(1.0 / R).epsilon(1.0e-10));
    }
}

TEST_CASE("hessian_curvature: scale-invariant under (psi1,psi2) -> "
         "c*(psi1,psi2), matching psi = R*phi not mattering")
{
    // Same circle configuration as above (R=2), scaled by an arbitrary
    // constant c on both fields simultaneously (as multiplying by the
    // background R(tau) would do) -- kappa must come out identical.
    const double R = 2.0;
    const double g1x = 0.0, g1y = 1.0, g1z = 0.0;
    const double g2x = 0.0, g2y = 0.0, g2z = -2.0 * R;
    const double h1xx = 0.0, h1yy = 0.0, h1zz = 0.0, h1xy = 0.0, h1yz = 0.0,
                h1xz = 0.0;
    const double h2xx = 2.0, h2yy = 0.0, h2zz = 2.0, h2xy = 0.0, h2yz = 0.0,
                h2xz = 0.0;

    const HessianCurvatureResult base = hessian_curvature(
        g1x, g1y, g1z, g2x, g2y, g2z, h1xx, h1yy, h1zz, h1xy, h1yz, h1xz,
        h2xx, h2yy, h2zz, h2xy, h2yz, h2xz);

    const double c = 5.3;
    const HessianCurvatureResult scaled = hessian_curvature(
        c * g1x, c * g1y, c * g1z, c * g2x, c * g2y, c * g2z, c * h1xx,
        c * h1yy, c * h1zz, c * h1xy, c * h1yz, c * h1xz, c * h2xx,
        c * h2yy, c * h2zz, c * h2xy, c * h2yz, c * h2xz);

    REQUIRE(base.valid);
    REQUIRE(scaled.valid);
    CHECK(scaled.kappa == doctest::Approx(base.kappa).epsilon(1.0e-9));
}

TEST_CASE("hessian_curvature: parallel gradients (no genuine 1D curve) "
         "give valid = false")
{
    const HessianCurvatureResult r = hessian_curvature(
        1.0, 2.0, -1.0, 3.0, 6.0, -3.0, // grad psi2 = 3 * grad psi1
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    CHECK_FALSE(r.valid);
}
