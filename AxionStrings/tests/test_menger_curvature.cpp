#include "doctest.h"

#include "../MengerCurvature.hpp"

#include <cmath>
#include <numbers>

namespace
{
constexpr double pi = std::numbers::pi;
}

TEST_CASE("menger_curvature: three collinear points give kappa = 0")
{
    CHECK(menger_curvature(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.5, 0.0, 0.0) ==
         doctest::Approx(0.0).epsilon(1.0e-12));
    CHECK(menger_curvature(0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 2.5, 5.0, 7.5) ==
         doctest::Approx(0.0).epsilon(1.0e-12));
}

TEST_CASE("menger_curvature: three points on a circle of radius R give "
         "kappa = 1/R")
{
    const double R = 3.7;
    const double theta1 = 0.0;
    const double theta2 = 2.0 * pi / 5.0;
    const double theta3 = 4.0 * pi / 5.0;
    const double x1 = R * std::cos(theta1), y1 = R * std::sin(theta1);
    const double x2 = R * std::cos(theta2), y2 = R * std::sin(theta2);
    const double x3 = R * std::cos(theta3), y3 = R * std::sin(theta3);

    CHECK(menger_curvature(x1, y1, 0.0, x2, y2, 0.0, x3, y3, 0.0) ==
         doctest::Approx(1.0 / R).epsilon(1.0e-10));
}

TEST_CASE("menger_curvature: is invariant under which point is treated as "
         "the middle one being tested from either side")
{
    const double kappa_forward =
        menger_curvature(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 0.0, 0.0);
    const double kappa_reversed =
        menger_curvature(2.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    CHECK(kappa_forward == doctest::Approx(kappa_reversed).epsilon(1.0e-12));
    CHECK(kappa_forward > 0.0);
}

TEST_CASE("menger_curvature: degenerate (coincident) points give 0, not NaN "
         "or a division by zero")
{
    CHECK(menger_curvature(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0) ==
         0.0);
}
