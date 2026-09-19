#include "doctest.h"

#include "../FlatBackground.hpp"

// Flat-space loop simulations (2026-09-19): R=1, R'=0, curvature=0
// identically, lambda=m_r^2 constant -- checked at several t, not just
// t=0, since nothing here should depend on t at all.

TEST_CASE("FlatBackground: R=1, R_prime=0, curvature=0 at every t")
{
    FlatBackground bkg(2.3);
    for (const double t : {0.0, 0.5, 1.0, 10.0, 100.0})
    {
        CHECK(bkg.R(t) == 1.0);
        CHECK(bkg.R_prime(t) == 0.0);
        CHECK(bkg.curvature_term_coeff(t) == 0.0);
    }
}

TEST_CASE("FlatBackground: lambda = m_r^2, constant in t")
{
    FlatBackground bkg(2.3);
    for (const double t : {0.0, 0.5, 1.0, 10.0, 100.0})
    {
        CHECK(bkg.lambda(t) == doctest::Approx(2.3 * 2.3));
    }
}

TEST_CASE("FlatBackground: default m_r=1, matching the project's lambda0=1 "
         "convention")
{
    FlatBackground bkg{};
    CHECK(bkg.m_r == 1.0);
    CHECK(bkg.lambda(0.0) == 1.0);
}
