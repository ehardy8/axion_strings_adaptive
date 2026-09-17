#include "doctest.h"

#include "../Velocity.hpp"

#include <cmath>

// milestone-1.md task 1.10 (velocities only): verifies the R-cancellation
// derivation in Velocity.hpp's header comment by directly comparing the
// literal psi/Pi formula (conventions.md sec.8) against the simplified,
// R-independent phi/phi_dot form, at several (R, tau, b_inv) points. If
// the derivation (in particular the beta = 1/b_inv reading, not stated
// explicitly in conventions.md) were wrong, these would not agree.

namespace
{
struct Point
{
    double psi1, psi2, Pi1, Pi2, R, tau, b_inv, m_r;
};
} // namespace

TEST_CASE("gamma_sq_v_sq: literal psi/Pi formula matches the simplified "
          "phi/phi_dot form (R cancels exactly)")
{
    const Point points[] = {
        {0.3, -0.2, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0},
        {0.3, -0.2, 0.05, 0.1, 2.0, 1.0, 1.0, 1.0},
        {0.3, -0.2, 0.05, 0.1, 2.0, 3.0, 1.0, 1.5},
        {-0.7, 0.4, -0.3, 0.2, 0.5, 2.0, 0.5, 0.8},
        {0.0, 0.9, 0.4, -0.1, 3.0, 4.0, 1.0, 1.0},
        {1.1, -1.1, 0.02, -0.02, 1.5, 0.7, 0.5, 2.0},
    };

    for (const auto &p : points)
    {
        // phi = psi/R, phi_dot = (Pi - psi/(b_inv tau))/R^2 (v=1).
        const double phi1     = p.psi1 / p.R;
        const double phi2     = p.psi2 / p.R;
        const double inv_bt   = 1.0 / (p.b_inv * p.tau);
        const double phi_dot1 = (p.Pi1 - inv_bt * p.psi1) / (p.R * p.R);
        const double phi_dot2 = (p.Pi2 - inv_bt * p.psi2) / (p.R * p.R);

        const double from_phi =
            gamma_sq_v_sq_from_phi(phi1, phi2, phi_dot1, phi_dot2, p.m_r);
        const double from_psi = gamma_sq_v_sq_from_psi_literal(
            p.psi1, p.psi2, p.Pi1, p.Pi2, p.R, p.tau, p.b_inv, p.m_r);

        CAPTURE(p.R);
        CAPTURE(p.tau);
        CAPTURE(p.b_inv);
        CHECK(from_phi == doctest::Approx(from_psi).epsilon(1.0e-9));
    }
}

TEST_CASE("gamma_sq_v_sq_from_phi is zero at rest (phi_dot = 0)")
{
    CHECK(gamma_sq_v_sq_from_phi(0.5, -0.3, 0.0, 0.0, 1.2) == 0.0);
}

TEST_CASE("gamma_sq_v_sq_from_phi is non-negative for the global string's "
          "e1 < 0 (checked over a range of phi, phi_dot)")
{
    for (double phi1 = -1.0; phi1 <= 1.0; phi1 += 0.3)
    {
        for (double phi_dot1 = -1.0; phi_dot1 <= 1.0; phi_dot1 += 0.3)
        {
            const double v = gamma_sq_v_sq_from_phi(phi1, 0.1, phi_dot1, 0.05, 1.0);
            CAPTURE(phi1);
            CAPTURE(phi_dot1);
            CHECK(v >= -1.0e-12);
        }
    }
}
