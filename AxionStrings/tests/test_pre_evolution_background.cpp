#include "doctest.h"

#include "../PreEvolutionBackground.hpp"

// conventions.md sec.7: pre-evolution's R(t)=R0(t/t0), lambda(t)~(t/t0)^-2
// are chosen so that R*m_r and m_r/H are both time-independent. Verifies
// the from-scratch derivation in PreEvolutionBackground.hpp (see its
// header comment) actually delivers that, at several gamma_pre values and
// several tau_pre.

TEST_CASE("m_r/H is exactly gamma_pre, independent of tau_pre")
{
    for (const double gamma_pre : {0.5, 1.0, 3.0})
    {
        PreEvolutionBackground bkg(gamma_pre);
        for (const double tau_pre : {-2.0, -0.3, 0.0, 0.7, 4.0})
        {
            CHECK(bkg.m_r_over_H_direct(tau_pre) ==
                  doctest::Approx(gamma_pre).epsilon(1.0e-9));
        }
    }
}

TEST_CASE("R * m_r is time-independent")
{
    PreEvolutionBackground bkg(2.0);

    const double R_mr_ref =
        bkg.R(0.0) * std::sqrt(bkg.lambda(0.0));

    for (const double tau_pre : {-1.5, 0.0, 0.6, 2.2})
    {
        const double R_mr = bkg.R(tau_pre) * std::sqrt(bkg.lambda(tau_pre));
        CHECK(R_mr == doctest::Approx(R_mr_ref).epsilon(1.0e-9));
    }
}

TEST_CASE("R_prime matches a numerical derivative of R")
{
    PreEvolutionBackground bkg(1.7);
    const double dtau = 1.0e-6;

    for (const double tau_pre : {-1.0, 0.0, 1.3})
    {
        const double numeric =
            (bkg.R(tau_pre + dtau) - bkg.R(tau_pre - dtau)) / (2.0 * dtau);
        CHECK(bkg.R_prime(tau_pre) ==
              doctest::Approx(numeric).epsilon(1.0e-6));
    }
}

TEST_CASE("R'' = curvature_term_coeff * R (the free-EOM identity the main "
          "Background satisfies too)")
{
    PreEvolutionBackground bkg(1.0);
    const double dtau = 1.0e-5;

    for (const double tau_pre : {-1.0, 0.0, 1.3})
    {
        const double R_second_numeric =
            (bkg.R(tau_pre + dtau) - 2.0 * bkg.R(tau_pre) +
             bkg.R(tau_pre - dtau)) /
            (dtau * dtau);
        const double expected =
            bkg.curvature_term_coeff(tau_pre) * bkg.R(tau_pre);
        CHECK(R_second_numeric == doctest::Approx(expected).epsilon(1.0e-4));
    }
}
