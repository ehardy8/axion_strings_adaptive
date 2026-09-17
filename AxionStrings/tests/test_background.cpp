#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../Background.hpp"

// milestone-1.md task 1.3, "done when": the three modes (c = 0, 1, 1+b_inv)
// all run, and m_r/H measured from the code matches the analytic expression
// in conventions.md sec.5 to round-off.
//
// H_over_mr_direct is computed from R(tau), R'(tau) and lambda(tau)
// (the quantities the RHS actually evaluates); H_over_mr_closed_form is the
// sec.5 shortcut (tau/tau0)^((c-a_inv)/(a_inv-1)). They are independent code
// paths -- agreement is a genuine check of the algebra, not circular.

namespace
{
constexpr double tol = 1.0e-12;

void check_closed_form_matches_direct(double a_inv, double c)
{
    CTauSchedule sched{};
    sched.c0 = c;
    Background bkg(a_inv, sched);

    for (const double tau : {0.3, 1.0, 2.7, 15.4})
    {
        const double direct      = bkg.H_over_mr_direct(tau);
        const double closed_form = bkg.H_over_mr_closed_form(tau);
        CAPTURE(a_inv);
        CAPTURE(c);
        CAPTURE(tau);
        CHECK(direct == doctest::Approx(closed_form).epsilon(tol));
    }
}
} // namespace

TEST_CASE("H/m_r closed form matches direct computation: physical (c=0), RD")
{
    // a_inv = 2, b_inv = 1 in radiation domination (conventions.md sec.2)
    check_closed_form_matches_direct(2.0, 0.0);
}

TEST_CASE("H/m_r closed form matches direct computation: fat string (c=1), RD")
{
    check_closed_form_matches_direct(2.0, 1.0);
}

TEST_CASE("H/m_r closed form matches direct computation: Moore (c=1+b_inv), RD")
{
    check_closed_form_matches_direct(2.0, 2.0); // b_inv = 1 in RD
}

TEST_CASE("R(tau) satisfies the power law and its own derivative")
{
    CTauSchedule sched{};
    sched.c0 = 0.0;
    Background bkg(2.0, sched);

    const double dtau = 1.0e-6;
    for (const double tau : {0.5, 1.0, 3.3})
    {
        const double numeric_deriv =
            (bkg.R(tau + dtau) - bkg.R(tau - dtau)) / (2.0 * dtau);
        CHECK(bkg.R_prime(tau) ==
              doctest::Approx(numeric_deriv).epsilon(1.0e-6));
    }
}

TEST_CASE("Homogeneous psi = R(tau) exactly solves the free EOM: "
          "R'' = curvature_term_coeff * R")
{
    // Cross-check used by AxionStringsLevel::initData: with psi1 = R(tau),
    // psi2 = 0, the potential term (|psi|^2 - R^2) vanishes identically, so
    // this configuration is an exact solution of the full nonlinear EOM iff
    // R'' = curvature_term_coeff(tau) * R(tau) for the analytic R(tau).
    CTauSchedule sched{};
    sched.c0 = 0.0;
    Background bkg(2.0, sched);

    const double dtau = 1.0e-5;
    for (const double tau : {0.5, 1.0, 3.3, 10.0})
    {
        const double R_second_numeric =
            (bkg.R(tau + dtau) - 2.0 * bkg.R(tau) + bkg.R(tau - dtau)) /
            (dtau * dtau);
        const double expected = bkg.curvature_term_coeff(tau) * bkg.R(tau);
        CHECK(R_second_numeric == doctest::Approx(expected).epsilon(1.0e-4));
    }
}

TEST_CASE("lambda(tau) is continuous across a c-switch")
{
    CTauSchedule sched{};
    sched.c0         = 1.0;
    sched.c1         = 2.0;
    sched.tau_switch = 5.0;
    sched.has_switch = true;
    Background bkg(2.0, sched);

    const double eps        = 1.0e-9;
    const double lambda_lhs = bkg.lambda(sched.tau_switch - eps);
    const double lambda_rhs = bkg.lambda(sched.tau_switch + eps);
    CHECK(lambda_lhs == doctest::Approx(lambda_rhs).epsilon(1.0e-6));
}

TEST_CASE("lambda(tau) reduces to the no-switch form before/without a switch")
{
    CTauSchedule sched{};
    sched.c0 = 1.5;
    Background bkg(2.0, sched);

    for (const double tau : {0.5, 1.0, 3.3})
    {
        CHECK(bkg.lambda(tau) ==
              doctest::Approx(std::pow(bkg.R(tau) / bkg.R0, -2.0 * sched.c0))
                  .epsilon(tol));
    }
}
