#include "doctest.h"

#include "../BoxPlan.hpp"
#include "../Background.hpp"

// milestone-1.md task 1.4: box planning of conventions.md sec.5.

TEST_CASE("General box plan reproduces N1, N2 at tau_f, for physical (c=0)")
{
    // Pick N, N1, N2 and check that plugging tau_f back into the Background
    // formulas reproduces H(tau_f)*L(tau_f) = N1 and 1/(dx(tau_f)*m_r(tau_f))
    // = N2, i.e. the box plan actually delivers what it was asked for.
    const int N       = 512;
    const double N1   = 1.0;
    const double N2   = 4.0;
    const double a_inv = 2.0; // radiation domination
    const double c      = 0.0; // physical

    const auto plan = compute_general_box_plan(N, N1, N2, a_inv, c);

    CTauSchedule sched{};
    sched.c0 = c;
    Background bkg(a_inv, sched);

    const double H         = bkg.R_prime(plan.tau_f) /
                      (bkg.R(plan.tau_f) * bkg.R(plan.tau_f));
    const double m_r       = std::sqrt(bkg.lambda(plan.tau_f));
    const double L_physical = plan.L_tilde * bkg.R(plan.tau_f);
    const double dx_physical = plan.delta_x * bkg.R(plan.tau_f);

    CHECK(H * L_physical == doctest::Approx(N1).epsilon(1.0e-9));
    CHECK(1.0 / (dx_physical * m_r) == doctest::Approx(N2).epsilon(1.0e-9));
}

TEST_CASE("General box plan reproduces N1, N2 at tau_f, for fat (c=1)")
{
    const int N        = 800;
    const double N1     = 1.5;
    const double N2     = 6.0;
    const double a_inv  = 2.0;
    const double c       = 1.0; // fat

    const auto plan = compute_general_box_plan(N, N1, N2, a_inv, c);

    CTauSchedule sched{};
    sched.c0 = c;
    Background bkg(a_inv, sched);

    const double H = bkg.R_prime(plan.tau_f) /
                      (bkg.R(plan.tau_f) * bkg.R(plan.tau_f));
    const double m_r        = std::sqrt(bkg.lambda(plan.tau_f));
    const double L_physical  = plan.L_tilde * bkg.R(plan.tau_f);
    const double dx_physical = plan.delta_x * bkg.R(plan.tau_f);

    CHECK(H * L_physical == doctest::Approx(N1).epsilon(1.0e-9));
    CHECK(1.0 / (dx_physical * m_r) == doctest::Approx(N2).epsilon(1.0e-9));
}

TEST_CASE("pre_evolution_L_tilde matches conventions.md sec.7 by hand")
{
    const double L_tilde_main = 10.0;
    const double a_inv        = 2.0;
    const double c            = 1.0;
    const double tau_i        = 4.0;

    // L_init = L_main * 1/(a_inv-1) * (tau_i/tau0)^((1-c)/(a_inv-1))
    //        = 10 * 1 * (4)^0 = 10, since c=1 makes the exponent 0.
    CHECK(pre_evolution_L_tilde(L_tilde_main, a_inv, c, tau_i) ==
          doctest::Approx(10.0).epsilon(1.0e-12));
}

TEST_CASE("pre_evolution_L_tilde exponent behaves for c != 1")
{
    const double L_tilde_main = 10.0;
    const double a_inv        = 2.0;
    const double c            = 0.0;
    const double tau_i        = 4.0;

    // exponent = (1-0)/(2-1) = 1, so L_init = 10 * 1 * 4^1 = 40
    CHECK(pre_evolution_L_tilde(L_tilde_main, a_inv, c, tau_i) ==
          doctest::Approx(40.0).epsilon(1.0e-12));
}

TEST_CASE("Moore dynamic range formula matches its algebraic inverse")
{
    // N = N2 * N1 * gamma * exp(D/2) <=> D = 2 log(N/(N2 N1 gamma))
    const double N1    = 1.0;
    const double N2    = 5.0;
    const double gamma = 3.0;
    const double D     = 6.0;

    const int N = static_cast<int>(std::round(N2 * N1 * gamma * std::exp(D / 2.0)));

    CHECK(moore_max_log_dynamic_range(N, N1, N2, gamma) ==
          doctest::Approx(D).epsilon(1.0e-3));
}
