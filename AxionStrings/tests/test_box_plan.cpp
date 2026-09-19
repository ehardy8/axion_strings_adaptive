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

TEST_CASE("AMR box plan: N_base/dx_base consistency and divisibility check")
{
    const int N_effective = 4096;
    const int max_level   = 4;
    const double N1 = 1.2, N2 = 8.0, a_inv = 2.0, c = 0.0; // physical

    const auto plan = compute_amr_box_plan(N_effective, N1, N2, a_inv, c, max_level);
    CHECK(plan.N_base == N_effective / (1 << max_level));
    CHECK(plan.dx_base == doctest::Approx(plan.dx_finest * (1 << max_level)));
    CHECK(static_cast<int>(plan.log_add.size()) == max_level);

    // Not evenly divisible -> N_base is the -1 sentinel, not a silently
    // wrong integer (e.g. from C++ truncating division).
    const auto bad_plan = compute_amr_box_plan(4097, N1, N2, a_inv, c, max_level);
    CHECK(bad_plan.N_base == -1);
}

TEST_CASE("AMR box plan: level thresholds are spaced by the derived "
          "Delta log(m_r/H), matching conventions.md sec.11's ln(4) at a_inv=2")
{
    const int N_effective = 4096;
    const int max_level   = 5;
    const double N1 = 1.0, N2 = 4.0, a_inv = 2.0, c = 0.0;

    const auto plan = compute_amr_box_plan(N_effective, N1, N2, a_inv, c, max_level);
    const double expected_delta = std::log(4.0); // a_inv=2, c=0 special case
    for (std::size_t i = 1; i < plan.log_add.size(); ++i)
    {
        CHECK(plan.log_add[i] - plan.log_add[i - 1] ==
             doctest::Approx(expected_delta).epsilon(1.0e-9));
    }

    // The finest level's own threshold sits exactly one level-spacing
    // before tau_f's own log(m_r/H) (tau_f is where *all* max_level levels
    // are needed, by compute_general_box_plan's construction).
    CTauSchedule sched{};
    sched.c0 = c;
    Background bkg(a_inv, sched);
    const double log_mr_over_h_tau_f = -std::log(bkg.H_over_mr_closed_form(plan.tau_f));
    CHECK(plan.log_add.back() ==
          doctest::Approx(log_mr_over_h_tau_f - expected_delta).epsilon(1.0e-9));
}

TEST_CASE("AMR box plan: level 1's threshold is exactly where the base "
          "grid alone first reaches the N2 target")
{
    // The physically meaningful boundary condition the whole schedule is
    // built from: at tau = tau_from_log_mr_over_h(log_add[0]), the *base*
    // grid (dx_base, no refinement at all) should give exactly N2 -- i.e.
    // level 1 is not needed a moment before this, and is needed from here
    // on. Checked directly against Background, not against
    // compute_amr_box_plan's own formula (which would be circular).
    const int N_effective = 8192;
    const int max_level   = 3;
    const double N1 = 1.0, N2 = 6.0, a_inv = 2.5, c = 0.0;

    const auto plan = compute_amr_box_plan(N_effective, N1, N2, a_inv, c, max_level);

    CTauSchedule sched{};
    sched.c0 = c;
    Background bkg(a_inv, sched);
    const double tau_1 = bkg.tau_from_log_mr_over_h(plan.log_add[0]);

    const double R_1   = bkg.R(tau_1);
    const double m_r_1 = std::sqrt(bkg.lambda(tau_1));
    const double N2_base_only = 1.0 / (R_1 * plan.dx_base * m_r_1);

    CHECK(N2_base_only == doctest::Approx(N2).epsilon(1.0e-6));
}

TEST_CASE("Moore box plan reproduces N2 at the switch and the achievable "
          "dynamic range at tau_end")
{
    // A fat (c0=1) background up to the switch -- R(tau)*m_r(tau) is
    // constant through the whole fat phase (independently re-derived in
    // this test, not assumed), so evaluating R_switch/m_r_switch at the
    // switch is representative of the whole phase.
    const int N        = 4000;
    const double N1     = 1.2;
    const double N2     = 8.0;
    const double a_inv  = 2.0;
    const double b_inv  = a_inv - 1.0;
    const double c0      = 1.0; // fat

    CTauSchedule sched{};
    sched.c0 = c0;
    Background bkg(a_inv, sched);

    const double tau_switch = 30.0;
    const double R_switch   = bkg.R(tau_switch);
    const double m_r_switch = std::sqrt(bkg.lambda(tau_switch));
    const double gamma      = 1.0 / bkg.H_over_mr_direct(tau_switch);

    const auto plan = compute_moore_box_plan(N, N1, N2, gamma, a_inv, b_inv,
                                             R_switch, m_r_switch, tau_switch);

    // N2 is reproduced exactly at the switch by construction.
    CHECK(1.0 / (N2 * plan.dx * R_switch * m_r_switch) ==
          doctest::Approx(1.0).epsilon(1.0e-12));

    // R(tau)*m_r(tau) is constant through the fat phase -- check this
    // independently at a different tau within the phase, confirming dx
    // (and hence N2) would be identical had it been evaluated there.
    const double tau_other   = 12.0;
    CHECK(bkg.R(tau_other) * std::sqrt(bkg.lambda(tau_other)) ==
          doctest::Approx(R_switch * m_r_switch).epsilon(1.0e-9));

    // H(tau) does not depend on c (only lambda does) -- compute H at
    // tau_switch and tau_end directly from R/R' (no lambda involved) and
    // check the ratio matches the achievable dynamic range D exactly.
    const double D_direct = std::log(bkg.H(tau_switch) / bkg.H(plan.tau_end));
    CHECK(D_direct == doctest::Approx(plan.D).epsilon(1.0e-9));
}
