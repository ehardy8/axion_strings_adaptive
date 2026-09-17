#ifndef BOXPLAN_HPP_
#define BOXPLAN_HPP_

// Box planning of docs/conventions.md sec.5: turning the two physical
// resolution targets (N1 = H(tau_f) L(tau_f), N2 = 1/(Delta(tau_f) m_r(tau_f)))
// specified at the final time into a comoving box size and stop time.
//
// AMReX/ParmParse-free like Background.hpp, for the same reason: pure host
// math, unit-testable standalone (see tests/test_box_plan.cpp).
//
// v = 1, lambda0 = 1, tau0 = 1 throughout (conventions.md sec.5 code units).

#include <cmath>

// General (non-Moore) box plan:
//   tau_f    = (N/(N1 N2))^((a_inv-1)/(a_inv-c)) tau0
//   L_tilde  = (N1/(R0 sqrt(lambda0))) (N/(N1 N2))^((a_inv-1)/(a_inv-c))
//
// Only valid for c != a_inv -- the exponent is singular exactly at c = a_inv
// (Moore), which is why Moore has its own box-planning formula below
// (conventions.md sec.5, "Moore-phase box planning"): during Moore, HL
// evolves as 1/tau rather than approaching a fixed target the same way, so
// this formula does not apply and must not be used to plan a Moore phase.
struct GeneralBoxPlan
{
    double tau_f{};
    double L_tilde{};
    double delta_x{};
    // dx / r_delta, r_delta = 3 (conventions.md sec.5). Informational only:
    // this bound is inherited from leapfrog and flagged there as "to be
    // re-derived for the RK integrator" -- not enforced.
    double delta_tau_leapfrog_bound{};
};

[[nodiscard]] inline GeneralBoxPlan
compute_general_box_plan(int N, double N1, double N2, double a_inv, double c)
{
    const double exponent = (a_inv - 1.0) / (a_inv - c);
    const double ratio    = static_cast<double>(N) / (N1 * N2);
    const double power    = std::pow(ratio, exponent);
    const double R0       = 1.0 / (a_inv - 1.0); // sec.5, lambda0 = 1

    GeneralBoxPlan plan{};
    plan.tau_f                     = power; // * tau0 = 1
    plan.L_tilde                   = (N1 / R0) * power;
    plan.delta_x                   = plan.L_tilde / N;
    plan.delta_tau_leapfrog_bound  = plan.delta_x / 3.0;
    return plan;
}

// Comoving box size for the pre-evolution stage (conventions.md sec.7):
//   L_tilde_init = L_tilde_main * (1/(a_inv-1)) * (tau_i/tau0)^((1-c)/(a_inv-1))
// c is the main run's c at its start (c0), tau0 = 1.
[[nodiscard]] inline double pre_evolution_L_tilde(double L_tilde_main,
                                                  double a_inv, double c,
                                                  double tau_i)
{
    constexpr double tau0 = 1.0;
    return L_tilde_main * (1.0 / (a_inv - 1.0)) *
           std::pow(tau_i / tau0, (1.0 - c) / (a_inv - 1.0));
}

// Moore-phase dynamic-range check (conventions.md sec.5): during c = 1+b_inv,
// HL falls as 1/tau while N2 only improves; the phase ends when HL reaches
// N1. gamma = m_r/H at the start of the Moore phase (fixed by the preceding
// fat phase, which holds N2 fixed from the very start -- conventions.md
// sec.5 point 1). Returns the maximum achievable log(H_start/H_end).
[[nodiscard]] inline double moore_max_log_dynamic_range(int N, double N1,
                                                        double N2,
                                                        double gamma)
{
    return 2.0 * std::log(static_cast<double>(N) / (N2 * N1 * gamma));
}

#endif // BOXPLAN_HPP_
