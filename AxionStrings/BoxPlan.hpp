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

// Full Moore-phase box plan (2026-09-19, with the user: automating what had
// been "geometry.prob_extent and evolution.stop_time are not derived in
// Moore mode -- set them by hand"). N (box size) and N1/N2 (Hubble-patch/
// core-resolution targets) are the inputs; dx, L_tilde, the achievable
// dynamic range and the run's stop time (tau_end) are all derived.
//
// dx: R(tau)*m_r(tau) is exactly constant through the whole fat phase (any
// a_inv) -- not merely "set at the start" but literally unchanging -- so
// choosing dx from N2 evaluated at the switch is equivalent to evaluating
// it anywhere in the fat phase, including the pre-evolution/main handoff.
// R_switch, m_r_switch are passed in as plain doubles (not a Background)
// to keep this file AMReX/Background-free, matching every other function
// here; the caller (AxionStringsParams.hpp) evaluates them.
//
// tau_end: H(tau) = R'(tau)/R(tau)^2 does not depend on c at all (only
// lambda does, per Background.hpp -- R(tau) is identical whichever side of
// a switch tau falls on), and R(tau) = R0(tau/tau0)^(1/b_inv) gives the
// closed form H(tau) ~ tau^(-a_inv/b_inv) unconditionally. The achievable
// dynamic range D = log(H_switch/H_end) (moore_max_log_dynamic_range,
// natural-log form used consistently here rather than the base-appropriate
// log(H0/H)_max notation in conventions.md, which differs only by the
// already-included factor of 2) therefore inverts to a closed-form
// tau_end = tau_switch * exp(D * b_inv / a_inv), with no separate ODE
// solve needed and no dependence on how long the preceding fat phase took.
struct MooreBoxPlan
{
    double dx{};
    double L_tilde{};
    double D{};       // achievable log(H_switch/H_end)
    double tau_end{};
};

[[nodiscard]] inline MooreBoxPlan
compute_moore_box_plan(int N, double N1, double N2, double gamma,
                       double a_inv, double b_inv, double R_switch,
                       double m_r_switch, double tau_switch)
{
    MooreBoxPlan plan{};
    plan.dx      = 1.0 / (N2 * R_switch * m_r_switch);
    plan.L_tilde = N * plan.dx;
    plan.D       = moore_max_log_dynamic_range(N, N1, N2, gamma);
    plan.tau_end = tau_switch * std::exp(plan.D * b_inv / a_inv);
    return plan;
}

#endif // BOXPLAN_HPP_
