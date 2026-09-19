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
#include <vector>

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

// AMR box plan (2026-09-19, with the user, milestone-2 "Phase 0"):
// conventions.md sec.5's own note -- "for AMR, N and N2 refer to the
// effective finest resolution, not the base grid" -- means
// compute_general_box_plan already gives the right L_tilde/tau_f/dx_finest
// when called with N = N_effective (the box's finest-level cell count).
// What's new here: deriving the *base* (level 0) grid size from
// N_effective/2^max_level, and the level-addition schedule -- the
// log(m_r/H) at which each level must come online so the finest active
// level's resolution never falls below the N2 target.
//
// Derivation (physical mode, c -- see the singularity note below):
// N2 at a fixed comoving spacing dx is N2(tau) = 1/(R(tau) dx m_r(tau));
// R(tau) m_r(tau) = R0 (tau/tau0)^((1-c)/b_inv) (Background.hpp's R(tau)
// and m_r(tau) = sqrt(lambda(tau)), lambda0=1), so at level 0's spacing
// dx_0, N2_0(tau) is a pure power law in tau, hence in x = m_r/H via
// x = tau^((a_inv-c)/b_inv). Level ell's spacing is dx_0/2^ell, and
// N2 scales as 1/dx, so N2_ell(tau) = 2^ell N2_0(tau) *identically*
// (a purely geometric statement, independent of c or m_r(tau)'s own time
// dependence). Requiring level ell to just reach the target N2 at the
// moment it is needed (2^(ell-1) N2_0 = N2, i.e. level ell-1's resolution
// has *just* degraded to the target) and solving for the corresponding x
// gives a power law in ell, whose consecutive spacing in log(x) works out
// to (checked numerically against the closed form and against
// conventions.md sec.11's ln(4) for a_inv=2, c=0 -- see
// tests/test_box_plan.cpp):
//   Delta log(m_r/H) per level = ln(2) * (a_inv - c) / (1 - c)
// Singular at c=1 (fat) -- consistent with conventions.md sec.11's own
// table entry that fat mode's comoving core width is exactly constant, so
// no further level is ever needed past the initial setup; this function
// is intended for physical mode (c=0), where the singularity is nowhere
// near (denominator = 1). Counting backward from tau_f (where, by
// compute_general_box_plan's own construction, N2 is reached using *all*
// max_level levels) by one level-spacing per level gives log_add[ell]
// (1-indexed: log_add[0] is level 1's own threshold).
struct AmrBoxPlan
{
    double L_tilde{};
    double dx_finest{};
    double dx_base{};
    int N_base{};   // -1 if N_effective is not evenly divisible by 2^max_level
    double tau_f{}; // from compute_general_box_plan(N_effective, ...)
    std::vector<double> log_add; // size max_level; log_add[ell-1] = level ell's threshold
};

[[nodiscard]] inline AmrBoxPlan
compute_amr_box_plan(int N_effective, double N1, double N2, double a_inv,
                     double c, int max_level)
{
    const GeneralBoxPlan general =
        compute_general_box_plan(N_effective, N1, N2, a_inv, c);

    AmrBoxPlan plan{};
    plan.L_tilde   = general.L_tilde;
    plan.dx_finest = general.delta_x;
    plan.tau_f     = general.tau_f;

    const int ratio = 1 << max_level;
    plan.N_base = (N_effective % ratio == 0) ? (N_effective / ratio) : -1;
    plan.dx_base = plan.dx_finest * static_cast<double>(ratio);

    const double b_inv = a_inv - 1.0;
    const double log_mr_over_h_at_tau_f =
        (a_inv - c) / b_inv * std::log(general.tau_f);
    const double delta_log_per_level = std::log(2.0) * (a_inv - c) / (1.0 - c);

    plan.log_add.resize(static_cast<std::size_t>(max_level));
    for (int ell = 1; ell <= max_level; ++ell)
    {
        plan.log_add[static_cast<std::size_t>(ell - 1)] =
            log_mr_over_h_at_tau_f -
            static_cast<double>(max_level - ell + 1) * delta_log_per_level;
    }
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
