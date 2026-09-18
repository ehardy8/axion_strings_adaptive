#ifndef BACKGROUND_HPP_
#define BACKGROUND_HPP_

// Analytic FRW background and c(tau) scheme of docs/conventions.md sec.2-5.
//
// Deliberately independent of AMReX/GRTeclyn: R(tau), lambda(tau) and their
// derived quantities are analytic functions of conformal time, evaluated on
// the host and passed into the RHS kernel as plain values (see
// AxionStringsRHS). Kept free of any framework dependency so it can be unit
// tested directly (see tests/test_background.cpp).
//
// Code units (sec.5): v = 1, lambda0 = 1, tau0 = 1/v = 1, R0 = 1/((a_inv-1)
// sqrt(lambda0)) = 1/b_inv. None of these are free parameters -- a_inv is
// the only physical input here.

#include <cmath>

// c(tau): a piecewise constant with at most one switch (sec.4). Re-anchoring
// continuity of lambda across the switch is handled by Background::lambda,
// not here -- this struct only answers "what is c now".
struct CTauSchedule
{
    double c0{0.0};
    double c1{0.0};
    double tau_switch{0.0};
    bool has_switch{false};

    [[nodiscard]] double c(double tau) const
    {
        return (has_switch && tau >= tau_switch) ? c1 : c0;
    }
};

class Background
{
  public:
    static constexpr double lambda0 = 1.0; // fixed by code-unit convention
    static constexpr double tau0     = 1.0; // = 1/v, v = 1

    double a_inv{2.0};
    double b_inv{1.0}; // = a_inv - 1
    double R0{1.0};    // = 1/b_inv
    CTauSchedule c_sched;

    Background() = default;

    Background(double a_a_inv, CTauSchedule a_c_sched)
        : a_inv(a_a_inv), b_inv(a_a_inv - 1.0), R0(1.0 / (a_a_inv - 1.0)),
          c_sched(a_c_sched)
    {
    }

    // R(tau) = R0 (tau/tau0)^(1/b_inv). Independent of c: the expansion
    // history is fixed (sec.2); only lambda responds to c (sec.4).
    [[nodiscard]] double R(double tau) const
    {
        return R0 * std::pow(tau / tau0, 1.0 / b_inv);
    }

    [[nodiscard]] double R_prime(double tau) const
    {
        return (R0 / (b_inv * tau0)) *
               std::pow(tau / tau0, 1.0 / b_inv - 1.0);
    }

    // [(1 - b_inv)/(b_inv^2 tau^2)] term in the equation of motion (sec.3).
    [[nodiscard]] double curvature_term_coeff(double tau) const
    {
        return (1.0 - b_inv) / (b_inv * b_inv * tau * tau);
    }

    // lambda(tau) = lambda0 (R/R0)^(-2c), re-anchored at the switch time so
    // that lambda (not dlambda/dtau) stays continuous (sec.4):
    //   lambda(tau) = lambda_s (R/R_s)^(-2 c_new)   for tau > tau_switch
    [[nodiscard]] double lambda(double tau) const
    {
        const double c_now = c_sched.c(tau);
        if (!c_sched.has_switch || tau < c_sched.tau_switch)
        {
            return lambda0 * std::pow(R(tau) / R0, -2.0 * c_now);
        }
        const double tau_s    = c_sched.tau_switch;
        const double R_s      = R(tau_s);
        const double lambda_s = lambda0 * std::pow(R_s / R0, -2.0 * c_sched.c0);
        return lambda_s * std::pow(R(tau) / R_s, -2.0 * c_now);
    }

    // m_r/H via the closed form of sec.5, H/m_r = (tau/tau0)^((c-a_inv)/
    // (a_inv-1)). This form assumes a single constant c all the way back to
    // tau0 -- it is only used (and only exact) in the no-switch case; see
    // H_over_mr_direct for the general, always-correct definition.
    [[nodiscard]] double H_over_mr_closed_form(double tau) const
    {
        const double c_now = c_sched.c(tau);
        return std::pow(tau / tau0, (c_now - a_inv) / (a_inv - 1.0));
    }

    // H/m_r computed directly from the analytic background: physical
    // H = Rdot/R, and cosmic dt = R dtau so Rdot = R'(tau)/R(tau), giving
    // H = R'(tau)/R(tau)^2; m_r = sqrt(lambda(tau)) (v = 1). Valid with or
    // without a c-switch, since it only uses R, R' and lambda directly.
    [[nodiscard]] double H_over_mr_direct(double tau) const
    {
        const double hubble = R_prime(tau) / (R(tau) * R(tau));
        const double m_r    = std::sqrt(lambda(tau));
        return hubble / m_r;
    }

    // Inverse of H_over_mr_closed_form (2026-09-18, with the user: a more
    // physical way to specify the main run's starting time than a raw
    // conformal time tau_i): given log(m_r/H) at the start, returns the tau
    // at which that ratio is first reached. Uses c_sched.c0 -- the value in
    // effect at/before the start -- rather than c_sched.c(tau), both because
    // that is what "the c at the start" means, and because
    // H_over_mr_closed_form's own no-switch caveat applies here too (any
    // switch is expected to happen after tau_i, not before it). Singular at
    // c0 = a_inv (Moore mode from the very start, not reached via a switch)
    // -- callers must check for that separately (AxionStringsParams::
    // read_tau_i does).
    [[nodiscard]] double tau_from_log_mr_over_h(double log_mr_over_h) const
    {
        const double c0 = c_sched.c0;
        return tau0 * std::exp(log_mr_over_h * (a_inv - 1.0) / (a_inv - c0));
    }
};

#endif // BACKGROUND_HPP_
