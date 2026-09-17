#ifndef PREEVOLUTIONBACKGROUND_HPP_
#define PREEVOLUTIONBACKGROUND_HPP_

// The pre-evolution background of conventions.md sec.7: R = R0(t/t0),
// lambda(t) = lambda(t0)(t/t0)^-2 in cosmic time t, chosen so that R*m_r and
// m_r/H are both time-independent -- this is what lets the network relax
// near the attractor before the main run's clock starts.
//
// This is NOT an instance of Background (task 1.3): R linear in cosmic time
// corresponds to a_inv=1, which is singular there (b_inv=a_inv-1=0, so
// R0=1/b_inv diverges). It is also not fat mode (c=1) under radiation
// domination -- checked: that gives m_r/H ~ tau, not constant. Derived
// independently (confirmed with the user, see docs/STATUS.md for the
// derivation): switching to pre-evolution's own conformal time tau_pre
// (dtau_pre = dt/R(t)) turns R(t)=R0(t/t0) into R(tau_pre) = R0 e^(alpha
// tau_pre), alpha = R0/t0. The same psi=R phi/v rescaling trick from sec.3
// (generic in R(tau), not tied to the power-law form) then gives an EOM
// with the same Laplacian/curvature/potential structure as the main RHS,
// with a *constant* curvature coefficient alpha^2 and
// lambda(tau_pre) = lambda_pre0 (R(tau_pre)/R0)^-2 (the c=1-shaped formula,
// with this new R). R*m_r and m_r/H are both exactly constant under this,
// matching sec.7's stated requirement (verified analytically and in
// tests/test_pre_evolution_background.cpp).
//
// R0 and alpha are pre-evolution's own time-gauge freedom (only ratios of
// R, and the elapsed tau_pre, are physical) -- fixed here to R0=alpha=1
// WLOG, leaving one physical input: gamma_pre, the constant m_r/H value to
// relax the network at. This means AxionStringsRHS is reused completely
// unchanged for pre-evolution -- only R(tau), lambda(tau) and the curvature
// coefficient differ from the main Background.
//
// AMReX-free like Background.hpp/BoxPlan.hpp -- pure math, unit-tested
// standalone.

#include <cmath>

class PreEvolutionBackground
{
  public:
    double gamma_pre{1.0}; // target (constant) m_r/H during pre-evolution

    PreEvolutionBackground() = default;
    explicit PreEvolutionBackground(double a_gamma_pre) : gamma_pre(a_gamma_pre)
    {
    }

    [[nodiscard]] double R(double tau_pre) const
    {
        return std::exp(tau_pre); // R0 = alpha = 1
    }

    [[nodiscard]] double R_prime(double tau_pre) const
    {
        return std::exp(tau_pre); // = R(tau_pre), since alpha = 1
    }

    // constant: alpha^2 = 1
    [[nodiscard]] double curvature_term_coeff(double /*tau_pre*/) const
    {
        return 1.0;
    }

    [[nodiscard]] double lambda(double tau_pre) const
    {
        const double r = R(tau_pre);
        return gamma_pre * gamma_pre / (r * r); // lambda_pre0 * (R/R0)^-2, R0=1
    }

    // Round-off cross-check (same pattern as Background::H_over_mr_direct,
    // computed from R, R' and lambda directly): should equal gamma_pre
    // exactly, independent of tau_pre.
    [[nodiscard]] double m_r_over_H_direct(double tau_pre) const
    {
        const double hubble = R_prime(tau_pre) / (R(tau_pre) * R(tau_pre));
        const double m_r    = std::sqrt(lambda(tau_pre));
        return m_r / hubble;
    }
};

#endif // PREEVOLUTIONBACKGROUND_HPP_
