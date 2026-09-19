#ifndef FLATBACKGROUND_HPP_
#define FLATBACKGROUND_HPP_

// Flat-space (no cosmological expansion) background, for loop simulations
// where the physics of interest is the string dynamics itself, not
// anything cosmological (2026-09-19, with the user).
//
// R(t) = 1 identically -- not a limit of the main Background class's
// power-law R(tau) = R0*(tau/tau0)^(1/b_inv) (that form cannot represent a
// *constant* R for any finite a_inv/b_inv; R=const would need b_inv -> Inf,
// which is not a value that can be passed through the existing formulas).
// Genuinely simpler than either Background or PreEvolutionBackground: with
// R'=0, the whole psi=R*phi/v rescaling is trivial (psi=phi/v, tau=t,
// Pi=phi_dot/v), and the curvature term (which arises entirely from R''/R
// in the sec.3 rescaling) vanishes identically, leaving the textbook flat
// -space complex-scalar EOM with a symmetry-breaking potential:
//   Pi_i' = laplacian(psi_i) - (lambda/2) psi_i (|psi|^2 - 1)
// lambda = m_r^2 is the only physical input (v=1 code units, as
// everywhere else in this project) -- defaults to 1, matching the
// existing lambda0=1 convention, so "loop radius relative to m_r" is
// simply the loop radius directly in code-length units when m_r=1.
//
// AMReX-free like Background.hpp/PreEvolutionBackground.hpp -- pure math,
// unit-tested standalone.

class FlatBackground
{
  public:
    double m_r{1.0}; // radial mass (sqrt(lambda)); v=1 code units

    FlatBackground() = default;
    explicit FlatBackground(double a_m_r) : m_r(a_m_r) {}

    [[nodiscard]] double R(double /*t*/) const { return 1.0; }

    [[nodiscard]] double R_prime(double /*t*/) const { return 0.0; }

    [[nodiscard]] double curvature_term_coeff(double /*t*/) const
    {
        return 0.0;
    }

    [[nodiscard]] double lambda(double /*t*/) const { return m_r * m_r; }
};

#endif // FLATBACKGROUND_HPP_
