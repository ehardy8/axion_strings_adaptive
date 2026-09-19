#ifndef ENERGY_HPP_
#define ENERGY_HPP_

// Total energy density: re-derived from scratch (2026-09-17, with the
// user) rather than taken literally from conventions.md sec.10, after
// finding what appears to be a factor-of-R^2 error in the kinetic/gradient
// terms there. Starting point, confirmed with the user: a canonically
// normalised COMPLEX scalar has NO 1/2 on its kinetic/gradient terms
// (unlike a real scalar) -- checked directly by matching
// L = A|phi_dot|^2 - A|grad phi|^2/R^2 - V(phi) against sec.3's given EOM
// (phi_ddot - grad^2(phi)/R^2 + (lambda/2)phi(|phi|^2-v^2) = 0) via
// d/dt(dL/d(phi_dot)*) = dV/dphi*, which forces A=1 exactly (a factor A=1/2
// would instead require the EOM's potential-derivative coefficient to be
// lambda, not lambda/2). So the physical energy density is
//   rho = |phi_dot|^2 + |grad phi|^2/R^2 + V(phi),   phi_dot = d(phi)/dt.
//
// Converting to psi = R phi/v and conformal time tau (Pi = d(psi)/d(tau)):
// phi = v psi/R, and using dt = R dtau and R'/R = 1/(b_inv tau) (sec.5):
//   d(phi)/d(tau) = (v/R)(Pi - psi/(b_inv tau))              [conformal]
//   phi_dot = [d(phi)/d(tau)]/R = (v/R^2)(Pi - psi/(b_inv tau))
// (checked three independent ways: direct algebra, this formula, and
// re-parametrising by t(tau) directly -- all agree). So
//   |phi_dot|^2 = (v^2/R^4)|Pi - psi/(b_inv tau)|^2,
// i.e. R^-4, not conventions.md's literal R^-2. The gradient term picks up
// the same extra 1/R^2 (grad_physical = grad_comoving/R, sec.3: "grad now
// the comoving gradient"), so it is R^-4 too. The potential term, once
// fully expanded, is identical either way:
//   V(phi) = (lambda/4)(|phi|^2-v^2)^2 = (lambda v^4)/(4R^4)(|psi|^2-R^2)^2
// -- conventions.md's literal formula and this derivation agree on the
// potential term exactly; only the kinetic/gradient terms differ, by R^2.
// v=1 in code units (sec.5).
//
// Pointwise (pre-averaging) formula only, AMReX-free like Background.hpp;
// the AMReX-side gradient computation and spatial reduction live in
// EnergyKernel.hpp.
//
// This is the *aggregate* formula only -- conventions.md sec.12 lists names
// for a full radial/axion/interaction decomposition (radial kinetic/
// gradient/mass, axion kinetic/gradient, interaction total/radial-side/
// axion-side, total field kinetic/potential/spatial) but does not give
// their formulas in terms of psi/Pi explicitly; that decomposition is not
// yet implemented here (see docs/STATUS.md).

#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
rho_tot_pointwise(double psi1, double psi2, double Pi1, double Pi2,
                  double grad_psi1_sq, double grad_psi2_sq, double R,
                  double lambda, double b_inv, double tau)
{
    const double c1 = Pi1 - psi1 / (b_inv * tau);
    const double c2 = Pi2 - psi2 / (b_inv * tau);
    const double kinetic  = c1 * c1 + c2 * c2;
    const double gradient = grad_psi1_sq + grad_psi2_sq;

    const double psi_sq_minus_R_sq = psi1 * psi1 + psi2 * psi2 - R * R;
    const double potential =
        0.25 * lambda * psi_sq_minus_R_sq * psi_sq_minus_R_sq;

    const double R2 = R * R;
    return (kinetic + gradient + potential) / (R2 * R2);
}

// Axion kinetic energy density, conventions.md sec.12's "Axion: kinetic
// 1/2 a_dot^2": a = f_a theta (f_a = sqrt(2) v = sqrt(2) in code units,
// sec.3), a_dot = f_a theta'/R (cosmic time from conformal: dt = R dtau),
// theta' = (psi1 Pi2 - Pi1 psi2)/|psi|^2 -- the same phase-rate identity
// Masking.hpp's "None" mode computes. Takes that already-computed theta'
// (masked_a_dot with MaskingScheme::None) rather than recomputing it, so
// there is exactly one place theta' is derived from psi/Pi (CLAUDE.md
// constraint 5's spirit, applied beyond masking itself).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
axion_kinetic_energy_pointwise(double theta_prime, double R)
{
    constexpr double f_a = 1.4142135623730951; // sqrt(2) v, v = 1
    const double a_dot   = f_a * theta_prime / R;
    return 0.5 * a_dot * a_dot;
}

// Axion gradient energy density, conventions.md sec.12's "Axion: gradient
// 1/2 (grad a)^2": a = f_a theta, and grad(a) (physical) = f_a
// grad(theta)/R (theta is a pure phase, so its physical gradient picks up
// only one power of 1/R relative to the comoving grad(theta) -- unlike
// psi's amplitude, which picks up two, sec.3: "grad now the comoving
// gradient"). grad(theta) in comoving direction d is
// (psi1 d_d(psi2) - psi2 d_d(psi1))/|psi|^2, the spatial analogue of the
// theta' identity used above; takes the already-summed comoving
// |grad(theta)|^2 rather than recomputing it component by component.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
axion_gradient_energy_pointwise(double grad_theta_sq_comoving, double R)
{
    constexpr double f_a = 1.4142135623730951; // sqrt(2) v, v = 1
    return 0.5 * f_a * f_a * grad_theta_sq_comoving / (R * R);
}

// Radial (Higgs/amplitude) energy density, conventions.md sec.12's
// "Radial: kinetic/gradient/mass" (2026-09-18, with the user: not yet
// implemented before now -- sec.12 names these but gives no psi/Pi
// formula). r = |phi| - v is the amplitude's own deviation from the vev;
// these three functions are its exact kinetic/gradient/potential energy,
// defined the same way axion_kinetic/gradient_energy_pointwise define the
// *phase*-only energy above: the physical energy carried by that one
// degree of freedom alone, not "total minus axion" (which would also
// include the interaction/cross terms conventions.md sec.12 lists as
// separate, un-requested observables).
//
// |phi| = v|psi|/R (v=1, code units), so d|phi|/dt (cosmic time) is the
// *radial* projection of phi_dot's numerator (c1,c2) = (Pi1,Pi2) -
// psi/(b_inv tau) onto psi's own direction psi_hat = psi/|psi|:
//   d|phi|/dt = (v/R^2) (psi1 c1 + psi2 c2)/|psi|
// -- exact, not linearised: the *orthogonal* (tangential) projection of
// the same (c1,c2) is exactly f_a*theta_prime (Masking.hpp's identity),
// confirmed algebraically (the (1/(b_inv tau))*psi_i terms cancel exactly
// in psi1 c2 - psi2 c1, leaving psi1 Pi2 - Pi1 psi2), so kinetic + radial
// = tangential exactly decomposes the full |phi_dot|^2 via Pythagoras --
// checked by expanding both sides. psi=0 is a measure-zero core point
// where "radial direction" is undefined; guarded to 0, the same fallback
// Masking.hpp/EnergyKernel.hpp already use for theta_prime there.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
radial_kinetic_energy_pointwise(double psi1, double psi2, double Pi1,
                                double Pi2, double R, double b_inv,
                                double tau)
{
    const double c1 = Pi1 - psi1 / (b_inv * tau);
    const double c2 = Pi2 - psi2 / (b_inv * tau);
    const double psi_sq = psi1 * psi1 + psi2 * psi2;
    if (psi_sq <= 0.0)
    {
        return 0.0;
    }
    const double radial_dot = (psi1 * c1 + psi2 * c2) / std::sqrt(psi_sq);
    const double R2 = R * R;
    return (radial_dot * radial_dot) / (R2 * R2);
}

// Radial gradient energy: |grad|phi||^2, physical. grad|psi| (comoving)
// follows the same per-direction orthogonal split as the kinetic term
// above -- an exact Pythagorean identity (expand both sides; the cross
// terms cancel), not an approximation:
//   (d_d psi1)^2 + (d_d psi2)^2 = (d_d|psi|)^2 + |psi|^2 (d_d theta)^2
// so |grad|psi||^2 = (grad_psi1_sq + grad_psi2_sq) - psi_sq*grad_theta_sq
// -- reuses the same intermediates EnergyKernel.hpp already computes for
// rho_tot/axion_gradient (grad_psi1_sq, grad_psi2_sq, the comoving
// grad_theta_sq), no separate per-direction loop needed. Well-defined at
// psi=0 with no extra guard: grad_theta_sq is already 0 there (its own
// psi_sq>0 guard, upstream), so this correctly reduces to the full
// (grad_psi1_sq+grad_psi2_sq) -- consistent, since at the exact core all
// gradient energy is "radial" by construction and axion_gradient is 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
radial_gradient_energy_pointwise(double grad_psi1_sq, double grad_psi2_sq,
                                 double psi_sq, double grad_theta_sq_comoving,
                                 double R)
{
    const double grad_mod_sq_comoving =
        grad_psi1_sq + grad_psi2_sq - psi_sq * grad_theta_sq_comoving;
    const double R2 = R * R;
    return grad_mod_sq_comoving / (R2 * R2);
}

// Radial (Higgs) mass/potential energy: the same (lambda/4)(|phi|^2-v^2)^2
// term already inside rho_tot_pointwise's "potential" local, exposed here
// standalone so it can be saved as its own diagnostic.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
radial_mass_energy_pointwise(double psi1, double psi2, double R,
                             double lambda)
{
    const double psi_sq_minus_R_sq = psi1 * psi1 + psi2 * psi2 - R * R;
    const double potential =
        0.25 * lambda * psi_sq_minus_R_sq * psi_sq_minus_R_sq;
    const double R2 = R * R;
    return potential / (R2 * R2);
}

#endif // ENERGY_HPP_
