#ifndef ENERGY_HPP_
#define ENERGY_HPP_

// Total energy density of conventions.md sec.10:
//   rho_tot = R^-2 v^2 <|psi_dot - psi/(b_inv tau)|^2 + |grad psi|^2
//             + (lambda v^2/4R^2)(|psi|^2-R^2)^2>
// with v=1 in code units (sec.5), "psi_dot" read as Pi = dpsi/dtau (the
// evolved state variable -- sec.3's derivation of this formula starts from
// phi' = (v/R)(Pi - psi/(b_inv tau)), so this is what the symbol must mean
// here). Pointwise (pre-averaging) formula only, AMReX-free like
// Background.hpp; the AMReX-side gradient computation and spatial
// reduction live in EnergyKernel.hpp.
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
    const double potential = (lambda / (4.0 * R * R)) * psi_sq_minus_R_sq *
                             psi_sq_minus_R_sq;

    return (kinetic + gradient + potential) / (R * R);
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

#endif // ENERGY_HPP_
