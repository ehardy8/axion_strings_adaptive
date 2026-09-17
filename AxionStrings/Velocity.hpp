#ifndef VELOCITY_HPP_
#define VELOCITY_HPP_

// String velocity estimator (conventions.md sec.8, adapted from
// arXiv:1707.05566; milestone-1.md task 1.10 -- velocities only here, not
// curvature or loops, both deferred as more involved per the user).
//
// Conventions.md's formula, in the psi/Pi variables:
//   gamma^2 v^2 = (1/(2R^4 m_r^2 c1^2)) |Pi - (beta/tau) psi|^2
//                 (1 - e1|psi|^2/(c1^3 R^2))
//                 - (e1/(4R^6 m_r^2 c1^5)) (Re(psi* Pi) - (beta/tau)|psi|^2)^2
// "beta" is not otherwise defined in the conventions snapshot available
// here; read as beta = 1/b_inv by analogy with the (1/(b_inv tau)) psi
// term that appears identically in the main EOM and rho_tot (sec.3/
// sec.10), and confirmed by the derivation below.
//
// Re-derived and verified (2026-09-17, with the user, applying the same
// scrutiny that caught the rho_tot R-power error): converting to the
// physical field phi = psi/v and its cosmic-time derivative
// phi_dot = (v/R^2)(Pi - psi/(b_inv tau)) (v=1; tasks 1.2/1.3/1.8's
// now-verified relations) makes every factor of R cancel exactly:
//   Pi - (beta/tau) psi = R^2 phi_dot                       (direct)
//   Re(psi* Pi) - (beta/tau)|psi|^2 = R^3 Re(phi* phi_dot)   (shown below)
// so the R^4 and R^6 prefactors divide out completely, leaving
//   gamma^2 v^2 = (1/(2 m_r^2 c1^2)) |phi_dot|^2 (1 - e1|phi|^2/c1^3)
//                 - (e1/(4 m_r^2 c1^5)) [Re(phi* phi_dot)]^2
// This R-independence is itself a strong consistency check -- a physical
// velocity should not depend on the comoving rescaling used to evolve the
// field -- and is verified directly in tests/test_velocity.cpp by
// comparing this form against the literal psi/Pi formula above at several
// (R, tau, b_inv) points.
//
// Derivation of the second identity: psi = R phi, Pi = psi' = R' phi + R
// phi' (conformal-time prime). Re(psi* Pi) = R R' |phi|^2 + R^2
// Re(phi* phi'). And (beta/tau)|psi|^2 = (1/(b_inv tau)) R^2|phi|^2 =
// R R'|phi|^2 (using R'/R = 1/(b_inv tau)). The R R'|phi|^2 terms cancel,
// leaving R^2 Re(phi* phi') = R^3 Re(phi* phi_dot) (phi' = R phi_dot).
//
// Global-string profile coefficients (conventions.md sec.8): c1=0.41222,
// e1=-0.025763 (d=0, unused here). Evaluated at the corners of each
// pierced plaquette and averaged over the network (sec.8) -- the AMReX-
// side reduction lives in VelocityKernel.hpp.
//
// AMReX-free like Background.hpp -- pure math, unit-tested standalone.

#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

namespace GlobalStringProfile
{
constexpr double c1 = 0.41222;
constexpr double e1 = -0.025763;
} // namespace GlobalStringProfile

// The physical-field form actually used by the kernel: phi1,phi2 = real/
// imaginary parts of the unrescaled field; phi_dot1,phi_dot2 = their
// cosmic-time derivatives; m_r = radial mass at the evaluation time.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
gamma_sq_v_sq_from_phi(double phi1, double phi2, double phi_dot1,
                       double phi_dot2, double m_r)
{
    using namespace GlobalStringProfile;
    const double phi_sq     = phi1 * phi1 + phi2 * phi2;
    const double phi_dot_sq = phi_dot1 * phi_dot1 + phi_dot2 * phi_dot2;
    const double re_phi_conj_phi_dot = phi1 * phi_dot1 + phi2 * phi_dot2;

    const double m_r_sq = m_r * m_r;
    const double c1_sq  = c1 * c1;
    const double term1  = (phi_dot_sq / (2.0 * m_r_sq * c1_sq)) *
                         (1.0 - e1 * phi_sq / (c1 * c1_sq));
    const double term2 = (e1 / (4.0 * m_r_sq * c1_sq * c1_sq * c1)) *
                        re_phi_conj_phi_dot * re_phi_conj_phi_dot;
    return term1 - term2;
}

// The literal psi/Pi form from conventions.md sec.8, for cross-checking
// gamma_sq_v_sq_from_phi against (tests/test_velocity.cpp) -- not used by
// the production kernel, which computes phi/phi_dot directly instead
// (avoiding recomputing the R-cancellation on every cell).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
gamma_sq_v_sq_from_psi_literal(double psi1, double psi2, double Pi1,
                               double Pi2, double R, double tau,
                               double b_inv, double m_r)
{
    using namespace GlobalStringProfile;
    const double beta_over_tau = 1.0 / (b_inv * tau);

    const double d1 = Pi1 - beta_over_tau * psi1;
    const double d2 = Pi2 - beta_over_tau * psi2;
    const double diff_sq = d1 * d1 + d2 * d2;

    const double psi_sq = psi1 * psi1 + psi2 * psi2;
    const double re_psi_conj_Pi = psi1 * Pi1 + psi2 * Pi2;
    const double re_term        = re_psi_conj_Pi - beta_over_tau * psi_sq;

    const double R2 = R * R;
    const double R4 = R2 * R2;
    const double R6 = R4 * R2;
    const double m_r_sq = m_r * m_r;
    const double c1_sq  = c1 * c1;

    const double term1 = (diff_sq / (2.0 * R4 * m_r_sq * c1_sq)) *
                        (1.0 - e1 * psi_sq / (c1 * c1_sq * R2));
    const double term2 = (e1 / (4.0 * R6 * m_r_sq * c1_sq * c1_sq * c1)) *
                        re_term * re_term;
    return term1 - term2;
}

#endif // VELOCITY_HPP_
