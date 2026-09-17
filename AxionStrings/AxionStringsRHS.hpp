#ifndef AXIONSTRINGSRHS_HPP_
#define AXIONSTRINGSRHS_HPP_

// GRTeclyn includes
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"

// Problem specific includes
#include "StateVariables.hpp"

// Complex-scalar RHS of conventions.md sec.3:
//
//   Pi_i' = laplacian(psi_i) + curvature_coeff * psi_i
//           - (lambda/2) * psi_i * (|psi|^2 - R^2)
//   psi_i' = Pi_i
//
// curvature_coeff, lambda and R_squared are analytic functions of conformal
// time (Background::curvature_term_coeff, Background::lambda, R(tau)^2);
// they are evaluated on the host from a_time and passed in here as plain
// values (milestone-1.md task 1.3, CLAUDE.md constraint 2: no ParmParse and
// no host-only calls inside the kernel).
//
// KO dissipation is omitted here by design: CLAUDE.md constraint 1 fixes
// sigma = 0 as the project default because it damps exactly the high-k
// modes q is measured from. If a nonzero sigma later proves necessary for
// stability it should be reintroduced deliberately, with its effect on the
// spectrum characterised (see conventions.md sec.6 and sec.14 decision log).
template <class deriv_t = FourthOrderDerivatives>
class AxionStringsRHS
{
  public:
    AxionStringsRHS(amrex::Real a_dx, amrex::Real a_curvature_coeff,
                    amrex::Real a_lambda, amrex::Real a_R_squared)
        : m_deriv(a_dx), m_curvature_coeff(a_curvature_coeff),
          m_lambda(a_lambda), m_R_squared(a_R_squared)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &rhs_state,
               const amrex::Array4<amrex::Real const> &state) const;

  private:
    deriv_t m_deriv;
    amrex::Real m_curvature_coeff;
    amrex::Real m_lambda;
    amrex::Real m_R_squared;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real laplacian(
        const amrex::Real *state_ptr_ijk, int comp,
        const amrex::Array1D<int, 0, AMREX_SPACEDIM> &strides,
        amrex::Long comp_stride) const;
};

#include "AxionStringsRHS.impl.hpp"

#endif
