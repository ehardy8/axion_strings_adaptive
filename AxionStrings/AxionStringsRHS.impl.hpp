#if !defined(AXIONSTRINGSRHS_HPP_)
#error "This file should only be included through AxionStringsRHS.hpp"
#endif

#ifndef AXIONSTRINGSRHS_IMPL_HPP_
#define AXIONSTRINGSRHS_IMPL_HPP_

#include "AxionStringsRHS.hpp"

template <class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
AxionStringsRHS<deriv_t>::laplacian(
    const amrex::Real *state_ptr_ijk, int comp,
    const amrex::Array1D<int, 0, AMREX_SPACEDIM> &strides,
    amrex::Long comp_stride) const
{
    amrex::Real result = 0.0;
    FOR (i)
    {
        result += m_deriv.diff2(state_ptr_ijk + comp * comp_stride, strides(i));
    }
    return result;
}

template <class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
AxionStringsRHS<deriv_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const auto *state_ptr_ijk = state.ptr(ix, iy, iz);
    const amrex::Long comp_stride = state.stride.a[2];

    amrex::Array1D<int, 0, AMREX_SPACEDIM> strides{
        AMREX_D_DECL(1, static_cast<int>(state.stride.a[0]),
                     static_cast<int>(state.stride.a[1]))};

    const amrex::Real lap_psi1 =
        laplacian(state_ptr_ijk, c_psi1, strides, comp_stride);
    const amrex::Real lap_psi2 =
        laplacian(state_ptr_ijk, c_psi2, strides, comp_stride);

    const auto cell = state.cellData(ix, iy, iz);
    const amrex::Real psi1 = cell[c_psi1];
    const amrex::Real psi2 = cell[c_psi2];

    const amrex::Real psi_sq_minus_R_sq =
        psi1 * psi1 + psi2 * psi2 - m_R_squared;
    const amrex::Real potential_coeff = 0.5 * m_lambda * psi_sq_minus_R_sq;

    auto rhs_cell = rhs_state.cellData(ix, iy, iz);

    rhs_cell[c_psi1] = cell[c_Pi1];
    rhs_cell[c_psi2] = cell[c_Pi2];

    rhs_cell[c_Pi1] =
        lap_psi1 + m_curvature_coeff * psi1 - potential_coeff * psi1;
    rhs_cell[c_Pi2] =
        lap_psi2 + m_curvature_coeff * psi2 - potential_coeff * psi2;
}

#endif // AXIONSTRINGSRHS_IMPL_HPP_
