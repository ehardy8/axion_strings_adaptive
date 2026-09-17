#ifndef MASKEDFIELDBUFFER_HPP_
#define MASKEDFIELDBUFFER_HPP_

// The one place (CLAUDE.md constraint 5) the axion time-derivative is
// masked before being written into the FFT buffer for a spectrum
// (milestone-1.md task 1.9): fills `buffer` with masked_a_dot
// (Masking.hpp) evaluated from the current state. Nothing else in the
// codebase should apply masking to a field destined for an FFT.

#include "Masking.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>

inline void fill_masked_a_dot_buffer(amrex::MultiFab &buffer,
                                     const amrex::MultiFab &state,
                                     const MaskingParams &params,
                                     amrex::Real R)
{
    const auto &state_arrs  = state.const_arrays();
    const auto &buffer_arrs = buffer.arrays();
    amrex::ParallelFor(
        buffer,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            const auto &a           = state_arrs[box_no];
            const amrex::Real psi1  = a(i, j, k, c_psi1);
            const amrex::Real psi2  = a(i, j, k, c_psi2);
            const amrex::Real Pi1   = a(i, j, k, c_Pi1);
            const amrex::Real Pi2   = a(i, j, k, c_Pi2);
            buffer_arrs[box_no](i, j, k) = static_cast<amrex::Real>(
                masked_a_dot(params, psi1, psi2, Pi1, Pi2, R));
        });
    amrex::Gpu::streamSynchronize();
}

#endif // MASKEDFIELDBUFFER_HPP_
