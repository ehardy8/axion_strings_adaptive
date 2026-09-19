#ifndef STRINGTAGGER_HPP_
#define STRINGTAGGER_HPP_

// Real tagging criterion (milestone-2 Phase 1, conventions.md sec.11,
// Buschmann et al. 2108.05368 sec. "AMR Simulation Framework"), replacing
// AxionStringsLevel's previous placeholder (FixedGridsTagger, a purely
// geometric fixed fraction of the domain -- see docs/STATUS.md).
//
// Primary criterion: pierced plaquettes at the cell's low-index corner,
// via StringFinder.hpp's plaquette_windings_at -- the exact same test the
// xi diagnostic uses, not a second copy (CLAUDE.md constraint 5).
//
// Secondary (optional) criterion: large gradients in psi, matching
// Buschmann et al.'s dx_ell^2 |laplacian(psi_i)| > threshold (their
// Laplacian, not a first-derivative gradient, despite "gradient
// criterion" being the common shorthand for it). Conventions.md sec.11 is
// explicit that this (and their coarse-level radial-mode criterion, not
// implemented here) are "to be assessed rather than copied" -- so this is
// wired up now, off by default via a threshold no physical configuration
// can cross, so it can be switched on with a single runtime parameter
// later without new code (2026-09-19, with the user).
//
// Buffering (conventions.md sec.11: the fastest string, v=c, must stay a
// full core width from any coarse-fine boundary between regrids) is not
// implemented here -- it is exactly what AMReX's own amr.n_error_buf
// already does (grows every tagged region by that many cells before
// clustering into grids), so it is a run-time parameter to set
// (Buschmann et al.: 11 cells, with regrid interval 0.2/2^level -- a
// starting point to test, not assumed correct for our string density),
// not code to write.
//
// `state` must already have at least 2 valid ghost cells (state.
// FillBoundary) before tagging: 1 for the plaquette test's (i+1,j+1,k+1)
// neighbours, 2 for FourthOrderDerivatives' Laplacian stencil.

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "StringFinder.hpp"

#include <AMReX_Array4.H>
#include <AMReX_TagBox.H>

#include <cmath>
#include <limits>

struct StringTaggerParams
{
    // Buschmann et al.'s value is 0.04; ours defaults to "off" (see header
    // comment) since conventions.md sec.11 says to assess this criterion
    // rather than copy it -- set explicitly, e.g. axion_strings.tagging.
    // gradient_threshold = 0.04, to turn it on.
    double gradient_threshold{std::numeric_limits<double>::max()};
};

class StringTagger
{
  public:
    AMREX_GPU_HOST_DEVICE
    StringTagger(amrex::Real a_dx, StringTaggerParams a_params)
        : m_deriv(a_dx), m_dx(a_dx), m_params(a_params)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
              const amrex::Array4<amrex::Real const> &state,
              const amrex::Array4<amrex::TagBox::TagType> &tags) const
    {
        const auto w = plaquette_windings_at(state, ix, iy, iz);
        if (w.w_xy != 0 || w.w_yz != 0 || w.w_zx != 0)
        {
            tags(ix, iy, iz) = amrex::TagBox::SET;
            return;
        }

        // Cheapest-first: only reach the Laplacian (2nd-derivative
        // stencil) for cells the plaquette test didn't already tag.
        const double lap1 = scalar_laplacian(state, ix, iy, iz, c_psi1);
        const double lap2 = scalar_laplacian(state, ix, iy, iz, c_psi2);
        const double dx2  = m_dx * m_dx;
        if (dx2 * std::abs(lap1) > m_params.gradient_threshold ||
            dx2 * std::abs(lap2) > m_params.gradient_threshold)
        {
            tags(ix, iy, iz) = amrex::TagBox::SET;
        }
    }

  private:
    FourthOrderDerivatives m_deriv;
    amrex::Real m_dx;
    StringTaggerParams m_params;

    // Matches AxionStringsRHS's own (private) laplacian() exactly -- the
    // same already-proven-correct low-level Array4 pointer/stride access
    // pattern, not a new one.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE double
    scalar_laplacian(const amrex::Array4<amrex::Real const> &state, int ix,
                     int iy, int iz, int comp) const
    {
        const auto *state_ptr_ijk     = state.ptr(ix, iy, iz);
        const amrex::Long comp_stride = state.stride.a[2];
        const amrex::Array1D<int, 0, AMREX_SPACEDIM> strides{AMREX_D_DECL(
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1]))};

        double result = 0.0;
        FOR (i)
        {
            result +=
                m_deriv.diff2(state_ptr_ijk + comp * comp_stride, strides(i));
        }
        return result;
    }
};

#endif // STRINGTAGGER_HPP_
