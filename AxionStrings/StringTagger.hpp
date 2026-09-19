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
// Third criterion, originally added for flat-space loop simulations
// (2026-09-19, with the user), and tried (same day) as an also-on-by-
// default option for the physical-mode network schedule too ("using the
// gradient of the radial mode instead of (or as well as) the winding,
// but... keep increasing the max refinement level based on when m_r
// Delta rises above the user specified value" -- i.e. this only ever
// changes which per-cell criterion tags cells within an already-
// permitted level, never AxionStringsLevel::tag_cells()'s own schedule-
// gating of *when* a level is allowed to exist): |grad(r)|^2, r=|psi| the
// radial (Higgs) amplitude mode -- the exact same comoving quantity
// Energy.hpp's radial_gradient_energy_pointwise computes (grad_psi1_sq +
// grad_psi2_sq - psi_sq*grad_theta_sq, the orthogonal radial/tangential
// decomposition of |grad(psi)|^2), reused here rather than re-derived.
// Motivation: a string element moving with local Lorentz factor Gamma
// has its core width contracted by 1/Gamma *in the lab/comoving frame*
// along its direction of motion, so the resolution this project's dx
// must supply scales with Gamma there -- but |grad(r)|, a lab-frame
// spatial derivative, already picks this up directly and automatically
// (a contracted core is *steeper* in lab-frame coordinates, regardless
// of why it's contracted), with no separate velocity estimate needed.
// For a stationary string this reduces to an ordinary (unboosted)
// resolution test on the core profile; it needs no boost-specific term
// at all.
//
// Off by default for network mode too, same as gradient_threshold above
// (2026-09-19, with the user, after testing): a physically-derived
// default (calibrated from an isolated static string's own peak core
// gradient, tied to BoxPlan.hpp's N2 target -- see AxionStringsParams::
// read_tagging_params()'s own comment for the full derivation) was tried
// against params_amr_validation_128.txt and tagged 100% of the domain at
// level 1's very first onset, not the ~37% the plaquette criterion alone
// gives there -- an isolated-core calibration does not describe the
// ambient radial-mode gradient in a real, dense Fourier-relaxed tangle.
// Kept alongside the plaquette criterion above (both can tag; neither
// replaces the other), available as an explicit opt-in via axion_strings.
// tagging.radial_gradient_threshold for whoever revisits the calibration.
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

    // Boost-aware radial-gradient criterion (see header comment) -- off by
    // default (same "genuinely off" convention as gradient_threshold
    // above), axion_strings.tagging.radial_gradient_threshold to enable.
    double radial_gradient_threshold{std::numeric_limits<double>::max()};
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

        const double dx2 = m_dx * m_dx;

        // Cheapest-first: only reach the Laplacian (2nd-derivative
        // stencil) for cells the plaquette test didn't already tag.
        const double lap1 = scalar_laplacian(state, ix, iy, iz, c_psi1);
        const double lap2 = scalar_laplacian(state, ix, iy, iz, c_psi2);
        if (dx2 * std::abs(lap1) > m_params.gradient_threshold ||
            dx2 * std::abs(lap2) > m_params.gradient_threshold)
        {
            tags(ix, iy, iz) = amrex::TagBox::SET;
            return;
        }

        // sqrt, not a squared threshold: squaring the default
        // std::numeric_limits<double>::max() "off" sentinel would itself
        // overflow (to +inf, which would still compare correctly, but
        // there is no reason to rely on that).
        if (m_dx * std::sqrt(radial_gradient_sq(state, ix, iy, iz)) >
           m_params.radial_gradient_threshold)
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

    // |grad(r)|^2, r=|psi|, in comoving units -- Energy.hpp's
    // radial_gradient_energy_pointwise's exact orthogonal decomposition
    // (grad_psi1_sq + grad_psi2_sq - psi_sq*grad_theta_sq), before its
    // /R^4 physical normalisation (not needed here: the tagger, like the
    // Laplacian criterion above, works directly in comoving/grid units).
    // Reuses the same low-level pointer/stride access as scalar_laplacian,
    // just diff1 (first derivative) instead of diff2.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE double
    radial_gradient_sq(const amrex::Array4<amrex::Real const> &state, int ix,
                       int iy, int iz) const
    {
        const double psi1   = state(ix, iy, iz, c_psi1);
        const double psi2   = state(ix, iy, iz, c_psi2);
        const double psi_sq = psi1 * psi1 + psi2 * psi2;

        const auto *state_ptr_ijk     = state.ptr(ix, iy, iz);
        const amrex::Long comp_stride = state.stride.a[2];
        const amrex::Array1D<int, 0, AMREX_SPACEDIM> strides{AMREX_D_DECL(
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1]))};

        double grad_psi1_sq  = 0.0;
        double grad_psi2_sq  = 0.0;
        double grad_theta_sq = 0.0;
        FOR (i)
        {
            const double d_psi1 = m_deriv.diff1(
                state_ptr_ijk + c_psi1 * comp_stride, strides(i));
            const double d_psi2 = m_deriv.diff1(
                state_ptr_ijk + c_psi2 * comp_stride, strides(i));
            grad_psi1_sq += d_psi1 * d_psi1;
            grad_psi2_sq += d_psi2 * d_psi2;
            const double dtheta_dir =
                (psi_sq > 0.0) ? (psi1 * d_psi2 - psi2 * d_psi1) / psi_sq
                              : 0.0;
            grad_theta_sq += dtheta_dir * dtheta_dir;
        }
        return grad_psi1_sq + grad_psi2_sq - psi_sq * grad_theta_sq;
    }
};

#endif // STRINGTAGGER_HPP_
