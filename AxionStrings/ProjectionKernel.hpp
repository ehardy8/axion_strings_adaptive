#ifndef PROJECTIONKERNEL_HPP_
#define PROJECTIONKERNEL_HPP_

// Line-of-sight projection of the total energy density onto the xy-plane,
// for optional visualisation snapshots (output follow-up, 2026-09-18, with
// the user: "make pictures of the simulation at the particular
// timeshots... this should be an option, not always done").
//
// Two deliberate choices, both per the user:
//  - The *unscreened* rho_tot is projected, not the screened one --
//    screening removes exactly the string cores this is meant to show.
//  - Each column takes the MAX of rho_tot along the line of sight (fixed
//    to the z-axis here), not a sum/average: a thin, high-density core
//    would otherwise be diluted by a long, mostly quiescent sight line,
//    whereas the max is exactly what makes a string visible against the
//    background.
//
// The string network's own xy-plaquette winding (PlaquetteWinding.hpp,
// StringFinder.hpp's convention: a pierced xy-plaquette marks a string
// segment threading roughly along z) is counted per column too, summed
// over the whole line of sight, as an independent cross-check overlay: a
// nonzero count should coincide with a peak in rho_tot_max at the same
// (i,j) if the energy/string-finding pipelines agree with each other.
//
// CPU-first (same decision as SpectrumKernel.hpp's shell-binning): a
// host-side MFIter loop, not a GPU kernel, since this runs once per output
// snapshot, not every substep.
//
// `state` must already have at least 2 valid ghost cells (state.
// FillBoundary) before calling this: FourthOrderDerivatives::diff1 needs 2
// for the energy gradient, and the plaquette winding needs 1 for its
// (i+1,j+1) corners.

#include "Energy.hpp"
#include "FourthOrderDerivatives.hpp"
#include "PlaquetteWinding.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

struct Projection
{
    int nx{0};
    int ny{0};
    // Row-major, index (i - domain.smallEnd(0)) * ny + (j - domain.
    // smallEnd(1)); both size nx*ny.
    std::vector<double> rho_tot_max_unscreened;
    std::vector<amrex::Long> string_hit_count;
};

[[nodiscard]] inline Projection
compute_energy_projection(const amrex::MultiFab &state, amrex::Real dx,
                          amrex::Real R, amrex::Real lambda,
                          amrex::Real R_prime_over_R,
                          const amrex::Box &domain)
{
    const auto lo = domain.smallEnd();
    const int nx  = domain.length(0);
    const int ny  = domain.length(1);

    Projection out{};
    out.nx = nx;
    out.ny = ny;
    out.rho_tot_max_unscreened.assign(
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny),
        -std::numeric_limits<double>::infinity());
    out.string_hit_count.assign(
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0);

    const FourthOrderDerivatives deriv(dx);

    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto &a        = state.const_array(mfi);
        const auto blo       = bx.smallEnd();
        const auto bhi       = bx.bigEnd();

        auto phase = [&](int ii, int jj, int kk) -> double
        {
            return std::atan2(a(ii, jj, kk, c_psi2), a(ii, jj, kk, c_psi1));
        };

        for (int k = blo[2]; k <= bhi[2]; ++k)
        {
            for (int j = blo[1]; j <= bhi[1]; ++j)
            {
                for (int i = blo[0]; i <= bhi[0]; ++i)
                {
                    const double psi1 = a(i, j, k, c_psi1);
                    const double psi2 = a(i, j, k, c_psi2);
                    const double Pi1  = a(i, j, k, c_Pi1);
                    const double Pi2  = a(i, j, k, c_Pi2);

                    const auto grad_psi1 = deriv.d1_scalar(i, j, k, a, c_psi1);
                    const auto grad_psi2 = deriv.d1_scalar(i, j, k, a, c_psi2);
                    double grad_psi1_sq  = 0.0;
                    double grad_psi2_sq  = 0.0;
                    FOR (dir)
                    {
                        grad_psi1_sq += grad_psi1(dir) * grad_psi1(dir);
                        grad_psi2_sq += grad_psi2(dir) * grad_psi2(dir);
                    }

                    const double rho = rho_tot_pointwise(
                        psi1, psi2, Pi1, Pi2, grad_psi1_sq, grad_psi2_sq, R,
                        lambda, R_prime_over_R);

                    const std::size_t idx =
                        static_cast<std::size_t>(i - lo[0]) *
                            static_cast<std::size_t>(ny) +
                        static_cast<std::size_t>(j - lo[1]);
                    out.rho_tot_max_unscreened[idx] =
                        std::max(out.rho_tot_max_unscreened[idx], rho);

                    const double th_000 = phase(i, j, k);
                    const double th_100 = phase(i + 1, j, k);
                    const double th_110 = phase(i + 1, j + 1, k);
                    const double th_010 = phase(i, j + 1, k);
                    const int w_xy =
                        plaquette_winding(th_000, th_100, th_110, th_010);
                    if (w_xy != 0)
                    {
                        ++out.string_hit_count[idx];
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(
        out.rho_tot_max_unscreened.data(),
        static_cast<int>(out.rho_tot_max_unscreened.size()));
    amrex::ParallelDescriptor::ReduceLongSum(
        out.string_hit_count.data(),
        static_cast<int>(out.string_hit_count.size()));
    return out;
}

#endif // PROJECTIONKERNEL_HPP_
