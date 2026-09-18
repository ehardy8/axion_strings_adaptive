#ifndef VELOCITYKERNEL_HPP_
#define VELOCITYKERNEL_HPP_

// AMReX-side reduction: evaluates gamma^2 v^2 (Velocity.hpp) at the
// corners of every pierced plaquette (conventions.md sec.8: "Evaluate at
// the corners of each pierced plaquette and average over the network"),
// reusing the exact same per-cell plaquette structure as StringFinder.hpp
// (xy, yz, zx planes anchored at the low-index corner, task 1.6).
//
// `state` must already have at least 1 valid ghost cell (state.
// FillBoundary) before calling this -- same requirement as StringFinder.hpp
// (same 7-corner stencil).

#include "PlaquetteWinding.hpp"
#include "StateVariables.hpp"
#include "Velocity.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

struct VelocityResult
{
    double sum_gamma_sq_v_sq{0.0};
    long count{0};
};

[[nodiscard]] inline VelocityResult
compute_velocity_at_pierced_corners(const amrex::MultiFab &state,
                                    amrex::Real R, amrex::Real tau,
                                    amrex::Real b_inv, amrex::Real m_r)
{
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<double, long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const auto &arrs    = state.const_arrays();
    const double inv_bt = 1.0 / (static_cast<double>(b_inv) *
                                 static_cast<double>(tau));
    const double R_d    = R;
    const double m_r_d  = m_r;

    reduce_op.eval(
        state, amrex::IntVect(0), reduce_data,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> ReduceTuple
        {
            const auto &a = arrs[box_no];

            auto phase = [&](int ii, int jj, int kk) -> double
            {
                return std::atan2(a(ii, jj, kk, c_psi2),
                                  a(ii, jj, kk, c_psi1));
            };
            auto gvv = [&](int ii, int jj, int kk) -> double
            {
                const double psi1 = a(ii, jj, kk, c_psi1);
                const double psi2 = a(ii, jj, kk, c_psi2);
                const double Pi1  = a(ii, jj, kk, c_Pi1);
                const double Pi2  = a(ii, jj, kk, c_Pi2);
                const double phi1 = psi1 / R_d;
                const double phi2 = psi2 / R_d;
                const double phi_dot1 =
                    (Pi1 - inv_bt * psi1) / (R_d * R_d);
                const double phi_dot2 =
                    (Pi2 - inv_bt * psi2) / (R_d * R_d);
                return gamma_sq_v_sq_from_phi(phi1, phi2, phi_dot1, phi_dot2,
                                             m_r_d);
            };

            const double th_000 = phase(i, j, k);
            const double th_100 = phase(i + 1, j, k);
            const double th_110 = phase(i + 1, j + 1, k);
            const double th_010 = phase(i, j + 1, k);
            const double th_001 = phase(i, j, k + 1);
            const double th_101 = phase(i + 1, j, k + 1);
            const double th_011 = phase(i, j + 1, k + 1);

            const int w_xy =
                plaquette_winding(th_000, th_100, th_110, th_010);
            const int w_yz =
                plaquette_winding(th_000, th_010, th_011, th_001);
            const int w_zx =
                plaquette_winding(th_000, th_001, th_101, th_100);

            double sum  = 0.0;
            long count = 0;

            if (w_xy != 0)
            {
                sum += gvv(i, j, k) + gvv(i + 1, j, k) +
                      gvv(i + 1, j + 1, k) + gvv(i, j + 1, k);
                count += 4;
            }
            if (w_yz != 0)
            {
                sum += gvv(i, j, k) + gvv(i, j + 1, k) +
                      gvv(i, j + 1, k + 1) + gvv(i, j, k + 1);
                count += 4;
            }
            if (w_zx != 0)
            {
                sum += gvv(i, j, k) + gvv(i, j, k + 1) +
                      gvv(i + 1, j, k + 1) + gvv(i + 1, j, k);
                count += 4;
            }

            return {sum, count};
        });

    const ReduceTuple result = reduce_data.value(reduce_op);
    VelocityResult out{};
    out.sum_gamma_sq_v_sq = amrex::get<0>(result);
    out.count             = amrex::get<1>(result);

    // amrex::ReduceOps only reduces within this rank's own boxes -- found
    // the hard way (2026-09-18): with a genuinely multi-rank run, every
    // rank silently reported only its own local sum/count as if they were
    // the whole domain's, without this. AMReX's own Reduce::Sum/Min/Max
    // free functions (AMReX_Reduce.H) never add this either -- it is
    // always the caller's job for a cross-rank total.
    amrex::ParallelDescriptor::ReduceRealSum(out.sum_gamma_sq_v_sq);
    amrex::ParallelDescriptor::ReduceLongSum(out.count);
    return out;
}

#endif // VELOCITYKERNEL_HPP_
