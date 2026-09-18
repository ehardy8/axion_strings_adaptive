#ifndef ENERGYKERNEL_HPP_
#define ENERGYKERNEL_HPP_

// AMReX-side energy reduction: computes the screened and unscreened
// spatial averages of rho_tot, axion kinetic energy and axion gradient
// energy (Energy.hpp), plus the point counts needed to reconstruct any
// derived combination afterwards -- deliberately not baking a single
// "string energy" formula in here (2026-09-17, with the user): the
// screened/unscreened *difference* only captures the energy diluted by
// the masked point-count fraction, not the core energy itself, and
// separating string energy from axion wave energy in a full network needs
// combining these raw averages in more than one way (e.g. core energy via
// n_total*unscreened - n_unmasked*screened, or -- for a static string --
// core energy plus the away-from-core (axion_gradient - axion_kinetic)
// combination, which cancels a propagating wave's equal kinetic/gradient
// contributions and leaves the string's own static long-range tail). All
// of that reconstruction happens downstream from the raw averages and
// counts saved here, every snapshot.
//
// "Screened" divides by the effective point count (conventions.md sec.10:
// "the effective number of lattice points averaged over must be multiplied
// by the same masking factor"), not the raw cell count; "unscreened" is
// the plain mean (MaskingScheme::None always has weight 1, so no separate
// code path is needed).
//
// `state` must already have at least 2 valid ghost cells (state.
// FillBoundary) before calling this -- FourthOrderDerivatives::diff1 uses a
// 2-cell-wide stencil.

#include "Energy.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Masking.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

struct TotalEnergyResult
{
    double rho_tot_unscreened{0.0};
    double rho_tot_screened{0.0};
    double rho_axion_kinetic_unscreened{0.0};
    double rho_axion_kinetic_screened{0.0};
    double rho_axion_gradient_unscreened{0.0};
    double rho_axion_gradient_screened{0.0};
    double n_total{0.0};
    double n_unmasked{0.0}; // sum of masking_weight over all cells
};

[[nodiscard]] inline TotalEnergyResult
compute_total_energy(const amrex::MultiFab &state, amrex::Real dx,
                     amrex::Real R, amrex::Real lambda, amrex::Real b_inv,
                     amrex::Real tau, const MaskingParams &screening,
                     long n_cells_total)
{
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                    amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
                    amrex::ReduceOpSum, amrex::ReduceOpSum>
        reduce_op;
    amrex::ReduceData<double, double, double, double, double, double, double>
        reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const auto &arrs = state.const_arrays();
    const FourthOrderDerivatives deriv(dx);

    reduce_op.eval(
        state, amrex::IntVect(0), reduce_data,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> ReduceTuple
        {
            const auto &a = arrs[box_no];
            const double psi1 = a(i, j, k, c_psi1);
            const double psi2 = a(i, j, k, c_psi2);
            const double Pi1  = a(i, j, k, c_Pi1);
            const double Pi2  = a(i, j, k, c_Pi2);
            const double psi_sq = psi1 * psi1 + psi2 * psi2;

            // psi_sq == 0.0 exactly is a genuine 0/0 for both quantities
            // below (grad_theta and theta_prime's numerators are also
            // exactly 0 there, since both are linear in psi1, psi2) --
            // string cores, not a numerical-precision edge case. Guarded
            // the same way Masking.hpp's masked_a_dot already is ("check
            // before dividing, so a 0/0 core is zeroed cleanly rather than
            // propagating a NaN through 0 * NaN") -- EnergyKernel.hpp's
            // masking_weight multiply happens *after* this, so it cannot
            // do that job on its own (found the hard way, 2026-09-18: a
            // genuinely random Fourier IC hits this on the very first
            // step, unlike the T1 test's two cores, deliberately offset
            // off-grid by construction).
            const auto grad_psi1 = deriv.d1_scalar(i, j, k, a, c_psi1);
            const auto grad_psi2 = deriv.d1_scalar(i, j, k, a, c_psi2);
            double grad_psi1_sq  = 0.0;
            double grad_psi2_sq  = 0.0;
            double grad_theta_sq = 0.0;
            FOR (dir)
            {
                grad_psi1_sq += grad_psi1(dir) * grad_psi1(dir);
                grad_psi2_sq += grad_psi2(dir) * grad_psi2(dir);
                const double dtheta_dir =
                    (psi_sq > 0.0)
                        ? (psi1 * grad_psi2(dir) - psi2 * grad_psi1(dir)) /
                              psi_sq
                        : 0.0;
                grad_theta_sq += dtheta_dir * dtheta_dir;
            }

            const double rho = rho_tot_pointwise(
                psi1, psi2, Pi1, Pi2, grad_psi1_sq, grad_psi2_sq, R, lambda,
                b_inv, tau);

            const double theta_prime =
                (psi_sq > 0.0) ? (psi1 * Pi2 - Pi1 * psi2) / psi_sq : 0.0;
            const double rho_a_kin =
                axion_kinetic_energy_pointwise(theta_prime, R);
            const double rho_a_grad =
                axion_gradient_energy_pointwise(grad_theta_sq, R);

            const double w = masking_weight(screening, psi1, psi2, R);

            return {rho,       rho * w,       rho_a_kin, rho_a_kin * w,
                   rho_a_grad, rho_a_grad * w, w};
        });

    const ReduceTuple result = reduce_data.value(reduce_op);
    // amrex::ReduceOps only reduces within this rank's own boxes -- found
    // the hard way (2026-09-18): with a genuinely multi-rank run, every
    // rank silently reported only its own local sums as if they were the
    // whole domain's, without this. AMReX's own Reduce::Sum/Min/Max free
    // functions (AMReX_Reduce.H) never add this either -- it is always
    // the caller's job for a cross-rank total.
    double sums[7] = {
        amrex::get<0>(result), amrex::get<1>(result), amrex::get<2>(result),
        amrex::get<3>(result), amrex::get<4>(result), amrex::get<5>(result),
        amrex::get<6>(result)};
    amrex::ParallelDescriptor::ReduceRealSum(sums, 7);
    const double sum_rho          = sums[0];
    const double sum_rho_w        = sums[1];
    const double sum_rho_a_kin    = sums[2];
    const double sum_rho_a_kin_w  = sums[3];
    const double sum_rho_a_grad   = sums[4];
    const double sum_rho_a_grad_w = sums[5];
    const double sum_w            = sums[6];

    const auto n_total_d = static_cast<double>(n_cells_total);

    TotalEnergyResult out{};
    out.n_total    = n_total_d;
    out.n_unmasked = sum_w;
    out.rho_tot_unscreened           = sum_rho / n_total_d;
    out.rho_axion_kinetic_unscreened = sum_rho_a_kin / n_total_d;
    out.rho_axion_gradient_unscreened = sum_rho_a_grad / n_total_d;
    if (sum_w > 0.0)
    {
        out.rho_tot_screened            = sum_rho_w / sum_w;
        out.rho_axion_kinetic_screened  = sum_rho_a_kin_w / sum_w;
        out.rho_axion_gradient_screened = sum_rho_a_grad_w / sum_w;
    }
    return out;
}

#endif // ENERGYKERNEL_HPP_
