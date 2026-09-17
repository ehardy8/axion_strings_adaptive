#ifndef ENERGYKERNEL_HPP_
#define ENERGYKERNEL_HPP_

// AMReX-side energy reduction: computes the screened and unscreened
// spatial averages of rho_tot (Energy.hpp) and the axion kinetic energy,
// via a single grid pass. "Screened" divides by the effective point count
// (conventions.md sec.10: "the effective number of lattice points averaged
// over must be multiplied by the same masking factor"), not the raw cell
// count; "unscreened" is the plain mean (MaskingScheme::None always has
// weight 1, so no separate code path is needed -- the masking module
// already supports f=1 as a first-class mode, milestone-1.md task 1.8's
// requirement of task 1.7).
//
// `state` must already have at least 2 valid ghost cells (state.
// FillBoundary) before calling this -- FourthOrderDerivatives::diff1 uses a
// 2-cell-wide stencil.

#include "Energy.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Masking.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

struct TotalEnergyResult
{
    double rho_tot_unscreened{0.0};
    double rho_tot_screened{0.0};
    double rho_axion_kinetic_unscreened{0.0};
    double rho_axion_kinetic_screened{0.0};
};

[[nodiscard]] inline TotalEnergyResult
compute_total_energy(const amrex::MultiFab &state, amrex::Real dx,
                     amrex::Real R, amrex::Real lambda, amrex::Real b_inv,
                     amrex::Real tau, const MaskingParams &screening,
                     long n_cells_total)
{
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                    amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum>
        reduce_op;
    amrex::ReduceData<double, double, double, double, double> reduce_data(
        reduce_op);
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
                psi1, psi2, Pi1, Pi2, grad_psi1_sq, grad_psi2_sq, R, lambda,
                b_inv, tau);

            const double theta_prime = (psi1 * Pi2 - Pi1 * psi2) /
                                      (psi1 * psi1 + psi2 * psi2);
            const double rho_a_kin =
                axion_kinetic_energy_pointwise(theta_prime, R);

            const double w = masking_weight(screening, psi1, psi2, R);

            return {rho, rho * w, rho_a_kin, rho_a_kin * w, w};
        });

    const ReduceTuple result = reduce_data.value(reduce_op);
    const double sum_rho          = amrex::get<0>(result);
    const double sum_rho_w        = amrex::get<1>(result);
    const double sum_rho_a_kin    = amrex::get<2>(result);
    const double sum_rho_a_kin_w  = amrex::get<3>(result);
    const double sum_w            = amrex::get<4>(result);

    TotalEnergyResult out{};
    out.rho_tot_unscreened          = sum_rho / static_cast<double>(n_cells_total);
    out.rho_axion_kinetic_unscreened = sum_rho_a_kin / static_cast<double>(n_cells_total);
    if (sum_w > 0.0)
    {
        out.rho_tot_screened           = sum_rho_w / sum_w;
        out.rho_axion_kinetic_screened = sum_rho_a_kin_w / sum_w;
    }
    return out;
}

#endif // ENERGYKERNEL_HPP_
