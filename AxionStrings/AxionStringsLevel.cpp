#include "AxionStringsLevel.hpp"
#include "AxionStringsParams.hpp"
#include "AxionStringsRHS.hpp"
#include "FixedGridsTagger.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"

void AxionStringsLevel::variableSetUp()
{
    BL_PROFILE("AxionStringsLevel::variableSetUp()");

    state_variable_set_up();

    s_background = AxionStringsParams::read_background();
    s_tau_i      = AxionStringsParams::read_tau_i();
}

void AxionStringsLevel::initData()
{
    BL_PROFILE("AxionStringsLevel::initData()");

    // Placeholder initial data (superseded by the Fourier-mode + pre-
    // evolution generator of milestone-1.md task 1.5): the homogeneous
    // configuration psi1 = R(tau), psi2 = 0. Since |psi|^2 - R^2 vanishes
    // identically here, this is an exact solution of the full nonlinear EOM
    // (see tests/test_background.cpp), giving a non-trivial but analytically
    // known trajectory to exercise the RHS and the c(tau)/lambda(tau)
    // machinery against.
    const amrex::Real tau  = s_tau_i;
    const amrex::Real R_i  = s_background.R(tau);
    const amrex::Real Pi_i = s_background.R_prime(tau);

    amrex::MultiFab &state_new = get_new_data(state_index);
    auto const &arrs           = state_new.arrays();

    amrex::ParallelFor(
        state_new, state_new.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            arrs[box_no](i, j, k, c_psi1) = R_i;
            arrs[box_no](i, j, k, c_psi2) = 0.0;
            arrs[box_no](i, j, k, c_Pi1)  = Pi_i;
            arrs[box_no](i, j, k, c_Pi2)  = 0.0;
        });
    amrex::Gpu::streamSynchronize();
}

void AxionStringsLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                         amrex::MultiFab &a_rhs,
                                         const amrex::Real a_time)
{
    BL_PROFILE("AxionStringsLevel::specific_eval_rhs()");

    const amrex::Real tau = s_tau_i + a_time;
    const amrex::Real curvature_coeff =
        s_background.curvature_term_coeff(tau);
    const amrex::Real lambda    = s_background.lambda(tau);
    const amrex::Real R_val     = s_background.R(tau);
    const amrex::Real R_squared = R_val * R_val;

    const auto dx                 = Geom().CellSize(0);
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    AxionStringsRHS<> rhs(dx, curvature_coeff, lambda, R_squared);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { rhs(ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]); });

    amrex::Gpu::streamSynchronize();
}

void AxionStringsLevel::specific_post_timestep()
{
    // Round-off cross-check for milestone-1 task 1.3's "done when": m_r/H
    // measured from the code (H_over_mr_direct, from R, R' and lambda as
    // the RHS actually evaluates them) against the conventions.md sec.5
    // closed form. Level 0 only to avoid duplicate prints (single-level
    // throughout milestone 1 in any case).
    if (Level() != 0)
    {
        return;
    }

    const amrex::Real tau = s_tau_i + get_state_data(state_index).curTime();
    const amrex::Real h_over_mr_direct = s_background.H_over_mr_direct(tau);
    const amrex::Real h_over_mr_closed =
        s_background.H_over_mr_closed_form(tau);

    amrex::Print() << "  [AxionStrings] tau = " << tau
                   << "  m_r/H (direct) = " << 1.0 / h_over_mr_direct
                   << "  m_r/H (closed form) = " << 1.0 / h_over_mr_closed
                   << "\n";
}

void AxionStringsLevel::tag_cells(amrex::TagBoxArray &tags,
                                 amrex::Real a_regrid_threshold)
{
    // Placeholder tagger: no refinement in milestone 1 (amr.max_level = 0
    // throughout), replaced by the string-based tagger in a later milestone.
    BL_PROFILE("AxionStringsLevel::tag_cells()");

    const auto &tag_arrs = tags.arrays();

    const amrex::Real dx    = Geom().CellSize(0);
    const int current_level = Level();

    FixedGridsTagger my_tagging_criterion{dx, current_level};

    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { my_tagging_criterion(ix, iy, iz, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}
