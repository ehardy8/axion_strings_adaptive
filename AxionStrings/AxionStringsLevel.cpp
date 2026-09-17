#include "AxionStringsLevel.hpp"
#include "AxionStringsParams.hpp"
#include "AxionStringsRHS.hpp"
#include "FixedGridsTagger.hpp"
#include "SmallDataIO.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"
#include "StringFinder.hpp"
#include "XiFormula.hpp"

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

    std::string ic_mode = "homogeneous";
    amrex::ParmParse("axion_strings").query("ic_mode", ic_mode);

    if (ic_mode == "straight_string_test")
    {
        // T1 (conventions.md sec.13/milestone-1.md task 1.6, plaquette-
        // detection slice only -- not the full core-profile/tension test,
        // which needs tasks 1.7/1.8 too): a straight string/antistring pair
        // along z, windings +1 and -1. A *single* isolated vortex cannot be
        // embedded in a periodic box (the net winding over a closed 2-torus
        // must vanish; a naive single-vortex field creates spurious extra
        // windings along the periodic seam where it fails to match itself)
        // -- a compensating antivortex, well separated from both the vortex
        // and the box edges, is what a periodic domain actually requires.
        // Offset by half a cell from grid points so no vertex sits exactly
        // on a (phase-undefined) core. The radial profile is a plain tanh,
        // not the exact equilibrium profile of sec.13 -- irrelevant here,
        // since plaquette detection depends only on the phase.
        const amrex::Real tau      = s_tau_i;
        const amrex::Real R_i      = s_background.R(tau);
        const amrex::Real m_r_test = std::sqrt(s_background.lambda(tau));

        std::array<amrex::Real, AMREX_SPACEDIM> center{};
        amrex::ParmParse().get("geometry.center", center);
        const auto dx      = Geom().CellSizeArray();
        const auto prob_lo = Geom().ProbLoArray();
        const amrex::Real L_x = Geom().ProbLength(0);
        const amrex::Real separation = 0.25 * L_x;
        const amrex::Real x1 = center[0] + 0.5 * dx[0] - separation;
        const amrex::Real y1 = center[1] + 0.5 * dx[1];
        const amrex::Real x2 = center[0] + 0.5 * dx[0] + separation;
        const amrex::Real y2 = center[1] + 0.5 * dx[1];

        amrex::MultiFab &state_new = get_new_data(state_index);
        auto const &arrs           = state_new.arrays();

        amrex::ParallelFor(
            state_new, state_new.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

                const amrex::Real rho1 =
                    std::sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
                const amrex::Real rho2 =
                    std::sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));
                const amrex::Real theta =
                    std::atan2(y - y1, x - x1) - std::atan2(y - y2, x - x2);
                const amrex::Real amp = R_i * std::tanh(m_r_test * rho1) *
                                        std::tanh(m_r_test * rho2);

                arrs[box_no](i, j, k, c_psi1) = amp * std::cos(theta);
                arrs[box_no](i, j, k, c_psi2) = amp * std::sin(theta);
                arrs[box_no](i, j, k, c_Pi1)  = 0.0;
                arrs[box_no](i, j, k, c_Pi2)  = 0.0;
            });
        amrex::Gpu::streamSynchronize();
        return;
    }

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
    // Level 0 only: single-level throughout milestone 1, and the network
    // scalars are a global (not per-level) diagnostic in any case.
    if (Level() != 0)
    {
        return;
    }

    amrex::MultiFab &state_new = get_new_data(state_index);
    // count_plaquettes reads the (i+1,j+1,k+1) neighbours of every valid
    // cell, so needs at least 1 valid ghost cell first.
    state_new.FillBoundary(Geom().periodicity());

    const PlaquetteCounts counts = count_plaquettes(state_new);

    const amrex::Real a_time_now = get_state_data(state_index).curTime();
    const amrex::Real tau        = s_tau_i + a_time_now;
    const double dx              = Geom().CellSize(0);
    const double L_tilde         = Geom().ProbLength(0);
    const double a_inv           = s_background.a_inv;

    // xi from pierced plaquettes (conventions.md sec.8); xi_weighted is the
    // winding-weighted variant (sec.12). Moore's trick needs an extra
    // H0^-2 factor here (sec.8) -- not yet implemented, see docs/STATUS.md.
    const double xi_plain = xi_from_plaquette_count(
        static_cast<double>(counts.n_p_plain), dx, L_tilde, a_inv, tau);
    const double xi_weighted = xi_from_plaquette_count(
        static_cast<double>(counts.n_p_weighted), dx, L_tilde, a_inv, tau);

    // Round-off cross-check from task 1.3: m_r/H measured from the code
    // (H_over_mr_direct, from R, R' and lambda as the RHS actually
    // evaluates them) against the conventions.md sec.5 closed form.
    const amrex::Real h_over_mr_direct = s_background.H_over_mr_direct(tau);
    const amrex::Real h_over_mr_closed =
        s_background.H_over_mr_closed_form(tau);
    const amrex::Real m_r_over_h = 1.0 / h_over_mr_direct;

    amrex::Print() << "  [AxionStrings] tau = " << tau
                   << "  N_p = " << counts.n_p_plain
                   << "  N_p_W = " << counts.n_p_weighted
                   << "  xi = " << xi_plain << "  xi_W = " << xi_weighted
                   << "  m_r/H (direct) = " << m_r_over_h
                   << "  m_r/H (closed form) = " << 1.0 / h_over_mr_closed
                   << "\n";

    const amrex::Real dt =
        a_time_now - get_state_data(state_index).prevTime();
    SmallDataIO network_scalars_file("network_scalars", dt, tau, 0.0,
                                     SmallDataIO::APPEND);
    if (!s_wrote_network_scalars_header)
    {
        network_scalars_file.write_header_line(
            {"N_p", "N_p_weighted", "xi", "xi_weighted", "m_r_over_H"});
        s_wrote_network_scalars_header = true;
    }
    const std::vector<amrex::Real> data_row{
        static_cast<amrex::Real>(counts.n_p_plain),
        static_cast<amrex::Real>(counts.n_p_weighted),
        static_cast<amrex::Real>(xi_plain),
        static_cast<amrex::Real>(xi_weighted), m_r_over_h};
    network_scalars_file.write_time_data_line(data_row);
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
