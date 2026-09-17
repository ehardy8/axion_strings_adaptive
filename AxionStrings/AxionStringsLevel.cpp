#include "AxionStringsLevel.hpp"
#include "AxionStringsParams.hpp"
#include "AxionStringsRHS.hpp"
#include "EnergyKernel.hpp"
#include "FixedGridsTagger.hpp"
#include "FourierIC.hpp"
#include "MaskedFieldBuffer.hpp"
#include "SmallDataIO.hpp"
#include "SpectrumKernel.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"
#include "StringFinder.hpp"
#include "XiFormula.hpp"

void AxionStringsLevel::variableSetUp()
{
    BL_PROFILE("AxionStringsLevel::variableSetUp()");

    state_variable_set_up();

    s_mode       = AxionStringsParams::read_mode();
    s_background = AxionStringsParams::read_background();
    s_tau_i      = AxionStringsParams::read_tau_i();

    if (s_mode == AxionStringsParams::Mode::PreEvolution)
    {
        s_pre_background = AxionStringsParams::read_pre_evolution_background();
        s_xi_target       = AxionStringsParams::read_xi_target();
        s_xi_cadence       = AxionStringsParams::read_xi_check_cadence();
        s_xi_check_interval = s_xi_cadence.coarse_interval;
    }
    else
    {
        s_energy_masking =
            AxionStringsParams::read_masking_params("axion_strings.masking");
    }
}

void AxionStringsLevel::initData()
{
    BL_PROFILE("AxionStringsLevel::initData()");

    std::string ic_mode = "homogeneous";
    amrex::ParmParse("axion_strings").query("ic_mode", ic_mode);

    if (ic_mode == "plane_wave_test")
    {
        // T3 (conventions.md sec.13/milestone-1.md task 1.9): a small-
        // amplitude axion plane wave, no strings, for the spectral
        // normalisation end-to-end check. psi1=R, Pi1=0 (no radial
        // excitation); Pi2 = C sin(2 pi p0 x/L) gives a known single-mode
        // theta' = (psi1 Pi2 - Pi1 psi2)/|psi|^2 = Pi2/R = (C/R) sin(...)
        // directly (psi2=0 keeps this exact, not just small-amplitude).
        const amrex::Real tau = s_tau_i;
        const amrex::Real R_i = s_background.R(tau);

        int p0 = 4;
        amrex::ParmParse("axion_strings").queryAdd("plane_wave_p0", p0);
        amrex::Real amplitude_C = 0.1;
        amrex::ParmParse("axion_strings")
            .queryAdd("plane_wave_amplitude", amplitude_C);

        const auto dx         = Geom().CellSizeArray();
        const auto prob_lo    = Geom().ProbLoArray();
        const amrex::Real L_x = Geom().ProbLength(0);
        const amrex::Real k0  = 2.0 * M_PI * p0 / L_x;

        amrex::MultiFab &state_new = get_new_data(state_index);
        auto const &arrs           = state_new.arrays();

        amrex::ParallelFor(
            state_new, state_new.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                arrs[box_no](i, j, k, c_psi1) = R_i;
                arrs[box_no](i, j, k, c_psi2) = 0.0;
                arrs[box_no](i, j, k, c_Pi1)  = 0.0;
                arrs[box_no](i, j, k, c_Pi2) =
                    amplitude_C * std::sin(k0 * x);
            });
        amrex::Gpu::streamSynchronize();

        // Check the spectrum of the exact IC just set, before any
        // evolution changes it (specific_post_timestep only fires after a
        // step -- by then R(tau) itself has moved on, so the field is no
        // longer this exact, by-hand-checkable configuration).
        bool compute_spectrum_flag = false;
        amrex::ParmParse("axion_strings")
            .queryAdd("compute_spectrum", compute_spectrum_flag);
        if (compute_spectrum_flag)
        {
            MaskingParams none_masking{};
            none_masking.scheme = MaskingScheme::None;
            amrex::MultiFab a_dot_buffer(state_new.boxArray(),
                                         state_new.DistributionMap(), 1, 0);
            fill_masked_a_dot_buffer(a_dot_buffer, state_new, none_masking,
                                     R_i);
            const Spectrum spectrum = compute_spectrum(a_dot_buffer, Geom());

            const double n_total_d =
                static_cast<double>(Geom().Domain().numPts());
            const double sum_sq = static_cast<double>(amrex::MultiFab::Dot(
                a_dot_buffer, 0, a_dot_buffer, 0, 1, 0));
            const double real_space_mean_sq = sum_sq / n_total_d;
            const double n_total_sq         = n_total_d * n_total_d;

            // S(p) = 4*pi<p^2|X_p|^2> is built from raw (un-Parseval-
            // normalised) |X_p|^2, same as full_cube/inscribed_sphere --
            // needs the same /n_total_sq. Approximate, not exact, even
            // after that: a single anisotropic mode shares its shell with
            // other (zero-power) lattice points at the same |p|, which
            // dilutes the shell average -- unlike the full-cube/inscribed-
            // sphere sums, which include every mode exactly once and so
            // match the real-space value exactly.
            double shell_integral = 0.0;
            for (double s : spectrum.shell_average)
            {
                shell_integral += s;
            }
            shell_integral /= n_total_sq;

            amrex::Print()
                << "  [T3] known amplitude_C = " << amplitude_C
                << ", p0 = " << p0
                << ": expected <a_dot^2> = " << 0.5 * amplitude_C * amplitude_C
                << "\n    measured <a_dot^2> real-space = "
                << real_space_mean_sq
                << "  Parseval (full cube) = "
                << spectrum.full_cube_energy / n_total_sq
                << "  Parseval (inscribed sphere) = "
                << spectrum.inscribed_sphere_energy / n_total_sq << "\n"
                << "    shell-binned integral (sum_s S(s)/n_total^2, "
                   "approximate) = "
                << shell_integral << "\n";
        }
        return;
    }

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
        // on a (phase-undefined) core.
        //
        // Radial profile g(rho) = rho_hat/sqrt(rho_hat^2+2), rho_hat =
        // m_r*rho: an approximation to sec.13's true equilibrium profile,
        // not exact (its near-core slope is 1/sqrt(2) =~ 0.707, not the
        // exact c1 = 0.41222), but with the *same functional form* at both
        // ends -- linear near the core, and g = 1/sqrt(1+2/rho_hat^2)
        // =~ 1 - 1/rho_hat^2 far away, matching sec.13's stated "1 - rho^-2"
        // asymptote exactly. This matters: an earlier attempt used tanh,
        // which decays to vacuum *exponentially* rather than as a power
        // law, and so is missing essentially all of the long-range
        // Goldstone tail that is responsible for the string's logarithmic
        // tension divergence in the first place -- for task 1.8's tension
        // test specifically, that tail is the point, not a detail. tanh
        // also turned out to be far enough from equilibrium that the field
        // "rang down" violently within a handful of steps (confirmed this
        // wasn't a deeper bug: the exact homogeneous solution stays at
        // exactly zero energy indefinitely under the same evolution code).
        const amrex::Real tau      = s_tau_i;
        const amrex::Real R_i      = s_background.R(tau);
        const amrex::Real m_r_test = std::sqrt(s_background.lambda(tau));
        const amrex::Real bkg_b_inv = s_background.b_inv;

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
                const amrex::Real rho1_hat = m_r_test * rho1;
                const amrex::Real rho2_hat = m_r_test * rho2;
                const amrex::Real g1 =
                    rho1_hat / std::sqrt(rho1_hat * rho1_hat + 2.0);
                const amrex::Real g2 =
                    rho2_hat / std::sqrt(rho2_hat * rho2_hat + 2.0);
                const amrex::Real amp  = R_i * g1 * g2;
                const amrex::Real psi1 = amp * std::cos(theta);
                const amrex::Real psi2 = amp * std::sin(theta);

                arrs[box_no](i, j, k, c_psi1) = psi1;
                arrs[box_no](i, j, k, c_psi2) = psi2;
                // Pi = (R'/R) psi, NOT 0: a genuinely static (non-radiating)
                // profile has a fixed *shape*, but psi's overall amplitude
                // must still track the growing background R(tau) -- exactly
                // like the homogeneous vacuum solution psi=R(tau) (tasks
                // 1.2/1.3). Setting Pi=0 here was a real bug: it left the
                // far field's |psi|/R decaying as R_i/R(tau) regardless of
                // any string, dropping below a typical masking threshold
                // within just 2-3 steps and collapsing the screened energy
                // to zero -- confirmed numerically (R_i/R(1.3) = 0.77 < 0.8)
                // to be exactly this effect, not a deeper bug (the exact
                // homogeneous solution stays at exactly zero energy
                // indefinitely under the same evolution code).
                const amrex::Real R_prime_over_R = 1.0 / (bkg_b_inv * tau);
                arrs[box_no](i, j, k, c_Pi1) = R_prime_over_R * psi1;
                arrs[box_no](i, j, k, c_Pi2) = R_prime_over_R * psi2;
            });
        amrex::Gpu::streamSynchronize();
        return;
    }

    const bool generating_pre_evolution_ic =
        (s_mode == AxionStringsParams::Mode::PreEvolution);

    if (generating_pre_evolution_ic || ic_mode == "fourier")
    {
        // Fourier-mode initial conditions (conventions.md sec.7): modes
        // occupied for |k| <= k_max, zero above, normalised to a target
        // mean-square variance. Pi1 = Pi2 = 0 (displacement-only ICs --
        // conventions.md does not specify a Pi variance, only psi's).
        const std::string prefix = generating_pre_evolution_ic
                                       ? "axion_strings.pre_evolution"
                                       : "axion_strings";
        const auto fourier_params =
            AxionStringsParams::read_fourier_ic_params(prefix);

        const amrex::Real m_r =
            generating_pre_evolution_ic
                ? std::sqrt(s_pre_background.lambda(0.0))
                : std::sqrt(s_background.lambda(s_tau_i));
        const amrex::Real dx = Geom().CellSize(0);
        const double k_max_cells = fourier_params.k_max_over_mr * m_r * dx;

        amrex::MultiFab &state_new = get_new_data(state_index);
        amrex::MultiFab psi1_mf(state_new.boxArray(),
                                state_new.DistributionMap(), 1, 0);
        amrex::MultiFab psi2_mf(state_new.boxArray(),
                                state_new.DistributionMap(), 1, 0);
        FourierIC::generate(psi1_mf, Geom(), k_max_cells,
                           fourier_params.mean_square_variance,
                           fourier_params.seed, 0);
        FourierIC::generate(psi2_mf, Geom(), k_max_cells,
                           fourier_params.mean_square_variance,
                           fourier_params.seed, 1);

        amrex::MultiFab::Copy(state_new, psi1_mf, 0, c_psi1, 1, 0);
        amrex::MultiFab::Copy(state_new, psi2_mf, 0, c_psi2, 1, 0);
        state_new.setVal(0.0, c_Pi1, 2, 0);
        state_new.FillBoundary(Geom().periodicity());
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

    amrex::Real curvature_coeff{};
    amrex::Real lambda{};
    amrex::Real R_squared{};
    if (s_mode == AxionStringsParams::Mode::PreEvolution)
    {
        // Pre-evolution's own clock starts at a_time = 0 = tau_pre.
        const amrex::Real tau_pre = a_time;
        curvature_coeff = s_pre_background.curvature_term_coeff(tau_pre);
        lambda           = s_pre_background.lambda(tau_pre);
        const amrex::Real R_val = s_pre_background.R(tau_pre);
        R_squared               = R_val * R_val;
    }
    else
    {
        const amrex::Real tau = s_tau_i + a_time;
        curvature_coeff        = s_background.curvature_term_coeff(tau);
        lambda                 = s_background.lambda(tau);
        const amrex::Real R_val = s_background.R(tau);
        R_squared               = R_val * R_val;
    }

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

    if (s_mode == AxionStringsParams::Mode::PreEvolution)
    {
        // xi-monitoring stopping loop (conventions.md sec.7, task 1.5).
        // Measuring xi means a plaquette count over the whole grid, which
        // is not cheap -- check on an adaptive cadence: coarse while xi is
        // far from the target, tightening as it gets close (user request).
        ++s_steps_since_xi_check;
        if (s_steps_since_xi_check < s_xi_check_interval)
        {
            return;
        }
        s_steps_since_xi_check = 0;

        amrex::MultiFab &pre_state = get_new_data(state_index);
        pre_state.FillBoundary(Geom().periodicity());
        const PlaquetteCounts counts = count_plaquettes(pre_state);

        const amrex::Real tau_pre = get_state_data(state_index).curTime();
        const double dx           = Geom().CellSize(0);
        const double L_tilde      = Geom().ProbLength(0);
        // xi is measured at tau = tau_i (the main run's start -- the
        // physical instant this pre-evolution state represents once
        // handed off), not tau_pre (confirmed with the user).
        const double xi = xi_from_plaquette_count(
            static_cast<double>(counts.n_p_plain), dx, L_tilde,
            s_background.a_inv, s_tau_i);
        const double ratio = xi / s_xi_target;

        amrex::Print() << "  [AxionStrings pre-evolution] tau_pre = "
                       << tau_pre << "  N_p = " << counts.n_p_plain
                       << "  xi(at tau_i) = " << xi
                       << "  target = " << s_xi_target
                       << "  ratio = " << ratio << "\n";

        // xi generally decreases during relaxation (user's observation):
        // stop the first time it drops to the target rather than waiting
        // for an exact match.
        if (ratio <= 1.0)
        {
            amrex::Print()
                << "  [AxionStrings pre-evolution] target xi reached -- "
                   "stopping (a checkpoint is written once the time-"
                   "stepping loop exits).\n";
            s_pre_evolution_target_reached = true;
            return;
        }

        if (ratio <= s_xi_cadence.fine_threshold)
        {
            s_xi_check_interval = s_xi_cadence.fine_interval;
        }
        else if (ratio <= s_xi_cadence.medium_threshold)
        {
            s_xi_check_interval = s_xi_cadence.medium_interval;
        }
        else
        {
            s_xi_check_interval = s_xi_cadence.coarse_interval;
        }
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

    // Energy (conventions.md sec.10, milestone-1.md task 1.8): rho_tot
    // (the sec.10 aggregate formula only -- the sec.12 radial/axion/
    // interaction split is not yet implemented, see docs/STATUS.md) and the
    // axion kinetic energy, both screened and unscreened. FourthOrder-
    // Derivatives::diff1 needs a 2-cell-wide stencil; the FillBoundary
    // above already covers it (state's ghost count is >=3 by default, for
    // the 4th-order Laplacian in the RHS).
    const amrex::Real lambda = s_background.lambda(tau);
    const long n_cells_total = Geom().Domain().numPts();
    const TotalEnergyResult energy =
        compute_total_energy(state_new, dx, s_background.R(tau), lambda,
                             s_background.b_inv, tau, s_energy_masking,
                             n_cells_total);

    // Core-energy diagnostics (2026-09-17, with the user): the screened/
    // unscreened *difference* of averages is diluted by the masked
    // point-count fraction (n_unmasked/n_total) and is not itself "the
    // string energy". The energy actually sitting in the masked (core)
    // cells is instead Sum_core rho = Sum_all rho - Sum_unmasked rho =
    // n_total*unscreened_avg - n_unmasked*screened_avg. For a *static*
    // string, adding the away-from-core (axion_gradient - axion_kinetic)
    // contribution approximately removes propagating axion wave energy
    // too: a free wave has equal kinetic/gradient energy on average, so
    // that difference cancels the wave's contribution and leaves only the
    // string's own static long-range tail. Printed here as a convenience,
    // interactive-use diagnostic only -- not saved to network_scalars.dat,
    // since the point of saving every raw screened/unscreened component
    // (below) is to let this and any other combination be reconstructed
    // afterwards without having picked one formula in advance.
    const double sum_core_rho_tot =
        energy.n_total * energy.rho_tot_unscreened -
        energy.n_unmasked * energy.rho_tot_screened;
    const double sum_unmasked_tail =
        energy.n_unmasked * (energy.rho_axion_gradient_screened -
                             energy.rho_axion_kinetic_screened);
    const double dx3 = dx * dx * dx;
    // 2 strings (vortex/antivortex pair), each spanning the full box in z.
    const double string_length_in_box = 2.0 * Geom().ProbLength(2);
    const double tension_core_only =
        sum_core_rho_tot * dx3 / string_length_in_box;
    const double tension_core_plus_tail =
        (sum_core_rho_tot + sum_unmasked_tail) * dx3 / string_length_in_box;

    amrex::Print() << "  [AxionStrings] tau = " << tau
                   << "  N_p = " << counts.n_p_plain
                   << "  N_p_W = " << counts.n_p_weighted
                   << "  xi = " << xi_plain << "  xi_W = " << xi_weighted
                   << "  m_r/H (direct) = " << m_r_over_h
                   << "  m_r/H (closed form) = " << 1.0 / h_over_mr_closed
                   << "\n"
                   << "    rho_tot: unscreened = " << energy.rho_tot_unscreened
                   << "  screened = " << energy.rho_tot_screened << "\n"
                   << "    rho_axion_kin: unscreened = "
                   << energy.rho_axion_kinetic_unscreened
                   << "  screened = " << energy.rho_axion_kinetic_screened
                   << "\n"
                   << "    rho_axion_grad: unscreened = "
                   << energy.rho_axion_gradient_unscreened
                   << "  screened = " << energy.rho_axion_gradient_screened
                   << "\n"
                   << "    n_total = " << energy.n_total
                   << "  n_unmasked = " << energy.n_unmasked
                   << "  tension (core only) = " << tension_core_only
                   << "  tension (core+tail) = " << tension_core_plus_tail
                   << "\n";

    const amrex::Real dt =
        a_time_now - get_state_data(state_index).prevTime();
    SmallDataIO network_scalars_file("network_scalars", dt, tau, 0.0,
                                     SmallDataIO::APPEND);
    if (!s_wrote_network_scalars_header)
    {
        network_scalars_file.write_header_line(
            {"N_p", "N_p_weighted", "xi", "xi_weighted", "m_r_over_H",
             "rho_tot_unscreened", "rho_tot_screened",
             "rho_axion_kin_unscreened", "rho_axion_kin_screened",
             "rho_axion_grad_unscreened", "rho_axion_grad_screened",
             "n_total", "n_unmasked"});
        s_wrote_network_scalars_header = true;
    }
    const std::vector<amrex::Real> data_row{
        static_cast<amrex::Real>(counts.n_p_plain),
        static_cast<amrex::Real>(counts.n_p_weighted),
        static_cast<amrex::Real>(xi_plain),
        static_cast<amrex::Real>(xi_weighted),
        m_r_over_h,
        static_cast<amrex::Real>(energy.rho_tot_unscreened),
        static_cast<amrex::Real>(energy.rho_tot_screened),
        static_cast<amrex::Real>(energy.rho_axion_kinetic_unscreened),
        static_cast<amrex::Real>(energy.rho_axion_kinetic_screened),
        static_cast<amrex::Real>(energy.rho_axion_gradient_unscreened),
        static_cast<amrex::Real>(energy.rho_axion_gradient_screened),
        static_cast<amrex::Real>(energy.n_total),
        static_cast<amrex::Real>(energy.n_unmasked)};
    network_scalars_file.write_time_data_line(data_row);

    // Spectrum (conventions.md sec.9/sec.12, milestone-1.md task 1.9):
    // opt-in (computing an FFT every snapshot is not free), gated behind
    // axion_strings.compute_spectrum. Masking is applied at exactly one
    // place -- MaskedFieldBuffer.hpp, filling the buffer handed to the FFT
    // -- per CLAUDE.md constraint 5.
    bool compute_spectrum_flag = false;
    amrex::ParmParse("axion_strings")
        .queryAdd("compute_spectrum", compute_spectrum_flag);
    if (compute_spectrum_flag)
    {
        amrex::MultiFab a_dot_buffer(state_new.boxArray(),
                                     state_new.DistributionMap(), 1, 0);
        fill_masked_a_dot_buffer(a_dot_buffer, state_new, s_energy_masking,
                                 s_background.R(tau));
        const Spectrum spectrum = compute_spectrum(a_dot_buffer, Geom());

        // Direct real-space mean square of the exact same buffer that was
        // fed to the FFT -- the honest way to get <(masked a_dot)^2>,
        // rather than trying to back it out of axion_kinetic_energy's own
        // f_a^2/2 prefactors (error-prone; avoided deliberately).
        const double sum_sq = static_cast<double>(amrex::MultiFab::Dot(
            a_dot_buffer, 0, a_dot_buffer, 0, 1, 0));
        const double real_space_mean_sq = sum_sq / energy.n_total;

        const double n_total_sq = energy.n_total * energy.n_total;
        // Parseval (conventions.md sec.9, verified with the user): a
        // spatial *average* of a squared real field equals
        // Sum_p|X_p|^2 / N_total^2 (one factor of N_total from the
        // average itself, one from Parseval for this "no normalisation,
        // round-trip = N_total" FFT convention) -- exact, independent of
        // how the shell-binning below groups modes.
        const double parseval_from_spectrum_full =
            spectrum.full_cube_energy / n_total_sq;
        const double parseval_from_spectrum_inscribed =
            spectrum.inscribed_sphere_energy / n_total_sq;

        // S(p) needs the same /n_total_sq as the two Parseval sums above
        // (it too is built from raw |X_p|^2); even then this is only
        // approximate, not exact, since a single shell can mix modes with
        // very different power (unlike the full-cube/inscribed-sphere
        // sums, which include every mode exactly once).
        double shell_integral = 0.0;
        for (double s : spectrum.shell_average)
        {
            shell_integral += s; // integral over |k|, bin width 1
        }
        shell_integral /= n_total_sq;

        amrex::Print()
            << "  [AxionStrings spectrum] <a_dot^2> real-space = "
            << real_space_mean_sq
            << "  Parseval (full cube) = " << parseval_from_spectrum_full
            << "  Parseval (inscribed sphere) = "
            << parseval_from_spectrum_inscribed << "\n"
            << "    full_cube/inscribed_sphere ratio = "
            << (spectrum.full_cube_energy / spectrum.inscribed_sphere_energy)
            << "  shell-binned integral (sum_s S(s)/n_total^2, approximate) = "
            << shell_integral << "\n";
    }
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

int AxionStringsLevel::okToContinue()
{
    if (s_mode == AxionStringsParams::Mode::PreEvolution &&
        s_pre_evolution_target_reached)
    {
        return 0;
    }
    return 1;
}

void AxionStringsLevel::specific_post_restart()
{
    // Pre-evolution -> main handoff (conventions.md sec.7): only on the one
    // restart command that performs it, flagged explicitly so a later,
    // ordinary restart of the main run's own progress does not re-apply
    // the rescale (see AxionStringsParams::read_restart_from_pre_evolution).
    if (s_mode != AxionStringsParams::Mode::Main ||
        !AxionStringsParams::read_restart_from_pre_evolution())
    {
        return;
    }

    // Pre-evolution's own clock (tau_pre) starts at 0, so whatever a_time
    // the checkpoint recorded *is* tau_pre at handoff.
    const amrex::Real tau_pre_end = get_state_data(state_index).curTime();

    const amrex::Real R_pre        = s_pre_background.R(tau_pre_end);
    const amrex::Real R_pre_prime  = s_pre_background.R_prime(tau_pre_end);
    const amrex::Real R_main       = s_background.R(s_tau_i);
    const amrex::Real R_main_prime = s_background.R_prime(s_tau_i);
    const amrex::Real kappa        = R_main / R_pre;
    const amrex::Real common_term  = R_main_prime - kappa * kappa * R_pre_prime;

    if (Level() == 0)
    {
        amrex::Print() << "  [AxionStrings] pre-evolution -> main handoff: "
                       << "tau_pre_end = " << tau_pre_end
                       << ", R_pre = " << R_pre
                       << ", R_main(tau_i) = " << R_main
                       << ", kappa = " << kappa << "\n";
    }

    // a_time continues counting up from wherever the checkpoint's clock
    // left off (it does not reset to 0), so s_tau_i must absorb that
    // offset for tau = s_tau_i + a_time to still equal tau_i right now.
    s_tau_i = s_tau_i - tau_pre_end;

    for (amrex::MultiFab *mf :
        {&get_new_data(state_index), &get_old_data(state_index)})
    {
        auto const &arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                auto &a                    = arrs[box_no];
                const amrex::Real psi1_pre = a(i, j, k, c_psi1);
                const amrex::Real psi2_pre = a(i, j, k, c_psi2);
                const amrex::Real Pi1_pre  = a(i, j, k, c_Pi1);
                const amrex::Real Pi2_pre  = a(i, j, k, c_Pi2);

                a(i, j, k, c_psi1) = kappa * psi1_pre;
                a(i, j, k, c_psi2) = kappa * psi2_pre;
                a(i, j, k, c_Pi1) =
                    kappa * kappa * Pi1_pre + (psi1_pre / R_pre) * common_term;
                a(i, j, k, c_Pi2) =
                    kappa * kappa * Pi2_pre + (psi2_pre / R_pre) * common_term;
            });
    }
    amrex::Gpu::streamSynchronize();
}
