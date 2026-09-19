#include "AxionStringsLevel.hpp"
#include "AxionStringsParams.hpp"
#include "AxionStringsRHS.hpp"
#include "EnergyKernel.hpp"
#include "FourierIC.hpp"
#include "MaskedFieldBuffer.hpp"
#include "ProjectionKernel.hpp"
#include "SmallDataIO.hpp"
#include "SpectrumKernel.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"
#include "StringFinder.hpp"
#include "StringTagger.hpp"
#include "VelocityKernel.hpp"
#include "XiFormula.hpp"

void AxionStringsLevel::variableSetUp()
{
    BL_PROFILE("AxionStringsLevel::variableSetUp()");

    state_variable_set_up();

    s_background = AxionStringsParams::read_background();
    s_tau_i      = AxionStringsParams::read_tau_i(s_background);

    // Read back AxionStringsParams::apply_box_plan's derived tau_f
    // (informational ParmParse injection; see its own comment) -- absent
    // in Moore mode, where okToContinue() falls back to evolution.
    // stop_time/max_steps instead.
    double tau_f{};
    s_has_tau_f =
        amrex::ParmParse("axion_strings").query("derived_tau_f", tau_f);
    s_tau_f = tau_f;

    // Read back the AMR level-addition schedule (absent for max_level=0
    // and Moore mode -- see AxionStringsLevel.hpp's own comment).
    s_has_level_schedule =
        amrex::ParmParse("axion_strings")
            .queryarr("derived_level_add_log_mr_over_h",
                     s_level_add_log_mr_over_h) != 0;

    // Every run needs these now, not just ones that skip relaxation --
    // there is no longer a separate, non-diagnostic relaxation-only mode.
    s_energy_masking =
        AxionStringsParams::read_masking_params("axion_strings.masking");
    s_output_cadence = AxionStringsParams::read_output_cadence();
    s_next_output_log_mr_over_h = s_output_cadence.first_log_mr_over_h;
    s_tagging_params = AxionStringsParams::read_tagging_params();

    std::string ic_mode = "homogeneous";
    amrex::ParmParse("axion_strings").queryAdd("ic_mode", ic_mode);
    if (ic_mode == "fourier_relaxed")
    {
        s_pre_background = AxionStringsParams::read_pre_evolution_background();
        s_xi_target       = AxionStringsParams::read_xi_target();
        s_xi_cadence       = AxionStringsParams::read_xi_check_cadence();
        s_xi_check_interval = s_xi_cadence.coarse_interval;
        s_phase = Phase::Relaxing;
    }
    else
    {
        s_phase = Phase::Evolving;
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

    const bool generating_relaxed_ic = (ic_mode == "fourier_relaxed");

    if (generating_relaxed_ic || ic_mode == "fourier")
    {
        // Fourier-mode initial conditions (conventions.md sec.7): modes
        // occupied for |k| <= k_max, zero above, normalised to a target
        // mean-square variance. Pi1 = Pi2 = 0 (displacement-only ICs --
        // conventions.md does not specify a Pi variance, only psi's).
        const std::string prefix = generating_relaxed_ic
                                       ? "axion_strings.pre_evolution"
                                       : "axion_strings";
        const auto fourier_params =
            AxionStringsParams::read_fourier_ic_params(prefix);

        const amrex::Real m_r =
            generating_relaxed_ic
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
    if (s_phase == Phase::Relaxing)
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

    if (s_phase == Phase::Relaxing)
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
        // physical instant this relaxed state represents once the
        // transition below occurs), not tau_pre (confirmed with the user).
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
        // transition the first time it drops to the target rather than
        // waiting for an exact match.
        if (ratio <= 1.0)
        {
            amrex::Print()
                << "  [AxionStrings pre-evolution] target xi reached -- "
                   "transitioning in place to the main evolution (no "
                   "restart/checkpoint involved).\n";
            apply_pre_evolution_to_main_rescale();
            s_phase = Phase::Evolving;
            return; // the main-mode diagnostics below start on the next call
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

    const amrex::Real a_time_now = get_state_data(state_index).curTime();
    const amrex::Real tau        = s_tau_i + a_time_now;

    // Output cadence (2026-09-18, with the user): the diagnostics below
    // (a full-grid plaquette count, energy/velocity reductions, an FFT)
    // are far too expensive to repeat every coarse step. m_r/H is
    // available for free from the analytic background (no grid pass), so
    // check it every step and only do the expensive work -- and only
    // write a snapshot -- once log(m_r/H) has reached the next threshold,
    // starting at axion_strings.output_first_log_mr_over_h and spaced by
    // axion_strings.output_delta_log_mr_over_h thereafter. The `while`
    // (not a single add) means a step that jumps past more than one
    // threshold still lands on the correct next one rather than drifting.
    //
    // log(m_r/H) itself is what's *reported* (network_scalars.dat's
    // m_r_over_H column) throughout, including during Moore, where it is
    // correctly constant -- but that means it cannot also drive the
    // cadence *trigger* once in Moore (see AxionStringsLevel.hpp's
    // s_moore_output_initialised comment): switch the trigger coordinate
    // to D_elapsed = log(H(tau_switch)/H(tau)) there instead, same spacing.
    const double log_mr_over_h =
        -std::log(static_cast<double>(s_background.H_over_mr_direct(tau)));
    const bool in_moore =
        s_background.c_sched.has_switch && tau >= s_background.c_sched.tau_switch;
    if (in_moore)
    {
        if (!s_moore_output_initialised)
        {
            s_moore_H_switch               = s_background.H(s_background.c_sched.tau_switch);
            s_next_output_moore_log_range  = 0.0;
            s_moore_output_initialised     = true;
        }
        const double d_elapsed =
            std::log(s_moore_H_switch / s_background.H(tau));
        if (d_elapsed < s_next_output_moore_log_range)
        {
            return;
        }
        while (s_next_output_moore_log_range <= d_elapsed)
        {
            s_next_output_moore_log_range += s_output_cadence.delta_log_mr_over_h;
        }
    }
    else
    {
        if (log_mr_over_h < s_next_output_log_mr_over_h)
        {
            return;
        }
        while (s_next_output_log_mr_over_h <= log_mr_over_h)
        {
            s_next_output_log_mr_over_h += s_output_cadence.delta_log_mr_over_h;
        }
    }

    // count_plaquettes reads the (i+1,j+1,k+1) neighbours of every valid
    // cell, so needs at least 1 valid ghost cell first.
    state_new.FillBoundary(Geom().periodicity());

    const PlaquetteCounts counts = count_plaquettes(state_new);

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

    // String velocities (conventions.md sec.8, milestone-1.md task 1.10 --
    // velocities only, not curvature or loops, both deferred as more
    // involved per the user). gamma^2 v^2 evaluated at the corners of
    // every pierced plaquette, averaged over the network.
    const amrex::Real m_r_now = std::sqrt(lambda);
    const VelocityResult velocity = compute_velocity_at_pierced_corners(
        state_new, s_background.R(tau), tau, s_background.b_inv, m_r_now);
    const double mean_gamma_sq_v_sq =
        (velocity.count > 0)
            ? velocity.sum_gamma_sq_v_sq / static_cast<double>(velocity.count)
            : 0.0;
    const double mean_gamma = std::sqrt(1.0 + mean_gamma_sq_v_sq);

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
    // string's own static long-range tail. Also saved to network_scalars
    // .dat below (2026-09-18, with the user: "we may as well save the
    // processed tension, even if it can be reconstructed") alongside every
    // raw screened/unscreened component, so both this and any other
    // combination remain reconstructable afterwards without having picked
    // one formula in advance.
    const double sum_core_rho_tot =
        energy.n_total * energy.rho_tot_unscreened -
        energy.n_unmasked * energy.rho_tot_screened;
    const double sum_unmasked_tail =
        energy.n_unmasked * (energy.rho_axion_gradient_screened -
                             energy.rho_axion_kinetic_screened);
    const double dx3 = dx * dx * dx;

    // Comoving/physical normalisation for tension (energy per unit length),
    // reworked 2026-09-18 with the user, deviating from a straight reading
    // of conventions.md sec.10 (documented here per the user's explicit
    // go-ahead to do so when the choice is written down).
    //
    // rho_tot_pointwise (Energy.hpp) is the genuine PHYSICAL energy
    // density (energy per physical volume) -- derived directly from the
    // canonically normalised phi Lagrangian, dt physical, grad physical.
    // Lattice points sit on a uniform comoving grid, so at fixed tau every
    // cell has the same physical volume (R dx)^3 and the arithmetic mean
    // over points IS the physical-volume average -- sum_core_rho_tot above
    // is dimensionless-count x rho, so the physical ENERGY in the core
    // cells is R(tau)^3 * dx^3 * sum_core_rho_tot, not dx^3 * sum_core_rho_tot:
    // the previous formula omitted this R^3.
    //
    // The previous string_length_in_box = 2*Geom().ProbLength(2) hardcoded
    // "exactly 2 straight strings spanning the box in z", correct only for
    // the T1 straight_string_test IC (and only as a *comoving* length -- it
    // also omitted converting to physical length, R*length_comoving). For a
    // general network (many loops of random orientation, as here) neither
    // assumption holds. Replaced with the same plaquette-count relation
    // XiFormula.hpp already uses to turn N_p into a length: each pierced
    // plaquette represents ~(2/3)*dx of string, the 2/3 being a *statistical*
    // correction for random orientation relative to the lattice (established
    // already, conventions.md sec.8) -- explicitly NOT exact for a single
    // deterministically axis-aligned string like T1's, so this new formula
    // is intended for network runs, not as a drop-in replacement for T1's
    // own (exact, deterministic-length) validation. T1's previously
    // recorded tension value (~3.72, docs/STATUS.md) was measured with the
    // old formula and is not expected to still hold with this one -- no
    // automated regression test depends on the exact number (checked:
    // tests/ has no tension test), so this is flagged for a follow-up
    // re-validation rather than chased now.
    //
    // Physical energy / physical length:
    //   mu = [R^3 dx^3 Sum_core(rho)] / [R * ell_comoving]
    //      = R^2 dx^3 Sum_core(rho) / ell_comoving,   ell_comoving = (2/3) N_p dx
    const double R_tau           = s_background.R(tau);
    const double ell_comoving    = (2.0 / 3.0) * counts.n_p_plain * dx;
    const double tension_core_only =
        (ell_comoving > 0.0)
            ? (R_tau * R_tau) * dx3 * sum_core_rho_tot / ell_comoving
            : 0.0;
    const double tension_core_plus_tail =
        (ell_comoving > 0.0)
            ? (R_tau * R_tau) * dx3 * (sum_core_rho_tot + sum_unmasked_tail) /
                  ell_comoving
            : 0.0;

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
                   << "    rho_radial_kin: unscreened = "
                   << energy.rho_radial_kinetic_unscreened
                   << "  screened = " << energy.rho_radial_kinetic_screened
                   << "\n"
                   << "    rho_radial_grad: unscreened = "
                   << energy.rho_radial_gradient_unscreened
                   << "  screened = " << energy.rho_radial_gradient_screened
                   << "\n"
                   << "    rho_radial_mass: unscreened = "
                   << energy.rho_radial_mass_unscreened
                   << "  screened = " << energy.rho_radial_mass_screened
                   << "\n"
                   << "    n_total = " << energy.n_total
                   << "  n_unmasked = " << energy.n_unmasked
                   << "  tension (core only) = " << tension_core_only
                   << "  tension (core+tail) = " << tension_core_plus_tail
                   << "\n"
                   << "    <gamma^2 v^2> = " << mean_gamma_sq_v_sq
                   << "  <gamma> = " << mean_gamma
                   << "  N_corners = " << velocity.count << "\n";

    // Restart-safe output bookkeeping (2026-09-18, fixing a bug flagged
    // earlier): a genuine restart is detected via amr.restart (not our own
    // tau_i-relative clock, since GRAmr::get_restart_time() -- set at
    // post_init/post_restart -- reports 0 even for a fresh run, and our
    // tau starts at tau_i != 0, so "> 0" alone cannot tell the two apart).
    // first_step true only for a brand-new run's very first snapshot,
    // which makes SmallDataIO rename any pre-existing file to
    // ".old.<random>" instead of silently appending to it (the original
    // bug); false on a real restart, which instead opens the existing
    // file for read+append so remove_duplicate_time_data() below can trim
    // any provisional rows from a run segment that is being redone.
    const bool is_restart =
        amrex::ParmParse("amr").countval("restart") > 0;
    const amrex::Real restart_time_tau =
        is_restart ? amrex::Real(s_tau_i + get_gr_amr_ptr()->get_restart_time())
                  : amrex::Real(0.0);
    // Passing tau itself as "dt" (rather than the numerical sub-step) is
    // deliberate: it guarantees SmallDataIO's restart-detection window
    // (m_time < m_restart_time + m_dt + eps) covers the first post-restart
    // snapshot regardless of how large the log(m_r/H)-spaced gap to it is,
    // while being a complete no-op on a fresh run (m_restart_time = 0
    // sentinel there).
    const bool first_network_scalars_step =
        !is_restart && !s_wrote_network_scalars_header;
    s_wrote_network_scalars_header = true;

    SmallDataIO network_scalars_file("network_scalars", tau, tau,
                                     restart_time_tau, SmallDataIO::APPEND,
                                     first_network_scalars_step);
    if (first_network_scalars_step)
    {
        network_scalars_file.write_header_line(
            {"N_p", "N_p_weighted", "xi", "xi_weighted", "m_r_over_H",
             "rho_tot_unscreened", "rho_tot_screened",
             "rho_axion_kin_unscreened", "rho_axion_kin_screened",
             "rho_axion_grad_unscreened", "rho_axion_grad_screened",
             "rho_radial_kin_unscreened", "rho_radial_kin_screened",
             "rho_radial_grad_unscreened", "rho_radial_grad_screened",
             "rho_radial_mass_unscreened", "rho_radial_mass_screened",
             "n_total", "n_unmasked", "mean_gamma_sq_v_sq", "mean_gamma",
             "n_velocity_corners", "tension_core_only",
             "tension_core_plus_tail"});
    }
    network_scalars_file.remove_duplicate_time_data();
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
        static_cast<amrex::Real>(energy.rho_radial_kinetic_unscreened),
        static_cast<amrex::Real>(energy.rho_radial_kinetic_screened),
        static_cast<amrex::Real>(energy.rho_radial_gradient_unscreened),
        static_cast<amrex::Real>(energy.rho_radial_gradient_screened),
        static_cast<amrex::Real>(energy.rho_radial_mass_unscreened),
        static_cast<amrex::Real>(energy.rho_radial_mass_screened),
        static_cast<amrex::Real>(energy.n_total),
        static_cast<amrex::Real>(energy.n_unmasked),
        static_cast<amrex::Real>(mean_gamma_sq_v_sq),
        static_cast<amrex::Real>(mean_gamma),
        static_cast<amrex::Real>(velocity.count),
        static_cast<amrex::Real>(tension_core_only),
        static_cast<amrex::Real>(tension_core_plus_tail)};
    network_scalars_file.write_time_data_line(data_row);

    // Spectrum (conventions.md sec.9/sec.12, milestone-1.md task 1.9):
    // opt-in (computing an FFT every snapshot is not free), gated behind
    // axion_strings.compute_spectrum. Masking is applied at exactly one
    // place -- MaskedFieldBuffer.hpp, filling the buffer handed to the FFT
    // -- per CLAUDE.md constraint 5.
    //
    // Both the screened (s_energy_masking, guaranteed scheme B by
    // AxionStringsParams::check_params() whenever this flag is on) and
    // unscreened (MaskingScheme::None, mirroring the plane_wave_test IC's
    // own check) spectra are computed and saved side by side (2026-09-18,
    // with the user: "output the unscreened spectrum to allow for
    // comparison and judging the impact of screening").
    bool compute_spectrum_flag = false;
    amrex::ParmParse("axion_strings")
        .queryAdd("compute_spectrum", compute_spectrum_flag);
    if (compute_spectrum_flag)
    {
        auto compute_and_report_spectrum =
            [&](const MaskingParams &masking, const char *label)
        {
            amrex::MultiFab a_dot_buffer(state_new.boxArray(),
                                         state_new.DistributionMap(), 1, 0);
            fill_masked_a_dot_buffer(a_dot_buffer, state_new, masking,
                                     s_background.R(tau));
            const Spectrum spectrum = compute_spectrum(a_dot_buffer, Geom());

            // Direct real-space mean square of the exact same buffer that
            // was fed to the FFT -- the honest way to get <(masked
            // a_dot)^2>, rather than trying to back it out of
            // axion_kinetic_energy's own f_a^2/2 prefactors (error-prone;
            // avoided deliberately). None's weight is 1 everywhere, so
            // energy.n_total (not n_unmasked) is the right divisor for
            // both cases.
            const double sum_sq = static_cast<double>(amrex::MultiFab::Dot(
                a_dot_buffer, 0, a_dot_buffer, 0, 1, 0));
            const double real_space_mean_sq = sum_sq / energy.n_total;

            const double n_total_sq = energy.n_total * energy.n_total;
            // Parseval (conventions.md sec.9, verified with the user): a
            // spatial *average* of a squared real field equals
            // Sum_p|X_p|^2 / N_total^2 (one factor of N_total from the
            // average itself, one from Parseval for this "no
            // normalisation, round-trip = N_total" FFT convention) --
            // exact, independent of how the shell-binning below groups
            // modes.
            const double parseval_full =
                spectrum.full_cube_energy / n_total_sq;
            const double parseval_inscribed =
                spectrum.inscribed_sphere_energy / n_total_sq;

            double shell_integral = 0.0;
            for (double s : spectrum.shell_average)
            {
                shell_integral += s; // integral over |k|, bin width 1
            }
            shell_integral /= n_total_sq;

            amrex::Print()
                << "  [AxionStrings spectrum, " << label
                << "] <a_dot^2> real-space = " << real_space_mean_sq
                << "  Parseval (full cube) = " << parseval_full
                << "  Parseval (inscribed sphere) = " << parseval_inscribed
                << "\n"
                << "    full_cube/inscribed_sphere ratio = "
                << (spectrum.full_cube_energy /
                   spectrum.inscribed_sphere_energy)
                << "  shell-binned integral (sum_s S(s)/n_total^2, "
                  "approximate) = "
                << shell_integral << "\n";

            return std::make_tuple(spectrum, real_space_mean_sq,
                                   parseval_full, parseval_inscribed);
        };

        MaskingParams none_masking{};
        none_masking.scheme = MaskingScheme::None;
        const auto [spectrum_screened, real_space_mean_sq_screened,
                   parseval_full_screened, parseval_inscribed_screened] =
            compute_and_report_spectrum(s_energy_masking, "screened");
        const auto [spectrum_unscreened, real_space_mean_sq_unscreened,
                   parseval_full_unscreened, parseval_inscribed_unscreened] =
            compute_and_report_spectrum(none_masking, "unscreened");

        // Persisted to axion_spectrum.dat (2026-09-18, with the user): the
        // shape S(p) is what a spectral index q actually gets measured
        // from, so it must be saved, not just cross-checked in the log.
        //
        // tau is repeated as the first column on every row (not just once
        // per block) rather than using SmallDataIO's blank-line block
        // separators, so that a single flat scan (network_scalars.dat's
        // own remove_duplicate_time_data(), which assumes the first
        // column is always time) is still exactly correct for trimming
        // provisional rows after a restart -- a blank line there would
        // break its std::stod parse.
        const bool first_spectrum_step =
            !is_restart && !s_wrote_spectrum_header;
        s_wrote_spectrum_header = true;

        SmallDataIO axion_spectrum_file("axion_spectrum", tau, tau,
                                        restart_time_tau, SmallDataIO::APPEND,
                                        first_spectrum_step);
        if (first_spectrum_step)
        {
            const std::vector<std::string> spectrum_header{
                "shell_average_screened", "shell_average_unscreened",
                "full_cube_energy_screened", "full_cube_energy_unscreened",
                "inscribed_sphere_energy_screened",
                "inscribed_sphere_energy_unscreened",
                "real_space_mean_sq_screened",
                "real_space_mean_sq_unscreened", "parseval_full_screened",
                "parseval_full_unscreened", "parseval_inscribed_screened",
                "parseval_inscribed_unscreened"};
            const std::vector<std::string> spectrum_pre_header{"tau",
                                                               "mode_index"};
            axion_spectrum_file.write_header_line(spectrum_header,
                                                  spectrum_pre_header);
        }
        axion_spectrum_file.remove_duplicate_time_data();
        // Both spectra share the same domain/binning, so their
        // shell_average vectors are always the same length.
        for (std::size_t mode_index = 0;
            mode_index < spectrum_screened.shell_average.size();
            ++mode_index)
        {
            const std::vector<amrex::Real> coords{
                tau, static_cast<amrex::Real>(mode_index)};
            const std::vector<amrex::Real> row{
                static_cast<amrex::Real>(
                    spectrum_screened.shell_average[mode_index]),
                static_cast<amrex::Real>(
                    spectrum_unscreened.shell_average[mode_index]),
                static_cast<amrex::Real>(spectrum_screened.full_cube_energy),
                static_cast<amrex::Real>(
                    spectrum_unscreened.full_cube_energy),
                static_cast<amrex::Real>(
                    spectrum_screened.inscribed_sphere_energy),
                static_cast<amrex::Real>(
                    spectrum_unscreened.inscribed_sphere_energy),
                static_cast<amrex::Real>(real_space_mean_sq_screened),
                static_cast<amrex::Real>(real_space_mean_sq_unscreened),
                static_cast<amrex::Real>(parseval_full_screened),
                static_cast<amrex::Real>(parseval_full_unscreened),
                static_cast<amrex::Real>(parseval_inscribed_screened),
                static_cast<amrex::Real>(parseval_inscribed_unscreened)};
            axion_spectrum_file.write_data_line(row, coords);
        }
    }

    // Optional 2D visualisation snapshot (output follow-up, 2026-09-18,
    // with the user): off by default, since it is only for spot-checking
    // the run by eye, not a routine diagnostic. See ProjectionKernel.hpp
    // for why it's a line-of-sight max of the *unscreened* rho_tot (fixed
    // to the z-axis) plus an overlaid xy-plaquette string-hit count.
    bool save_projection_flag = false;
    amrex::ParmParse("axion_strings")
        .queryAdd("save_projection", save_projection_flag);
    if (save_projection_flag)
    {
        const Projection projection = compute_energy_projection(
            state_new, dx, s_background.R(tau), lambda, s_background.b_inv,
            tau, Geom().Domain());

        const bool first_projection_step =
            !is_restart && !s_wrote_projection_header;
        s_wrote_projection_header = true;

        SmallDataIO axion_projection_file(
            "axion_projection", tau, tau, restart_time_tau,
            SmallDataIO::APPEND, first_projection_step);
        if (first_projection_step)
        {
            const std::vector<std::string> projection_header{
                "rho_tot_max_unscreened", "string_hit_count"};
            const std::vector<std::string> projection_pre_header{
                "tau", "i", "j"};
            axion_projection_file.write_header_line(projection_header,
                                                     projection_pre_header);
        }
        axion_projection_file.remove_duplicate_time_data();
        for (int i = 0; i < projection.nx; ++i)
        {
            for (int j = 0; j < projection.ny; ++j)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(projection.ny) +
                    static_cast<std::size_t>(j);
                const std::vector<amrex::Real> coords{
                    tau, static_cast<amrex::Real>(i),
                    static_cast<amrex::Real>(j)};
                const std::vector<amrex::Real> row{
                    static_cast<amrex::Real>(
                        projection.rho_tot_max_unscreened[idx]),
                    static_cast<amrex::Real>(projection.string_hit_count[idx])};
                axion_projection_file.write_data_line(row, coords);
            }
        }
    }
}

void AxionStringsLevel::tag_cells(amrex::TagBoxArray &tags,
                                 amrex::Real a_regrid_threshold)
{
    // Real (string-based) tagger, milestone-2 Phase 1 (StringTagger.hpp):
    // replaces the geometric FixedGridsTagger placeholder used throughout
    // milestone 1.
    BL_PROFILE("AxionStringsLevel::tag_cells()");

    // No refinement during relaxation (2026-09-19, with the user): a_time
    // there is pre-evolution's own local clock (tau_pre), not tau = s_tau_i
    // + a_time -- computing log(m_r/H) from it via s_background would be
    // meaningless (wrong background entirely). Relaxation is an artificial
    // process anyway, not real cosmological evolution, so "no refinement
    // needed yet" is also the physically sensible answer here, matching
    // the schedule-gating below's own spirit for the very start of a run.
    if (s_phase == Phase::Relaxing)
    {
        return;
    }

    // Schedule-gating: only allow tagging to create the next level once
    // log(m_r/H) has actually reached its threshold (BoxPlan.hpp::
    // compute_amr_box_plan) -- otherwise the current level's own
    // resolution is, by construction, still adequate, and tagging anyway
    // would refine long before it's needed. The current level's own index
    // (0-based) is exactly the index into the schedule for the *next*
    // level's threshold (s_level_add_log_mr_over_h[ell-1] is level ell's
    // threshold). Absent for max_level=0 and Moore mode (guarded together
    // at parameter-read time), in which case tagging is never gated here.
    const int current_level = Level();
    if (s_has_level_schedule)
    {
        if (current_level >=
            static_cast<int>(s_level_add_log_mr_over_h.size()))
        {
            return; // already at the finest permitted level
        }
        const amrex::Real a_time = get_state_data(state_index).curTime();
        const double tau         = s_tau_i + static_cast<double>(a_time);
        const double log_mr_over_h =
            -std::log(static_cast<double>(s_background.H_over_mr_direct(tau)));
        if (log_mr_over_h <
            s_level_add_log_mr_over_h[static_cast<std::size_t>(current_level)])
        {
            return; // this level's own resolution is still adequate
        }
    }

    // FillBoundary: the plaquette test reads (i+1,j+1,k+1) neighbours, and
    // the optional gradient criterion's Laplacian needs
    // FourthOrderDerivatives' usual 2-cell-wide stencil -- both covered by
    // the state's default ghost count (>=3), same precondition as every
    // other per-cell diagnostic pass in this file.
    amrex::MultiFab &state_new = get_new_data(state_index);
    state_new.FillBoundary(Geom().periodicity());

    const auto &tag_arrs   = tags.arrays();
    const auto &state_arrs = state_new.const_arrays();
    const amrex::Real dx   = Geom().CellSize(0);

    StringTagger tagger{dx, s_tagging_params};

    amrex::ParallelFor(
        tags, [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { tagger(ix, iy, iz, state_arrs[box_no], tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}

int AxionStringsLevel::okToContinue()
{
    if (s_phase == Phase::Relaxing)
    {
        // Never stops on its own here: the transition to Evolving (once
        // the xi-monitoring loop's target is reached) happens in place
        // inside specific_post_timestep(), not by ending the run.
        return 1;
    }
    if (!s_has_tau_f)
    {
        // Moore mode: apply_box_plan does not derive tau_f there (see its
        // own comment) -- fall back to evolution.stop_time/max_steps, set
        // by hand (docs/STATUS.md's known Moore-mode limitation).
        return 1;
    }
    const amrex::Real tau_now = s_tau_i + get_state_data(state_index).curTime();
    return (tau_now < s_tau_f) ? 1 : 0;
}

void AxionStringsLevel::apply_pre_evolution_to_main_rescale()
{
    // Pre-evolution -> main handoff (conventions.md sec.7), applied in
    // place within the same run (2026-09-18, with the user: replacing the
    // old restart-based handoff entirely -- each pre-evolution run is only
    // ever used for one main run, so no ensemble-reuse value is lost, and
    // this avoids a second cluster job submission and an intermediate
    // full-grid checkpoint on a memory-constrained cluster). Math is
    // unchanged from the old specific_post_restart().

    // Pre-evolution's own clock (tau_pre) starts at 0, so the a_time this
    // level has reached *is* tau_pre at the moment of transition.
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

    // a_time keeps counting up through the transition (it never resets to
    // 0), so s_tau_i must absorb tau_pre_end for tau = s_tau_i + a_time to
    // still equal tau_i right now, and to correctly track the main
    // evolution's tau from here on.
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
