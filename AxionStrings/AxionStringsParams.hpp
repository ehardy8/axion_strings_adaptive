#ifndef AXIONSTRINGSPARAMS_HPP_
#define AXIONSTRINGSPARAMS_HPP_

// Reads the axion_strings.* parameter block (conventions.md sec.4-5) into a
// Background. Kept separate from Background.hpp so that header stays free
// of any AMReX/GRTeclyn dependency (see its own comment).
//
// CLAUDE.md constraint 2: these are only ever called once at start-up
// (SimulationParameters::check_params and AxionStringsLevel::variableSetUp),
// never from specific_eval_rhs.

#include "Background.hpp"
#include "BoxPlan.hpp"
#include "GRParmParse.hpp"
#include "Masking.hpp"
#include "PreEvolutionBackground.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace AxionStringsParams
{

enum class Mode
{
    Main,
    PreEvolution
};

// axion_strings.mode: "main" (default) or "pre_evolution" (conventions.md
// sec.7). Pre-evolution and the main run are separate executions
// communicating through a checkpoint file, not two branches of a single
// run -- this selects which one this execution is.
inline Mode read_mode()
{
    GRParmParse pp("axion_strings");
    std::string mode = "main";
    pp.queryAdd("mode", mode);
    if (mode == "main")
    {
        return Mode::Main;
    }
    if (mode == "pre_evolution")
    {
        return Mode::PreEvolution;
    }
    pp.error("mode", "must be \"main\" or \"pre_evolution\"");
    return Mode::Main; // unreachable; pp.error aborts
}

inline Background read_background()
{
    GRParmParse pp("axion_strings");

    double a_inv = 2.0; // radiation domination (conventions.md sec.2)
    pp.get("a_inv", a_inv);
    if (a_inv <= 1.0)
    {
        pp.error("a_inv", "must be > 1");
    }

    CTauSchedule sched{};
    pp.get("c0", sched.c0);

    const bool has_c1         = pp.contains("c1");
    const bool has_tau_switch = pp.contains("tau_switch");
    if (has_c1 != has_tau_switch)
    {
        pp.error("c1", "c1 and tau_switch must both be set, or both left "
                       "unset (conventions.md sec.4: at most one switch)");
    }
    sched.has_switch = has_c1;
    if (sched.has_switch)
    {
        pp.get("c1", sched.c1);
        pp.get("tau_switch", sched.tau_switch);
        if (sched.tau_switch <= 0.0)
        {
            pp.error("tau_switch", "must be > 0");
        }
    }

    return {a_inv, sched};
}

// Conformal time at which the AMReX clock (a_time = 0) begins. tau0 = 1 is
// the fixed reference point where m_r = H (conventions.md sec.5); tau_i is
// an independent run parameter picking where the simulation actually starts.
inline double read_tau_i()
{
    GRParmParse pp("axion_strings");
    double tau_i = 1.0;
    pp.get("tau_i", tau_i);
    if (tau_i <= 0.0)
    {
        pp.error("tau_i", "must be > 0 (conformal time must stay positive)");
    }
    return tau_i;
}

// Pre-evolution's own background (conventions.md sec.7, milestone-1.md
// task 1.5) -- see PreEvolutionBackground.hpp for the derivation. Only
// meaningful when read_mode() == Mode::PreEvolution.
inline PreEvolutionBackground read_pre_evolution_background()
{
    GRParmParse pp("axion_strings.pre_evolution");
    double gamma_pre = 1.0;
    pp.get("gamma", gamma_pre);
    if (gamma_pre <= 0.0)
    {
        pp.error("gamma", "must be > 0");
    }
    return PreEvolutionBackground(gamma_pre);
}

// Target xi at the *main run's* start (tau_i, not necessarily tau0 -- an
// independent input in general): pre-evolution stops once the measured xi
// (on its own grid) first drops to this value (sec.8's xi<->N_p relation
// inverted at tau=tau_i, via XiFormula.hpp). Confirmed with the user: this
// is a direct run input, not derived from an external (xi0,H_c,H0) table.
inline double read_xi_target()
{
    GRParmParse pp("axion_strings.pre_evolution");
    double xi_target{};
    pp.get("xi_target", xi_target);
    if (xi_target <= 0.0)
    {
        pp.error("xi_target", "must be > 0");
    }
    return xi_target;
}

struct FourierICParams
{
    double k_max_over_mr{};
    double mean_square_variance{};
    std::uint64_t seed{};
};

// The two Fourier-mode IC generation parameters of conventions.md sec.7,
// plus the RNG seed (CLAUDE.md constraint 7: recorded in the checkpoint and
// output metadata, and reproducibility depends on it being fixed).
inline FourierICParams read_fourier_ic_params(const std::string &prefix)
{
    GRParmParse pp(prefix);
    FourierICParams params{};
    pp.get("k_max_over_mr", params.k_max_over_mr);
    pp.get("mean_square_variance", params.mean_square_variance);
    if (params.k_max_over_mr <= 0.0)
    {
        pp.error("k_max_over_mr", "must be > 0");
    }
    if (params.mean_square_variance < 0.0)
    {
        pp.error("mean_square_variance", "must be >= 0");
    }
    long long seed = 0;
    pp.get("seed", seed);
    params.seed = static_cast<std::uint64_t>(seed);
    return params;
}

// Set (only on the specific restart command that performs the pre-
// evolution -> main handoff) to trigger AxionStringsLevel::
// specific_post_restart's rescale. A later restart of the main run's own
// progress (e.g. resuming after a crash) must NOT re-apply it, so this is
// not inferred from amr.restart alone -- it needs to be requested
// explicitly, each time.
inline bool read_restart_from_pre_evolution()
{
    GRParmParse pp("axion_strings");
    bool restart_from_pre_evolution = false;
    pp.queryAdd("restart_from_pre_evolution", restart_from_pre_evolution);
    return restart_from_pre_evolution;
}

struct XiCheckCadence
{
    long coarse_interval{50};
    long medium_interval{10};
    long fine_interval{2};
    double medium_threshold{1.5};
    double fine_threshold{1.1};
};

// The adaptive xi-check cadence requested by the user: check infrequently
// while xi is far from the target (measuring it is slow), more often once
// within medium_threshold x of it, and more often still within
// fine_threshold x.
inline XiCheckCadence read_xi_check_cadence()
{
    GRParmParse pp("axion_strings.pre_evolution");
    XiCheckCadence cadence{};
    pp.queryAdd("check_interval_coarse", cadence.coarse_interval);
    pp.queryAdd("check_interval_medium", cadence.medium_interval);
    pp.queryAdd("check_interval_fine", cadence.fine_interval);
    pp.queryAdd("check_threshold_medium", cadence.medium_threshold);
    pp.queryAdd("check_threshold_fine", cadence.fine_threshold);
    if (cadence.coarse_interval < 1 || cadence.medium_interval < 1 ||
        cadence.fine_interval < 1)
    {
        pp.error("check_interval_coarse", "all check intervals must be >= 1");
    }
    if (!(cadence.fine_threshold > 1.0 &&
         cadence.fine_threshold < cadence.medium_threshold))
    {
        pp.error("check_threshold_medium",
                "must have 1 < check_threshold_fine < check_threshold_medium");
    }
    return cadence;
}

struct OutputCadence
{
    double first_log_mr_over_h{};
    double delta_log_mr_over_h{};
};

// Main-run diagnostic output cadence (2026-09-18, with the user): the full
// per-snapshot diagnostics (plaquette count, energy/velocity reductions,
// the spectrum FFT) are too expensive to repeat every coarse step, so they
// only run at snapshots spaced by delta_log_mr_over_h in log(m_r/H),
// starting from the first snapshot at first_log_mr_over_h -- both required
// (no default value is fixed in conventions.md, and a silent default here
// would itself be a physics-adjacent numerical choice).
inline OutputCadence read_output_cadence()
{
    GRParmParse pp("axion_strings");
    OutputCadence cadence{};
    pp.get("output_first_log_mr_over_h", cadence.first_log_mr_over_h);
    pp.get("output_delta_log_mr_over_h", cadence.delta_log_mr_over_h);
    if (cadence.delta_log_mr_over_h <= 0.0)
    {
        pp.error("output_delta_log_mr_over_h", "must be > 0");
    }
    return cadence;
}

// Masking (conventions.md sec.10, milestone-1.md task 1.7). The threshold
// is a runtime parameter, never a compile-time constant (CLAUDE.md
// constraint 4) -- it will be scanned; conventions.md sec.10/sec.14 record
// 0.8, 0.9 and 0.95 all having been tried without reconciliation, so no
// value here should be read as settled.
inline MaskingParams read_masking_params(const std::string &prefix)
{
    GRParmParse pp(prefix);
    MaskingParams params{};

    std::string scheme = "A";
    pp.queryAdd("scheme", scheme);
    if (scheme == "none" || scheme == "None")
    {
        params.scheme = MaskingScheme::None;
    }
    else if (scheme == "A")
    {
        params.scheme = MaskingScheme::A;
    }
    else if (scheme == "B")
    {
        params.scheme = MaskingScheme::B;
    }
    else
    {
        pp.error("scheme", "must be \"A\", \"B\" or \"none\"");
    }

    pp.queryAdd("threshold", params.threshold);
    if (params.scheme == MaskingScheme::B &&
        (params.threshold <= 0.0 || params.threshold >= 1.0))
    {
        pp.error("threshold", "must be in (0, 1) for scheme B");
    }

    return params;
}

// Box planning (conventions.md sec.5, milestone-1.md task 1.4). Skipped
// entirely if axion_strings.N is absent, so ad hoc/manual grid setups (e.g.
// the tasks 1.2/1.3 smoke tests) are unaffected. When present, N, N1, N2
// derive L_tilde/tau_f (and hence geometry.prob_extent/evolution.stop_time)
// for a constant-c run, or -- whenever a switch is configured, since the
// only supported use of a switch is the fat->Moore protocol (sec.4) -- run
// the Moore-phase dynamic-range check instead (sec.5, "Moore-phase box
// planning"): the general formula is singular at c = a_inv and does not
// apply there.
//
// Already-set geometry.prob_extent/amr.n_cell/evolution.stop_time are
// cross-checked rather than overwritten, so a run can also be configured
// entirely by hand; a mismatch aborts before the run starts.
inline void apply_box_plan(const Background &background, double tau_i,
                          Mode mode)
{
    GRParmParse pp("axion_strings");
    if (!pp.contains("N"))
    {
        return;
    }

    int N{};
    double N1{};
    double N2{};
    pp.get("N", N);
    pp.get("N1", N1);
    pp.get("N2", N2);
    if (N <= 0 || N1 <= 0.0 || N2 <= 0.0)
    {
        pp.error("N", "N, N1 and N2 must all be > 0");
    }

    const double a_inv = background.a_inv;
    constexpr double tol = 1.0e-9;
    const bool moore_mode = background.c_sched.has_switch ||
                            std::abs(background.c_sched.c0 - a_inv) < tol;

    amrex::Print() << "AxionStrings box planning (conventions.md sec.5):\n"
                   << "  a_inv = " << a_inv
                   << ", b_inv = " << background.b_inv
                   << ", R0 = " << background.R0
                   << ", tau0 = " << Background::tau0 << "\n"
                   << "  N = " << N << ", N1 = " << N1 << ", N2 = " << N2
                   << "\n";

    GRParmParse amr_pp("amr");
    int max_level = 0;
    amr_pp.queryAdd("max_level", max_level);
    if (max_level == 0)
    {
        if (amr_pp.contains("n_cell"))
        {
            std::array<int, AMREX_SPACEDIM> n_cell{};
            amr_pp.get("n_cell", n_cell);
            for (int n : n_cell)
            {
                if (n != N)
                {
                    amr_pp.error(
                        "n_cell",
                        "does not match axion_strings.N -- either remove "
                        "amr.n_cell to let it be derived, or fix "
                        "axion_strings.N to match");
                }
            }
        }
        else
        {
            amr_pp.addarr("n_cell", std::vector<int>{N, N, N});
            amrex::Print() << "  -> amr.n_cell set to " << N << " " << N
                           << " " << N << "\n";
        }
    }
    else
    {
        amrex::Print()
            << "  amr.max_level = " << max_level
            << " > 0: axion_strings.N is the effective finest resolution "
               "(conventions.md sec.5/sec.11), not amr.n_cell -- set "
               "amr.n_cell by hand.\n";
    }

    if (moore_mode)
    {
        double gamma{};
        double target_log_range{};
        pp.get("gamma", gamma);
        pp.get("target_log_range", target_log_range);
        if (gamma <= 0.0)
        {
            pp.error("gamma", "must be > 0");
        }

        const double max_log_range =
            moore_max_log_dynamic_range(N, N1, N2, gamma);

        amrex::Print()
            << "  Moore mode (c = a_inv, reached "
            << (background.c_sched.has_switch ? "via a switch" : "from tau_i")
            << "): gamma = " << gamma
            << ", requested log range = " << target_log_range
            << ", achievable log range = " << max_log_range << "\n"
            << "  geometry.prob_extent and evolution.stop_time are not "
               "derived in Moore mode -- set them by hand (sec.5 Moore-phase "
               "box planning is a dynamic-range check here, not a closed "
               "form for L_tilde/tau_f).\n";

        if (max_log_range < target_log_range)
        {
            pp.error(
                "target_log_range",
                "cannot be reached with the requested N, N1, N2, gamma (see "
                "achievable log range printed above)");
        }
        return;
    }

    const double c   = background.c_sched.c0;
    const auto plan  = compute_general_box_plan(N, N1, N2, a_inv, c);

    amrex::Print()
        << "  tau_f = " << plan.tau_f << " (from tau_i = " << tau_i << ")\n"
        << "  L_tilde (main) = " << plan.L_tilde
        << ", delta_x (main) = " << plan.delta_x << "\n"
        << "  delta_tau (leapfrog bound, dx/3) = "
        << plan.delta_tau_leapfrog_bound
        << " -- provisional, not yet re-derived for the RK integrator "
           "(sec.5/sec.6)\n";

    if (plan.tau_f <= tau_i)
    {
        pp.error("N", "the requested N, N1, N2, c configuration gives tau_f "
                      "<= tau_i; cannot evolve forward from tau_i");
    }

    // Decision (see docs/STATUS.md): pre-evolution uses the *same*
    // geometry.prob_extent as the main run, not the sec.7 L_tilde_init
    // formula's (generally different) value. Reasons: (1) AMReX's restart
    // mechanism does not resize the domain -- verified directly -- so a
    // differently-sized pre-evolution box could not actually be handed off
    // via amr.restart without a custom regridding step nobody has designed
    // yet; (2) the sec.7 formula only optimises pre-evolution's own
    // resolution/box-size trade-off (it is computed without reference to
    // the actual, a priori unknown, stopping time), while
    // specific_post_restart's rescale (kappa, derived from the psi=R phi/v
    // chain rule) already corrects the physical normalisation exactly,
    // whatever box size was used -- so using L_tilde_main throughout costs
    // nothing beyond a possibly-suboptimal (not incorrect) choice of
    // pre-evolution resolution. L_tilde_init is still computed and printed
    // here for reference/comparison, deliberately not acted on.
    const double L_tilde = plan.L_tilde;
    if (mode == Mode::PreEvolution)
    {
        const double L_tilde_init_fyi =
            pre_evolution_L_tilde(plan.L_tilde, a_inv, c, tau_i);
        amrex::Print()
            << "  L_tilde (pre-evolution, sec.7 formula, FYI only) = "
            << L_tilde_init_fyi
            << " -- NOT used; this run uses L_tilde_main = " << L_tilde
            << " instead (see docs/STATUS.md)\n";
    }

    GRParmParse geom_pp("geometry");
    if (geom_pp.contains("prob_extent"))
    {
        std::array<double, AMREX_SPACEDIM> prob_extent{};
        geom_pp.get("prob_extent", prob_extent);
        for (double L : prob_extent)
        {
            if (std::abs(L - L_tilde) > 1.0e-6 * L_tilde)
            {
                geom_pp.error(
                    "prob_extent",
                    "does not match the box-planning-derived L_tilde -- "
                    "either remove geometry.prob_extent to let it be "
                    "derived, or fix axion_strings.N/N1/N2 to match");
            }
        }
    }
    else
    {
        geom_pp.addarr("prob_extent",
                       std::vector<double>{L_tilde, L_tilde, L_tilde});
        amrex::Print() << "  -> geometry.prob_extent set to " << L_tilde
                       << " " << L_tilde << " " << L_tilde << "\n";
    }

    if (mode == Mode::PreEvolution)
    {
        // Pre-evolution does not stop at a fixed evolution.stop_time -- it
        // stops when the measured xi first drops to the target (sec.7/
        // task 1.5's xi-monitoring loop). Leave evolution.stop_time (and
        // max_steps) as a user-set safety cap.
        return;
    }

    const double stop_time = plan.tau_f - tau_i;
    GRParmParse evolution_pp("evolution");
    if (evolution_pp.contains("stop_time"))
    {
        double existing_stop_time{};
        evolution_pp.get("stop_time", existing_stop_time);
        if (std::abs(existing_stop_time - stop_time) > 1.0e-6 * stop_time)
        {
            evolution_pp.error(
                "stop_time",
                "does not match tau_f - tau_i derived from axion_strings.N/"
                "N1/N2 -- either remove evolution.stop_time to let it be "
                "derived, or fix the box-planning inputs to match");
        }
    }
    else
    {
        evolution_pp.add("stop_time", stop_time);
        amrex::Print() << "  -> evolution.stop_time set to " << stop_time
                       << "\n";
    }
}

inline void check_params()
{
    const Mode mode              = read_mode();
    const Background background = read_background();
    const double tau_i          = read_tau_i();
    apply_box_plan(background, tau_i, mode);

    if (mode == Mode::PreEvolution)
    {
        read_pre_evolution_background();
        read_xi_target();
        read_fourier_ic_params("axion_strings.pre_evolution");
        read_xi_check_cadence();

        GRParmParse amr_pp("amr");
        int check_int = -1;
        amr_pp.queryAdd("check_int", check_int);
        if (check_int < 0)
        {
            amr_pp.error(
                "check_int",
                "must be >= 0 for a pre_evolution run -- its whole purpose "
                "is to hand off to the main run through a checkpoint "
                "(conventions.md sec.7)");
        }
    }
    else
    {
        std::string ic_mode = "homogeneous";
        GRParmParse("axion_strings").queryAdd("ic_mode", ic_mode);
        if (ic_mode == "fourier")
        {
            // Main run, generating Fourier-mode ICs directly at tau_i --
            // the user's "skip pre-evolution" option (conventions.md
            // sec.7).
            read_fourier_ic_params("axion_strings");
        }
        // ic_mode == "homogeneous"/"straight_string_test": ad hoc/manual
        // configs (tasks 1.2/1.3/1.6 smoke tests), or a restart (either
        // the pre-evolution handoff via restart_from_pre_evolution, or an
        // ordinary resumption of the main run's own progress) -- nothing
        // further to validate here.

        read_output_cadence();

        // The persisted spectrum (axion_spectrum.dat) is only useful if it
        // is genuinely screened -- scheme "A" has no top-hat at all, and
        // "none" is the raw unscreened field -- so require scheme B
        // whenever the routine per-snapshot spectrum output is on. (The
        // plane_wave_test IC's own one-shot spectrum check in initData()
        // is a separate, self-contained identity check that hardcodes
        // MaskingScheme::None regardless of this setting, so it is
        // unaffected.)
        bool compute_spectrum_flag = false;
        GRParmParse spectrum_pp("axion_strings");
        spectrum_pp.queryAdd("compute_spectrum", compute_spectrum_flag);
        if (compute_spectrum_flag &&
            read_masking_params("axion_strings.masking").scheme !=
                MaskingScheme::B)
        {
            spectrum_pp.error(
                "compute_spectrum",
                "requires axion_strings.masking.scheme = B -- the "
                "spectrum saved to axion_spectrum.dat must be the "
                "genuinely screened field");
        }
    }
}

} // namespace AxionStringsParams

#endif // AXIONSTRINGSPARAMS_HPP_
