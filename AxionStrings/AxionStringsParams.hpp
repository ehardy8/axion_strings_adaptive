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

    // axion_strings.log_mr_over_h_switch (2026-09-19, with the user): the
    // fat->Moore protocol's switch time, specified physically rather than
    // as a raw tau_switch -- mirrors log_mr_over_h_i's own rationale
    // exactly. Resolved here, before the raw c1/tau_switch reading below,
    // using only a_inv/sched.c0 (already read above) via a throwaway
    // Background -- tau_from_log_mr_over_h only ever touches those two
    // fields (its own doc comment), so this sidesteps needing the *final*
    // Background (which needs the switch resolved first) to resolve the
    // switch. c1 is set to a_inv automatically: reaching Moore is the only
    // supported use of a switch (see apply_box_plan's own comment), so
    // there is no other value for the user to sensibly choose here.
    const bool has_log_switch = pp.contains("log_mr_over_h_switch");
    const bool has_c1         = pp.contains("c1");
    const bool has_tau_switch = pp.contains("tau_switch");
    if (has_log_switch && (has_c1 || has_tau_switch))
    {
        pp.error("log_mr_over_h_switch",
                 "cannot be combined with c1/tau_switch -- pick one way to "
                 "specify the switch");
    }
    if (has_log_switch)
    {
        constexpr double tol = 1.0e-9;
        if (std::abs(a_inv - sched.c0) < tol)
        {
            pp.error("log_mr_over_h_switch",
                     "cannot derive tau_switch when axion_strings.c0 = "
                     "a_inv already -- there is no fat phase to switch "
                     "from");
        }
        double log_mr_over_h_switch{};
        pp.get("log_mr_over_h_switch", log_mr_over_h_switch);

        const Background pre_switch_background{a_inv, CTauSchedule{sched.c0}};
        sched.has_switch = true;
        sched.c1         = a_inv;
        sched.tau_switch =
            pre_switch_background.tau_from_log_mr_over_h(log_mr_over_h_switch);
    }
    else
    {
        if (has_c1 != has_tau_switch)
        {
            pp.error("c1", "c1 and tau_switch must both be set, or both "
                           "left unset (conventions.md sec.4: at most one "
                           "switch)");
        }
        sched.has_switch = has_c1;
        if (sched.has_switch)
        {
            pp.get("c1", sched.c1);
            pp.get("tau_switch", sched.tau_switch);
        }
    }
    if (sched.has_switch && sched.tau_switch <= 0.0)
    {
        pp.error("tau_switch", "must be > 0");
    }

    return {a_inv, sched};
}

// Conformal time at which the AMReX clock (a_time = 0) begins, derived from
// the (more physical, 2026-09-18 with the user) run parameter
// axion_strings.log_mr_over_h_i: log(m_r/H) at the start. tau0 = 1 is the
// fixed reference point where m_r = H (conventions.md sec.5); tau_i is
// wherever that ratio is instead. See Background::tau_from_log_mr_over_h
// for the inversion and its no-switch-at-the-start caveat.
inline double read_tau_i(const Background &background)
{
    GRParmParse pp("axion_strings");
    double log_mr_over_h_i{};
    pp.get("log_mr_over_h_i", log_mr_over_h_i);

    constexpr double tol = 1.0e-9;
    if (std::abs(background.a_inv - background.c_sched.c0) < tol)
    {
        pp.error(
            "log_mr_over_h_i",
            "cannot derive tau_i when axion_strings.c0 = a_inv (Moore mode "
            "from the very start) -- H_over_mr_closed_form's exponent is "
            "singular there; reach Moore only via a switch (axion_strings."
            "c1/tau_switch), not as the starting c0");
    }

    return background.tau_from_log_mr_over_h(log_mr_over_h_i);
}

// Pre-evolution's own background (conventions.md sec.7, milestone-1.md
// task 1.5) -- see PreEvolutionBackground.hpp for the derivation. Only
// meaningful when axion_strings.ic_mode == "fourier_relaxed".
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
// derive L_tilde/tau_f (and hence geometry.prob_extent) for a constant-c
// run, or -- whenever a switch is configured, since the only supported use
// of a switch is the fat->Moore protocol (sec.4) -- run the Moore-phase
// dynamic-range check instead (sec.5, "Moore-phase box planning"): the
// general formula is singular at c = a_inv and does not apply there.
//
// Already-set geometry.prob_extent/amr.n_cell are cross-checked rather
// than overwritten, so a run can also be configured entirely by hand; a
// mismatch aborts before the run starts.
//
// evolution.stop_time is deliberately NOT derived here (2026-09-18, with
// the user: replacing the old restart-based pre-evolution -> main handoff
// with a single continuous run -- see AxionStringsLevel::okToContinue()).
// tau_f (the derived stop condition) is instead injected as
// axion_strings.derived_tau_f, purely informational/for AxionStringsLevel
// ::variableSetUp() to read back, since okToContinue() is now the
// authoritative stop condition for the general (non-Moore) box plan --
// comparing the live tau against this value directly, correctly spanning
// both an optional relaxation phase and the main evolution without caring
// how long relaxation took. evolution.stop_time is instead forced to -1
// (unlimited) if not already set, purely to defeat GRTeclyn's own
// BaseParameterChecker, which otherwise silently defaults it to 1.0 --
// found the hard way, see docs/STATUS.md. Moore mode is the one case that
// still needs evolution.stop_time/max_steps set by hand (unaffected by
// any of this): its box plan does not derive tau_f, so okToContinue()
// cannot use it there either.
inline void apply_box_plan(const Background &background, double tau_i)
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

    if (moore_mode && !background.c_sched.has_switch)
    {
        // c0 = a_inv from the very start -- no preceding fat phase, so the
        // "dx from N2 via R*m_r = const through the fat phase" derivation
        // below does not apply (there is no fat phase). Unchanged from
        // before: a dynamic-range check only, geometry/stop_time by hand.
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
            << "  Moore mode (c = a_inv from tau_i, no fat phase): gamma = "
            << gamma << ", requested log range = " << target_log_range
            << ", achievable log range = " << max_log_range << "\n"
            << "  geometry.prob_extent and evolution.stop_time are not "
               "derived here -- set them by hand (sec.5 Moore-phase box "
               "planning is a dynamic-range check only when there is no "
               "preceding fat phase to derive dx from).\n";

        if (max_log_range < target_log_range)
        {
            pp.error(
                "target_log_range",
                "cannot be reached with the requested N, N1, N2, gamma (see "
                "achievable log range printed above)");
        }
        return;
    }

    if (moore_mode)
    {
        // Reached via a fat->Moore switch (2026-09-19, with the user):
        // fully derived, mirroring the non-Moore branch below. gamma is
        // *derived* from tau_switch, not a separate manual input -- the
        // previous interface let gamma and tau_switch silently disagree,
        // since nothing checked them against each other.
        const double tau_switch = background.c_sched.tau_switch;
        const double gamma      = 1.0 / background.H_over_mr_direct(tau_switch);
        const double R_switch   = background.R(tau_switch);
        const double m_r_switch = std::sqrt(background.lambda(tau_switch));

        const auto plan = compute_moore_box_plan(
            N, N1, N2, gamma, a_inv, background.b_inv, R_switch, m_r_switch,
            tau_switch);

        amrex::Print()
            << "  Moore mode (c = a_inv, reached via a switch at "
               "tau_switch = "
            << tau_switch << ", gamma = m_r/H there = " << gamma << ")\n"
            << "  L_tilde (from N2 at the switch) = " << plan.L_tilde
            << ", delta_x = " << plan.dx << "\n"
            << "  achievable Moore dynamic range D = log(H_switch/H_end) = "
            << plan.D << " -> tau_end = " << plan.tau_end << "\n";

        double target_log_range{};
        if (pp.queryAdd("target_log_range", target_log_range) &&
            plan.D < target_log_range)
        {
            pp.error(
                "target_log_range",
                "cannot be reached with the requested N, N1, N2 (see "
                "achievable log range printed above) -- increase N or "
                "reduce N1/N2/the switch's log(m_r/H)");
        }

        GRParmParse geom_pp("geometry");
        if (geom_pp.contains("prob_extent"))
        {
            std::array<double, AMREX_SPACEDIM> prob_extent{};
            geom_pp.get("prob_extent", prob_extent);
            for (double L : prob_extent)
            {
                if (std::abs(L - plan.L_tilde) > 1.0e-6 * plan.L_tilde)
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
            geom_pp.addarr(
                "prob_extent",
                std::vector<double>{plan.L_tilde, plan.L_tilde, plan.L_tilde});
            amrex::Print() << "  -> geometry.prob_extent set to "
                           << plan.L_tilde << " " << plan.L_tilde << " "
                           << plan.L_tilde << "\n";
        }

        pp.add("derived_tau_f", plan.tau_end);

        GRParmParse evolution_pp("evolution");
        if (!evolution_pp.contains("stop_time"))
        {
            evolution_pp.add("stop_time", -1.0);
            amrex::Print()
                << "  -> evolution.stop_time set to -1 (unlimited): "
                   "AxionStringsLevel::okToContinue() is the authoritative "
                   "stop condition now, comparing the live tau against "
                   "axion_strings.derived_tau_f = "
                << plan.tau_end << "\n";
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

    // FYI-only cross-reference to sec.7's L_tilde_init formula for an
    // optional relaxation phase, computed here purely for comparison (see
    // docs/STATUS.md): this run always uses L_tilde_main instead, since
    // pre-evolution and the main evolution now share one live grid rather
    // than communicating through a checkpoint, so no separate box size is
    // ever actually needed.
    const double L_tilde = plan.L_tilde;
    const double L_tilde_init_fyi =
        pre_evolution_L_tilde(plan.L_tilde, a_inv, c, tau_i);
    amrex::Print() << "  L_tilde (sec.7 pre-evolution formula, FYI only) = "
                   << L_tilde_init_fyi
                   << " -- NOT used; this run uses L_tilde_main = " << L_tilde
                   << " throughout (see docs/STATUS.md)\n";

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

    pp.add("derived_tau_f", plan.tau_f);

    GRParmParse evolution_pp("evolution");
    if (!evolution_pp.contains("stop_time"))
    {
        evolution_pp.add("stop_time", -1.0);
        amrex::Print()
            << "  -> evolution.stop_time set to -1 (unlimited): "
               "AxionStringsLevel::okToContinue() is the authoritative stop "
               "condition now, comparing the live tau against "
               "axion_strings.derived_tau_f = "
            << plan.tau_f << "\n";
    }
}

// axion_strings.ic_mode = "fourier_relaxed" (2026-09-18, with the user,
// replacing the old two-execution restart-based handoff entirely: each
// pre-evolution run was only ever used for one main run anyway -- no
// ensemble-reuse value lost -- and a single continuous run avoids both a
// second cluster job submission and an intermediate full-grid checkpoint):
// generate a Fourier-mode IC in the pre-evolution time-gauge, relax it via
// the xi-monitoring loop, then AxionStringsLevel::specific_post_timestep
// transitions in place (rescale + switch phase) to the main evolution,
// all within the same process, no restart involved. See
// AxionStringsLevel::Phase.
inline void check_params()
{
    const Background background = read_background();
    const double tau_i          = read_tau_i(background);
    apply_box_plan(background, tau_i);

    std::string ic_mode = "homogeneous";
    GRParmParse("axion_strings").queryAdd("ic_mode", ic_mode);
    if (ic_mode == "fourier_relaxed")
    {
        read_pre_evolution_background();
        read_xi_target();
        read_fourier_ic_params("axion_strings.pre_evolution");
        read_xi_check_cadence();
    }
    else if (ic_mode == "fourier")
    {
        // Generating Fourier-mode ICs directly at tau_i -- the user's
        // "skip relaxation entirely" option (conventions.md sec.7).
        read_fourier_ic_params("axion_strings");
    }
    // ic_mode == "homogeneous"/"straight_string_test"/"plane_wave_test":
    // ad hoc/manual configs (tasks 1.2/1.3/1.6/1.9 smoke tests) -- nothing
    // further to validate for the IC itself.

    read_output_cadence();

    // The persisted spectrum (axion_spectrum.dat) is only useful if it is
    // genuinely screened -- scheme "A" has no top-hat at all, and "none"
    // is the raw unscreened field -- so require scheme B whenever the
    // routine per-snapshot spectrum output is on. (The plane_wave_test
    // IC's own one-shot spectrum check in initData() is a separate, self-
    // contained identity check that hardcodes MaskingScheme::None
    // regardless of this setting, so it is unaffected.)
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

} // namespace AxionStringsParams

#endif // AXIONSTRINGSPARAMS_HPP_
