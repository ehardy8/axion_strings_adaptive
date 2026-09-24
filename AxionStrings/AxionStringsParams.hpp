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
#include "CurvatureKernel.hpp"
#include "GRParmParse.hpp"
#include "Masking.hpp"
#include "PreEvolutionBackground.hpp"
#include "StringTagger.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
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

// Handoff smoothing (2026-09-23, with the user): apply_pre_evolution_to_
// main_rescale()'s R-rescale matches psi/Pi exactly across the pre-
// evolution -> main handoff, but lambda and curvature_term_coeff
// themselves still jump discontinuously between the pre-evolution and
// main schedules (they generically sit on different (a_inv,c) curves) --
// confirmed by direct simulation to excite a measurable core-breathing
// transient (oscillatory rho_radial_kin/mass in network_scalars.dat, a
// bump in axion_spectrum.dat that migrates from the old core scale
// towards the new one over several tau) that the usual masking scheme
// does not remove. On by default (3 periods) as of 2026-09-23: tested at
// n=3 and n=6 with no further improvement from the longer window (n=3
// already captures essentially all of the benefit smoothing alone can
// buy -- docs/STATUS.md has the full comparison), so 3 is kept as the
// standing default rather than left opt-in; every fourier_relaxed config
// that doesn't set this explicitly now gets the smoothed handoff. Set to
// 0 to recover the old instantaneous-jump behaviour. When > 0,
// AxionStringsLevel blends lambda/curvature_term_coeff from their frozen
// pre-evolution values at the handoff to the main schedule's own values
// via a smoothstep, over this many main-schedule core-oscillation periods
// (2*pi/sqrt(lambda_main(tau_i)) -- a physical, resolution-independent
// timescale, not a bare tau window).
inline double read_handoff_transition_n_periods()
{
    GRParmParse pp("axion_strings.pre_evolution");
    double n_periods = 3.0;
    pp.queryAdd("handoff_transition_n_periods", n_periods);
    if (n_periods < 0.0)
    {
        pp.error("handoff_transition_n_periods", "must be >= 0");
    }
    return n_periods;
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

// Real (string-based) tagging criterion, milestone-2 Phase 1
// (StringTagger.hpp/StringTaggerParams, conventions.md sec.11). The
// primary (plaquette-based) criterion has no parameters of its own; the
// secondary gradient criterion (Buschmann et al.: dx_ell^2|laplacian(
// psi_i)| > threshold, their value 0.04) is optional -- defaults to
// effectively off (a threshold no field configuration can cross), since
// conventions.md sec.11 says to assess this criterion rather than copy
// it. queryAdd, not get: this is meant to be left unset for most runs.
//
// Radial-gradient criterion (2026-09-19, with the user: "using the
// gradient of the radial mode instead of (or as well as) the winding" for
// the physical-mode network schedule, keeping the existing log(m_r/H)
// level-*addition* schedule -- when a new level is *allowed to exist* --
// completely unchanged; only which per-cell criterion is used to tag
// cells within an already-permitted level changes here). Off by default
// (same "genuinely off" convention as gradient_threshold above), axion_
// strings.tagging.radial_gradient_threshold to enable explicitly.
//
// A physically-derived default was tried and rejected the same day, kept
// here as a documented, ruled-out data point rather than silently
// forgotten: near a global string's core, |psi| = R(tau) |phi|/v rises as
// R(tau) m_r(tau) c1 rho / sqrt(2) (conventions.md sec.8/13's exact
// equilibrium near-core slope c1 = 0.41222, and the extra 1/sqrt(2) from
// phi's own v/sqrt(2) normalisation there), so its comoving spatial
// gradient peaks at exactly R(tau) m_r(tau) c1/sqrt(2) at an isolated,
// static core. Since BoxPlan.hpp's own N2(tau) = 1/(R(tau) dx m_r(tau)),
// dx * |grad(r)|_peak = c1/(sqrt(2) N2(tau)) identically at that core --
// no separate normalisation by R(tau) needed, it cancels exactly. Setting
// the threshold to this value at N2 = the run's own target was meant to
// reproduce "the tagged condition for stationary strings" while
// automatically widening for a degraded-resolution or boosted core.
// **Tested directly against params_amr_validation_128.txt (the same
// config Phase 1 validated at ~37% tagged at level 1's first onset with
// the plaquette criterion alone) and rejected**: with this derived
// threshold (0.2915 at N2=1), level 1 came online at exactly the
// predicted TIME=5.091 (schedule-gating itself unaffected, as intended),
// but tagged 100% of the domain immediately, not ~37% -- the ambient
// radial-mode gradient in a real, dense Fourier-relaxed tangle (as
// opposed to a single isolated static string, which is all the
// derivation above modelled) is evidently comparable to or above this
// threshold almost everywhere, not just at cores. The derivation isn't
// necessarily wrong on its own terms, but a single-isolated-core
// calibration doesn't transfer to the network's actual field statistics
// without further work -- if this criterion is revisited for network use,
// treat this as the starting point to explain, not re-derive from scratch.
inline StringTaggerParams read_tagging_params()
{
    GRParmParse pp("axion_strings.tagging");
    StringTaggerParams params{};
    pp.queryAdd("gradient_threshold", params.gradient_threshold);
    if (params.gradient_threshold <= 0.0)
    {
        pp.error("gradient_threshold", "must be > 0");
    }
    pp.queryAdd("radial_gradient_threshold", params.radial_gradient_threshold);
    if (params.radial_gradient_threshold <= 0.0)
    {
        pp.error("radial_gradient_threshold", "must be > 0");
    }
    pp.queryAdd("force_full_refinement", params.force_full_refinement);
    return params;
}

// Regrid-frequency/buffer safety margin (2026-09-19, with the user).
// AMReX only re-tags cells every amr.regrid_int level-native steps, so
// the amr.n_error_buf cell buffer grown around every tagged cell must be
// wide enough that the fastest possible string segment (v=1, these
// units) cannot cross out of the refined-plus-buffer region before the
// next regrid check -- otherwise it can briefly sit unrefined right at a
// coarse-fine boundary, exactly the systematic conventions.md sec.11's
// "level-timing" section warns about. GRAmrLevel::ComputeDt sets
// dt_level = evolution.dt_multiplier * dx_level exactly (checked
// directly in GRTeclyn's source, not assumed), so a v=1 signal crosses
// dt_multiplier cells every level-native step, hence dt_multiplier *
// amr.regrid_int cells over one whole regrid interval -- independent of
// level, since dt_level and dx_level scale together under subcycling.
//
// Rather than tuning amr.n_error_buf and amr.regrid_int independently by
// hand against each other (as every params file did before this), two
// axion_strings.tagging.* inputs express the same choice directly and
// keep them consistent automatically. Deliberately *not* an empirical
// optimum search over their trade-off right now (2026-09-19, with the
// user: "we will ultimately run on a much bigger grid on a cluster so
// the optimal might change dramatically" -- today's small-grid numbers
// would not transfer anyway); this just makes the safety-preserving
// relationship between them explicit and scannable later.
//
// axion_strings.tagging.regrid_interval_steps (int, no default -- this
// whole mechanism is off unless set, so every existing hand-tuned params
// file that already sets amr.n_error_buf/amr.regrid_int directly is
// completely unaffected): how the fixed v=1 safety requirement is
// *distributed* between regrid time and buffer width -- sets amr.
// regrid_int directly, and raising it trades more-frequent regridding
// for a correspondingly larger required buffer (derived below); 1 keeps
// the buffer at the historical bare minimum (regrid every step).
//
// axion_strings.tagging.buffer_safety_factor (double, default 1.0, must
// be >= 1.0): multiplies the bare-minimum derived buffer for extra
// margin beyond the v=1 bound above -- e.g. covering AMReX's own box-
// clustering rounding the buffered region to whole grids/blocking_factor
// multiples in a way that could shave cells off one side, or simply
// wanting headroom before trusting a new configuration. Only meaningful
// once regrid_interval_steps has opted into this mechanism.
//
// amr.n_error_buf, if not already set by hand, is derived as
// ceil(buffer_safety_factor * regrid_interval_steps * dt_multiplier) and
// injected; if already set, it is only checked to be at *least* this
// (more buffer is always safe, just costs more volume -- unlike amr.
// regrid_int below, this is not required to match exactly). amr.
// regrid_int, if not already set, is set to regrid_interval_steps; if
// already set by hand too, it is cross-checked to match exactly (same
// "cross-checked, not overwritten" pattern as every other derived AMR
// parameter in this file, e.g. amr.n_cell/geometry.prob_extent in
// apply_box_plan above).
inline void apply_regrid_buffer_policy()
{
    GRParmParse amr_pp("amr");
    int max_level = 0;
    amr_pp.queryAdd("max_level", max_level);
    if (max_level <= 0)
    {
        return; // no regridding happens; nothing to derive
    }

    GRParmParse tag_pp("axion_strings.tagging");
    int regrid_interval_steps = 0; // 0 = mechanism not requested
    tag_pp.queryAdd("regrid_interval_steps", regrid_interval_steps);
    if (regrid_interval_steps == 0)
    {
        return; // opt-in only -- leaves hand-set amr.n_error_buf/
                // regrid_int (if any) completely alone
    }
    if (regrid_interval_steps < 0)
    {
        tag_pp.error("regrid_interval_steps", "must be > 0");
    }

    double buffer_safety_factor = 1.0;
    tag_pp.queryAdd("buffer_safety_factor", buffer_safety_factor);
    if (buffer_safety_factor < 1.0)
    {
        tag_pp.error("buffer_safety_factor", "must be >= 1.0");
    }

    GRParmParse evolution_pp("evolution");
    double dt_multiplier = 0.0;
    evolution_pp.get("dt_multiplier", dt_multiplier);

    const double n_error_buf_min =
        buffer_safety_factor * regrid_interval_steps * dt_multiplier;
    const int n_error_buf_derived =
        std::max(1, static_cast<int>(std::ceil(n_error_buf_min)));

    amrex::Print()
        << "Regrid buffer policy (conventions.md sec.11): "
           "regrid_interval_steps = "
        << regrid_interval_steps
        << ", buffer_safety_factor = " << buffer_safety_factor
        << ", evolution.dt_multiplier = " << dt_multiplier << "\n"
        << "  -> minimum buffer to keep a v=1 string inside the refined "
           "region between regrids = "
        << n_error_buf_derived << " cells\n";

    if (amr_pp.contains("n_error_buf"))
    {
        int n_error_buf_set = 0;
        amr_pp.get("n_error_buf", n_error_buf_set);
        if (n_error_buf_set < n_error_buf_derived)
        {
            amr_pp.error(
                "n_error_buf",
                "is smaller than the derived safety minimum for the "
                "requested axion_strings.tagging.regrid_interval_steps/"
                "buffer_safety_factor -- either remove amr.n_error_buf to "
                "let it be derived, or raise it to at least the minimum "
                "printed above");
        }
    }
    else
    {
        amr_pp.add("n_error_buf", n_error_buf_derived);
        amrex::Print() << "  -> amr.n_error_buf set to " << n_error_buf_derived
                       << "\n";
    }

    if (amr_pp.contains("regrid_int"))
    {
        int regrid_int_set = 0;
        amr_pp.get("regrid_int", regrid_int_set);
        if (regrid_int_set != regrid_interval_steps)
        {
            amr_pp.error(
                "regrid_int",
                "does not match axion_strings.tagging.regrid_interval_"
                "steps -- either remove amr.regrid_int to let it be "
                "derived, or fix axion_strings.tagging.regrid_interval_"
                "steps to match");
        }
    }
    else
    {
        amr_pp.add("regrid_int", regrid_interval_steps);
        amrex::Print() << "  -> amr.regrid_int set to " << regrid_interval_steps
                       << "\n";
    }
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

// Local string curvature (2026-09-26, with the user, CurvatureKernel.hpp):
// a distribution of curvatures along the network, binned by the
// dimensionless kappa/m_r (curvature relative to the string's own core
// scale). Off by default (axion_strings.compute_curvature), same "opt-in,
// real extra cost" convention as compute_spectrum. The bin range is a
// runtime parameter, never hardcoded (CLAUDE.md constraint 4's spirit,
// same reasoning as the masking threshold) -- kappa/m_r > 1 corresponds
// to a radius of curvature tighter than the core width itself, not
// physically meaningful for a string of finite width, so the default
// upper edge is set just past that; the lower edge is a practical floor
// below which curvature is indistinguishable from lattice-scale noise on
// an already-straight segment.
inline CurvatureParams read_curvature_params()
{
    GRParmParse pp("axion_strings.curvature");
    CurvatureParams params{};
    pp.queryAdd("n_bins", params.n_bins);
    pp.queryAdd("log10_kappa_over_mr_min", params.log10_kappa_over_mr_min);
    pp.queryAdd("log10_kappa_over_mr_max", params.log10_kappa_over_mr_max);
    if (params.n_bins <= 0)
    {
        pp.error("n_bins", "must be > 0");
    }
    if (params.log10_kappa_over_mr_max <= params.log10_kappa_over_mr_min)
    {
        pp.error("log10_kappa_over_mr_max",
                "must be > log10_kappa_over_mr_min");
    }

    // Two independent ways to build "the position/direction of the
    // string" for the curvature measurement (CurvatureKernel.hpp's own
    // header comment has the full story) -- kept switchable rather than
    // picking one, so they can be cross-checked against each other.
    // Defaults to tangent_vector: the one validated first (a live network
    // test's histogram went from 0/720 bins ever populated to a smooth,
    // stable unimodal distribution across 30 snapshots, 2026-09-27).
    std::string method = "tangent_vector";
    pp.queryAdd("method", method);
    if (method == "tangent_vector")
    {
        params.method = CurvatureMethod::TangentVector;
    }
    else if (method == "interpolated_position")
    {
        params.method = CurvatureMethod::InterpolatedPosition;
    }
    else if (method == "hessian_analytic")
    {
        params.method = CurvatureMethod::HessianAnalytic;
    }
    else
    {
        pp.error("method", "must be \"tangent_vector\", "
                           "\"interpolated_position\" or "
                           "\"hessian_analytic\"");
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

// axion_strings.pre_evolution.gamma default (2026-09-22, with the user):
// StringTagger.hpp's tag_cells() refines nothing during Phase::Relaxing,
// and the schedule-gated tagger only brings a finer level online once
// log(m_r/H) actually crosses *that level's* threshold -- so as long as
// axion_strings.log_mr_over_h_i is below the level-1 threshold (checked
// below, not just assumed), the main run genuinely starts with only the
// coarsest (level-0) grid present, and relaxing at anything other than
// dx_base would hand off a field resolved at the wrong resolution for the
// grid it is about to sit on. gamma = 1/dx_base is therefore not a
// separate choice from the box plan above, just an unstated consequence
// of it -- so derive it the same way as amr.n_cell/geometry.prob_extent
// (derive-if-absent, cross-check-if-present) rather than requiring every
// ic_mode=fourier_relaxed parameter file to compute 1/dx_base by hand,
// which is exactly the bug class docs/STATUS.md's params_full_test_1024
// .txt entry documents (a stale/mismatched gamma, caught only because the
// run's own startup printout was read carefully). Unlike amr.n_cell/
// geometry.prob_extent above, an explicit mismatch here only *warns*, not
// errors (2026-09-22, with the user): relaxing at a resolution other than
// dx_base is a legitimate thing to want deliberately, not only ever a
// mistake -- the warning still catches the accidental/stale case (the
// original motivation), it just no longer refuses to run over it.
inline void apply_pre_evolution_gamma_default(const std::string &ic_mode,
                                              double dx_base)
{
    if (ic_mode != "fourier_relaxed")
    {
        return; // no relaxation phase -- gamma is not used at all
    }

    const double gamma_expected = 1.0 / dx_base;
    GRParmParse pre_pp("axion_strings.pre_evolution");
    if (pre_pp.contains("gamma"))
    {
        double gamma_set{};
        pre_pp.get("gamma", gamma_set);
        if (std::abs(gamma_set - gamma_expected) > 1.0e-6 * gamma_expected)
        {
            // Warning, not error (2026-09-22, with the user): a deliberate
            // choice to relax at a resolution other than dx_base is a
            // legitimate thing to want, not just a mistake -- this still
            // surfaces every mismatch (including the accidental,
            // forgot-to-update-it kind this mechanism exists to catch),
            // it just no longer refuses to run over it.
            pre_pp.warning(
                "gamma",
                "does not match the box-planning-derived 1/dx_base (level "
                "0) -- using your value as set. If this wasn't deliberate "
                "(e.g. a stale value left over from a different N/N1/N2/"
                "amr.max_level), remove axion_strings.pre_evolution.gamma "
                "to let it be derived instead (see the box plan printed "
                "above for the dx_base this run actually uses)");
        }
    }
    else
    {
        pre_pp.add("gamma", gamma_expected);
        amrex::Print()
            << "  -> axion_strings.pre_evolution.gamma set to "
            << gamma_expected
            << " (= 1/dx_base): relaxes at the coarsest/level-0 "
               "resolution, matching what the main run's own grid looks "
               "like the moment it starts\n";
    }
}

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

    // Moore should run single-level (conventions.md sec.11: "AMR pays peak
    // cost for no asymptotic saving" there -- Moore's comoving core width
    // grows, so refinement demand is maximal at the start and only falls;
    // there is no phase where AMR is cheaper than just running at the
    // resolution Moore needs throughout). Checked before either Moore
    // branch below, both of which assume amr.n_cell/geometry.prob_extent
    // are single-level quantities.
    if (moore_mode && max_level > 0)
    {
        pp.error("N", "amr.max_level > 0 is not supported reaching Moore "
                      "(conventions.md sec.11: Moore should run "
                      "single-level) -- set amr.max_level = 0 for a "
                      "fat->Moore run");
    }

    // amr.n_cell derivation is deferred past this point: which grid size it
    // derives from (N directly for a single-level run, or N's role as the
    // *effective finest* resolution for AMR -- conventions.md sec.5/sec.11)
    // depends on max_level, resolved together with the rest of box planning
    // below rather than duplicated here.

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

    const double c = background.c_sched.c0;

    std::string ic_mode = "homogeneous";
    GRParmParse("axion_strings").queryAdd("ic_mode", ic_mode);

    if (max_level > 0)
    {
        // AMR box plan (2026-09-19, milestone 2 "Phase 0"/wiring):
        // BoxPlan.hpp's compute_amr_box_plan. N is the *effective finest*
        // resolution (conventions.md sec.5/sec.11 -- not amr.n_cell, which
        // is the coarser base-level grid derived from it here), N1/N2 the
        // same Hubble-patch/core-resolution targets as ever, now understood
        // to be achieved at tau_f using the finest level.
        const auto plan =
            compute_amr_box_plan(N, N1, N2, a_inv, c, max_level);

        amrex::Print()
            << "  AMR: N (effective, finest) = " << N << ", amr.max_level = "
            << max_level << "\n"
            << "  tau_f = " << plan.tau_f << " (from tau_i = " << tau_i
            << ")\n"
            << "  L_tilde = " << plan.L_tilde
            << ", dx_finest = " << plan.dx_finest
            << ", dx_base (level 0) = " << plan.dx_base << "\n";

        if (plan.N_base <= 0)
        {
            pp.error("N", "N is not evenly divisible by 2^amr.max_level -- "
                          "choose an effective resolution N that is a clean "
                          "multiple of 2^max_level");
        }
        if (plan.tau_f <= tau_i)
        {
            pp.error("N", "the requested N, N1, N2, c configuration gives "
                          "tau_f <= tau_i; cannot evolve forward from tau_i");
        }

        amrex::Print() << "  Level-addition schedule (log(m_r/H) at which "
                          "each level must be active):\n";
        for (int ell = 1; ell <= max_level; ++ell)
        {
            amrex::Print()
                << "    level " << ell << ": log(m_r/H) = "
                << plan.log_add[static_cast<std::size_t>(ell - 1)] << "\n";
        }

        if (ic_mode == "fourier_relaxed" &&
            !GRParmParse("axion_strings.pre_evolution").contains("gamma"))
        {
            // apply_pre_evolution_gamma_default's whole premise (only the
            // coarsest level exists when the main run starts) requires
            // log_mr_over_h_i to sit below level 1's own threshold --
            // checked here rather than just assumed. Only when gamma is
            // left to the default, though (2026-09-24, fixing a real bug
            // the user hit): this premise simply does not apply once the
            // user has set axion_strings.pre_evolution.gamma explicitly --
            // that IS "the resolution you actually want to relax at" the
            // warning below itself suggests as the fix, and an earlier
            // version of this check fired unconditionally regardless, so
            // following that exact suggested fix could not actually avoid
            // it.
            //
            // Warning, not error (2026-09-24, with the user -- relaxing
            // this from the abort it used to be): proceeding with the
            // same default (1/dx_base, "assume no refinement yet") is a
            // legitimate deliberate choice here too, same as an explicit-
            // but-mismatched gamma already only warns about just below --
            // the field ends up relaxed at the coarsest level's
            // resolution even though a finer level is already active at
            // tau_i, which is not necessarily wrong (relaxation itself
            // never refines, regardless -- StringTagger.hpp's tag_cells()
            // is a no-op during Phase::Relaxing -- so the relaxed field
            // was always going to be coarsest-level-resolution; the only
            // question is whether that resolution is now coarser, at
            // tau_i, than what log_mr_over_h_i has already committed the
            // main run to needing).
            const double log_mr_over_h_i =
                -std::log(background.H_over_mr_direct(tau_i));
            if (log_mr_over_h_i >= plan.log_add[0])
            {
                amrex::Print()
                    << "  WARNING: axion_strings.log_mr_over_h_i = "
                    << log_mr_over_h_i
                    << " is already at or past the level-1 threshold "
                       "printed above -- the main run will start with "
                       "level 1 already active, but axion_strings.pre_"
                       "evolution.gamma was left unset, so it defaults to "
                       "1/dx_base, assuming only the coarsest level is "
                       "present at tau_i. Proceeding with that default "
                       "anyway (relaxing at the coarsest level's "
                       "resolution regardless of levels already active at "
                       "tau_i) -- set axion_strings.pre_evolution.gamma "
                       "explicitly to relax at a different resolution "
                       "instead.\n";
            }
        }
        apply_pre_evolution_gamma_default(ic_mode, plan.dx_base);

        if (amr_pp.contains("n_cell"))
        {
            std::array<int, AMREX_SPACEDIM> n_cell{};
            amr_pp.get("n_cell", n_cell);
            for (int n : n_cell)
            {
                if (n != plan.N_base)
                {
                    amr_pp.error(
                        "n_cell",
                        "does not match the box-planning-derived base-level "
                        "grid size (N/2^max_level) -- either remove "
                        "amr.n_cell to let it be derived, or fix "
                        "axion_strings.N/N2/amr.max_level to match");
                }
            }
        }
        else
        {
            amr_pp.addarr(
                "n_cell",
                std::vector<int>{plan.N_base, plan.N_base, plan.N_base});
            amrex::Print() << "  -> amr.n_cell set to " << plan.N_base << " "
                           << plan.N_base << " " << plan.N_base << "\n";
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

        pp.add("derived_tau_f", plan.tau_f);
        // Consumed by the tagger (milestone-2 Phase 1, not yet
        // implemented): which level must be active at the current tau.
        pp.addarr("derived_level_add_log_mr_over_h", plan.log_add);

        GRParmParse evolution_pp("evolution");
        if (!evolution_pp.contains("stop_time"))
        {
            evolution_pp.add("stop_time", -1.0);
            amrex::Print()
                << "  -> evolution.stop_time set to -1 (unlimited): "
                   "AxionStringsLevel::okToContinue() is the authoritative "
                   "stop condition now, comparing the live tau against "
                   "axion_strings.derived_tau_f = "
                << plan.tau_f << "\n";
        }
        return;
    }

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
        amrex::Print() << "  -> amr.n_cell set to " << N << " " << N << " "
                       << N << "\n";
    }

    const auto plan = compute_general_box_plan(N, N1, N2, a_inv, c);

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

    // max_level == 0: no refinement at all, so trivially "only the
    // coarsest level is present" -- no schedule check needed, unlike the
    // AMR branch above.
    apply_pre_evolution_gamma_default(ic_mode, plan.delta_x);

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
    // Flat-space (no cosmological expansion) loop simulations (2026-09-19,
    // with the user): box planning, the log(m_r/H) output cadence and the
    // AMR level schedule all assume an expanding background (Background's
    // own required inputs, e.g. log_mr_over_h_i, are not physically
    // meaningful here) and simply do not apply. geometry.prob_extent/
    // amr.n_cell are set directly (like the ad hoc tasks 1.2/1.3 smoke
    // test configs), and the stop condition is evolution.stop_time/
    // max_steps, set by hand -- there is no tau_f to derive.
    bool flat_space = false;
    GRParmParse("axion_strings").queryAdd("flat_space", flat_space);
    if (flat_space)
    {
        read_masking_params("axion_strings.masking");
        read_tagging_params();
        apply_regrid_buffer_policy();
        return;
    }

    const Background background = read_background();
    const double tau_i          = read_tau_i(background);
    apply_box_plan(background, tau_i);
    apply_regrid_buffer_policy();

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
