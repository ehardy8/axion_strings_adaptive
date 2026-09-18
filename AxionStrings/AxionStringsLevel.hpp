#ifndef AXIONSTRINGSLEVEL_HPP_
#define AXIONSTRINGSLEVEL_HPP_

#include "AxionStringsParams.hpp"
#include "Background.hpp"
#include "DefaultLevelBld.hpp"
#include "GRAmrLevel.hpp"
#include "Masking.hpp"
#include "PreEvolutionBackground.hpp"

class AxionStringsLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    //! Define data descriptors, and cache the axion_strings.* parameters.
    static void variableSetUp();

    //! Initialize data at problem start-up.
    void initData() override;

    //! Advance this level for one step
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           const amrex::Real a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    void specific_update_ode(amrex::MultiFab &a_soln) override {};

    // to do post each time step on every level
    void specific_post_timestep() override;

    //! Error estimation for regridding.
    void tag_cells(amrex::TagBoxArray &tags,
                   amrex::Real a_regrid_threshold) override;

    // The authoritative stop condition (Main_AxionStrings.cpp's loop, via
    // amr.okToContinue()): stops once tau reaches the box-planning-derived
    // tau_f, correctly spanning an optional relaxation phase (Phase::
    // Relaxing) and the main evolution (Phase::Evolving) without caring how
    // long relaxation took, since the transition between them is in place,
    // not a restart. Falls back to leaving evolution.stop_time/max_steps as
    // the stop condition (never returns 0) when box planning could not
    // derive tau_f (Moore mode -- see AxionStringsParams::apply_box_plan).
    int okToContinue() override;

    // Which phase of a single continuous run this level is currently in
    // (2026-09-18, with the user: replacing the old restart-based pre-
    // evolution -> main handoff -- each pre-evolution run was only ever
    // used for one main run anyway, so nothing is lost, and this avoids a
    // second cluster job submission and an intermediate full-grid
    // checkpoint). Every run starts in Evolving unless axion_strings.
    // ic_mode == "fourier_relaxed", which starts in Relaxing; specific_post
    // _timestep() transitions Relaxing -> Evolving in place (rescale, via
    // apply_pre_evolution_to_main_rescale()) once the xi-monitoring loop's
    // target is reached, never restarting the process.
    enum class Phase
    {
        Relaxing,
        Evolving
    };
    inline static Phase s_phase{Phase::Evolving};

    // Background and c(tau) schedule, cached once in variableSetUp() from
    // the axion_strings.* parameters (conventions.md sec.4-5). Conformal
    // time is tau = s_tau_i + a_time, since the AMReX clock always starts
    // at a_time = 0 and never resets (including across the in-place
    // Relaxing -> Evolving transition, which instead adjusts s_tau_i -- see
    // apply_pre_evolution_to_main_rescale()).
    inline static Background s_background{};
    inline static amrex::Real s_tau_i{1.0};

    // Box-planning-derived final tau (AxionStringsParams::apply_box_plan's
    // axion_strings.derived_tau_f, read back here), and whether it was
    // actually available (absent in Moore mode -- see okToContinue()).
    inline static double s_tau_f{0.0};
    inline static bool s_has_tau_f{false};

    // Relaxation phase only (conventions.md sec.7).
    inline static PreEvolutionBackground s_pre_background{};
    inline static double s_xi_target{0.0};
    inline static AxionStringsParams::XiCheckCadence s_xi_cadence{};
    inline static long s_steps_since_xi_check{0};
    inline static long s_xi_check_interval{50};

    // Whether each output file's header has been written yet (fresh-run
    // bookkeeping only -- a genuine restart never re-writes a header, see
    // specific_post_timestep).
    inline static bool s_wrote_network_scalars_header{false};
    inline static bool s_wrote_spectrum_header{false};
    inline static bool s_wrote_projection_header{false};

    // Screening used for the energy diagnostics (task 1.8) -- a runtime
    // parameter, never a compile-time constant (CLAUDE.md constraint 4).
    inline static MaskingParams s_energy_masking{};

    // Output cadence (2026-09-18, with the user): the full per-snapshot
    // diagnostics only run once log(m_r/H) has advanced past
    // s_next_output_log_mr_over_h, which starts at
    // axion_strings.output_first_log_mr_over_h and is bumped by
    // axion_strings.output_delta_log_mr_over_h each time it is crossed.
    inline static AxionStringsParams::OutputCadence s_output_cadence{};
    inline static double s_next_output_log_mr_over_h{0.0};

  private:

    AxionStringsLevel &getLevel(int lev)
    {
        return dynamic_cast<AxionStringsLevel &>(parent->getLevel(lev));
    }

    // The old specific_post_restart's handoff rescale (conventions.md
    // sec.7), now applied in place from specific_post_timestep() the
    // moment the xi-monitoring loop's target is reached, rather than via
    // an external restart. Derivation unchanged: kappa = R_main(tau_i) /
    // R_pre(tau_pre_end) from the psi = R phi/v chain rule.
    void apply_pre_evolution_to_main_rescale();
};

#endif /* AXIONSTRINGSLEVEL_HPP_ */
