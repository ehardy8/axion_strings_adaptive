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

    // Signals the pre-evolution -> main run loop (Main_AxionStrings.cpp) to
    // stop once the xi-monitoring loop's target is reached (task 1.5).
    int okToContinue() override;

    // Pre-evolution -> main handoff (conventions.md sec.7): rescales the
    // restored state from pre-evolution's normalisation into the main
    // run's, when axion_strings.restart_from_pre_evolution is set.
    void specific_post_restart() override;

    // Background and c(tau) schedule, cached once in variableSetUp() from
    // the axion_strings.* parameters (conventions.md sec.4-5). Conformal
    // time is tau = s_tau_i + a_time, since the AMReX clock always starts
    // at a_time = 0 (or, after a restart, at whatever time the checkpoint
    // recorded).
    inline static Background s_background{};
    inline static amrex::Real s_tau_i{1.0};
    inline static AxionStringsParams::Mode s_mode{
        AxionStringsParams::Mode::Main};

    // Pre-evolution only (conventions.md sec.7).
    inline static PreEvolutionBackground s_pre_background{};
    inline static double s_xi_target{0.0};
    inline static AxionStringsParams::XiCheckCadence s_xi_cadence{};
    inline static long s_steps_since_xi_check{0};
    inline static long s_xi_check_interval{50};
    inline static bool s_pre_evolution_target_reached{false};

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
};

#endif /* AXIONSTRINGSLEVEL_HPP_ */
