#ifndef AXIONSTRINGSLEVEL_HPP_
#define AXIONSTRINGSLEVEL_HPP_

#include "Background.hpp"
#include "DefaultLevelBld.hpp"
#include "GRAmrLevel.hpp"

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

    // Background and c(tau) schedule, cached once in variableSetUp() from
    // the axion_strings.* parameters (conventions.md sec.4-5). Conformal
    // time is tau = s_tau_i + a_time, since the AMReX clock always starts
    // at a_time = 0.
    inline static Background s_background{};
    inline static amrex::Real s_tau_i{1.0};

    // Whether the network_scalars.dat header has been written yet.
    inline static bool s_wrote_network_scalars_header{false};

  private:

    AxionStringsLevel &getLevel(int lev)
    {
        return dynamic_cast<AxionStringsLevel &>(parent->getLevel(lev));
    }
};

#endif /* AXIONSTRINGSLEVEL_HPP_ */
