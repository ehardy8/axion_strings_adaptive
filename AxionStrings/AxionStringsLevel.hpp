#ifndef AXIONSTRINGSLEVEL_HPP_
#define AXIONSTRINGSLEVEL_HPP_

#include "AxionStringsRHS.hpp"
#include "DefaultLevelBld.hpp"
#include "DerivedVariables.hpp"
#include "GRAmrLevel.hpp"

class AxionStringsLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    //! Define data descriptors.
    static void variableSetUp();

    //! Initialize data at problem start-up.
    void initData() override;

    //! Advance this level for one step
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           const amrex::Real a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    void specific_update_ode(amrex::MultiFab &a_soln) override {};

    // to do post each time step on every level
    void specific_post_timestep() override {};

    //! Error estimation for regridding.
    void tag_cells(amrex::TagBoxArray &tags,
                   amrex::Real a_regrid_threshold) override;

    template <class model_t>
    void eval_model_specific_rhs(amrex::MultiFab &a_soln,
                                 amrex::MultiFab &a_rhs);

  private:

    AxionStringsLevel &getLevel(int lev)
    {
        return dynamic_cast<AxionStringsLevel &>(parent->getLevel(lev));
    }
};

#endif /* AXIONSTRINGSLEVEL_HPP_ */
