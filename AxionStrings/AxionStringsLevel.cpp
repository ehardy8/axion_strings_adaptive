#include "AxionStringsLevel.hpp"
#include "AxionStringsRHS.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateTypes.hpp"
#include <numeric>

void AxionStringsLevel::variableSetUp()
{
    BL_PROFILE("AxionStringsLevel::variableSetUp()");

    // Set up the state variables
    state_variable_set_up();

    // The first two derived variables calculate the analytic solution
    // for phi and Pi (only meaningful while the Wave placeholder model is
    // in place -- see task 1.1 note in Wave.hpp).

    const std::string &comp_type                = {"analytic_soln"};
    const amrex::Vector<std::string> comp_names = {"phi_analytic",
                                                   "Pi_analytic"};

    int ncomp_analytic{static_cast<int>(comp_names.size())};

    derive_lst.add(
        comp_type, amrex::IndexType::TheCellType(), ncomp_analytic, comp_names,
        calc_analytic_solution, [=](const amrex::Box &box) { return box; },
        &amrex::cell_quartic_interp);
    derive_lst.addComponent("analytic_soln", desc_lst, state_index, 0, 1);

    const int ncomp_rho{1};
    const int num_ghosts_rho{2};

    derive_lst.add(
        "rho", amrex::IndexType::TheCellType(), ncomp_rho,
        calc_energy_density<Wave>, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts_rho); },
        &amrex::cell_quartic_interp);

    derive_lst.addComponent("rho", desc_lst, state_index, 0, NUM_VARS);
}

void AxionStringsLevel::initData()
{
    BL_PROFILE("AxionStringsLevel::initData()");

    std::array<amrex::Real, AMREX_SPACEDIM> center{};

    amrex::ParmParse pp;
    pp.get("geometry.center", center);

    amrex::MultiFab &state_new = get_new_data(state_index);

    int dcomp{0};
    const amrex::Real current_time{0.0};

    calc_analytic_mf_3d<Wave>(state_new, dcomp, geom, current_time);
}

void AxionStringsLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                         amrex::MultiFab &a_rhs,
                                         const amrex::Real a_time)
{
    BL_PROFILE("AxionStringsLevel::specific_eval_rhs()");

    eval_model_specific_rhs<Wave>(a_soln, a_rhs);

    amrex::Gpu::streamSynchronize();
}

template <class model_t>
void AxionStringsLevel::eval_model_specific_rhs(amrex::MultiFab &a_soln,
                                               amrex::MultiFab &a_rhs)
{
    const auto dx                 = Geom().CellSize(0);
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    model_t my_model;
    AxionStringsRHS rhs(dx, my_model);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { rhs(ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]); });
}

void AxionStringsLevel::tag_cells(amrex::TagBoxArray &tags,
                                 amrex::Real a_regrid_threshold)
{
    BL_PROFILE("AxionStringsLevel::tag_cells()");

    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrs   = tags.arrays();
    const auto &state_arrs = state_new.arrays();

    const amrex::Real dx    = Geom().CellSize(0);
    const int current_level = Level();

    FixedGridsTagger my_tagging_criterion{dx, current_level};

    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { my_tagging_criterion(ix, iy, iz, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}
