#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "AxionStringsParams.hpp"

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        // AxionStringsParams::check_params() runs first: it may inject
        // geometry.prob_extent/amr.n_cell/evolution.stop_time from the
        // axion_strings.* box-planning block (task 1.4), and
        // BaseParameterChecker::check_params() below requires those to
        // already be set (it reads geometry.prob_extent with a bare .get()).
        AxionStringsParams::check_params();
        BaseParameterChecker::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
