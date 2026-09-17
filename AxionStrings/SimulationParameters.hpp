#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "Wave.hpp"

// Milestone-1 task 1.1 placeholder: only the Wave model is wired up, to
// establish the build/run baseline. Superseded by the complex-scalar
// parameter block of conventions.md sec.5 in task 1.4.
class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();
        Wave::params_t::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
