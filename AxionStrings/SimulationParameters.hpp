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
        BaseParameterChecker::check_params();
        AxionStringsParams::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
