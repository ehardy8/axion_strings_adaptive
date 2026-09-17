#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"

// assign an enum to each variable
enum
{
    c_phi,
    c_Pi,

    NUM_VARS,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names{"phi", "Pi"};

static const std::array<BCParity, NUM_VARS> parities = {BCParity::even,
                                                         BCParity::even};
static const std::array<amrex::Real, NUM_VARS> asymptotic_values{};
// Periodic boundary conditions are used throughout, so parities aren't
// actually used, but they must be defined.

} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
