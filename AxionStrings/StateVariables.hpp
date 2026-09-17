#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"

// Complex scalar psi = psi1 + i psi2, with Pi_i = psi_i' (conventions.md
// sec.3). psi = R(tau) phi / v is the rescaled field actually evolved.
enum
{
    c_psi1,
    c_psi2,
    c_Pi1,
    c_Pi2,

    NUM_VARS,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names{"psi1", "psi2", "Pi1", "Pi2"};

static const std::array<BCParity, NUM_VARS> parities = {
    BCParity::even, BCParity::even, BCParity::even, BCParity::even};
static const std::array<amrex::Real, NUM_VARS> asymptotic_values{};
// Periodic boundary conditions are used throughout, so parities aren't
// actually used, but they must be defined.

} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
