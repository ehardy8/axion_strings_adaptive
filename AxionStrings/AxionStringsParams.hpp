#ifndef AXIONSTRINGSPARAMS_HPP_
#define AXIONSTRINGSPARAMS_HPP_

// Reads the axion_strings.* parameter block (conventions.md sec.4-5) into a
// Background. Kept separate from Background.hpp so that header stays free
// of any AMReX/GRTeclyn dependency (see its own comment).
//
// CLAUDE.md constraint 2: these are only ever called once at start-up
// (SimulationParameters::check_params and AxionStringsLevel::variableSetUp),
// never from specific_eval_rhs.

#include "Background.hpp"
#include "GRParmParse.hpp"

namespace AxionStringsParams
{

inline Background read_background()
{
    GRParmParse pp("axion_strings");

    double a_inv = 2.0; // radiation domination (conventions.md sec.2)
    pp.get("a_inv", a_inv);
    if (a_inv <= 1.0)
    {
        pp.error("a_inv", "must be > 1");
    }

    CTauSchedule sched{};
    pp.get("c0", sched.c0);

    const bool has_c1         = pp.contains("c1");
    const bool has_tau_switch = pp.contains("tau_switch");
    if (has_c1 != has_tau_switch)
    {
        pp.error("c1", "c1 and tau_switch must both be set, or both left "
                       "unset (conventions.md sec.4: at most one switch)");
    }
    sched.has_switch = has_c1;
    if (sched.has_switch)
    {
        pp.get("c1", sched.c1);
        pp.get("tau_switch", sched.tau_switch);
        if (sched.tau_switch <= 0.0)
        {
            pp.error("tau_switch", "must be > 0");
        }
    }

    return {a_inv, sched};
}

// Conformal time at which the AMReX clock (a_time = 0) begins. tau0 = 1 is
// the fixed reference point where m_r = H (conventions.md sec.5); tau_i is
// an independent run parameter picking where the simulation actually starts.
inline double read_tau_i()
{
    GRParmParse pp("axion_strings");
    double tau_i = 1.0;
    pp.get("tau_i", tau_i);
    if (tau_i <= 0.0)
    {
        pp.error("tau_i", "must be > 0 (conformal time must stay positive)");
    }
    return tau_i;
}

inline void check_params()
{
    read_background();
    read_tau_i();
}

} // namespace AxionStringsParams

#endif // AXIONSTRINGSPARAMS_HPP_
