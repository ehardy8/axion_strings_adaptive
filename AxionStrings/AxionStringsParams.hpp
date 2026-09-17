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
#include "BoxPlan.hpp"
#include "GRParmParse.hpp"

#include <array>
#include <cmath>
#include <vector>

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

// Box planning (conventions.md sec.5, milestone-1.md task 1.4). Skipped
// entirely if axion_strings.N is absent, so ad hoc/manual grid setups (e.g.
// the tasks 1.2/1.3 smoke tests) are unaffected. When present, N, N1, N2
// derive L_tilde/tau_f (and hence geometry.prob_extent/evolution.stop_time)
// for a constant-c run, or -- whenever a switch is configured, since the
// only supported use of a switch is the fat->Moore protocol (sec.4) -- run
// the Moore-phase dynamic-range check instead (sec.5, "Moore-phase box
// planning"): the general formula is singular at c = a_inv and does not
// apply there.
//
// Already-set geometry.prob_extent/amr.n_cell/evolution.stop_time are
// cross-checked rather than overwritten, so a run can also be configured
// entirely by hand; a mismatch aborts before the run starts.
inline void apply_box_plan(const Background &background, double tau_i)
{
    GRParmParse pp("axion_strings");
    if (!pp.contains("N"))
    {
        return;
    }

    int N{};
    double N1{};
    double N2{};
    pp.get("N", N);
    pp.get("N1", N1);
    pp.get("N2", N2);
    if (N <= 0 || N1 <= 0.0 || N2 <= 0.0)
    {
        pp.error("N", "N, N1 and N2 must all be > 0");
    }

    const double a_inv = background.a_inv;
    constexpr double tol = 1.0e-9;
    const bool moore_mode = background.c_sched.has_switch ||
                            std::abs(background.c_sched.c0 - a_inv) < tol;

    amrex::Print() << "AxionStrings box planning (conventions.md sec.5):\n"
                   << "  a_inv = " << a_inv
                   << ", b_inv = " << background.b_inv
                   << ", R0 = " << background.R0
                   << ", tau0 = " << Background::tau0 << "\n"
                   << "  N = " << N << ", N1 = " << N1 << ", N2 = " << N2
                   << "\n";

    GRParmParse amr_pp("amr");
    int max_level = 0;
    amr_pp.queryAdd("max_level", max_level);
    if (max_level == 0)
    {
        if (amr_pp.contains("n_cell"))
        {
            std::array<int, AMREX_SPACEDIM> n_cell{};
            amr_pp.get("n_cell", n_cell);
            for (int n : n_cell)
            {
                if (n != N)
                {
                    amr_pp.error(
                        "n_cell",
                        "does not match axion_strings.N -- either remove "
                        "amr.n_cell to let it be derived, or fix "
                        "axion_strings.N to match");
                }
            }
        }
        else
        {
            amr_pp.addarr("n_cell", std::vector<int>{N, N, N});
            amrex::Print() << "  -> amr.n_cell set to " << N << " " << N
                           << " " << N << "\n";
        }
    }
    else
    {
        amrex::Print()
            << "  amr.max_level = " << max_level
            << " > 0: axion_strings.N is the effective finest resolution "
               "(conventions.md sec.5/sec.11), not amr.n_cell -- set "
               "amr.n_cell by hand.\n";
    }

    if (moore_mode)
    {
        double gamma{};
        double target_log_range{};
        pp.get("gamma", gamma);
        pp.get("target_log_range", target_log_range);
        if (gamma <= 0.0)
        {
            pp.error("gamma", "must be > 0");
        }

        const double max_log_range =
            moore_max_log_dynamic_range(N, N1, N2, gamma);

        amrex::Print()
            << "  Moore mode (c = a_inv, reached "
            << (background.c_sched.has_switch ? "via a switch" : "from tau_i")
            << "): gamma = " << gamma
            << ", requested log range = " << target_log_range
            << ", achievable log range = " << max_log_range << "\n"
            << "  geometry.prob_extent and evolution.stop_time are not "
               "derived in Moore mode -- set them by hand (sec.5 Moore-phase "
               "box planning is a dynamic-range check here, not a closed "
               "form for L_tilde/tau_f).\n";

        if (max_log_range < target_log_range)
        {
            pp.error(
                "target_log_range",
                "cannot be reached with the requested N, N1, N2, gamma (see "
                "achievable log range printed above)");
        }
        return;
    }

    const double c   = background.c_sched.c0;
    const auto plan  = compute_general_box_plan(N, N1, N2, a_inv, c);

    amrex::Print()
        << "  tau_f = " << plan.tau_f << " (from tau_i = " << tau_i << ")\n"
        << "  L_tilde = " << plan.L_tilde << ", delta_x = " << plan.delta_x
        << "\n"
        << "  delta_tau (leapfrog bound, dx/3) = "
        << plan.delta_tau_leapfrog_bound
        << " -- provisional, not yet re-derived for the RK integrator "
           "(sec.5/sec.6)\n";

    if (plan.tau_f <= tau_i)
    {
        pp.error("N", "the requested N, N1, N2, c configuration gives tau_f "
                      "<= tau_i; cannot evolve forward from tau_i");
    }

    GRParmParse geom_pp("geometry");
    if (geom_pp.contains("prob_extent"))
    {
        std::array<double, AMREX_SPACEDIM> prob_extent{};
        geom_pp.get("prob_extent", prob_extent);
        for (double L : prob_extent)
        {
            if (std::abs(L - plan.L_tilde) > 1.0e-6 * plan.L_tilde)
            {
                geom_pp.error(
                    "prob_extent",
                    "does not match the box-planning-derived L_tilde -- "
                    "either remove geometry.prob_extent to let it be "
                    "derived, or fix axion_strings.N/N1/N2 to match");
            }
        }
    }
    else
    {
        geom_pp.addarr(
            "prob_extent",
            std::vector<double>{plan.L_tilde, plan.L_tilde, plan.L_tilde});
        amrex::Print() << "  -> geometry.prob_extent set to " << plan.L_tilde
                       << " " << plan.L_tilde << " " << plan.L_tilde << "\n";
    }

    const double stop_time = plan.tau_f - tau_i;
    GRParmParse evolution_pp("evolution");
    if (evolution_pp.contains("stop_time"))
    {
        double existing_stop_time{};
        evolution_pp.get("stop_time", existing_stop_time);
        if (std::abs(existing_stop_time - stop_time) > 1.0e-6 * stop_time)
        {
            evolution_pp.error(
                "stop_time",
                "does not match tau_f - tau_i derived from axion_strings.N/"
                "N1/N2 -- either remove evolution.stop_time to let it be "
                "derived, or fix the box-planning inputs to match");
        }
    }
    else
    {
        evolution_pp.add("stop_time", stop_time);
        amrex::Print() << "  -> evolution.stop_time set to " << stop_time
                       << "\n";
    }
}

inline void check_params()
{
    const Background background = read_background();
    const double tau_i          = read_tau_i();
    apply_box_plan(background, tau_i);
}

} // namespace AxionStringsParams

#endif // AXIONSTRINGSPARAMS_HPP_
