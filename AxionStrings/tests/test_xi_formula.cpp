#include "doctest.h"

#include "../XiFormula.hpp"

// milestone-1.md task 1.5/1.6, conventions.md sec.8: xi <-> N_p relation.

TEST_CASE("plaquette_count_for_target_xi inverts xi_from_plaquette_count")
{
    const double N_p_in = 12345.0;
    const double dx      = 0.05;
    const double N        = 64;
    const double L_tilde   = dx * N;
    const double a_inv     = 2.0;
    const double tau       = 3.7;

    const double xi =
        xi_from_plaquette_count(N_p_in, dx, L_tilde, a_inv, tau);
    const double N_p_out =
        plaquette_count_for_target_xi(xi, static_cast<int>(N), L_tilde,
                                      a_inv, tau);

    CHECK(N_p_out == doctest::Approx(N_p_in).epsilon(1.0e-9));
}

TEST_CASE("target N_p scales as expected with tau and xi_target")
{
    const int N          = 128;
    const double L_tilde  = 2.0;
    const double a_inv    = 2.0;

    const double N_p_at_tau1 =
        plaquette_count_for_target_xi(1.0, N, L_tilde, a_inv, 1.0);
    const double N_p_at_tau2 =
        plaquette_count_for_target_xi(1.0, N, L_tilde, a_inv, 2.0);

    // N_p_target ~ 1/tau^2 at fixed xi_target (matches xi ~ N_p * tau^2)
    CHECK(N_p_at_tau1 / N_p_at_tau2 == doctest::Approx(4.0).epsilon(1.0e-9));

    // linear in xi_target at fixed tau
    const double N_p_double_xi =
        plaquette_count_for_target_xi(2.0, N, L_tilde, a_inv, 1.0);
    CHECK(N_p_double_xi == doctest::Approx(2.0 * N_p_at_tau1).epsilon(1.0e-9));
}
