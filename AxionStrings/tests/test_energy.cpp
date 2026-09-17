#include "doctest.h"

#include "../Background.hpp"
#include "../Energy.hpp"

// A check that holds regardless of how docs/STATUS.md's open rho_tot
// normalisation question resolves: the homogeneous vacuum solution
// psi = R(tau), Pi = R'(tau) (tasks 1.2/1.3: an exact solution of the free
// EOM, since Pi - psi/(b_inv tau) = R' - R/(b_inv tau) = 0 identically, as
// R'/R = 1/(b_inv tau) exactly) must give exactly zero energy density --
// uniform, at the potential minimum, at rest relative to the background.

TEST_CASE("rho_tot_pointwise is exactly zero for the homogeneous vacuum "
          "solution psi=R(tau), Pi=R'(tau)")
{
    CTauSchedule sched{};
    sched.c0 = 0.0;
    Background bkg(2.0, sched);

    for (const double tau : {0.5, 1.0, 3.3, 10.0})
    {
        const double R      = bkg.R(tau);
        const double Rprime = bkg.R_prime(tau);
        const double lambda = bkg.lambda(tau);

        const double rho = rho_tot_pointwise(
            /*psi1=*/R, /*psi2=*/0.0, /*Pi1=*/Rprime, /*Pi2=*/0.0,
            /*grad_psi1_sq=*/0.0, /*grad_psi2_sq=*/0.0, R, lambda,
            bkg.b_inv, tau);

        CHECK(rho == doctest::Approx(0.0).epsilon(1.0e-10));
    }
}

TEST_CASE("axion_kinetic_energy_pointwise is exactly zero when theta' = 0")
{
    CHECK(axion_kinetic_energy_pointwise(0.0, 2.0) == 0.0);
}

TEST_CASE("axion_gradient_energy_pointwise is exactly zero when grad(theta) "
          "= 0")
{
    CHECK(axion_gradient_energy_pointwise(0.0, 2.0) == 0.0);
}

TEST_CASE("axion_gradient_energy_pointwise scales as 1/R^2 at fixed "
          "comoving grad(theta)")
{
    const double grad_theta_sq = 1.3;
    const double e1            = axion_gradient_energy_pointwise(grad_theta_sq, 1.0);
    const double e2            = axion_gradient_energy_pointwise(grad_theta_sq, 2.0);
    CHECK(e1 / e2 == doctest::Approx(4.0).epsilon(1.0e-12));
}
