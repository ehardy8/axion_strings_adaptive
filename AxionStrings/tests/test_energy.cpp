#include "doctest.h"

#include "../Background.hpp"
#include "../Energy.hpp"

#include <array>

// A check that holds regardless of how docs/STATUS.md's open rho_tot
// normalisation question resolves: the homogeneous vacuum solution
// psi = R(tau), Pi = R'(tau) (tasks 1.2/1.3: an exact solution of the free
// EOM, since Pi - (R'/R) psi = R' - (R'/R) R = 0 identically) must give
// exactly zero energy density -- uniform, at the potential minimum, at
// rest relative to the background.

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
            /*R_prime_over_R=*/1.0 / (bkg.b_inv * tau));

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

// Radial energy: the kinetic and gradient functions are defined as an
// *exact* orthogonal (radial/tangential) decomposition of the full
// |phi_dot|^2 and |grad phi|^2 -- not an approximation -- so the direct
// correctness check is that radial + tangential reconstructs the full,
// independently-computed quantity exactly, for arbitrary (not just
// vacuum) psi/Pi. The tangential piece is theta_prime^2 * psi_sq / R^4
// (NOT axion_kinetic_energy_pointwise, which uses the *linearised*,
// vacuum-amplitude f_a = sqrt(2)*v prefactor -- the two only coincide
// when psi_sq = R^2, i.e. at the vacuum, checked separately below).
TEST_CASE("radial_kinetic_energy_pointwise + tangential = exact |phi_dot|^2, "
          "for arbitrary (non-vacuum) psi/Pi")
{
    const double b_inv = 1.7;
    const double tau   = 2.3;
    const double R     = 1.9;

    for (const auto &[psi1, psi2, Pi1, Pi2] :
        {std::array<double, 4>{0.6, -0.3, 1.1, 0.4},
         std::array<double, 4>{2.5, 1.8, -0.7, 2.2},
         std::array<double, 4>{-1.2, 0.9, 0.05, -1.6}})
    {
        const double c1 = Pi1 - psi1 / (b_inv * tau);
        const double c2 = Pi2 - psi2 / (b_inv * tau);
        const double psi_sq = psi1 * psi1 + psi2 * psi2;
        const double exact_full_kinetic = (c1 * c1 + c2 * c2) / (R * R * R * R);

        const double theta_prime = (psi1 * Pi2 - Pi1 * psi2) / psi_sq;
        const double tangential_kinetic =
            theta_prime * theta_prime * psi_sq / (R * R * R * R);

        const double radial_kinetic = radial_kinetic_energy_pointwise(
            psi1, psi2, Pi1, Pi2, R, 1.0 / (b_inv * tau));

        CHECK(radial_kinetic + tangential_kinetic ==
             doctest::Approx(exact_full_kinetic).epsilon(1.0e-10));
    }
}

TEST_CASE("radial_kinetic_energy_pointwise is exactly zero at psi=(0,0) "
          "(guarded, not a 0/0 NaN)")
{
    CHECK(radial_kinetic_energy_pointwise(0.0, 0.0, 1.0, -0.4, 1.5,
                                          1.0 / (1.7 * 2.3)) == 0.0);
}

TEST_CASE("radial_gradient_energy_pointwise + tangential = exact "
          "|grad phi|^2, for arbitrary (non-vacuum) psi and grad(psi)")
{
    const double R = 1.9;

    // (psi1, psi2, grad_psi1_sq, grad_psi2_sq, grad_theta_sq_comoving) --
    // grad_theta_sq chosen independently of psi/grad_psi_sq (it is itself
    // an exact function of psi and grad(psi) in the real kernel, but this
    // function only needs the *identity* to hold given consistent inputs,
    // which the derivation guarantees for any grad_psi1/grad_psi2 pair
    // -- checked here by constructing grad_theta_sq_comoving from an
    // explicit grad_psi1, grad_psi2 rather than picking it freely).
    const std::array<std::array<double, 4>, 3> cases{
        {{0.6, -0.3, 1.1, 0.4}, {2.5, 1.8, -0.7, 2.2}, {-1.2, 0.9, 0.05, -1.6}}};
    for (const auto &c : cases)
    {
        const double psi1 = c[0];
        const double psi2 = c[1];
        const double grad_psi1 = c[2]; // single-direction stand-in
        const double grad_psi2 = c[3];
        const double psi_sq = psi1 * psi1 + psi2 * psi2;
        const double grad_psi1_sq = grad_psi1 * grad_psi1;
        const double grad_psi2_sq = grad_psi2 * grad_psi2;
        const double grad_theta_dir =
            (psi1 * grad_psi2 - psi2 * grad_psi1) / psi_sq;
        const double grad_theta_sq = grad_theta_dir * grad_theta_dir;

        const double exact_full_gradient =
            (grad_psi1_sq + grad_psi2_sq) / (R * R * R * R);
        const double tangential_gradient =
            grad_theta_sq * psi_sq / (R * R * R * R);

        const double radial_gradient = radial_gradient_energy_pointwise(
            grad_psi1_sq, grad_psi2_sq, psi_sq, grad_theta_sq, R);

        CHECK(radial_gradient + tangential_gradient ==
             doctest::Approx(exact_full_gradient).epsilon(1.0e-10));
    }
}

TEST_CASE("radial_gradient_energy_pointwise: all gradient energy is "
          "radial at psi_sq=0 (grad_theta_sq is 0 there too, upstream "
          "guard)")
{
    const double grad_psi1_sq = 0.7;
    const double grad_psi2_sq = 1.3;
    const double R            = 1.5;
    CHECK(radial_gradient_energy_pointwise(grad_psi1_sq, grad_psi2_sq, 0.0,
                                           0.0, R) ==
         doctest::Approx((grad_psi1_sq + grad_psi2_sq) / (R * R * R * R)));
}

TEST_CASE("radial_mass_energy_pointwise is exactly zero at the vacuum "
          "(|psi|=R) and matches the direct potential formula off it")
{
    const double R      = 1.5;
    const double lambda = 0.8;
    CHECK(radial_mass_energy_pointwise(R, 0.0, R, lambda) == 0.0);

    const double psi1 = 2.1;
    const double psi2 = -0.4;
    const double psi_sq_minus_R_sq = psi1 * psi1 + psi2 * psi2 - R * R;
    const double expected = 0.25 * lambda * psi_sq_minus_R_sq *
                            psi_sq_minus_R_sq / (R * R * R * R);
    CHECK(radial_mass_energy_pointwise(psi1, psi2, R, lambda) ==
         doctest::Approx(expected));
}
