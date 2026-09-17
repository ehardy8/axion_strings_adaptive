#ifndef XIFORMULA_HPP_
#define XIFORMULA_HPP_

// The xi <-> plaquette-count relation of docs/conventions.md sec.8:
//
//   xi = (2/3) N_p (dx / L_tilde^3) v^2 ((a_inv-1)^3/a_inv^2) (tau/tau0)^2
//
// Shared between task 1.6 (the string finder, which measures N_p and reports
// xi) and task 1.5 (pre-evolution, which needs the inverse: how many
// plaquettes give a target xi). v = 1, tau0 = 1 (sec.5 code units).
//
// AMReX-free like Background.hpp/BoxPlan.hpp -- pure math, unit-tested
// standalone.
//
// dx and L_tilde are the comoving lattice spacing and box size of whichever
// grid the plaquettes are counted on -- the pre-evolution grid in task 1.5
// (dx_pre = L_tilde_pre / N), or the main run's grid once task 1.6 lands.
// tau is the conformal time at which the count is taken; task 1.5 evaluates
// this at tau_i (the main run's start time -- sec.5 note: not necessarily
// tau0), since that is the physical instant the handed-off state represents.

#include <cmath>

[[nodiscard]] inline double xi_from_plaquette_count(double N_p, double dx,
                                                     double L_tilde,
                                                     double a_inv, double tau)
{
    const double a_inv_minus_1 = a_inv - 1.0;
    return (2.0 / 3.0) * N_p * (dx / (L_tilde * L_tilde * L_tilde)) *
           (a_inv_minus_1 * a_inv_minus_1 * a_inv_minus_1 / (a_inv * a_inv)) *
           (tau * tau);
}

// Inverse of the above, specialised to dx = L_tilde / N (N grid points per
// side): the plaquette count needed to reach xi_target at time tau.
[[nodiscard]] inline double
plaquette_count_for_target_xi(double xi_target, int N, double L_tilde,
                              double a_inv, double tau)
{
    const double a_inv_minus_1 = a_inv - 1.0;
    return 3.0 * N * xi_target * L_tilde * L_tilde * a_inv * a_inv /
           (2.0 * a_inv_minus_1 * a_inv_minus_1 * a_inv_minus_1 * tau * tau);
}

#endif // XIFORMULA_HPP_
