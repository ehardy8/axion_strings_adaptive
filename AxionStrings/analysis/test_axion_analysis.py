"""Tests for axion_analysis.py. Plain assert-based (no pytest dependency,
matching the project's lightweight doctest-suite style for the C++ side) --
run directly with `python3 test_axion_analysis.py`.
"""
import numpy as np
from scipy import integrate

from axion_analysis import (
    Background, L_tilde_general, comoving_g, detect_background,
    instantaneous_emission, k_over_H, nearest_snapshot_for_delta_log,
    v3_drho_dk,
)


def test_t_cosmic_matches_quadrature():
    """t_cosmic's closed form vs. direct numerical quadrature of R(tau)."""
    for a_inv, c0 in [(2.0, 1.0), (2.0, 0.0), (3.0, 0.0), (1.5, 1.0)]:
        bg = Background(a_inv, c0)
        tau_a, tau_b = 5.0, 12.3
        expected, _ = integrate.quad(bg.R, tau_a, tau_b)
        got = bg.t_cosmic(tau_b) - bg.t_cosmic(tau_a)
        assert abs(got - expected) / expected < 1e-8, (a_inv, c0, got, expected)
    print("PASS: t_cosmic matches quadrature of R(tau)")


def test_v3_drho_dk_and_comoving_g_consistency():
    """comoving_g(S, R, ...) == R^3 * v3_drho_dk(S, R, ...) by construction."""
    S = np.array([0.3, 1.7, 4.2])
    R, L_tilde, N = 3.5, 256.0, 128
    lhs = comoving_g(S, R, L_tilde, N)
    rhs = R**3 * v3_drho_dk(S, R, L_tilde, N)
    assert np.allclose(lhs, rhs)
    print("PASS: comoving_g == R^3 * v3_drho_dk")


def test_k_over_h_scales_correctly():
    p = np.array([1.0, 2.0, 10.0])
    R, H, L_tilde = 4.0, 0.1, 256.0
    x = k_over_H(p, R, H, L_tilde)
    # linear in p, and doubling R (at fixed H, L_tilde) halves x (k/H ~ 1/R
    # through the physical box length L = R*L_tilde).
    x2 = k_over_H(p, 2 * R, H, L_tilde)
    assert np.allclose(x2, x / 2.0)
    assert np.allclose(x / p, x[0] / p[0])  # exactly linear in p
    print("PASS: k_over_H scales as expected with p and R")


def test_background_lam_continuous_and_frozen_across_a_switch():
    """lam(tau)'s re-anchored switch branch (mirrors Background.hpp's
    CTauSchedule) must be continuous at tau_switch, and log_mr_over_h must
    freeze at gamma from the switch onward -- checked against independent
    hand computation, not just internal self-consistency."""
    a_inv, c0, tau_switch = 2.0, 1.0, 20.0
    bg = Background(a_inv, c0, c1=a_inv, tau_switch=tau_switch)

    eps = 1.0e-6
    lam_before = bg.lam(tau_switch - eps)
    lam_after = bg.lam(tau_switch + eps)
    assert abs(lam_before - lam_after) / lam_before < 1.0e-4, (lam_before, lam_after)

    gamma = tau_switch  # exact identity for c0=1 (fat), derived by hand
    for tau in [tau_switch, tau_switch * 1.5, tau_switch * 3.0]:
        got = bg.log_mr_over_h(tau)
        assert abs(got - np.log(gamma)) < 1.0e-8, (tau, got, np.log(gamma))

    # Before the switch, log_mr_over_h must match the plain no-switch
    # background exactly (independent construction, not the same object).
    bg_no_switch = Background(a_inv, c0)
    for tau in [1.0, 5.0, tau_switch - 0.001]:
        got, expected = bg.log_mr_over_h(tau), bg_no_switch.log_mr_over_h(tau)
        assert abs(got - expected) < 1.0e-9 * max(abs(expected), 1.0), (tau, got, expected)
    print("PASS: Background.lam is continuous at the switch and frozen after it")


def test_detect_background_recovers_a_known_switch():
    """Synthesise network_scalars-like (tau, m_r_over_H) rows for a real
    fat->Moore run (growing then frozen at gamma), and check
    detect_background recovers the exact tau_switch used to build it --
    an independent construction (built from the *true* switch-aware
    Background), not a restatement of detect_background's own logic."""
    a_inv, c0 = 2.0, 1.0
    true_tau_switch = 37.5
    bg_true = Background(a_inv, c0, c1=a_inv, tau_switch=true_tau_switch)

    tau = np.concatenate([np.linspace(5.0, true_tau_switch, 10),
                          np.linspace(true_tau_switch + 2.0, 80.0, 10)])
    m_r_over_h = bg_true.m_r(tau) / bg_true.H(tau)

    network_scalars = {"tau": tau, "m_r_over_H": m_r_over_h}
    bg_detected = detect_background(network_scalars, a_inv, c0)

    assert bg_detected.has_switch
    assert abs(bg_detected.tau_switch - true_tau_switch) / true_tau_switch < 1.0e-6
    print(
        f"PASS: detect_background recovers tau_switch={bg_detected.tau_switch:.4f} "
        f"(true={true_tau_switch})"
    )


def test_detect_background_no_switch_when_no_plateau():
    a_inv, c0 = 2.0, 0.0
    bg_true = Background(a_inv, c0)
    tau = np.linspace(2.0, 50.0, 15)
    m_r_over_h = bg_true.m_r(tau) / bg_true.H(tau)
    network_scalars = {"tau": tau, "m_r_over_H": m_r_over_h}
    bg_detected = detect_background(network_scalars, a_inv, c0)
    assert not bg_detected.has_switch
    print("PASS: detect_background finds no switch when there is none")


def test_nearest_snapshot_for_delta_log():
    bg = Background(2.0, 1.0)  # R(tau) = tau, log(m_r/H) = log(tau)
    available_taus = np.exp(np.arange(3.0, 6.0, 0.1))  # log spaced by 0.1
    tau_ref = available_taus[-1]
    tau1, achieved = nearest_snapshot_for_delta_log(
        bg, available_taus, tau_ref, delta_log=0.25
    )
    assert abs(achieved - 0.25) < 0.06, achieved  # snaps to nearest 0.1 grid
    assert tau1 in available_taus
    print(f"PASS: nearest_snapshot_for_delta_log (achieved Delta log={achieved:.4f})")


def test_instantaneous_emission_recovers_known_F():
    """End-to-end check of eq.(33)'s implementation: construct drho_a/dk by
    directly, numerically integrating the *defining* relation (Fleury &
    Moore 1806.04677 eq.23, at fixed comoving mode -- see
    instantaneous_emission's docstring) for a known, simple F_test(x) and
    Gamma(t)/H(t) = const, then check instantaneous_emission recovers
    F_test to good accuracy from two closely-spaced snapshots. This is an
    independent construction (direct quadrature of the defining integral),
    not a restatement of instantaneous_emission's own formula, so it
    actually exercises the R^3-weighting/comoving-matching/Delta-t logic
    rather than just checking self-consistency.
    """
    a_inv, c0 = 2.0, 1.0  # fat string: R(tau) = tau
    bg = Background(a_inv, c0)
    N = 256
    L_tilde = L_tilde_general(N, a_inv, c0)

    x0, sigma = 8.0, 3.0

    def F_test(x):
        return np.exp(-0.5 * ((x - x0) / sigma) ** 2)

    # Normalise F_test over the actual x-grid used below, so the recovered
    # F should match on a like-for-like footing (same discretisation).
    p_grid = np.arange(1, 200)
    tau_probe = 60.0  # arbitrary; only used to build a representative x-grid
    x_grid = k_over_H(p_grid, bg.R(tau_probe), bg.H(tau_probe), L_tilde)
    dx_grid = x_grid[1] - x_grid[0]
    norm = np.sum(F_test(x_grid)) * dx_grid
    F_test_normalised = lambda x: F_test(x) / norm  # noqa: E731

    GAMMA_OVER_H = 1.0  # constant Gamma(t)/H(t), code units

    def integrand(tau_prime, p):
        R_p = bg.R(tau_prime)
        H_p = bg.H(tau_prime)
        kappa = 2.0 * np.pi * p / L_tilde
        x_p = kappa / (R_p * H_p)
        # dt' = R' dtau' (cosmic from conformal); integrand of eq.(23) at
        # fixed kappa is (Gamma'/H') R'^3 F(x'), and dt'/dtau' = R'.
        return GAMMA_OVER_H * R_p**3 * F_test_normalised(x_p) * R_p

    tau0 = 1.0
    tau2 = 60.0
    p_test = np.arange(1, 60)  # away from the very lowest/highest modes

    def max_rel_err(tau1):
        S1 = np.empty_like(p_test, dtype=float)
        S2 = np.empty_like(p_test, dtype=float)
        for i, p in enumerate(p_test):
            g1, _ = integrate.quad(integrand, tau0, tau1, args=(p,), limit=200)
            g2, _ = integrate.quad(integrand, tau0, tau2, args=(p,), limit=200)
            # Convert G(kappa, tau) back to shell_average S(p, tau):
            # G = S*f_a^2*L_tilde*R^2/(2*pi*N^6)  (comoving_g o v3_drho_dk)
            S1[i] = g1 * 2.0 * np.pi * N**6 / (2.0 * L_tilde * bg.R(tau1) ** 2)
            S2[i] = g2 * 2.0 * np.pi * N**6 / (2.0 * L_tilde * bg.R(tau2) ** 2)

        x, F, valid = instantaneous_emission(
            bg, N, L_tilde, tau1, p_test, S1, tau2, p_test, S2
        )
        expected = F_test_normalised(x)
        # Compare only where F_test itself is not negligible: a numerical
        # derivative over a *finite* Delta log necessarily smears the
        # recovered F across a range of x rather than sampling one instant,
        # and F_test's Gaussian tails fall off so fast that even a tiny
        # smear blows up the *relative* error there -- not a bug (checked
        # below: the error shrinks with Delta log, the correct diagnostic
        # for "is this converging to the right continuum answer").
        mask = valid & (expected > 1e-2 * expected.max())
        assert mask.sum() > 10, "too few usable points for a meaningful check"
        return np.max(np.abs(F[mask] - expected[mask]) / expected[mask])

    # Convergence check: a real bug (wrong power of R, wrong Delta t, wrong
    # normalisation) would not shrink away as Delta log -> 0; finite-
    # difference truncation error does, and should shrink roughly linearly
    # (eq.33's derivative is a first-order-accurate finite difference).
    err_coarse = max_rel_err(tau1=55.0)   # Delta log(m_r/H) = 0.087
    err_fine = max_rel_err(tau1=59.5)     # Delta log(m_r/H) = 0.0084, ~10x smaller
    assert err_fine < err_coarse, (err_coarse, err_fine)
    assert err_fine < 0.15, f"error {err_fine:.4f} still too large at small Delta log"
    print(
        f"PASS: instantaneous_emission converges to a known F(x) as Delta log "
        f"shrinks (max rel err {err_coarse:.3f} -> {err_fine:.3f})"
    )


if __name__ == "__main__":
    test_t_cosmic_matches_quadrature()
    test_v3_drho_dk_and_comoving_g_consistency()
    test_k_over_h_scales_correctly()
    test_background_lam_continuous_and_frozen_across_a_switch()
    test_detect_background_recovers_a_known_switch()
    test_detect_background_no_switch_when_no_plateau()
    test_nearest_snapshot_for_delta_log()
    test_instantaneous_emission_recovers_known_F()
    print("\nAll tests passed.")
