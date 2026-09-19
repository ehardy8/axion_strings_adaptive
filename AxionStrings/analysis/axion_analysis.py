"""Shared, reusable analysis utilities for AxionStrings output.

Consolidates what had been re-derived ad hoc in several one-off scratchpad
scripts this project (Background(tau), the v^-3 drho_a/dk unit conversion,
k/H) into one place, and adds the instantaneous-emission-spectrum
extraction F(k/H, m_r/H) of Fleury & Moore, 1806.04677 sec.4.2.1 eq.(33)
(2026-09-19, with the user: "all these types of plots are going to have to
be made many times over the course of this project, so let's have a nice
setup").

A fat->Moore switch (2026-09-19, with the user) is supported via the
optional c1/tau_switch fields, mirroring AxionStrings/Background.hpp's
CTauSchedule exactly (re-anchored lambda(tau) at the switch) -- see
detect_background() to build one of these automatically from a run's own
network_scalars.dat rather than needing the switch parameters passed in by
hand (easy to get wrong/stale relative to the actual run).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Background(tau) -- mirrors AxionStrings/Background.hpp exactly, including
# CTauSchedule's re-anchored lambda(tau) across an optional single switch.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Background:
    a_inv: float
    c0: float
    c1: Optional[float] = None
    tau_switch: Optional[float] = None

    @property
    def has_switch(self) -> bool:
        return self.tau_switch is not None

    @property
    def b_inv(self) -> float:
        return self.a_inv - 1.0

    @property
    def R0(self) -> float:
        return 1.0 / self.b_inv

    def R(self, tau):
        return self.R0 * (np.asarray(tau) / 1.0) ** (1.0 / self.b_inv)

    def R_prime(self, tau):
        return (self.R0 / (self.b_inv * 1.0)) * (np.asarray(tau) / 1.0) ** (
            1.0 / self.b_inv - 1.0
        )

    def lam(self, tau):
        tau = np.asarray(tau, dtype=float)
        pre_switch = 1.0 * (self.R(tau) / self.R0) ** (-2.0 * self.c0)
        if not self.has_switch:
            return pre_switch
        R_s = self.R(self.tau_switch)
        lambda_s = 1.0 * (R_s / self.R0) ** (-2.0 * self.c0)
        post_switch = lambda_s * (self.R(tau) / R_s) ** (-2.0 * self.c1)
        return np.where(tau < self.tau_switch, pre_switch, post_switch)

    def H(self, tau):
        return self.R_prime(tau) / self.R(tau) ** 2

    def m_r(self, tau):
        return np.sqrt(self.lam(tau))

    def log_mr_over_h(self, tau):
        return np.log(self.m_r(tau) / self.H(tau))

    def tau_from_log_mr_over_h(self, log_mr_over_h):
        """Inverse of the no-switch H_over_mr closed form -- only valid
        for resolving a time *before* c_sched.c0 stops applying (mirrors
        Background::tau_from_log_mr_over_h's own caveat)."""
        return np.exp(log_mr_over_h * (self.a_inv - 1.0) / (self.a_inv - self.c0))

    def t_cosmic(self, tau):
        """Cosmic time t = integral of R(tau') dtau', in closed form for
        this power-law R(tau) = R0 (tau/tau0)^(1/b_inv):
            t(tau) = tau^(a_inv/b_inv) / a_inv + const.
        The additive constant is arbitrary (fixed here so t_cosmic(1)=0)
        and irrelevant for any use of this function, since only time
        *differences* are ever physically meaningful (e.g. the eq.(33)
        finite difference below) -- verified against direct numerical
        quadrature of R(tau) in test_axion_analysis.py. Not switch-aware
        (matches every use so far -- the eq.33 finite difference is always
        taken within a single phase); extend if that ever changes.
        """
        tau = np.asarray(tau)
        return (tau ** (self.a_inv / self.b_inv) - 1.0) / self.a_inv


def detect_background(network_scalars, a_inv, c0):
    """Build a Background from a1_inv/c0 plus, if network_scalars shows a
    plateau in m_r_over_H (the unmistakable signature of a fat->Moore
    switch -- m_r/H is frozen by construction from the switch onward),
    the c1=a_inv/tau_switch that reproduces it -- rather than needing the
    switch parameters passed in by hand and kept in sync with the actual
    run. tau_switch is solved exactly from the plateau's own value (not
    just estimated): the plateaued m_r_over_H is exactly what the *no*
    switch formula would have given at the true tau_switch, so inverting
    that formula (tau_from_log_mr_over_h) recovers it exactly, for any
    c0 -- not just the c0=1 fat case where tau_switch happens to equal
    gamma.
    """
    tau = np.asarray(network_scalars["tau"])
    m_r_over_h = np.asarray(network_scalars["m_r_over_H"])
    if len(tau) < 2:
        return Background(a_inv, c0)

    diffs = np.abs(np.diff(m_r_over_h))
    plateaued = diffs < 1.0e-9 * np.maximum(np.abs(m_r_over_h[:-1]), 1.0)
    if not np.any(plateaued):
        return Background(a_inv, c0)

    i = int(np.argmax(plateaued))  # first index where m_r_over_H[i]==m_r_over_H[i+1]
    gamma = float(m_r_over_h[i])
    no_switch = Background(a_inv, c0)
    tau_switch = float(no_switch.tau_from_log_mr_over_h(np.log(gamma)))
    return Background(a_inv, c0, c1=a_inv, tau_switch=tau_switch)


def L_tilde_general(N: int, a_inv: float, c0: float) -> float:
    """Comoving box length for a general (non-Moore) box plan (BoxPlan.hpp),
    N1=N2=1 (this project's production convention throughout)."""
    R0 = (a_inv - 1.0) ** -1
    return (1.0 / R0) * float(N) ** ((a_inv - 1.0) / (a_inv - c0))


F_A = np.sqrt(2.0)  # f_a = sqrt(2) v, v = 1 in code units (Energy.hpp)


# ---------------------------------------------------------------------------
# File loading -- named fields instead of magic column indices, so a future
# change to network_scalars.dat's column layout (it has already changed
# once this project, when the radial energy components were added) can't
# silently desync a script that indexes by position.
# ---------------------------------------------------------------------------
_NETWORK_COLUMNS = [
    "tau", "N_p", "N_p_weighted", "xi", "xi_weighted", "m_r_over_H",
    "rho_tot_unscreened", "rho_tot_screened",
    "rho_axion_kin_unscreened", "rho_axion_kin_screened",
    "rho_axion_grad_unscreened", "rho_axion_grad_screened",
    "rho_radial_kin_unscreened", "rho_radial_kin_screened",
    "rho_radial_grad_unscreened", "rho_radial_grad_screened",
    "rho_radial_mass_unscreened", "rho_radial_mass_screened",
    "n_total", "n_unmasked", "mean_gamma_sq_v_sq", "mean_gamma",
    "n_velocity_corners", "tension_core_only", "tension_core_plus_tail",
]

_SPECTRUM_COLUMNS = [
    "tau", "mode_index",
    "shell_average_screened", "shell_average_unscreened",
    "full_cube_energy_screened", "full_cube_energy_unscreened",
    "inscribed_sphere_energy_screened", "inscribed_sphere_energy_unscreened",
    "real_space_mean_sq_screened", "real_space_mean_sq_unscreened",
    "parseval_full_screened", "parseval_full_unscreened",
    "parseval_inscribed_screened", "parseval_inscribed_unscreened",
]


def _load_named(path, expected_columns):
    data = np.loadtxt(path, comments="#")
    if data.shape[1] != len(expected_columns):
        raise ValueError(
            f"{path}: expected {len(expected_columns)} columns "
            f"{expected_columns}, got {data.shape[1]} -- header layout has "
            "likely changed; update the column list in axion_analysis.py."
        )
    return {name: data[:, i] for i, name in enumerate(expected_columns)}


def load_network_scalars(path):
    return _load_named(path, _NETWORK_COLUMNS)


def load_axion_spectrum(path):
    return _load_named(path, _SPECTRUM_COLUMNS)


# ---------------------------------------------------------------------------
# Spectrum unit conversion (2026-09-19 fix -- see docs/STATUS.md): the FFT
# buffer holds theta' (comoving dtheta/dtau), not the physical a_dot =
# f_a*theta'/R, so converting the raw shell_average S(p) to the physical,
# dimensionless v^-3 drho_a/dk needs an explicit (f_a/R)^2 factor that
# conventions.md sec.10's documented rescale omits (and gets R's power
# wrong besides). Verified against the independently-validated
# rho_axion_kinetic_screened (network_scalars.dat): integral(drho_a/dk)/
# (2*rho_axion_kinetic_screened) -> 1.00 +/- 0.02 for log(m_r/H) >~ 4.
# ---------------------------------------------------------------------------
def v3_drho_dk(S, R, L_tilde, N, f_a=F_A):
    """v^-3 drho_a/dk(p), from the raw shell_average S(p)."""
    return S * f_a**2 * L_tilde / (2.0 * np.pi * R * N**6)


def k_over_H(p, R, H, L_tilde):
    """k/H for comoving mode index p, at the physical box length R*L_tilde."""
    L_phys = R * L_tilde
    return p * (2.0 * np.pi / (L_phys * H))


def comoving_g(S, R, L_tilde, N, f_a=F_A):
    """R^3 * (v^-3 drho_a/dk) -- the R^3-weighted comoving spectral energy
    density tracked at fixed comoving mode number p (equivalently, fixed
    comoving wavenumber kappa = k*R) across time, per Fleury & Moore
    1806.04677 eq.(33)'s derivation (see instantaneous_emission's
    docstring for the full derivation)."""
    return R**3 * v3_drho_dk(S, R, L_tilde, N, f_a)


# ---------------------------------------------------------------------------
# Instantaneous emission spectrum F(k/H, m_r/H), 1806.04677 eq.(33):
#   F(k/H, m_r/H) = (A/R^3) d/dt [R^3 drho_a/dk],  A = H/Gamma fixed by
#   normalising integral(F dx) = 1 (eq.22), x = k/H.
# ---------------------------------------------------------------------------
def instantaneous_emission(bg: Background, N, L_tilde,
                           tau1, p1, S1, tau2, p2, S2, f_a=F_A):
    """Extract F(k/H) from two spectrum snapshots (tau1 < tau2, same run).

    Derivation (not in the paper, done here to make the implementation
    checkable): the paper's eq.(23) gives drho_a/dk[k,t] as a time integral
    over the *redshifted* momentum k' = k*R(t)/R(t'). Substituting the
    comoving wavenumber kappa = k*R(t) (constant for a fixed Fourier mode
    -- k' = kappa/R(t') is then just that same mode's *own* physical
    momentum at t'), the eq.(23) integral becomes, at fixed kappa,

        G(kappa, t) = integral_t0^t dt' (Gamma'/H') R'^3 F(kappa/(R'H'), ...)

    whose integrand no longer depends on the upper limit t, so
    d(G)/dt|_kappa = (Gamma(t)/H(t)) R(t)^3 F(k/H, m_r/H) exactly, with
    k = kappa/R(t) the *current* physical momentum -- i.e. eq.(33), with
    the derivative taken *at fixed comoving mode index p* (kappa = 2*pi*p/
    L_tilde, L_tilde fixed for a run), not at fixed physical k. This is
    exactly "redshift the earlier snapshot's R^3-weighted spectrum forward
    before differencing" -- comoving_g() already does the R^3 weighting;
    matching by mode index p (not k or k/H) *is* the redshift-consistent
    comparison, so no additional interpolation in k is needed or correct.

    Normalises the raw derivative to unit integral over x=k/H directly
    (mirrors the paper's own approach -- "the factor A=H/Gamma is fixed by
    requiring that F is normalised to 1" -- so neither Gamma(t)'s often
    complicated closed form (eq.17) nor any overall constant in
    v3_drho_dk's own definition need to be right for F's *shape*; only the
    R^3 weighting, the p-matching, and the cosmic-time Delta t matter).

    Parameters
    ----------
    tau1, p1, S1 : earlier snapshot's tau and (mode_index, shell_average)
        arrays (screened, typically -- axion emission, not string cores).
    tau2, p2, S2 : later snapshot, same run, same quantities.

    Returns
    -------
    x : k/H, evaluated at tau2 (the pair's later/reference time).
    F : normalised instantaneous emission spectrum at each x.
    valid : boolean mask, False where the raw (unnormalised) derivative
        came out negative -- expected occasionally near the core scale
        (1806.04677: "interactions with radial modes induce small
        oscillations... subject to fluctuations at frequencies near the
        core"), not necessarily a bug; plot/interpret these with care
        rather than silently dropping or trusting them.
    """
    if not (np.array_equal(p1, p2)):
        # Match on the common mode indices rather than assuming identical
        # arrays -- the two snapshots can in principle have different
        # inscribed-sphere cutoffs if N ever changes mid-project.
        common = np.intersect1d(p1, p2)
        i1 = np.searchsorted(p1, common)
        i2 = np.searchsorted(p2, common)
        p1, S1 = p1[i1], S1[i1]
        p2, S2 = p2[i2], S2[i2]
        p = common
    else:
        p = p1

    R1, R2 = bg.R(tau1), bg.R(tau2)
    g1 = comoving_g(S1, R1, L_tilde, N, f_a)
    g2 = comoving_g(S2, R2, L_tilde, N, f_a)

    dt = bg.t_cosmic(tau2) - bg.t_cosmic(tau1)
    if dt <= 0:
        raise ValueError(f"tau2 ({tau2}) must be later than tau1 ({tau1})")

    F_raw = (g2 - g1) / dt
    valid = F_raw > 0.0

    H2 = bg.H(tau2)
    x = k_over_H(p, R2, H2, L_tilde)

    # Normalise to integral(F dx) = 1 -- x is uniformly spaced in p (dx =
    # 2*pi/(L_phys*H) is the same constant for every mode at this
    # snapshot), so a plain Riemann sum is exact for this grid.
    dx = x[1] - x[0] if len(x) > 1 else 1.0
    area = np.sum(np.where(valid, F_raw, 0.0)) * dx
    if area <= 0:
        raise ValueError(
            "Instantaneous emission spectrum has non-positive total area -- "
            "Delta log is probably too small (dominated by fluctuation "
            "noise); try a larger spacing between tau1 and tau2."
        )
    F = F_raw / area

    return x, F, valid


def nearest_snapshot_for_delta_log(bg: Background, spectrum_taus, tau_ref,
                                   delta_log):
    """Pick the snapshot tau1 (from the available, discrete spectrum_taus)
    whose log(m_r/H) is closest to log(m_r/H)(tau_ref) - delta_log. Returns
    (tau1, achieved_delta_log) -- the achieved value can differ from the
    requested one since snapshots only exist at the run's own output
    cadence; always report which was actually used, don't assume it."""
    target_log = bg.log_mr_over_h(tau_ref) - delta_log
    taus = np.asarray(sorted(set(spectrum_taus)))
    logs = bg.log_mr_over_h(taus)
    idx = np.argmin(np.abs(logs - target_log))
    tau1 = taus[idx]
    achieved = bg.log_mr_over_h(tau_ref) - logs[idx]
    return tau1, achieved
