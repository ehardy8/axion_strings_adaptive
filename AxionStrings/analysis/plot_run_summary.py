#!/usr/bin/env python3
"""Standard per-run summary plots for AxionStrings output: xi, energies
(total/axion/radial, raw and H^2 f_a^2-normalised), string tension, and the
axion power spectrum vs k/H (raw and H f_a^2-normalised). Consolidates what
had been a one-off scratchpad script (2026-09-18/19) into the reusable
`analysis/` setup (2026-09-19, with the user: "let's get the full spectrum,
energy and tension plots... into a proper reusable python script as well"),
built on axion_analysis.py rather than re-deriving Background(tau)/unit
conversions locally again.

Usage:
    python3 plot_run_summary.py NETWORK_FILE SPECTRUM_FILE OUT_DIR \
        A_INV C0 N [TITLE_TAG]

See plot_instantaneous_emission.py for the separate F(k/H, m_r/H)
instantaneous-emission-spectrum plot.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from axion_analysis import (
    F_A, Background, L_tilde_general, k_over_H, load_axion_spectrum,
    load_network_scalars, v3_drho_dk,
)

NETWORK_FILE = sys.argv[1] if len(sys.argv) > 1 else "network_scalars.dat"
SPECTRUM_FILE = sys.argv[2] if len(sys.argv) > 2 else "axion_spectrum.dat"
OUT_DIR = sys.argv[3] if len(sys.argv) > 3 else "."
A_INV = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
C0 = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
N = int(sys.argv[6]) if len(sys.argv) > 6 else 256
TITLE_TAG = sys.argv[7] if len(sys.argv) > 7 else f"N={N}"

bg = Background(A_INV, C0)
L_TILDE = L_tilde_general(N, A_INV, C0)

net = load_network_scalars(NETWORK_FILE)
tau = net["tau"]
log_mr_over_h = bg.log_mr_over_h(tau)

# Cross-check: m_r/H from the file vs. this module's own Background(tau).
mrh_check = bg.m_r(tau) / bg.H(tau)
print(
    "m_r/H cross-check (file vs. Background(tau)): max relative diff = "
    f"{np.max(np.abs(mrh_check - net['m_r_over_H']) / net['m_r_over_H']):.3e}"
)
print(f"L_tilde (recomputed) = {L_TILDE:.6g}")

# ---------------------------------------------------------------------------
# xi vs log(m_r/H)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(log_mr_over_h, net["xi"], "o-", color="#2b6cb0", markersize=4, linewidth=1.5)
ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.6, label=r"$\xi=1$ (attractor)")
ax.set_xlabel(r"$\log(m_r/H)$")
ax.set_ylabel(r"$\xi$")
ax.set_title(f"String network density vs. log(m_r/H)\n{TITLE_TAG}, screened (scheme B)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/xi_vs_log_mr_over_h.png", dpi=150)
print(f"Wrote {OUT_DIR}/xi_vs_log_mr_over_h.png")


def make_energy_panels(ax_tot, ax_axion, ax_radial, normalise):
    scale = (bg.H(tau) ** 2 * F_A**2) if normalise else 1.0
    unit = r"$H^2 f_a^2$" if normalise else ""

    ax_tot.plot(log_mr_over_h, net["rho_tot_unscreened"] / scale, "o-", label=r"$\rho_{tot}$ unscreened", color="#c53030")
    ax_tot.plot(log_mr_over_h, net["rho_tot_screened"] / scale, "s--", label=r"$\rho_{tot}$ screened", color="#e53e3e")
    ax_tot.set_yscale("log")
    ax_tot.set_xlabel(r"$\log(m_r/H)$")
    ax_tot.set_ylabel(rf"$\rho$ [{unit}]" if normalise else r"$\rho$")
    ax_tot.set_title("Total energy density")
    ax_tot.legend(fontsize=8)
    ax_tot.grid(True, alpha=0.3)

    ax_axion.plot(log_mr_over_h, net["rho_axion_kin_unscreened"] / scale, "o-", label="axion kin. unscreened", color="#2b6cb0")
    ax_axion.plot(log_mr_over_h, net["rho_axion_kin_screened"] / scale, "s--", label="axion kin. screened", color="#4299e1")
    ax_axion.plot(log_mr_over_h, net["rho_axion_grad_unscreened"] / scale, "o-", label="axion grad. unscreened", color="#2f855a")
    ax_axion.plot(log_mr_over_h, net["rho_axion_grad_screened"] / scale, "s--", label="axion grad. screened", color="#48bb78")
    ax_axion.set_yscale("log")
    ax_axion.set_xlabel(r"$\log(m_r/H)$")
    ax_axion.set_ylabel(rf"$\rho$ [{unit}]" if normalise else r"$\rho$")
    ax_axion.set_title("Axion (phase) kinetic / gradient energy")
    ax_axion.legend(fontsize=7)
    ax_axion.grid(True, alpha=0.3)

    ax_radial.plot(log_mr_over_h, net["rho_radial_kin_unscreened"] / scale, "o-", label="radial kin. unscreened", color="#b7791f")
    ax_radial.plot(log_mr_over_h, net["rho_radial_kin_screened"] / scale, "s--", label="radial kin. screened", color="#ecc94b")
    ax_radial.plot(log_mr_over_h, net["rho_radial_grad_unscreened"] / scale, "o-", label="radial grad. unscreened", color="#805ad5")
    ax_radial.plot(log_mr_over_h, net["rho_radial_grad_screened"] / scale, "s--", label="radial grad. screened", color="#b794f4")
    ax_radial.plot(log_mr_over_h, net["rho_radial_mass_unscreened"] / scale, "o-", label="radial mass unscreened", color="#c53030")
    ax_radial.plot(log_mr_over_h, net["rho_radial_mass_screened"] / scale, "s--", label="radial mass screened", color="#fc8181")
    ax_radial.set_yscale("log")
    ax_radial.set_xlabel(r"$\log(m_r/H)$")
    ax_radial.set_ylabel(rf"$\rho$ [{unit}]" if normalise else r"$\rho$")
    ax_radial.set_title("Radial (Higgs) kinetic / gradient / mass energy")
    ax_radial.legend(fontsize=7)
    ax_radial.grid(True, alpha=0.3)


fig, axes = plt.subplots(1, 3, figsize=(19, 5))
make_energy_panels(*axes, normalise=False)
fig.suptitle(f"{TITLE_TAG} -- energies (raw code units)")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/energies_vs_log_mr_over_h.png", dpi=150)
print(f"Wrote {OUT_DIR}/energies_vs_log_mr_over_h.png")

fig, axes = plt.subplots(1, 3, figsize=(19, 5))
make_energy_panels(*axes, normalise=True)
fig.suptitle(f"{TITLE_TAG} -- energies normalised by $H^2 f_a^2$")
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/energies_vs_log_mr_over_h_normalized.png", dpi=150)
print(f"Wrote {OUT_DIR}/energies_vs_log_mr_over_h_normalized.png")

# ---------------------------------------------------------------------------
# String tension -- left in raw (v^2) units; H^2 f_a^2 has different
# dimensions (mass^4) to tension (mass^2, energy/length), so the energy
# normalisation above does not apply here.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(log_mr_over_h, net["tension_core_only"], "o-", label="core only", color="#805ad5")
ax.plot(log_mr_over_h, net["tension_core_plus_tail"], "s--", label="core + tail", color="#b794f4")
ax.set_xlabel(r"$\log(m_r/H)$")
ax.set_ylabel(r"tension $\mu$")
ax.set_title(f"{TITLE_TAG} -- string tension vs. log(m_r/H)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/tension_vs_log_mr_over_h.png", dpi=150)
print(
    f"Wrote {OUT_DIR}/tension_vs_log_mr_over_h.png (core-only range: "
    f"[{np.nanmin(net['tension_core_only']):.3g}, {np.nanmax(net['tension_core_only']):.3g}])"
)

# ---------------------------------------------------------------------------
# Axion spectrum at several times, x-axis k/H, y-axis v^-3 drho_a/dk (raw
# and H f_a^2-normalised) -- see axion_analysis.v3_drho_dk/docs/STATUS.md
# for the 2026-09-19 unit-conversion fix this builds on.
# ---------------------------------------------------------------------------
spec = load_axion_spectrum(SPECTRUM_FILE)
spec_tau, mode_index, shell_screened = spec["tau"], spec["mode_index"], spec["shell_average_screened"]

unique_tau = np.unique(spec_tau)
n_pick = min(6, len(unique_tau))
picks = unique_tau[np.linspace(0, len(unique_tau) - 1, n_pick).round().astype(int)]


def make_spectrum_plot(normalise):
    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.get_cmap("viridis")
    for idx, t in enumerate(picks):
        mask = spec_tau == t
        p, S = mode_index[mask], shell_screened[mask]
        Rt, Ht = bg.R(t), bg.H(t)
        x = k_over_H(p, Rt, Ht, L_TILDE)
        y = v3_drho_dk(S, Rt, L_TILDE, N)
        if normalise:
            # f_a^2 cancels analytically -- see docs/STATUS.md; the whole
            # point of this normalisation is to remove the (otherwise
            # free) f_a dependence from the plotted spectral shape.
            y = S * L_TILDE / (2.0 * np.pi * Rt * Ht * N**6)
        good = (p > 0) & (S > 0)
        ax.loglog(x[good], y[good], "-", color=cmap(idx / max(1, n_pick - 1)),
                  label=f"log(m_r/H)={bg.log_mr_over_h(t):.2f}", alpha=0.85)
    ax.set_xlabel(r"$k/H$")
    ylabel = r"$(\partial\rho_a/\partial k) / (H f_a^2)$" if normalise else r"$v^{-3}\,\partial\rho_a/\partial k$"
    ax.set_ylabel(ylabel)
    suffix = "normalised by $H f_a^2$" if normalise else "raw code units"
    ax.set_title(f"Axion (screened) power spectrum at several times\n{TITLE_TAG}, {suffix}")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    suffix_path = "_normalized" if normalise else ""
    out_path = f"{OUT_DIR}/axion_spectrum_vs_k_over_H{suffix_path}.png"
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path} ({n_pick} snapshots: {picks})")


make_spectrum_plot(normalise=False)
make_spectrum_plot(normalise=True)
