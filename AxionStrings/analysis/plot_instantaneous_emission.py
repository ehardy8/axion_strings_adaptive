#!/usr/bin/env python3
"""Plot the instantaneous axion emission spectrum F(k/H, m_r/H), Fleury &
Moore 1806.04677 eq.(33) (2026-09-19, with the user).

F is extracted by comparing the R^3-weighted comoving axion spectrum at
pairs of snapshots separated by DELTA_LOG in log(m_r/H) -- see
axion_analysis.instantaneous_emission's docstring for the derivation and
docs/STATUS.md for the (corrected, previously wrong) drho_a/dk unit
conversion this builds on.

Usage:
    python3 plot_instantaneous_emission.py NETWORK_FILE SPECTRUM_FILE \
        OUT_DIR A_INV C0 N [DELTA_LOG] [TITLE_TAG]

DELTA_LOG (default 0.2) is the target spacing between the two snapshots
used for each finite difference -- deliberately a single, easy-to-change
number (per the user: "leave this as a choice that can easily be
modified"), not hardcoded inline in the plotting logic below.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from axion_analysis import (
    Background, L_tilde_general, instantaneous_emission, load_axion_spectrum,
    load_network_scalars, nearest_snapshot_for_delta_log,
)

NETWORK_FILE = sys.argv[1] if len(sys.argv) > 1 else "network_scalars.dat"
SPECTRUM_FILE = sys.argv[2] if len(sys.argv) > 2 else "axion_spectrum.dat"
OUT_DIR = sys.argv[3] if len(sys.argv) > 3 else "."
A_INV = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
C0 = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
N = int(sys.argv[6]) if len(sys.argv) > 6 else 256
DELTA_LOG = float(sys.argv[7]) if len(sys.argv) > 7 else 0.2
TITLE_TAG = sys.argv[8] if len(sys.argv) > 8 else f"N={N}"

bg = Background(A_INV, C0)
L_tilde = L_tilde_general(N, A_INV, C0)

spec = load_axion_spectrum(SPECTRUM_FILE)
spec_taus = np.unique(spec["tau"])

# Reference (later) snapshots to show F at: spread across the run, skipping
# the very first one or two (right at/near the pre-evolution handoff, which
# earlier plots showed to be a transient outlier, and which -- being the
# start of the series -- has no earlier snapshot far enough back anyway).
usable_taus = spec_taus[2:]
n_pick = min(6, len(usable_taus))
targets = usable_taus[np.linspace(0, len(usable_taus) - 1, n_pick).round().astype(int)]

fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis")
n_dropped_total = 0
n_points_total = 0

for idx, tau2 in enumerate(targets):
    tau1, achieved_delta_log = nearest_snapshot_for_delta_log(
        bg, spec_taus, tau2, DELTA_LOG
    )
    if tau1 == tau2:
        print(f"Skipping log(m_r/H)={bg.log_mr_over_h(tau2):.2f}: no earlier "
              f"snapshot available for Delta log={DELTA_LOG}")
        continue

    mask1 = spec["tau"] == tau1
    mask2 = spec["tau"] == tau2
    p1, S1 = spec["mode_index"][mask1], spec["shell_average_screened"][mask1]
    p2, S2 = spec["mode_index"][mask2], spec["shell_average_screened"][mask2]

    x, F, valid = instantaneous_emission(bg, N, L_tilde, tau1, p1, S1, tau2, p2, S2)
    n_dropped_total += np.sum(~valid)
    n_points_total += len(valid)

    label = f"log(m_r/H)={bg.log_mr_over_h(tau2):.2f} (Delta log={achieved_delta_log:.3f})"
    ax.loglog(x[valid], F[valid], "-", color=cmap(idx / max(1, n_pick - 1)),
              label=label, alpha=0.85)
    if np.any(~valid):
        # Show but visually de-emphasise the negative-derivative points
        # (1806.04677: "subject to fluctuations at frequencies near the
        # core") -- plotted as their absolute value, dotted, not connected
        # to the main curve, so they read as flagged rather than trusted.
        ax.loglog(x[~valid], np.abs(F[~valid]), ":", color=cmap(idx / max(1, n_pick - 1)),
                  alpha=0.3, linewidth=0.8)

    m_r_over_H_2 = bg.m_r(tau2) / bg.H(tau2)
    ax.axvline(m_r_over_H_2 / 2.0, color=cmap(idx / max(1, n_pick - 1)),
               linestyle="--", alpha=0.15, linewidth=1)

ax.set_xlabel(r"$k/H$")
ax.set_ylabel(r"$F(k/H,\, m_r/H)$")
ax.set_title(
    f"Instantaneous axion emission spectrum\n{TITLE_TAG}, screened, "
    f"target $\\Delta\\log(m_r/H)$={DELTA_LOG}\n"
    "(dotted = negative raw derivative, fluctuation-dominated -- 1806.04677 sec.4.2.1;\n"
    "dashed vertical lines = $m_r/2H$ UV cutoff, matching colour)"
)
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
out_path = f"{OUT_DIR}/instantaneous_emission_F.png"
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
if n_points_total:
    print(f"Negative (flagged) points: {n_dropped_total}/{n_points_total} "
         f"({100 * n_dropped_total / n_points_total:.1f}%)")
