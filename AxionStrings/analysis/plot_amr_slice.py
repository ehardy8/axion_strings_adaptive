#!/usr/bin/env python3
"""Plot an AxionStrings AMReX plotfile's energy density with the AMR grid
structure overlaid, so refinement can be inspected visually (2026-09-19,
requested by the user: "a picture of the energy density at a fairly late
timeshot... overlay the grid... to show refinement"). Two modes:

- "slice": a single-plane cut (fast, but only shows the grid at that one
  plane -- a string segment elsewhere along the line of sight, and the
  refinement following it, would not appear at all).
- "projection" (default, per the user's follow-up: "more like the
  projection we had before, with the grid corresponding to the point
  where the projected density originates"): a max-intensity projection
  along the line of sight, matching `ProjectionKernel.hpp`'s own "max of
  unscreened rho_tot per column" convention. `annotate_grids()` on a yt
  ProjectionPlot draws each grid box's (x,y) *footprint*, i.e. the union
  of refined regions anywhere along each column's full line of sight --
  not literally "the one cell the column's max came from" (yt has no
  built-in per-pixel provenance for that), but the natural and correct
  generalisation: it shows where finer resolution was available anywhere
  that could have contributed to the projected value at that point.

Needs `yt` (not part of the project's normal analysis/ dependencies,
since it is only used for this one-off AMR-structure visualisation, not
the routine per-run diagnostics in axion_analysis.py/plot_run_summary.py).
Install into a virtualenv, e.g.:

    python3 -m venv .venv-yt && source .venv-yt/bin/activate
    pip install yt

Usage:
    python3 plot_amr_slice.py PLOTFILE_DIR OUT_PNG \
        A_INV C0 TAU_I_POSTHANDOFF [AXIS] [MODE] [TITLE_TAG]

PLOTFILE_DIR is an AMReX plotfile directory (e.g. "plots/plt00316"), as
written by `amr.plot_int`. TAU_I_POSTHANDOFF is `s_tau_i` after any
pre-evolution handoff (`tau_i - tau_pre_end`; 0 if the run never went
through Phase::Relaxing) -- read it from the run's own
"pre-evolution -> main handoff" log line, or pass `tau_i` directly for a
run with no relaxation phase. AXIS defaults to "z"; MODE is "projection"
(default) or "slice".

Only a_inv=2, c0=0 (physical mode, no switch) is implemented for the
R(tau)/lambda(tau) background needed to compute rho_tot -- this is a
visualisation tool for a specific test, not (yet) a general Background
reimplementation like axion_analysis.py's.
"""
import sys

import numpy as np
import yt

PLOTFILE_DIR = sys.argv[1]
OUT_PNG = sys.argv[2]
A_INV = float(sys.argv[3])
C0 = float(sys.argv[4])
TAU_I_POSTHANDOFF = float(sys.argv[5])
AXIS = sys.argv[6] if len(sys.argv) > 6 else "z"
MODE = sys.argv[7] if len(sys.argv) > 7 else "projection"
TITLE_TAG = sys.argv[8] if len(sys.argv) > 8 else ""

if MODE not in ("projection", "slice"):
    raise ValueError(f"MODE must be 'projection' or 'slice', got {MODE!r}")

if A_INV != 2.0 or C0 != 0.0:
    raise NotImplementedError(
        "R(tau)/lambda(tau) below are hardcoded to physical mode "
        "(a_inv=2, c0=0); generalise via axion_analysis.Background if a "
        "different mode is needed."
    )

ds = yt.load(PLOTFILE_DIR)
ds.force_periodicity()  # the domain is periodic; needed for gradient ghost zones

a_time = ds.current_time.to_value("code_time")
tau = a_time + TAU_I_POSTHANDOFF
log_mr_over_h = 2.0 * np.log(tau)
print(f"a_time = {a_time}, tau = {tau}, log(m_r/H) = {log_mr_over_h}")

# Physical mode (a_inv=2, c0=0): R(tau) = tau, lambda(tau) = 1 identically
# (Background.R/Background.lam in axion_analysis.py, specialised).
R = tau
lam = 1.0
b_inv = A_INV - 1.0

ds.add_gradient_fields(("boxlib", "psi1"))
ds.add_gradient_fields(("boxlib", "psi2"))


def _rho_tot(field, data):
    """Energy.hpp's rho_tot_pointwise, reimplemented from the raw psi/Pi
    plotfile fields -- not screened/masked (Masking.hpp is not applied
    here; this is the same 'unscreened' pass used elsewhere so string
    cores are visible, not removed)."""
    psi1 = data["boxlib", "psi1"].d
    psi2 = data["boxlib", "psi2"].d
    Pi1 = data["boxlib", "Pi1"].d
    Pi2 = data["boxlib", "Pi2"].d
    c1 = Pi1 - psi1 / (b_inv * tau)
    c2 = Pi2 - psi2 / (b_inv * tau)
    kinetic = c1**2 + c2**2
    grad_psi1_sq = (
        data["boxlib", "psi1_gradient_x"].d ** 2
        + data["boxlib", "psi1_gradient_y"].d ** 2
        + data["boxlib", "psi1_gradient_z"].d ** 2
    )
    grad_psi2_sq = (
        data["boxlib", "psi2_gradient_x"].d ** 2
        + data["boxlib", "psi2_gradient_y"].d ** 2
        + data["boxlib", "psi2_gradient_z"].d ** 2
    )
    gradient = grad_psi1_sq + grad_psi2_sq
    psi_sq_minus_R_sq = psi1**2 + psi2**2 - R**2
    potential = 0.25 * lam * psi_sq_minus_R_sq**2
    return (kinetic + gradient + potential) / R**4


ds.add_field(
    ("gas", "rho_tot"),
    function=_rho_tot,
    sampling_type="local",
    units="",
    take_log=True,
)

common_kwargs = dict(center=ds.domain_center, width=ds.domain_width[0])
if MODE == "projection":
    pw = yt.ProjectionPlot(ds, AXIS, ("gas", "rho_tot"), method="max", **common_kwargs)
    mode_desc = f"max-intensity projection along {AXIS}"
else:
    pw = yt.SlicePlot(ds, AXIS, ("gas", "rho_tot"), **common_kwargs)
    mode_desc = f"{AXIS}-slice"

pw.set_cmap(("gas", "rho_tot"), "inferno")
pw.set_colorbar_label(("gas", "rho_tot"), r"$\rho_{\rm tot}$ (unscreened, code units)")
# periodic=False: annotate_grids()'s default tiles neighbouring periodic
# images of the grid structure for wraparound visualisation -- confusing
# here, since it makes it look like the plotted *data* extends beyond the
# real domain too (it does not; only the box outlines were being tiled).
# On a ProjectionPlot this draws each grid's projected (x,y) footprint,
# i.e. the union over the line of sight -- see module docstring.
pw.annotate_grids(edgecolors="cyan", linewidth=0.6, periodic=False)
tag = f", {TITLE_TAG}" if TITLE_TAG else ""
pw.annotate_title(
    f"rho_tot (unscreened), {mode_desc}, tau={tau:.3f}, "
    f"log(m_r/H)={log_mr_over_h:.3f}{tag}"
)
pw.save(OUT_PNG)
print(f"saved {OUT_PNG}")
