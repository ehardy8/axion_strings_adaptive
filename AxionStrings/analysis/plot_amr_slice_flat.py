#!/usr/bin/env python3
"""Flat-space variant of plot_amr_slice.py (2026-09-19), for the
four-string-collision N_p-dropout investigation: plot_amr_slice.py
hardcodes the FRW physical-mode background (R(tau)=tau, lam=1) and
explicitly refuses any other mode. Flat space (axion_strings.flat_space=1)
has R=1, lambda=m_r^2, R_prime_over_R=0 identically -- see
FlatBackground.hpp -- so rho_tot_pointwise's kinetic term is just
Pi1^2+Pi2^2 with no b_inv/tau subtraction, and there is no log(m_r/H) to
report (just a_time itself).

Usage:
    python3 plot_amr_slice_flat.py PLOTFILE_DIR OUT_PNG M_R [AXIS] [MODE] [TITLE_TAG]
"""
import sys

import numpy as np
import yt

PLOTFILE_DIR = sys.argv[1]
OUT_PNG = sys.argv[2]
M_R = float(sys.argv[3])
AXIS = sys.argv[4] if len(sys.argv) > 4 else "z"
MODE = sys.argv[5] if len(sys.argv) > 5 else "projection"
TITLE_TAG = sys.argv[6] if len(sys.argv) > 6 else ""

if MODE not in ("projection", "slice"):
    raise ValueError(f"MODE must be 'projection' or 'slice', got {MODE!r}")

ds = yt.load(PLOTFILE_DIR)
ds.force_periodicity()

a_time = ds.current_time.to_value("code_time")
print(f"a_time = {a_time}")

lam = M_R * M_R

ds.add_gradient_fields(("boxlib", "psi1"))
ds.add_gradient_fields(("boxlib", "psi2"))


def _rho_tot(field, data):
    """Energy.hpp's rho_tot_pointwise at R=1, R_prime_over_R=0 (flat
    space): kinetic is just Pi1^2+Pi2^2, no R^4 rescale."""
    psi1 = data["boxlib", "psi1"].d
    psi2 = data["boxlib", "psi2"].d
    Pi1 = data["boxlib", "Pi1"].d
    Pi2 = data["boxlib", "Pi2"].d
    kinetic = Pi1**2 + Pi2**2
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
    psi_sq_minus_1 = psi1**2 + psi2**2 - 1.0
    potential = 0.25 * lam * psi_sq_minus_1**2
    return kinetic + gradient + potential


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
pw.annotate_grids(edgecolors="cyan", linewidth=0.6, periodic=False)
tag = f", {TITLE_TAG}" if TITLE_TAG else ""
pw.annotate_title(f"rho_tot (unscreened), {mode_desc}, t={a_time:.3f}{tag}")
pw.save(OUT_PNG)
print(f"saved {OUT_PNG}")
