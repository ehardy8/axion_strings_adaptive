# Analysis scripts

Reusable, tested Python utilities for turning `network_scalars.dat`/
`axion_spectrum.dat`/`axion_projection.dat` into plots. Consolidates the
`Background(tau)` reimplementation, file loading, and unit conversions that
had been copy-pasted into one-off scratchpad scripts (2026-09-19) -- new
plots should build on `axion_analysis.py` rather than re-deriving these
from scratch again.

## Layout

- `axion_analysis.py` -- `Background` (mirrors `Background.hpp`), file
  loaders with named columns (not positional indices), the `v3_drho_dk`/
  `k_over_H` spectrum unit conversions (see `docs/STATUS.md` for the
  2026-09-19 bugfix these encode), and `instantaneous_emission` (the
  F(k/H, m_r/H) extraction of Fleury & Moore 1806.04677 eq.33).
- `test_axion_analysis.py` -- run with `python3 test_axion_analysis.py`
  (needs `scipy`, only for the tests' independent cross-check integrals).
  No pytest dependency, matching the C++ side's lightweight doctest style.
  Re-run this after touching `axion_analysis.py` -- the spectrum unit
  conversion in particular has a documented history of looking-plausible-
  but-wrong (see `docs/STATUS.md` and `[[unverified_formula_was_actually_wrong]]`
  in the session memory).
- `plot_instantaneous_emission.py` -- CLI driver for the F(k/H) plot.
  `python3 plot_instantaneous_emission.py NETWORK_FILE SPECTRUM_FILE OUT_DIR
  A_INV C0 N [DELTA_LOG] [TITLE_TAG]`.
- `plot_run_summary.py` -- CLI driver for the standard per-run plot set:
  `xi`, energies (total/axion/radial, raw and `H^2 f_a^2`-normalised),
  string tension, and the axion spectrum vs `k/H` (raw and
  `H f_a^2`-normalised). `python3 plot_run_summary.py NETWORK_FILE
  SPECTRUM_FILE OUT_DIR A_INV C0 N [TITLE_TAG]`. Supersedes the earlier
  one-off scratchpad version of this script (2026-09-18/19) -- verified
  pixel-identical output on both the `N=256` fat and `N=384` physical runs
  before that script was retired.

## Adding a new plot

Import from `axion_analysis` rather than reimplementing `Background`,
`v3_drho_dk`, or `k_over_H` locally -- every prior duplication of these
(there have been several) has eventually diverged or had a units bug. If a
new derived quantity is genuinely reusable (not a one-off), add it to
`axion_analysis.py` with a docstring and a test, following
`instantaneous_emission`'s pattern: derive it, then verify it against an
*independent* construction (a different quadrature, a different code path,
or an already-validated quantity from `network_scalars.dat`) -- not just a
restatement of the same formula, which would validate arithmetic but not
physics.
