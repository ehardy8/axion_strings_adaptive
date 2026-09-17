# Milestone-1 status

Tracks progress against `milestone-1.md` task by task. Updated as work lands.

| Task | Status | Notes |
|---|---|---|
| 1.1 Skeleton from KleinGordon | **Done** | `AxionStrings/` example, out-of-tree against `GRTeclyn`/`amrex` submodules. Builds and runs the unmodified Wave model end to end (8 steps, plotfiles + checkpoint). See README for the macOS/homebrew-LLVM `PATH` workaround needed to build. |
| 1.2 Complex state | **Done** | `NUM_VARS=4` (`psi1,psi2,Pi1,Pi2`) in `StateVariables.hpp`. |
| 1.3 RHS with c(tau) | **Done** | `Background.hpp` (analytic `R(tau)`, `lambda(tau)` with switch re-anchoring, AMReX-free) + `AxionStringsRHS` (complex EOM, no `model_t` template -- one physics, `c(tau)` only enters via `lambda(tau)`). `R(tau)`/`lambda(tau)`/curvature coeff computed on host in `specific_eval_rhs` and passed in as plain values (no ParmParse in the kernel). KO dissipation omitted entirely (not just `sigma=0`) per CLAUDE.md constraint 1. Verified: all three modes (c=0,1,1+b_inv) run; `m_r/H` from the code (`H_over_mr_direct`, using R/R'/lambda) matches the sec.5 closed form to round-off, in both a standalone doctest suite (`AxionStrings/tests/test_background.cpp`, AMReX-free) and at runtime via a `specific_post_timestep` diagnostic. Qualitative behaviour also matches the sec.4 table: m_r/H grows (c=0), grows more slowly (c=1), constant (c=2=Moore in RD). Placeholder initial data: homogeneous `psi1=R(tau)`, an exact solution of the free EOM (real ICs are task 1.5). |
| 1.4 Parameters and box planning | Not started | |
| 1.5 Initial conditions and restarts | Not started | |
| 1.6 String finder and xi | Not started | |
| 1.7 Masking | Not started | |
| 1.8 Energies | Not started | |
| 1.9 Spectra | Not started | |
| 1.10 Remaining diagnostics | Not started | |
| T1 static straight string | Not started | |
| T2 mask unit test | Not started | |
| T3 plane wave | Not started | |
| T4 collapsing loop | Not started | |
| T5 float vs double | Not started | |

## Known issues / open questions

- `Background::H_over_mr_closed_form` only equals `H_over_mr_direct` in the
  no-switch case (verified analytically and by unit test). Across a c(tau)
  switch, the sec.5 closed form `(tau/tau0)^((c-a_inv)/(a_inv-1))` does not
  hold as-is with the original `tau0` (conventions.md sec.12 flags "a
  separate branch for constant-log runs" for the related `log(m_r/H)`).
  `H_over_mr_direct` (from `R`, `R'`, `lambda` directly) is unaffected and
  always correct; nothing in the current code relies on the closed form
  across a switch. Worth resolving properly before the fat->Moore protocol
  (task 1.4+) needs a closed-form `log(m_r/H)` post-switch.
