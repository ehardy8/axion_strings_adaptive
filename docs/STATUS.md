# Milestone-1 status

Tracks progress against `milestone-1.md` task by task. Updated as work lands.

| Task | Status | Notes |
|---|---|---|
| 1.1 Skeleton from KleinGordon | **Done** | `AxionStrings/` example, out-of-tree against `GRTeclyn`/`amrex` submodules. Builds and runs the unmodified Wave model end to end (8 steps, plotfiles + checkpoint). See README for the macOS/homebrew-LLVM `PATH` workaround needed to build. |
| 1.2 Complex state | **Done** | `NUM_VARS=4` (`psi1,psi2,Pi1,Pi2`) in `StateVariables.hpp`. |
| 1.3 RHS with c(tau) | **Done** | `Background.hpp` (analytic `R(tau)`, `lambda(tau)` with switch re-anchoring, AMReX-free) + `AxionStringsRHS` (complex EOM, no `model_t` template -- one physics, `c(tau)` only enters via `lambda(tau)`). `R(tau)`/`lambda(tau)`/curvature coeff computed on host in `specific_eval_rhs` and passed in as plain values (no ParmParse in the kernel). KO dissipation omitted entirely (not just `sigma=0`) per CLAUDE.md constraint 1. Verified: all three modes (c=0,1,1+b_inv) run; `m_r/H` from the code (`H_over_mr_direct`, using R/R'/lambda) matches the sec.5 closed form to round-off, in both a standalone doctest suite (`AxionStrings/tests/test_background.cpp`, AMReX-free) and at runtime via a `specific_post_timestep` diagnostic. Qualitative behaviour also matches the sec.4 table: m_r/H grows (c=0), grows more slowly (c=1), constant (c=2=Moore in RD). Placeholder initial data: homogeneous `psi1=R(tau)`, an exact solution of the free EOM (real ICs are task 1.5). |
| 1.4 Parameters and box planning | **Done** | `BoxPlan.hpp` (AMReX-free, unit-tested): general box plan (`tau_f`, `L_tilde`, `delta_x` from `N,N1,N2,a_inv,c`, sec.5) plus the Moore-phase dynamic-range check (`log(H0/H)_max` from `N,N1,N2,gamma`). `AxionStringsParams::apply_box_plan` wires this into `SimulationParameters::check_params()` (runs first, before `BaseParameterChecker`, since it may need to inject `geometry.prob_extent`/`amr.n_cell`/`evolution.stop_time` before those are validated elsewhere): if `axion_strings.N` is absent, box planning is skipped entirely (ad hoc/manual configs, e.g. the 1.2/1.3 smoke tests, are unaffected); if present, derived quantities are echoed at startup and injected into ParmParse (so they land in `parameters_and_version.txt` too) when not already user-set, or cross-checked (abort on mismatch) when they are. Any configured switch is treated as the fat->Moore protocol (sec.4's only described use) and routed to the dynamic-range check instead of the general formula, which is singular at `c = a_inv`; a single constant `c0 = a_inv` (Moore run with no switch) is also detected and routed there. Verified by running all of: no-box-planning (unchanged), auto-derived geometry/n_cell/stop_time, a deliberate mismatch (aborts), Moore mode with an achievable and an unachievable target dynamic range (accepts/aborts correctly), and the switch-triggered Moore path. |
| 1.5 Initial conditions and restarts | **Done, with one open architectural gap (see below)** | `XiFormula.hpp`: sec.8 xi<->N_p relation. `PreEvolutionBackground.hpp` (user-confirmed derivation): exponential `R(tau_pre)`, constant curvature coefficient, `c=1`-shaped `lambda`; `R*m_r` and `m_r/H` both exactly constant (unit-tested). `FourierIC.hpp`: Gaussian random field via `amrex::FFT::R2C` (needs `USE_FFT=TRUE` + `FFTW_HOME` in the `GNUmakefile` -- see README), modes occupied for `|k|<=k_max` and zero above (and at `k=0`, since the target is a *variance*), normalised to a target mean-square variance by generate-then-measure-then-rescale (sidesteps needing the FFT library's own normalisation convention). Mode amplitudes come from a hash of `(seed, component, i, j, k)`, not an RNG stream -- reproducible at fixed `(seed, parameters, rank count)` per CLAUDE.md constraint 7 (checked: *not* also independent of rank count -- 1 vs 2 ranks differ, though each is internally reproducible; this was never required). `ic_mode=fourier` (`AxionStringsLevel::initData`) generates directly at `tau_i`, the user's requested skip-pre-evolution option; `axion_strings.mode=pre_evolution` always uses it in the pre-evolution frame. The xi-monitoring stopping loop (`specific_post_timestep`, pre-evolution mode) implements the user's adaptive cadence exactly: coarse checks (measuring xi is a full plaquette count, not cheap) tightening to medium within `check_threshold_medium` (default 1.5x) of target, fine within `check_threshold_fine` (default 1.1x), stopping the first time xi drops to the target (xi generally decreases -- user's observation) via an `okToContinue()` override; the existing post-loop checkpoint write in `Main_AxionStrings.cpp` then fires unchanged. Handoff rescale (`specific_post_restart`, gated on the one-shot `axion_strings.restart_from_pre_evolution` flag, kept separate from an ordinary resume so it isn't reapplied) derived from the `psi=R phi/v` chain rule across the two different time parametrisations; verified **N_p is preserved exactly** across the rescale (164 both immediately before and after handoff in the end-to-end test below), a strong correctness check since a positive real rescale cannot change phase/winding. Full pipeline (pre-evolution run -> adaptive-cadence stop -> checkpoint -> main restart -> rescale -> continued evolution) run end to end and produced sensible diagnostics throughout. |
| 1.6 String finder and xi | **Done** | `PlaquetteWinding.hpp` (dual-purpose: falls back to plain host macros when compiled standalone for its unit test, picks up AMReX's real GPU macros when included from the kernel file -- one definition, not inlined at multiple call sites, per CLAUDE.md constraint 5's spirit) implements the sec.8 winding test, unit-tested including higher-winding and branch-cut-wraparound cases. `StringFinder.hpp`'s `count_plaquettes` reduces this over the 3 plaquettes/cell (xy,yz,zx, low-index corner) via `amrex::ReduceOps`, giving plain and winding-weighted counts. Wired into `specific_post_timestep` (xi via `XiFormula.hpp`, both plain and weighted, alongside the existing m_r/H diagnostic) and also written to a `network_scalars.dat` time series via GRTeclyn's `SmallDataIO`. T1 (plaquette-detection/xi-normalisation slice; the full core-profile/tension T1 needs tasks 1.7/1.8 too) verified with a `straight_string_test` IC mode: **a single isolated vortex cannot be embedded in a periodic box** (net winding over a closed 2-torus must vanish; discovered this the hard way -- a naive single-vortex field gave exactly double the expected count, from spurious extra windings along the periodic seam where the field failed to match itself) -- fixed with a compensating vortex/antivortex pair. Confirmed the resulting N_p exactly equals `2 * N_z` (one winding of each sign per z-layer, zero yz/zx contamination) at two different resolutions/`c`/`tau_i`, and xi matches `XiFormula.hpp`'s formula by hand-calculation. Moore's trick's extra `H0^-2` factor (sec.8) is not yet implemented -- flagged below. |
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

- Box planning in Moore mode (`apply_box_plan`) only runs the sec.5
  dynamic-range check; it does not derive `geometry.prob_extent` or
  `evolution.stop_time` the way the general (non-Moore) path does, since
  milestone-1.md task 1.4 asks specifically for "the Moore-phase
  dynamic-range formula as a startup check", not a closed form for
  `L_tilde`/`tau_f` (sec.5 doesn't give one -- the general formula is
  singular at `c = a_inv`). A production fat->Moore run currently needs
  `geometry.prob_extent`/`evolution.stop_time` set by hand.

- Moore's trick needs an extra `H0^-2` factor in the xi formula (conventions.md
  sec.8, "source document eq. 84") that `XiFormula.hpp`/`StringFinder.hpp`
  don't yet implement -- not needed for the single-level, non-Moore T1
  validation done so far, but must be added before xi is trusted for a
  Moore-mode run.

- **Open architectural gap in the pre-evolution -> main handoff (task 1.5):**
  `AxionStringsLevel::specific_post_restart`'s rescale is verified correct
  (N_p preserved exactly across handoff, in an end-to-end test), but it
  relies on plain `amr.restart`, and **AMReX's restart mechanism does not
  resize the domain** -- tested directly: giving the main run a different
  `geometry.prob_extent` than the pre-evolution checkpoint's is silently
  ignored, and the checkpoint's own (smaller) domain is kept. Conventions.md
  sec.7's `L_tilde_init` formula (`BoxPlan::pre_evolution_L_tilde`) generally
  gives a *different* comoving box size for pre-evolution than for the main
  run (not necessarily smaller -- e.g. for physical/c=0 with `tau_i > tau0`,
  `L_init = L_main * tau_i`, i.e. *larger*), so this will bite in practice,
  not just in principle. A real fix needs a custom procedure that reads the
  pre-evolution checkpoint's raw field data and places it onto the main
  run's own (differently sized) grid -- e.g. periodic tiling/cropping --
  rather than a plain AMReX restart. Not yet designed or implemented;
  everything tested so far used matching box sizes for both executions.

- Pre-evolution's own diagnostics reuse the general `network_scalars.dat`/
  `specific_post_timestep` machinery only partially -- the xi-monitoring
  loop prints to stdout but does not (yet) also write a `network_scalars.dat`-
  style time series for the pre-evolution phase. Low priority (pre-evolution
  is a short-lived warm-up execution, not itself a headline result) but
  worth adding for completeness/debugging.

- `Main_AxionStrings.cpp`'s time-stepping loop compares `amr.levelSteps(0)`
  (an absolute step count, carried over from a checkpoint on restart)
  against `evolution.max_steps`. Discovered while testing the handoff: a
  restart needs `max_steps` set to (checkpoint's step count) + (desired
  further steps), not just "how many more steps" -- inherited from the
  original KleinGordon example, not something task 1.5 changed, but worth
  knowing when configuring a real restart.

- Pre-evolution's EOM (conventions.md sec.7) is not an instance of the main
  `Background` class -- `R = R0(t/t0)` (linear in *cosmic* time) is a
  genuinely distinct, singular case (`a_inv=1`, where `b_inv=0` and
  `R0=1/b_inv` diverges in the main formulas), not reducible to fat mode
  under radiation domination (checked: fat/RD gives `m_r/H ~ tau`, not
  constant). Derived independently (not yet reviewed/confirmed by the user):
  switching to pre-evolution's own conformal time `tau_pre` (`dtau_pre =
  dt/R(t)`) turns `R(t)=R0(t/t0)` into `R(tau_pre) = R0 e^(alpha*tau_pre)`,
  `alpha = R0/t0`. The same `psi=R phi/v` rescaling trick from sec.3 (generic
  in `R(tau)`, not tied to the power-law form) then gives an EOM with the
  *same* Laplacian/curvature/potential structure as the main RHS, just with
  a constant curvature coefficient `alpha^2` (instead of
  `(1-b_inv)/(b_inv^2 tau^2)`) and `lambda(tau_pre) = lambda_pre0
  (R(tau_pre)/R0)^-2` (the `c=1`-shaped formula, with the new `R`). Verified
  `R*m_r` and `m_r/H` are both exactly constant under this, matching sec.7's
  stated requirement. Since `R0`/`alpha` are pre-evolution's own time-gauge
  freedom, `R0=alpha=1` WLOG, leaving one physical input: `gamma_pre`, the
  constant `m_r/H` to relax the network at. This means `AxionStringsRHS`
  can be reused unchanged for pre-evolution -- only a new
  `PreEvolutionBackground` (same 3-quantity interface as `Background`) is
  needed. Not yet implemented pending confirmation this derivation is right.
