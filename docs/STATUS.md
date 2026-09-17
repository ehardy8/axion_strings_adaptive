# Milestone-1 status

Tracks progress against `milestone-1.md` task by task. Updated as work lands.

| Task | Status | Notes |
|---|---|---|
| 1.1 Skeleton from KleinGordon | **Done** | `AxionStrings/` example, out-of-tree against `GRTeclyn`/`amrex` submodules. Builds and runs the unmodified Wave model end to end (8 steps, plotfiles + checkpoint). See README for the macOS/homebrew-LLVM `PATH` workaround needed to build. |
| 1.2 Complex state | **Done** | `NUM_VARS=4` (`psi1,psi2,Pi1,Pi2`) in `StateVariables.hpp`. |
| 1.3 RHS with c(tau) | **Done** | `Background.hpp` (analytic `R(tau)`, `lambda(tau)` with switch re-anchoring, AMReX-free) + `AxionStringsRHS` (complex EOM, no `model_t` template -- one physics, `c(tau)` only enters via `lambda(tau)`). `R(tau)`/`lambda(tau)`/curvature coeff computed on host in `specific_eval_rhs` and passed in as plain values (no ParmParse in the kernel). KO dissipation omitted entirely (not just `sigma=0`) per CLAUDE.md constraint 1. Verified: all three modes (c=0,1,1+b_inv) run; `m_r/H` from the code (`H_over_mr_direct`, using R/R'/lambda) matches the sec.5 closed form to round-off, in both a standalone doctest suite (`AxionStrings/tests/test_background.cpp`, AMReX-free) and at runtime via a `specific_post_timestep` diagnostic. Qualitative behaviour also matches the sec.4 table: m_r/H grows (c=0), grows more slowly (c=1), constant (c=2=Moore in RD). Placeholder initial data: homogeneous `psi1=R(tau)`, an exact solution of the free EOM (real ICs are task 1.5). |
| 1.4 Parameters and box planning | **Done** | `BoxPlan.hpp` (AMReX-free, unit-tested): general box plan (`tau_f`, `L_tilde`, `delta_x` from `N,N1,N2,a_inv,c`, sec.5) plus the Moore-phase dynamic-range check (`log(H0/H)_max` from `N,N1,N2,gamma`). `AxionStringsParams::apply_box_plan` wires this into `SimulationParameters::check_params()` (runs first, before `BaseParameterChecker`, since it may need to inject `geometry.prob_extent`/`amr.n_cell`/`evolution.stop_time` before those are validated elsewhere): if `axion_strings.N` is absent, box planning is skipped entirely (ad hoc/manual configs, e.g. the 1.2/1.3 smoke tests, are unaffected); if present, derived quantities are echoed at startup and injected into ParmParse (so they land in `parameters_and_version.txt` too) when not already user-set, or cross-checked (abort on mismatch) when they are. Any configured switch is treated as the fat->Moore protocol (sec.4's only described use) and routed to the dynamic-range check instead of the general formula, which is singular at `c = a_inv`; a single constant `c0 = a_inv` (Moore run with no switch) is also detected and routed there. Verified by running all of: no-box-planning (unchanged), auto-derived geometry/n_cell/stop_time, a deliberate mismatch (aborts), Moore mode with an achievable and an unachievable target dynamic range (accepts/aborts correctly), and the switch-triggered Moore path. |
| 1.5 Initial conditions and restarts | **Done** | `XiFormula.hpp`: sec.8 xi<->N_p relation. `PreEvolutionBackground.hpp` (user-confirmed derivation): exponential `R(tau_pre)`, constant curvature coefficient, `c=1`-shaped `lambda`; `R*m_r` and `m_r/H` both exactly constant (unit-tested). `FourierIC.hpp`: Gaussian random field via `amrex::FFT::R2C` (needs `USE_FFT=TRUE` + `FFTW_HOME` in the `GNUmakefile` -- see README), modes occupied for `|k|<=k_max` and zero above (and at `k=0`, since the target is a *variance*), normalised to a target mean-square variance by generate-then-measure-then-rescale (sidesteps needing the FFT library's own normalisation convention). Mode amplitudes come from a hash of `(seed, component, i, j, k)`, not an RNG stream -- reproducible at fixed `(seed, parameters, rank count)` per CLAUDE.md constraint 7 (checked: *not* also independent of rank count -- 1 vs 2 ranks differ, though each is internally reproducible; this was never required). `ic_mode=fourier` (`AxionStringsLevel::initData`) generates directly at `tau_i`, the user's requested skip-pre-evolution option; `axion_strings.mode=pre_evolution` always uses it in the pre-evolution frame. The xi-monitoring stopping loop (`specific_post_timestep`, pre-evolution mode) implements the user's adaptive cadence exactly: coarse checks (measuring xi is a full plaquette count, not cheap) tightening to medium within `check_threshold_medium` (default 1.5x) of target, fine within `check_threshold_fine` (default 1.1x), stopping the first time xi drops to the target (xi generally decreases -- user's observation) via an `okToContinue()` override; the existing post-loop checkpoint write in `Main_AxionStrings.cpp` then fires unchanged. Handoff rescale (`specific_post_restart`, gated on the one-shot `axion_strings.restart_from_pre_evolution` flag, kept separate from an ordinary resume so it isn't reapplied) derived from the `psi=R phi/v` chain rule across the two different time parametrisations; verified **N_p is preserved exactly** across the rescale (164 both immediately before and after handoff in the end-to-end test below), a strong correctness check since a positive real rescale cannot change phase/winding. Full pipeline (pre-evolution run -> adaptive-cadence stop -> checkpoint -> main restart -> rescale -> continued evolution) run end to end, including via `axion_strings.N/N1/N2` box planning in both parameter files, and produced sensible diagnostics throughout. Pre-evolution deliberately uses the same box as the main run rather than sec.7's `L_tilde_init` formula -- see the decision recorded below. |
| 1.6 String finder and xi | **Done** | `PlaquetteWinding.hpp` (dual-purpose: falls back to plain host macros when compiled standalone for its unit test, picks up AMReX's real GPU macros when included from the kernel file -- one definition, not inlined at multiple call sites, per CLAUDE.md constraint 5's spirit) implements the sec.8 winding test, unit-tested including higher-winding and branch-cut-wraparound cases. `StringFinder.hpp`'s `count_plaquettes` reduces this over the 3 plaquettes/cell (xy,yz,zx, low-index corner) via `amrex::ReduceOps`, giving plain and winding-weighted counts. Wired into `specific_post_timestep` (xi via `XiFormula.hpp`, both plain and weighted, alongside the existing m_r/H diagnostic) and also written to a `network_scalars.dat` time series via GRTeclyn's `SmallDataIO`. T1 (plaquette-detection/xi-normalisation slice; the full core-profile/tension T1 needs tasks 1.7/1.8 too) verified with a `straight_string_test` IC mode: **a single isolated vortex cannot be embedded in a periodic box** (net winding over a closed 2-torus must vanish; discovered this the hard way -- a naive single-vortex field gave exactly double the expected count, from spurious extra windings along the periodic seam where the field failed to match itself) -- fixed with a compensating vortex/antivortex pair. Confirmed the resulting N_p exactly equals `2 * N_z` (one winding of each sign per z-layer, zero yz/zx contamination) at two different resolutions/`c`/`tau_i`, and xi matches `XiFormula.hpp`'s formula by hand-calculation. Moore's trick's extra `H0^-2` factor (sec.8) is not yet implemented -- flagged below. |
| 1.7 Masking | **Done** | `Masking.hpp` (dual-purpose like `PlaquetteWinding.hpp`: AMReX-free for its unit test, picks up real GPU macros when included from a kernel file -- one definition, applied at exactly one place, per CLAUDE.md constraint 5). Three modes behind one interface (`masked_a_dot`, `masking_weight`): `None` (f=1, no top-hat -- the unscreened pass task 1.8 needs), `A` (bare numerator `psi1 Pi2 - Pi1 psi2`, no division, no top-hat), `B` (the same divided by `|psi|^2`, zeroed below a **runtime** threshold on `|psi|/R` -- CLAUDE.md constraint 4; `AxionStringsParams::read_masking_params` reads it, not yet called from anywhere since no energy/spectrum pipeline exists yet to call it). `masking_weight` gives the effective-point-count factor sec.10 requires for unbiased averaging (1 for None/A, the 0/1 top-hat for B). T2 verified as two doctest suites (`tests/test_masking.cpp`), deliberately not blind to "a static string does not test the mask" (sec.13): (1) direct array-level test across the threshold boundary, confirming scheme B is exactly zero below it and bit-exact-equal to the unmasked formula at/above it, plus `None`/`A` behave as specified even deep inside a core; (2) an analytic boosted-string configuration (`psi = R tanh(m rho) e^{i theta}` translated at constant velocity, `Pi` from the exact chain rule, no numerical differentiation) confirming the masked-out region tracks the *moving* core rather than a fixed location, and that the unmasked phase rate is genuinely non-zero away from it (not just trivially zero as in the static case) -- caught two real bugs while building this: a `0/0` NaN at the exact core before the top-hat check short-circuited it, and a first draft of the boosted-string test that sampled only points on the string's own symmetry axis, where the phase-rotation term vanishes trivially regardless of masking. |
| 1.8 Energies | In progress -- normalisation fixed, T1 tension not yet clean | `Energy.hpp`/`EnergyKernel.hpp` implement `rho_tot` and the axion kinetic energy, both screened/unscreened via the sec.7 masking module (`MaskingScheme::None` gives the required unscreened pass "for free"). Wired into `specific_post_timestep`/`network_scalars.dat`. **`rho_tot`'s normalisation was wrong and has been corrected** (2026-09-17, with the user): conventions.md sec.10's formula, read literally, has the kinetic/gradient terms scaling as `R^-2`; re-deriving from scratch (confirmed with the user: a canonically-normalised *complex* scalar has no 1/2 on kinetic/gradient, checked by matching `L=A\|phi_dot\|^2-...` against sec.3's given EOM, which forces `A=1`) gives `R^-4` instead -- verified three independent ways (direct algebra, the derived formula, and reparametrising via `t(tau)` directly). The potential term is unaffected (matches either way once fully expanded). Full derivation is in `Energy.hpp`'s header comment. Separately, found and fixed a real bug in the T1 test IC that the user correctly flagged as "quite a bad sign": `Pi` was set to `0` everywhere, so the far field's `\|psi\|/R` decayed as `R_i/R(tau)` regardless of any string, dropping below the masking threshold everywhere within 2-3 steps and collapsing the screened energy to exactly 0 (confirmed numerically: `R_i/R(1.3)=0.769<0.8`, matching exactly when the collapse occurred) -- fixed by setting `Pi=(R'/R)*psi`, matching how the homogeneous vacuum solution tracks the expanding background. Also swapped the radial profile from `tanh` to `rho_hat/sqrt(rho_hat^2+2)`, which has the same near-core linear behaviour but the *correct* `1-rho_hat^-2` far-field power-law falloff sec.13 specifies (tanh decays exponentially, missing essentially all of the long-range Goldstone tail responsible for the string's logarithmic tension divergence in the first place). Confirmed via a control (the exact homogeneous solution stays at exactly zero energy indefinitely under the same code) that neither bug was a deeper problem with the evolution/energy pipeline. After both fixes: screened energy no longer collapses, and tension values are the right order of magnitude, but the L-dependence still is not cleanly logarithmic (L=8 remains negative; the implied slope between L=16/32/64 is not constant) -- the analytic 2-string profile is still evidently not close enough to the true equilibrium for a quantitatively clean measurement. Not yet resolved: likely needs either a numerically-relaxed (gradient-flow) profile or a more careful measurement protocol; paused to check in with the user before pursuing further. |
| 1.9 Spectra | Not started | |
| 1.10 Remaining diagnostics | Not started | |
| T1 static straight string | Not started | |
| T2 mask unit test | Not started | |
| T3 plane wave | Not started | |
| T4 collapsing loop | Not started | |
| T5 float vs double | Not started | |

## Known issues / open questions

- **Resolved (2026-09-17, with the user): `rho_tot`'s normalisation
  (task 1.8).** Conventions.md sec.10's formula, read literally, has the
  kinetic/gradient terms scaling as `R^-2`. Confirmed with the user that a
  canonically-normalised *complex* scalar has no 1/2 on kinetic/gradient
  terms (checked directly: matching `L = A\|phi_dot\|^2 - A\|grad phi\|^2/R^2
  - V(phi)` against sec.3's given EOM via `d/dt(dL/d(phi_dot)*) = dV/dphi*`
  forces `A=1` exactly). Re-deriving `rho_tot` from this via `phi = v*psi/R`
  and `Pi = d(psi)/d(tau)` gives kinetic/gradient terms scaling as `R^-4`,
  not `R^-2` -- verified three independent ways (direct algebra, the derived
  formula, and reparametrising via `t(tau)` directly). The potential term is
  unaffected (matches conventions.md's literal formula exactly, once fully
  expanded, either way). `Energy.hpp` now implements the re-derived (`R^-4`)
  formula; its header comment has the full derivation. Not yet revisited:
  whether conventions.md's own `rho_tot` (the `R^-2` form) is deliberately a
  *different*, redshift-tracking comoving quantity for some other purpose
  (sec.11's "total box energy should redshift correctly"), distinct from the
  instantaneous physical energy density `Energy.hpp` now computes -- if a
  future diagnostic needs that comoving-tracked quantity specifically, this
  will need a second look.

- **T1 tension test: normalisation and a real IC bug both fixed, but the
  measured tension's L-dependence is still not cleanly logarithmic.**
  Two bugs found and fixed here (2026-09-17): (1) the `rho_tot` normalisation
  above; (2) the `straight_string_test` IC set `Pi=0` everywhere, so the far
  field's `\|psi\|/R` decayed as `R_i/R(tau)` regardless of any string,
  dropping below the masking threshold *everywhere* within 2-3 steps and
  collapsing the screened energy to exactly 0 -- confirmed numerically
  (`R_i/R(1.3) = 0.769 < 0.8`, exactly matching when the collapse occurred)
  and confirmed via a control (the exact homogeneous solution stays at
  exactly zero energy indefinitely under the same evolution code) that this
  was the IC's bug, not a deeper problem. Fixed by setting
  `Pi = (R'/R) * psi`, matching how the homogeneous vacuum solution tracks
  the expanding background. Also swapped the radial profile from `tanh` to
  `rho_hat/sqrt(rho_hat^2+2)` (same near-core linear behaviour, but the
  correct `1-rho_hat^-2` far-field falloff sec.13 specifies -- tanh decays
  exponentially, missing the long-range Goldstone tail the logarithmic
  tension divergence actually comes from). After both fixes, screened energy
  no longer collapses and tension values are the right order of magnitude,
  but: L=8 is still negative, and the implied slope between L=16/32/64 is
  not constant (not cleanly `~ln(L)`). The two-vortex analytic profile
  (product of two single-string ansatze) is evidently still not close enough
  to the true equilibrium for a quantitatively clean measurement -- likely
  needs either a numerically-relaxed (gradient-flow) profile or a more
  careful measurement protocol. Paused here to check in with the user before
  pursuing further.

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

- **Decision (2026-09-17, discussed with the user): pre-evolution uses the
  same `geometry.prob_extent`/`amr.n_cell` as the main run, not
  conventions.md sec.7's `L_tilde_init` formula.** Background: AMReX's
  restart mechanism does not resize the domain (verified directly -- giving
  the main run a different `geometry.prob_extent` than the checkpoint's is
  silently ignored), so a differently-sized pre-evolution box could not
  actually be handed off via `amr.restart` without a custom regridding step.
  Resolved by recognising that `L_tilde_init` only optimises pre-evolution's
  own resolution/box-size trade-off -- it is computed from `a_inv, c, tau_i,
  L_tilde_main` alone, *without* reference to pre-evolution's actual (a
  priori unknown) stopping time, so it was never doing exact physical-size
  matching in the first place; `specific_post_restart`'s rescale (`kappa`,
  from the `psi=R phi/v` chain rule) already corrects the normalisation
  exactly, whatever box size was used, at whatever time pre-evolution
  actually stops. Using `L_tilde_main` throughout costs nothing beyond a
  possibly-suboptimal (not incorrect) pre-evolution resolution, and both
  executions get the same box "for free" by construction when both read the
  same `axion_strings.N/N1/N2/a_inv/c0/tau_i`.
  `BoxPlan::pre_evolution_L_tilde` is still computed and printed (FYI only,
  explicitly marked "NOT used") in pre-evolution mode, for reference. If a
  future run wants pre-evolution's resolution deliberately optimised via the
  sec.7 formula after all, that still needs the custom read-checkpoint-onto-
  a-new-grid procedure described in the previous version of this note --
  not implemented.
  Verified end to end using `axion_strings.N/N1/N2` box planning (not just
  hand-matched `geometry.prob_extent`) in both the pre-evolution and main
  parameter files: handoff succeeds, `kappa` comes out sane (~0.997 for the
  small test case used).

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
