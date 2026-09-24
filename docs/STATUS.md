# Milestone-1 status

Tracks progress against `milestone-1.md` task by task. Updated as work lands.

| Task | Status | Notes |
|---|---|---|
| 1.1 Skeleton from KleinGordon | **Done** | `AxionStrings/` example, out-of-tree against `GRTeclyn`/`amrex` submodules. Builds and runs the unmodified Wave model end to end (8 steps, plotfiles + checkpoint). See README for the macOS/homebrew-LLVM `PATH` workaround needed to build. |
| 1.2 Complex state | **Done** | `NUM_VARS=4` (`psi1,psi2,Pi1,Pi2`) in `StateVariables.hpp`. |
| 1.3 RHS with c(tau) | **Done** | `Background.hpp` (analytic `R(tau)`, `lambda(tau)` with switch re-anchoring, AMReX-free) + `AxionStringsRHS` (complex EOM, no `model_t` template -- one physics, `c(tau)` only enters via `lambda(tau)`). `R(tau)`/`lambda(tau)`/curvature coeff computed on host in `specific_eval_rhs` and passed in as plain values (no ParmParse in the kernel). KO dissipation omitted entirely (not just `sigma=0`) per CLAUDE.md constraint 1. Verified: all three modes (c=0,1,1+b_inv) run; `m_r/H` from the code (`H_over_mr_direct`, using R/R'/lambda) matches the sec.5 closed form to round-off, in both a standalone doctest suite (`AxionStrings/tests/test_background.cpp`, AMReX-free) and at runtime via a `specific_post_timestep` diagnostic. Qualitative behaviour also matches the sec.4 table: m_r/H grows (c=0), grows more slowly (c=1), constant (c=2=Moore in RD). Placeholder initial data: homogeneous `psi1=R(tau)`, an exact solution of the free EOM (real ICs are task 1.5). |
| 1.4 Parameters and box planning | **Done** | `BoxPlan.hpp` (AMReX-free, unit-tested): general box plan (`tau_f`, `L_tilde`, `delta_x` from `N,N1,N2,a_inv,c`, sec.5) plus the Moore-phase dynamic-range check (`log(H0/H)_max` from `N,N1,N2,gamma`). `AxionStringsParams::apply_box_plan` wires this into `SimulationParameters::check_params()` (runs first, before `BaseParameterChecker`, since it may need to inject `geometry.prob_extent`/`amr.n_cell`/`evolution.stop_time` before those are validated elsewhere): if `axion_strings.N` is absent, box planning is skipped entirely (ad hoc/manual configs, e.g. the 1.2/1.3 smoke tests, are unaffected); if present, derived quantities are echoed at startup and injected into ParmParse (so they land in `parameters_and_version.txt` too) when not already user-set, or cross-checked (abort on mismatch) when they are. Any configured switch is treated as the fat->Moore protocol (sec.4's only described use) and routed to the dynamic-range check instead of the general formula, which is singular at `c = a_inv`; a single constant `c0 = a_inv` (Moore run with no switch) is also detected and routed there. Verified by running all of: no-box-planning (unchanged), auto-derived geometry/n_cell/stop_time, a deliberate mismatch (aborts), Moore mode with an achievable and an unachievable target dynamic range (accepts/aborts correctly), and the switch-triggered Moore path. |
| 1.5 Initial conditions and restarts | **Done** | `XiFormula.hpp`: sec.8 xi<->N_p relation. `PreEvolutionBackground.hpp` (user-confirmed derivation): exponential `R(tau_pre)`, constant curvature coefficient, `c=1`-shaped `lambda`; `R*m_r` and `m_r/H` both exactly constant (unit-tested). `FourierIC.hpp`: Gaussian random field via `amrex::FFT::R2C` (needs `USE_FFT=TRUE` + `FFTW_HOME` in the `GNUmakefile` -- see README), modes occupied for `|k|<=k_max` and zero above (and at `k=0`, since the target is a *variance*), normalised to a target mean-square variance by generate-then-measure-then-rescale (sidesteps needing the FFT library's own normalisation convention). Mode amplitudes come from a hash of `(seed, component, i, j, k)`, not an RNG stream -- reproducible at fixed `(seed, parameters, rank count)` per CLAUDE.md constraint 7 (checked: *not* also independent of rank count -- 1 vs 2 ranks differ, though each is internally reproducible; this was never required). `ic_mode=fourier` (`AxionStringsLevel::initData`) generates directly at `tau_i`, the user's requested skip-pre-evolution option; `axion_strings.mode=pre_evolution` always uses it in the pre-evolution frame. The xi-monitoring stopping loop (`specific_post_timestep`, pre-evolution mode) implements the user's adaptive cadence exactly: coarse checks (measuring xi is a full plaquette count, not cheap) tightening to medium within `check_threshold_medium` (default 1.5x) of target, fine within `check_threshold_fine` (default 1.1x), stopping the first time xi drops to the target (xi generally decreases -- user's observation) via an `okToContinue()` override; the existing post-loop checkpoint write in `Main_AxionStrings.cpp` then fires unchanged. Handoff rescale (`specific_post_restart`, gated on the one-shot `axion_strings.restart_from_pre_evolution` flag, kept separate from an ordinary resume so it isn't reapplied) derived from the `psi=R phi/v` chain rule across the two different time parametrisations; verified **N_p is preserved exactly** across the rescale (164 both immediately before and after handoff in the end-to-end test below), a strong correctness check since a positive real rescale cannot change phase/winding. Full pipeline (pre-evolution run -> adaptive-cadence stop -> checkpoint -> main restart -> rescale -> continued evolution) run end to end, including via `axion_strings.N/N1/N2` box planning in both parameter files, and produced sensible diagnostics throughout. Pre-evolution deliberately uses the same box as the main run rather than sec.7's `L_tilde_init` formula -- see the decision recorded below. **Superseded (2026-09-18) -- see the "single continuous run" entry further down**: `axion_strings.mode`/`Mode::PreEvolution`, `restart_from_pre_evolution` and `specific_post_restart` described here have all been removed; the mechanism (xi-monitoring cadence, handoff rescale math, "same box throughout" decision) is unchanged, just applied in place within one run instead of via an external checkpoint/restart. |
| 1.6 String finder and xi | **Done** | `PlaquetteWinding.hpp` (dual-purpose: falls back to plain host macros when compiled standalone for its unit test, picks up AMReX's real GPU macros when included from the kernel file -- one definition, not inlined at multiple call sites, per CLAUDE.md constraint 5's spirit) implements the sec.8 winding test, unit-tested including higher-winding and branch-cut-wraparound cases. `StringFinder.hpp`'s `count_plaquettes` reduces this over the 3 plaquettes/cell (xy,yz,zx, low-index corner) via `amrex::ReduceOps`, giving plain and winding-weighted counts. Wired into `specific_post_timestep` (xi via `XiFormula.hpp`, both plain and weighted, alongside the existing m_r/H diagnostic) and also written to a `network_scalars.dat` time series via GRTeclyn's `SmallDataIO`. T1 (plaquette-detection/xi-normalisation slice; the full core-profile/tension T1 needs tasks 1.7/1.8 too) verified with a `straight_string_test` IC mode: **a single isolated vortex cannot be embedded in a periodic box** (net winding over a closed 2-torus must vanish; discovered this the hard way -- a naive single-vortex field gave exactly double the expected count, from spurious extra windings along the periodic seam where the field failed to match itself) -- fixed with a compensating vortex/antivortex pair. Confirmed the resulting N_p exactly equals `2 * N_z` (one winding of each sign per z-layer, zero yz/zx contamination) at two different resolutions/`c`/`tau_i`, and xi matches `XiFormula.hpp`'s formula by hand-calculation. Moore's trick's extra `H0^-2` factor (sec.8) is not yet implemented -- flagged below. |
| 1.7 Masking | **Done** | `Masking.hpp` (dual-purpose like `PlaquetteWinding.hpp`: AMReX-free for its unit test, picks up real GPU macros when included from a kernel file -- one definition, applied at exactly one place, per CLAUDE.md constraint 5). Three modes behind one interface (`masked_a_dot`, `masking_weight`): `None` (f=1, no top-hat -- the unscreened pass task 1.8 needs), `A` (bare numerator `psi1 Pi2 - Pi1 psi2`, no division, no top-hat), `B` (the same divided by `|psi|^2`, zeroed below a **runtime** threshold on `|psi|/R` -- CLAUDE.md constraint 4; `AxionStringsParams::read_masking_params` reads it, not yet called from anywhere since no energy/spectrum pipeline exists yet to call it). `masking_weight` gives the effective-point-count factor sec.10 requires for unbiased averaging (1 for None/A, the 0/1 top-hat for B). T2 verified as two doctest suites (`tests/test_masking.cpp`), deliberately not blind to "a static string does not test the mask" (sec.13): (1) direct array-level test across the threshold boundary, confirming scheme B is exactly zero below it and bit-exact-equal to the unmasked formula at/above it, plus `None`/`A` behave as specified even deep inside a core; (2) an analytic boosted-string configuration (`psi = R tanh(m rho) e^{i theta}` translated at constant velocity, `Pi` from the exact chain rule, no numerical differentiation) confirming the masked-out region tracks the *moving* core rather than a fixed location, and that the unmasked phase rate is genuinely non-zero away from it (not just trivially zero as in the static case) -- caught two real bugs while building this: a `0/0` NaN at the exact core before the top-hat check short-circuited it, and a first draft of the boosted-string test that sampled only points on the string's own symmetry axis, where the phase-rotation term vanishes trivially regardless of masking. |
| 1.8 Energies | Mostly done -- core-energy extraction validated; the log(L) tail growth is a known, deferred rough edge | `Energy.hpp`/`EnergyKernel.hpp` implement `rho_tot`, axion kinetic and axion gradient energy, each screened/unscreened via the sec.7 masking module (`MaskingScheme::None` gives the unscreened pass "for free"). **`rho_tot`'s normalisation was wrong and has been corrected** (2026-09-17, with the user): conventions.md sec.10's formula, read literally, has kinetic/gradient terms scaling as `R^-2`; re-deriving from scratch (a canonically-normalised *complex* scalar has no 1/2 on kinetic/gradient -- checked by matching `L=A\|phi_dot\|^2-...` against sec.3's given EOM, which forces `A=1`) gives `R^-4` instead, verified three independent ways. The potential term is unaffected (matches either way once expanded). Full derivation in `Energy.hpp`'s header comment. Found and fixed a real T1 IC bug the user correctly flagged as "quite a bad sign": `Pi=0` everywhere let the far field's `\|psi\|/R` decay as `R_i/R(tau)` regardless of any string, collapsing the screened energy to 0 within 2-3 steps (confirmed numerically: `R_i/R(1.3)=0.769<0.8`, exactly matching when it broke) -- fixed via `Pi=(R'/R)*psi`, matching the homogeneous solution's background-tracking. Also swapped the radial profile from `tanh` to `rho_hat/sqrt(rho_hat^2+2)` for the correct `1-rho_hat^-2` far-field falloff (tanh decays exponentially, missing the long-range Goldstone tail responsible for the log tension divergence). A control (exact homogeneous solution stays at exactly zero energy indefinitely) confirmed neither bug was a deeper pipeline problem. **Also reworked how tension is extracted** (2026-09-17, with the user): the naive screened/unscreened *difference* of averages is diluted by the masked point-count fraction and is not itself "the string energy" -- switched to explicitly summing the core contribution (`n_total*unscreened - n_unmasked*screened`), which **converges cleanly to a constant (~3.72) independent of L for L>=32** (checked L=16..256 at `N_2=1`, per the user's guidance that `m_r*dx=1` is their production convention, not the far-finer resolutions tried initially) -- exactly the expected behaviour for energy localised at the core. All raw components (screened+unscreened rho_tot/axion_kinetic/axion_gradient, plus `n_total`/`n_unmasked`) are saved to `network_scalars.dat` every snapshot rather than baking in one derived formula, per the user: any combination (core-only, core+tail, or others) can be reconstructed afterwards. A second, "core+tail" measure (core plus the away-from-core `axion_gradient - axion_kinetic`, intended to cancel a propagating axion wave's equal kinetic/gradient contribution and isolate a static string's long-range tail) is implemented as a convenience diagnostic but its L-growth is **not** clean log(L) (closer to `L^0.6`); found that the axion kinetic term is exactly zero in this static test (`Pi` is always parallel to `psi` by construction, so `theta'=0` identically), so this measure is really just "core + total unmasked gradient energy" with no wave actually being subtracted -- the discrepancy is suspected to be either the two-vortex product-profile's amplitude between the cores not matching a true relaxed dipole (which the gradient identity is sensitive to via its `1/\|psi\|^2` factor), or `N_2=1` being too coarse right at the cores where `\|grad(theta)\|^2` is largest. Left as a known, explicitly deferred rough edge (user: "keep this in the back of our minds, but press on") rather than pursued further now. |
| 1.9 Spectra | Core machinery done (T3 passes); physical-unit rescaling and radial-mode spectra not yet wired up | `SpectrumKernel.hpp`'s `compute_spectrum` builds on `amrex::FFT::R2C`. Verified its round-trip normalisation matches conventions.md sec.9's stated convention exactly ("no normalisation, forward+backward multiplies by N^3"), so AMReX's raw `forward()` output can be used directly as `X~_p`, no correction factor needed (unlike `FourierIC.hpp`'s generation direction). R2C stores `kx>=0` only (Hermitian symmetry); every `kx` plane except `0` and `Nx/2` has an unstored conjugate partner at `Nx-kx` and is counted twice, to recover sums over the *full* spectrum. Binning follows sec.12 exactly (`modeIndex=floor(p_mod+0.5)`, `4*pi*p_mod^2\|X_p\|^2` per shell, dropped outside the inscribed sphere). `MaskedFieldBuffer.hpp` applies masking at exactly the one place CLAUDE.md constraint 5 asks for -- filling the buffer handed to the FFT, reusing `Masking.hpp` unchanged. **Resolved an `N^6` question with the user**: task 1.8's real-space/spectral cross-check formula (`(1/2N^6)*Sum\|a_dot(p)\|^2`) looked like it should have `N^3` (standard Parseval for this convention); redoing it carefully showed `N^6` is correct -- a spatial *average* already divides by `N^3` once, and Parseval divides by `N^3` again, so `<a_dot^2> = (1/N^6)*Sum_p\|a_dot(p)\|^2` exactly. T3 (a known-amplitude, known-mode plane wave, checked immediately after `initData` rather than after any evolution step, since `R(tau)` and hence the field itself moves on after even one step) **passes exactly**: real-space `<a_dot^2>`, Parseval from the full-cube sum, and Parseval from the inscribed-sphere sum all equal the hand-computed expected value `0.5*C^2` to machine precision, at both a low mode (`p0=4`) and one that aliases down from outside the naive Nyquist range (`p0=20` on a 32-point grid aliases to `p0'=12`, a correct consequence of `kx` being capped at `Nx/2` by the sampling theorem, not a bug). The shell-binned integral (`Sum_s S(s)/N_total^2`) is only *approximate* even after the same `N_total^2` normalisation (caught this: it was missing from an early diagnostic print) -- a single anisotropic mode shares its shell with other same-radius lattice points that carry zero power, diluting the angular average; unlike the full-cube/inscribed-sphere sums (which include every mode exactly once, so match the real-space value exactly), this is inherently a cross-check, not an identity. Not yet done: rescaling the raw `(p, S(p))` output to physical `(k/H, v^-3 d(rho)/dk)` units (sec.10/sec.12's `(2pi/(LH), (1/2pi) R L_tilde/N^6)` factors -- not yet re-derived/verified the way `rho_tot` and this cross-check were); radial-mode (`r`) spectra (sec.10: "follow identically" from the axion case, not yet implemented); wiring the axion spectrum into the routine per-snapshot diagnostics (currently opt-in via `axion_strings.compute_spectrum`, used so far only for the T3 IC-mode check and as an ad hoc `specific_post_timestep` addition). |
| 1.10 Remaining diagnostics | Velocities and curvature done; loop-finding still deferred | `Velocity.hpp` implements conventions.md sec.8's `gamma^2 v^2` estimator (global-string profile coefficients `c1=0.41222`, `e1=-0.025763`), evaluated at the corners of every pierced plaquette and averaged over the network. Sec.8's literal formula is written in terms of `psi`/`Pi` and a symbol `beta` that is not otherwise defined in the conventions snapshot available here; read as `beta = 1/b_inv` by analogy with the identical `(1/(b_inv tau))*psi` term in the main EOM/`rho_tot`, and **verified rather than assumed**: converting to the physical field `phi=psi/v` and its cosmic-time derivative makes every factor of `R` in the literal formula cancel exactly (shown algebraically in the header comment, and confirmed to 1e-9 numerically in `tests/test_velocity.cpp` across 6 diverse `(R,tau,b_inv)` points) -- a physical velocity should not depend on the comoving rescaling used to evolve the field, so the clean cancellation is itself a strong consistency check on the `beta` reading. The production kernel (`VelocityKernel.hpp`'s `compute_velocity_at_pierced_corners`) works directly in the R-independent `phi`/`phi_dot` form (reusing `StringFinder.hpp`'s exact 7-corner/3-plaquette structure) rather than recomputing the cancellation on every cell; the literal `psi`/`Pi` form is kept only as the cross-check. Wired into `specific_post_timestep` alongside the existing xi/energy/spectrum diagnostics (`mean_gamma_sq_v_sq`, `mean_gamma`, `n_velocity_corners` -- printed and saved to `network_scalars.dat`); the zero-strings case (`N_corners=0`) is handled explicitly rather than dividing by zero. Smoke-tested against both `straight_string_test` (T1: small, finite, non-negative `<gamma^2 v^2>` throughout a run of a near-static string, as expected) and the default homogeneous IC (no strings: exactly `<gamma^2 v^2>=0`, `<gamma>=1`, `N_corners=0`). Loop-finding (T4's other half) is still not implemented. **Curvature** (2026-09-26/28, with the user; opt-in via `axion_strings.compute_curvature`, method selectable via `axion_strings.curvature.method`, params file `params_curvature_smoke_test.txt`): `CurvatureKernel.hpp` bins a distribution of local string curvature (`kappa/m_r`, dimensionless) per snapshot to `curvature_distribution.dat`, via three independent, runtime-switchable methods rather than one trusted blindly -- this measurement turned out subtler than it first looked. A first version (three-point Menger curvature on raw face *centres*) was bimodal by construction: on a well-resolved lattice a chord between adjacent pierced faces is forced into exactly straight or bent-by-the-lattice's-own-turning-angle, never anything between (confirmed live: 0 of 720 histogram bins ever populated across a real network run). `tangent_vector` (Frenet `\|dT/ds\|` from `grad(psi1) x grad(psi2)`, `TangentField.hpp`) and `interpolated_position` (Menger curvature on `PlaquetteCrossing.hpp`'s sub-grid bilinear zero-crossing positions, not face centres) both fix this independently; a spot check reconstructing real string segments (`CurvatureChainDebug.hpp`, opt-in debug tool) confirmed the fix tracks visibly-real bends and straight stretches, and surfaced a genuine numerical caveat specific to `interpolated_position`'s tail (near-coincident consecutive crossings inflating kappa on a near-degenerate triangle) that is documented rather than silently patched. `hessian_analytic` (`HessianCurvature.hpp`: analytic single-point curvature of the curve `{psi1=0} intersect {psi2=0}` from psi1/psi2's gradients and Hessians, verified against a circle of known radius) is a third, structurally independent construction with no chain or chord at all; its finite-difference stencil was upgraded to 4th order (reusing GRTeclyn's own `FourthOrderDerivatives`, already the project's standard elsewhere), measurably tightening its agreement with `interpolated_position` (mean disagreement across a full network run: 0.037 -> 0.012 dex, most of the improvement at late times/low statistics). Composite/multi-level (Stage 2, 2026-09-27): `compute_curvature_histogram` takes an optional coverage mask (`count_plaquettes`/`compute_velocity_at_pierced_corners`'s own convention), driven once per AMR level by a composite loop in `AxionStringsLevel.cpp` and merged, so refined regions contribute their own finer-scale bends instead of being silently measured at the coarsest resolution -- verified against a real multi-level run (genuine regridding via the string tagger): measurement totals track the independently-trusted `N_p` composite total throughout, including exactly when a finer level activates. Output has both `plain_count` and `length_weighted` per bin -- **only `length_weighted` is safe to use across levels**: its per-measurement weight is a genuine physical chord/cell length (correctly smaller for finer-level measurements, mirroring how `N_p_plain` is converted to a length before being combined across levels elsewhere in this file), while `plain_count` is a raw measurement count that overweights whatever a finer level covers, since a fixed physical length gets chopped into proportionally more (smaller) counts there. See "Known issues" below for an open, unresolved question about how the measured curvature compares to the naive Hubble-scale expectation. |
| Output infrastructure (cadence, restart-safety, persistence) | **Done** | Follow-up to task 1.10, prompted by the user asking what output infrastructure existed and how it should work for a multi-run ensemble (2026-09-18). Three changes: (1) **Full-grid outputs off by default** -- `params_test.txt` now sets `amr.plot_int = amr.check_int = -1`; `Mode::PreEvolution` already hard-requires `check_int >= 0` for its handoff checkpoint (`AxionStringsParams::check_params`), so that path is unaffected. (2) **A real output cadence**, replacing "every coarse step": the full per-snapshot diagnostics (plaquette count, energy/velocity reductions, the spectrum FFT) are gated behind a cheap check of `log(m_r/H)` (free from the analytic `Background`, no grid pass) against `s_next_output_log_mr_over_h`, which starts at the new required `axion_strings.output_first_log_mr_over_h` and is bumped by `axion_strings.output_delta_log_mr_over_h` (also required -- no default fixed here, since that is a physics-adjacent numerical choice) each time it is crossed, via a `while` loop so a step that jumps past more than one threshold still lands correctly rather than drifting. (3) **Two real bugs fixed, not just flagged**: `network_scalars.dat` used `SmallDataIO`'s "old" constructor, whose `first_step = (time == dt)` heuristic assumes time starts at 0 -- ours starts at `tau_i` -- so a fresh run silently *appended* to (rather than renaming-old and overwriting) a pre-existing file, and the hardcoded `restart_time = 0.0` meant a genuine restart never deduplicated the redone tail. Fixed by detecting a real restart via `amr.restart` (not our own tau-relative clock, which can't distinguish "fresh start" from "restart" since `GRAmr::get_restart_time()` is 0 in both cases) and passing the actual restart time (`tau_i + get_restart_time()`) plus a correct `first_step`; `remove_duplicate_time_data()` is now called before every write (a verified no-op on a fresh run). The spectrum -- what a spectral index `q` actually gets measured from -- was previously only ever `amrex::Print()`'d, never saved; it is now persisted to a new `axion_spectrum.dat`, one row per mode index per snapshot with `tau` repeated as the literal first column (not blank-line-separated blocks, so the same flat restart-dedup logic applies exactly and safely). `AxionStringsParams::check_params()` now aborts at start-up if `axion_strings.compute_spectrum` is on without `axion_strings.masking.scheme = B`, since scheme A has no top-hat at all -- the persisted spectrum must be the genuinely screened field (the `plane_wave_test` IC's own one-shot spectrum check in `initData()` is unaffected: it hardcodes `MaskingScheme::None` locally regardless of this setting, since it is testing the FFT/Parseval identity, not screening). `tension_core_only`/`tension_core_plus_tail` are now also saved as columns (previously print-only) per the user: "we may as well save the processed tension, even if it can be reconstructed". Verified end-to-end against `straight_string_test`: (a) a coarse cadence (`delta=0.5`) demonstrably skips most coarse steps (3 snapshots out of 10, vs. 1 per step before); (b) a fresh rerun over pre-existing output files renames them to `.old.<random>` instead of mixing data; (c) restarting from an *earlier* checkpoint than the last snapshot and re-running forward past an already-recorded snapshot's tau produces no duplicate row in either file (checked both `network_scalars.dat` and the full 17-row `axion_spectrum.dat` block). One accepted wrinkle: a genuine restart always takes one extra, off-cadence snapshot exactly at the restart point (since `s_next_output_log_mr_over_h` resets to `first_log_mr_over_h` on every fresh process and immediately fires), which is harmless (no duplication, just one bonus sample) but not persisted across restarts -- left as is rather than adding checkpoint-side state for what is now an off-by-default, edge-case path. |
| Output infrastructure: unscreened spectrum + 2D visualisation | **Done** | Two follow-up asks (2026-09-18). (1) **Unscreened spectrum alongside the screened one**, "to allow for comparison and judging the impact of screening": `axion_spectrum.dat`'s columns are now `_screened`/`_unscreened` pairs (e.g. `shell_average_screened`, `shell_average_unscreened`), computed via two calls to `fill_masked_a_dot_buffer`/`compute_spectrum` per snapshot (screened using `s_energy_masking`, guaranteed scheme B; unscreened hardcoding `MaskingScheme::None`, mirroring the `plane_wave_test` IC's own check) -- doubles the FFT cost of an already opt-in diagnostic. Sanity-checked on `straight_string_test`: unscreened `<a_dot^2>` is consistently higher than screened, as expected once the cores are no longer masked out. (2) **Optional 2D visualisation snapshots**, off by default (`axion_strings.save_projection`), for spot-checking a run by eye rather than routine diagnostics: `ProjectionKernel.hpp`'s `compute_energy_projection` does a host-side (CPU-first, same call as `SpectrumKernel.hpp`'s binning: once per snapshot, not every substep) `MFIter` line-of-sight projection onto the xy-plane (z fixed as the line of sight), taking the **max** of the **unscreened** `rho_tot` per column -- max rather than sum/average so a thin core isn't diluted by a long quiet sight line; unscreened because screening would remove exactly what this is meant to show -- plus an overlaid `string_hit_count` (the xy-plaquette winding, `PlaquetteWinding.hpp`, summed over the line of sight) as an independent cross-check that the energy and string-finding pipelines agree on where the strings are. Saved to a new `axion_projection.dat`, one row per `(i,j)` pixel per snapshot (`tau` repeated as the literal first column, same restart-dedup trick as `axion_spectrum.dat`) -- a flat `(tau, i, j, value...)` table, straightforward to pivot into a 2D array for a coloured plot in any plotting tool. Verified end to end on `straight_string_test` (32x32 grid, 1024 rows/snapshot as expected): `string_hit_count` correctly picks out both vortex cores at their expected grid locations (count = `N_z` at each, one winding per z-layer, zero elsewhere). **Found a real bug this way, not just a visualisation quirk**: the projected energy's *global* max sat away from either core, in a smooth ridge on the far side of the box from both strings, exceeding the core's own value -- correctly flagged by the user as implausible ("the gradient energy density far from the core is much smaller... there should be no way for this to dominate"), and investigated rather than dismissed. **Root cause, confirmed both analytically (independent Python re-evaluation of the IC formula, no AMReX involved) and via instrumented debug prints of the actual running field**: `straight_string_test`'s phase field `theta = atan2(y-y1,x-x1) - atan2(y-y2,x-x2)` (`AxionStringsLevel.cpp`'s `initData`) is evaluated from literal, non-periodic `(x,y)` coordinates and is simply **not periodic** on the torus -- nothing in the formula knows the domain wraps at `y=Ly`. Both vortices sit at the same `y=y1=y2`, so the seam `y=0<->y=Ly` (diametrically opposite them) is exactly where this shows up worst: `psi2` jumps from `-0.81` to `+0.78` between two *physically adjacent* periodic grid points (a real discontinuity, confirmed present already in the pure, unevolved IC at `tau_i` via a from-scratch Python re-implementation of the exact formula -- not introduced by evolution or by any bug in the new energy/projection code, both of which were also directly ruled out by inspecting the actual ghost-cell values, which wrap correctly). That fake discontinuity then feeds a large, spurious gradient-energy contribution into the evolved solution from step 1 onward. **This is a distinct defect from the previously-documented (task 1.8) "product profile doesn't cleanly cancel between cores" tail issue** -- my first attempt at explaining this finding conflated the two; this one is specifically about the *phase* construction's non-periodicity, not the amplitude profile, and is large enough to be a real correctness problem for T1 (tension, and now core-location-by-energy), not a minor rough edge. Not yet fixed -- flagged to the user for a decision on the right periodic-aware construction (e.g. summing periodic images of each vortex's phase, or a different unwrapping) before relying on `straight_string_test` results, including the existing tension numbers, further. |
| Input parameterization: `log(m_r/H_i)` replaces `tau_i` | **Done** | Prompted by comparing our inputs against the user's old fixed-grid code (`main_global.cpp`): its `Hmri` argument (`H/m_r` at the start) is more physical than a raw conformal time, and the user asked to switch to it (2026-09-18). `Background::tau_from_log_mr_over_h` (new) is the exact inverse of the existing `H_over_mr_closed_form`, using `c_sched.c0` (the value in effect at/before the start -- same convention `apply_box_plan`/`pre_evolution_L_tilde` already use, and consistent with `H_over_mr_closed_form`'s own documented "only exact in the no-switch case" caveat, since any switch is expected to happen after `tau_i`, not before it). Singular at `c0 = a_inv` (Moore mode from the very start rather than reached via a switch) -- `AxionStringsParams::read_tau_i` (now takes a `const Background&`) checks for this explicitly and aborts with a clear message rather than dividing by zero silently. `axion_strings.tau_i` is gone; `axion_strings.log_mr_over_h_i` is the new required input. Verified: a new round-trip unit test (`test_background.cpp`, "`tau_from_log_mr_over_h` is the exact inverse of `H_over_mr_closed_form`") across several `(a_inv, c0, log_mr_over_h)` combinations, plus `tau_from_log_mr_over_h(0) == tau0`; `params_test.txt` updated to `log_mr_over_h_i = 0.0` (algebraically equivalent to the old `tau_i = 1.0` for its `a_inv=2, c0=0` configuration) and confirmed to reproduce byte-identical `tau`/`m_r/H` output at runtime; both new failure paths (`c0 = a_inv`, and the parameter simply missing) checked to abort with clear, specific messages rather than a generic ParmParse error. While reading the old code's full argument list to build this mapping, also matched up the rest of our `axion_strings.*` parameters against its `argv[]` positions for the user's reference (`N`/`N1`/`N2` for grid/core/Hubble resolution, `pre_evolution.gamma`/`seed` for pre-evolution, `dt_multiplier` for the old `R_delta`, `save_projection` for `domovie`) -- GW evolution, loop extraction, backreaction, the second "systematics" grid, and Moore's *extra* numerical trick (beyond just running with `c=1+b_inv`) remain unimplemented, matching what was already known/deferred. |
| Single continuous run: pre-evolution -> main handoff without a restart | **Done** | Replaces the restart-based handoff entirely (2026-09-18, with the user): a cluster running the requested `N=512` production job might not have enough memory/scheduling headroom for two separate job submissions bridged by a full-grid checkpoint, and since each pre-evolution relaxation is only ever used for one main-run realisation (ensembles need different ICs per member), nothing was actually lost by removing the restart path -- "keep the workflow simple and direct." Mechanism: `axion_strings.mode`/`Mode` (Main vs PreEvolution) is gone entirely; a run is now driven purely by `axion_strings.ic_mode`, with a new value `"fourier_relaxed"` that starts the run in a new internal-only `AxionStringsLevel::Phase::Relaxing` state (`s_phase`, not user-facing) -- everything else defaults straight into `Phase::Evolving`. `specific_post_timestep()`'s xi-monitoring loop, on reaching `xi_target`, now calls the extracted `apply_pre_evolution_to_main_rescale()` (identical rescale math to the old `specific_post_restart`, just invoked in place on the live `MultiFab` instead of after an external restart) and flips `s_phase` to `Evolving` directly, rather than stopping the run for a checkpoint. `okToContinue()` is now the authoritative stop condition throughout (previously it did nothing in Main mode, relying entirely on `Main_AxionStrings.cpp`'s `evolution.stop_time` check): it compares the live `tau` against a newly-exposed `axion_strings.derived_tau_f` (box planning's `tau_f`, informational ParmParse injection, read back in `variableSetUp()`), correctly spanning both phases without needing to know in advance how long relaxation takes, since the transition is in-place rather than a restart with its own separate clock. `evolution.stop_time` is now forced to `-1` (unlimited) by box planning whenever unset, both because it is no longer needed (`okToContinue()` replaces it) and to defeat a real trap found in the previous entry's restart-based version: GRTeclyn's own `BaseParameterChecker` silently defaults `evolution.stop_time` to `1.0` if left unset. `amr.check_int >= 0` is no longer required for any run (checkpointing is fully optional everywhere now, crash-recovery only). **Known, accepted limitation, not addressed here** (user: "there are likely to be complications for the Moore protocol, but I don't think these will be insurmountable"): box planning's Moore-mode branch (`c0 = a_inv`, or any configured switch) still does not derive `tau_f` at all (pre-existing limitation, unrelated to this change -- a fat->Moore run already needed `geometry.prob_extent`/`evolution.stop_time` set by hand before this); `okToContinue()` detects this (`s_has_tau_f == false`) and falls back to always returning 1, i.e. `evolution.stop_time`/`max_steps` set by hand are still the stop condition for that case, unchanged from before. Verified end to end at small scale (`N=16`, a lax `xi_target` for a fast test) before replacing the two production param files: relaxation runs, the transition fires and prints the same handoff diagnostics as before (`tau_pre_end`, `kappa`, etc.), the main evolution continues seamlessly in the same process, `network_scalars.dat` starts recording immediately once the output cadence's first threshold is crossed post-transition, and the run stops itself, with no user-set step count, at exactly the box-planning-derived `tau_f` (confirmed by hand: `s_tau_i` after the transition plus the final `a_time` reproduces `tau_f` to the last printed digit). `params_pre_evolution.txt`/`params_main.txt` replaced by a single `params_production.txt`. |
| **CRITICAL BUG (fixed 2026-09-18): every multi-rank diagnostic was silently wrong** | **Fixed** | Found while running the first genuinely multi-rank (`mpirun -n 4`) job of this whole project -- every earlier validation this session (T1, T3, box planning, the continuous-run refactor) used `-n 1`, which hid this completely. Root cause: `amrex::ReduceOps`/`ReduceData::value()` only reduces *within the calling rank's own boxes* -- confirmed directly from `amrex/Src/Base/AMReX_Reduce.H`: AMReX's own `Reduce::Sum`/`Min`/`Max` free-function wrappers never follow their `reduce_data.value(reduce_op)` with an `amrex::ParallelDescriptor::Reduce*Sum` call either, so a cross-rank total is always the *caller's* responsibility, not something the API does automatically. Three of our four `amrex::ReduceOps`-based kernels were missing this entirely: `StringFinder.hpp`'s `count_plaquettes` (`N_p`, hence every `xi`), `EnergyKernel.hpp`'s `compute_total_energy` (`rho_tot`, axion kinetic/gradient, `n_unmasked`, hence tension), and `VelocityKernel.hpp`'s `compute_velocity_at_pierced_corners` (`<gamma^2 v^2>`, `N_corners`) -- only `SpectrumKernel.hpp` and `ProjectionKernel.hpp` (both written with an explicit host-side loop + `ParallelDescriptor::ReduceRealSum`/`ReduceRealMax`/`ReduceLongSum` already) were unaffected. Confirmed the actual failure mode directly with a controlled test (`N=32`, 1 vs 2 vs 4 ranks, otherwise identical): `N_p`, `rho_tot_unscreened`, and `tension_core_only` each scaled down as almost exactly `1/n_ranks` -- each rank was reporting only its own local partial sum as if it were the whole domain's. Averaged (ratio-of-two-equally-wrong-sums) quantities like `rho_tot_screened` looked deceptively plausible despite being built from the same broken sums, which is what makes this class of bug dangerous -- CLAUDE.md's testing-discipline section calls out exactly this risk ("a bug that leaves a diagnostic plausible but wrong gets *harder* to spot as the ensemble grows"). Fixed by adding the missing `amrex::ParallelDescriptor::ReduceRealSum`/`ReduceLongSum` calls (batched into one 7-wide call in `EnergyKernel.hpp`) to all three kernels, matching the pattern already correct elsewhere. Verified: the same 1/2/4-rank test now gives bit-for-bit identical results at every rank count. **Practical implication**: any earlier local run in this project used `-n 1` and was therefore unaffected by this specific bug (its results stand), but this means the codebase had never actually been exercised at `n_ranks > 1` until now -- worth remembering if something *else* rank-count-dependent turns up later. |
| **Tension normalisation fix: physical energy/length, general (not just T1) string length** | **Fixed** | Found while producing tension/energy/spectrum plots for an `N=256` fat-string screened run (2026-09-18, with the user). `tension_core_only`/`tension_core_plus_tail` (`AxionStringsLevel.cpp`) divided the core energy by `string_length_in_box = 2*Geom().ProbLength(2)`, hardcoded for T1's exactly-2-straight-strings-spanning-z geometry -- silently wrong for any general network (many loops, random orientation), where it just tracks `N_p` diluting rather than measuring a real tension. Separately, dimensional analysis turned up a second, independent bug: `rho_tot_pointwise` (`Energy.hpp`) is the genuine *physical* energy density (energy/physical-volume, derived from the canonically normalised `phi` Lagrangian), so the physical energy in the core cells is `R(tau)^3 * dx^3 * Sum_core(rho)`, not `dx^3 * Sum_core(rho)` -- the old formula omitted the `R(tau)^3` entirely and divided by a *comoving* length, missing another factor of `R(tau)`; net effect, exactly one missing factor of `R(tau)^2`. This went unnoticed because T1's own validation (the `~3.72`, `L`-independent result recorded under task 1.8 above) was evaluated at `R(tau) approx 1`, where the missing factor is invisible. **Fix** (deliberately deviating from a literal reading of conventions.md sec.10's presentation, per the user's explicit go-ahead to do so as long as the choice is documented here): `mu = R(tau)^2 * dx^3 * Sum_core(rho) / ell_comoving`, with `ell_comoving = (2/3) * N_p * dx` -- the same plaquette-count -> length relation `XiFormula.hpp` already uses for `xi`, i.e. a *statistical* average valid for a randomly-oriented network, explicitly **not** exact for T1's single deterministically axis-aligned string (T1 used the exact `2*L_z` for that reason, and does not have an automated regression test tying it to a specific tension value, so this change does not silently break CI -- but T1's previously recorded `~3.72` figure is now stale and would need re-deriving against the new formula if that specific check matters again). Verified on the `N=256` screened re-run: `tension_core_only` now sits in a roughly flat `O(7-9)` band after the initial transient (the core-localised contribution, no longer diluting to zero over the run), while `tension_core_plus_tail` grows slowly and roughly logarithmically with `log(m_r/H)` (`~13` to `~23` over `log(m_r/H) in [3, 5.5]`) -- qualitatively the expected `mu ~ mu_core + pi v^2 ln(...)` shape for a global string's long-range Goldstone tail, a much more physically sensible result than the old formula's monotonic ~300x decay to zero. Unit tests (39/39, `tests/`) unaffected (no test currently covers `AxionStringsLevel.cpp`'s tension block directly). |
| **Spectrum physical-unit rescale (`k/H`, `v^-3 drho_a/dk`): re-derived from scratch and verified against the energy pipeline; the documented conventions.md formula was wrong** | **Done -- and corrected upstream in conventions.md sec.10 (2026-09-19, with the user)**, with a dated note and decision-log entry (sec.14) | Task 1.9's own status line above already flagged this rescale as not yet re-derived/verified -- correctly, as it turned out. First pass (2026-09-18) fixed only the *x*-axis: conventions.md sec.5's `L` (physical, `= R(tau)*L_tilde`) vs `L_tilde` (comoving) are distinct symbols, and an initial script used the comoving one in both the `k/H` and amplitude factors, putting every spectrum's peak 1-2 decades too high in `k/H`. That fix was necessary but not sufficient. **Second, more serious bug found 2026-09-19** while adding the `H f_a^2`-normalised version (user: "I'm not convinced by the normalisation... please check"): `MaskedFieldBuffer.hpp` fills the FFT buffer with `theta'` (`d(theta)/d(tau)`, comoving conformal-time phase rate, via `masked_a_dot`), *not* the physical axion time-derivative `a_dot = f_a*theta'/R` -- despite the code/diagnostics calling it "a_dot" throughout, which is a misnomer that obscured the gap. Converting `theta'` to `a_dot` needs an explicit `(f_a/R)^2` factor; conventions.md sec.10's documented rescale `(1/2pi) R L_tilde/N^6` has `R` to the wrong power (`+1`, not `-1`) and no `f_a` dependence at all -- confirmed genuinely wrong, not just unverified, via a direct numerical cross-check against the independently-validated `rho_axion_kinetic_screened` (`network_scalars.dat`): `integral(drho_a/dk) / (2*rho_axion_kinetic_screened)` drifted from `8.7` to `122` across a run using the documented formula (wrong, growing `R`-dependence), a clear signature of a missing/wrong power of `R`, not an accidentally-right constant. **Re-derived from scratch**, tracking the `4*pi*p^2` shell-binning convention (`SpectrumKernel.hpp`) through to a genuine `d^3k` integral, and using the already-established, already-trusted identity `Sum_s S(s) ~= inscribed_sphere_energy` (documented as an *approximate* cross-check, task 1.9's own entry above) as the anchor: `v^-3 drho_a/dk (p) = S(p) * f_a^2 * L_tilde / (2*pi*R*N^6)` (replacing the old `S(p)*(1/2pi)*R*L_tilde/N^6`). Verified: the same integral-vs-`rho_axion_kinetic_screened` ratio now sits at `1.00 +/- 0.02` for `log(m_r/H) >~ 4` on both the `N=256` fat and `N=384` physical runs, degrading only for the very first (transient, right-at-handoff) snapshot -- consistent with, not worse than, the pre-existing "approximate, not exact" caveat on the shell-binning identity itself. A useful side confirmation: dividing this corrected quantity by `H f_a^2` makes the `f_a^2` cancel analytically, exactly the expected/intended feature of that normalisation (removing the free `f_a` dependence to compare spectral *shape* across otherwise-unrelated axion-mass choices) -- the old, wrong formula did not have this property. (`axion_spectrum.dat`'s raw `p`/`shell_average_*` columns were unaffected and still correct throughout -- only the downstream rescale was wrong.) **conventions.md sec.10 corrected accordingly** (2026-09-19, with the user's go-ahead to amend it directly): the amplitude factor now reads `f_a^2 L_tilde/(2*pi*R*N^6)`, with a dated note explaining the error and a decision-log entry (sec.14) -- this project's living document and this repo's own formula are no longer diverged. |
| **Radial (Higgs) energy components: kinetic, gradient, mass** | **Done** | Task 1.8's own status line above flagged this decomposition as not yet implemented (conventions.md sec.12 names it but gives no psi/Pi formula). Requested and implemented 2026-09-18/19, with the user. `Energy.hpp` gains three new pointwise functions, defined analogously to the existing axion (phase-only) kinetic/gradient functions -- the exact physical energy of the `\|phi\|` amplitude degree of freedom alone, not "total minus axion" (which would also pull in cross/interaction terms conventions.md separately names and which were not requested). Derived and unit-tested as an *exact* orthogonal (radial/tangential) decomposition of the full kinetic and gradient energy -- not an approximation: `d\|phi\|/dt` and `grad\|phi\|` are the projections of `phi_dot`/`grad(phi)` onto `psi`'s own direction, with the orthogonal (tangential) projections reproducing `theta_prime`/`grad(theta)` exactly (the `(1/(b_inv tau))psi_i` terms cancel algebraically in the cross product) -- 5 new doctest cases in `tests/test_energy.cpp` check this decomposition reconstructs the independently-computed full kinetic/gradient energy to `1e-10`, for several non-vacuum `(psi,Pi)` configurations, plus the `psi=0` guarded-core case and the mass term's vacuum/off-vacuum values. `EnergyKernel.hpp`'s `TotalEnergyResult`/`compute_total_energy` extended to a 13-wide reduction (was 7) computing screened+unscreened radial kinetic/gradient/mass alongside the existing quantities (screened was the explicit ask; unscreened added for consistency with every other energy component already saved both ways). Wired into `network_scalars.dat` (6 new columns) and the per-snapshot `amrex::Print()`. Verified end to end on both the `N=256` fat-string and `N=384` physical-string screened runs below: radial kinetic and mass track each other closely (equipartition, as expected for an oscillating massive mode) with gradient somewhat lower, all three screened consistently below unscreened, no NaN/Inf. |
| **Physical-unit plot normalisation: energies in `H^2 f_a^2`, spectrum in `H f_a^2`** | Done, as an analysis-script step | Requested 2026-09-19: the natural units for these observables are `H^2 f_a^2` (energies) and `H f_a^2` (`v^-3 drho_a/dk`), `f_a = sqrt(2) v = sqrt(2)` in code units. Implemented in the scratchpad plotting script (not the C++ output -- consistent with the project's existing "save raw components, rescale downstream" pattern already used for `xi`/tension/spectrum): each energy/spectrum figure is now produced twice, raw and rescaled. Tension is *not* rescaled by these units -- it has different dimensions (`v^2`/mass^2, vs. `H^2 f_a^2`'s mass^4) so the requested rescale does not apply to it. Implementing the spectrum's `H f_a^2` version is what surfaced the pre-existing `drho_a/dk` unit bug documented in the row above -- the energy-side `H^2 f_a^2` normalisation had no such issue (it rescales `network_scalars.dat`'s already-correct, already-validated `rho_*` columns directly, nothing new to derive). |
| **Reusable analysis setup + instantaneous emission spectrum F(k/H, m_r/H)** | **Done** | Requested 2026-09-19 ("all these types of plots are going to have to be made many times... let's have a nice setup") together with a request to plot Fleury & Moore 1806.04677's `F` (the shape of the axion emission spectrum at a single instant, sec.4.2.1 eq.33). New `AxionStrings/analysis/` (in-repo, not scratchpad): `axion_analysis.py` consolidates `Background(tau)`, named-column file loaders (`network_scalars.dat`'s column layout has already changed once this project; positional indexing was an accident waiting to happen), and the (2026-09-19-corrected) spectrum unit conversions into one place, plus a new `instantaneous_emission()`. **Derivation** (not spelled out in the paper -- done here to make the implementation checkable): substituting the comoving wavenumber `kappa = k*R(t)` into the paper's eq.(23) collapses its time-integral's dependence on the upper limit, giving `d/dt[R^3 drho_a/dk]|_kappa = (Gamma/H) R^3 F(k/H, m_r/H)` -- eq.(33), with the derivative taken **at fixed comoving mode index `p`** (`kappa = 2*pi*p/L_tilde`, `L_tilde` fixed per run), not at fixed physical `k` -- i.e. "redshift the earlier snapshot forward" *is* comparing the same `p` across both snapshots after `R^3`-weighting, not an interpolation in `k`. Normalises to `integral(F dx)=1` directly (mirroring the paper's own approach) rather than computing `Gamma(t)` independently (eq.17's general form is complicated; the normalisation-to-1 requirement sidesteps needing it, for `F`'s shape specifically). Cosmic-time `Delta t` (the finite difference is in `t`, not `tau`) uses a closed-form `t_cosmic(tau) = tau^(a_inv/b_inv)/a_inv`, verified against direct numerical quadrature of `R(tau)`. **Tested against an independent construction, not just self-consistency**: `test_axion_analysis.py` forward-integrates a *known* `F_test(x)` and constant `Gamma/H` through the defining eq.(23) integral (via `scipy.integrate.quad`, a completely separate code path) to synthesise `S(p)` at two close snapshots, then checks `instantaneous_emission` recovers `F_test` -- as a *convergence* check (error shrinks as `Delta log -> 0`, from `0.75` to `0.06` over a `10x` reduction), not a single-tolerance check, since a steep Gaussian test function's tails are expected to show a large but shrinking finite-difference smear even with a correct implementation (checked directly: this smearing is a real, expected numerical-differentiation artefact, not a bug -- confirmed by varying `Delta log` and watching the error track it). Verified end to end on both the `N=256` fat and `N=384` physical runs: fat-string `F` clearly shows the paper's IR peak at `k/H~5-15` and the different-time curves cluster together (scaling behaviour) for `log(m_r/H)>~4`, matching Figure 11's fat-string panel qualitatively; physical-string `F` is visibly noisier (`18%` of points have a negative raw derivative, vs `5%` for fat), matching the paper's own observation that the physical case is harder to extract cleanly. Negative-derivative points (1806.04677 sec.4.2.1: "subject to fluctuations at frequencies near the core") are flagged (dotted, de-emphasised) rather than silently dropped or trusted. Not yet done: fitting the power-law index `q` from `F` (paper's Appendix E) -- the `F` extraction itself was the immediate ask; fitting `q` is a natural next step using the same `instantaneous_emission()` output. **Follow-up (same day)**: the standard per-run plot set (`xi`, energies raw/normalised, tension, spectrum raw/normalised -- previously a one-off scratchpad script) was likewise consolidated into `AxionStrings/analysis/plot_run_summary.py`, built on the same `axion_analysis.py` utilities rather than its own copy of `Background(tau)`/unit conversions. Verified pixel-identical output (`max pixel diff = 0`) against the scratchpad script's already-validated plots on both the `N=256` fat and `N=384` physical runs before treating the scratchpad version as retired. |
| **N=384 physical-string (`c0=0`) test run: `log(m_r/H_i)=3`, `xi_target=1`** | **Done** | Requested 2026-09-19 as a cross-check of the tension/energy pipeline against a second evolution scheme. Box planning for physical mode differs qualitatively from the fat-string demo runs: `m_r` is constant (`lambda(tau)=lambda0` always), and `N2=1` only fixes resolution at the *final* time, giving a much finer comoving grid (`dx=1/sqrt(N)=0.05103` for `N=384`, vs fat mode's `dx=1`) -- `pre_evolution.gamma` generalised accordingly to `1/dx = sqrt(N) = 19.596` (previously hardcoded to `1` for the fat-mode case, where `dx` happens to equal 1). **Performance note, not yet resolved**: this run took approximately 8 hours (overnight) at `~15-20s/coarse-step`, far slower than the `N=256` fat run's `~0.2-0.3s/step` -- more than the `3.375x` cell-count ratio predicts. Per-rank CPU usage was only `11-16%` during a live check, suggesting an MPI/memory-bandwidth bottleneck (`N=384` with `max_grid_size=64` gives 216 small boxes over 4 ranks, vs `N=256`'s 64 boxes) rather than a raw compute limit -- flagged to the user mid-run rather than assumed; not investigated further since the run finished successfully overnight. Worth profiling (e.g. `max_grid_size=128`) before relying on `N=384`+-scale runs routinely. **Physics results, verified clean (no NaN/Inf, exact stop at `tau_f`)**: `xi` stays in `[0.75, 1.05]` throughout -- tighter around the attractor than the fat-string run's `[0.65, 1.18]`. `tension_core_only` (`~7.5-8`) and `tension_core_plus_tail` (`~17->24`, growing roughly logarithmically with `log(m_r/H)`) closely match the `N=256` fat-string run's values -- a strong cross-check that the corrected tension normalisation measures something physical, independent of the `c(tau)` evolution scheme used to evolve the network (tension is a property of the string profile itself, not of the numerical scheme). One isolated single-snapshot spike in `rho_axion_grad_unscreened` around `log(m_r/H)~4.3` (visible in the raw-units energy plot) was not investigated further -- plausibly a genuine transient (e.g. a small loop event) rather than a bug, but not confirmed either way. |
| **Fat->Moore switch: derived box planning, auto-stop, and a Moore-phase output cadence fix** | **Done** | Requested 2026-09-19 ("start with the fat string system... evolve until the log(m_r/H) that will be the constant value... switch to the Moore approach... box size chosen to optimise time spent in Moore, with N2 an input"). Design confirmed with the user first: `N` (box size) is the input, achievable dynamic range `D` is reported (not targeted); `N1` stays an explicit input (a genuine statistics-vs-dynamic-range tradeoff, not something to default silently). **Key physics, worked out before writing any code**: `R(tau)*m_r(tau)` is *exactly* constant through the whole fat phase (any `a_inv`), not merely "set at the start" -- so `dx = 1/(N2*R(tau_switch)*m_r(tau_switch))` is representative of the entire fat phase, evaluated anywhere. Separately, `H(tau) = R'(tau)/R(tau)^2` does not depend on `c` at all (only `lambda` does), giving a closed form for where Moore should stop: `tau_end = tau_switch * exp(D*b_inv/a_inv)`, `D` from the existing `moore_max_log_dynamic_range`. Both are new `BoxPlan.hpp` additions (`compute_moore_box_plan`), unit-tested by plugging the derived `tau_end` back into `Background` directly (mirroring `test_box_plan.cpp`'s existing style) and confirming `H(tau_switch)/H(tau_end)` reproduces `D` and `R*m_r` really is constant elsewhere in the fat phase, not just at the switch. New `Background::H(tau)` (previously only `H_over_mr_*` existed) backs this, tested against `H_over_mr_direct(tau)*m_r(tau)` on both sides of a switch. **`AxionStringsParams.hpp` changes**: `read_background()` gains `axion_strings.log_mr_over_h_switch` as an alternative to raw `c1`/`tau_switch` (mirroring `log_mr_over_h_i`'s own rationale) -- resolved *before* `tau_switch` is finalised (using only `a_inv`/`c0`, already read at that point) via a throwaway `Background`, avoiding the chicken-and-egg problem of needing the *final* `Background` (which needs the switch resolved) to resolve the switch; `c1` is set to `a_inv` automatically, since reaching Moore is the only supported use of a switch. `apply_box_plan()`'s Moore branch (previously a validation-only check, "geometry.prob_extent and evolution.stop_time... set them by hand") now derives `gamma` from `tau_switch` (removing a previous silent-inconsistency risk between separately-supplied `gamma`/`tau_switch`), `geometry.prob_extent`, and injects `derived_tau_f = tau_end` via the *same* key the existing (non-Moore) continuous-run mechanism already reads back -- so `AxionStringsLevel::okToContinue()` needed **zero changes**. The `c0 = a_inv` from the very start (no switch) case is unchanged, since there is no fat phase for the `R*m_r=const` derivation to apply to. **Caught by the smoke test, not anticipated in the design**: the existing output cadence (`specific_post_timestep`) is driven entirely by `log(m_r/H)` crossing thresholds -- which is frozen by construction throughout Moore, so it fired exactly once (at the switch) and never again, silently defeating the entire point of extending time in Moore (more statistics at a fixed `log(m_r/H)`). Fixed by cadencing on `D_elapsed(tau) = log(H(tau_switch)/H(tau))` instead once `tau >= tau_switch` (same natural quantity `BoxPlan.hpp` uses for the achievable range), same spacing parameter, no new input; the *reported* `log(m_r/H)` in `network_scalars.dat` is unaffected -- correctly constant throughout Moore, only the cadence *trigger* changes coordinate. Verified end to end at `N=24` (`params_moore_smoke_test.txt`): printed `tau_switch`/`gamma`/`L_tilde`/`D`/`tau_end` all matched hand computation exactly; `m_r/H` tracked `tau` exactly through the fat phase then froze at `gamma` from the switch onward; ~27 snapshots recorded throughout Moore (previously 0) up to the derived stop, no NaN/Inf, clean exit. 46/46 unit tests pass (7 new). Separately confirmed with the user: `N` (lattice size) and `N1` (Hubble patches at the final time) are already required, independent inputs for standard (non-Moore) box-planned runs via the pre-existing `compute_general_box_plan` -- unchanged by this work. **Follow-up, full-scale `N=256` test (2026-09-19)**: `params_moore_256.txt`, switch at `log(m_r/H)=5`, reproduced every hand-computed box-planning number exactly (`tau_switch=148.413`, `D=1.090`, `tau_end=256`), froze `m_r/H` at the switch, recorded ~10 Moore-phase snapshots, no NaN/Inf, clean exit. **Found and fixed a second real bug this surfaced**, in the *analysis* tooling rather than the simulation: `AxionStrings/analysis/axion_analysis.py`'s `Background` (the Python reimplementation used for all plots) had no notion of a switch at all -- for every post-switch row it silently used the still-growing no-switch formula, an initial `m_r/H` cross-check mismatch of 65%. Fixed by extending `Background` with the same re-anchored `lambda(tau)` `CTauSchedule` has, plus a new `detect_background()` that recovers `tau_switch` *exactly* from the plateau already visible in a run's own `network_scalars.dat` (inverting the no-switch closed form at the plateau's own value) rather than needing the switch parameters passed in by hand and kept in sync -- both `plot_run_summary.py` and `plot_instantaneous_emission.py` now use it. 3 new Python tests (8/8 passing). **Presentation limitation above, fixed for the `xi` plot (2026-09-19)**: `plot_run_summary.py`'s x-axis is `log(m_r/H)`, which is exactly the quantity Moore freezes, so every Moore-phase point landed at the same x-coordinate (a correct but uninformative vertical cluster). Per the user ("let's plot in terms of log(tau) instead... for exactly the reason you identified"): the `xi` panel now plots against `log(tau)` whenever `bg.has_switch`, with a marked vertical line at the switch -- the fat portion is pixel-identical to before (`log(tau) == log(m_r/H)` exactly through the fat phase, the same `gamma = tau_switch` identity noted above), and the Moore portion now spreads out legibly from `log(tau_switch)` to `log(tau_end)`. The other panels (energies, tension, spectrum) still use `log(m_r/H)` and have the same vertical-cluster property during Moore -- not yet addressed, and `plot_instantaneous_emission.py`'s `F(k/H)` extraction is similarly not yet Moore-aware (it detects a switch and restricts itself to the fat-phase snapshots rather than mispairing Moore ones, but does not yet extract `F` *from* the Moore phase itself). |
| T1 static straight string | Not started | |
| T2 mask unit test | Not started | |
| T3 plane wave | Not started | |
| T4 collapsing loop | Not started | |
| T5 float vs double | Not started (build default changed, smoke-tested only) | 2026-09-19, with the user, ahead of moving to cluster-scale runs: `AxionStrings/GNUmakefile` now defaults to `PRECISION=FLOAT` (override with `make PRECISION=DOUBLE ...`; both configs coexist, precision-suffixed executable names). Rationale and risk are conventions.md sec.14's own existing decision-log entry, not new: float roughly halves state memory/`FillBoundary`/`FillPatch` volume, but headroom in `\|psi\|^2 - R^2` erodes with `log(m_r/H)` (~0.46 less reach than double); q/spectrum are the least-exposed observable, radial energy the early-warning one. Smoke-tested only so far, not the formal T5 comparison: the float build ran both `params_amr_validation_128.txt` (40 steps) and the full `params_circular_loop_test.txt` cleanly (no NaN/Inf), and the loop run's `loop_scalars.dat` matched the earlier double-precision run to ~7 significant figures throughout (e.g. `total_energy` at `t=10`: `3438.363` float vs `3438.365` double) -- consistent with ordinary float roundoff, not a bug, but not a substitute for actually running T5 (multiple `log(m_r/H)` values, comparing xi/spectrum/radial energy/total-energy conservation) before trusting float for a real production measurement. |

# Milestone 2: Adaptive Mesh Refinement

Requested 2026-09-19, using Buschmann et al. 2108.05368 ("Dark Matter from Axion Strings with
Adaptive Mesh Refinement") as a reference for what a real implementation of this involves.
Scope confirmed against `docs/conventions.md` sec.11 before starting: **physical mode (`c0=0`)
is the target** (the only mode where refinement demand and the falling cost of supplying it
move together); fat is restart-only viable; **Moore stays single-level**, unaffected by any of
this. A codebase survey found the tagging *plumbing* already exists (GRTeclyn's
`GRAmrLevel::errorEst()` -> `tag_cells()`) but the criterion is a geometric placeholder
(`AxionStringsLevel::tag_cells()`, explicitly commented as such); all five diagnostic kernels
(`StringFinder`, `EnergyKernel`, `VelocityKernel`, `SpectrumKernel`, `ProjectionKernel`) are
single-level only, and `specific_post_timestep()` hard-gates to `Level()==0`.

Planned phases: 0) AMR box planning (this entry). 1) Real (plaquette-based) tagging. 2) Composite
(cross-level) diagnostics, using **field-theoretic (not refinement-level) masking** for the
energy/xi observables -- confirmed with the user, deviating from Buschmann et al.'s own choice of
masking by finest-refinement-level, since our masking scheme (conventions.md sec.10) is already a
first-class, runtime-scannable convention and switching conventions between single-level and AMR
runs would itself be a systematic. 3) Spectrum on the hierarchy (average fine data down, FFT the
coarse level only -- already "settled" policy, sec.11). 4) Interface validation, including the
level-timing systematic conventions.md sec.11 calls "the central AMR risk." 5) Production
readiness (multi-level restart, performance).

| Phase | Status | Notes |
|---|---|---|
| **0. AMR box planning** | **Done** | New `BoxPlan.hpp::compute_amr_box_plan(N_effective, N1, N2, a_inv, c, max_level)`. Conventions.md sec.5's own note -- "for AMR, N and N2 refer to the effective finest resolution, not the base grid" -- means the existing `compute_general_box_plan` already gives the right `L_tilde`/`tau_f`/`dx_finest` when called with `N=N_effective`; what's new is deriving the *base* (level 0) grid size (`N_effective/2^max_level`, with a `-1` sentinel -- not a silently-wrong truncated integer -- if not evenly divisible) and the **level-addition schedule**: the `log(m_r/H)` at which each level must come online so the current finest level's resolution never drops below the `N2` target. **Derivation** (worked out from scratch, then verified numerically against a direct simulation of the underlying power law before writing any C++): `N2` at fixed comoving spacing scales purely geometrically with level (`N2_level_ell = 2^ell * N2_level_0`, independent of `c` or how `m_r(tau)` itself evolves -- a consequence of `dx_ell = dx_0/2^ell` alone), while `N2_level_0(tau)` is a power law in `tau` set by `R(tau)*m_r(tau) = R0*tau^((1-c)/b_inv)`. Converting to `x = m_r/H` and combining, the level spacing works out to `Delta log(m_r/H) = ln(2)*(a_inv-c)/(1-c)` -- singular at `c=1` (fat), matching conventions.md's own observation that fat mode's comoving core width is exactly constant, so no further level is ever needed there; for `a_inv=2, c=0` this reduces to `ln(4) ~= 1.386`, exactly reproducing conventions.md sec.11's number (previously stated only for that specific case) and Buschmann et al.'s own empirical level-addition spacings (`2.6, 3.9, 5.3, 6.7`, differences `1.3-1.4`). Thresholds are counted backward from `tau_f` (where, by `compute_general_box_plan`'s own construction, `N2` is reached using *all* `max_level` levels), one level-spacing per level. **Tested against independent constructions, not just the formula restated**: (1) numerically verified the `Delta log` formula against a from-scratch simulation of `N2_0(tau)`'s power law before any C++ was written; (2) `tests/test_box_plan.cpp` checks `N_base`/`dx_base`/divisibility, that consecutive `log_add` thresholds are spaced by the derived `Delta log` (including a **non-`a_inv=2`** case, `a_inv=2.5`, to exercise the general formula rather than only its special case), that the finest level's own threshold sits exactly one spacing before `tau_f`'s `log(m_r/H)`, and -- the actual physically meaningful boundary condition, checked directly against `Background`, not against `compute_amr_box_plan`'s own formula (which would be circular) -- that the **base grid alone**, evaluated at `tau_from_log_mr_over_h(log_add[0])`, gives exactly `N2` in units the user chose: `N1`/`N2`/`max_level` as explicit inputs, `L_tilde`/`dx_base`/`N_base`/the level schedule all derived. 49/49 unit tests passing (3 new). **Wired up (same day)**: `apply_box_plan()` now branches on `amr.max_level` -- `>0` uses `compute_amr_box_plan`, injecting `amr.n_cell` from `N_base` (cross-checked, not overwritten, if already set -- same pattern as every other derived quantity here), `geometry.prob_extent` from `L_tilde`, `derived_tau_f`, and a new `axion_strings.derived_level_add_log_mr_over_h` array for Phase 1's tagger to consume later; `max_level==0` keeps the exact prior single-level behaviour. Added a guard: `moore_mode && max_level>0` now aborts with a clear message, since conventions.md sec.11 is explicit that Moore should stay single-level. **Smoke-tested with GRTeclyn's real AMR machinery running for the first time in this project** (`params_amr_smoke_test.txt`, `N=64` effective, `max_level=2`, using the *existing placeholder* `tag_cells()` since Phase 1 hasn't replaced it yet): startup printout matched hand computation exactly (`N_base=16`, `L_tilde=8`, level thresholds `1.386`/`2.773`); the run then genuinely regridded at `lbase=0` and `lbase=1`, advanced real level-1 and level-2 steps (correct `ref_ratio=2` subcycled `dt`), and exited clean with no NaN/Inf -- confirms the box-planning-derived grid/geometry are consumable by GRTeclyn's actual multi-level driver (quartic interpolation, composite FillPatch), not just internally self-consistent numbers. Phase 1 (replacing the placeholder tagger with the real plaquette-based criterion) is next. |
| **1. Real (string-based) tagging** | **Done, with real findings flagged for follow-up** | New `StringTagger.hpp`, replacing `AxionStringsLevel::tag_cells()`'s geometric `FixedGridsTagger` placeholder entirely. **Primary criterion**: pierced plaquettes at the cell's low-index corner -- refactored `StringFinder.hpp` to expose the per-cell test (`plaquette_windings_at`) as a shared function first, so the tagger and the `xi` diagnostic call the *exact same* code (CLAUDE.md constraint 5), not a second hand-copied version that could silently drift. **Secondary criterion, implemented now rather than deferred** (per the user: "we might need to turn the gradient trigger on at some point, so it might be worth implementing this now"): Buschmann et al.'s `dx_ell^2|laplacian(psi_i)| > threshold` (their value 0.04), using the *same* low-level Laplacian access pattern already proven correct in `AxionStringsRHS`'s own EOM kernel, not a new one. New `axion_strings.tagging.gradient_threshold` (optional, `queryAdd` not `get`) defaults to `std::numeric_limits<double>::max()` -- genuinely off, not just "a large number that happens not to matter" -- so it can be switched on later (e.g. `= 0.04`) with no new code. Buffering (conventions.md sec.11: strings must stay a core width from any coarse-fine boundary between regrids) needed no new code either -- it is exactly AMReX's own `amr.n_error_buf`, already a project parameter, just previously always `1`. **Smoke-tested on real (not placeholder) refinement for the first time**, `N=128` effective/`max_level=2`/physical mode: first attempt used a raw, unrelaxed Fourier IC and found ~100% tagging at every level -- diagnosed as the IC's own deliberately-noisy lattice defects spuriously piercing nearly every plaquette (exactly the effect pre-evolution exists to remove, not a tagger bug), so redone with `ic_mode=fourier_relaxed`. Also found and fixed a smoke-test grid-configuration mistake, not a tagger bug: with `amr.max_grid_size` equal to the base grid's own size, level 0 is a single box, so *any* tagged cell forces refinement of the entire domain -- fixed by shrinking `max_grid_size` well below `N_base` so AMReX has room to make localised patches. **Critical gap found from the first smoke test, and fixed**: the findings below from that first pass (level 0->1 pinned near 100%; late-time runaway) turned out to trace back to `tag_cells()` never having consulted Phase 0's level-addition schedule at all -- it tagged unconditionally on the plaquette/gradient tests regardless of whether refinement was actually *due yet* at the current `log(m_r/H)`. Fixed by adding, at the top of `tag_cells()`: an early return during `Phase::Relaxing` (pre-evolution's `a_time` is `PreEvolutionBackground`'s own local clock, not `tau = s_tau_i + a_time`, so `log(m_r/H)` from the main `Background` would be meaningless there), then a check that `log(m_r/H)` at the *current level's* time has actually reached `s_level_add_log_mr_over_h[current_level]` (read back from Phase 0's `derived_level_add_log_mr_over_h` via `queryarr`) before permitting any tagging that would create the next level. **Gating verified at small scale first** (`params_tagger_smoke_test.txt`, `N_base=32`): a temporary debug print inside the gated branch, cross-checked against the log's own `"lbase = 1"` regrid trace lines, confirmed level 1 tagging starts at exactly `a_time=4.879` (`log(m_r/H)=2.192`, just past the `2.079` threshold) -- matching the hand-computed threshold almost exactly. (One methodology lesson from this check, worth remembering: AMReX's periodic `"Level N ... % of domain"` grid-summary print is not tied to every regrid event, so its absence between two log lines does not mean no refinement happened in between -- an apparent "level 1 not until step 186" discrepancy was purely this, not a bug, resolved by checking `"lbase = N"` trace lines instead.) **Then re-validated at the scale the user asked for** (`params_amr_validation_128.txt`: `N_base=128`, `max_level=2`, physical mode, `log_mr_over_h_i=2.0`): level 1 came online at `TIME=5.091` (`log(m_r/H)=3.477` vs. the predicted `3.466`, 0.3% off); level 2 at `TIME=10.748` (`log(m_r/H)=4.859` vs. predicted `4.852`, 0.15% off) -- confirming Phase 0's derivation and Phase 1's gating together, end-to-end, at a properly-resourced scale, and directly bearing out the user's own framing: "at the start of simulations no refinement is needed, but we should successively start refining more as log(m_r/H) increases." **A further, non-obvious finding from this run**: the *schedule* (which levels exist) increases monotonically as required, but the *refined fraction within an active level* does not -- level 1 covers 37% of its domain the moment it turns on and falls steadily to ~8% by `tau_f`, as the string network dilutes (`N_p_weighted` falls from ~22700 to ~420 over the run, per `network_scalars.dat`) and correspondingly fewer of the `blocking_factor=16` blocks get touched by a string. This is expected scaling-regime dilution, not a gating defect, and should not be conflated with the (correctly monotonic) level-addition schedule. Level-0-only `xi` (uncorrected for masking, since Phase 2's composite diagnostics don't exist yet) also stayed bounded this time (~0.44 rising to ~0.62, settling back to ~0.53) rather than diverging -- in contrast to the earlier `N_base=32` test, consistent with that earlier runaway having been a symptom of the missing gate plus an under-resourced base grid, both now addressed. One concrete cost data point: average evolution speed dropped from ~1520 to ~224 code units/h the step level 2 turned on. Gradient criterion still untouched by this validation (default off); a dedicated sensitivity study remains follow-up work, as does Phase 2 (composite, cross-level `xi`/energies with field-theoretic masking, needed to check `xi` properly rather than level-0-only). 49/49 unit tests still passing (unaffected -- verification here was by real multi-level runs). |
| **2. Composite (cross-level) diagnostics** | **Done** | Until now, `specific_post_timestep()` computed `xi`/energies/velocities from `Level()==0`'s own data only -- once refinement existed (Phase 1), that data was simply *ignored*, not approximated: the interesting, string-dense refined region contributed nothing. Fixed by making all three reductions genuinely composite, using the *existing* field-based mask unchanged throughout (the user's "field-theoretic masking" choice from the planning discussion, as opposed to Buschmann et al.'s own choice of masking by finest-refinement-level -- a cell's mask weight still comes from its own field value, at whatever level it lives on, never from which level it happens to be). **Coverage**: a new `gather_composite_levels()` (`AxionStringsLevel.cpp`, anonymous namespace) walks `0..parent->finestLevel()`, building each level's own `amrex::iMultiFab` coverage mask via `amrex::makeFineMask` (1 = keep, 0 = covered by a finer level) so every physical point is counted exactly once, at whichever level actually covers it -- not zero times (the old behaviour) and not twice. **Ghost cells**: each level's diagnostic copy is filled via `amrex::AmrLevel::FillPatch` (genuine composite, quartic-interpolated across a coarse-fine boundary -- the same interpolation the evolution itself uses), not a plain `FillBoundary`, which only exchanges same-level/periodic neighbours and would leave cells right at a coarse-fine interface stale -- precisely the kind of resolution-correlated bias conventions.md sec.11 warns AMR can introduce if not handled carefully. **Combining levels correctly**: a coarser cell/plaquette/corner represents proportionally more comoving volume or string length than a finer one, so levels are combined weighted by `dx_level` (plaquette length, `StringFinder.hpp`'s `count_plaquettes`; velocity-corner length-average, `VelocityKernel.hpp`) or `dx_level^3` (energy volume, `EnergyKernel.hpp`'s new `compute_total_energy_sums`/`compute_composite_total_energy`) -- never by raw cell/corner *count*, which would silently overweight finer levels (more, smaller cells for the same physical region). `xi`'s own formula is exactly linear in comoving string length (`XiFormula.hpp`'s existing `(2/3) N_p dx` term), so summing lengths across levels first and converting once (new `xi_from_comoving_length`, with the old `xi_from_plaquette_count` now a thin single-level wrapper around it, unit-tested to be identical) is *exact*, not an approximation -- proved as its own doctest case (`test_xi_formula.cpp`), not just asserted. `TotalEnergyResult.n_total`/`n_unmasked` change *meaning* under this refactor, from a raw point count (valid only because every cell shared one `dx` pre-Phase-2) to a comoving *volume* (count times `dx^3`, summed across whichever levels contributed) -- single-level callers see the identical numeric value either way (`n_cells*dx^3` vs `n_cells`), so the downstream tension formula (`AxionStringsLevel.cpp`) was adjusted to stop applying its own separate `dx^3` (the product is unchanged, only how it is factored) rather than double-counting the volume. The spectrum stays deliberately level-0-only (established policy, conventions.md sec.11/task 1.9) -- its own real-space/Parseval cross-checks were still dividing by `energy.n_total`, which would have silently become the wrong (whole-hierarchy) denominator once Phase 2 landed, so this was repointed at a new, explicitly level-0-only `level0_n_cells` instead. **Verified**: 51/51 unit tests passing (2 new, both on `xi_from_comoving_length`'s exactness). Real multi-level runs (`params_phase2_smoke_test.txt`, a deliberately raw/unrelaxed Fourier IC chosen specifically to trigger refinement almost immediately, since correctness plumbing -- not physical realism -- was what needed exercising quickly) confirmed a strong, independent invariant at every snapshot, across 2-level and 3-level hierarchies alike, under both default (scheme A) and screened (scheme B) masking: composite `n_total` equals `L_tilde^3` to full double precision (`1448.154688`, matching `11.313708...^3` exactly) -- i.e. the coverage masks partition the whole comoving domain with no gaps and no double-counting. With scheme B enabled, `n_unmasked` correctly differs from `n_total` and evolves sensibly (grows as the network dilutes); tension stays finite and well-behaved across level transitions. Not yet done: a real (non-throwaway, physically relaxed) end-to-end re-validation of `xi`/tension against the earlier Phase 1 physical-mode run now that they are composite (the `params_amr_validation_128.txt` numbers quoted in Phase 1's own entry above predate this fix and are level-0-only); Phase 3 (spectrum on the hierarchy) is unaffected and still pending separately. |
| **3. Spectrum on the hierarchy** | **Done** | The plan agreed when this milestone started ("average fine data down, FFT the coarse level only") is what this implements -- distinct from Phase 2, which is about *which cells contribute at all*; Phase 3 is about *what value* a level-0-resolution cell should hold once something finer exists above it. Confirmed first that this genuinely was still needed after Phase 2: nothing in this codebase's post-timestep flow keeps level 0's own stored state in sync with a finer level's more accurate solution in a refined region (grepped the AMR driver for an existing `average_down`/sync step -- none exists), so level 0's own data there is a real, independently-timestepped, *worse* solution, not a stale-but-equivalent one. New `collapse_to_level0()` (`AxionStringsLevel.cpp`, anonymous namespace, alongside `gather_composite_levels`) cascades `amrex::average_down` from the finest level down to level 0, one adjacent pair at a time (needed for `max_level=2`: level 2 into a working copy of level 1 first, *then* that corrected level 1 into level 0, not a single direct ratio-4 average, since level 1's own valid data in a level-2-covered cell is itself superseded and must be corrected first). Verified from AMReX's own `average_down` implementation, not assumed, that it only overwrites the *covered* subset of the destination (an internal `ParallelCopy` from a temporary matching just the fine footprint) -- so starting from a copy of `gather_composite_levels`'s own per-level state and cascading downward leaves every uncovered cell exactly as it already was, and only touches what actually has better data available. The result feeds the *existing*, otherwise-unchanged single-level spectrum pipeline (`fill_masked_a_dot_buffer` + `compute_spectrum`) -- masking still happens at exactly one place (CLAUDE.md constraint 5), only the field it operates on has changed. `energy.n_total`/`level0_n_cells` normalisation (already split apart in Phase 2, since energy's own `n_total` became a whole-hierarchy volume) needed no further change -- the spectrum's own point count is still level 0's, just better-informed data at each of those points now. **Verified**: 51/51 unit tests unaffected (this is pure AMReX-hierarchy plumbing, no standalone-testable math changed). A temporary debug comparison (`amrex::MultiFab::Dot`, removed after use) against a real 3-level run confirmed both halves of the expected signature: at `finest_level=0` the corrected and raw fields are *bit-identical* (`collapse_to_level0` is a true no-op with one level -- zero regression risk for every existing single-level physical/fat/Moore validation), and once levels 1/2 come online the corrected field's sum-of-squares is consistently *smaller* than the raw level-0-only field's (`|corrected|^2 < |raw|^2` at every checked snapshot) -- exactly the sign expected of a genuine averaging-down (a plain mean reduces variance; Jensen's inequality), not a no-op or a sign error. Real multi-level runs with `compute_spectrum=1` (screened, scheme B, and unscreened side by side) completed cleanly through a 3-level hierarchy with no NaN/Inf in any reported spectrum quantity (one *expected*, unrelated `0/0 = nan` in a log-only cosmetic ratio was traced to a fully-masked buffer at a very early, single-level snapshot -- reproduces identically with or without this change, not a Phase 3 artifact). Not yet done: comparing an actual measured spectral shape/`q` before and after this fix on a real (non-throwaway) physical-mode run, to quantify how much the pre-Phase-3 level-0-only spectrum was actually biased in practice. |

## Known issues / open questions

- **Open (2026-09-27/28, with the user): the measured curvature distribution's
  peak/mean does not track `2*pi*H(tau)*R(tau)/m_r(tau)`** -- the correct
  comoving-frame horizon-scale benchmark (derivation and cross-check against
  `Background.hpp` in `CurvatureKernel.hpp`'s header comment; two earlier
  frame-conversion mistakes were made and caught in this same investigation
  before landing on this formula, so treat any *new* use of it with the same
  scrutiny). Tracked across a full single-level network run
  (`interpolated_position`, `N=256`): the benchmark declines steadily
  (`log10(kappa/m_r) = +0.36` at `m_r/H~7` down to `-0.22` at `m_r/H~110`),
  but the measured peak stays pinned near `kappa/m_r ~ 0.25-0.4` throughout,
  and the mean plateaus around `tau~5.5` rather than continuing to decline --
  the gap between measured and benchmark grows from `~0.05` dex to `~0.47`
  dex over the run, i.e. the curvature scale becomes *increasingly*
  decoupled from the horizon scale as the run deepens, not less. A
  multi-baseline check (`CurvatureChainDebug.hpp`, walking wider chord
  separations along real reconstructed chains) confirmed mean curvature
  does decrease with a wider measurement baseline, consistent with real
  sub-horizon small-scale structure on the string (a genuine, well-
  documented phenomenon in the cosmic/axion-string literature) rather than
  a plotting or unit bug -- but the chains tested so far (tens of cells)
  are far too short to reach a baseline comparable to the horizon itself
  (hundreds to thousands of cells at these resolutions), so whether the
  curvature spectrum eventually *does* approach the benchmark at a large
  enough baseline is untested. Not yet understood whether the plateau is
  physical (small-scale structure saturates) or a resolution artifact of
  the nearest-neighbour/local measurement methods themselves.

- **Open, investigated (2026-09-19, with the user: "genuinely investigate
  and worry about... before proceeding"): a coarser base grid
  (`dx_base ~ 0.25`) makes the plaquette tagger fire on ~100% of the domain
  at level 1's very first onset, instead of the ~25-40% seen at
  `dx_base <= 0.177` (Phase 1's own validated runs).** First hit trying
  `max_level=3` at `N_base=128` (`params_full_test_1024.txt`) to extend
  Milestone 2's dynamic range beyond the validated `max_level=2` depth --
  `L_tilde`/`tau_f` scale as `sqrt(N_effective)` (`compute_general_box_
  plan`), not independently of it, so extending the hierarchy at fixed
  `N_base` necessarily coarsens `dx_base`. Five controlled tests, each
  changing exactly one variable, to isolate the cause:
  1. **Hierarchy depth ruled out.** Reproduced the identical ~100% tagging
     at `max_level=2` (not 3) by matching `params_full_test_1024.txt`'s
     `(dx_base, level-1 threshold)` exactly (`N_base=64`, `N_effective=
     256`, same `log_mr_over_h_i=2.0`) -- both give `level-1 threshold =
     ln(N_base) - max_level*ln(2) = 4*ln(2) = 2.773` by construction, and
     both showed 100%. Not a depth-3-specific bug.
  2. **Schedule-derived resolution margin: matters, but not sufficient.**
     The margin between `log_mr_over_h_i` and level 1's own threshold sets
     `N2_0(tau_i)` -- the base grid's points-per-core-width at the moment
     of handoff, `N2_0(tau) = 1/(dx_base * tau)` for `a_inv=2,c=0` (a clean
     closed form, cross-checked against the schedule's own power-law
     derivation and found identical). Raising `N2_0(tau_i)` from 1.47 to
     2.43 (lowering `log_mr_over_h_i` to 1.0, same `dx_base=0.25`) only
     partially helped: 100% -> 74%, still far from healthy.
  3. **Margin alone does not predict this.** A genuinely healthy run
     (`N_base=160`, `dx_base=0.158`, `N2_0(tau_i)=2.33`) and test 2's
     `N2_0(tau_i)=2.43` (marginally *better*) gave wildly different
     outcomes (~27% vs 74%) despite near-identical schedule-derived
     margins -- ruling out `N2_0(tau_i)` as a sufficient predictor on its
     own, whatever role it plays.
  4. **`pre_evolution.gamma`'s absolute value ruled out.** `gamma_pre` sets
     the fat-mode relaxation's own `m_r/H` -- with no explicit dissipation
     in this code (`sigma=0` by default, CLAUDE.md constraint 1), a lower
     `gamma_pre` means fewer oscillation cycles per Hubble time and weaker
     Hubble-friction damping of non-string defects, a plausible independent
     mechanism. Tested by decoupling `gamma_pre` from `dx_base`: held
     `dx_base=0.25` fixed, raised `gamma_pre` to 6.0 (matching validated
     runs) while lowering `k_max_over_mr` to 42.667 to hold the Fourier
     IC's own bandwidth (`k_max_cells = k_max_over_mr * gamma_pre *
     dx_base`) fixed at 64 -- isolating gamma_pre's dynamical role from any
     change in injected-noise bandwidth. Still ~100%. Ruled out as the sole
     cause.
  5. **A unit mismatch via `pre_evolution_L_tilde` ruled out.** That sec.7
     formula (a *different*, generally larger comoving box for the
     relaxation phase) is confirmed unused in this code path -- read the
     source directly: `apply_box_plan` computes it and prints it
     explicitly labelled "FYI only... NOT used", since pre-evolution and
     the main run share one live grid (`Geom()`) rather than communicating
     through a checkpoint. Not the mechanism.
  **Left standing**: the one factor common to every failing configuration
  and absent from every healthy one, across all five tests, is `dx_base`
  itself -- every failure used `dx_base=0.25`; every success used
  `dx_base <= 0.177`. Leading untested hypothesis: a coarser `dx_base`
  changes the *relative* numerical weight of the gradient-energy term
  (`~1/dx^2` via `FourthOrderDerivatives`) against the (`dx`-independent)
  mass term in the relaxation EOM, changing how the *same* initial noise
  pattern (confirmed statistically identical in grid-index space between
  the `max_level=2` and `max_level=3` `N_base=128` cases -- same seed, same
  `N_base`, same `k_max_cells`) evolves during relaxation. Not yet
  confirmed. Also worth weighing: `N2=1` is already flagged elsewhere in
  this document (Phase 1 entry) as "the project's established
  aggressive-minimum convention" -- `dx_base=0.25` may simply be the first
  configuration to fall on the wrong side of a margin that was already
  known to be thin. **Practical guidance until resolved**: treat
  `dx_base >~ 0.2` as suspect; the two genuinely validated configurations
  both used `dx_base <= 0.177`. Full test files and reasoning kept in
  `params_full_test_1024.txt`'s own header comment.

  **Follow-up (2026-09-19, same day, prompted by the user: "could a
  stronger cut on the max initial k improve the situation?")** -- yes,
  substantially, though the mechanism turned out more subtle than "less
  initial noise." A fine-grained trace (`axion_strings.pre_evolution.
  check_interval_coarse=1`, checking every step instead of the default
  50-step cadence) of the failing `dx_base=0.25`, `k_max_over_mr=64`
  config showed relaxation starting at `xi ~ 14-18` (`N_p ~ 1.5` million
  pierced plaquettes out of 2.1 million cells -- i.e. `k_max_over_mr=64`
  saturates the `N_base=128` Fourier grid's own Nyquist limit, `N_base/2
  =64` in FFT-index units, populating essentially every representable
  mode) and undergoing a genuine, real ~50x collapse in `N_p` over ~50
  light-crossing times before crossing the target -- not an artificially
  truncated relaxation as first suspected (the cadence-tightening logic,
  `AxionStringsLevel.cpp`'s `s_xi_check_interval` schedule, checked out
  correctly; the default 50-step cadence was just too coarse to *see* the
  trajectory, not too coarse to *let it happen*).
  Lowering `k_max_over_mr` while holding everything else fixed
  (`N_base=128`, `max_level=3`, `dx_base=0.25`, `log_mr_over_h_i=2.0`):
  `64 -> 100%` tagged at level-1 onset; `16 -> 91%`, declining slowly;
  `4 -> 47%`, declining to `35%` within the same window (comparable to the
  genuinely validated configs' ~25-40%); `0.5 -> N_p=0` at the very first
  check -- **no strings ever formed at all** (too little power to trigger
  the Kibble mechanism), confirming a real floor below which this knob
  cannot be pushed. So there is a working window, and the user's physical
  intuition ("no structure below the string core scale, `k <~ few * m_r`")
  is directionally right and empirically confirmed as the single most
  effective lever found in this investigation.
  What it is *not*, however, is simply "less initial mess": a matching
  fine-grained trace of the genuinely healthy `N_base=160` config (`dx_base
  =0.158`, same `k_max_over_mr=64`) started at `xi ~ 25-33` -- *higher*
  than the failing config's ~14-18 -- yet still relaxed to a clean, ~27%-
  tagged handoff. Both configurations start from an extremely dense,
  double-digit-`xi` tangle; only one collapses cleanly. This is consistent
  with (not yet proof of) the standing `dx_base`-dependent-EOM hypothesis
  above: a large `xi` built from *coherent*, large-scale windings (more of
  which populate at a properly band-limited `k_max_cells` relative to that
  grid's own Nyquist -- `N_base=160`'s `k_max_cells=64` is a real,
  enforced cutoff at 80% of Nyquist; `N_base=128`'s is not a cutoff at all,
  since it coincides with Nyquist exactly) collapses via genuine loop
  annihilation, while a large `xi` with a substantial *incoherent*,
  lattice-noise admixture (from populating literally every representable
  mode) may not collapse the same way even once the raw plaquette count
  crosses the target -- consistent with lower `k_max_over_mr` helping by
  removing exactly that admixture, independent of the starting `xi` value
  itself. Not fully proven; would need direct inspection of the field's
  own spectral content pre- and post-relaxation to confirm.
  **Updated practical guidance**: prefer `k_max_over_mr` a handful of
  units at most (this session's tests bracket a working window around
  `2-8`, not yet narrowed further), and treat `k_max_over_mr` approaching
  or exceeding `N_base/2` (the Fourier IC grid's own Nyquist limit) as a
  second, independent red flag alongside `dx_base >~ 0.2` above.

  **Second follow-up, same day: does refinement fire prematurely, before
  the schedule says it should (the user: "the main run starts at an
  m_r/H such that refinement probably isn't needed")?** Checked directly
  rather than inferred: `grid_places()`'s own `new finest: 0` prints
  persist through every regrid check (`amr.regrid_int=2`, so every 2
  steps) both through all of relaxation *and* for ~26 steps after
  handoff -- level 1 is not created until `log(m_r/H) ~ 2.78`, matching
  the analytic threshold (`2.773`) almost exactly. **The schedule-gating
  is not the bug and is not firing early.** But `log_mr_over_h_i=2.0`
  sits only `Delta_log=0.773` below that threshold (`tau_i=2.718` vs.
  `Delta_tau~1.3` to cross it) -- refinement becomes due very early in
  the run's own life, not after a long, gentle settling period, so
  whatever quality problem the field has gets exposed almost immediately.

  **Third follow-up, same day: does widening that margin (a *larger*
  `log_mr_over_h_i`, shrinking the time spent evolving un-refined on an
  increasingly under-resolved coarse grid, `N2_0(tau) = 1/(dx_base*tau)`
  decreasing monotonically) help, tested alongside a smaller
  `k_max_over_mr`?** Tested `k_max_over_mr=8` (between the previously-
  bracketed 4 and 16) and `log_mr_over_h_i=2.5` (margin 0.773 -> 0.273),
  both together and in isolation, same `N_base=128/max_level=3/dx_base
  =0.25` box throughout:
  - `k_max_over_mr=8` alone (`log_i=2.0` unchanged): 87% at onset,
    declining only to 70% within the test window -- a much steeper
    falloff between `k_max_over_mr=4` (47%->35%, healthy) and `=8` than
    between `8` and `16` (87% vs 91%). The earlier "working window of
    2-8" was too generous; the real transition sits close to 4.
  - `log_mr_over_h_i=2.5` alone (`k_max_over_mr=64` unchanged): 100%,
    persistently, across 9 consecutive regrid checks -- *not even
    beginning to decline*, arguably worse than the original `log_i=2.0`
    baseline (which did eventually start declining slowly).
  - Both together (`k_max_over_mr=8`, `log_i=2.5`): 100%, persistently.
  **Conclusion: raising `log_mr_over_h_i` does not help, and the evidence
  points the other way.** `N2_0(tau_i)` is a direct, monotonically
  *decreasing* function of `tau_i` alone (`=1/(dx_base*tau_i)`) --
  raising `log_i` unavoidably starts the main run at an already-worse
  absolute resolution (`N2_0(tau_i)=1.146` at `log_i=2.5` vs. `1.471` at
  `log_i=2.0`), and this appears to dominate over any benefit from
  spending less time exposed to further coarse-grid aging before
  refinement rescues it. The "two competing effects" framing from the
  previous discussion was itself the error: there is no real settling
  benefit to trade against, since this code has no dissipation
  (`sigma=0`) for "ordinary evolution" to lean on in the first place --
  only the resolution-degradation effect is real, and it argues for the
  *opposite* of a larger `log_mr_over_h_i`.
  **Revised guidance**: don't touch `log_mr_over_h_i` as a lever for this
  problem -- keep it set by the physics question being asked, not by
  this investigation. Narrow the `k_max_over_mr` search below 8, closer
  to (but confirmed above) the 4 value that worked; 4-5 is the best
  currently-known value for this box configuration, not the wider "2-8"
  bracket stated above.

  **Correction from the user, same day**: a larger `log_mr_over_h_i` may
  in fact be *better*, not worse -- at later times there are fewer total
  strings (network dilution), so more physical relaxation toward the
  attractor has already happened by handoff, which can outweigh the
  naive `N2_0(tau)` resolution-margin argument above. More fundamentally:
  this project cares about *late-time*, near-attractor behavior, not
  precise control of IC-generation nuisance parameters -- moderate
  variation in exactly how the initial tangle is generated is acceptable
  as long as the handed-off state starts close to the scaling attractor.
  This investigation is left here at a practically-adopted value
  (`k_max_over_mr ~ 4-5`) rather than further chased to a fully-resolved
  first-principles mechanism, consistent with that philosophy -- not
  because the remaining questions (why `dx_base` and `log_mr_over_h_i`
  affect things the way they empirically did here) are uninteresting, but
  because pinning them down further is not needed to proceed.

- **Tried and rejected (2026-09-19, with the user): a physically-derived
  default for `StringTagger`'s radial-gradient criterion (`axion_strings.
  tagging.radial_gradient_threshold`), calibrated from an isolated static
  string's own peak core gradient tied to the box-planning `N2` target
  (full derivation in `AxionStringsParams.hpp`'s `read_tagging_params()`
  comment), also produces ~100% tagging at level 1's first onset** -- a
  distinct mechanism from the `dx_base`-dependent plaquette-tagger issue
  just above (this hit `params_amr_validation_128.txt`, `N_base=128`,
  `dx_base=0.17678`, one of the *validated-healthy* configurations for the
  plaquette criterion alone, ~37% at onset there). Level 1 came online at
  exactly the predicted `TIME=5.091` (the log(m_r/H) *schedule-gating*
  itself is unaffected, as intended -- this criterion only changes which
  cells get tagged once a level is already due), but with this derived
  threshold active it tagged the *entire* domain immediately. Conclusion:
  a single isolated static core's peak gradient does not describe the
  ambient radial-mode gradient level in a real, dense Fourier-relaxed
  tangle -- the two are evidently comparable almost everywhere, not just
  at cores. Reverted to off by default (same convention as the existing
  `gradient_threshold`/Buschmann-Laplacian criterion), available as an
  explicit opt-in for anyone who wants to revisit the calibration --
  treat this finding as the starting point, not a reason to re-derive
  from scratch. The plaquette criterion remains the network schedule's
  only currently-tagging-by-default criterion.

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

- **Accepted, not fixed (2026-09-18, with the user): `straight_string_test`'s
  phase field is not periodic on the torus, and this is not just a cosmetic
  issue.** Found via the new `axion_strings.save_projection` visualisation
  (task: output infrastructure follow-up) -- the projected energy's global
  max sat away from either string core, which the user correctly flagged as
  physically implausible rather than accepting at face value. Root cause,
  confirmed independently of any AMReX/simulation code (a from-scratch
  Python re-evaluation of the exact IC formula) and via instrumented debug
  prints of the real running field (which also ruled out a ghost-cell bug --
  periodic wrapping was confirmed correct): `theta = atan2(y-y1,x-x1) -
  atan2(y-y2,x-x2)` is evaluated from literal, non-periodic `(x,y)`, so it
  is simply not periodic in `y` when both vortices share the same `y`
  (the `y=0<->y=L_y` seam, diametrically opposite them, is where this shows
  up worst). Tried the user's proposed fix (summing periodic images of the
  phase): **does not converge** -- a square truncation window plateaus at
  one residual, a circular window at a *different* one, the classic
  signature of a conditionally-convergent lattice sum (the same
  mathematical issue Ewald summation exists to solve for periodic dipole
  lattices; growing the truncation window doesn't help). Quantified how the
  residual scales with box size at fixed string separation and fixed grid
  resolution instead: the local psi-jump shrinks like `1/L`, so the local
  peak energy *density* and domain average shrink like `1/L^2` (fast enough
  that by `L/separation ~ 8` the seam's peak density already drops below
  the core's own, resolving the "impossible" visual finding) -- but the
  *total, summed* spurious energy (the quantity that actually feeds
  `tension_core_only`, since that is built from `n_total * average`, i.e. a
  sum) converges toward a roughly constant floor rather than vanishing, so
  growing the box does not make the T1 tension number itself trustworthy at
  high precision. **Decision: live with it.** A real string network (the
  actual production case) never uses this closed-form two-vortex construction
  at all -- ICs there come from `FourierIC.hpp`'s Gaussian random field, not
  an analytic vortex ansatz, so this specific non-periodicity mechanism does
  not apply there. Left unfixed rather than pursuing a rigorous periodic
  Green's function or an IC-relaxation pass, both real options if this ever
  needs revisiting. **Flagged to keep in mind**: if an unexplained energy
  excess, gradient anomaly, or "value doesn't peak where it should" symptom
  ever resurfaces in a *different* context (in particular anything using an
  analytic closed-form phase/winding construction rather than the Fourier-
  mode network ICs), this class of bug -- a smooth-looking formula that
  quietly assumes an infinite, non-periodic domain -- is worth checking for
  again before assuming it is something new.

- `Background::H_over_mr_closed_form` only equals `H_over_mr_direct` in the
  no-switch case (verified analytically and by unit test). Across a c(tau)
  switch, the sec.5 closed form `(tau/tau0)^((c-a_inv)/(a_inv-1))` does not
  hold as-is with the original `tau0` (conventions.md sec.12 flags "a
  separate branch for constant-log runs" for the related `log(m_r/H)`).
  `H_over_mr_direct` (from `R`, `R'`, `lambda` directly) is unaffected and
  always correct; nothing in the current code relies on the closed form
  across a switch. Worth resolving properly before the fat->Moore protocol
  (task 1.4+) needs a closed-form `log(m_r/H)` post-switch.

- **Trap for `Mode::PreEvolution` configs: `evolution.stop_time` silently
  defaults to `1.0` if left unset.** `apply_box_plan` deliberately does not
  set it in pre-evolution mode (it stops via the xi-monitoring loop
  instead, per sec.7) -- but GRTeclyn's own `BaseParameterChecker` (outside
  our code, standard for every run) independently `queryAdd`s a `1.0`
  default for it regardless. Found while validating
  `params_pre_evolution.txt` (2026-09-18): a run with `stop_time` left
  unset silently stopped at `t=1.0` rather than actually waiting for
  `xi_target`, with no warning. Fix: set `evolution.stop_time = -1`
  (unlimited, same convention as `evolution.max_steps = -1`) explicitly in
  any pre-evolution config -- now done in `params_pre_evolution.txt`,
  worth remembering for any future one.

- Box planning in Moore mode (`apply_box_plan`) only runs the sec.5
  dynamic-range check; it does not derive `geometry.prob_extent` or
  `evolution.stop_time` the way the general (non-Moore) path does, since
  milestone-1.md task 1.4 asks specifically for "the Moore-phase
  dynamic-range formula as a startup check", not a closed form for
  `L_tilde`/`tau_f` (sec.5 doesn't give one -- the general formula is
  singular at `c = a_inv`). A production fat->Moore run currently needs
  `geometry.prob_extent`/`evolution.stop_time` set by hand.

- **Resolved (2026-09-19, with the user): the supposed Moore-mode `H0^-2` xi
  factor does not apply and was a documentation error.** Previously flagged
  here as "must be added before xi is trusted for a Moore-mode run." Queried
  by the user ("the definition of xi is straightforward... could there just
  be a mistake in the documentation?") and checked against first principles:
  `R(tau)` (hence `t(tau)`) is independent of the `c(tau)` schedule entirely
  (only `lambda(tau)` depends on `c` -- the same identity the Moore
  box-planning/switch work already relies on), and `XiFormula.hpp`'s formula
  is built entirely from `R(tau)`/`t(tau)` and lattice quantities, so it
  already gives the correct physical xi in Moore mode unmodified -- there is
  no physical mechanism for a Moore-specific correction. The user then
  checked the actual source document directly and confirmed the cited factor
  ("source document eq. 84") is outdated and does not apply here.
  `conventions.md` sec.8 corrected accordingly (with a dated note and a
  decision-log entry, sec.14) -- no code change needed, since `XiFormula.hpp`
  never implemented the (spurious) factor in the first place.

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

- **Open, parked (2026-09-19, with the user): both flat-space loop ICs
  (`circular_loop` and `four_string_collision`) have a hard periodic-
  boundary phase discontinuity, a recurrence of the `straight_string_test`
  periodicity bug in a new context where it actually matters.** Found
  while investigating why both ICs' `N_p`/`string_length` (in the new
  `loop_scalars.dat` diagnostics) collapse to exactly 0 for a stretch of
  several light-crossing times before either intended collision/collapse
  physics has had time to occur, then (for `circular_loop`) partially
  recover. Confirmed by direct evaluation of each IC formula at a periodic
  seam (not just by inspecting the evolved run): for `four_string_collision`,
  `psi2` at `y=-10` vs. the periodically-identical `y=+10` differs by
  `~1.85` against a full amplitude of `~0.92` (essentially a full sign
  flip); for `circular_loop`, `psi2` at `z=-8` vs. `z=+8` differs by
  `~1.97` against amplitude `~0.98`, deep in the "vacuum" region far from
  the ring. Mechanism in both cases: a winding angle computed as
  `atan2(u, w)` from a *raw, unwrapped* coordinate `u` (respectively
  `gamma*(y +/- approach)` and `z`) that saturates to `+/-pi/2` for large
  `|u|` -- the two periodic images of the same point land on opposite
  sides of that sign-dependent saturation, so the phase (not just a small
  amplitude tail) jumps by order 1 exactly at the seam. Visually confirmed
  too: projections/slices of `rho_tot` show a bright artifact band at the
  affected boundary already at `t=0` (before any evolution), which then
  grows and radiates inward over a few light-crossing times, visibly
  smearing the intended string cores together well before `N_p` hits 0 --
  i.e. genuine numerical destruction of the winding, not a plaquette-
  detector/resolution artifact (the user's first hypothesis, checked and
  ruled out this way).
  **`circular_loop`'s own in-code comment, which claimed this construction
  was exempt ("not the straight_string_test's own periodicity bug... this
  is a smooth, quantitatively small amplitude residual"), was wrong** --
  written without checking, corrected once actually verified numerically.
  Both loop ICs share the same underlying issue as `straight_string_test`
  (`docs/STATUS.md`'s and the project memory's existing entry on that),
  but the earlier "live with it" decision for that case relied on real
  network ICs coming from `FourierIC` instead -- an escape that does not
  apply here, since the loop-collision/collapse study *is* built entirely
  from this class of closed-form construction.
  Real fixes considered but not yet chosen (**deliberately parked, not
  decided**): (1) a short damped relaxation pass before the real, measured
  evolution starts, reusing the existing `Phase::Relaxing` machinery and
  consistent with the project's own IC-precision philosophy (exact IC
  generation doesn't need to be right, only late-time/attractor behaviour
  does); (2) push the pair/ring further from the affected periodic
  boundary and explicitly quantify the residual -- known from the
  `straight_string_test` investigation to shrink the *local* artifact but
  not the *total* spurious energy, which plateaus at a floor rather than
  vanishing, so on its own this is not rigorous for anything needing
  quantitative energy/tension precision, only qualitative dynamics; (3) an
  exactly periodic closed-form construction (elliptic/theta-function-based
  multi-vortex solution) -- rigorous but a substantially bigger
  implementation effort. **Do not trust any energy/tension number, or
  fine dynamical detail, from either loop IC until this is resolved**;
  qualitative large-scale behaviour away from the affected boundary (e.g.
  `circular_loop`'s initial collapse trend) is less affected but still not
  fully quantified.

- **Confirmed (2026-09-19, with the user, ahead of moving to cluster-scale
  runs): a wave whose wavelength is well resolved on a fine level but
  poorly resolved on the coarser level it borders loses a large fraction
  of its energy crossing that coarse-fine boundary -- genuinely
  dissipated, not mostly reflected.** Requested as a small, standalone
  test before trusting AMR for the real network runs: "does it bounce
  back or is the energy lost from the system?" New `ic_mode =
  wave_packet_test` (`AxionStringsLevel.cpp`), a localised Gaussian
  wave packet in the phase (`psi2`) direction -- exactly massless/
  dispersionless in the linearised theory around the flat-space vacuum
  (group velocity 1 for every `k`, no curvature term), so any shape
  change, reflection or energy loss seen numerically is a genuine
  discretisation effect, not physics. Run on a **static, non-regridding**
  two-level grid (`amr.initial_grid_file` for a fixed level-1 box,
  `amr.regrid_int = -1` to disable all further regridding) so a single,
  motionless coarse-fine boundary could be studied cleanly -- `dx_fine=
  0.05` (8 points/wavelength, well resolved) one side, `dx_coarse=0.1`
  (4 points/wavelength, poorly resolved) the other, `wavelength=0.4`,
  `params_wave_boundary_test.txt`.
  **Result**: the packet crosses the boundary essentially undistorted in
  shape up to the interface, then over a short window (`t~4.4-8.4`, right
  where the packet actually overlaps the boundary) the domain's total
  energy (the existing composite `EnergyKernel.hpp` diagnostic, exact,
  not approximated) drops sharply from `0.1108` to `0.0813` -- **27.6%
  gone by `t=20`**, most of it in that one crossing window, with only a
  slow, much smaller (~1-1.5%) continued decline afterward. A 1D `psi2`
  profile along the propagation axis (`analysis/wave_boundary_test/
  wave_profiles.png`) shows a clear, coherent **transmitted** packet
  continuing into the coarse region at the correct speed (though visibly
  degrading further as it continues to propagate there), and only a
  *tiny* **reflected** disturbance left behind in the fine region --
  nothing like a specular bounce.
  **Controls, to isolate the mechanism** (`params_wave_control_uniform_
  coarse.txt`/`_fine.txt`, same packet, same duration, single level, no
  boundary at all): a uniform grid at the *coarse* resolution alone
  (same 4 points/wavelength) retains **99.15%** of the energy over the
  same `t=20`; a uniform grid at the *fine* resolution retains **99.96%**
  (`analysis/wave_boundary_test/wave_energy_comparison.png`). So poor
  resolution on its own costs ~1%, not ~28% -- the large loss is
  specifically an artefact of the coarse-fine *interface* (its
  interpolation/ghost-cell exchange cannot represent content the coarse
  side has no basis functions for, and that content is mostly discarded
  rather than reflected or conservatively transferred), not simply a
  consequence of running part of the domain at coarse resolution.
  A rough, illustrative region-by-region energy split (`analysis/
  wave_boundary_test/wave_energy_budget.png`, yt's own 2nd-order
  gradients -- qualitative only, not the code's exact 4th-order value)
  is consistent with the above: ~64% of the original energy ends up in
  the coherent transmitted packet, <1% in the reflected remnant, and the
  rest genuinely unaccounted for in either region.
  **Implication for the real network runs**: this is a real, sizeable
  systematic specifically tied to the coarse-fine boundary itself, for
  any field content marginally resolved on the coarser side of a
  refinement jump -- directly the kind of effect conventions.md sec.11's
  "level-timing systematic" section already flags as "the central AMR
  risk," now with a concrete, quantified magnitude behind it (~28% of
  the energy in a marginally-resolved mode, in one boundary crossing).
  Not yet investigated: whether the buffer-width/regrid-frequency inputs
  added the same day (`AxionStringsParams.hpp`'s `apply_regrid_buffer_
  policy`) reduce this (a wider buffer keeps genuinely under-resolved
  content farther from any coarse-fine boundary at the moment of
  regridding, but does not by itself change how the boundary treats
  content that does reach it), nor whether string-network energy content
  near the tagging threshold is typically well- or poorly-resolved on the
  coarser side in practice. Flagged here rather than acted on immediately
  -- a real, now-quantified risk to weigh before relying on deep AMR
  hierarchies for precision energy/spectrum measurements close to a
  level boundary.

  **Follow-up (2026-09-20, with the user): a wavelength ladder across
  three resolution regimes and two successive boundaries, plus spectra.**
  Requested after the single-boundary result above, since 4 points/
  wavelength (that test's "coarse" side) turned out to still basically
  represent the wave, not fail outright: "wavepackets dominated by
  different wavelengths (fits comfortably in the coarser grid, marginal
  on the coarser grid, solidly too small ... but fits on the finest
  grid)". Built a **static three-level** grid (`params_wave_boundary_
  3level.txt`, `wave_test_grids_3level.txt`), levels placed away from the
  periodic seam (first attempt put the finest level flush against it,
  producing a NaN within 2 steps once its periodically-wrapped ghost
  cells needed level-1 data from a region only level 0 covers -- an
  8-unit level-1 buffer on each side fixed it): finest `dx=0.025` |
  intermediate `dx=0.05` | coarse `dx=0.1`, boundaries at `x=32` and
  `x=40`. Same packet construction as above, three wavelengths chosen so
  points/wavelength *at the coarsest level* is 8 ("comfortable"), 4
  ("marginal", matching the single-boundary test above), or 2 ("severe",
  the Nyquist limit) -- giving 2x and 4x that at the intermediate/finest
  levels automatically (fixed ref_ratio=2), so each case crosses a
  *different pair* of points-per-wavelength transitions: comfortable
  32->16->8, marginal 16->8->4, severe 8->4->2.
  **Energy result, and it cross-validates cleanly**: plotting the exact
  moment of each crossing shows the loss depends on the *local*
  points/wavelength transition there, not on which case or which
  boundary -- confirmed by it recurring at consistent values across
  independent occurrences: a 16->8 crossing costs ~2.6-2.8% (comfortable's
  2nd boundary: 99.8%->97.2%; marginal's 1st: 99.8%->97.4%); an 8->4
  crossing costs ~27-28% (marginal's 2nd: 97.4%->70.0%; severe's 1st:
  100%->72.1%, matching the single-boundary test's 27.6% closely). The
  one genuine surprise: severe's 2nd crossing (4->2, all the way to
  Nyquist) cost only a further ~1.6% (72.1%->70.5%) -- far less than a
  second 8->4-sized hit would suggest, i.e. the loss is not monotonically
  worsening per crossing once a wave has already been degraded by an
  earlier one. Not explained -- flagged as observed, not modelled.
  **Spectra** (`analysis/wave_boundary_test/wave3_spectra.png`, psi2
  power spectrum, linear-interpolated onto a uniform finest-resolution
  grid before FFT -- checked this doesn't bias the result by repeating
  with nearest-neighbour resampling instead: same peak location either
  way): comfortable's spectrum is unchanged early to late, essentially
  down to the FFT noise floor at 3x the fundamental (relative power
  `2e-12`, against `2e-27` at t=0 -- i.e. present but utterly negligible).
  Marginal and severe both develop a **real, resolution-induced 3rd-
  harmonic peak** absent at t=0 (`~1e-30` there) and absent in the
  comfortable case: `1.6%` of peak power for marginal, `3.0%` for severe,
  confirmed at the same `k` under both resampling methods (not an
  analysis artefact) and confirmed absent from the initial condition (not
  something already present in the IC). Mechanism not pinned down
  (candidates: genuine cubic-nonlinearity harmonic generation, amplified
  by poor resolution changing the fundamental/harmonic's relative
  survival; or an artefact of the coarse-fine prolongation operator
  itself injecting spurious high-k content) -- flagged as a further,
  distinct symptom of the same under-resolved-boundary problem, not
  further chased. Profile grid (`analysis/wave_boundary_test/wave3_
  profiles.png`) shows the qualitative counterpart: severe's crossing
  visibly leaves the most prominent reflected remnant of the three cases
  (small but clearly a separate, slowly-separating packet, unlike
  marginal or comfortable), consistent with reflection becoming
  proportionally more relevant as resolution worsens, alongside the
  dominant dissipation channel already established above.

  **Major correction (2026-09-20, same day, with the user): most of the
  "energy loss" above is not dissipation at all -- it is substantially a
  measurement artefact of the energy diagnostic's own finite-difference
  gradient stencil, and it is largely reversible.** Requested: normalise
  the spectrum as `(1/rho_0) d(rho)/d(log k)` (so its area under the
  curve reads off the retained energy fraction directly), move the two
  boundaries 1.5x further apart, and run the marginal/severe cases much
  longer (`t=200`, ~2 periodic round trips) to see whether reflected
  energy bounces back and forth. It does more than bounce: **the total
  energy trace does not settle to a lower plateau at all -- it oscillates
  with large amplitude** (marginal: ~66%-95% of the initial value,
  severe: ~64%-95%), repeating over multiple periods with no sign of
  decay toward a fixed asymptote (`analysis/wave_boundary_test/
  wave3v2_energy_full.png`). Splitting the energy into its kinetic
  (`Pi2^2`) and gradient (`|grad(psi2)|^2`) parts separately
  (`wave3v2_kinetic_vs_gradient.png`) shows why: **gradient energy swings
  far more violently (33%-95%) than kinetic (75%-100%)**, and both
  recover substantially whenever the wave packet re-enters a well-
  resolved region (visible directly in the profile grids,
  `wave3v2_profiles_marginal.png`/`_severe.png`, as the packet
  refocusing after a spread-out interval on coarser grids). This matches
  a computable, resolution-dependent property of the project's own
  4th-order first-derivative stencil (`FourthOrderDerivatives::diff1`,
  used for the gradient-energy term): its Fourier transfer function
  `k_eff(k)` is *not* `k` away from the continuum limit, so `(k_eff/k)^2`
  -- the fraction of the true gradient energy the stencil actually
  reports -- is `0.98` at 16 points/wavelength, `0.72` at 4 points/
  wavelength, and exactly `0` at 2 points/wavelength (the stencil's
  transfer function has an exact null at the Nyquist wavenumber,
  independent of phase) -- these numbers bracket the observed retained
  fractions well. **This does not mean there is no genuine loss** (the
  3rd-harmonic generation documented above is real, irreversible mode
  transfer, confirmed absent from the initial condition; the energy
  trace's peaks also decay slowly run-to-run, e.g. marginal's peaks go
  ~95% -> recovering to ~89% by `t=200`, not fully back to 100%) -- but
  it means **the single-boundary result quoted earlier in this entry
  (27.6% "lost") should not be read as 27.6% of the wave's energy being
  destroyed**; a large fraction of that number is the diagnostic
  temporarily failing to see gradient energy that is still physically
  present, reversible once the content returns to good resolution.
  **Practical implication, not yet acted on**: any energy-conservation
  check in this project (including the ones already used for network
  runs) can show large *apparent* non-conservation purely from field
  content moving across a resolution change, with no real energy loss
  behind it -- this needs to be kept in mind before treating an energy
  discrepancy near a refinement boundary as evidence of a bug or of real
  dissipation. Spectra with the requested normalisation
  (`wave3v2_spectra_normalized.png`, four times per case: initial,
  between the two boundaries, a low-energy point, and a recovered/high-
  energy point) show this directly -- the area under each curve tracks
  the (reversible) energy swing, and the harmonic peak's relative
  prominence is visibly smaller at the recovered point than at the low
  point for both marginal and severe, consistent with genuine partial
  recovery, not just an amplitude rescaling of a fixed shape.

  **Further correction (2026-09-20, same day, with the user, pushing
  back): the "measurement artefact" framing above is right for the
  *amplitude* of the energy oscillation, but wrong about *why* the
  severe case's energy comes back -- the user was right to be
  unconvinced, and correctly guessed the real mechanism (reflection, not
  mismeasured transmission).** Checked directly by region-splitting the
  (exact, kinetic+gradient, `yt`-integrated) energy into West (`x<48`,
  the original/finest side), Middle (`48<x<60`) and East (`x>60`, the
  coarse "transmitted" side) over the full `t=200` severe run
  (`analysis/wave_boundary_test/severe_region_split.png`): **East never
  exceeds 0.22% of the initial total energy, at any time in the entire
  run.** Essentially none of the wave ever gets past the second boundary
  -- what looked like a "transmitted packet" in the earlier profile plot
  (a large-*amplitude* feature in the coarse region) carries negligible
  real energy; the actual dynamics is the wave crossing boundary 1
  freely (West <-> Middle, both directions, repeatedly) while being
  *almost totally reflected* at boundary 2, trapped oscillating between
  West and Middle for the whole run -- this, not measurement bias in the
  East region, is what the earlier energy-oscillation plot was showing.
  **Mechanism, verified numerically**: severe's boundary-2 crossing is
  exactly 2 points/wavelength -- the discrete Nyquist limit. The
  project's own Laplacian stencil (`FourthOrderDerivatives::diff2`, used
  in `AxionStringsRHS`) gives a discrete dispersion relation whose *group
  velocity* (`d(omega)/dk`, `omega=k_eff` for the massless mode) is a
  generic, stencil-independent zero at the Brillouin-zone boundary
  (`k*dx=pi`) by lattice symmetry -- checked directly for this stencil:
  `d(k_eff)/d(k dx) -> 0` as `k dx -> pi`, vs. a substantial `0.87` (close
  to the continuum value 1) at `k dx = pi/2` (marginal's own coarsest
  crossing, 4 points/wavelength). A wave with zero group velocity on the
  grid it would need to enter simply cannot propagate there -- so it
  reflects, essentially completely, exactly as observed. Marginal's own
  coarsest crossing (4 points/wavelength) is *not* at this cutoff, so it
  should behave differently (real, if reduced, transmission) -- not yet
  checked with the same region-split, flagged as the natural next step
  if this is worth confirming. **Revised bottom line**: the reversible
  "measurement bias" explanation above is real and relevant to *how much*
  the trapped wave's measured energy swings as it moves between West's
  good resolution and Middle's poorer one, but the *dominant* effect,
  at least for a crossing this severe, is a genuine physical cutoff --
  near-total reflection at a resolution jump landing exactly on the
  discrete Nyquist limit -- not a diagnostic seeing-and-not-seeing
  problem. Do not generalise the "it's just measurement bias, real
  energy is fine" reading from this entry without checking the region
  split first.

  **Same-day follow-up, user asking specifically "what happens between
  t=50 and t=100 -- I'd expect energy conserved, so how does it
  change?"**: the "West" bin above (`x<48`) was itself not uniform, and
  lumping it together hid the answer. `x<48` contains, left to right: the
  coarse-left buffer `[0,12)` (level 0 only, *the same `dx=0.1` as the
  main coarse region*), an intermediate corridor `[12,24)`, and the
  finest zone `[24,48)`. Tracing the actual `psi2` profile through this
  window (`analysis/wave_boundary_test/severe_t50_100_detail.png`) shows
  the reflected wave drifting left from the finest zone, crossing `x=24`
  (finest/intermediate, an easy 8->4 points/wavelength transition) into
  the intermediate corridor, then visibly scattering right at **`x=12`
  -- which is *also* a 4->2 points/wavelength transition, the same
  Nyquist cutoff as boundary 2** -- before reflecting back right and
  recrossing into the finest zone by `t=100`. **So there are two Nyquist
  walls, at `x=12` and `x=60`, not one -- the wave is trapped in a
  48-wide cavity between them**, with `x=24` and `x=48` (both ordinary,
  non-Nyquist crossings) transparent in between. The `t=50-100` energy
  dip-and-recovery is exactly the reversible gradient-stencil measurement
  bias from above, operating at this *second*, previously unmarked
  internal boundary as the trapped wave transits the poorly-resolved
  corridor near `x=12` and back -- not a new or additional real loss, just
  the same mechanism recurring somewhere the original three-bin
  West/Middle/East split didn't distinguish. Lesson for any future
  version of this test: when defining "regions" to track energy in a
  multi-level static grid, split by *actual dx*, not by the two headline
  boundaries -- any region spanning more than one resolution internally
  will hide this kind of substructure.

  **Semi-analytic derivation (2026-09-21, with the user, asked to confirm
  the recovery is understood, not just observed)**: derived the reversible
  measurement bias directly from `FourthOrderDerivatives::diff1`'s own
  weights (`weight_far=1/12`, `weight_near=2/3`) rather than fitting it,
  and used it to explain two things the region-split account above didn't:
  why the *uniform*-grid controls (original single-boundary test) showed
  ~99% retention despite sitting at the same "marginal" 4 points/wavelength
  as the boundary test's coarse side, and why the severe run's troughs
  (~62-70%) go well below the naive single-mode prediction for its
  intermediate zone.
  - For a pure mode `cos(kx)`, `diff1` measures an effective wavenumber
    `k_eff` with `k_eff*dx = 2*(weight_near*sin(k dx) - weight_far*sin(2 k
    dx))`; gradient energy is suppressed by the exact, phase-independent
    factor `T(N) = (k_eff/k)^2` for `N = lambda/dx` points per wavelength.
    Verified numerically (apply the real stencil to the actual IC formula,
    compare to the continuum-exact gradient computed by finite-differencing
    a 200x finer sample): `T(8)=0.976`, `T(4)=0.721`, `T(2)=0.000` (exact
    null at Nyquist) -- matches the closed form to machine precision.
  - `Pi2`'s own energy (point samples, no derivative operator) is exactly
    *unbiased* for `N>=3` -- a discrete-orthogonality fact (the sampled
    mean square of `cos(kx+phi)` over any `N>=3` uniformly spaced points
    per period is exactly `1/2` for *any* phase `phi`) confirmed
    numerically to 4 decimal places at `N=4` and `N=8`. **This is why the
    uniform-grid controls looked fine**: `total_energy`'s reported
    percentage is `E(t)/E(0)`, and on a uniform grid both `E(0)` and every
    later `E(t)` carry the *same* resolution's bias, which cancels in the
    ratio -- a uniform grid can never reveal this effect via its own
    energy-conservation trace, only a resolution *change* between the
    numerator and denominator can, which is exactly what the boundary
    tests are.
  - Assuming equipartition (exact for this IC: kinetic and gradient
    energy are equal by construction, `Pi2=-dpsi2/dx` at `t=0`) and
    `Pi2` unbiased, the predicted measured fraction after a pure `N=8->N=4`
    move is `(1+T(4))/(1+T(8)) = 87%` -- reproduced almost exactly by
    `analysis/wave_boundary_test/nyquist_measurement_bias_semianalytic.png`
    (right panel), and the severe run's *plateaus* (95%, 90%, 85%,
    decreasing slightly on each successive bounce) sit close under this
    model's trivial 100% self-consistency check, the residual gap being
    the already-established small genuine loss (harmonic generation) on
    top of the reversible piece.
  - But the *troughs* (~62-70%, and noticeably jittery/noisy rather than
    smooth) go well below the 87% "purely N=4" prediction. Exactly at
    `N=2` (Nyquist), the discrete-orthogonality argument above breaks
    down -- the sampled mean square of `cos(kx+phi)` at `N=2` is
    `cos^2(phi)`, genuinely phase-dependent, so **even the kinetic term
    can nearly vanish** depending on the wave's exact alignment with the
    grid at that instant (checked numerically for the severe packet's
    actual phase at `x=12`/`x=60`: measured kinetic and gradient both
    collapsed to ~1-3% of true). The severe wave doesn't sit purely at
    `N=4`; it grazes these `N=2` Nyquist walls each bounce, so both terms
    can transiently collapse together -- explaining both the extra depth
    of the troughs beyond the single-mode `N=4` estimate and their jitter
    (a phase-alignment effect, not noise in the usual sense).
  - **Answer to "are we sure the energy increase on re-entering the
    finest region is understood": yes** -- it is the same reversible
    stencil bias recalibrating back toward its `N=8` value, not energy
    being created; the analytic model predicts the recovery scale
    correctly, and the leftover few-percent gap plus the progressively
    lower plateaus track the independently-established genuine harmonic
    generation, not a new unexplained effect.

**`axion_strings.pre_evolution.gamma` now defaults, rather than being a
required manual input (2026-09-22, with the user).** Every existing
`ic_mode = fourier_relaxed` parameter file had to compute `1/dx_base` by
hand and enter it as `gamma` -- exactly the bug class the
`params_full_test_1024.txt` entry above documents (a stale/mismatched
`gamma`, caught only because the run's own startup printout was read
carefully). Prompted by the user asking, while discussing that entry,
whether the pre-evolution's resolution should match the main run's
resolution at `tau_i` -- yes, it should, and checking `tag_cells()`
confirmed why the existing convention (`gamma = 1/dx_base`) was already
correct in every validated example: refinement is off entirely during
`Phase::Relaxing`, and the schedule-gated tagger only brings level 1
online once `log(m_r/H)` crosses *its own* threshold, so as long as
`log_mr_over_h_i` sits below that threshold, the main run genuinely
starts with only the coarsest grid present.

New `apply_pre_evolution_gamma_default()` (`AxionStringsParams.hpp`),
called from `apply_box_plan()` in both the AMR and uniform-grid branches
(guarded to `ic_mode == "fourier_relaxed"`, and left untouched for Moore
mode, which already derives/requires `gamma` its own way): derives
`gamma = 1/dx_base` and injects it via `GRParmParse::add` if the user
hasn't set it, or cross-checks an already-set value against that
derivation and aborts on mismatch -- the same derive-if-absent,
cross-check-if-present pattern already used for `amr.n_cell`/
`geometry.prob_extent`/`evolution.stop_time` in the same function, now
extended to this parameter too. Also added a startup check (AMR branch
only) that `log_mr_over_h_i` is actually below the level-1 threshold --
the load-bearing assumption behind "the coarsest level is the only one
present at `tau_i`" -- and aborts with a clear message rather than
silently deriving a `gamma` that under-resolves the string cores relative
to the grid they're about to sit on, if it isn't. Deliberately *not* the
more general "derive from whichever level is actually active at `tau_i`"
version discussed first -- the user asked to assume only the coarsest
level is ever present at the start, matching every existing example, and
the new safety check turns a violation of that assumption into a startup
error instead of a silent one.

**Verified**: all 51 existing unit tests still pass (unaffected --
`BoxPlan.hpp`/`Background.hpp` themselves are untouched, only
`AxionStringsParams.hpp`'s own injection logic changed). Smoke-tested all
four paths directly against the real executable on
`params_amr_validation_128.txt` (AMR) and `params_demo_384_physical.txt`
(uniform grid, `max_level=0`): `gamma` omitted derives to the exact
existing hand-set value in both (`5.656854249` vs `5.65685`;
`19.59591794` vs `19.596`); `gamma` already correct runs with no error;
`gamma` deliberately wrong aborts with the cross-check message; and
`log_mr_over_h_i` pushed past the level-1 threshold (with `gamma` unset)
aborts with the new schedule-violation message, rather than silently
under-resolving the relaxed field.

**Follow-up, same day: the mismatch check is a warning, not an abort
(2026-09-22, with the user, asked directly "can gamma still be
overridden by the user if they want?").** The original cross-check
aborted on any mismatch, which technically blocked a *deliberate* choice
to relax at some other resolution, not just an accidental stale value --
the user asked for exactly this to be possible, so
`apply_pre_evolution_gamma_default` now calls `GRParmParse::warning`
(GRTeclyn's existing, already-wired-up non-fatal diagnostic --
`SetupFunctions.hpp`'s `mainSetup` already surfaces
`GRParmParse::warnings_issued()` at startup) instead of `.error()` on a
mismatch, printing the same explanation but letting the run proceed with
the user's value. The `log_mr_over_h_i`-vs-level-1-threshold safety check
added above is unchanged and still aborts -- that one guards an assumption
the *default's derivation itself* depends on (nothing to override there,
since violating it doesn't correspond to a valid alternative choice, just
a wrong default). Re-verified: rebuild clean, all 51 unit tests still
pass, and the deliberately-wrong-`gamma` case now prints
`Warning from parameter axion_strings.pre_evolution.gamma = 3.0: ...`
(matching the exact format of the codebase's other parameter warnings,
e.g. `evolution.stop_time`'s) and completes its steps normally rather
than aborting.

**Multi-threading (`USE_OMP`) enabled, after the user asked whether extra
`--cpus-per-task` on a SLURM job actually get used (2026-09-22).**

Audited every parallel loop in `AxionStrings/` (not just spot-checked)
before touching anything, since a naive "just turn OpenMP on" is exactly
the kind of change that can silently corrupt output via a data race
without ever crashing. Two things made this audit tractable:
- `grep -rl MFIter *.hpp *.cpp` found only two files with a hand-written
  MFIter loop (`ProjectionKernel.hpp`, `SpectrumKernel.hpp`) -- every
  other per-cell kernel in the project (`specific_eval_rhs`, the tagger,
  every `initData` IC branch, the pre-evolution->main rescale) already
  goes through `amrex::ParallelFor(some_multifab, ...)`, and every
  reduction (`EnergyKernel.hpp`, `VelocityKernel.hpp`, `StringFinder.hpp`)
  through `amrex::ReduceOps`.
- Read AMReX's own CPU implementation of both
  (`amrex/Src/Base/AMReX_MFParallelForC.H`'s `ParallelFor_doit`,
  `amrex/Src/Base/AMReX_Reduce.H`) rather than assuming: both already
  wrap their internal `MFIter` loop in `#ifdef AMREX_USE_OMP #pragma omp
  parallel #endif`, tiled per-box. So every one of those call sites was
  *already* thread-safe and gets real multi-threading for free the moment
  `USE_OMP=TRUE` is set at build time -- no AxionStrings-specific kernel
  code needed to change for the main per-step cost (the RHS evaluation)
  or any of the routine diagnostics.

The two hand-written exceptions were a genuine hazard, not a formality:
`ProjectionKernel.hpp`'s `compute_energy_projection` and
`SpectrumKernel.hpp`'s `compute_spectrum` both accumulate into shared
arrays/scalars indexed by something *other* than which box a thread is
processing (`(i,j)` column and k-shell respectively) -- multiple boxes on
one rank routinely land in the same output slot (a very ordinary
consequence of AMReX's own domain decomposition), so wrapping either
loop's existing `#pragma omp parallel` would race on the max-update,
`++`, and `+=` accumulations. Fixing that properly needs per-thread
partial buffers merged afterward; both files' own header comments already
document them as running "once per output snapshot, not every substep",
i.e. not the cost `USE_OMP` exists for, so the right call was to leave
both loops deliberately single-threaded and document *why* directly next
to the loop (not just in this entry, where it would be easy to miss),
rather than either race silently or spend the complexity budget on a path
that was never the bottleneck.

`AxionStrings/GNUmakefile` now documents `USE_OMP` (default `FALSE`,
matching GRTeclyn's own default and every executable built so far) as an
explicit, available `make USE_OMP=TRUE ...` flag, with the reasoning
above summarised right there for whoever next needs it.

**Verified**: all 51 unit tests pass and the existing non-OMP build
(`COMP=llvm`, this project's own Mac toolchain) still compiles and links
cleanly after the `ProjectionKernel.hpp`/`SpectrumKernel.hpp`/
`GNUmakefile` edits (comments only for the first two -- no behaviour
change there either way).

**Not verified locally, flagged honestly rather than glossed over**:
actually compiling and running a `USE_OMP=TRUE` build. Six attempts
across three toolchains on this Mac all failed for reasons specific to
*this machine*, not the code: Apple's system clang (what `mpicxx` wraps
by default here) doesn't support `-fopenmp` at all; Homebrew's LLVM
supports it but its libc++ `<math.h>` conflicts with the Apple SDK
(`FP_INFINITE`/`FP_NORMAL` undeclared) -- the exact, already-documented
issue this file's own toolchain note works around by putting `/usr/bin`
first, which is precisely what removes OpenMP support; and Homebrew's GCC
16 hits its own SDK header search-path issue (`wchar.h`/`stdlib.h` not
found) that a couple of quick `CPATH`/sysroot attempts didn't resolve.
None of these are expected to occur on a real Linux cluster (`COMP=gnu`
against a native system GCC, which is what the cluster guide already
recommends) -- but this means the `USE_OMP=TRUE` path has only been
verified by careful reading of AMReX's own source, not by an actual
multi-thread run compared against a single-thread one on this codebase.
**First thing to do on a real cluster**: build with `USE_OMP=TRUE`, run a
short test at `OMP_NUM_THREADS=1` and again at `OMP_NUM_THREADS=4` (or
similar) on the same seed, and diff `network_scalars.dat` -- if the two
aren't identical (or at least statistically indistinguishable for a
`fourier_relaxed` run with genuine floating-point reduction-order
sensitivity), something in this audit missed a case.

**Follow-up, same day: a worked cluster-scale example parameter file
(2026-09-22, with the user, asked to make sure every option discussed in
the cluster guide -- `axion_strings.tagging.regrid_interval_steps`
specifically named -- is genuinely settable, plus a clear example for a
"moderately big cluster run").**

New `AxionStrings/params_cluster_512base_2level.txt`: `N_base=512`,
`amr.max_level=2` (`N_effective=2048`) -- a deliberate step up from every
base grid tried so far in this project (previously largest was
`N_base=160`, `params_full_test_640.txt`). Uses the
`regrid_interval_steps`/`buffer_safety_factor` mechanism (not hand-set
`amr.n_error_buf`/`regrid_int`) with `buffer_safety_factor=1.5` for the
extra margin a not-yet-tried configuration warrants, leaves
`pre_evolution.gamma` unset to demonstrate the new auto-derivation, and
turns on every optional diagnostic (masking scheme B, spectrum,
projection) so the file doubles as a complete worked example of every
category of setting the cluster guide documents.

Re-verified every parameter the guide's parameter tables reference
against its actual `pp.get`/`pp.queryAdd` call site (not just the ones
touched recently) -- all confirmed genuinely wired, including the native
AMReX-read ones (`amr.blocking_factor`, `regrid_int`, `n_error_buf`,
`max_grid_size`, all read directly by `amrex::AmrMesh`/`amrex::Amr`, not
by this project's own code) and `evolution.dt_multiplier` (read
unconditionally by `GRAmrLevel::ComputeDt`, not only when the regrid-
buffer-policy mechanism happens to also read it).

**Verified**: smoke-tested the new file directly (not just written and
assumed correct) -- `evolution.max_steps=0` first, confirming every
box-plan number in the file's own header comment (`L_tilde`, `dx_base`,
`dx_finest`, both level thresholds, the auto-derived `gamma`) matches the
program's actual startup printout exactly; then 3 real steps at 8 MPI
ranks (the $512^3$ level-0 grid, matching the earlier fixed-grid speed
benchmark's cell count almost exactly, so ~45s/step here was expected,
not a red flag) with no NaN or crash. The `axion_strings.save_projection`
"unused ParmParse variable" warning that showed up in that 3-step run is
expected, not a bug: `ic_mode=fourier_relaxed` starts in
`Phase::Relaxing`, and the diagnostic-output code path that reads
`save_projection` is gated to `Phase::Evolving`, which 3 steps of
relaxation is nowhere near reaching -- already independently confirmed by
reading the call site directly, not just inferred from this run's
behaviour.

- **Investigated and fixed (2026-09-23, with the user, following up on
  the user's own prior-work observation: "strings are produced with
  excited core modes that lead to strong oscillations in the spectrum
  close to mr and obscure the physics"): the pre-evolution -> main handoff
  excites a real, measured core-breathing transient, traced to a
  discontinuous jump in `lambda`/`curvature_term_coeff` at the switch --
  now smoothed by default.**

Direct simulation (a small, fast physical-string probe, `N=256`, single
level, `fourier_relaxed` IC, fine output cadence spanning the handoff)
found two independent, converging signatures right after
`apply_pre_evolution_to_main_rescale()`: (1) `rho_radial_kin_unscreened`
and `rho_radial_mass_unscreened` in `network_scalars.dat` rise and fall
out of phase with each other (kinetic peaking at `tau=2.93`, mass peaking
later at `tau=3.49`) -- the phase lag of an underdamped oscillator
exchanging kinetic/potential energy, not a monotonically-relaxing
profile; (2) `axion_spectrum.dat`'s unscreened spectrum shows a bump that
starts near mode index `p~16` and migrates down towards `p~7-9` over the
probed window, still not settled. Both survive masking (the *screened*
columns show the same shape, just damped in amplitude).

**Root cause**: for a physical run (`c0=0`), `Background::lambda(tau)=1`
identically (constant comoving core mass), but `PreEvolutionBackground::
lambda(tau_pre)=gamma_pre^2/R_pre(tau_pre)^2` *decreases* through
relaxation -- using this probe's own numbers (`gamma_pre=16`, auto-
derived; `R_pre(tau_pre_end)=4.158`), the comoving core mass jumps from
`3.85` just before handoff to `1.0` just after: a discontinuous ~4x
change in the string core's natural width, landing at the exact instant
`apply_pre_evolution_to_main_rescale()` runs. That rescale is an exact
kinematic identity for `psi`/`Pi` (the `psi=R phi` chain rule) but has no
mechanism to fix up `lambda` itself, which is what actually sets the
core's static profile -- a core equilibrated for the pre-evolution mass
rings when suddenly sitting in a background that wants a different one.
`curvature_term_coeff` jumps too, for an unrelated reason (pre-evolution's
own schedule has it constant at `1.0`; the main schedule's own formula,
`(1-b_inv)/(b_inv^2 tau^2)`, is identically `0` for `a_inv=2`).

**Fix**: `axion_strings.pre_evolution.handoff_transition_n_periods`, read
in `AxionStringsParams::read_handoff_transition_n_periods()`, freezes
`lambda`/`curvature_term_coeff` at their pre-evolution values the instant
of handoff (`apply_pre_evolution_to_main_rescale()`) and blends them to
the main schedule's own values via a smoothstep (`3t^2-2t^3` -- value-
continuous and zero-slope at both ends, so the smoothing itself adds no
new kick) in `specific_eval_rhs()`'s main branch, over this many main-
schedule core-oscillation periods (`2*pi/sqrt(lambda_main(tau_i))` -- a
physical, resolution-independent timescale, not a bare `tau` window).

**Verified** (same seed, bit-identical relaxation trajectory up to
handoff in every comparison, since nothing before the handoff changed):
at `n=3`, `rho_radial_mass` peak dropped from `0.093` to `0.016` (~6x),
the spectrum bump peak from `1.08e14` to `3.27e13` (~3.3x), and the
`xi(tau)` trajectory changed from a spurious overshoot (`0.98 -> 1.43`)
to a smooth monotonic decline (`0.98 -> 0.79`) -- `xi`'s overshoot was
never separately targeted, so its disappearance is independent evidence
the mechanism, not just its two originally-flagged symptoms, was
correctly diagnosed. Doubling the window to `n=6` gave **no further
improvement** (rad_mass/rad_kin differ from `n=3` by <0.2% throughout the
probed window) -- `n=3` already captures essentially all of the benefit
smoothing alone can buy; the residual that's left (still decaying, not
zero) is most likely genuine pre-existing formation/annihilation
radiation baked into the field before handoff even happens, not a
handoff-smoothness artifact, and would need a different lever (more
relaxation time before handoff, or just trusting diagnostics only after
some settling margin post-handoff) if it needs reducing further. A
follow-up test of that different lever -- doubling `gamma_pre` to buy
more relaxation time -- was tried and made things **worse**, not better:
`gamma_pre` also directly inflates the size of the jump being smoothed
(`lambda_start` rose from `14.8` to `24.1` despite the longer relaxation),
so `rad_kin` peak rose ~39% even with `n=3` on top; this is a confounded
knob (changes both the mismatch size and the relaxation time at once)
and was reverted -- **the default `gamma_pre` (auto-derived, unset in
every example file) is kept**.

**Decision**: `handoff_transition_n_periods` now **defaults to `3.0`**
(previously would have defaulted to `0`, i.e. opt-in) -- every
`fourier_relaxed` config that doesn't set this explicitly gets the
smoothed handoff automatically; set to `0` to recover the old
instantaneous-jump behaviour. `AxionStrings/params_cluster_512base_
2level.txt` and `docs/cluster_guide/cluster_getting_started.tex`/`.pdf`
updated to document the new default (Section~5.3 "The handoff from
relaxation to the main run" in the guide). 54/54 unit tests passing
throughout (pure `AxionStringsLevel.cpp`/`AxionStringsParams.hpp` logic,
not exercised by the standalone doctest suite).

**Not yet pursued, considered and set aside for now**: an alternative
design that avoids the ξ-based stopping criterion entirely (stop
relaxation at the analytically-solved `tau_pre_end` where `lambda_pre`
exactly crosses `lambda_main(tau_i)`, a closed form needing no empirical
fit, making `lambda` exactly continuous by construction rather than
smoothed after the fact) was discussed with the user but not implemented
-- the smoothstep fix above was judged sufficient once verified. Revisit
if a future config needs `xi` at handoff decoupled from what the smoothed-
jump default happens to produce.

- **Added and tested (2026-09-24, with the user): `axion_strings.tagging.
  force_full_refinement`, a refinement-systematics control -- does AMR
  itself bias any observable, as opposed to just trading cost for
  resolution?**

Prompted by an earlier finding the same day (`N_p` jumping discontinuously
right at a level's onset in an AMR run, see the `tension_core_only` drift
entry above) -- the user wanted a direct way to test whether refinement
*itself* biases measured observables, not just resolution. New `bool
StringTaggerParams::force_full_refinement` (default `false`,
`axion_strings.tagging.force_full_refinement`): when set, `tag_cells()`
tags every cell at every level unconditionally, bypassing both the
string-based criteria and the schedule-gating that normally only permits
a level once `log(m_r/H)` crosses its own threshold. Placed in
`AxionStringsLevel::tag_cells()` right after the existing `Phase::
Relaxing` early return (relaxation never refines regardless of this flag,
so both an ordinary and a fully-refined run share bit-identical relaxed
ICs) and before the schedule-gating block. Uses the same per-cell
`amrex::TagBox::SET` pattern `StringTagger`'s own kernel already uses,
not a new API. Documented in `params_cluster_512base_2level.txt`'s
reference block.

**Verified the mechanism directly** (`N=128`, `max_level=2`, physical,
otherwise matching the validated `params_tagger_smoke_test.txt` config):
confirmed via the run's own regrid trace that with the flag set, level 1
*and* level 2 both come online at "100% of domain" within 2-3 regrid
cycles of the main run starting (`TIME=3.5355` at handoff; level 1 by
`TIME=3.606`, level 2 by `TIME=3.606` too) -- not literally the first
step, since AMReX's regrid can only add one level per cycle, but
effectively immediate. The ordinarily-tagged run, same seed, instead
brings level 1 online later (`TIME=4.808`) and only partially (~40-43%
of the domain), exactly as expected from the string-based criteria.

**Then ran the actual A/B comparison this feature exists for**: same
seed/IC, `k_max_over_mr=16`, `xi_target=1.0`, masking scheme B, one run
with ordinary adaptive tagging and one with `force_full_refinement=1`.
Compared `xi` and `tension_core_only` across `tau=1.72` to `11.05` --
spanning well before, during, and after the adaptive run's own level-1
(`tau~4.8`) and level-2 onsets. **Result**: `xi` and `tension_core_only`
agree to within ~0-3% at every single snapshot, with no discontinuity or
growing divergence at either refinement transition; from `tau~6.1`
onward the two runs' `N_p` become numerically identical and the
differences drop to exactly 0.0% (the adaptive tagger has independently
reached ~100% coverage by then too, at this small scale/short duration).
`N_p` itself (the raw plaquette count, not length-normalised) differs by
up to ~4x throughout the early/mid range -- expected and not a concern:
the fully-refined run has far more fine cells everywhere, so it detects
more raw windings, but `xi`/`tension_core_only` correctly length-weight
per level (`ell_comoving = (2/3) * Sum_l N_p_l * dx_l`), so the *physical*
quantities converge even though the raw count doesn't.

**Caveats, since this was a first, quick check, not a rigorous
validation**: small grid (`N=128`), short duration, one seed. The
adaptive run's own tagged fraction converges to ~100% (matching
`force_full_refinement` exactly by construction) well before the run
ends at this scale, so the most informative window is really the earlier
part (`tau<6`) where adaptive coverage is genuinely partial (~40%) yet
still tracks the fully-refined reference to within a few percent -- a
bigger box, run for longer, would keep that partial-coverage regime
informative for a larger fraction of the run, and multiple seeds would
give a real error bar on the residual differences rather than a single
noisy trace. Not done as part of this first pass. Not a production
setting either way: refining everything discards AMR's entire cost
advantage, so `force_full_refinement` should only ever be used for this
kind of systematics check, never for an ensemble run.

- **Bug found and fixed (2026-09-25, with the user): `axion_strings.masking.scheme = A` was a silent
  no-op for every energy diagnostic -- `screened` and `unscreened` came out numerically identical.**

Found while writing the network-evolution physics reference PDF: the user asked directly whether
scheme A screens energies by multiplying by `|psi|/R`, and checking `Masking.hpp::masking_weight()`
against that showed it returned `1.0` unconditionally for anything other than scheme B --
`rho_tot_screened`, both axion energies, and all three radial energies were therefore identical to
their own unscreened values under scheme A, for every run that has ever used it. This did not affect
scheme B (the default, and the only scheme `compute_spectrum` allows) or scheme None.

**Root cause and fix**: `masking_weight()` only ever implemented scheme B's hard top-hat; scheme A
fell through to the same `return 1.0` as None, rather than the smooth weight `initialMD/conventions.md`
actually specifies for it, `f = (1 + r/f_a)^2` (`r = |phi| - v` the radial deviation from vacuum, the
same `r` `Energy.hpp`'s `radial_*_energy_pointwise` functions use; `f_a = sqrt(2) v`). Implemented
exactly as specified -- in terms of the already-computed `mod = |psi|/R`, `r = mod - 1` -- giving `w`
that varies smoothly from `~0.086` at the exact core (`|phi|=0`) up to `1` in the far field, no hard
cutoff. `masked_a_dot()`'s own scheme-A branch (the bare numerator, no division) was already correct
and is conventions.md's *other*, independently-stated form of scheme A -- specific to the axion-kinetic
spectral quantity, algebraically distinct from `w` above (a `|psi|^2`-type weighting on that one
quantity, not literally `w * theta_prime`) -- left unchanged.

**Verified**: `tests/test_masking.cpp`'s existing T2a case asserted `masking_weight(scheme=A, ...) ==
1.0` for a point deep inside a core -- itself a symptom of the bug, now replaced with a check against
the correct formula (54/54 unit tests passing). Confirmed live in a real `fourier_relaxed` run
(`masking.scheme=A`): `n_unmasked` is now a genuinely continuous value (`3835.89` out of `n_total=
4096`, not an integer/exact match) rather than always equal to `n_total`, and screened energies now
visibly differ from unscreened (e.g. `rho_axion_kin`: unscreened `1.172`, screened `0.429`, at one
snapshot) rather than being identical.

**Practical impact**: since `compute_spectrum` requires scheme B, and every validated example file in
this project uses either B (with spectrum on) or leaves masking unset (default `A` in `MaskingParams`,
but never exercised for energies before this fix since nothing reads `rho_*_screened` under scheme A in
any checked-in analysis unless a user explicitly ran with `masking.scheme=A` and looked at the
screened energy columns specifically) -- no checked-in validated result is known to depend on the buggy
behaviour, but any past ad hoc run that *did* set `masking.scheme=A` to inspect screened energies would
have silently gotten the unscreened values back. `docs/physics_reference/network_evolution_physics.tex`'s
masking section was rewritten at the same time to describe both mechanisms (the energy weight `w` and
the spectrum-only `masked_a_dot` construction) separately and correctly.
