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
| 1.8 Energies | Mostly done -- core-energy extraction validated; the log(L) tail growth is a known, deferred rough edge | `Energy.hpp`/`EnergyKernel.hpp` implement `rho_tot`, axion kinetic and axion gradient energy, each screened/unscreened via the sec.7 masking module (`MaskingScheme::None` gives the unscreened pass "for free"). **`rho_tot`'s normalisation was wrong and has been corrected** (2026-09-17, with the user): conventions.md sec.10's formula, read literally, has kinetic/gradient terms scaling as `R^-2`; re-deriving from scratch (a canonically-normalised *complex* scalar has no 1/2 on kinetic/gradient -- checked by matching `L=A\|phi_dot\|^2-...` against sec.3's given EOM, which forces `A=1`) gives `R^-4` instead, verified three independent ways. The potential term is unaffected (matches either way once expanded). Full derivation in `Energy.hpp`'s header comment. Found and fixed a real T1 IC bug the user correctly flagged as "quite a bad sign": `Pi=0` everywhere let the far field's `\|psi\|/R` decay as `R_i/R(tau)` regardless of any string, collapsing the screened energy to 0 within 2-3 steps (confirmed numerically: `R_i/R(1.3)=0.769<0.8`, exactly matching when it broke) -- fixed via `Pi=(R'/R)*psi`, matching the homogeneous solution's background-tracking. Also swapped the radial profile from `tanh` to `rho_hat/sqrt(rho_hat^2+2)` for the correct `1-rho_hat^-2` far-field falloff (tanh decays exponentially, missing the long-range Goldstone tail responsible for the log tension divergence). A control (exact homogeneous solution stays at exactly zero energy indefinitely) confirmed neither bug was a deeper pipeline problem. **Also reworked how tension is extracted** (2026-09-17, with the user): the naive screened/unscreened *difference* of averages is diluted by the masked point-count fraction and is not itself "the string energy" -- switched to explicitly summing the core contribution (`n_total*unscreened - n_unmasked*screened`), which **converges cleanly to a constant (~3.72) independent of L for L>=32** (checked L=16..256 at `N_2=1`, per the user's guidance that `m_r*dx=1` is their production convention, not the far-finer resolutions tried initially) -- exactly the expected behaviour for energy localised at the core. All raw components (screened+unscreened rho_tot/axion_kinetic/axion_gradient, plus `n_total`/`n_unmasked`) are saved to `network_scalars.dat` every snapshot rather than baking in one derived formula, per the user: any combination (core-only, core+tail, or others) can be reconstructed afterwards. A second, "core+tail" measure (core plus the away-from-core `axion_gradient - axion_kinetic`, intended to cancel a propagating axion wave's equal kinetic/gradient contribution and isolate a static string's long-range tail) is implemented as a convenience diagnostic but its L-growth is **not** clean log(L) (closer to `L^0.6`); found that the axion kinetic term is exactly zero in this static test (`Pi` is always parallel to `psi` by construction, so `theta'=0` identically), so this measure is really just "core + total unmasked gradient energy" with no wave actually being subtracted -- the discrepancy is suspected to be either the two-vortex product-profile's amplitude between the cores not matching a true relaxed dipole (which the gradient identity is sensitive to via its `1/\|psi\|^2` factor), or `N_2=1` being too coarse right at the cores where `\|grad(theta)\|^2` is largest. Left as a known, explicitly deferred rough edge (user: "keep this in the back of our minds, but press on") rather than pursued further now. |
| 1.9 Spectra | Core machinery done (T3 passes); physical-unit rescaling and radial-mode spectra not yet wired up | `SpectrumKernel.hpp`'s `compute_spectrum` builds on `amrex::FFT::R2C`. Verified its round-trip normalisation matches conventions.md sec.9's stated convention exactly ("no normalisation, forward+backward multiplies by N^3"), so AMReX's raw `forward()` output can be used directly as `X~_p`, no correction factor needed (unlike `FourierIC.hpp`'s generation direction). R2C stores `kx>=0` only (Hermitian symmetry); every `kx` plane except `0` and `Nx/2` has an unstored conjugate partner at `Nx-kx` and is counted twice, to recover sums over the *full* spectrum. Binning follows sec.12 exactly (`modeIndex=floor(p_mod+0.5)`, `4*pi*p_mod^2\|X_p\|^2` per shell, dropped outside the inscribed sphere). `MaskedFieldBuffer.hpp` applies masking at exactly the one place CLAUDE.md constraint 5 asks for -- filling the buffer handed to the FFT, reusing `Masking.hpp` unchanged. **Resolved an `N^6` question with the user**: task 1.8's real-space/spectral cross-check formula (`(1/2N^6)*Sum\|a_dot(p)\|^2`) looked like it should have `N^3` (standard Parseval for this convention); redoing it carefully showed `N^6` is correct -- a spatial *average* already divides by `N^3` once, and Parseval divides by `N^3` again, so `<a_dot^2> = (1/N^6)*Sum_p\|a_dot(p)\|^2` exactly. T3 (a known-amplitude, known-mode plane wave, checked immediately after `initData` rather than after any evolution step, since `R(tau)` and hence the field itself moves on after even one step) **passes exactly**: real-space `<a_dot^2>`, Parseval from the full-cube sum, and Parseval from the inscribed-sphere sum all equal the hand-computed expected value `0.5*C^2` to machine precision, at both a low mode (`p0=4`) and one that aliases down from outside the naive Nyquist range (`p0=20` on a 32-point grid aliases to `p0'=12`, a correct consequence of `kx` being capped at `Nx/2` by the sampling theorem, not a bug). The shell-binned integral (`Sum_s S(s)/N_total^2`) is only *approximate* even after the same `N_total^2` normalisation (caught this: it was missing from an early diagnostic print) -- a single anisotropic mode shares its shell with other same-radius lattice points that carry zero power, diluting the angular average; unlike the full-cube/inscribed-sphere sums (which include every mode exactly once, so match the real-space value exactly), this is inherently a cross-check, not an identity. Not yet done: rescaling the raw `(p, S(p))` output to physical `(k/H, v^-3 d(rho)/dk)` units (sec.10/sec.12's `(2pi/(LH), (1/2pi) R L_tilde/N^6)` factors -- not yet re-derived/verified the way `rho_tot` and this cross-check were); radial-mode (`r`) spectra (sec.10: "follow identically" from the axion case, not yet implemented); wiring the axion spectrum into the routine per-snapshot diagnostics (currently opt-in via `axion_strings.compute_spectrum`, used so far only for the T3 IC-mode check and as an ad hoc `specific_post_timestep` addition). |
| 1.10 Remaining diagnostics | Velocities done (per user: "implement only the string velocities not the loops or curvatures"); curvature and loop-finding explicitly out of scope for now | `Velocity.hpp` implements conventions.md sec.8's `gamma^2 v^2` estimator (global-string profile coefficients `c1=0.41222`, `e1=-0.025763`), evaluated at the corners of every pierced plaquette and averaged over the network. Sec.8's literal formula is written in terms of `psi`/`Pi` and a symbol `beta` that is not otherwise defined in the conventions snapshot available here; read as `beta = 1/b_inv` by analogy with the identical `(1/(b_inv tau))*psi` term in the main EOM/`rho_tot`, and **verified rather than assumed**: converting to the physical field `phi=psi/v` and its cosmic-time derivative makes every factor of `R` in the literal formula cancel exactly (shown algebraically in the header comment, and confirmed to 1e-9 numerically in `tests/test_velocity.cpp` across 6 diverse `(R,tau,b_inv)` points) -- a physical velocity should not depend on the comoving rescaling used to evolve the field, so the clean cancellation is itself a strong consistency check on the `beta` reading. The production kernel (`VelocityKernel.hpp`'s `compute_velocity_at_pierced_corners`) works directly in the R-independent `phi`/`phi_dot` form (reusing `StringFinder.hpp`'s exact 7-corner/3-plaquette structure) rather than recomputing the cancellation on every cell; the literal `psi`/`Pi` form is kept only as the cross-check. Wired into `specific_post_timestep` alongside the existing xi/energy/spectrum diagnostics (`mean_gamma_sq_v_sq`, `mean_gamma`, `n_velocity_corners` -- printed and saved to `network_scalars.dat`); the zero-strings case (`N_corners=0`) is handled explicitly rather than dividing by zero. Smoke-tested against both `straight_string_test` (T1: small, finite, non-negative `<gamma^2 v^2>` throughout a run of a near-static string, as expected) and the default homogeneous IC (no strings: exactly `<gamma^2 v^2>=0`, `<gamma>=1`, `N_corners=0`). Not implemented, per explicit user instruction to defer as "more tricky": curvature and loop-finding (T4). |
| Output infrastructure (cadence, restart-safety, persistence) | **Done** | Follow-up to task 1.10, prompted by the user asking what output infrastructure existed and how it should work for a multi-run ensemble (2026-09-18). Three changes: (1) **Full-grid outputs off by default** -- `params_test.txt` now sets `amr.plot_int = amr.check_int = -1`; `Mode::PreEvolution` already hard-requires `check_int >= 0` for its handoff checkpoint (`AxionStringsParams::check_params`), so that path is unaffected. (2) **A real output cadence**, replacing "every coarse step": the full per-snapshot diagnostics (plaquette count, energy/velocity reductions, the spectrum FFT) are gated behind a cheap check of `log(m_r/H)` (free from the analytic `Background`, no grid pass) against `s_next_output_log_mr_over_h`, which starts at the new required `axion_strings.output_first_log_mr_over_h` and is bumped by `axion_strings.output_delta_log_mr_over_h` (also required -- no default fixed here, since that is a physics-adjacent numerical choice) each time it is crossed, via a `while` loop so a step that jumps past more than one threshold still lands correctly rather than drifting. (3) **Two real bugs fixed, not just flagged**: `network_scalars.dat` used `SmallDataIO`'s "old" constructor, whose `first_step = (time == dt)` heuristic assumes time starts at 0 -- ours starts at `tau_i` -- so a fresh run silently *appended* to (rather than renaming-old and overwriting) a pre-existing file, and the hardcoded `restart_time = 0.0` meant a genuine restart never deduplicated the redone tail. Fixed by detecting a real restart via `amr.restart` (not our own tau-relative clock, which can't distinguish "fresh start" from "restart" since `GRAmr::get_restart_time()` is 0 in both cases) and passing the actual restart time (`tau_i + get_restart_time()`) plus a correct `first_step`; `remove_duplicate_time_data()` is now called before every write (a verified no-op on a fresh run). The spectrum -- what a spectral index `q` actually gets measured from -- was previously only ever `amrex::Print()`'d, never saved; it is now persisted to a new `axion_spectrum.dat`, one row per mode index per snapshot with `tau` repeated as the literal first column (not blank-line-separated blocks, so the same flat restart-dedup logic applies exactly and safely). `AxionStringsParams::check_params()` now aborts at start-up if `axion_strings.compute_spectrum` is on without `axion_strings.masking.scheme = B`, since scheme A has no top-hat at all -- the persisted spectrum must be the genuinely screened field (the `plane_wave_test` IC's own one-shot spectrum check in `initData()` is unaffected: it hardcodes `MaskingScheme::None` locally regardless of this setting, since it is testing the FFT/Parseval identity, not screening). `tension_core_only`/`tension_core_plus_tail` are now also saved as columns (previously print-only) per the user: "we may as well save the processed tension, even if it can be reconstructed". Verified end-to-end against `straight_string_test`: (a) a coarse cadence (`delta=0.5`) demonstrably skips most coarse steps (3 snapshots out of 10, vs. 1 per step before); (b) a fresh rerun over pre-existing output files renames them to `.old.<random>` instead of mixing data; (c) restarting from an *earlier* checkpoint than the last snapshot and re-running forward past an already-recorded snapshot's tau produces no duplicate row in either file (checked both `network_scalars.dat` and the full 17-row `axion_spectrum.dat` block). One accepted wrinkle: a genuine restart always takes one extra, off-cadence snapshot exactly at the restart point (since `s_next_output_log_mr_over_h` resets to `first_log_mr_over_h` on every fresh process and immediately fires), which is harmless (no duplication, just one bonus sample) but not persisted across restarts -- left as is rather than adding checkpoint-side state for what is now an off-by-default, edge-case path. |
| Output infrastructure: unscreened spectrum + 2D visualisation | **Done** | Two follow-up asks (2026-09-18). (1) **Unscreened spectrum alongside the screened one**, "to allow for comparison and judging the impact of screening": `axion_spectrum.dat`'s columns are now `_screened`/`_unscreened` pairs (e.g. `shell_average_screened`, `shell_average_unscreened`), computed via two calls to `fill_masked_a_dot_buffer`/`compute_spectrum` per snapshot (screened using `s_energy_masking`, guaranteed scheme B; unscreened hardcoding `MaskingScheme::None`, mirroring the `plane_wave_test` IC's own check) -- doubles the FFT cost of an already opt-in diagnostic. Sanity-checked on `straight_string_test`: unscreened `<a_dot^2>` is consistently higher than screened, as expected once the cores are no longer masked out. (2) **Optional 2D visualisation snapshots**, off by default (`axion_strings.save_projection`), for spot-checking a run by eye rather than routine diagnostics: `ProjectionKernel.hpp`'s `compute_energy_projection` does a host-side (CPU-first, same call as `SpectrumKernel.hpp`'s binning: once per snapshot, not every substep) `MFIter` line-of-sight projection onto the xy-plane (z fixed as the line of sight), taking the **max** of the **unscreened** `rho_tot` per column -- max rather than sum/average so a thin core isn't diluted by a long quiet sight line; unscreened because screening would remove exactly what this is meant to show -- plus an overlaid `string_hit_count` (the xy-plaquette winding, `PlaquetteWinding.hpp`, summed over the line of sight) as an independent cross-check that the energy and string-finding pipelines agree on where the strings are. Saved to a new `axion_projection.dat`, one row per `(i,j)` pixel per snapshot (`tau` repeated as the literal first column, same restart-dedup trick as `axion_spectrum.dat`) -- a flat `(tau, i, j, value...)` table, straightforward to pivot into a 2D array for a coloured plot in any plotting tool. Verified end to end on `straight_string_test` (32x32 grid, 1024 rows/snapshot as expected): `string_hit_count` correctly picks out both vortex cores at their expected grid locations (count = `N_z` at each, one winding per z-layer, zero elsewhere). **Found a real bug this way, not just a visualisation quirk**: the projected energy's *global* max sat away from either core, in a smooth ridge on the far side of the box from both strings, exceeding the core's own value -- correctly flagged by the user as implausible ("the gradient energy density far from the core is much smaller... there should be no way for this to dominate"), and investigated rather than dismissed. **Root cause, confirmed both analytically (independent Python re-evaluation of the IC formula, no AMReX involved) and via instrumented debug prints of the actual running field**: `straight_string_test`'s phase field `theta = atan2(y-y1,x-x1) - atan2(y-y2,x-x2)` (`AxionStringsLevel.cpp`'s `initData`) is evaluated from literal, non-periodic `(x,y)` coordinates and is simply **not periodic** on the torus -- nothing in the formula knows the domain wraps at `y=Ly`. Both vortices sit at the same `y=y1=y2`, so the seam `y=0<->y=Ly` (diametrically opposite them) is exactly where this shows up worst: `psi2` jumps from `-0.81` to `+0.78` between two *physically adjacent* periodic grid points (a real discontinuity, confirmed present already in the pure, unevolved IC at `tau_i` via a from-scratch Python re-implementation of the exact formula -- not introduced by evolution or by any bug in the new energy/projection code, both of which were also directly ruled out by inspecting the actual ghost-cell values, which wrap correctly). That fake discontinuity then feeds a large, spurious gradient-energy contribution into the evolved solution from step 1 onward. **This is a distinct defect from the previously-documented (task 1.8) "product profile doesn't cleanly cancel between cores" tail issue** -- my first attempt at explaining this finding conflated the two; this one is specifically about the *phase* construction's non-periodicity, not the amplitude profile, and is large enough to be a real correctness problem for T1 (tension, and now core-location-by-energy), not a minor rough edge. Not yet fixed -- flagged to the user for a decision on the right periodic-aware construction (e.g. summing periodic images of each vortex's phase, or a different unwrapping) before relying on `straight_string_test` results, including the existing tension numbers, further. |
| Input parameterization: `log(m_r/H_i)` replaces `tau_i` | **Done** | Prompted by comparing our inputs against the user's old fixed-grid code (`main_global.cpp`): its `Hmri` argument (`H/m_r` at the start) is more physical than a raw conformal time, and the user asked to switch to it (2026-09-18). `Background::tau_from_log_mr_over_h` (new) is the exact inverse of the existing `H_over_mr_closed_form`, using `c_sched.c0` (the value in effect at/before the start -- same convention `apply_box_plan`/`pre_evolution_L_tilde` already use, and consistent with `H_over_mr_closed_form`'s own documented "only exact in the no-switch case" caveat, since any switch is expected to happen after `tau_i`, not before it). Singular at `c0 = a_inv` (Moore mode from the very start rather than reached via a switch) -- `AxionStringsParams::read_tau_i` (now takes a `const Background&`) checks for this explicitly and aborts with a clear message rather than dividing by zero silently. `axion_strings.tau_i` is gone; `axion_strings.log_mr_over_h_i` is the new required input. Verified: a new round-trip unit test (`test_background.cpp`, "`tau_from_log_mr_over_h` is the exact inverse of `H_over_mr_closed_form`") across several `(a_inv, c0, log_mr_over_h)` combinations, plus `tau_from_log_mr_over_h(0) == tau0`; `params_test.txt` updated to `log_mr_over_h_i = 0.0` (algebraically equivalent to the old `tau_i = 1.0` for its `a_inv=2, c0=0` configuration) and confirmed to reproduce byte-identical `tau`/`m_r/H` output at runtime; both new failure paths (`c0 = a_inv`, and the parameter simply missing) checked to abort with clear, specific messages rather than a generic ParmParse error. While reading the old code's full argument list to build this mapping, also matched up the rest of our `axion_strings.*` parameters against its `argv[]` positions for the user's reference (`N`/`N1`/`N2` for grid/core/Hubble resolution, `pre_evolution.gamma`/`seed` for pre-evolution, `dt_multiplier` for the old `R_delta`, `save_projection` for `domovie`) -- GW evolution, loop extraction, backreaction, the second "systematics" grid, and Moore's *extra* numerical trick (beyond just running with `c=1+b_inv`) remain unimplemented, matching what was already known/deferred. |
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
