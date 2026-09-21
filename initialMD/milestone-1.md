# Milestone 1 — single-level port

**Goal:** a complex scalar on a fixed FRW background in GRTeclyn, running on one level, whose
diagnostics reproduce the existing fixed-grid code.

No refinement in this milestone. `amr.max_level = 0` throughout. Getting the physics and the
diagnostics right on a uniform grid is a prerequisite for every AMR question that follows, and
the single-level code is also the production configuration for the Moore runs (conventions
§11), so it is not throwaway scaffolding.

Read `../CLAUDE.md` and `conventions.md` first. Section references below are to
`conventions.md`.

---

## Tasks

### 1.1 — Skeleton from the KleinGordon example

Copy `Examples/KleinGordon` to a new example directory and strip it to a working build. Keep
the structure: `StateVariables.hpp`, an RHS header plus `.impl.hpp`, a `Level` class, a
`SimulationParameters`, a `Main_*.cpp`, a `GNUmakefile`.

*Done when:* it builds and runs the unmodified wave model end to end under the new name.

### 1.2 — Complex state

`NUM_VARS = 4`: `ψ₁, ψ₂, Π₁, Π₂` (§3). Names in `StateVariables::names`, parities set (unused
under periodic BCs but must be defined).

### 1.3 — RHS with c(τ)

Implement (§3):

```
Π_i' = ∇²ψ_i + [(1 − b_inv)/(b_inv² τ²)] ψ_i − (λ(τ)/2) v² ψ_i (|ψ|² − R²(τ))
ψ_i' = Π_i
```

with `λ(τ) = λ₀ (R/R₀)^(−2c)` and `c` supplied by a small time function supporting one
mid-run switch (§4). `R(τ)` and `λ(τ)` are analytic; compute them from `a_time` **outside** the
kernel and pass plain values in.

Constraints from `CLAUDE.md`: no `ParmParse` in the RHS; `σ = 0` for dissipation; kernel body
GPU-clean.

*Done when:* the three modes (c = 0, 1, 1+b_inv) all run, and `m_r/H` measured from the code
matches the analytic expression in §5 to round-off.

### 1.4 — Parameters and box planning

A parameter block implementing §5: inputs `N`, `N₁`, `N₂`, `a_inv`, `c` (and switch time, if
used), `λ₀`. Derive `τ₀`, `τ_f`, `L̃`, `R₀` and the time step. Echo every derived quantity at
startup and write them to the output metadata.

Include the Moore-phase dynamic-range formula as a startup check: if the requested
configuration cannot achieve the requested range, say so before burning the allocation.

### 1.5 — Initial conditions and restarts

Generate initial conditions in-code (§7), controlled by two parameters:

- `k_max/m_r` — occupy Fourier modes for `|k| ≤ k_max`, zero above. This sets the initial
  string density.
- the mean-square variance of the field.

Then the pre-evolution stage: `R = R₀(t/t₀)` and `λ(t) = λ(t₀)(t/t₀)^(−2)`, run until the
target ξ is reached, then rescale and hand off to the main run. Pre-evolution and main run
communicate through a checkpoint, so restart must work — including a change of `c` at restart,
which the fat→Moore protocol needs.

Record seed, parameters and code commit in checkpoints and output metadata. Bit-reproducibility
from `(seed, parameters, commit)` at fixed rank count is a requirement, not an aspiration.

### 1.6 — String finder and ξ

Plaquette winding per §8, counting each of the three planes per cell. Both plain and
winding-weighted counts. ξ from the §8 formula, with the Moore branch.

*Done when:* T1 passes — a straight string of known length gives the expected plaquette count
and the expected ξ.

### 1.7 — Masking, behind one interface

A single masking module applied at exactly one point: the field written into the FFT buffer
(`CLAUDE.md` constraint 5). Both schemes (§10):

- A: `ȧ_scr = ψ₁Π₂ − Π₁ψ₂`
- B: the same divided by `|ψ|²`, with a top-hat at a **runtime** threshold

The effective point count for averaging must carry the same factor (§10).

*Done when:* T2 passes — the array-level unit test, plus a boosted string.

### 1.8 — Energies, screened and unscreened

The components listed in §12, computed **both ways**. The screened/unscreened difference is
what gives the string energy and hence the tension, so the unscreened pass is a required
observable, not a debugging aid. The masking module from 1.7 must therefore support `f ≡ 1` as
a first-class mode, and the effective-point-count correction applies to the screened averages
only.

Implement the real-space vs spectral cross-check for the axion kinetic energy (§10) and run it
every snapshot; it is cheap and catches normalisation errors immediately.

*Done when:* a straight string of known length reproduces the expected tension from the
screened/unscreened difference, including its logarithmic dependence on box size (T1).

### 1.9 — Spectra

FFT with the conventions of §9. Binning and normalisation per §12. Preserve the two built-in
cross-checks: full-cube vs inscribed-sphere energy, and 1D spectral energy vs the real-space
average.

*Done when:* T3 passes — a known plane wave integrates to its known energy through the whole
chain.

### 1.10 — Remaining diagnostics

Velocities (§8, global coefficients `c₁ = 0.41222`, `e₁ = −0.025763`), curvature with `s`
recorded, loop distribution. These are needed for the comparison but are not on the critical
path for q.

---

## Acceptance

**Deterministic (must pass, in CI):**

| Test | Guards |
|---|---|
| T1 static straight string | profile, tension from the screened/unscreened difference, plaquette detection, ξ normalisation |
| T2 mask unit test + boosted string | the entire `ȧ` path; the velocity estimator |
| T3 plane wave | stencil and integrator convergence order; spectral normalisation end-to-end |
| T4 collapsing loop | dynamics against Nambu–Goto; loop finder |

**Statistical (against the fixed-grid code, matched parameters):**

Initial conditions are generated independently in each code, so this comparison is statistical
over ensembles — not field-by-field. Round-off-level agreement is demanded of the deterministic
tests above instead.

- ξ(log), energy components and spectra agree within the ensemble scatter
- residual drift attributable to the leapfrog → RK change, demonstrated rather than assumed
- same convergence order under refinement of `m_r Δ`
- the screened/unscreened difference reproduces the published tension behaviour

**T5, float vs double:** run at several values of log, not one (§13, §14).

---

## Explicitly out of scope here

Refinement of any kind; tagging; interface systematics; the level-timing study; gravitational
waves; the improved (order a⁴) Laplacian.

---

## Questions to bring back rather than decide

- Any discrepancy between these definitions and the **production** fixed-grid code, especially
  the masking threshold (§10, open).
- Any case where matching the old code would require changing a convention in §1–§12.
- Any systematic that turns out to correlate with `log(m_r/H)`.
