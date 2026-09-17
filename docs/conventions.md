# Axion String AMR: Conventions

Snapshot of the living conventions document, rev 25, 2026-09-17.
Maintained outside the repo; propose changes rather than editing here.

---

## 1. Scope and status

Conventions for simulations of **global (axion) strings** in an expanding universe, for a new
AMR code built on [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn) (AMReX-based;
`Examples/KleinGordon` is the starting template).

Derived from the existing *Conventions for Abelian Higgs string simulations*, reduced to the
global case. Everything gauged is deliberately dropped: link variables, plaquette field
strengths, the Gauss constraint, the `Θ = R^c θ` rescaling, the compact formulation, and the
type I / type II distinction. Gravitational waves are **out of scope**, which removes the
`U_ij` evolution and the TT projection — the most AMR-hostile operation in the source document.

What remains: a single complex scalar with a Mexican-hat potential on a fixed FLRW background,
plus the diagnostics needed to measure the network and the axion emission spectrum.

**Physics target.** A high-statistics AMR study of the scaling regime — ensembles over initial
conditions, controlled refinement-boundary systematics, and defensible statistical
uncertainties on ξ and on the spectral index q and its dependence on log(m_r/H).

**Status.** Draft. Items marked *open* are unsettled; §14 records what has been settled and why.

---

## 2. Background cosmology

FLRW with scale factor `R(t) ∝ t^(1/a_inv)`, where `a_inv` is left as a code parameter. The
Hubble rate is `H = Ṙ/R = 1/(a_inv t)`.

Conformal time `τ = ∫ dt'/R(t')`, in which the metric is `g_μν = R² η_μν` with
`η_μν = (−1,1,1,1)`. Writing `R(τ) ∝ τ^(1/b_inv)` with

```
b_inv = a_inv − 1
```

gives conformal Hubble `H = 1/(b_inv τ)` and physical `H ∝ τ^(−1−1/b_inv)`.

| Era | a_inv | b_inv |
|---|---|---|
| Radiation domination | 2 | 1 |
| Matter domination | 3/2 | 1/2 |

**Decision:** `a_inv` stays general in the code. The `(1 − b_inv)/(b_inv² τ²)` term in the
equation of motion vanishes in radiation domination but is retained, since it costs nothing and
keeps matter domination available. All production runs assume radiation domination unless
stated otherwise.

The background is analytic and fixed: no metric is evolved. GRTeclyn's CCZ4 machinery is
unused, and the `KleinGordon` example — which assumes a flat non-dynamical background — is the
correct base.

---

## 3. Fields, potential and equation of motion

A complex scalar `φ` with

```
V(φ) = (λ/4) (|φ|² − v²)²,      m_r² = λ v²
```

Decomposing `φ = (1/√2)(√2 v + r) exp(i a /(√2 v))` gives the radial mode `r` of mass `m_r` and
the axion `a` of period `2π f_a`.

### Relation to the published convention

The *Axions from Strings* papers write `V = (m_r²/2 f_a²)(|φ|² − f_a²/2)²`. The two agree under

```
v = f_a/√2 ,      λ = 2 m_r²/f_a²
```

which is consistent with `m_r² = λ v²`. **This mapping must be applied whenever comparing to
published figures** — it is the most likely place for a silent factor-of-2 to enter.

### Equation of motion

In cosmic time,

```
φ̈ + 3H φ̇ − ∇²φ/R² + (λ/2) φ (|φ|² − v²) = 0
```

The field actually evolved is the rescaled `ψ = R(τ) φ / v`, which cancels the Hubble friction
entirely. In conformal time, with `∇` now the comoving gradient:

```
ψ'' − ∇²ψ − [(1 − b_inv)/(b_inv² τ²)] ψ + (λ/2) v² ψ (|ψ|² − R²) = 0
```

Primes are derivatives with respect to conformal time. The third term vanishes in radiation
domination (`b_inv = 1`) but is kept (§2).

**Note on the port.** Because the friction term is absent after rescaling, the RHS differs from
GRTeclyn's `KleinGordonRHS` only by the complex doubling of the state and the time-dependent
factors `λ(τ)` and `R²(τ)` in the potential term. Both are analytic functions of time available
from `a_time`. The state vector is

```
ψ₁, ψ₂, Π₁, Π₂     (NUM_VARS = 4)
```

with `Π_i = ψ_i'`.

---

## 4. The c(τ) scheme

All evolution modes are one parameter:

```
λ(τ) = λ₀ (R/R₀)^(−2c)
```

so that `m_r = √λ v ∝ R^(−c)`.

| c | Mode | Comoving core width `∝ R^(c−1)` | m_r/H |
|---|---|---|---|
| 0 | Physical | shrinks as 1/R | grows |
| 1 | Fat string | constant | grows, more slowly |
| 1 + b_inv (= 2 in RD) | Moore | grows as R | **constant** |

**In the global case `c` enters only through `λ(τ)`.** The `2cH θ'` term and the `Θ = R^c θ`
rescaling in the source document are gauge-sector artifacts and have no analogue here. This is
a genuine simplification, not merely a deletion: it is what makes a mid-run change of `c` clean.

### Switching c mid-run

`λ` appears algebraically in the equation of motion, so only `dλ/dτ` jumps at a switch — and
the EOM never sees it. Continuity of `λ` is imposed by re-anchoring at the switch time `τ_s`:

```
λ(τ) = λ_s (R/R_s)^(−2 c_new)      for τ > τ_s
```

The intended use is the **fat → Moore protocol**: run `c = 1` until the target log(m_r/H) is
reached, then switch to `c = 1 + b_inv` so the system is already near the correct attractor
when the constant-log phase begins.

### Implementation note

`c` is a small time-dependent function — a piecewise constant with at most one switch — not a
compile-time constant. The derived quantity `m_r/H` follows from

```
H/m_r = (τ/τ₀)^((c − a_inv)/(a_inv − 1))
```

**Keep this light.** The physical-string AMR programme is the point of the code; fat and Moore
are options exercised through `c(τ)` and should not shape the surrounding structure or the
naming. One function, one switch time, one parameter block.

---

## 5. Code units, parameters and box planning

Dimensional quantities are normalised by `v`, i.e. **`v = 1` in code units**. Quantities at the
time when `m_r = H` carry subscript 0. `δx` is dimensionful, so `R` is dimensionless.

```
t₀ = 1/(a_inv √λ₀ v)
R₀ = 1/((a_inv − 1) √λ₀)        so that   τ₀ = 1/v
```

`λ₀ = 1` in the code (its value does not affect the physics).

### Derived relations

```
t/t₀   = (τ/τ₀)^(a_inv/(a_inv−1))
R      = R₀ (t/t₀)^(1/a_inv) = R₀ (τ/τ₀)^(1/(a_inv−1))
m_r    = √λ₀ v (τ/τ₀)^(−c/(a_inv−1))
H/v    = √λ₀ (τ/τ₀)^(−a_inv/(a_inv−1))
H/m_r  = (τ/τ₀)^((c − a_inv)/(a_inv−1))
```

### Inputs

The two physical resolution parameters are specified **at the final time**:

- `N₁ = H(τ_f) L(τ_f)` — number of Hubble patches in the box at the end
- `N₂ = (Δ(τ_f) m_r(τ_f))⁻¹` — lattice points per string core at the end

with `Δ = L/N` the physical lattice spacing and `δx = L̃/N` the comoving one. Then

```
τ_f  = (N/(N₁N₂))^((a_inv−1)/(a_inv−c)) τ₀
L̃ v = (N₁/(R₀√λ₀)) (N/(N₁N₂))^((a_inv−1)/(a_inv−c))
```

For AMR, `N` and `N₂` refer to the **effective** finest resolution, not the base grid — see §11.

### Moore-phase box planning

During `c = 1 + b_inv` the box is fixed in comoving terms while `H⁻¹` grows faster, so
`HL ∝ 1/τ` while `N₂` only improves. The run ends when `HL` falls to `N₁`, giving (radiation
domination)

```
log(H₀/H)_max = 2 log[ N / (N₂ N₁ γ) ],     γ = m_r/H
```

equivalently `N = N₂ N₁ γ e^(D/2)` for dynamic range `D`.

Two consequences worth stating plainly:

1. The binding parameter is `N₂` **at the switch**, and the preceding fat phase holds `N₂`
   fixed — so it is set at the very start. **The fat pre-evolution and the Moore phase must be
   budgeted together, not sequentially.**
2. `N ∝ γ` at fixed dynamic range, so memory grows as `e^(3 log)`. Moore is a **moderate-log
   cross-check, realistically capped around log 5–6**. (Consistent with the source document's
   table: 1200 ok for 4.4, 1500 for 5.1, ≳2000 for 5.8.)

Also from the source document: what sets the initial ξ in the main run is `ξ_in × (m_r/H)²`.

### Time step

```
δτ ≲ δx / r_Δ ,     r_Δ ≳ 3
```

to be re-derived for the RK integrator (§6) rather than inherited from leapfrog.

---

## 6. Discretisation and time integration

**Spatial derivatives:** fourth order, using GRTeclyn's `FourthOrderDerivatives`. The published
work notes that second order gives similar results for string length, but fourth order is
retained as the default.

**Time integration:** method of lines with Runge–Kutta, as provided by GRTeclyn, with
subcycling in time on refined levels. This replaces the second-order leapfrog of the existing
code.

*Consequence:* step-size systematics must be re-derived, not carried across. Some drift
relative to the old code is expected on integrator grounds alone; a deterministic test with an
analytic reference (§13) is what distinguishes expected drift from a bug.

**Kreiss–Oliger dissipation: default `σ = 0`.**

This deserves emphasis. GRTeclyn's `KleinGordonRHS` calls `m_deriv.add_dissipation(...)` as a
matter of course. KO dissipation damps precisely the high-k modes whose spectral index is the
headline observable. Inherited silently it would be a systematic sitting directly on top of the
result. If a non-zero `σ` proves necessary for stability, its effect on the spectrum must be
characterised as carefully as the lattice spacing, and reported.

**Open — the improved (order a⁴) Laplacian.** The source document's §7 improvement (following
[1406.1688](https://arxiv.org/abs/1406.1688)) is tuned to the lattice dispersion relation of a
*uniform* grid. At a refinement boundary the interpolation carries its own dispersion error, so
the improvement may not deliver what it does on a fixed grid. Combined with the existing note
that it "seems to make the evolution worse", the decision is to **test this late rather than
port it early**.

*Terminology:* this is the **Moore improvement** (a stencil); §4's constant-log evolution is
**Moore's trick** (a scheme). They are unrelated. Use the full names throughout to avoid
confusion.

---

## 7. Initial conditions

### Generation

Initial conditions are **generated in the code**, not loaded from a field configuration file.
Two parameters control them:

- `k_max/m_r` — Fourier modes are occupied for `|k| ≤ k_max` and zero above. Larger `k_max`
  gives shorter-wavelength initial fluctuations and more initial strings, so this is the handle
  on initial string density.
- the **mean-square variance** of the field, set to a specified value.

A pre-evolution stage then follows, with `R = R₀(t/t₀)` and `λ(t) = λ(t₀)(t/t₀)^(−2)` so that
`R m_r` and `m_r/H` are both time-independent, run until the network reaches the desired ξ. The
field is then rescaled back to `φ` and handed to the main run.

The pre-evolution is what puts the network near the attractor before the main clock starts;
without it the system spends part of the run relaxing and releasing energy that contaminates
the spectrum. It is the direct analogue of "method (b)" in the published appendix. Comoving box
size for the pre-evolution:

```
L̃_init = L̃_main (1/(a_inv − 1)) (τ_i/τ₀)^((1−c)/(a_inv−1))
```

Initial conditions are specified by the triple `(ξ₀, H_c, H₀)`; the source document's §19 tables
of tested values for physical, fat and Moore carry over and should be re-tested rather than
assumed, since the integrator has changed.

### Restart protocol

Pre-evolution and main run are separate executions communicating through a **checkpoint file**.
This is required anyway for the fat → Moore protocol (§4) and for the fat-string AMR
configuration (§11), and it gives cheap ensemble management for free.

### Seeds and reproducibility

Every run records the RNG seed and the full parameter set in its checkpoint and in its output
metadata. A run must be **bit-reproducible** from `(seed, parameters, code commit)` on the same
number of ranks. This is not optional polish: the whole programme is an ensemble study, and a
systematic that cannot be reproduced cannot be measured.

### Milestone-1 validation path

Because initial conditions are generated independently in each code, the comparison against the
existing fixed-grid code is **statistical**: ensembles at matched parameters, compared on
ξ(log), the energy components and the spectra. Field-by-field agreement is not available and is
not sought.

The normalisations that a statistical comparison cannot pin are pinned instead by the
deterministic tests of §13, which is where round-off-level agreement is demanded.

---

## 8. String identification, ξ and velocities

### Plaquette winding

In the global case the gauge-invariant winding of the source document reduces to the plain
axion phase difference. A plaquette is pierced if traversing its four vertices accumulates a net
`2π` in `γ(x) = arg(φ)`, with each consecutive difference reduced to `(−π, π]`.

The published algorithm (1806.04677 App. A.2) is the reference implementation and is already
the community standard — Buschmann et al. adopt it as their primary AMR tagging criterion.

### String length per Hubble

`ξ = lim ℓ t²/V`, obtained from the number of pierced plaquettes `N_p`:

```
ξ = (2/3) N_p (δx / L̃³) v² ((a_inv − 1)³/a_inv²) (τ/τ₀)²
```

The `2/3` corrects for random orientation of string relative to the lattice. For Moore's trick
there is an extra factor of `H₀^(−2)` (source document eq. 84).

The same expression inverted gives the number of plaquettes that must be pierced to reach a
target `ξ_i` when setting initial conditions.

> **AMR caveat.** `δx` appears explicitly, so on a hierarchy the count must use `δx_ℓ` per level
> with covered coarse cells masked, or strings are double-counted. The `2/3` factor is derived
> for a uniform lattice; whether it survives near a coarse–fine interface is to be **tested
> directly**, not assumed. See §11.

### Velocities

Adapting [1707.05566](https://arxiv.org/abs/1707.05566), in code units:

```
γ²v² = (1/(2R⁴m_r²c₁²)) |ψ̇ − (β/τ)ψ|² (1 − e₁|ψ|²/(c₁³R²))
       − (e₁/(4R⁶m_r²c₁⁵)) (Re(ψ*ψ̇) − (β/τ)|ψ|²)²
```

For global strings the profile coefficients are the `e = 0` row of the source document's table:

| | c₁ | d | e₁ |
|---|---|---|---|
| Global | 0.41222(1) | 0 | −0.025763(1) |

Evaluate at the corners of each pierced plaquette and average over the network.

### Loop distribution and curvature

Both carry over. Curvature uses the non-uniform stencil of
[2412.08699](https://arxiv.org/abs/2412.08699) on segments of `s = 4n+1` points, currently
`s = 9`. **Note the known sensitivity:** `⟨κ⟩` varies by roughly a factor of 2 between `s = 9`
and `s = 13`, so `s` must be reported alongside any curvature result and held fixed across an
ensemble. Curvature output is normalised to the Hubble length.

---

## 9. Fourier transform conventions

Forward transform as implemented:

```
X̃_p = Σ_{n=0}^{N−1} X_n exp(−2πi (p·n)/N)
```

with `p_i = 0, …, N−1`. **No normalisation factor** is applied on either direction, so a
forward followed by a backward transform multiplies every entry by `N³`.

Memory layout: positive frequencies counting up in the first half, then most-negative to
least-negative in the second half. Largest stored frequency `N/2 − 1`.

Continuum convention — all factors of `2π` on the momentum integral, `k = 2πp/L`:

```
f(x) = ∫ d³k/(2π)³ exp(ik·x) f̃(k)
Σ_{p=0}^{N−1} → (L/2π) ∫ dk
(L/N) X̃_p → f̃(k)
```

These are the conventions every normalisation in §10 depends on. Any FFT library substituted
for the current one must be checked against all three relations, not just the first.

---

## 10. Energies, masking and spectra

### Total energy

```
ρ_tot = R⁻² v² ⟨ |ψ̇ − ψ/(b_inv τ)|² + |∇ψ|² + (λv²/4R²)(|ψ|² − R²)² ⟩
```

Spatial averages are averages over lattice points.

**Every energy component is computed both screened and unscreened.** The string tension is
obtained from the *difference* of the two: the energy removed by screening is the string
contribution, and `ρ_tot = ρ_s + ρ_a + ρ_r` is how the split is defined. So an unscreened pass
is not a debugging convenience — it is required for one of the headline observables, and the
masking module must support `f ≡ 1` as a first-class mode. Note that the effective-point-count
correction below applies to the screened averages only.

### Axion energy and spectrum

With `ρ_a = ⟨ȧ²⟩` (note the factor of 2 convention) and `∫ d|k| ∂ρ_a/∂|k| = ρ_a`:

```
∂ρ_a/∂|k| = (|k|²/(2πL)³) ∫ dΩ_k |ã̇(k)|²
```

From the discrete transform, `getSpectrum` returns `(p, 4π ⟨|p|² |X_p|²⟩)` averaged over modes
with `|p|` within ±0.5 of `p`. To obtain `(k/H, v⁻³ ∂ρ_a/∂k)` rescale by

```
( 2π/(L H) ,  (1/2π) R L̃ / N⁶ )
```

from `∂ρ/∂k = (2L/N⁶) ⟨p² |X̃_p|²⟩`.

**Cross-check to implement:** `ρ_a,kin` computed as `⟨½ȧ²⟩` over lattice sites must agree with
`(1/2N⁶) Σ |ã̇(p)|²`. Cheap, and catches normalisation errors immediately.

Radial mode spectra follow identically from `ρ_r = ½ṙ² + ½(∇r)² + ½m_r²r²`.

### Masking — both schemes, compared

Strings must be screened out before transforming. **Two schemes are implemented and the
difference between them is measured, not assumed:**

| Scheme | Form | Provenance |
|---|---|---|
| A — smooth | `f = (1 + r/f_a)²`, i.e. `ȧ_scr = ψ₁ψ̇₂ − ψ̇₁ψ₂` | published papers; also Buschmann et al. fiducial |
| B — top-hat | `f = 1` where `|ψ|²` exceeds `x`, else 0 | source document §23 |

A screening distance `d_s` in units of `m_r⁻¹` parameterises scheme B's effective reach and
should be reported.

**The masking threshold is a runtime parameter, not a compile-time constant — deliberately, so
it can be scanned.** The existing fixed-grid code uses a top-hat at `|ψ|/R < 0.8`; §23 of the
source conventions records testing `x = 0.9, 0.95` and notes `mod < 0.75` as good. These are not
reconciled, and the value may also differ between that code and the production version. **To be
tested rather than settled now:** run the threshold over its plausible range and report the
induced spread in q alongside the statistical error.

**For unbiased averaging, the effective number of lattice points averaged over must be
multiplied by the same masking factor.** This is an easy omission and produces a smooth,
plausible, wrong answer.

The A/B difference is expected at the percent level — which is exactly the scale at which q's
log-dependence lives. That is why both are needed.

### Extracting q

The instantaneous emission spectrum `F(k)` is obtained by finite-differencing the spectra in
time with the appropriate redshift factors, and `q` from a fit over a momentum range bounded
below by finite-volume effects and above by core-scale contamination. **Both bounds must be
recorded with every quoted q**, since the two groups' disagreement is partly a disagreement
about fit ranges.

---

## 11. AMR conventions and open questions

### What each level buys

In comoving coordinates with a fixed box, physical spacing grows as `R ∝ (m_r/H)^(1/2)`, so
points per core fall as `(m_r/H)^(−1/2)` and **each refinement level buys a factor of 4 in
`m_r/H`**:

```
Δlog per level = ln 4 ≈ 1.39
```

Cross-check: Buschmann et al. add levels at log ≈ 2.6, 3.9, 5.3, 6.7 — spacings of 1.3–1.4, as
predicted. Their final gap of 2.0 to log 8.7 is resolution being allowed to degrade before
paying for the last level.

### Refinement demand by mode

The refined volume fraction is (string length per comoving volume) × (comoving core width)²:

```
f_refined ≈ 4ξ (H/m_r)²
```

with a buffer and block-structured granularity inflating this by perhaps two orders of
magnitude.

| Mode | Comoving core width | Demand over time | AMR verdict |
|---|---|---|---|
| Physical (c=0) | shrinks | grows, and cheaply — demand rises exactly as supply gets cheap | **the target** |
| Fat (c=1) | constant | binding at the *start*; f ~ 1 below log ≈ 2 | viable **only as a restart** at log ≈ 4–5 |
| Moore (c=1+b_inv) | grows | maximal at the start, decreasing | **run single-level**; AMR pays peak cost for no asymptotic saving |

The anticorrelation in the physical case — demand growing at exactly the rate the cost of
supplying it falls — is the whole reason AMR suits this problem. Neither other mode has it.

*Practical consequence:* the Moore programme needs nothing beyond milestone 1
(`amr.max_level = 0`).

### Tagging

Primary criterion: pierced plaquettes (§8), tagging the low-index corner. Secondary criteria to
be assessed rather than copied — Buschmann et al. add a gradient criterion `Δx_ℓ²∇²ψ > 0.04`
and a coarse-level radial-mode criterion, both with admittedly phenomenological thresholds.

Buffer sizing: enough that the fastest string (`v = c`) stays a full core width from any
coarse–fine boundary between regrids. Their choice was 11 cells with regrid interval
`Δη_ℓ = 0.2/2^ℓ`; adopt as a starting point and **test**, since it was tuned for a different
string density.

### The level-timing systematic — *the central AMR risk*

Effective resolution changes **discontinuously during the run**, at particular values of log.
The disputed observable is not q but `dq/dlog`. Any resolution-dependent bias in q is therefore
imprinted as a spurious contribution to its trend — and level-addition times are by
construction perfectly correlated with the fit abscissa.

This failure mode **does not exist on a fixed grid**, where resolution degrades smoothly and
monotonically. It is created by the AMR.

**Test:** run identical initial conditions with levels added at different logs (and with levels
added earlier than needed, so resolution is never binding), and ask whether the inferred
`dq/dlog` moves. Cheap, decisive, and plausibly the source of the disagreement in the
literature.

### Diagnostics on the hierarchy

- **ξ on a hierarchy.** Level-dependent `δx_ℓ`, covered coarse cells masked. Does the `2/3`
  orientation factor survive near an interface?
- **Spectra on a hierarchy.** Settled — average fine data down and FFT the coarse level,
  following the existing single-FFT conventions of §9–10. Then the coarse Nyquist bounds the
  usable fit range — that bound must be verified against a uniform run and the fit range
  reported with every quoted q. Level-decomposed or composite spectra are rejected: a bespoke
  estimator would risk introducing AMR-dependent structure into the very observable in dispute,
  and the coarse grid carries ample modes to fit the power law.
- **Energy accounting across interfaces.** The cleanest global diagnostic: total box energy
  should redshift correctly, and any leak localised at coarse–fine boundaries shows up here
  before it reaches the spectrum.
- **Interface effects on evolution.** Probe with a straight string and a collapsing loop placed
  deliberately across a refinement boundary, where the answer is known.
- **Interpolation order.** GRTeclyn registers derived variables with
  `amrex::cell_quartic_interp`. Interpolation choice sets spurious reflection of high-k
  radiation at interfaces — a parameter to study, not accept.

---

## 12. Observables inventory

Everything the new code must reproduce, as extracted from the fixed-grid implementation.
Gravitational-wave outputs are omitted (§1).

### Network scalars

| Quantity | Definition / notes |
|---|---|
| `log(m_r/H)` | the abscissa for everything; `−((c − b_inv − 1)/b_inv) log(τ/τ₀)`, with a separate branch for constant-log runs |
| ξ | from pierced plaquettes, §8. Each cell contributes up to 3 (one per plane) |
| ξ_W | winding-weighted variant; most strings are singly wound, so the two agree closely — a divergence is a diagnostic in itself |
| γ²v² distribution | 500 bins, γ ∈ [1,10], evaluated at corners of pierced plaquettes and averaged over the network |
| κ distribution | 200 bins, normalised to the Hubble length, segments of `s = 4n+1` points. **Report `s`** (§8) |
| Loop distribution | per loop: start, end, length, closed flag |

### Energy components

All screened, all spatially averaged. **The effective point count must carry the same masking
factor** (§10).

- **Radial:** kinetic `½ṙ²`, gradient `½(∇r)²`, mass `½m_r²r²`, and their sum
- **Axion:** kinetic `½ȧ²`, gradient `½(∇a)²`
- **Interaction:** total, plus separate radial-side and axion-side parts
- **Total field:** kinetic, potential, spatial, and their sum
- **Unscreened counterparts of all of the above**, from which the string energy and hence the
  tension follow by difference

The split matters because the redshifting of each component is what converts a total-energy
time series into an emission rate.

### Spectra

For each of the axion and radial sectors, a component index selects which quantity is
transformed:

| Component | Transformed field | Extra factor |
|---|---|---|
| kinetic | masked `ȧ` (resp. `ṙ`) | — |
| gradient | masked `a` (resp. `r`) | `(2π/RL)² |p|²` |

Note the gradient spectrum is built by transforming the *field* and multiplying by `k²`, not by
transforming its gradient. The two differ by lattice dispersion at high k, so **the convention
must be stated with any published spectrum**.

**Binning.** `modeIndex = floor(p_mod + 0.5)` where `p_mod` is the mode modulus; accumulate
`4π p_mod² |X_p|²` into the shell; divide by the shell count at the end. Modes outside the
inscribed sphere (`modeIndex ≥ N/2`) are dropped.

**Normalisation** to `(k/H, v⁻³ ∂ρ/∂k)`: multiply by `(2π/(LH), (1/2π) R L̃ / N⁶)` — §10.

**Built-in cross-check.** The existing code reports the ratio of full-cube to inscribed-sphere
energy, and the 1D energy from the spectrum against the real-space average. Both should be
preserved; they are cheap and catch normalisation errors on the spot.

### Derived

- `F(k)`: instantaneous emission spectrum, by finite-differencing spectra in time with the
  appropriate redshift factors
- `q`: power-law index of `F(k)`, fit over a recorded momentum range. **The IR and UV bounds of
  the fit are reported with every quoted q** (§10)

### Open

Masking threshold, pending the production version (§10). The definitions above are taken from
the non-production code, but the observables themselves are unchanged between versions.

---

## 13. Acceptance tests

Deterministic tests with known answers. **These are not a supplement to statistical validation
against the old code — they are the only thing that catches a class of error statistics cannot
see.** A bug that leaves a diagnostic plausible but wrong gets *harder* to spot as the ensemble
grows, because the error bars shrink around the wrong answer.

All of these run at the single-level stage, and each has an AMR variant (below).

### T1 — Static straight string

Infinite string along z, `φ = (v/√2) g(m_r ρ) e^(iθ)`, with `g(ρ) = c₁ρ + O(ρ³)` at small ρ and
`1 − ρ⁻² + O(ρ⁻⁴)` at large ρ; `c₁ = 0.41222` for the global case (§8).

Pins: core profile and its equilibrium; the tension via the energy integral, including the
expected logarithmic divergence with box size; plaquette detection — exactly the expected set of
plaquettes should be flagged; the ξ normalisation, since a known string length must give a
known ξ.

### T2 — Mask unit test

**A static string does *not* test the mask.** With `φ` time-independent, `ψ = Rφ/v` gives
`ψ' = (R'/R)ψ`, so

```
ψ₁ψ'₂ − ψ'₁ψ₂ = 0   identically
```

and the masked `ȧ` is zero everywhere regardless of whether the mask is applied at the right
lattice site. Any test built on a static configuration is blind to exactly the failure mode we
care about.

Two tests that are not blind to it:

- **Direct unit test.** Assert that the buffer written for the FFT is zero at precisely the set
  of points where the mask predicate holds, and bit-identical to the unmasked field elsewhere.
  A pure array comparison, no physics required, and it catches an index slip immediately.
- **Boosted straight string.** A string moving with known velocity has `θ̇ ≠ 0` near the core,
  so correct masking must zero a known moving cylinder. Doubles as the test of the `γ²v²`
  estimator (§8) against the imposed boost.

### T3 — Plane wave / single travelling wave

A small-amplitude wave in the axion field at known `k` and amplitude, no strings.

Pins: convergence order of the spatial stencil and of the time integrator; lattice dispersion;
and the **spectral normalisation end-to-end** — a known amplitude at a known mode must integrate
to the known total energy through the full binning and rescaling chain of §12.

### T4 — Collapsing circular loop

Pins: dynamics against the Nambu–Goto expectation for the radius; energy loss into radiation;
the loop finder and its closed-loop flag; and, at collapse, the radial-mode channel.

### T5 — Float versus double

Identical initial conditions, single level, at several values of log. Compare ξ, the axion
spectrum, the radial energy and total energy conservation. **Run at more than one log** — the
cancellation in `|ψ|² − R²` degrades with log, so a pass at low log certifies nothing at high
log (§14).

### AMR variants

Each of T1–T4 repeated with the object placed deliberately **across a refinement boundary**, and
with the boundary moving through it. These are the tests that quantify interface systematics,
and they are only meaningful because the fixed-grid answer is known.

### CI

T1–T3 should run in continuous integration from before any network evolution exists. They are
cheap, and they guard the `ȧ` path that every headline result depends on.

---

## 14. Decision log

Dated record of settled decisions and their reasons. Append; do not rewrite.

**2026-09-17**

- **Base code: GRTeclyn, `Examples/KleinGordon` template.** It solves a scalar field on a flat
  non-dynamical background, so "turning GR off" means not using CCZ4 rather than disabling it.
  Inherits AMReX AMR, fourth-order derivatives, periodic BCs, GPU-ready kernels and
  derived-variable I/O.
- **Gravitational waves out of scope.** Removes the `U_ij` evolution and, more importantly, the
  TT projection's global FFT — the most AMR-hostile operation in the source conventions.
- **Gauge sector dropped.** Global strings only. No Gauss constraint, no link variables, no
  type I/II.
- **`a_inv` kept general**, with the `(1−b_inv)/(b_inv²τ²)` term retained even though it
  vanishes in radiation domination. Revisit only if specialising proves materially simpler.
- **Both masking schemes implemented** (§10), with their difference measured rather than
  assumed.
- **RK / method of lines accepted** over leapfrog, to avoid reimplementing GRTeclyn's level and
  subcycling machinery.
- **KO dissipation defaults to zero** (§6) — it damps the modes the headline observable is
  measured from.
- **Moore's trick included as an option via `c(τ)`, not as a structural feature.** It is
  supporting evidence, not the main result; the code should read as a physical-string AMR code.
- **Moore and fat-string phases run single-level.** Their refinement demand is front-loaded, so
  AMR costs peak effort for no asymptotic saving.
- **Milestone 1 = single-level port validated against the existing fixed-grid code**, using
  initial data read from file rather than regenerated. Statistical agreement is the working
  criterion, supplemented by one deterministic test with an analytic answer to pin
  normalisations.
- **CPU first, GPU later** — but all kernels kept GPU-clean from the first commit, since the
  retrofit is painful and the discipline is free.
- **q is extracted from the coarse grid only.** The fit range sits well below the core scale, so
  those modes are resolved on the base level regardless of refinement. Level-decomposed or
  composite spectral estimators are rejected as too dangerous — they would let the mesh
  structure enter the disputed observable.
- **Masking threshold left as a scannable runtime parameter** (§10), with its effect on q to be
  measured later. The code and the source conventions currently disagree on the value, and the
  production code may differ from both.
- **Scheme A and scheme B are nested.** The fixed-grid code's
  `adt = (ψ₁ψ̇₂ − ψ̇₁ψ₂)/|ψ|²`; the published `ȧ_scr` is the bare numerator. So A/B is a
  one-line switch — divide and apply the top-hat, or do neither — which makes the comparison
  essentially free.
- **Grid precision: float to be tested, double the fallback.** GRTeclyn uses `amrex::Real` on
  the evolution path, so precision is a build flag. The specific risk is the cancellation in
  `|ψ|² − R²`: since `R² ∝ m_r/H`, float headroom erodes as log grows, so **a test at low log
  does not certify high log**. The axion phase `atan2(ψ₂, ψ₁)` is unaffected by this (both
  components are O(R)), so the spectrum and q are the least exposed quantities and the radial
  energy is the early-warning diagnostic. Reductions, energy sums and spectral accumulation stay
  in double regardless. Float buys ensemble size (×2 at fixed N), not reach (≈0.46 in log).
- **Energies are computed screened *and* unscreened.** The string tension comes from the
  difference, so the unscreened pass is a required observable rather than a diagnostic extra
  (§10, §12).
- **No file-based initial data.** Initial conditions are generated in-code from two parameters —
  occupied Fourier modes up to `k_max/m_r`, and a specified mean-square variance — followed by
  the pre-evolution stage (§7). This supersedes the earlier plan to load a configuration written
  by the fixed-grid code, and makes the milestone-1 comparison statistical rather than
  field-by-field.

### Corrections

- **2026-09-17.** An earlier estimate that each refinement level buys `ln 2 ≈ 0.69` in log was
  wrong; the correct figure is `ln 4 ≈ 1.39` (§11). Six levels therefore buy ≈ 8.3 in log, not
  ≈ 4.2 — the reach is considerably better than first stated.
- **2026-09-17.** An earlier claim that a static straight string would expose a misaligned
  screening mask was wrong. For any time-independent `φ`, the masked `ȧ` vanishes identically,
  so a static configuration is blind to that failure mode. A boosted string or a direct
  array-level unit test is required instead (§13, T2).
