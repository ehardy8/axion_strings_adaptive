# CLAUDE.md

Standing brief for Claude Code sessions on this repository. Read this first, then
`docs/conventions.md` before touching physics or diagnostics.

## What this project is

An adaptive-mesh-refinement code for simulating **global (axion) strings** in an expanding
universe, built on [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn) (AMReX-based).

The scientific goal is **not** simply to reach higher `log(m_r/H)`. It is a high-statistics
study: ensembles over initial conditions, controlled refinement-boundary systematics, and
defensible statistical uncertainties on the string density `ξ` and on the axion emission
spectral index `q` and its dependence on `log(m_r/H)`.

This matters for how you write code. **AMR here is a variance-reduction tool.** A uniform-grid
run at fixed `HL` and `m_r Δ` costs roughly `e^(4 log)`, so ensembles become unaffordable
exactly where the physics is interesting. AMR buys back the budget to run *many* realisations,
not one heroic one. Cheap restarts, strict reproducibility and compact per-run diagnostics are
therefore first-class requirements, not polish.

## Starting point

`Examples/KleinGordon` in GRTeclyn is the template, not `Examples/ScalarField` and not the
CCZ4 machinery. KleinGordon solves a scalar field on a flat, non-dynamical background — so
"turning GR off" means *not using* CCZ4 rather than disabling it.

What it gives us: a `model_t`-templated potential, `FourthOrderDerivatives`, periodic
boundaries, GPU-ready `amrex::ParallelFor` kernels, and the derived-variable I/O system.

What we must add: complex field (4 state variables), the `c(τ)` time-dependence, a
string-based tagger, file-based initial data, and all diagnostics.

## Hard constraints

These are decisions already made, with reasons. Do not quietly reverse them.

1. **Kreiss–Oliger dissipation defaults to `σ = 0`.** The KleinGordon RHS calls
   `m_deriv.add_dissipation(...)` as a matter of course. KO dissipation damps precisely the
   high-k modes that `q` is measured from. If non-zero `σ` proves necessary for stability, it
   must be characterised as carefully as the lattice spacing and reported.

2. **Kernels stay GPU-clean from the first commit.** No host-only calls inside a
   `ParallelFor` body; resolve all parameters into plain values before the launch. Production
   is CPU-first, GPU later — but the retrofit is painful and the discipline is free.
   Corollary: do not call `amrex::ParmParse` inside `specific_eval_rhs` (the KleinGordon
   example does this every step; hoist it out).

3. **Reductions, energy sums and spectral accumulation stay in double**, whatever precision
   the grids use.

4. **The masking threshold is a runtime parameter**, never a compile-time constant. It will be
   scanned.

5. **Masking is applied at one place** — the field written into the FFT buffer — behind a
   small interface. The reference fixed-grid code inlines `loc_screen_fact(...)` at four
   separate call sites, and an index slip in one of them went unnoticed precisely because
   there was no single place where it was obviously wrong.

6. **Fit ranges are recorded with every quoted `q`**, and the masking scheme with every
   quoted spectrum. Part of the disagreement in the literature is a disagreement about these.

7. **Reproducibility.** A run must be bit-reproducible from `(seed, parameters, code commit)`
   on the same rank count. Seed and full parameter set go into the checkpoint and the output
   metadata.

## Design shape

- **State:** `ψ₁, ψ₂, Π₁, Π₂` (`NUM_VARS = 4`), with `Π_i = ψ_i'`.
- **One RHS**, templated on the model. Physical, fat and Moore evolution differ *only* through
  `c(τ)` entering `λ(τ)`. Do not build separate code paths for them.
- **Keep Moore's trick light.** It is supporting evidence, not the main result. It must not
  shape the surrounding structure or the naming — this should read as a physical-string AMR
  code that happens to support other modes.
- **Two run modes from one binary:** single-level (`amr.max_level = 0`) for Moore runs and for
  milestone-1 validation; AMR for the physical-string ensemble.
- **Restarts are a core feature**, not an afterthought. Needed for the fat→Moore protocol, for
  fat-string AMR (which is viable only as a restart), and for ensemble management.

## Testing discipline

Deterministic tests with known answers are **not a supplement** to statistical validation —
they are the only thing that catches errors statistics cannot see. A bug that leaves a
diagnostic plausible but wrong gets *harder* to spot as the ensemble grows, because the error
bars shrink around the wrong answer.

`docs/conventions.md` §13 specifies the tests. T1–T3 should run in CI from before any network
evolution exists.

One non-obvious point, worth internalising: **a static string does not test the mask.** For
time-independent `φ`, the masked `ȧ` vanishes identically. Use the direct array-level unit
test, or a boosted string.

## Where things live

| | |
|---|---|
| `docs/conventions.md` | physics, code units, observable definitions, decision log. Canonical. |
| `docs/milestone-1.md` | the current work plan and its acceptance criteria |
| `CLAUDE.md` | this file — standing constraints |

`docs/conventions.md` is a **snapshot** of a living document maintained outside the repo. Treat
it as read-mostly: propose changes rather than editing it directly, so the two do not diverge.
Implementation notes, build instructions, test status and known issues belong in the repo and
are yours to maintain.

## When to ask rather than decide

Ask about: anything that changes a physics convention, a normalisation, or an observable
definition; anything that would introduce a systematic correlated with `log(m_r/H)`; and any
choice of numerical parameter that is not already fixed in `docs/conventions.md`.

Decide freely: code organisation, build configuration, test scaffolding, I/O plumbing,
performance work that does not change results.
