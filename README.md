# Axion Strings AMR

An adaptive-mesh-refinement code for simulating global (axion) strings in an
expanding universe, built on [GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn)
(AMReX-based). See [`CLAUDE.md`](CLAUDE.md) and [`docs/conventions.md`](docs/conventions.md)
for the physics/engineering conventions, and [`docs/milestone-1.md`](docs/milestone-1.md)
for the current work plan. [`docs/STATUS.md`](docs/STATUS.md) tracks task-by-task progress.

## Layout

- `GRTeclyn/` — git submodule, upstream engine. Not modified.
- `amrex/` — git submodule, AMReX (GRTeclyn's dependency), sibling to `GRTeclyn/`
  as its build system expects.
- `AxionStrings/` — our simulation code, an out-of-tree "example" that points
  at `GRTeclyn/` and `amrex/` rather than living inside the submodule.
- `docs/` — local copies of the canonical conventions/milestone docs (see
  `CLAUDE.md` for how these relate to the externally maintained originals).

## Building

```bash
git submodule update --init --recursive
cd AxionStrings
make COMP=llvm -j4
```

### macOS toolchain note

Homebrew's LLVM (`/opt/homebrew/opt/llvm/bin`) is commonly first on `PATH` on
this machine, and its `clang`/`clang++` conflict with the macOS SDK headers
when driven through Open MPI's `mpicxx` wrapper (`FP_INFINITE`/`FP_NORMAL`
etc. undeclared in `<math.h>`). The fix is to build with a `PATH` that puts
`/usr/bin` (Apple clang) ahead of the homebrew LLVM directory, so `mpicxx`
resolves to Apple clang instead:

```bash
env PATH="/usr/bin:/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/bin:/usr/sbin:/sbin" \
    make COMP=llvm -j4
```

`COMP=llvm` is correct here (not `COMP=gnu`) because the underlying compiler
is clang-family either way.

## Running

```bash
env PATH="/usr/bin:/opt/homebrew/bin:/opt/homebrew/sbin:/usr/local/bin:/bin:/usr/sbin:/sbin" \
    mpirun -n 2 AxionStrings/AxionStrings3d.llvm.MPI.ex AxionStrings/params_test.txt
```

Pick the evolution mode with `axion_strings.c0` (0 = physical, 1 = fat
string, `1 + b_inv` = Moore; `b_inv = 1` in radiation domination, so `c0=2`):

```bash
... AxionStrings/params_test.txt axion_strings.c0=2.0
```

## Unit tests

`Background.hpp`, `BoxPlan.hpp`, `XiFormula.hpp`, `PlaquetteWinding.hpp`,
`PreEvolutionBackground.hpp`, `Masking.hpp` and `Energy.hpp` are
deliberately independent of AMReX, so their tests build and run in seconds
with no MPI/AMReX dependency:

```bash
cd AxionStrings/tests
clang++ -std=c++20 -O1 -I ../../GRTeclyn/external/doctest \
    test_background.cpp test_box_plan.cpp test_xi_formula.cpp \
    test_plaquette_winding.cpp test_pre_evolution_background.cpp \
    test_masking.cpp test_energy.cpp \
    -o /tmp/test_all
/tmp/test_all
```

## Fourier-mode initial conditions (task 1.5)

`FourierIC.hpp` uses `amrex::FFT::R2C`, which needs FFTW and AMReX's FFT
module. The `GNUmakefile` already sets `USE_FFT := TRUE` and
`FFTW_HOME ?= /opt/homebrew/opt/fftw` (adjust `FFTW_HOME` if FFTW lives
somewhere else on your machine; install it with `brew install fftw` if
needed).
