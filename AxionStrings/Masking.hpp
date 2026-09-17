#ifndef MASKING_HPP_
#define MASKING_HPP_

// Screening/masking of conventions.md sec.10, milestone-1.md task 1.7:
// strings must be screened out of the field written into the FFT buffer
// before transforming. Applied at exactly one place -- this module -- not
// inlined at several call sites (CLAUDE.md constraint 5: an index slip in
// one of several inlined copies is what went unnoticed in the reference
// code). The threshold is a runtime parameter, never a compile-time
// constant (CLAUDE.md constraint 4), so it can be scanned.
//
// Two schemes, nested (conventions.md sec.14 decision log): scheme A's
// a_dot_scr is the bare numerator (psi1 Pi2 - Pi1 psi2); scheme B is the
// same divided by |psi|^2 (the genuine dOmega/dtau, i.e. the axion phase
// rate), with a top-hat that zeros it below a runtime threshold on
// |psi|/R. A third, "None" mode (f = 1 identically, no top-hat) is a
// first-class option: task 1.8 needs it for the *unscreened* pass the
// string tension comes from (the difference of the screened and
// unscreened energies).
//
// gamma(x) = arg(psi) = arg(phi) exactly (psi = R phi/v, R > 0 real
// everywhere -- sec.3), so, like the string finder, this operates directly
// on the stored state; no rescaling needed beyond dividing by R for the
// top-hat's |psi|/R.
//
// Falls back to plain host-only macros when AMReX headers haven't already
// defined these (compiled standalone for its unit test); picks up AMReX's
// real GPU-decorated versions when included from a kernel file -- same
// pattern as PlaquetteWinding.hpp.
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

enum class MaskingScheme
{
    None, // f = 1, no top-hat: the required unscreened pass (task 1.8)
    A,    // smooth: bare numerator, no division, no top-hat
    B     // top-hat: numerator / |psi|^2, zeroed below the threshold
};

struct MaskingParams
{
    MaskingScheme scheme{MaskingScheme::A};
    double threshold{0.8}; // scheme B only: top-hat on |psi|/R (runtime --
                           // to be scanned, not settled; conventions.md
                           // sec.10/sec.14 record 0.8, 0.9, 0.95 all having
                           // been tried without reconciliation)
};

// Whether this lattice point counts towards the effective point count used
// for unbiased spatial averaging (conventions.md sec.10: "the effective
// number of lattice points averaged over must be multiplied by the same
// masking factor"). 1 for None/A (no top-hat); 0 or 1 for B.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
masking_weight(const MaskingParams &params, double psi1, double psi2,
              double R)
{
    if (params.scheme != MaskingScheme::B)
    {
        return 1.0;
    }
    const double mod = std::sqrt(psi1 * psi1 + psi2 * psi2) / R;
    return (mod >= params.threshold) ? 1.0 : 0.0;
}

// The masked axion time-derivative written into the FFT buffer.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
masked_a_dot(const MaskingParams &params, double psi1, double psi2,
            double Pi1, double Pi2, double R)
{
    const double raw = psi1 * Pi2 - Pi1 * psi2;
    if (params.scheme == MaskingScheme::A)
    {
        return raw;
    }
    // Scheme B: check the top-hat *before* dividing, so an exactly-zero
    // core (psi = 0, a 0/0 in the raw/|psi|^2 below) is zeroed cleanly
    // rather than propagating a NaN through 0 * NaN.
    if (params.scheme == MaskingScheme::B &&
        masking_weight(params, psi1, psi2, R) == 0.0)
    {
        return 0.0;
    }
    // None and B (above threshold) both evaluate the genuine phase rate
    // raw/|psi|^2 (conventions.md sec.14: the fixed-grid code's
    // adt = .../|psi|^2).
    const double psi_sq = psi1 * psi1 + psi2 * psi2;
    return raw / psi_sq;
}

#endif // MASKING_HPP_
