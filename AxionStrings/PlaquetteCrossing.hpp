#ifndef PLAQUETTECROSSING_HPP_
#define PLAQUETTECROSSING_HPP_

// Sub-grid string position (2026-09-27, with the user): where a pierced
// plaquette's own bilinear interpolant of (psi1, psi2) actually vanishes,
// rather than just the plaquette's centre. This is the standard vortex-
// core-finding technique (bilinear interpolation of the two real field
// components, solved simultaneously for their common zero) -- an
// independent second way to build "the position of the string", to be
// used as a cross-check against TangentField.hpp's gradient-based
// tangent (CurvatureKernel.hpp switches between the two via
// axion_strings.curvature.method).
//
// Corners are always passed in the same traversal order as
// PlaquetteWinding.hpp's plaquette_winding (theta1..theta4, i.e. the
// (u,v) unit-square corners (0,0),(1,0),(1,1),(0,1)) -- the two
// components f (psi1) and g (psi2), each bilinear in (u,v) on that unit
// square:
//   f(u,v) = f00 + (f10-f00) u + (f01-f00) v + (f00-f10-f01+f11) uv
// (and likewise for g). Their common zero is found by eliminating v
// between f(u,v)=0 and g(u,v)=0, leaving a quadratic in u.
//
// A genuinely pierced plaquette (winding +-1) has, generically, exactly
// one such zero inside the unit square; 0 or 2 in-range roots are
// treated as "not found" here (ambiguous/degenerate) rather than guessed
// at -- CurvatureKernel.hpp counts these as excluded measurements, same
// discipline as its other topology checks.
//
// Dual-purpose like MengerCurvature.hpp/TangentField.hpp/
// PlaquetteWinding.hpp: falls back to plain host-only macros when
// compiled standalone for its unit test
// (tests/test_plaquette_crossing.cpp).
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

struct BilinearZero
{
    double u{0.0};
    double v{0.0};
    bool valid{false};
};

// f/g corners in (00,10,11,01) order, matching plaquette_winding's own
// traversal convention exactly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE BilinearZero
find_bilinear_zero(double f00, double f10, double f11, double f01,
                   double g00, double g10, double g11, double g01)
{
    const double a0 = f00;
    const double a1 = f10 - f00;
    const double a2 = f01 - f00;
    const double a3 = f00 - f10 - f01 + f11;

    const double b0 = g00;
    const double b1 = g10 - g00;
    const double b2 = g01 - g00;
    const double b3 = g00 - g10 - g01 + g11;

    // Eliminating v between f(u,v)=0 and g(u,v)=0 leaves A u^2+B u+C = 0.
    const double A = b1 * a3 - b3 * a1;
    const double B = b0 * a3 + b1 * a2 - b2 * a1 - b3 * a0;
    const double C = b0 * a2 - b2 * a0;

    constexpr double eps = 1.0e-14;
    double roots[2];
    int n_roots = 0;
    if (std::abs(A) < eps)
    {
        if (std::abs(B) < eps)
        {
            return BilinearZero{}; // f,g proportional in u -- no isolated
                                   // solution from this elimination.
        }
        roots[0] = -C / B;
        n_roots  = 1;
    }
    else
    {
        const double disc = B * B - 4.0 * A * C;
        if (disc < 0.0)
        {
            return BilinearZero{};
        }
        const double sq = std::sqrt(disc);
        roots[0]        = (-B + sq) / (2.0 * A);
        roots[1]        = (-B - sq) / (2.0 * A);
        n_roots         = 2;
    }

    BilinearZero found{};
    int n_accepted = 0;
    for (int r = 0; r < n_roots; ++r)
    {
        const double u = roots[r];
        if (u < 0.0 || u > 1.0)
        {
            continue;
        }
        // v from whichever of f=0/g=0 (as a linear equation in v at this
        // u) is better conditioned, to avoid amplifying error near a
        // small denominator.
        const double den_f = a2 + a3 * u;
        const double den_g = b2 + b3 * u;
        double v{};
        if (std::abs(den_f) >= std::abs(den_g))
        {
            if (std::abs(den_f) < eps)
            {
                continue;
            }
            v = -(a0 + a1 * u) / den_f;
        }
        else
        {
            if (std::abs(den_g) < eps)
            {
                continue;
            }
            v = -(b0 + b1 * u) / den_g;
        }
        if (v < 0.0 || v > 1.0)
        {
            continue;
        }
        found.u     = u;
        found.v     = v;
        found.valid = true;
        ++n_accepted;
    }
    if (n_accepted != 1)
    {
        return BilinearZero{}; // 0 -> no crossing found; >1 -> ambiguous.
    }
    return found;
}

#endif // PLAQUETTECROSSING_HPP_
