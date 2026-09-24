#ifndef MENGERCURVATURE_HPP_
#define MENGERCURVATURE_HPP_

// Three-point Menger (circumradius) curvature: kappa = 4*Area/(d12 d23
// d13), coordinate-free, well-defined for any three non-collinear points,
// equal to 1/R for three points on a circle of radius R.
//
// Re-added 2026-09-27 (with the user) for the "interpolated position"
// curvature method (CurvatureKernel.hpp/PlaquetteCrossing.hpp): the first
// version of the curvature diagnostic used this on raw face *centres* and
// was bimodal by construction (see CurvatureKernel.hpp's header comment
// for the full story) -- fed sub-grid-interpolated positions instead, it
// is a genuine second, independent way to build "the position of the
// string" for cross-checking against the gradient-based tangent method.
//
// Dual-purpose like PlaquetteWinding.hpp: falls back to plain host-only
// macros when compiled standalone for its unit test
// (tests/test_menger_curvature.cpp), picks up AMReX's real GPU-decorated
// versions when included from CurvatureKernel.hpp -- one definition, not
// inlined at two call sites (CLAUDE.md constraint 5's spirit).
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

// 0 if the three points are degenerate (coincident, or one pairwise
// distance is exactly zero), rather than dividing by zero.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
menger_curvature(double x1, double y1, double z1, double x2, double y2,
                 double z2, double x3, double y3, double z3)
{
    const double d12 =
        std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) +
                 (z2 - z1) * (z2 - z1));
    const double d23 =
        std::sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) +
                 (z3 - z2) * (z3 - z2));
    const double d13 =
        std::sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) +
                 (z3 - z1) * (z3 - z1));
    if (d12 <= 0.0 || d23 <= 0.0 || d13 <= 0.0)
    {
        return 0.0;
    }
    // Area from the cross product of two edge vectors, magnitude/2.
    const double ux = x2 - x1;
    const double uy = y2 - y1;
    const double uz = z2 - z1;
    const double vx = x3 - x1;
    const double vy = y3 - y1;
    const double vz = z3 - z1;
    const double cx = uy * vz - uz * vy;
    const double cy = uz * vx - ux * vz;
    const double cz = ux * vy - uy * vx;
    const double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
    return 4.0 * area / (d12 * d23 * d13);
}

#endif // MENGERCURVATURE_HPP_
