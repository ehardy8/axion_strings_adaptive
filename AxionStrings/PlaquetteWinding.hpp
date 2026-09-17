#ifndef PLAQUETTEWINDING_HPP_
#define PLAQUETTEWINDING_HPP_

// Plaquette winding of docs/conventions.md sec.8: a plaquette is pierced if
// traversing its four vertices accumulates a net 2*pi in gamma(x)=arg(phi),
// with each consecutive difference reduced to (-pi, pi].
//
// Falls back to plain host-only macros when AMReX headers haven't already
// defined these (i.e. when this header is compiled standalone for
// tests/test_plaquette_winding.cpp), and picks up AMReX's real GPU-decorated
// versions automatically when they have (i.e. when included from
// StringFinder.hpp for the actual kernel). Keeps this the single place the
// winding logic is defined, in both contexts, per CLAUDE.md constraint 5's
// spirit (one place, not inlined at several call sites).
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double wrap_angle_diff(double delta)
{
    constexpr double pi     = 3.14159265358979323846;
    constexpr double two_pi = 2.0 * pi;
    delta                   = std::fmod(delta, two_pi);
    if (delta <= -pi)
    {
        delta += two_pi;
    }
    else if (delta > pi)
    {
        delta -= two_pi;
    }
    return delta;
}

// Corners passed in traversal order around the plaquette loop. Returns the
// (generically +-1, occasionally higher) winding number; 0 if not pierced.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
plaquette_winding(double theta1, double theta2, double theta3, double theta4)
{
    constexpr double two_pi = 2.0 * 3.14159265358979323846;
    const double total      = wrap_angle_diff(theta2 - theta1) +
                         wrap_angle_diff(theta3 - theta2) +
                         wrap_angle_diff(theta4 - theta3) +
                         wrap_angle_diff(theta1 - theta4);
    return static_cast<int>(std::lround(total / two_pi));
}

#endif // PLAQUETTEWINDING_HPP_
