#ifndef TANGENTFIELD_HPP_
#define TANGENTFIELD_HPP_

// The string's local tangent direction, estimated directly from the
// continuous field rather than from which lattice axis a pierced face
// happens to sit on (2026-09-26/27, with the user -- see
// CurvatureKernel.hpp's own header comment for why the face-position-only
// approach was replaced: on a well-resolved lattice, chords between
// adjacent pierced faces can only ever be dead straight or bent by
// exactly the lattice's own minimum turning angle, so any curvature built
// only from face *positions* is forced into two delta-function-like
// values, never a real distribution).
//
// T ~ grad(psi1) x grad(psi2), unit-normalised: this is (up to
// normalisation) the same topological current density whose discretised
// flux through a plaquette gives the winding number
// (StringFinder.hpp::plaquette_winding) -- but evaluated here from actual
// field gradients, so it can point in any direction, not just along a
// lattice axis. psi = R(tau) phi / v (R a spatially uniform scale factor
// at fixed time) so grad(psi) is parallel to grad(phi); the cross product
// of two such gradients only needs a direction, so R cancels in the
// unit-normalisation and does not need to be supplied here.
//
// Dual-purpose like MengerCurvature.hpp/PlaquetteWinding.hpp: falls back
// to plain host-only macros when compiled standalone for its unit test
// (tests/test_tangent_field.cpp), picks up AMReX's real GPU-decorated
// versions when included from CurvatureKernel.hpp.
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

struct TangentVector
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
    bool valid{false}; // false if the gradients were (numerically) parallel
                       // or vanishing, so no direction is defined here.
};

// grad_psi1 x grad_psi2, unit-normalised. Not a per-cell/per-face lookup --
// takes the two gradient vectors directly, so it is testable with plain
// numbers and reusable regardless of how the gradients were estimated
// (CurvatureKernel.hpp uses centred finite differences on the grid).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TangentVector
tangent_from_gradients(double dpsi1_dx, double dpsi1_dy, double dpsi1_dz,
                       double dpsi2_dx, double dpsi2_dy, double dpsi2_dz)
{
    const double tx = dpsi1_dy * dpsi2_dz - dpsi1_dz * dpsi2_dy;
    const double ty = dpsi1_dz * dpsi2_dx - dpsi1_dx * dpsi2_dz;
    const double tz = dpsi1_dx * dpsi2_dy - dpsi1_dy * dpsi2_dx;

    const double mag = std::sqrt(tx * tx + ty * ty + tz * tz);
    TangentVector out{};
    if (mag <= 0.0)
    {
        return out; // valid = false
    }
    out.x     = tx / mag;
    out.y     = ty / mag;
    out.z     = tz / mag;
    out.valid = true;
    return out;
}

#endif // TANGENTFIELD_HPP_
