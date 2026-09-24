#ifndef HESSIANCURVATURE_HPP_
#define HESSIANCURVATURE_HPP_

// Third, structurally independent curvature method (2026-09-27, with the
// user): the analytic curvature of the implicitly-defined curve
// {psi1(x)=0} intersect {psi2(x)=0}, from the gradients and Hessians of
// psi1, psi2 at a single point -- no connectivity graph, no chain of
// neighbouring measurements, no chord/finite-baseline of any kind. Where
// TangentField.hpp needs a second point to form |dT/ds|, and
// PlaquetteCrossing.hpp needs a chain of three interpolated positions,
// this needs only local derivatives at the one point being measured.
//
// Derivation. Let x(s) be an arc-length parametrisation of the curve, so
// psi1(x(s)) = psi2(x(s)) = 0 identically. Differentiating once:
//   grad(psi_a) . T = 0        (a = 1,2; T = dx/ds, the unit tangent)
// so T is orthogonal to both gradients -- consistent with
// TangentField.hpp's T ~ grad(psi1) x grad(psi2), and grad(psi1),
// grad(psi2) span the 2D plane normal to T (assuming they are not
// parallel, i.e. the curve is genuinely 1D there). Differentiating again:
//   (Hess(psi_a) T) . T + grad(psi_a) . T' = 0,   T' = dT/ds = kappa*N
// so, writing the curvature vector kappa*N in that normal plane as
// kappa*N = a*grad(psi1) + b*grad(psi2), the two equations above become a
// 2x2 linear system in (a, b):
//   [ g11  g12 ] [a]   [ -T.H1.T ]
//   [ g12  g22 ] [b] = [ -T.H2.T ]
// with g_ab = grad(psi_a).grad(psi_b), H_a = Hess(psi_a). Solving and
// forming kappa*N = a*grad(psi1) + b*grad(psi2) gives kappa = |kappa*N|.
//
// Verified analytically (tests/test_hessian_curvature.cpp) against a
// circle of radius R written as an explicit intersection of two level
// sets (psi1 = y, psi2 = x^2 + (z-R)^2 - R^2, both zero exactly on a
// circle of radius R through the origin, tangent to the z-axis there) --
// this formula gives exactly kappa = 1/R at the origin, evaluated purely
// from psi1, psi2's own analytic gradients/Hessians with no reference to
// R itself being "the curve's radius" anywhere in the code. Also scale-
// invariant under (psi1,psi2) -> c*(psi1,psi2) for any nonzero constant c
// (both g_ab and T.H.T scale as c^2 and c respectively in a way that
// cancels exactly in (a,b), and a,b then scale as 1/c against the grad
// terms) -- so, like TangentField.hpp, it does not matter whether this is
// fed psi = R phi or phi directly.
//
// Dual-purpose like TangentField.hpp/MengerCurvature.hpp/
// PlaquetteWinding.hpp: falls back to plain host-only macros when
// compiled standalone for its unit test
// (tests/test_hessian_curvature.cpp).
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

#include <cmath>

struct HessianCurvatureResult
{
    double kappa{0.0};
    bool valid{false}; // false if grad(psi1), grad(psi2) are (numerically)
                       // parallel or vanishing -- no 1D curve/no normal
                       // plane defined here, same degeneracy
                       // TangentField.hpp guards against.
};

// h1xx..h1xz / h2xx..h2xz are the independent entries of the (symmetric)
// Hessians of psi1, psi2: xx, yy, zz, xy, yz, xz.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE HessianCurvatureResult
hessian_curvature(double g1x, double g1y, double g1z, double g2x,
                  double g2y, double g2z, double h1xx, double h1yy,
                  double h1zz, double h1xy, double h1yz, double h1xz,
                  double h2xx, double h2yy, double h2zz, double h2xy,
                  double h2yz, double h2xz)
{
    const double tx0 = g1y * g2z - g1z * g2y;
    const double ty0 = g1z * g2x - g1x * g2z;
    const double tz0 = g1x * g2y - g1y * g2x;
    const double tmag = std::sqrt(tx0 * tx0 + ty0 * ty0 + tz0 * tz0);
    if (tmag <= 0.0)
    {
        return HessianCurvatureResult{};
    }
    const double Tx = tx0 / tmag;
    const double Ty = ty0 / tmag;
    const double Tz = tz0 / tmag;

    auto quad_form = [&](double hxx, double hyy, double hzz, double hxy,
                         double hyz, double hxz)
    {
        return Tx * Tx * hxx + Ty * Ty * hyy + Tz * Tz * hzz +
              2.0 * Tx * Ty * hxy + 2.0 * Ty * Tz * hyz +
              2.0 * Tx * Tz * hxz;
    };
    const double Tt_H1_T =
        quad_form(h1xx, h1yy, h1zz, h1xy, h1yz, h1xz);
    const double Tt_H2_T =
        quad_form(h2xx, h2yy, h2zz, h2xy, h2yz, h2xz);

    const double g11 = g1x * g1x + g1y * g1y + g1z * g1z;
    const double g22 = g2x * g2x + g2y * g2y + g2z * g2z;
    const double g12 = g1x * g2x + g1y * g2y + g1z * g2z;
    const double det = g11 * g22 - g12 * g12;

    constexpr double eps = 1.0e-14;
    if (std::abs(det) < eps)
    {
        return HessianCurvatureResult{};
    }

    const double rhs1 = -Tt_H1_T;
    const double rhs2 = -Tt_H2_T;
    const double a = (rhs1 * g22 - rhs2 * g12) / det;
    const double b = (rhs2 * g11 - rhs1 * g12) / det;

    const double Kx = a * g1x + b * g2x;
    const double Ky = a * g1y + b * g2y;
    const double Kz = a * g1z + b * g2z;

    HessianCurvatureResult out{};
    out.kappa = std::sqrt(Kx * Kx + Ky * Ky + Kz * Kz);
    out.valid = true;
    return out;
}

#endif // HESSIANCURVATURE_HPP_
