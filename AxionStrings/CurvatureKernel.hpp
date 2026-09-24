#ifndef CURVATUREKERNEL_HPP_
#define CURVATUREKERNEL_HPP_

// Local string curvature (2026-09-26/27, with the user; composite/multi-
// level extension, Stage 2, 2026-09-27): a distribution of curvatures
// along the network, binned by kappa/m_r (dimensionless -- curvature
// relative to the string's own core scale, so bins are directly
// comparable across resolution and time), saved per snapshot. Operates
// across the whole AMR hierarchy (AxionStringsLevel.cpp's composite
// driver calls compute_curvature_histogram once per level with that
// level's own coverage mask, and merges the results) -- not just level 0,
// so refined regions actually contribute their own finer-scale bends
// instead of being silently measured at the coarsest available
// resolution.
//
// The connectivity (the user's own observation, 2026-09-25): a pierced
// face is shared by exactly two cells, and each of those cells has --
// generically -- exactly one *other* pierced face, since a string cannot
// simply end inside a cell (conservation of the winding number: net
// piercing flux through a cell's boundary is zero unless a string
// endpoint sits inside it, which does not happen for closed loops or a
// percolating network). So a cell with exactly two pierced faces among
// its six has a well-defined local "chord" connecting those two face
// centres, and each pierced face's two neighbours along the string (one
// via each of the two cells sharing it) fall straight out of this, with
// no explicit tracing or connectivity graph needed -- every lookup is a
// fixed, local read of a neighbouring cell's own six faces.
//
// Curvature itself is *not* computed from these face positions (a first
// version did, via three-point Menger curvature -- see git history for
// MengerCurvature.hpp, removed here). On a well-resolved lattice
// (m_r*dx << 1, needed for correctness anyway) a chord between two
// adjacent pierced faces is forced into exactly one of two shapes: dead
// straight (string enters/exits a cell on opposite faces -> the three
// points are *exactly* collinear -> kappa = 0 exactly), or a single-cell
// turn (adjacent faces of the same cell -> kappa pinned at ~1/dx,
// regardless of the network's actual physical bend). A live 15-minute
// network test confirmed this exhaustively: 0 of 720 histogram bins
// across 30 snapshots ever got a single count -- literally every
// measurement fell in "below range" (kappa exactly 0) or "above range"
// (kappa/m_r clustering at 1/(m_r*dx) = 16 for this test's resolution).
// This is not a resolution/bin-range problem; it is inherent to building
// curvature only out of *which lattice axis* a face sits on.
//
// Fix, two independent ways (the user's own suggestions, refined
// together, 2026-09-27) -- kept as a runtime-switchable
// axion_strings.curvature.method, not a single fixed choice, precisely
// *because* this is subtle enough that the two should be cross-checked
// against each other rather than trusted on their own:
//
// "tangent_vector": estimate the string's local tangent direction from
// the *continuous* field instead of face positions, via TangentField.hpp's
// grad(psi1) x grad(psi2) (the same topological current whose discretised
// flux through a plaquette gives the winding number,
// StringFinder.hpp::plaquette_winding -- but from actual centred finite
// differences, so it can point in any direction, not just along a lattice
// axis). Curvature is then the discrete Frenet-Serret |dT/ds|: evaluate
// the unit tangent at the two cells straddling a pierced face (cell A =
// the face's canonical owner, cell B = the neighbour across it -- exactly
// the two cells the connectivity above already identifies) and take
//   kappa = |T_B - T_A| / (arc length between where T_A and T_B are
//                          themselves sampled)
// T_A and T_B are each evaluated from centred differences at a cell
// centre, i.e. at (approximately) the midpoint of their own chord -- T_A
// near the midpoint of P1->P2 (the owner cell's chord), T_B near the
// midpoint of P2->P3 (the neighbour's chord). Their true separation is
// therefore *one* chord length, not the full P1->P2->P3 arc: an earlier
// version divided by the two-chord arc length and underestimated kappa by
// a factor of 2 (caught 2026-09-24 via a sister project's port of this
// method, confirmed here both algebraically and against an exact circle --
// see the fix note at the call site). The three face *positions* are still
// needed, only to get that one-chord arc length (P1 = the owner cell's
// other pierced face, P2 = the shared face, P3 = the neighbour's other
// pierced face) -- T itself does not depend on face positions at all, so
// it is free to vary continuously even between immediately-adjacent
// pierced faces.
//
// "interpolated_position": instead of fixing the bimodality by dropping
// face positions, fix the positions themselves -- PlaquetteCrossing.hpp
// bilinearly interpolates (psi1, psi2) across each pierced plaquette and
// solves for their common zero, the standard sub-grid vortex-core-finding
// technique, giving a genuine (non-lattice-quantised) position for the
// string on that face. Feed those three interpolated positions (P1, P2,
// P3, same three faces as above) into the original three-point Menger
// curvature (MengerCurvature.hpp) -- the same formula the first, buggy
// version used on face *centres*, now on positions that can fall anywhere
// within the plaquette rather than only at its middle.
//
// These are two different constructions of the *same* quantity (a
// tangent-turning-rate vs a circumradius through interpolated positions),
// so agreement between them is real cross-validation, not just running the
// same computation twice -- with tangent_vector's normalisation fixed
// (2026-09-24, see above), they do agree; see the note near the end of
// this comment. A spot check reconstructing actual chains of interpolated
// positions
// (CurvatureChainDebug.hpp) confirmed both track visibly-real bends and
// straight stretches, but also turned up a genuine caveat specific to
// interpolated_position: consecutive interpolated crossings can land
// almost coincident when the string passes close to a cell edge, and the
// Menger circumradius blows up on that near-degenerate triangle -- the
// top outliers in a real run's tail all had a chord length a few percent
// of dx. Left unguarded deliberately (2026-09-27, with the user) rather
// than adding a minimum-chord-length exclusion, so this is a known,
// documented property of that method's tail, not a silently-patched one.
//
// "hessian_analytic": a third, structurally different construction
// (HessianCurvature.hpp) -- the analytic curvature of the curve
// {psi1=0} intersect {psi2=0} from the gradients and Hessians of psi1,
// psi2 at a *single point*, with no chain, no chord, no second
// measurement at all (verified against a circle of known radius
// analytically, see that file's own header comment). Evaluated once per
// qualifying *cell* (not per face -- there is no face-specific
// information in this formula at all). This turns out to give almost the
// *same* measurement count as the other two methods, not half of it as a
// naive count might suggest: get_cell_piercing's count, summed over every
// cell, double-counts each distinct pierced face exactly once from each
// of its two adjacent cells, so (since almost every nonzero-count cell
// has count exactly 2, the generic case) the number of qualifying cells
// is itself ~N_p -- the same total the other two methods reach via one
// canonical-owner measurement per distinct face. Confirmed directly on
// real network data: 51104 cell-based measurements vs 51722/49823
// face-based, at N_p=51802, all within a few percent of each other.
//
// A live comparison (2026-09-27, with the user) across a full network
// test found hessian_analytic tracking interpolated_position closely --
// mean log10(kappa/m_r) within ~0.02-0.05 at nearly every snapshot, only
// drifting to ~0.1 late in the run as statistics thin out -- while
// tangent_vector sat on its own, consistently offset ~0.3-0.4 dex lower
// throughout. At the time this was read as evidence that tangent_vector
// measures a genuinely different quantity. It does not: the offset was a
// plain factor-of-2 normalisation bug (kappa = |dT| / arc_length instead
// of |dT| / (arc_length/2), see the fix note above and at the
// TangentVector call site, fixed 2026-09-24) -- log10(2) = 0.30 dex,
// exactly the low end of the observed offset, and the remainder is
// consistent with the method's own discretisation error relative to the
// other two. Once fixed, all three are expected to (and should be
// re-checked to) track each other; the earlier "genuinely different
// observables, don't expect agreement" framing above should be read in
// that light. Two structurally independent constructions (one from
// sub-grid positions, one from purely local derivatives) landing on
// essentially the same number is real cross-validation, not a
// coincidence of both being "the same formula in disguise".
//
// Deduplication (tangent_vector/interpolated_position only --
// hessian_analytic measures once per cell, so this does not apply to
// it): each pierced face is shared by two cells, and a naive "check
// every cell's two faces" pass would compute (and bin) every face twice,
// from either side, with the identical value -- wasteful, and wrong for
// a distribution. Avoided via the *same* low-index-corner convention
// already used everywhere else in this project for exactly this reason
// (StringFinder.hpp's own docstring): a FaceId's own (i,j,k) is *always*
// its canonical low-index-owner cell, by construction, so a face is only
// ever measured when the cell currently being processed equals that
// owner -- deferred to the appropriate neighbour otherwise, never both.
//
// `state` must already have at least 3 valid ghost cells (state.
// FillBoundary) before calling this: computing a single cell's own
// six-face piercing needs plaquette_windings_at at up to (i+1,j+1,k+1)
// relative to it (StringFinder.hpp's own 7-corner stencil), and finding a
// pierced face's *other* neighbour needs that same six-face computation
// repeated one cell further out -- matching StringTagger.hpp's own
// stated ghost requirement, already the project's standing default. The
// tangent-vector evaluation only needs +-1 relative to each of the two
// cells it is evaluated at, comfortably within the same 3-ghost-cell
// budget.

#include "FourthOrderDerivatives.hpp"
#include "HessianCurvature.hpp"
#include "MengerCurvature.hpp"
#include "PlaquetteCrossing.hpp"
#include "StringFinder.hpp"
#include "TangentField.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <vector>

// A face's canonical identity: axis = 0,1,2 for the x-normal (yz), y-normal
// (zx), z-normal (xy) orientations respectively (matching StringFinder.hpp's
// w_yz/w_zx/w_xy ordering); (i,j,k) is *always* that face's low-index-owner
// cell (StringFinder.hpp/StringTagger.hpp's own convention), regardless of
// which cell is asking about it.
struct StringFaceId
{
    int axis{0};
    int i{0};
    int j{0};
    int k{0};

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator==(const StringFaceId &other) const
    {
        return axis == other.axis && i == other.i && j == other.j &&
              k == other.k;
    }
};

// The (up to) two pierced faces of a single cell, plus the total count
// actually found (checking all six, not stopping early at two, so a
// pathological >2 case is correctly counted rather than silently
// truncated).
struct CellPiercing
{
    int count{0};
    StringFaceId face[2]{};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE CellPiercing
get_cell_piercing(const amrex::Array4<amrex::Real const> &state, int i, int j,
                  int k)
{
    const auto low  = plaquette_windings_at(state, i, j, k);
    const auto hix  = plaquette_windings_at(state, i + 1, j, k).w_yz;
    const auto hiy  = plaquette_windings_at(state, i, j + 1, k).w_zx;
    const auto hiz  = plaquette_windings_at(state, i, j, k + 1).w_xy;

    const bool nonzero[6] = {low.w_yz != 0, low.w_zx != 0, low.w_xy != 0,
                            hix != 0,       hiy != 0,       hiz != 0};
    const StringFaceId ids[6] = {{0, i, j, k},     {1, i, j, k},
                                {2, i, j, k},     {0, i + 1, j, k},
                                {1, i, j + 1, k}, {2, i, j, k + 1}};

    CellPiercing result{};
    for (int m = 0; m < 6; ++m)
    {
        if (nonzero[m])
        {
            if (result.count < 2)
            {
                result.face[result.count] = ids[m];
            }
            ++result.count;
        }
    }
    return result;
}

// Physical (well, comoving -- this operates entirely in comoving grid
// units, same convention as every other tagging/finding kernel here)
// coordinates of a face's own centre.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_center(const StringFaceId &f, const amrex::GpuArray<amrex::Real, 3> &dx,
           const amrex::GpuArray<amrex::Real, 3> &prob_lo, double &x,
           double &y, double &z)
{
    x = prob_lo[0] + dx[0] * (f.i + (f.axis == 0 ? 0.0 : 0.5));
    y = prob_lo[1] + dx[1] * (f.j + (f.axis == 1 ? 0.0 : 0.5));
    z = prob_lo[2] + dx[2] * (f.k + (f.axis == 2 ? 0.0 : 0.5));
}

// The cell on the *other* side of a face from its own canonical owner
// (owner - 1 along the face's own axis, matching how a face's "high" side
// is the neighbour's own "low" side everywhere else in this file).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
other_cell_across_face(const StringFaceId &f, int &oi, int &oj, int &ok)
{
    oi = f.i;
    oj = f.j;
    ok = f.k;
    if (f.axis == 0)
    {
        oi -= 1;
    }
    else if (f.axis == 1)
    {
        oj -= 1;
    }
    else
    {
        ok -= 1;
    }
}

// The string's local tangent direction at a cell, from centred finite
// differences of the field itself (TangentField.hpp's header comment) --
// unlike a face position, this can point in any direction, not just along
// a lattice axis. Needs +-1 valid cell around (i,j,k) in every direction.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TangentVector
tangent_field_at_cell(const amrex::Array4<amrex::Real const> &state,
                      const amrex::GpuArray<amrex::Real, 3> &dx, int i,
                      int j, int k)
{
    const double dpsi1_dx = (static_cast<double>(state(i + 1, j, k, c_psi1)) -
                             static_cast<double>(state(i - 1, j, k, c_psi1))) /
                            (2.0 * dx[0]);
    const double dpsi1_dy = (static_cast<double>(state(i, j + 1, k, c_psi1)) -
                             static_cast<double>(state(i, j - 1, k, c_psi1))) /
                            (2.0 * dx[1]);
    const double dpsi1_dz = (static_cast<double>(state(i, j, k + 1, c_psi1)) -
                             static_cast<double>(state(i, j, k - 1, c_psi1))) /
                            (2.0 * dx[2]);
    const double dpsi2_dx = (static_cast<double>(state(i + 1, j, k, c_psi2)) -
                             static_cast<double>(state(i - 1, j, k, c_psi2))) /
                            (2.0 * dx[0]);
    const double dpsi2_dy = (static_cast<double>(state(i, j + 1, k, c_psi2)) -
                             static_cast<double>(state(i, j - 1, k, c_psi2))) /
                            (2.0 * dx[1]);
    const double dpsi2_dz = (static_cast<double>(state(i, j, k + 1, c_psi2)) -
                             static_cast<double>(state(i, j, k - 1, c_psi2))) /
                            (2.0 * dx[2]);
    return tangent_from_gradients(dpsi1_dx, dpsi1_dy, dpsi1_dz, dpsi2_dx,
                                  dpsi2_dy, dpsi2_dz);
}

// Sub-grid position of a pierced face: bilinearly interpolates (psi1,
// psi2) across the plaquette and finds their common zero
// (PlaquetteCrossing.hpp), using each axis's own corner order (matching
// plaquette_windings_at's w_yz/w_zx/w_xy conventions exactly, so this
// finds a zero of the *same* plaquette that was tested for piercing).
// Returns false (leaving x,y,z untouched) if no unique interior zero was
// found (PlaquetteCrossing.hpp's own header comment).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
interpolated_crossing_position(const amrex::Array4<amrex::Real const> &state,
                               const amrex::GpuArray<amrex::Real, 3> &dx,
                               const amrex::GpuArray<amrex::Real, 3> &prob_lo,
                               const StringFaceId &f, double &x, double &y,
                               double &z)
{
    const int i = f.i;
    const int j = f.j;
    const int k = f.k;

    auto psi1 = [&](int ii, int jj, int kk)
    { return static_cast<double>(state(ii, jj, kk, c_psi1)); };
    auto psi2 = [&](int ii, int jj, int kk)
    { return static_cast<double>(state(ii, jj, kk, c_psi2)); };

    BilinearZero zero{};
    if (f.axis == 0) // w_yz: corners (j,k)->(j+1,k)->(j+1,k+1)->(j,k+1)
    {
        zero = find_bilinear_zero(
            psi1(i, j, k), psi1(i, j + 1, k), psi1(i, j + 1, k + 1),
            psi1(i, j, k + 1), psi2(i, j, k), psi2(i, j + 1, k),
            psi2(i, j + 1, k + 1), psi2(i, j, k + 1));
        if (!zero.valid)
        {
            return false;
        }
        x = prob_lo[0] + dx[0] * i;
        y = prob_lo[1] + dx[1] * (j + zero.u);
        z = prob_lo[2] + dx[2] * (k + zero.v);
        return true;
    }
    if (f.axis == 1) // w_zx: corners (k,i)->(k+1,i)->(k+1,i+1)->(k,i+1)
    {
        zero = find_bilinear_zero(
            psi1(i, j, k), psi1(i, j, k + 1), psi1(i + 1, j, k + 1),
            psi1(i + 1, j, k), psi2(i, j, k), psi2(i, j, k + 1),
            psi2(i + 1, j, k + 1), psi2(i + 1, j, k));
        if (!zero.valid)
        {
            return false;
        }
        x = prob_lo[0] + dx[0] * (i + zero.v);
        y = prob_lo[1] + dx[1] * j;
        z = prob_lo[2] + dx[2] * (k + zero.u);
        return true;
    }
    // axis == 2, w_xy: corners (i,j)->(i+1,j)->(i+1,j+1)->(i,j+1)
    zero = find_bilinear_zero(
        psi1(i, j, k), psi1(i + 1, j, k), psi1(i + 1, j + 1, k),
        psi1(i, j + 1, k), psi2(i, j, k), psi2(i + 1, j, k),
        psi2(i + 1, j + 1, k), psi2(i, j + 1, k));
    if (!zero.valid)
    {
        return false;
    }
    x = prob_lo[0] + dx[0] * (i + zero.u);
    y = prob_lo[1] + dx[1] * (j + zero.v);
    z = prob_lo[2] + dx[2] * k;
    return true;
}

// Analytic single-point curvature (HessianCurvature.hpp) at a cell
// centre, from psi1, psi2's gradients *and* Hessians. 4th-order accurate
// (2026-09-27, with the user -- upgraded from an initial hand-rolled
// 2nd-order version): reuses GRTeclyn's own FourthOrderDerivatives
// (Source/Grids/FourthOrderDerivatives.hpp, already the project's
// standard stencil everywhere else -- CLAUDE.md's "KleinGordon is the
// template") rather than re-deriving the mixed-partial stencil by hand;
// its d2_scalar's mixed_diff2 is exactly the tensor product of two
// 4th-order 1D first-derivative stencils, the standard construction for
// a 4th-order mixed second derivative. Since this method has no chain or
// chord (unlike tangent_vector/interpolated_position), the stencil
// accuracy directly *is* the method's accuracy -- there is no other
// discretisation error to be dominated by. Needs +-2 valid cells around
// (i,j,k) in every direction (comfortably inside this file's standing
// 3-ghost-cell requirement); assumes an isotropic grid (dx[0]==dx[1]==
// dx[2]), true everywhere in this project (geometry.coarsest_dx is a
// single value) -- FourthOrderDerivatives itself only takes one dx.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE HessianCurvatureResult
hessian_curvature_at_cell(const amrex::Array4<amrex::Real const> &state,
                          const amrex::GpuArray<amrex::Real, 3> &dx, int i,
                          int j, int k)
{
    const FourthOrderDerivatives derivs(dx[0]);
    const Tensor::Rank1 grad1     = derivs.d1_scalar(i, j, k, state, c_psi1);
    const Tensor::Rank1 grad2     = derivs.d1_scalar(i, j, k, state, c_psi2);
    const Tensor::Sym12Rank2 hess1 = derivs.d2_scalar(i, j, k, state, c_psi1);
    const Tensor::Sym12Rank2 hess2 = derivs.d2_scalar(i, j, k, state, c_psi2);

    return hessian_curvature(
        grad1(0), grad1(1), grad1(2), grad2(0), grad2(1), grad2(2),
        hess1(0, 0), hess1(1, 1), hess1(2, 2), hess1(0, 1), hess1(1, 2),
        hess1(0, 2), hess2(0, 0), hess2(1, 1), hess2(2, 2), hess2(0, 1),
        hess2(1, 2), hess2(0, 2));
}

enum class CurvatureMethod
{
    TangentVector,       // TangentField.hpp: |dT/ds| from field gradients
    InterpolatedPosition, // MengerCurvature.hpp on PlaquetteCrossing.hpp's
                         // sub-grid zero-crossing positions
    HessianAnalytic      // HessianCurvature.hpp: single-point analytic
                         // curvature from gradients + Hessians
};

struct CurvatureParams
{
    int n_bins{24};
    double log10_kappa_over_mr_min{-3.0};
    double log10_kappa_over_mr_max{1.0};
    CurvatureMethod method{CurvatureMethod::TangentVector};
};

struct CurvatureHistogram
{
    std::vector<double> bin_center_log10;   // log10(kappa/m_r) at bin centre
    std::vector<double> plain_count;        // one measurement = one count
    std::vector<double> length_weighted;    // weighted by (d12+d23)/2

    long n_measurements{0};        // accepted, binned (in-range) measurements
    long n_below_range{0};         // accepted but below the histogram's range
    long n_above_range{0};         // accepted but above the histogram's range
    long n_excluded_not_two{0};    // this cell itself has != 2 pierced faces
    long n_excluded_neighbor{0};   // the neighbour needed for the 3rd point
                                   // itself has != 2 pierced faces
    long n_excluded_topology{0};   // neighbour's own faces did not include
                                   // the shared face -- should not happen if
                                   // plaquette_windings_at is self-consistent
                                   // (see this file's header comment); kept
                                   // as a defensive check, not assumed away
    long n_excluded_degenerate_tangent{0}; // tangent_vector method only:
                                           // grad(psi1) and grad(psi2) were
                                           // (numerically) parallel or
                                           // vanishing at cell A or cell B,
                                           // so no tangent direction is
                                           // defined there -- expected to be
                                           // rare, since a pierced face
                                           // implies a genuine phase winding
                                           // nearby.
    long n_excluded_no_crossing{0}; // interpolated_position method only:
                                    // PlaquetteCrossing.hpp found 0 or >1
                                    // interior zeros on one of the three
                                    // faces -- expected to be rare for the
                                    // same reason as above.
    long n_excluded_degenerate_hessian{0}; // hessian_analytic method
                                           // only: grad(psi1), grad(psi2)
                                           // (numerically) parallel or
                                           // vanishing at the cell centre.
};

// Per-level driver (2026-09-27, with the user: composite/multi-level
// extension, Stage 2 of the plan agreed when this diagnostic was first
// designed). `state` is *some* level's own state -- level 0's for a
// single-level run, or any one level of a composite hierarchy -- already
// ghost-filled with at least 3 valid ghost cells (a plain FillBoundary
// for a single level; AxionStringsLevel.cpp's gather_composite_levels
// for a real hierarchy, whose composite FillPatch correctly interpolates
// ghost cells across a coarse-fine boundary rather than leaving them
// stale -- the same concern StringFinder.hpp's own composite path
// already guards against). `m_r` normalises the histogram's bin axis to
// the dimensionless kappa/m_r.
//
// `mask`, when non-null (StringFinder.hpp::count_plaquettes's own
// convention exactly), is a per-cell coverage mask on `state`'s own
// BoxArray/DistributionMapping (1 = include, 0 = skip): a cell's own
// piercing is skipped entirely when masked, deferring that physical
// region to whichever finer level actually covers it, so the same face
// is never measured at two resolutions. The mask is checked *only* at
// each cell the outer loop is centred on (never at a neighbour reached
// during that cell's own chain/gradient lookups) -- sufficient to
// prevent double-counting (every face's canonical owner is a single,
// specific cell, and that owner's own mask value is what decides whether
// this level claims it), and it avoids needing to index the mask
// iMultiFab (itself built with zero ghost cells, unlike `state`) at an
// arbitrary neighbour that may sit in a different box or ghost region --
// exactly the granularity StringFinder.hpp/VelocityKernel.hpp's own mask
// checks already use, since neither has a neighbour concept either.
// AxionStringsLevel.cpp's composite driver calls this once per level and
// sums the resulting histograms; for a single-level run (mask always
// null) it reduces to exactly the original level-0-only behaviour.
[[nodiscard]] inline CurvatureHistogram
compute_curvature_histogram(const amrex::MultiFab &state,
                            const amrex::Geometry &geom, double m_r,
                            const CurvatureParams &params,
                            const amrex::iMultiFab *mask = nullptr)
{
    CurvatureHistogram hist{};
    hist.bin_center_log10.resize(static_cast<std::size_t>(params.n_bins));
    hist.plain_count.assign(static_cast<std::size_t>(params.n_bins), 0.0);
    hist.length_weighted.assign(static_cast<std::size_t>(params.n_bins), 0.0);
    const double bin_width =
        (params.log10_kappa_over_mr_max - params.log10_kappa_over_mr_min) /
        static_cast<double>(params.n_bins);
    for (int b = 0; b < params.n_bins; ++b)
    {
        hist.bin_center_log10[static_cast<std::size_t>(b)] =
            params.log10_kappa_over_mr_min + (b + 0.5) * bin_width;
    }

    const auto dx      = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    long n_measurements = 0, n_below = 0, n_above = 0, n_not_two = 0,
        n_neighbor = 0, n_topology = 0, n_degenerate_tangent = 0,
        n_no_crossing = 0, n_degenerate_hessian = 0;
    std::vector<double> plain_local(static_cast<std::size_t>(params.n_bins),
                                    0.0);
    std::vector<double> weighted_local(
        static_cast<std::size_t>(params.n_bins), 0.0);

    // Shared by all three methods: log-bin kappa/m_r (or count it as
    // below/above range), given the per-measurement weight already
    // computed by whichever method is active.
    auto bin_measurement = [&](double kappa_over_mr, double weight)
    {
        if (kappa_over_mr <= 0.0)
        {
            // A perfectly straight local segment (kappa=0) has no
            // well-defined log; not an error, just outside a log-spaced
            // histogram's domain.
            ++n_below;
            return;
        }
        const double log10_k = std::log10(kappa_over_mr);
        if (log10_k < params.log10_kappa_over_mr_min)
        {
            ++n_below;
            return;
        }
        if (log10_k >= params.log10_kappa_over_mr_max)
        {
            ++n_above;
            return;
        }
        int bin = static_cast<int>(
            (log10_k - params.log10_kappa_over_mr_min) / bin_width);
        bin = amrex::Clamp(bin, 0, params.n_bins - 1);
        plain_local[static_cast<std::size_t>(bin)] += 1.0;
        weighted_local[static_cast<std::size_t>(bin)] += weight;
        ++n_measurements;
    };

    // Host-side loop, not a GPU ParallelFor (same reasoning as
    // SpectrumKernel.hpp's own shell-binning: this is a histogram
    // accumulation into a small, shared set of bins, not an independent
    // per-cell reduction -- computed once per snapshot, not every
    // substep, so the CPU cost here is not the bottleneck USE_OMP is for).
    const bool has_mask = (mask != nullptr);
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx  = mfi.validbox();
        const auto &state_arr = state.const_array(mfi);
        const auto mask_arr =
            has_mask ? mask->const_array(mfi) : amrex::Array4<int const>{};
        const auto lo = bx.smallEnd();
        const auto hi = bx.bigEnd();

        for (int k = lo[2]; k <= hi[2]; ++k)
        {
            for (int j = lo[1]; j <= hi[1]; ++j)
            {
                for (int i = lo[0]; i <= hi[0]; ++i)
                {
                    if (has_mask && mask_arr(i, j, k) == 0)
                    {
                        // Covered by a finer level -- that level's own
                        // pass owns this physical region, not this one.
                        continue;
                    }
                    const CellPiercing here =
                        get_cell_piercing(state_arr, i, j, k);
                    if (here.count != 2)
                    {
                        if (here.count > 0)
                        {
                            ++n_not_two;
                        }
                        continue;
                    }

                    if (params.method == CurvatureMethod::HessianAnalytic)
                    {
                        // No face/neighbour bookkeeping needed at all --
                        // this method's value depends only on this cell's
                        // own local derivatives (this file's header
                        // comment). Weighted by dx[0] (isotropic grid, one
                        // cell's worth of string length), since there is
                        // no chord length for a single-point measurement.
                        const HessianCurvatureResult hc =
                            hessian_curvature_at_cell(state_arr, dx, i, j,
                                                      k);
                        if (!hc.valid)
                        {
                            ++n_degenerate_hessian;
                            continue;
                        }
                        const double kappa_over_mr =
                            (m_r > 0.0) ? hc.kappa / m_r : 0.0;
                        bin_measurement(kappa_over_mr, dx[0]);
                        continue;
                    }

                    for (int f = 0; f < 2; ++f)
                    {
                        const StringFaceId &face = here.face[f];
                        if (!(face.i == i && face.j == j && face.k == k))
                        {
                            // Not this cell's own face to measure -- its
                            // canonical owner is a neighbour, deferred to
                            // when that neighbour is processed.
                            continue;
                        }
                        const StringFaceId &other_of_this =
                            here.face[1 - f];

                        int ni{}, nj{}, nk{};
                        other_cell_across_face(face, ni, nj, nk);
                        const CellPiercing neighbor =
                            get_cell_piercing(state_arr, ni, nj, nk);
                        if (neighbor.count != 2)
                        {
                            ++n_neighbor;
                            continue;
                        }
                        StringFaceId far_face{};
                        bool found_far = false;
                        if (neighbor.face[0] == face)
                        {
                            far_face  = neighbor.face[1];
                            found_far = true;
                        }
                        else if (neighbor.face[1] == face)
                        {
                            far_face  = neighbor.face[0];
                            found_far = true;
                        }
                        if (!found_far)
                        {
                            ++n_topology;
                            continue;
                        }

                        double kappa{};
                        double weight{};

                        if (params.method == CurvatureMethod::TangentVector)
                        {
                            // Face positions are only needed for the arc
                            // length (ds) below -- the curvature value
                            // itself comes from the tangent field, not
                            // these positions.
                            double x1{}, y1{}, z1{}, x2{}, y2{}, z2{}, x3{},
                                y3{}, z3{};
                            face_center(other_of_this, dx, prob_lo, x1, y1,
                                       z1);
                            face_center(face, dx, prob_lo, x2, y2, z2);
                            face_center(far_face, dx, prob_lo, x3, y3, z3);
                            const double d12 = std::sqrt(
                                (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) +
                                (z2 - z1) * (z2 - z1));
                            const double d23 = std::sqrt(
                                (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) +
                                (z3 - z2) * (z3 - z2));
                            const double arc_length = d12 + d23;

                            // Cell A = (i,j,k), the current cell (owner of
                            // `face`); cell B = (ni,nj,nk), the neighbour
                            // across it -- exactly the two cells the
                            // connectivity graph above already identifies as
                            // straddling this face.
                            const TangentVector t_a = tangent_field_at_cell(
                                state_arr, dx, i, j, k);
                            const TangentVector t_b = tangent_field_at_cell(
                                state_arr, dx, ni, nj, nk);
                            if (!t_a.valid || !t_b.valid || arc_length <= 0.0)
                            {
                                ++n_degenerate_tangent;
                                continue;
                            }

                            // T_a and T_b are each sampled at (approx.) the
                            // midpoint of their own chord -- T_a near the
                            // midpoint of other_of_this->face, T_b near the
                            // midpoint of face->far_face. Their true
                            // separation is therefore half of arc_length
                            // (one chord, not two), not arc_length itself:
                            // dividing by the full two-chord arc length
                            // under-estimated kappa by a factor of 2 (caught
                            // 2026-09-24 via cross-check with a sister
                            // project's port of this method; confirmed
                            // algebraically and numerically against an
                            // exact circle: kappa*R=0.5 with the old
                            // denominator, 1.0 with this one).
                            weight = 0.5 * (d12 + d23);
                            const double dtx = t_b.x - t_a.x;
                            const double dty = t_b.y - t_a.y;
                            const double dtz = t_b.z - t_a.z;
                            kappa = std::sqrt(dtx * dtx + dty * dty +
                                              dtz * dtz) /
                                   weight;
                        }
                        else // InterpolatedPosition
                        {
                            double x1{}, y1{}, z1{}, x2{}, y2{}, z2{}, x3{},
                                y3{}, z3{};
                            const bool ok1 = interpolated_crossing_position(
                                state_arr, dx, prob_lo, other_of_this, x1, y1,
                                z1);
                            const bool ok2 = interpolated_crossing_position(
                                state_arr, dx, prob_lo, face, x2, y2, z2);
                            const bool ok3 = interpolated_crossing_position(
                                state_arr, dx, prob_lo, far_face, x3, y3, z3);
                            if (!ok1 || !ok2 || !ok3)
                            {
                                ++n_no_crossing;
                                continue;
                            }
                            kappa = menger_curvature(x1, y1, z1, x2, y2, z2,
                                                     x3, y3, z3);
                            const double d12 = std::sqrt(
                                (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) +
                                (z2 - z1) * (z2 - z1));
                            const double d23 = std::sqrt(
                                (x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) +
                                (z3 - z2) * (z3 - z2));
                            weight = 0.5 * (d12 + d23);
                        }

                        const double kappa_over_mr =
                            (m_r > 0.0) ? kappa / m_r : 0.0;
                        bin_measurement(kappa_over_mr, weight);
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(plain_local.data(),
                                             params.n_bins);
    amrex::ParallelDescriptor::ReduceRealSum(weighted_local.data(),
                                             params.n_bins);
    amrex::ParallelDescriptor::ReduceLongSum(n_measurements);
    amrex::ParallelDescriptor::ReduceLongSum(n_below);
    amrex::ParallelDescriptor::ReduceLongSum(n_above);
    amrex::ParallelDescriptor::ReduceLongSum(n_not_two);
    amrex::ParallelDescriptor::ReduceLongSum(n_neighbor);
    amrex::ParallelDescriptor::ReduceLongSum(n_topology);
    amrex::ParallelDescriptor::ReduceLongSum(n_degenerate_tangent);
    amrex::ParallelDescriptor::ReduceLongSum(n_no_crossing);
    amrex::ParallelDescriptor::ReduceLongSum(n_degenerate_hessian);

    hist.plain_count      = plain_local;
    hist.length_weighted  = weighted_local;
    hist.n_measurements   = n_measurements;
    hist.n_below_range    = n_below;
    hist.n_above_range    = n_above;
    hist.n_excluded_not_two  = n_not_two;
    hist.n_excluded_neighbor = n_neighbor;
    hist.n_excluded_topology = n_topology;
    hist.n_excluded_degenerate_tangent = n_degenerate_tangent;
    hist.n_excluded_no_crossing        = n_no_crossing;
    hist.n_excluded_degenerate_hessian = n_degenerate_hessian;
    return hist;
}

// Combines one level's histogram into a running composite total: every
// field is a plain sum (bin-for-bin for the two vectors), valid because
// every level in a composite run shares the same CurvatureParams (bin
// range/count, method) -- `acc` and `part` therefore always have
// identically-sized/positioned bins. `acc` starts as a default-
// constructed CurvatureHistogram (empty bins, n_bins taken from the
// first `part` merged in).
inline void
merge_curvature_histogram(CurvatureHistogram &acc,
                          const CurvatureHistogram &part)
{
    if (acc.bin_center_log10.empty())
    {
        acc.bin_center_log10 = part.bin_center_log10;
        acc.plain_count.assign(part.plain_count.size(), 0.0);
        acc.length_weighted.assign(part.length_weighted.size(), 0.0);
    }
    for (std::size_t b = 0; b < part.plain_count.size(); ++b)
    {
        acc.plain_count[b]     += part.plain_count[b];
        acc.length_weighted[b] += part.length_weighted[b];
    }
    acc.n_measurements   += part.n_measurements;
    acc.n_below_range    += part.n_below_range;
    acc.n_above_range    += part.n_above_range;
    acc.n_excluded_not_two  += part.n_excluded_not_two;
    acc.n_excluded_neighbor += part.n_excluded_neighbor;
    acc.n_excluded_topology += part.n_excluded_topology;
    acc.n_excluded_degenerate_tangent += part.n_excluded_degenerate_tangent;
    acc.n_excluded_no_crossing        += part.n_excluded_no_crossing;
    acc.n_excluded_degenerate_hessian += part.n_excluded_degenerate_hessian;
}

#endif // CURVATUREKERNEL_HPP_
