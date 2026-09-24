#ifndef CURVATURECHAINDEBUG_HPP_
#define CURVATURECHAINDEBUG_HPP_

// One-off visual validation tool (2026-09-27, with the user): reconstruct
// a handful of actual, extended string segments from
// PlaquetteCrossing.hpp's interpolated crossing positions -- walking the
// *same* connectivity graph CurvatureKernel.hpp's main histogram loop
// uses, just continued for many steps instead of stopping at 3 points --
// and dump each segment's points plus the interpolated_position method's
// per-point kappa, so they can be plotted and checked by eye: do visibly
// tight bends actually get high kappa, and visibly straight stretches
// actually get kappa near 0? This is not a third independent curvature
// formula (it reuses menger_curvature + interpolated_crossing_position
// verbatim) -- it is a check that the numbers already being computed
// correspond to something geometrically sensible, not just an internally
// self-consistent formula.
//
// Not wired into the production histogram path at all; opt-in
// (axion_strings.curvature.debug_dump_chains) and intended to be run by
// hand for a spot check, not as a routine diagnostic.

#include "CurvatureKernel.hpp"

#include <array>
#include <fstream>
#include <string>
#include <vector>

// Walk outward from `start_cell`'s first pierced face, through successive
// cells that each have exactly two pierced faces (the same criterion
// CurvatureKernel.hpp's main loop uses), for up to `max_steps` hops.
// Stops early at a cell without exactly two pierced faces (a genuine end
// of a locally well-defined single strand), if a face's interpolated
// crossing cannot be found, or if the walk would step outside
// `safe_box` -- get_cell_piercing/interpolated_crossing_position read up
// to 2 cells beyond whichever cell they are called at (StringFinder.hpp's
// 7-corner stencil, applied twice), and unlike the main histogram loop
// (which only ever takes two such hops from a FillBoundary'd MultiFab
// with the documented 3-ghost-cell margin), an unbounded walk here has no
// such guarantee -- unchecked, it reads out of the FAB's allocated
// memory within a few dozen steps on real network data (found the hard
// way, 2026-09-27: a second call at a later snapshot, further from the
// domain's low corner, segfaulted where the first call happened not to).
// Returns the sequence of interpolated positions along the path (empty
// if even the starting face has no valid crossing).
[[nodiscard]] inline std::vector<std::array<double, 3>>
walk_curvature_chain(const amrex::Array4<amrex::Real const> &state,
                     const amrex::GpuArray<amrex::Real, 3> &dx,
                     const amrex::GpuArray<amrex::Real, 3> &prob_lo,
                     const amrex::Box &safe_box, int start_i, int start_j,
                     int start_k, int max_steps)
{
    std::vector<std::array<double, 3>> pts;

    if (!safe_box.contains(amrex::IntVect(start_i, start_j, start_k)))
    {
        return pts;
    }
    const CellPiercing start = get_cell_piercing(state, start_i, start_j,
                                                 start_k);
    if (start.count != 2)
    {
        return pts;
    }

    double x0{}, y0{}, z0{};
    if (!interpolated_crossing_position(state, dx, prob_lo, start.face[0],
                                        x0, y0, z0))
    {
        return pts;
    }
    pts.push_back({x0, y0, z0});

    StringFaceId current_face = start.face[0];
    int ci = start_i, cj = start_j, ck = start_k;

    for (int step = 0; step < max_steps; ++step)
    {
        if (!safe_box.contains(amrex::IntVect(ci, cj, ck)))
        {
            break;
        }
        const CellPiercing piercing = get_cell_piercing(state, ci, cj, ck);
        if (piercing.count != 2)
        {
            break;
        }
        StringFaceId other{};
        if (piercing.face[0] == current_face)
        {
            other = piercing.face[1];
        }
        else if (piercing.face[1] == current_face)
        {
            other = piercing.face[0];
        }
        else
        {
            break; // topology inconsistency -- stop rather than guess
        }

        double xo{}, yo{}, zo{};
        if (!interpolated_crossing_position(state, dx, prob_lo, other, xo,
                                            yo, zo))
        {
            break;
        }
        pts.push_back({xo, yo, zo});

        // Step into whichever of the face's two adjacent cells is *not*
        // the one we are currently in.
        const int cand1i = other.i, cand1j = other.j, cand1k = other.k;
        int cand2i{}, cand2j{}, cand2k{};
        other_cell_across_face(other, cand2i, cand2j, cand2k);
        if (cand1i == ci && cand1j == cj && cand1k == ck)
        {
            ci = cand2i;
            cj = cand2j;
            ck = cand2k;
        }
        else
        {
            ci = cand1i;
            cj = cand1j;
            ck = cand1k;
        }
        current_face = other;
    }

    return pts;
}

// Finds up to `n_chains` well-separated starting cells (every
// `stride`-th pierced-and-well-defined cell encountered in MFIter order,
// so chains are spread across the domain rather than clustered) and
// dumps each one's walked chain to `path`: columns chain_id, step,
// x, y, z, kappa, kappa_over_mr. Single-rank use only (IOProcessor-gated,
// no cross-rank gather) -- a one-off debug tool, not production I/O.
inline void
dump_curvature_chains(const amrex::MultiFab &state,
                      const amrex::Geometry &geom, double m_r, int n_chains,
                      int max_steps, int stride, const std::string &path)
{
    if (!amrex::ParallelDescriptor::IOProcessor())
    {
        return;
    }
    const auto dx      = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    std::ofstream out(path);
    out << "# chain_id step x y z kappa kappa_over_mr\n";

    int n_dumped = 0;
    long n_seen  = 0;
    for (amrex::MFIter mfi(state); mfi.isValid() && n_dumped < n_chains;
        ++mfi)
    {
        const amrex::Box &bx  = mfi.validbox();
        const auto &state_arr = state.const_array(mfi);
        // Shrunk by 3 (this file's own header comment): the margin
        // documented everywhere else in this codebase that calls
        // get_cell_piercing/interpolated_crossing_position, generously
        // covering their up-to-2-cells-ahead stencils.
        const amrex::Box safe_box = amrex::grow(mfi.fabbox(), -3);
        const auto lo              = bx.smallEnd();
        const auto hi              = bx.bigEnd();
        for (int k = lo[2]; k <= hi[2] && n_dumped < n_chains; ++k)
        {
            for (int j = lo[1]; j <= hi[1] && n_dumped < n_chains; ++j)
            {
                for (int i = lo[0]; i <= hi[0] && n_dumped < n_chains; ++i)
                {
                    const CellPiercing here =
                        get_cell_piercing(state_arr, i, j, k);
                    if (here.count != 2)
                    {
                        continue;
                    }
                    ++n_seen;
                    if ((n_seen - 1) % stride != 0)
                    {
                        continue;
                    }

                    const auto pts = walk_curvature_chain(
                        state_arr, dx, prob_lo, safe_box, i, j, k,
                        max_steps);
                    if (pts.size() < 3)
                    {
                        continue;
                    }
                    for (std::size_t p = 1; p + 1 < pts.size(); ++p)
                    {
                        const double kappa = menger_curvature(
                            pts[p - 1][0], pts[p - 1][1], pts[p - 1][2],
                            pts[p][0], pts[p][1], pts[p][2], pts[p + 1][0],
                            pts[p + 1][1], pts[p + 1][2]);
                        const double kappa_over_mr =
                            (m_r > 0.0) ? kappa / m_r : 0.0;
                        out << n_dumped << " " << p << " " << pts[p][0]
                            << " " << pts[p][1] << " " << pts[p][2] << " "
                            << kappa << " " << kappa_over_mr << "\n";
                    }
                    ++n_dumped;
                }
            }
        }
    }
}

#endif // CURVATURECHAINDEBUG_HPP_
