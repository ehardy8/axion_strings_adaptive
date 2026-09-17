#ifndef STRINGFINDER_HPP_
#define STRINGFINDER_HPP_

// Plaquette counting of conventions.md sec.8/sec.12: for each cell, the 3
// plaquettes (xy, yz, zx) anchored at its low-index corner (sec.11's tagging
// convention). n_p_plain counts pierced plaquettes; n_p_weighted sums
// |winding| over them -- sec.12: "most strings are singly wound, so the two
// agree closely; a divergence is a diagnostic in itself".
//
// The phase gamma(x) = arg(psi) equals arg(phi) exactly, since psi = R
// phi / v with R > 0 real everywhere (sec.3) -- so this operates directly on
// the stored state, no R-rescaling needed.
//
// `state` must already have at least 1 valid ghost cell (state.FillBoundary)
// before calling count_plaquettes.

#include "PlaquetteWinding.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <cmath>

struct PlaquetteCounts
{
    long n_p_plain{0};
    long n_p_weighted{0};
};

[[nodiscard]] inline PlaquetteCounts
count_plaquettes(const amrex::MultiFab &state)
{
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<long, long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const auto &arrs = state.const_arrays();

    reduce_op.eval(
        state, amrex::IntVect(0), reduce_data,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> ReduceTuple
        {
            const auto &a = arrs[box_no];

            auto phase = [&](int ii, int jj, int kk) -> double
            {
                return std::atan2(static_cast<double>(a(ii, jj, kk, c_psi2)),
                                  static_cast<double>(a(ii, jj, kk, c_psi1)));
            };

            const double th_000 = phase(i, j, k);
            const double th_100 = phase(i + 1, j, k);
            const double th_110 = phase(i + 1, j + 1, k);
            const double th_010 = phase(i, j + 1, k);
            const double th_001 = phase(i, j, k + 1);
            const double th_101 = phase(i + 1, j, k + 1);
            const double th_011 = phase(i, j + 1, k + 1);

            const int w_xy =
                plaquette_winding(th_000, th_100, th_110, th_010);
            const int w_yz =
                plaquette_winding(th_000, th_010, th_011, th_001);
            const int w_zx =
                plaquette_winding(th_000, th_001, th_101, th_100);

            const long plain = (w_xy != 0 ? 1 : 0) + (w_yz != 0 ? 1 : 0) +
                              (w_zx != 0 ? 1 : 0);
            const long weighted =
                std::abs(w_xy) + std::abs(w_yz) + std::abs(w_zx);

            return {plain, weighted};
        });

    const ReduceTuple result = reduce_data.value(reduce_op);
    PlaquetteCounts counts{};
    counts.n_p_plain    = amrex::get<0>(result);
    counts.n_p_weighted = amrex::get<1>(result);
    return counts;
}

#endif // STRINGFINDER_HPP_
