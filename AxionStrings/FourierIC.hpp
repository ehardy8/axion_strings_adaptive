#ifndef FOURIERIC_HPP_
#define FOURIERIC_HPP_

// Generates a real Gaussian random field with a target mean-square variance,
// Fourier modes occupied for |k| <= k_max and zero above (conventions.md
// sec.7's two initial-condition parameters). Deterministic per
// (seed, mode index): each independent spectral mode's amplitude and phase
// are drawn from a hash of (seed, component, i, j, k), not from an RNG
// stream whose state depends on execution order -- this is what makes
// generation reproducible at fixed (seed, parameters, rank count), per
// CLAUDE.md constraint 7. Checked empirically that it is *not* also
// independent of rank count (1 vs 2 ranks gives a different field, though
// each is internally reproducible) -- amrex::FFT::R2C's internal spectral
// layout evidently depends on the domain decomposition in a way this
// mode-indexed hash does not undo. Rank-count independence was not
// required and is not claimed; only fixed-rank-count reproducibility is.

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>

#include <cmath>
#include <cstdint>

namespace FourierIC
{

// SplitMix64-style hash: cheap, well-mixed, deterministic.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint64_t
hash_u64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    x = x ^ (x >> 31);
    return x;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
uniform_from_hash(std::uint64_t h)
{
    // top 53 bits -> a value in (0, 1)
    return (static_cast<double>(h >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

// Deterministic complex Gaussian amplitude (unit variance per real and
// imaginary part) for spectral mode (comp, i, j, k) under the given seed.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gaussian_mode(std::uint64_t seed, int comp, int i, int j, int k, double &re,
             double &im)
{
    const std::uint64_t key =
        seed ^ (static_cast<std::uint64_t>(comp) << 60) ^
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 40) ^
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)) << 20) ^
        static_cast<std::uint64_t>(static_cast<std::uint32_t>(k));
    const std::uint64_t h1 = hash_u64(key);
    const std::uint64_t h2 = hash_u64(key ^ 0x1234567890ABCDEFULL);

    const double u1 = amrex::max(uniform_from_hash(h1), 1.0e-300);
    const double u2 = uniform_from_hash(h2);

    const double r     = std::sqrt(-2.0 * std::log(u1));
    const double theta = 2.0 * 3.14159265358979323846 * u2;
    re                 = r * std::cos(theta);
    im                 = r * std::sin(theta);
}

// Wraps a 0-based FFT array index into a signed wavenumber: conventions.md
// sec.9's memory layout (positive frequencies counting up in the first
// half, then most-negative to least-negative in the second half).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int wrapped_freq(int idx, int n)
{
    return (idx <= n / 2) ? idx : idx - n;
}

// Fills `field` (real, matching `geom`'s domain) with a Gaussian random
// field occupying modes |k| <= k_max_cells (in units of 1/dx -- the caller
// converts from k_max/m_r), zero above and at k=0 (a target *variance* means
// zero-mean fluctuations). Rescales in real space afterwards to hit
// target_mean_square exactly, sidestepping the need to track the FFT
// library's own normalisation convention through the round trip.
inline void generate(amrex::MultiFab &field, const amrex::Geometry &geom,
                    double k_max_cells, double target_mean_square,
                    std::uint64_t seed, int comp_seed_offset)
{
    const amrex::Box domain = geom.Domain();
    amrex::FFT::R2C<amrex::Real> r2c(domain);

    amrex::MultiFab zero_in(field.boxArray(), field.DistributionMap(), 1, 0);
    zero_in.setVal(0.0);

    const std::uint64_t full_seed =
        seed ^ (static_cast<std::uint64_t>(comp_seed_offset) << 32);
    const int ny           = domain.length(1);
    const int nz           = domain.length(2);
    const double k_max_sq  = k_max_cells * k_max_cells;

    auto fill_mode = [=] AMREX_GPU_DEVICE(int i, int j, int k,
                                         auto &spectral_data)
    {
        const int kx = i; // R2C stores kx >= 0 only
        const int ky = wrapped_freq(j, ny);
        const int kz = wrapped_freq(k, nz);
        const double k_sq =
            static_cast<double>(kx * kx + ky * ky + kz * kz);
        if (k_sq > k_max_sq || k_sq == 0.0)
        {
            spectral_data = amrex::GpuComplex<amrex::Real>(0.0, 0.0);
            return;
        }
        double re{};
        double im{};
        gaussian_mode(full_seed, 0, i, j, k, re, im);
        spectral_data = amrex::GpuComplex<amrex::Real>(
            static_cast<amrex::Real>(re), static_cast<amrex::Real>(im));
    };

    r2c.forwardThenBackward(zero_in, field, fill_mode);

    const double sum_sq =
        static_cast<double>(amrex::MultiFab::Dot(field, 0, field, 0, 1, 0));
    const double n_cells = static_cast<double>(domain.numPts());
    const double current_mean_square = sum_sq / n_cells;
    if (current_mean_square > 0.0)
    {
        const auto rescale = static_cast<amrex::Real>(
            std::sqrt(target_mean_square / current_mean_square));
        field.mult(rescale, 0, 1);
    }
}

} // namespace FourierIC

#endif // FOURIERIC_HPP_
