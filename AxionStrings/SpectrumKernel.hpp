#ifndef SPECTRUMKERNEL_HPP_
#define SPECTRUMKERNEL_HPP_

// Power spectrum of a real scalar field, conventions.md sec.9/sec.12
// (milestone-1.md task 1.9). Built on amrex::FFT::R2C, whose forward
// transform's round-trip normalisation (verified against sec.9: "no
// normalisation factor... forward followed by backward multiplies every
// entry by N^3") matches conventions.md's stated convention exactly, so
// its raw output can be used directly as X~_p with no correction factor
// (unlike FourierIC.hpp's generation direction, which sidesteps needing to
// know this by measuring and rescaling in real space instead).
//
// R2C only stores kx >= 0 (Hermitian symmetry for a real field); ky, kz
// range over the full domain with the standard wraparound (FourierIC.hpp's
// wrapped_freq, sec.9's memory layout). Every kx plane except kx=0 and
// kx=Nx/2 (self-conjugate) has an implicit conjugate partner at Nx-kx that
// R2C does not store -- it must be counted twice to recover sums over the
// *full* spectrum that sec.9/sec.12's Parseval-based formulas assume.
//
// Binning follows sec.12 exactly: modeIndex = floor(p_mod + 0.5); shells
// accumulate 4*pi*p_mod^2*|X_p|^2 and are divided by shell count at the
// end (getSpectrum's raw, un-rescaled output); modes with modeIndex >=
// N/2 (outside the inscribed sphere) are dropped from the binned spectrum,
// though both the full-cube and inscribed-sphere *totals* are kept
// (sec.12: "the existing code reports the ratio of full-cube to inscribed-
// sphere energy").
//
// CPU-first (project decision log): binning is done via a host-side loop
// over the (typically modest-sized) spectral MultiFab rather than a GPU
// scatter-reduction; fine for a diagnostic computed once per snapshot, not
// every substep.

#include "FourierIC.hpp" // wrapped_freq

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <vector>

struct Spectrum
{
    std::vector<double> p;             // shell index (dimensionless mode number)
    std::vector<double> shell_average; // getSpectrum's raw 4*pi<p^2|X_p|^2>
    std::vector<long> shell_count;

    double full_cube_energy{0.0};       // Sum over ALL modes of |X_p|^2
    double inscribed_sphere_energy{0.0}; // Sum over modes with p_mod < N/2
};

[[nodiscard]] inline Spectrum
compute_spectrum(const amrex::MultiFab &field, const amrex::Geometry &geom)
{
    const amrex::Box domain = geom.Domain();
    amrex::FFT::R2C<amrex::Real> r2c(domain);

    const auto [spectral_ba, spectral_dm] = r2c.getSpectralDataLayout();
    typename amrex::FFT::R2C<amrex::Real>::cMF spectral(spectral_ba,
                                                        spectral_dm, 1, 0);

    r2c.forward(field, spectral);

    const int nx    = domain.length(0);
    const int ny    = domain.length(1);
    const int nz    = domain.length(2);
    const int n_max = std::max({nx, ny, nz});
    const int n_shells = n_max / 2 + 1;

    std::vector<double> shell_sum(n_shells, 0.0);
    std::vector<long> shell_count(n_shells, 0);
    double full_cube_energy       = 0.0;
    double inscribed_sphere_energy = 0.0;

    // Deliberately NOT #pragma omp parallel (2026-09-22, with the user,
    // implementing USE_OMP): shell_sum/shell_count are indexed by
    // spherical k-shell, not by box, so boxes covering different parts
    // of k-space routinely land in the *same* shell -- and full_cube_
    // energy/inscribed_sphere_energy are plain scalar accumulators.
    // Parallelising this MFIter loop naively would race on all four.
    // Same call as ProjectionKernel.hpp's identical case: this file's
    // own header comment already documents it as "computed once per
    // snapshot, not every substep", so it's not the bottleneck USE_OMP
    // is for, and a correct parallel version needs thread-local partial
    // reductions merged afterward, not just a pragma.
    for (amrex::MFIter mfi(spectral); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto &arr      = spectral.const_array(mfi);
        const auto lo        = bx.smallEnd();
        const auto hi        = bx.bigEnd();

        for (int k = lo[2]; k <= hi[2]; ++k)
        {
            const int kz = FourierIC::wrapped_freq(k, nz);
            for (int j = lo[1]; j <= hi[1]; ++j)
            {
                const int ky = FourierIC::wrapped_freq(j, ny);
                for (int i = lo[0]; i <= hi[0]; ++i)
                {
                    const int kx = i; // R2C stores kx >= 0 only
                    const double p_mod_sq =
                        static_cast<double>(kx * kx + ky * ky + kz * kz);
                    const double p_mod = std::sqrt(p_mod_sq);

                    const auto &c = arr(i, j, k, 0);
                    const double re = static_cast<double>(c.real());
                    const double im = static_cast<double>(c.imag());
                    const double mag_sq = re * re + im * im;

                    // kx=0 and kx=Nx/2 planes are self-conjugate; every
                    // other kx plane's conjugate (at Nx-kx) is not stored
                    // and must be counted here too.
                    const int weight = (kx == 0 || kx == nx / 2) ? 1 : 2;

                    full_cube_energy += weight * mag_sq;
                    if (p_mod < n_max / 2.0)
                    {
                        inscribed_sphere_energy += weight * mag_sq;
                        const auto shell =
                            static_cast<int>(std::floor(p_mod + 0.5));
                        if (shell < n_shells)
                        {
                            shell_sum[shell] +=
                                weight * 4.0 * M_PI * p_mod_sq * mag_sq;
                            shell_count[shell] += weight;
                        }
                    }
                }
            }
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(shell_sum.data(), n_shells);
    amrex::ParallelDescriptor::ReduceLongSum(shell_count.data(), n_shells);
    amrex::ParallelDescriptor::ReduceRealSum(full_cube_energy);
    amrex::ParallelDescriptor::ReduceRealSum(inscribed_sphere_energy);

    Spectrum result{};
    result.full_cube_energy       = full_cube_energy;
    result.inscribed_sphere_energy = inscribed_sphere_energy;
    for (int s = 0; s < n_shells; ++s)
    {
        if (shell_count[s] > 0)
        {
            result.p.push_back(s);
            result.shell_average.push_back(shell_sum[s] /
                                           static_cast<double>(shell_count[s]));
            result.shell_count.push_back(shell_count[s]);
        }
    }
    return result;
}

#endif // SPECTRUMKERNEL_HPP_
