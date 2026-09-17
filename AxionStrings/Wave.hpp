#ifndef WAVE_HPP_
#define WAVE_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// GRTeclyn includes
#include "GRParmParse.hpp"
#include "StateVariables.hpp"

// Placeholder model carried over unmodified from the GRTeclyn KleinGordon
// example (milestone-1 task 1.1). It is replaced by the complex axion-string
// scalar in task 1.2/1.3.
class Wave
{
  public:

    struct params_t
    {
        amrex::Real k_r{1.0};
        amrex::Real scalar_mass{0.0};

        static void check_params()
        {
            GRParmParse wave_pp("wave");
            amrex::Real k_r{1.0};
            wave_pp.queryAdd("wave_vector", k_r);
            amrex::Real scalar_mass{0.0};
            wave_pp.queryAdd("scalar_mass", scalar_mass);
        }

        void fill_params()
        {
            GRParmParse wave_pp("wave");
            wave_pp.get("wave_vector", k_r);
            wave_pp.get("scalar_mass", scalar_mass);
        }
    };

    params_t m_params;
    amrex::Real m_t0{0.0};

    Wave()
    {
        amrex::ParmParse pp;
        pp.query("axion_strings.initial_time", m_t0);

        m_params.fill_params();
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    {
        amrex::Real omega = m_params.k_r;

        amrex::Real rr2 = x * x + y * y + z * z;

        return std::cos(m_params.k_r * rr2 - omega * (t + m_t0));
    };

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    {
        amrex::Real omega = m_params.k_r;

        amrex::Real rr2 = x * x + y * y + z * z;

        return omega * std::sin(m_params.k_r * rr2 - omega * (t + m_t0));
    };

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const amrex::Real &phi) const
    {
        V_of_phi =
            0.5 * m_params.scalar_mass * m_params.scalar_mass * phi * phi;

        dVdphi = m_params.scalar_mass * m_params.scalar_mass * phi;
    }
};

#endif // WAVE_HPP_
