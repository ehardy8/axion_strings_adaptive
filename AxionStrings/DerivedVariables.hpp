#ifndef DERIVEDVARIABLES_HPP_
#define DERIVEDVARIABLES_HPP_

// AMReX includes
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

// GRTeclyn includes
#include "Coordinates.hpp"

// AxionStrings includes
#include "AxionStringsRHS.hpp"
#include "Wave.hpp"

AMREX_FORCE_INLINE void
calc_analytic_solution(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                       const amrex::MultiFab & /*mf_in*/,
                       const amrex::Geometry &geom, const amrex::Real time,
                       const int * /*bcomp*/, int /*scomp*/);

template <typename model_t>
AMREX_FORCE_INLINE void
calc_energy_density(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                    const amrex::MultiFab &mf_in,
                    const amrex::Geometry /*&geom*/, const amrex::Real /*time*/,
                    const int * /*bcomp*/, int /*scomp*/);

template <typename model_t>
AMREX_FORCE_INLINE void calc_analytic_mf_3d(amrex::MultiFab &mf_out, int dcomp,
                                            const amrex::Geometry &geom,
                                            const amrex::Real time);

// This needs to be a .hpp file because there is a template definition inside
// The template arguments are expanded inline by the compiler in the .hpp file.
#include "DerivedVariables.impl.hpp"

#endif // DERIVEDVARIABLES_HPP_
