GRTECLYN_HOME := $(realpath ../GRTeclyn)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

EBASE = AxionStrings

AMREX_HOME ?= $(realpath ../amrex)

# Fourier-mode initial conditions (task 1.5) use amrex::FFT::R2C.
USE_FFT := TRUE
FFTW_HOME ?= /opt/homebrew/opt/fftw

# Grid precision (2026-09-19, with the user): default to float ahead of
# moving to bigger, cluster-scale runs -- amrex::Real on the evolution
# path is a build flag, and this is a bandwidth-bound explicit stencil
# code, so float roughly halves both the state's memory footprint and
# FillBoundary/FillPatch communication volume (the "stepping forward and
# regridding" cost the user asked about). This is not a novel choice:
# conventions.md sec.14's own decision log already flags "grid precision:
# float to be tested, double the fallback" and specifies the risk
# precisely (cancellation in |psi|^2 - R^2 erodes float's headroom as
# log(m_r/H) grows, costing about Delta_log ~ 0.46 of reach relative to
# double -- q/the spectrum are the least-exposed observables, radial
# energy is the early-warning diagnostic; conventions.md sec.13's T5 test
# is the formal check of this, not yet run). Reductions/energy sums/
# spectral accumulation stay in double regardless (CLAUDE.md constraint
# 3) -- this flag only affects per-cell field storage and stencil work.
# Override with `make PRECISION=DOUBLE ...` to go back (command-line
# make variables always win over this default) -- the executable name
# picks up a .FLOAT/.DOUBLE tag either way, so both builds can coexist.
PRECISION ?= FLOAT

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source

# Name all directories for relevant source and header files
# Don't forget the dot! (i.e., example specific files are in current directory)
src_dirs := $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Tagging \
            $(realpath .)

# Find all source and header files in the src_dirs
CEXE_sources += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.impl.hpp)))

VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules
