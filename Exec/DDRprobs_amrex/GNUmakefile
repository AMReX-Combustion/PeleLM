TOP = ../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

# Enable flctin namelist
#USE_FLCTS       = TRUE

#USE_PARTICLES   = FALSE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = TRUE

# Enable time-dependent forcing term
#USE_TDEPFORCING = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = DRM19
#Chemistry_Model = CH42STEP
#Chemistry_Model = DODECANE_LU
#Chemistry_Model = WANGDODECANE
#fEXE_sources += wang_dodecane_phi_1000.f
fEXE_sources += drm19Soln_seed_1.00.f
#fEXE_sources += 2stepSoln1.00.f
#fEXE_sources += gri30ion_phi_1.00_Le1.f

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_$(DIM)d.F90 Prob_F.F90
CEXE_headers += 
FEXE_headers += probdata.H Prob_F.H
#f90EXE_sources += pmf_generic.f90

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
