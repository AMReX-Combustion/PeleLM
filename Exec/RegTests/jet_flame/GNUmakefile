PELELM_HOME = ../../..
SUBMODS ?= ${PELELM_HOME}/Submodules
AMREX_HOME         ?= ${SUBMODS}/amrex
IAMR_HOME          ?= ${SUBMODS}/IAMR
PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics
AMREX_HYDRO_HOME   ?= ${SUBMODS}/AMReX-Hydro

#
# Build configuration
#

# AMREX options
DIM             = 3

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# MISC options
DEBUG           = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# PeleLM options
DO_2S_CONVERGENCE=FALSE

USE_TURBINFLOW  = TRUE

ifeq ($(DO_2S_CONVERGENCE),TRUE)
	Chemistry_Model = chem-CH4-2step
	F90EXE_sources += 2step_Phi0.8_refined_CubicInterp.F90
	DEFINES+=-DINTERP_PMF_AS_POINT
else
	Chemistry_Model = drm19
endif

# This sets the EOS model to one in $(PELE_PHYSICS_HOME)/Eos
Eos_Model   := Fuego

# This sets the transport to one in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Simple

Blocs   := .
CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
