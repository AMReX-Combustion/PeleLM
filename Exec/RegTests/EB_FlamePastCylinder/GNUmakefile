PELELM_HOME = ../../..
SUBMODS = ${PELELM_HOME}/Submodules
AMREX_HOME         ?= ${SUBMODS}/amrex
IAMR_HOME          ?= ${SUBMODS}/IAMR
PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics
AMREX_HYDRO_HOME   ?= ${SUBMODS}/AMReX-Hydro

#
# Build configuration
#

# AMREX options
DIM             = 2
USE_EB          = TRUE

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# MISC options
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = LiDryer

# This sets the EOS to one in $(PELE_PHYSICS_HOME)/Eos
Eos_Model     := Fuego

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
