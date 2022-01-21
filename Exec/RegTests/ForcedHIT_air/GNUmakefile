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
DIM             = 3

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

USE_TURBULENT_FORCING = TRUE
USE_VELOCITY          = FALSE
USE_FAST_FORCE        = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = air
# Chemistry_Model = TPRF_98species_gpu

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

ifeq ($(USE_CUDA),TRUE)
CXXFLAGS += -Xptxas --disable-optimizer-constants
endif

ifeq (${USE_VELOCITY}, TRUE)
  #AMRVIS_DIR defines the directory in which we will find pAmrvis (e.g. DataServices, AmrData and FABUTIL)
  AMRVIS_DIR = $(AMREX_HOME)/Src/Extern/amrdata
endif
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
