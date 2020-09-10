TOP = ../../../../
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR_dev
PELELM_HOME ?= ${TOP}/PeleLM_dev

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = TRUE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

#USE_EB = TRUE

TINY_PROFILE = FALSE
#USE_WBAR = TRUE

DEFINES+=-DALWAYS_NEW_J

#######################
USE_SUNDIALS_PP = TRUE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  # use KLU sparse features -- only useful if CVODE is used on CPU
  USE_KLU_PP = FALSE
  # Activates GPU version
  USE_CUDA_SUNDIALS_PP = FALSE
endif
#######################

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
#Transport_dir := EGLib
Transport_dir := Constant

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
