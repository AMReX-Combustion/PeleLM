TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
IAMR_HOME          ?= ${TOP}/IAMR
PELELM_HOME        ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics

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

# CVODE
USE_SUNDIALS_PP = FALSE
USE_KLU_PP      = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Null

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
