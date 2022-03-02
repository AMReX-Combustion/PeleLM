TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELM_HOME        ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro


# AMReX
DIM             = 3
DEBUG           = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLM

# PelePhysics
Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
