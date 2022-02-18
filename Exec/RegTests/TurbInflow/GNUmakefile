PELELM_HOME = ../../..
SUBMODS = ${PELELM_HOME}/Submodules
AMREX_HOME         ?= ${SUBMODS}/amrex
IAMR_HOME          ?= ${SUBMODS}/IAMR
PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics
AMREX_HYDRO_HOME   ?= ${SUBMODS}/AMReX-Hydro

DEBUG	= FALSE
#DEBUG	= TRUE
DIM	= 3
COMP    = gcc
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_OMP   = FALSE

EBASE     = main
#EBASE     = turbToFab

#NEEDS_FORT = FALSE
NEEDS_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_FORT), TRUE)
  F90EXE_sources += $(EBASE)_nd.F90
endif

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
include $(AMREX_HOME)/Src/Base/Make.package
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
include $(AMREX_HOME)/Tools/C_util/Make.package
vpathdir += $(AMREX_HOME)/Tools/C_util

INCLUDE_LOCATIONS += $(PELE_PHYSICS_HOME)/Utility/TurbInflow/
include $(PELE_PHYSICS_HOME)/Utility/TurbInflow/Make.package
vpathdir += $(PELE_PHYSICS_HOME)/Utility/TurbInflow

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
