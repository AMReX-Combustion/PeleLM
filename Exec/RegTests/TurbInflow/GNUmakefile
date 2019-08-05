AMREX_HOME ?= ../../../../amrex
IAMR_HOME ?= ../../../../IAMR

DEBUG	= FALSE
DIM	= 3
COMP    = gcc
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_OMP   = TRUE

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
F90EXE_sources += $(EBASE)_nd.F90

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
include $(AMREX_HOME)/Src/Base/Make.package
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(IAMR_HOME)/Source/TurbInflow
include $(IAMR_HOME)/Source/TurbInflow/Make.package
vpathdir += $(IAMR_HOME)/Source/TurbInflow

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
