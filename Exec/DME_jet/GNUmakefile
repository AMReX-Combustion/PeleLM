TOP = ../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# User-settable variables.
#
PRECISION       = DOUBLE
DEBUG           = FALSE
DIM             = 2
COMP            = gnu
USE_OMP         = FALSE
USE_MPI         = TRUE
VERBOSE         = FALSE


# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable time-dependent forcing term
#USE_TDEPFORCING = FALSE

# Choose model (from list below), and pmf file
CHEMISTRY_MODEL = LUDME
fEXE_sources += LuDME_0700.f

Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F PROB_F.F
CEXE_headers += 
FEXE_headers += probdata.H PROB_F.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

#DEFINES += -DDO_JBB_HACK
DEFINES += -DTRIGGER_NEW_J
#DEFINES += -DALWAYS_NEW_J
#DEFINES += -DSOLN_IS_1D

#DEFINES += -DDO_JBB_HACK_TEMP
#DEFINES += -DHACK_TEMP_MIN=273.d0
#DEFINES += -DHACK_TEMP_MAX=1900.d0
