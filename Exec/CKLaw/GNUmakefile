TOP = ../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = TRUE
DIM    	          = 2
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = TRUE

# Enable flctin namelist
USE_FLCTS         = FALSE

CHEMISTRY_MODEL= BURKEMOD
fEXE_sources += pmf_BurkeMOD_040_10atm.f

Blocs   := .

CEXE_sources +=
FEXE_sources += Prob_$(DIM)d.F
CEXE_headers += 
FEXE_headers += probdata.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
