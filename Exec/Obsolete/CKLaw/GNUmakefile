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
COMP              = gnu
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = FALSE

Chemistry_Model= BURKEMOD
fEXE_sources += pmf_BurkeMOD_040_10atm.f

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_$(DIM)d.F90 Prob_F.F90
CEXE_headers += 
FEXE_headers += probdata.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
