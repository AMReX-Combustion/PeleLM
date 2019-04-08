TOP = ../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
DIM    	          = 3
COMP              = gnu
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = FALSE

Chemistry_Model = CHEMH
fEXE_sources += H-flame_0500.f

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_F.F90 Prob_$(DIM)d.F90
CEXE_headers += 
FEXE_headers += probdata.H MyProb_F.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM


