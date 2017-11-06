TOP = ../../..
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
DIM    	          = 2
COMP              = gnu
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = FALSE

Chemistry_Model = DRM19
fEXE_sources += drm19Soln_seed_0.125.f

Blocs   := .

CEXE_sources +=
FEXE_sources += Prob_$(DIM)d.F PROB_F.F
CEXE_headers += 
FEXE_headers += probdata.H MyPROB_F.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM


