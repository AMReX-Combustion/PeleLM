PELELM_HOME = ../../..
SUBMODS = ${PELELM_HOME}/Submodules
AMREX_HOME         ?= ${SUBMODS}/amrex
IAMR_HOME          ?= ${SUBMODS}/IAMR
PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics
AMREX_HYDRO_HOME   ?= ${SUBMODS}/AMReX-Hydro

#
# Build configuration
#

DIM             = 3
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE

TINY_PROFILE = FALSE

DEFINES+=-DALWAYS_NEW_J

# This allows GetForce to call a Fortran routine
DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE


#######################
USE_SUNDIALS_3x4x = FALSE

#######################

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
