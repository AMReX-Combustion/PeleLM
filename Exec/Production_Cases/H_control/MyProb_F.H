#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define SDIM BL_SPACEDIM

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>
#include <Prob_F.H>
#include <ChemDriver_F.H>
#include <PeleLM_F.H>
#if defined(BL_DO_FLCT)
#include <infl_frc.H>
#endif

#define BL_PROB_PREMIXED_FIXED_INFLOW 'premixed_fixed_inlet'
#define BL_PROB_PREMIXED_CONTROLLED_INFLOW 'premixed_controlled_inlet'
#define BL_PROB_PREMIXED_FREE 'premixed_free'
#define BL_PROB_DIFFUSION 'diffusion_two-stream'
#define BL_PROB_UNDEFINED 'undefined_problem'
#define BL_PROB_FUELBLOB 'fuelblob'
#define BL_PROB_CHAMBER 'chamber'


#define BL_XLO 0
#define BL_YLO 1
#define BL_XHI 2
#define BL_YHI 3
#define BL_ZLO 4
#define BL_ZHI 5
#define BL_INTERIOR 6

#define BL_FUELPIPE 1
#define BL_AMBIENT  2
#define BL_COFLOW   3
#define BL_STICK    4
#define BL_WALL     5
#define BL_VOLUME   6
#define BL_PIPEEND  7

