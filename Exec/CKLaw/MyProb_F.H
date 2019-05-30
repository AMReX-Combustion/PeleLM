#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#define SDIM BL_SPACEDIM

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <Prob_F.H>
#include <AMReX_ArrayLim.H>
#include <ChemDriver_F.H>
#include <PeleLM_F.H>
#if defined(BL_DO_FLCT)
#include <infl_frc.H>
#endif

#define XLO 0
#define YLO 1
#define XHI 2
#define YHI 3
#define ZLO 4
#define ZHI 5

#define M_FUEL  1
#define M_CO    2
#define M_WALL  3

      
#define BL_FUELPIPE 1
#define BL_COFLOW   2
#define BL_STICK    3
#define BL_WALL     4
#define BL_AMBIENT  5
#define BL_VOLUME   6
#define BL_PIPEEND  7
