#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>
#include <mechanism.H>
#include <NavierStokesBase.H>

// Xvel, Yvel, Zvel and Density are already defined in NavierStokesBase.H

#define DEF_first_spec Density + 1
#define DEF_RhoH DEF_first_spec + NUM_SPECIES
#define DEF_Temp DEF_RhoH + 1
#define DEF_RhoRT DEF_Temp + 1
#define DEF_first_passive DEF_RhoRT + 1
#ifdef PELELM_USE_SOOT
#define DEF_NUM_SOOT_VARS NUM_SOOT_MOMENTS + 1
#define DEF_first_soot DEF_first_passive
#else
#define DEF_NUM_SOOT_VARS 0
#endif
#define DEF_NUM_PASSIVE DEF_NUM_SOOT_VARS
#define DEF_NUM_STATE AMREX_SPACEDIM + NUM_SPECIES + 4 + DEF_NUM_PASSIVE
#define DEF_NUM_SCALARS  DEF_NUM_STATE - Density

#endif
