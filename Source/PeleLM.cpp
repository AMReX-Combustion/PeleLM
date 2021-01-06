//
// "Divu_Type" means S, where divergence U = S
// "RhoYchemProd_Type" means -omega_l/rho, i.e., the mass rate of decrease of species l due
//             to kinetics divided by rho
//
//
#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <vector>
#include <sys/stat.h>

#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ErrorList.H>
#include <PeleLM.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_SPACE.H>
#include <AMReX_Interpolater.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Utility.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <NS_util.H>
#include <AMReX_GpuContainers.H>

#include <PPHYS_CONSTANTS.H>
#include <PeleLM_K.H>
#include <pelelm_prob.H>
#include <pelelm_prob_parm.H>
#include <PeleLM_parm.H>
#include <pmf_data.H>

#if defined(AMREX_USE_NEWMECH) || defined(AMREX_USE_VELOCITY)
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#include <reactor.h>

#include <Prob_F.H>
#include <NAVIERSTOKES_F.H>
#include <DERIVE_F.H>

#include <iamr_godunov.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EBAmrUtil.H>
#include <iamr_mol.H>
#endif

#include <AMReX_buildInfo.H>
//fixme, for writesingle level plotfile
#include<AMReX_PlotFileUtil.H>

using namespace amrex;

static Box stripBox; // used for debugging

#ifdef AMREX_USE_FLOAT
const Real Real_MIN = FLT_MIN;
const Real Real_MAX = FLT_MAX;
#else
const Real Real_MIN = DBL_MIN;
const Real Real_MAX = DBL_MAX;
#endif

const int  GEOM_GROW   = 1;
const int  PRESS_GROW  = 1;
const int  DIVU_GROW   = 1;
const Real bogus_value =  1.e20;
const int  DQRAD_GROW  = 1;
const int  YDOT_GROW   = 1;
const int  HYPF_GROW   = 1;
const int  LinOp_grow  = 1;

static Real              typical_RhoH_value_default = -1.e10;
static const std::string typical_values_filename("typical_values.fab");

namespace
{
  bool initialized = false;
}
//
// Set all default values in Initialize()!!!
//
namespace
{
  std::set<std::string> ShowMF_Sets;
  std::string           ShowMF_Dir;
  bool                  ShowMF_Verbose;
  bool                  ShowMF_Check_Nans;
  FABio::Format         ShowMF_Fab_Format;
  bool                  do_not_use_funccount;
  Real                  crse_dt;
  int                   chem_box_chop_threshold;
  int                   num_deltaT_iters_MAX;
  Real                  deltaT_norm_max;
  int                   num_forkjoin_tasks;
  bool                  forkjoin_verbose;
}

Real PeleLM::p_amb_old;
Real PeleLM::p_amb_new;
Real PeleLM::dp0dt;
Real PeleLM::thetabar;
Real PeleLM::dpdt_factor;
int  PeleLM::closed_chamber;
int  PeleLM::num_divu_iters;
int  PeleLM::init_once_done;
int  PeleLM::do_OT_radiation;
int  PeleLM::do_heat_sink;
int  PeleLM::RhoH;
int  PeleLM::do_diffuse_sync;
int  PeleLM::do_reflux_visc;
int  PeleLM::RhoYdot_Type;
int  PeleLM::FuncCount_Type;
int  PeleLM::divu_ceiling;
Real PeleLM::divu_dt_factor;
Real PeleLM::min_rho_divu_ceiling;
int  PeleLM::have_rhort;
int  PeleLM::RhoRT;
int  PeleLM::first_spec;
int  PeleLM::last_spec;
Vector<std::string>  PeleLM::spec_names;
int  PeleLM::floor_species;
int  PeleLM::do_set_rho_to_species_sum;
Real PeleLM::prandtl;
Real PeleLM::schmidt;
Real PeleLM::constant_thick_val;
Array<Real, 4> PeleLM::Beta_mix;
Array<Real, NUM_SPECIES> PeleLM::spec_Bilger_fact;
Real PeleLM::Zfu;
Real PeleLM::Zox;
bool PeleLM::mixture_fraction_ready;
int  PeleLM::unity_Le;
int  PeleLM::use_wbar;
Real PeleLM::htt_tempmin;
Real PeleLM::htt_tempmax;
Real PeleLM::htt_hmixTYP;
int  PeleLM::zeroBndryVisc;
int  PeleLM::do_check_divudt;
int  PeleLM::hack_nochem;
int  PeleLM::hack_nospecdiff;
int  PeleLM::hack_noavgdivu;
int  PeleLM::use_tranlib;
Real PeleLM::trac_diff_coef;
bool PeleLM::plot_reactions;
bool PeleLM::plot_consumption;
bool PeleLM::plot_heat_release;
int  PeleLM::ncells_chem;
bool PeleLM::use_typ_vals_chem = 0;
Real PeleLM::relative_tol_chem = 1.0e-10;
Real PeleLM::absolute_tol_chem = 1.0e-10;
static bool plot_rhoydot;
int  PeleLM::nGrowAdvForcing=1;
bool PeleLM::avg_down_chem;
int  PeleLM::reset_typical_vals_int=-1;
Real PeleLM::typical_Y_val_min=1.e-10;
std::map<std::string,Real> PeleLM::typical_values_FileVals;

std::unique_ptr<ProbParm> PeleLM::prob_parm;
std::unique_ptr<ACParm> PeleLM::ac_parm;

std::string                                PeleLM::turbFile;
std::map<std::string, Vector<std::string> > PeleLM::auxDiag_names;

// Active control defaults
bool         PeleLM::ctrl_active = 0;
bool         PeleLM::ctrl_use_temp = 0;
Real         PeleLM::ctrl_tauControl = -1.0;
Real         PeleLM::ctrl_cfix = 0.0;
Real         PeleLM::ctrl_coftOld = -1.0;
Real         PeleLM::ctrl_sest = 0.0;
Real         PeleLM::ctrl_corr = 1.0;
Real         PeleLM::ctrl_V_in = 0.0;
Real         PeleLM::ctrl_V_in_old = 0.0;
Real         PeleLM::ctrl_changeMax = 1.0;
Real         PeleLM::ctrl_tBase = 0.0;
Real         PeleLM::ctrl_dV = 0.0;
Real         PeleLM::ctrl_scale = 1.0;
Real         PeleLM::ctrl_zBase = 0.0;
Real         PeleLM::ctrl_h = -1.0;
Real         PeleLM::ctrl_velMax = 10000.0;
Real         PeleLM::ctrl_temperature = -1.0;
int          PeleLM::ctrl_verbose = 0;
int          PeleLM::ctrl_NavgPts = 3;
int          PeleLM::ctrl_flameDir = AMREX_SPACEDIM-1;
int          PeleLM::ctrl_pseudoGravity = 0;
int          PeleLM::ctrl_nfilled = -1;
std::string  PeleLM::ctrl_AChistory = "AC_History.dat";
Vector<Real> PeleLM::ctrl_time_pts;
Vector<Real> PeleLM::ctrl_velo_pts;
Vector<Real> PeleLM::ctrl_cntl_pts;

Vector<Real> PeleLM::typical_values;

int PeleLM::sdc_iterMAX;
int PeleLM::num_mac_sync_iter;
int PeleLM::deltaT_verbose = 0;
int PeleLM::deltaT_crashOnConvFail = 1;

int PeleLM::mHtoTiterMAX;
Vector<amrex::Real> PeleLM::mTmpData;

Vector<AMRErrorTag> PeleLM::errtags;

static
std::string
to_upper (const std::string& s)
{
  std::string rtn = s;
  for (unsigned int i=0; i<rtn.length(); i++)
  {
    rtn[i] = toupper(rtn[i]);
  }
  return rtn;
}

static std::map<std::string,FABio::Format> ShowMF_Fab_Format_map;


#ifdef AMREX_PARTICLES
namespace
{
  bool do_curvature_sample = false;
}

void
PeleLM::read_particle_params ()
{
  ParmParse ppht("ht");
  ppht.query("do_curvature_sample", do_curvature_sample);
}

int
PeleLM::timestamp_num_extras ()
{
  return do_curvature_sample ? 1 : 0;
}

void
PeleLM::timestamp_add_extras (int lev,
				    Real time,
				    MultiFab& mf)
{
  if (do_curvature_sample)
  {
    AmrLevel& amr_level = parent->getLevel(lev);
    int cComp = mf.nComp()-1;

    amr_level.derive("mean_progress_curvature", time, mf, cComp);

    mf.setBndry(0,cComp,1);
  }
}

#endif /*AMREX_PARTICLES*/

void
PeleLM::compute_rhohmix (Real      time,
                         MultiFab& rhohmix,
                         int       dComp)
{
   const Real strt_time = ParallelDescriptor::second();
   const TimeLevel whichTime = which_time(State_Type,time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
   const MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type): get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(rhohmix, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = S.array(mfi,Density);
         auto const& rhoY    = S.array(mfi,first_spec);
         auto const& T       = S.array(mfi,Temp);
         auto const& rhoHm   = rhohmix.array(mfi,dComp);

         amrex::ParallelFor(bx, [rho, rhoY, T, rhoHm]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
         });
      }
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::compute_rhohmix(): lev: " << level << ", time: " << run_time << '\n';
   }
}

int 
PeleLM::getSpeciesIdx(const std::string& spName)
{
    for (int i=0; i<NUM_SPECIES; i++) {
        if (spName == spec_names[i]) {
            return i;
        }
    }
    return -1;
}

void
PeleLM::Initialize ()
{
  if (initialized) return;

  PeleLM::Initialize_specific();
  
  NavierStokesBase::Initialize();

  //
  // Set all default values here!!!
  //
  ShowMF_Fab_Format_map["ASCII"] = FABio::FAB_ASCII;
  ShowMF_Fab_Format_map["IEEE"] = FABio::FAB_IEEE;
  ShowMF_Fab_Format_map["NATIVE"] = FABio::FAB_NATIVE;
  ShowMF_Fab_Format_map["8BIT"] = FABio::FAB_8BIT;
  ShowMF_Fab_Format_map["IEEE_32"] = FABio::FAB_IEEE_32;
  ShowMF_Verbose          = true;
  ShowMF_Check_Nans       = true;
  ShowMF_Fab_Format       = ShowMF_Fab_Format_map["ASCII"];
  do_not_use_funccount    = false;
  crse_dt                 = -1;
  chem_box_chop_threshold = -1;

  PeleLM::p_amb_old                 = -1.0;
  PeleLM::p_amb_new                 = -1.0;
  PeleLM::num_divu_iters            = 1;
  PeleLM::init_once_done            = 0;
  PeleLM::do_OT_radiation           = 0;
  PeleLM::do_heat_sink              = 0;
  PeleLM::RhoH                      = -1;
  PeleLM::do_diffuse_sync           = 1;
  PeleLM::do_reflux_visc            = 1;
  PeleLM::RhoYdot_Type                 = -1;
  PeleLM::FuncCount_Type            = -1;
  PeleLM::divu_ceiling              = 0;
  PeleLM::divu_dt_factor            = .5;
  PeleLM::min_rho_divu_ceiling      = -1.e20;
  PeleLM::have_rhort                = 0;
  PeleLM::RhoRT                     = -1;
  PeleLM::first_spec                = -1;
  PeleLM::last_spec                 = -2;
  PeleLM::floor_species             = 1;
  PeleLM::do_set_rho_to_species_sum = 1;
  PeleLM::prandtl                   = .7;
  PeleLM::schmidt                   = .7;
  PeleLM::constant_thick_val        = -1;

  PeleLM::Beta_mix = {0};
  PeleLM::spec_Bilger_fact = {0};
  PeleLM::Zfu = -1;
  PeleLM::Zox = -1;
  PeleLM::mixture_fraction_ready    = false;
  PeleLM::unity_Le                  = 0;
  PeleLM::use_wbar                  = 1;
  PeleLM::htt_tempmin               = 298.0;
  PeleLM::htt_tempmax               = 40000.;
  PeleLM::htt_hmixTYP               = -1.;
  PeleLM::zeroBndryVisc             = 0;
  PeleLM::do_check_divudt           = 1;
  PeleLM::hack_nochem               = 0;
  PeleLM::hack_nospecdiff           = 0;
  PeleLM::hack_noavgdivu            = 0;
  PeleLM::use_tranlib               = 0;
  PeleLM::trac_diff_coef            = 0.0;
  PeleLM::turbFile                  = "";
  PeleLM::plot_reactions            = false;
  PeleLM::plot_consumption          = true;
  PeleLM::plot_heat_release         = true;
  plot_rhoydot                            = false;
  PeleLM::avg_down_chem             = false;
  PeleLM::reset_typical_vals_int    = -1;
  PeleLM::typical_values_FileVals.clear();

  PeleLM::sdc_iterMAX               = 1;
  PeleLM::num_mac_sync_iter         = 1;
  PeleLM::mHtoTiterMAX              = 20;
  PeleLM::ncells_chem               = 1;

  ParmParse pp("ns");
  pp.query("do_diffuse_sync",do_diffuse_sync);
  AMREX_ASSERT(do_diffuse_sync == 0 || do_diffuse_sync == 1);
  pp.query("do_reflux_visc",do_reflux_visc);
  AMREX_ASSERT(do_reflux_visc == 0 || do_reflux_visc == 1);

  verbose = pp.contains("v");

  pp.query("divu_ceiling",divu_ceiling);
  AMREX_ASSERT(divu_ceiling >= 0 && divu_ceiling <= 3);
  pp.query("divu_dt_factor",divu_dt_factor);
  AMREX_ASSERT(divu_dt_factor>0 && divu_dt_factor <= 1.0);
  pp.query("min_rho_divu_ceiling",min_rho_divu_ceiling);
  if (divu_ceiling) AMREX_ASSERT(min_rho_divu_ceiling >= 0.0);

  pp.query("htt_tempmin",htt_tempmin);
  pp.query("htt_tempmax",htt_tempmax);

  pp.query("floor_species",floor_species);
  AMREX_ASSERT(floor_species == 0 || floor_species == 1);

  pp.query("do_set_rho_to_species_sum",do_set_rho_to_species_sum);

  pp.query("num_divu_iters",num_divu_iters);

  pp.query("do_not_use_funccount",do_not_use_funccount);

  pp.query("schmidt",schmidt);
  pp.query("prandtl",prandtl);
  pp.query("unity_Le",unity_Le);
  unity_Le = unity_Le ? 1 : 0;
  if (unity_Le)
  {
    schmidt = prandtl;
    if (verbose) amrex::Print() << "PeleLM::read_params: Le=1, setting Sc = Pr" << '\n';
  }

  pp.query("use_wbar",use_wbar);
  pp.query("sdc_iterMAX",sdc_iterMAX);
  pp.query("num_mac_sync_iter",num_mac_sync_iter);

  pp.query("thickening_factor",constant_thick_val);
  if (constant_thick_val != -1)
  {
    if (verbose)
      amrex::Print() << "PeleLM::read_params: using a constant thickening factor = " 
		     << constant_thick_val << '\n';
  }

  pp.query("hack_nochem",hack_nochem);
  pp.query("hack_nospecdiff",hack_nospecdiff);
  pp.query("hack_noavgdivu",hack_noavgdivu);
  pp.query("do_check_divudt",do_check_divudt);
  pp.query("avg_down_chem",avg_down_chem);
  pp.query("reset_typical_vals_int",reset_typical_vals_int);

  pp.query("do_OT_radiation",do_OT_radiation);
  do_OT_radiation = (do_OT_radiation ? 1 : 0);
  pp.query("do_heat_sink",do_heat_sink);
  do_heat_sink = (do_heat_sink ? 1 : 0);

  pp.query("use_tranlib",use_tranlib);
  if (use_tranlib == 1) {
    if (verbose) amrex::Print() << "PeleLM::read_params: Using Tranlib transport " << '\n';
  }
  else {
    if (verbose) amrex::Print() << "PeleLM::read_params: Using EGLib transport " << '\n';
  }

  pp.query("turbFile",turbFile);

  pp.query("zeroBndryVisc",zeroBndryVisc);
  //
  // Read in scalar value and use it as tracer.
  //
  pp.query("scal_diff_coefs",trac_diff_coef);

  for (int i = 0; i < visc_coef.size(); i++)
    visc_coef[i] = bogus_value;

  // Get some useful amr inputs
  ParmParse ppa("amr");

  // Useful for debugging
  ParmParse pproot;
  if (int nsv=pproot.countval("ShowMF_Sets"))
  {
    Vector<std::string> ShowMF_set_names(nsv);
    pproot.getarr("ShowMF_Sets",ShowMF_set_names);
    for (int i=0; i<nsv; ++i) {
      ShowMF_Sets.insert(ShowMF_set_names[i]);
    }
    ShowMF_Dir="."; pproot.query("ShowMF_Dir",ShowMF_Dir);
    pproot.query("ShowMF_Verbose",ShowMF_Verbose);
    pproot.query("ShowMF_Check_Nans",ShowMF_Check_Nans);
    std::string format="NATIVE"; pproot.query("ShowMF_Fab_Format",format);
    if (ShowMF_Fab_Format_map.count(to_upper(format)) == 0) {
      amrex::Abort("Unknown FABio format label");
    }
    ShowMF_Fab_Format = ShowMF_Fab_Format_map[format];

    if (ShowMF_Verbose>0 && ShowMF_set_names.size()>0) {
      amrex::Print() << "   ******************************  Debug: ShowMF_Sets: ";
      for (int i=0; i<ShowMF_set_names.size(); ++i) {
	amrex::Print() << ShowMF_set_names[i] << " ";
      }
      amrex::Print() << '\n';
    }

  }

#ifdef AMREX_PARTICLES
  read_particle_params();
#endif

  if (verbose)
  {
    amrex::Print() << "\nDumping ParmParse table:\n \n";

    if (ParallelDescriptor::IOProcessor()) {
      ParmParse::dumpTable(std::cout);
    }

    amrex::Print() << "\n... done dumping ParmParse table.\n" << '\n';
  }

  amrex::ExecOnFinalize(PeleLM::Finalize);

  initialized = true;
}

void
PeleLM::Initialize_specific ()
{

    num_deltaT_iters_MAX    = 10;
    deltaT_norm_max         = 1.e-10;
    num_forkjoin_tasks      = 1;
    forkjoin_verbose        = false;
  
    ParmParse pplm("peleLM");
    
    pplm.query("num_forkjoin_tasks",num_forkjoin_tasks);
    pplm.query("forkjoin_verbose",forkjoin_verbose);
    pplm.query("num_deltaT_iters_MAX",num_deltaT_iters_MAX);
    pplm.query("deltaT_norm_max",deltaT_norm_max);
    pplm.query("deltaT_verbose",deltaT_verbose);
    pplm.query("deltaT_convFailCrash",deltaT_crashOnConvFail);

    pplm.query("use_typ_vals_chem",use_typ_vals_chem);
    pplm.query("relative_tol_chem",relative_tol_chem);
    pplm.query("absolute_tol_chem",absolute_tol_chem);
    
    // Get boundary conditions
    Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
    Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
    pplm.getarr("lo_bc",lo_bc_char,0,AMREX_SPACEDIM);
    pplm.getarr("hi_bc",hi_bc_char,0,AMREX_SPACEDIM);


    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    bool flag_closed_chamber = false;
    for (int dir = 0; dir<AMREX_SPACEDIM; dir++){
      if (!lo_bc_char[dir].compare("Interior")){
        lo_bc[dir] = 0;
      } else if (!lo_bc_char[dir].compare("Inflow")){
        lo_bc[dir] = 1;
        flag_closed_chamber = true;
      } else if (!lo_bc_char[dir].compare("Outflow")){
        lo_bc[dir] = 2;
        flag_closed_chamber = true;
      } else if (!lo_bc_char[dir].compare("Symmetry")){
        lo_bc[dir] = 3;
      } else if (!lo_bc_char[dir].compare("SlipWallAdiab")){
        lo_bc[dir] = 4;
      } else if (!lo_bc_char[dir].compare("NoSlipWallAdiab")){
        lo_bc[dir] = 5;
      } else if (!lo_bc_char[dir].compare("SlipWallIsotherm")){
        lo_bc[dir] = 6;
      } else if (!lo_bc_char[dir].compare("NoSlipWallIsotherm")){
        lo_bc[dir] = 7;
      } else {
        amrex::Abort("Wrong boundary condition word in lo_bc, please use: Interior, Inflow, Outflow, "
                     "Symmetry, SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm");
      }

      if (!hi_bc_char[dir].compare("Interior")){
        hi_bc[dir] = 0;
      } else if (!hi_bc_char[dir].compare("Inflow")){
        hi_bc[dir] = 1;
        flag_closed_chamber = true;
      } else if (!hi_bc_char[dir].compare("Outflow")){
        hi_bc[dir] = 2;
        flag_closed_chamber = true;
      } else if (!hi_bc_char[dir].compare("Symmetry")){
        hi_bc[dir] = 3;
      } else if (!hi_bc_char[dir].compare("SlipWallAdiab")){
        hi_bc[dir] = 4;
      } else if (!hi_bc_char[dir].compare("NoSlipWallAdiab")){
        hi_bc[dir] = 5;
      } else if (!hi_bc_char[dir].compare("SlipWallIsotherm")){
        hi_bc[dir] = 6;
      } else if (!hi_bc_char[dir].compare("NoSlipWallIsotherm")){
        hi_bc[dir] = 7;
      } else {
        amrex::Abort("Wrong boundary condition word in hi_bc, please use: Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
      }
    }

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }
  
    read_geometry();
    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (DefaultGeometry().isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (DefaultGeometry().isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "PeleLM::variableSetUp:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Abort("PeleLM::Initialize()");
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "PeleLM::variableSetUp:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Abort("PeleLM::Initialize()");
                }
            } 
        }
    }

    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (!DefaultGeometry().isPeriodic(dir))
            {
              if (lo_bc[dir] == Interior)
              {
                  std::cerr << "PeleLM::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("PeleLM::Initialize()");
              }
              if (hi_bc[dir] == Interior)
              {
                  std::cerr << "PeleLM::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("PeleLM::Initialize()");
              }
            }
        }
    } 

   PeleLM::closed_chamber  = 1;
   if (flag_closed_chamber){
      PeleLM::closed_chamber  = 0;
   }

   PeleLM::dpdt_factor = 1.0;
   pplm.query("dpdt_factor",dpdt_factor);

}

void
PeleLM::Finalize ()
{
  initialized = false;
}

static
Box
getStrip(const Geometry& geom)
{
  const Box& box = geom.Domain();
  IntVect be = box.bigEnd();
  IntVect se = box.smallEnd();
  se[0] = (int) 0.5*(se[0]+be[0]);
  be[0] = se[0];
  return Box(se,be);
}

void
showMFsub(const std::string&   mySet,
          const MultiFab&      mf,
          const Box&           box,
          const std::string&   name,
          int                  lev = -1,
          int                  iter = -1) // Default value = no append 2nd integer
{
  if (ShowMF_Sets.count(mySet)>0)
  {
    const FABio::Format saved_format = FArrayBox::getFormat();
    FArrayBox::setFormat(ShowMF_Fab_Format);
    std::string DebugDir(ShowMF_Dir);
    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(DebugDir, 0755))
        amrex::CreateDirectoryFailed(DebugDir);
    ParallelDescriptor::Barrier();

    std::string junkname = name;
    if (lev>=0) {
      junkname = amrex::Concatenate(junkname+"_",lev,1);
    }
    if (iter>=0) {
      junkname = amrex::Concatenate(junkname+"_",iter,1);
    }
    junkname = DebugDir + "/" + junkname;

    if (ShowMF_Verbose>0) {
      amrex::Print() << "   ******************************  Debug: writing " 
		     << junkname << '\n';
    }

    FArrayBox sub(box,mf.nComp());

    mf.copyTo(sub,0,0,mf.nComp(),0);

    if (ShowMF_Check_Nans)
    {
      AMREX_ASSERT(!sub.contains_nan<RunOn::Host>(box,0,mf.nComp()));
    }
    std::ofstream os;
    os.precision(15);
    os.open(junkname.c_str());
    sub.writeOn(os);
    os.close();
    FArrayBox::setFormat(saved_format);
  }
}

void
showMF(const std::string&   mySet,
       const MultiFab&      mf,
       const std::string&   name,
       int                  lev = -1,
       int                  iter = -1, // Default value = no append 2nd integer
       int                  step = -1) // Default value = no append 3nd integer
{
  if (ShowMF_Sets.count(mySet)>0)
  {
    const FABio::Format saved_format = FArrayBox::getFormat();
    FArrayBox::setFormat(ShowMF_Fab_Format);
    std::string DebugDir(ShowMF_Dir);
    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(DebugDir, 0755))
        amrex::CreateDirectoryFailed(DebugDir);
    ParallelDescriptor::Barrier();

    std::string junkname = name;
    if (lev>=0) {
      junkname = amrex::Concatenate(junkname+"_",lev,1);
    }
    if (iter>=0) {
      junkname = amrex::Concatenate(junkname+"_",iter,1);
    }
    if (step>=0) {
      junkname = amrex::Concatenate(junkname+"_",step,5);
    }
    junkname = DebugDir + "/" + junkname;

    if (ShowMF_Verbose>0) {
      amrex::Print() << "   ******************************  Debug: writing " 
		     << junkname << '\n';
    }

#if 0
    if (ShowMF_Check_Nans)
    {
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
      {
        //                AMREX_ASSERT(!mf[mfi].contains_nan(mfi.validbox(),0,mf.nComp()));
      }
    }
#endif
    VisMF::Write(mf,junkname);
    FArrayBox::setFormat(saved_format);
  }
}

PeleLM::FPLoc 
PeleLM::fpi_phys_loc (int p_bc)
{
  //
  // Location of data that FillPatchIterator returns at physical boundaries
  //
  if (p_bc == EXT_DIR || p_bc == HOEXTRAP || p_bc == FOEXTRAP)
  {
    return HT_Edge;
  }
  return HT_Center;
}

LM_Error_Value::LM_Error_Value (Real _min_time, Real _max_time, int _max_level)
    : lmef(0), value(0), min_time(_min_time), max_time(_max_time), max_level(_max_level)
{
}

LM_Error_Value::LM_Error_Value (LMEF _lmef,
                                Real _value, Real _min_time,
                                Real _max_time, int _max_level)
    : lmef(_lmef), lmef_box(0), value(_value), min_time(_min_time), max_time(_max_time), max_level(_max_level)
{
}

LM_Error_Value::LM_Error_Value (LMEF_BOX _lmef_box, const amrex::RealBox& _box, amrex::Real _min_time,
                                amrex::Real _max_time, int _max_level)
    : lmef(0), lmef_box(_lmef_box), min_time(_min_time), max_time(_max_time), box(_box), max_level(_max_level)
{
}

void
LM_Error_Value::tagCells(int* tag, const int* tlo, const int* thi,
                         const int* tagval, const int* clearval,
                         const Real* data, const int* dlo, const int* dhi,
                         const int* lo, const int* hi, const int* nvar,
                         const int* domain_lo, const int* domain_hi,
                         const Real* dx, const Real* xlo,
                         const Real* prob_lo, const Real* time,
                         const int* level) const
{
    AMREX_ASSERT(lmef);

    bool max_level_applies = ( (max_level < 0) || ( (max_level >= 0) && (*level < max_level) ) );
    bool valid_time_range = (min_time >= 0) && (max_time >= 0) && (min_time <= max_time);
    bool in_valid_time_range = valid_time_range && (*time >= min_time ) && (*time <= max_time);
    bool one_side_lo = !valid_time_range && (min_time >= 0) && (*time >= min_time);
    bool one_side_hi = !valid_time_range && (max_time >= 0) && (*time <= max_time);
    bool time_window_applies = in_valid_time_range || one_side_lo || one_side_hi || !valid_time_range;
 
    if (max_level_applies && time_window_applies)
    {
      lmef(tag, AMREX_ARLIM_ANYD(tlo), AMREX_ARLIM_ANYD(thi),
           tagval, clearval,
           data, AMREX_ARLIM_ANYD(dlo), AMREX_ARLIM_ANYD(dhi),
           lo, hi, nvar,
           domain_lo, domain_hi, dx, xlo, prob_lo, time, level, &value);
    }
}

void
LM_Error_Value::tagCells1(int* tag, const int* tlo, const int* thi,
                          const int* tagval, const int* clearval,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level) const
{
    AMREX_ASSERT(lmef_box);

    bool max_level_applies = ( (max_level < 0) || ( (max_level >= 0) && (*level < max_level) ) );
    bool valid_time_range = (min_time >= 0) && (max_time >= 0) && (min_time <= max_time);
    bool in_valid_time_range = valid_time_range && (*time >= min_time ) && (*time <= max_time);
    bool one_side_lo = !valid_time_range && (min_time >= 0) && (*time >= min_time);
    bool one_side_hi = !valid_time_range && (max_time >= 0) && (*time <= max_time);
    bool time_window_applies = in_valid_time_range || one_side_lo || one_side_hi || !valid_time_range;
 
    if (max_level_applies && time_window_applies)
    {
      lmef_box(tag, AMREX_ARLIM_ANYD(tlo), AMREX_ARLIM_ANYD(thi),
               tagval, clearval,
               AMREX_ZFILL(box.lo()), AMREX_ZFILL(box.hi()),
               lo, hi, domain_lo, domain_hi, dx, xlo, prob_lo, time, level);
    }
}

void
PeleLM::variableCleanUp ()
{

   NavierStokesBase::variableCleanUp();
   ShowMF_Sets.clear();
   auxDiag_names.clear();
   typical_values.clear();
   prob_parm.reset();
   ac_parm.reset();
}

PeleLM::PeleLM ()
{
   if (!init_once_done)
      init_once();

   if (!do_temp)
      amrex::Abort("do_temp MUST be true");

   if (!have_divu)
      amrex::Abort("have_divu MUST be true");

   // p_amb_old and p_amb_new contain the old-time and new-time
   // pressure at level 0.  For closed chamber problems they change over time.
   // set p_amb_old and new if they haven't been set yet
   // to the value in mod_Fvar_def.F90 set in PROB_F.F90
   // only the coarse level advance and the level 0-1 mac_sync 
   // can modify these later
   if (p_amb_old == -1.0)
   {
      p_amb_old = prob_parm->P_mean;
   }
   if (p_amb_new == -1.0)
   {
      p_amb_new = prob_parm->P_mean;
   }

   updateFluxReg = false;

   EdgeState              = 0;
   EdgeFlux               = 0;
   SpecDiffusionFluxn     = 0;
   SpecDiffusionFluxnp1   = 0;
   if (use_wbar) {
      SpecDiffusionFluxWbar  = 0;
   }
}

PeleLM::PeleLM (Amr&            papa,
                int             lev,
                const Geometry& level_geom,
                const BoxArray& bl,
                const DistributionMapping& dm,
                Real            time)
  :
  NavierStokesBase(papa,lev,level_geom,bl,dm,time)
{
  if (!init_once_done)
    init_once();

  if (!do_temp)
    amrex::Abort("do_temp MUST be true");

  if (!have_divu)
    amrex::Abort("have_divu MUST be true");

  // p_amb_old and p_amb_new contain the old-time and new-time
  // pressure at level 0.  For closed chamber problems they change over time.
  // set p_amb_old and new if they haven't been set yet
  // to the value in mod_Fvar_def.F90 set in PROB_F.F90
  // only the coarse level advance and the level 0-1 mac_sync 
  // can modify these later
  if (p_amb_old == -1.0)
  {
    p_amb_old = prob_parm->P_mean;
  }
  if (p_amb_new == -1.0)
  {
    p_amb_new = prob_parm->P_mean;
  }

  updateFluxReg = false;

  define_data();
  
}

PeleLM::~PeleLM ()
{
  ;
}

void
PeleLM::define_data ()
{
   const int nGrow       = 0;
#ifdef AMREX_USE_EB
   const int nGrowEdges  = 2; // We need 2 growth cells for the redistribution when using MOL EB
#else
   const int nGrowEdges  = 0; 
#endif
   const int nEdgeStates = desc_lst[State_Type].nComp();

   mTmpData.resize(mHtoTiterMAX);

   raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nEdgeStates, nGrowEdges)});
   EdgeState = raii_fbs.back()->get();

   raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nEdgeStates, nGrowEdges)});
   EdgeFlux  = raii_fbs.back()->get();
     
   if (NUM_SPECIES>0 && !unity_Le)
   {
     raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, NUM_SPECIES+3, nGrow)});
     SpecDiffusionFluxn   = raii_fbs.back()->get();

     raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, NUM_SPECIES+3, nGrow)});
     SpecDiffusionFluxnp1 = raii_fbs.back()->get();

#ifndef NDEBUG
     for ( int i = 0; i<AMREX_SPACEDIM; i++) {
       SpecDiffusionFluxnp1[i]->setVal(1.2345e40);
       SpecDiffusionFluxn[i]->setVal(1.2345e40);
     }
#endif
     
     if (use_wbar) {
        raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, NUM_SPECIES, nGrow)});
        SpecDiffusionFluxWbar = raii_fbs.back()->get();
     }
   }

   for (const auto& kv : auxDiag_names)
   {
      auxDiag[kv.first] = std::unique_ptr<MultiFab>(new MultiFab(grids,dmap,kv.second.size(),0));
      auxDiag[kv.first]->setVal(0.0);
   }

   // HACK for debugging
   if (level==0)
      stripBox = getStrip(geom);
}

void
PeleLM::init_once ()
{
   //
   // Computes the static variables unique to PeleLM.
   // Check that (some) things are set up correctly.
   //
   int dummy_State_Type;

   int tTemp = -1;
   int have_temp = isStateVariable("temp", dummy_State_Type, tTemp);
   AMREX_ALWAYS_ASSERT(tTemp == Temp);

   have_temp = have_temp && State_Type == dummy_State_Type;
   int tRhoH = -1;
   have_temp = have_temp && isStateVariable("rhoh", dummy_State_Type, tRhoH);
   AMREX_ALWAYS_ASSERT(tRhoH == RhoH);
   have_temp = have_temp && State_Type == dummy_State_Type;

   int tRhoRT = -1;
   have_rhort = isStateVariable("RhoRT", dummy_State_Type, tRhoRT);
   AMREX_ALWAYS_ASSERT(tRhoRT == RhoRT);
   have_rhort = have_rhort && State_Type == dummy_State_Type;

   if (!have_temp)
     amrex::Abort("PeleLM::init_once(): RhoH & Temp must both be the state");
     
   if (Temp < RhoH)
     amrex::Abort("PeleLM::init_once(): must have RhoH < Temp");
   //
   // Temperature must be non-conservative, rho*h must be conservative.
   //
   if (advectionType[Temp] == Conservative)
     amrex::Abort("PeleLM::init_once(): Temp must be non-conservative");

   if (advectionType[RhoH] != Conservative)
     amrex::Abort("PeleLM::init_once(): RhoH must be conservative");
   //
   // Species checks.
   //
   AMREX_ASSERT(Temp > RhoH && RhoH > Density);
   //
   // Here we want to count relative to Density instead of relative
   // to RhoH, so we can put RhoH after the species instead of before.  
   // This logic should work in both cases.
   last_spec  = first_spec + NUM_SPECIES - 1;
     
   for (int i = first_spec; i <= last_spec; i++)
      if (advectionType[i] != Conservative)
         amrex::Error("PeleLM::init_once: species must be conservative");
     
   int diffuse_spec = is_diffusive[first_spec];
   for (int i = first_spec+1; i <= last_spec; i++)
      if (is_diffusive[i] != diffuse_spec)
         amrex::Error("PeleLM::init_once: Le != 1; diffuse");
   //
   // Load integer pointers into Fortran common, reqd for proper ICs.
   //
   const int density = (int)Density;

   //
   // make space for typical values
   //
   typical_values.resize(NUM_STATE,-1); // -ve means don't use for anything
   typical_values[RhoH] = typical_RhoH_value_default;

   Vector<std::string> specNames;
   EOS::speciesNames(specNames);
   ParmParse pp("ht");
   for (int i=0; i<NUM_SPECIES; ++i) {
     const std::string ppStr = std::string("typValY_") + specNames[i];
     if (pp.countval(ppStr.c_str())>0) {
       pp.get(ppStr.c_str(),typical_values_FileVals[specNames[i]]);
     }
   }
   Vector<std::string> otherKeys = {"Temp", "RhoH", "Vel"};
   for (int i=0; i<otherKeys.size(); ++i) {
     const std::string ppStr(std::string("typVal_")+otherKeys[i]);
     if (pp.countval(ppStr.c_str())>0) {
       pp.get(ppStr.c_str(),typical_values_FileVals[otherKeys[i]]);
     }
   }

   //
   // Chemistry.
   //
   int ydot_good = RhoYdot_Type >= 0 && RhoYdot_Type <  desc_lst.size()
                                     && RhoYdot_Type != Divu_Type
                                     && RhoYdot_Type != State_Type;


   //
   // This is the minimum number of boxes per MPI proc that I want
   // when chopping up chemistry work.
   //
   pp.query("chem_box_chop_threshold",chem_box_chop_threshold);

   if (chem_box_chop_threshold <= 0)
   {
#ifdef AMREX_USE_OMP
     chem_box_chop_threshold = 8;
#else
     chem_box_chop_threshold = 4;
#endif
   }

   if (!ydot_good)
      amrex::Error("PeleLM::init_once(): need RhoYdot_Type if do_chemistry");

   if (desc_lst[RhoYdot_Type].nComp() < NUM_SPECIES)
      amrex::Error("PeleLM::init_once(): RhoYdot_Type needs NUM_SPECIES components");
   //
   // Enforce Le = 1, unless !unity_Le
   //
   if (unity_Le && (schmidt != prandtl) )
   {
      if (verbose)
        amrex::Print() << "**************** WARNING ***************\n"
                       << "PeleLM::init_once() : currently must have"
                       << "equal Schmidt and Prandtl numbers unless !unity_Le.\n"
                       << "Setting Schmidt = Prandtl\n"
                       << "**************** WARNING ***************\n";
      schmidt = prandtl;
   }
   //
   // We are done.
   //
   num_state_type = desc_lst.size();

   if (verbose)
      amrex::Print() << "PeleLM::init_once(): num_state_type = " << num_state_type << '\n';

   pp.query("plot_reactions",plot_reactions);
   if (plot_reactions)
   {
      auxDiag_names["REACTIONS"].resize(NUM_REACTIONS);
      amrex::Print() << "NUM_REACTIONS = "<< NUM_REACTIONS << '\n';
      for (int i = 0; i < auxDiag_names["REACTIONS"].size(); ++i)
         auxDiag_names["REACTIONS"][i] = amrex::Concatenate("R",i+1);
      amrex::Print() << "***** Make sure to increase amr.regrid_int !!!!!" << '\n';
   }

   pp.query("plot_consumption",plot_consumption);
   pp.query("plot_auxDiags",plot_consumption); // This is for backward comptibility - FIXME
   if (plot_consumption && consumptionName.size()>0)
   {
      auxDiag_names["CONSUMPTION"].resize(consumptionName.size());
      for (int j=0; j<consumptionName.size(); ++j)
      {
         auxDiag_names["CONSUMPTION"][j] = consumptionName[j] + "_ConsumptionRate";
      }
   }

   pp.query("plot_heat_release",plot_heat_release);
   if (plot_heat_release)
   {
      auxDiag_names["HEATRELEASE"].resize(1);
      auxDiag_names["HEATRELEASE"][0] = "HeatRelease";
   }

   init_once_done = 1;
}

void
PeleLM::restart (Amr&          papa,
                 std::istream& is,
                 bool          bReadSpecial)
{

   NavierStokesBase::restart(papa,is,bReadSpecial);

   define_data();

   bool running_sdc_from_strang_chk = false;

   if (running_sdc_from_strang_chk)
   {
      MultiFab& rYdot_old = get_old_data(RhoYdot_Type);
      MultiFab& rYdot_new = get_new_data(RhoYdot_Type);
      MultiFab& S_old = get_old_data(State_Type);
      MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(rYdot_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& R_old   = rYdot_old.array(mfi);
         auto const& R_new   = rYdot_new.array(mfi);
         auto const& rho_old = S_old.array(mfi,Density);
         auto const& rho_new = S_new.array(mfi,Density);

         amrex::ParallelFor(bx, NUM_SPECIES, [rho_old, rho_new, R_old, R_new]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            R_old(i,j,k,n) *= rho_old(i,j,k);
            R_new(i,j,k,n) *= rho_new(i,j,k);
         });
      }
   }

   // Deal with typical values
   set_typical_values(true);

   if (closed_chamber) {
      std::string line;
      std::string file=papa.theRestartFile();

      std::string File(file + "/PAMB");
      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream isp(fileCharPtrString, std::istringstream::in);

      // read in title line
      std::getline(isp, line);

      isp >> p_amb_old;
      p_amb_new = p_amb_old;
  }

  //
  // Initialize mixture fraction data.
  //
  init_mixture_fraction();

  //
  // Initialize active control
  //
  initActiveControl();

}

void
PeleLM::init_mixture_fraction()
{
      // Get default fuel and oxy tank composition
      Vector<std::string> specNames;
      EOS::speciesNames(specNames);
      amrex::Real YF[NUM_SPECIES], YO[NUM_SPECIES];
      for (int i=0; i<NUM_SPECIES; ++i) {
         YF[i] = 0.0;
         YO[i] = 0.0;  
         if (!specNames[i].compare("O2"))  YO[i] = 0.233;
         if (!specNames[i].compare("N2"))  YO[i] = 0.767;
         if (!specNames[i].compare(fuelName)) YF[i] = 1.0;
      }

      // Overwrite with user-defined value if provided in input file
      ParmParse pp("peleLM");
      std::string MFformat;
      int hasUserMF = pp.contains("mixtureFraction.format");
      if ( hasUserMF ) {
         pp.query("mixtureFraction.format", MFformat);
         if ( !MFformat.compare("Cantera")) {             // use a Cantera-like format with <SpeciesName>:<Value>, default in 0.0
            std::string MFCompoType;
            pp.query("mixtureFraction.type", MFCompoType);
            Vector<std::string> compositionIn;
            int entryCount = pp.countval("mixtureFraction.oxidTank");
            compositionIn.resize(entryCount);
            pp.getarr("mixtureFraction.oxidTank",compositionIn,0,entryCount);
            parseComposition(compositionIn, MFCompoType, YO);
            entryCount = pp.countval("mixtureFraction.fuelTank");
            compositionIn.resize(entryCount);
            pp.getarr("mixtureFraction.fuelTank",compositionIn,0,entryCount);
            parseComposition(compositionIn, MFCompoType, YF);
         } else if ( !MFformat.compare("RealList")) {     // use a list of Real. MUST contains an entry for each species in the mixture
            std::string MFCompoType;
            pp.query("mixtureFraction.type", MFCompoType);
            if ( !MFCompoType.compare("mass") ) {
               int entryCount = pp.countval("mixtureFraction.oxidTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               Vector<amrex::Real> compositionIn(NUM_SPECIES);
               pp.getarr("mixtureFraction.oxidTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  YO[i] = compositionIn[i];
               }
               entryCount = pp.countval("mixtureFraction.fuelTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               pp.getarr("mixtureFraction.fuelTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  YF[i] = compositionIn[i];
               }
            } else if ( !MFCompoType.compare("mole") ) {
               amrex::Real XF[NUM_SPECIES], XO[NUM_SPECIES];
               int entryCount = pp.countval("mixtureFraction.oxidTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               Vector<amrex::Real> compositionIn(NUM_SPECIES);
               pp.getarr("mixtureFraction.oxidTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  XO[i] = compositionIn[i];
               }
               entryCount = pp.countval("mixtureFraction.fuelTank");
               AMREX_ALWAYS_ASSERT(entryCount==NUM_SPECIES);
               pp.getarr("mixtureFraction.fuelTank",compositionIn,0,NUM_SPECIES);
               for (int i=0; i<NUM_SPECIES; ++i) {
                  XF[i] = compositionIn[i];
               }
               EOS::X2Y(XO, YO);
               EOS::X2Y(XF, YF);
            } else {
               Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
            }
         } else {
            Abort("Unknown mixtureFraction.format ! Should be 'Cantera' or 'RealList'");
         }
      }

      // Only interested in CHON -in that order. Compute Bilger weights
      amrex::Real atwCHON[4] = {0.0};
      EOS::atomic_weightsCHON(atwCHON);
      Beta_mix[0] = ( atwCHON[0] != 0.0 ) ? 2.0/atwCHON[0] : 0.0;
      Beta_mix[1] = ( atwCHON[1] != 0.0 ) ? 1.0/(2.0*atwCHON[1]) : 0.0;
      Beta_mix[2] = ( atwCHON[2] != 0.0 ) ? -1.0/atwCHON[2] : 0.0;
      Beta_mix[3] = 0.0;

      // Compute each species weight for the Bilger formulation based on elemental compo
      // Only interested in CHON -in that order.
      int ecompCHON[NUM_SPECIES*4];
      EOS::element_compositionCHON(ecompCHON);
      amrex::Gpu::DeviceVector<amrex::Real> mwt_v(NUM_SPECIES);
      amrex::Real* mwt_d = mwt_v.data();
      Box dumbx({AMREX_D_DECL(0,0,0)},{AMREX_D_DECL(0,0,0)});
      amrex::ParallelFor(dumbx, [mwt_d]
      AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      { 
         EOS::molecular_weight(mwt_d);
      });
      amrex::Real mwt[NUM_SPECIES];
#if AMREX_USE_GPU
      amrex::Gpu::dtoh_memcpy(mwt,mwt_d,sizeof(amrex::Real)*NUM_SPECIES);
#else
      std::memcpy(mwt,mwt_d,sizeof(amrex::Real)*NUM_SPECIES);
#endif
      Zfu = 0.0;
      Zox = 0.0;
      for (int i=0; i<NUM_SPECIES; ++i) {
         spec_Bilger_fact[i] = 0.0;
         for (int k = 0; k < 4; k++) {
            spec_Bilger_fact[i] += Beta_mix[k] * (ecompCHON[i*4 + k]*atwCHON[k]/mwt[i]);
         }
         Zfu += spec_Bilger_fact[i]*YF[i];
         Zox += spec_Bilger_fact[i]*YO[i];
      }

      mixture_fraction_ready = true;
}

void
PeleLM::set_typical_values(bool is_restart)
{
  if (level==0)
  {
    const int nComp = typical_values.size();

    if (is_restart) 
    {
      AMREX_ASSERT(nComp==NUM_STATE);

      for (int i=0; i<nComp; ++i)
        typical_values[i] = -1;
            
      if (ParallelDescriptor::IOProcessor())
      {
        const std::string tvfile = parent->theRestartFile() + "/" + typical_values_filename;
        std::ifstream tvis;
        tvis.open(tvfile.c_str(),std::ios::in|std::ios::binary);
                
        if (tvis.good())
        {
          FArrayBox tvfab;
          tvfab.readFrom(tvis);
          if (tvfab.nComp() != typical_values.size())
            amrex::Abort("Typical values file has wrong number of components");
          for (int i=0; i<typical_values.size(); ++i)
            typical_values[i] = tvfab.dataPtr()[i];
        }
      }

      ParallelDescriptor::ReduceRealMax(typical_values.dataPtr(),nComp); //FIXME: better way?
    }
    else  
    {
      Vector<const MultiFab*> Slevs(parent->finestLevel()+1);
      
      for (int lev = 0; lev <= parent->finestLevel(); lev++)
      {
        Slevs[lev] = &(getLevel(lev).get_new_data(State_Type));
      }

      auto scaleMax = VectorMax(Slevs,FabArrayBase::mfiter_tile_size,Density,NUM_STATE-AMREX_SPACEDIM,0);
      auto scaleMin = VectorMin(Slevs,FabArrayBase::mfiter_tile_size,Density,NUM_STATE-AMREX_SPACEDIM,0);      
      auto velMaxV = VectorMaxAbs(Slevs,FabArrayBase::mfiter_tile_size,0,AMREX_SPACEDIM,0);

      auto velMax = *max_element(std::begin(velMaxV), std::end(velMaxV));

      for (int i=0; i<NUM_STATE - AMREX_SPACEDIM; ++i) {
        typical_values[i + AMREX_SPACEDIM] = std::abs(scaleMax[i] - scaleMin[i]);
		  if ( typical_values[i + AMREX_SPACEDIM] <= 0.0)
			  typical_values[i + AMREX_SPACEDIM] = 0.5*std::abs(scaleMax[i] + scaleMin[i]);
      }

      for (int i=0; i<AMREX_SPACEDIM; ++i) {
        typical_values[i] = velMax;
      }

      AMREX_ALWAYS_ASSERT(typical_values[Density] > 0);
      typical_values[RhoH] = typical_values[RhoH] / typical_values[Density];
      for (int i=0; i<NUM_SPECIES; ++i) {
        typical_values[first_spec + i] = std::max(typical_values[first_spec + i]/typical_values[Density],
                                                  typical_Y_val_min);
      }
    }
    //
    // If typVals specified in inputs, these take precedence componentwise.
    //
    for (std::map<std::string,Real>::const_iterator it=typical_values_FileVals.begin(), 
           End=typical_values_FileVals.end();
         it!=End;
         ++it)
    {
      const int idx = getSpeciesIdx(it->first);
      if (idx>=0)
      {
        typical_values[first_spec+idx] = it->second;
      }
      else
      {
        if (it->first == "Temp")
          typical_values[Temp] = it->second;
        else if (it->first == "RhoH")
          typical_values[RhoH] = it->second;
        else if (it->first == "Vel")
        {
          for (int d=0; d<AMREX_SPACEDIM; ++d)
            typical_values[d] = std::max(it->second,typical_Y_val_min);
        }
      }
    }

    amrex::Print() << "Typical vals: " << '\n';
    amrex::Print() << "\tVelocity: ";
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
      amrex::Print() << typical_values[i] << ' ';
    }
    amrex::Print() << '\n';
    amrex::Print() << "\tDensity: " << typical_values[Density] << '\n';
    amrex::Print() << "\tTemp:    " << typical_values[Temp]    << '\n';
    amrex::Print() << "\tRhoH:    " << typical_values[RhoH]    << '\n';
    for (int i=0; i<NUM_SPECIES; ++i)
      {
        amrex::Print() << "\tY_" << spec_names[i] << ": " << typical_values[first_spec+i] <<'\n';
      }

#ifdef USE_SUNDIALS_PP
    if (use_typ_vals_chem) {
      amrex::Print() << "Using typical values for the absolute tolerances of the ode solver\n";
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
      {
      Vector<Real> typical_values_chem;
      typical_values_chem.resize(NUM_SPECIES+1);
      for (int i=0; i<NUM_SPECIES; ++i) {
	      typical_values_chem[i] = typical_values[first_spec+i] * typical_values[Density];
      }
      typical_values_chem[NUM_SPECIES] = typical_values[Temp];
      SetTypValsODE(typical_values_chem);
#ifndef AMREX_USE_GPU
      ReSetTolODE();
#endif  
      }
    }  
#endif

  }
}

void
PeleLM::reset_typical_values (const MultiFab& S)
{
  //
  // NOTE: Assumes that this level has valid data everywhere.
  //
  const int nComp = typical_values.size();

  AMREX_ASSERT(nComp == S.nComp());

  for (int i=0; i<nComp; ++i)
  {
    const Real thisMax = S.max(i);
    const Real thisMin = S.min(i);
    const Real newVal = std::abs(thisMax - thisMin);
    if (newVal > 0)
    {
      if ( (i>=first_spec && i<=last_spec) )
        typical_values[i] = newVal / typical_values[Density];
      else
        typical_values[i] = newVal;
    }
  }
  //
  // If typVals specified in inputs, these take precedence componentwise.
  //
  if (parent->levelSteps(0) == 0)
  {
    for (std::map<std::string,Real>::const_iterator it=typical_values_FileVals.begin(), 
           End=typical_values_FileVals.end();
         it!=End;
         ++it)
    {
      const int idx = getSpeciesIdx(it->first);
      if (idx>=0)
      {
        typical_values[first_spec+idx] = it->second;
      }
      else
      {
        if (it->first == "Temp")
          typical_values[Temp] = it->second;
        else if (it->first == "RhoH")
          typical_values[RhoH] = it->second;
        else if (it->first == "Vel")
        {
          for (int d=0; d<AMREX_SPACEDIM; ++d)
            typical_values[d] = std::max(it->second,typical_Y_val_min);
        }
      }
    }
  }

  amrex::Print() << "New typical vals: " << '\n';
  amrex::Print() << "\tVelocity: ";
  for (int i=0; i<AMREX_SPACEDIM; ++i) {
    amrex::Print() << typical_values[i] << ' ';
  }
  amrex::Print() << '\n';
  amrex::Print() << "\tDensity:  " << typical_values[Density] << '\n';
  amrex::Print() << "\tTemp:     " << typical_values[Temp]    << '\n';
  amrex::Print() << "\tRhoH:     " << typical_values[RhoH]    << '\n';
  for (int i=0; i<NUM_SPECIES; ++i)
    {
      amrex::Print() << "\tY_" << spec_names[i] << ": " << typical_values[first_spec+i] <<'\n';
    }
}

Real
PeleLM::estTimeStep ()
{
   Real estdt = NavierStokesBase::estTimeStep();

   if (fixed_dt > 0.0 || !divu_ceiling)
     //
     // The n-s function did the right thing in this case.
     //
     return estdt;

   const Real strt_time = ParallelDescriptor::second();

   Real ns_estdt = estdt;
   Real divu_dt = 1.0e20;

   const int   nGrow   = 1;
   const Real  cur_time = state[State_Type].curTime();
   std::unique_ptr<MultiFab> DivU (getDivCond(0,cur_time));
   int  divu_check_flag = divu_ceiling;
   Real divu_dt_fac     = divu_dt_factor;
   Real rho_min         = min_rho_divu_ceiling;
   const auto& dxinv    = geom.InvCellSizeArray();

   FillPatchIterator U_fpi(*this,*DivU,nGrow,cur_time,State_Type,Xvel,AMREX_SPACEDIM);
   MultiFab& Umf=U_fpi.get_mf();

   if ( divu_ceiling == 1 ) {
      divu_dt = amrex::ReduceMin(rho_ctime, *DivU, 0,
                                 [divu_check_flag,divu_dt_fac,rho_min,dxinv]
      AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& rho,
                                            Array4<Real const> const& divu ) noexcept -> Real
      {   
         using namespace amrex::literals;
         const auto lo = amrex::lbound(bx);
         const auto hi = amrex::ubound(bx);
         amrex::Real dt = 1.e37_rt;
         for       (int k = lo.z; k <= hi.z; ++k) {
            for    (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  Real dtcell = est_divu_dt_1(i, j, k, divu_check_flag, divu_dt_fac, rho_min, dxinv,
                                              rho, divu );
                  dt = amrex::min(dt,dtcell);
               }
            }
         }
         return dt;
      });
   } else if ( divu_ceiling == 2 ) {
      divu_dt = amrex::ReduceMin(rho_ctime, Umf, *DivU, 0,
                                 [divu_check_flag,divu_dt_fac,rho_min,dxinv]
      AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& rho,
                                            Array4<Real const> const& vel,
                                            Array4<Real const> const& divu ) noexcept -> Real
      {   
         using namespace amrex::literals;
         const auto lo = amrex::lbound(bx);
         const auto hi = amrex::ubound(bx);
         amrex::Real dt = 1.e37_rt;
         for       (int k = lo.z; k <= hi.z; ++k) {
            for    (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  Real dtcell = est_divu_dt_2(i, j, k, divu_check_flag, divu_dt_fac, rho_min, dxinv,
                                              rho, vel, divu );
                  dt = amrex::min(dt,dtcell);
               }
            }
         }
         return dt;
      });
   } else if ( divu_ceiling == 3 ) {
      amrex::Abort("divu_ceiling == 3 currently not available. If amrex::ReduceMin has been updated then this should be fixed");
   } else {
      amrex::Abort("Unknown divu_ceiling value. It should be between 0 and 3.");
   }

   ParallelDescriptor::ReduceRealMin(divu_dt);

   if (verbose)
   {
      amrex::Print() << "PeleLM::estTimeStep(): estdt, divu_dt = "
                     << estdt << ", " << divu_dt << '\n';
   }

   estdt = std::min(estdt, divu_dt);

   if (estdt < ns_estdt && verbose)
      amrex::Print() << "PeleLM::estTimeStep(): timestep reduced from "
                     << ns_estdt << " to " << estdt << '\n';

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::estTimeStep(): lev: " << level
                     << ", time: " << run_time << '\n';
   }

   return estdt;
}

void
PeleLM::checkTimeStep (Real dt)
{
   if (fixed_dt > 0.0 || !divu_ceiling)
      return;

   const int   nGrow    = 1;
   const Real  cur_time  = state[State_Type].curTime();
   MultiFab*   DivU      = getDivCond(0,cur_time);

   FillPatchIterator U_fpi(*this,*DivU,nGrow,cur_time,State_Type,Xvel,AMREX_SPACEDIM);
   MultiFab& Umf=U_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(Umf,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      const Box&  bx       = mfi.tilebox();
      auto const& rho      = rho_ctime.array(mfi);
      auto const& vel      = Umf.array(mfi);
      auto const& divu     = DivU->array(mfi);
      auto const& vol      = volume.const_array(mfi);
      AMREX_D_TERM(auto const& areax = (area[0]).const_array(mfi);,
                   auto const& areay = (area[1]).const_array(mfi);,
                   auto const& areaz = (area[2]).const_array(mfi););
      int  divu_check_flag = divu_ceiling;
      Real divu_dt_fac     = divu_dt_factor;
      Real rho_min         = min_rho_divu_ceiling;
      const auto dxinv     = geom.InvCellSizeArray();

      amrex::ParallelFor(bx, [rho, vel, divu, vol, AMREX_D_DECL(areax,areay,areaz),
                              divu_check_flag, divu_dt_fac, rho_min, dxinv, dt]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         check_divu_dt(i, j, k, divu_check_flag, divu_dt_fac, rho_min, dxinv,
                       rho, vel, divu, vol, AMREX_D_DECL(areax,areay,areaz), dt); 
      });
   }
   delete DivU;
}

void
PeleLM::setTimeLevel (Real time,
                      Real dt_old,
                      Real dt_new)
{
   NavierStokesBase::setTimeLevel(time, dt_old, dt_new);

   state[RhoYdot_Type].setTimeLevel(time,dt_old,dt_new);

   state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

//
// This (minus the NEWMECH stuff) is copied from NavierStokes.cpp
//

void
PeleLM::initData ()
{
  //
  // Initialize the state and the pressure.
  //
  int         ns       = NUM_STATE - AMREX_SPACEDIM;
  const Real* dx       = geom.CellSize();
  MultiFab&   S_new    = get_new_data(State_Type);
  MultiFab&   P_new    = get_new_data(Press_Type);
  const Real  cur_time = state[State_Type].curTime();

#ifdef AMREX_USE_NEWMECH
  //
  // This code has a few drawbacks.  It assumes that the physical
  // domain size of the current problem is the same as that of the
  // one that generated the pltfile.  It also assumes that the pltfile
  // has at least as many levels as does the current problem.  If
  // either of these are false this code is likely to core dump.
  //
  ParmParse pp("ht");

  std::string pltfile;
  pp.query("pltfile", pltfile);
  if (pltfile.empty())
    amrex::Abort("You must specify `pltfile'");
  if (verbose)
    amrex::Print() << "initData: reading data from: " << pltfile << '\n';

  DataServices::SetBatchMode();
  FileType fileType(NEWPLT);
  DataServices dataServices(pltfile, fileType);

  if (!dataServices.AmrDataOk())
    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
  AmrData&                  amrData     = dataServices.AmrDataRef();
  Vector<std::string> names;
  EOS::speciesNames(names);
  Vector<std::string>        plotnames   = amrData.PlotVarNames();

  int idT = -1, idX = -1;
  for (int i = 0; i < plotnames.size(); ++i)
  {
    if (plotnames[i] == "temp")       idT = i;
    if (plotnames[i] == "x_velocity") idX = i;
  }
  //
  // In the plotfile the mass fractions directly follow the velocities.
  //
  int idSpec = idX + AMREX_SPACEDIM;

  for (int i = 0; i < AMREX_SPACEDIM; i++)
  {
    amrData.FillVar(S_new, level, plotnames[idX+i], Xvel+i);
    amrData.FlushGrids(idX+i);
  }
  amrData.FillVar(S_new, level, plotnames[idT], Temp);
  amrData.FlushGrids(idT);

  for (int i = 0; i < NUM_SPECIES; i++)
  {
    amrData.FillVar(S_new, level, plotnames[idSpec+i], first_spec+i);
    amrData.FlushGrids(idSpec+i);
  }

  if (verbose) amrex::Print() << "initData: finished init from pltfile" << '\n';
#endif

  const auto geomdata = geom.data();
  S_new.setVal(0.0);
  P_new.setVal(0.0);

  ProbParm const* lprobparm = prob_parm.get();
  PmfData const* lpmfdata = pmf_data_g;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
      const Box& box = mfi.validbox();
      auto sfab = S_new.array(mfi);

      amrex::ParallelFor(box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
#ifdef AMREX_USE_NEWMECH
        amrex::Abort("USE_NEWMECH feature no longer working and has to be fixed/redone");
#else
        pelelm_initdata(i, j, k, sfab, geomdata, *lprobparm, lpmfdata);
#endif
      });
  }

  showMFsub("1D",S_new,stripBox,"1D_S",level);
 
// Here we save a reference state vector to apply it later to covered cells
// in order to avoid non-physical values after diffusion solves
// First we have to put Pnew in S_new so as to not impose NaNs for covered cells
  MultiFab::Copy(S_new,P_new,0,RhoRT,1,1);

  //
  // Initialize GradP
  //
  computeGradP(state[Press_Type].curTime());

#ifdef AMREX_USE_EB
  set_body_state(S_new);
#endif
  

#ifdef AMREX_USE_VELOCITY
  //
  // We want to add the velocity from the supplied plotfile
  // to what we already put into the velocity field via FORT_INITDATA.
  //
  // This code has a few drawbacks.  It assumes that the physical
  // domain size of the current problem is the same as that of the
  // one that generated the pltfile.  It also assumes that the pltfile
  // has at least as many levels (with the same refinement ratios) as does
  // the current problem.  If either of these are false this code is
  // likely to core dump.
  //
  ParmParse pp("ht");

  std::string velocity_plotfile;
  pp.query("velocity_plotfile", velocity_plotfile);

  if (!velocity_plotfile.empty())
  {
    if (verbose)
      amrex::Print() << "initData: reading data from: " << velocity_plotfile << '\n';

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(velocity_plotfile, fileType);

    if (!dataServices.AmrDataOk())
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData&                  amrData   = dataServices.AmrDataRef();
    Vector<std::string>        plotnames = amrData.PlotVarNames();

    if (amrData.FinestLevel() < level)
      amrex::Abort("initData: not enough levels in plotfile");

    if (amrData.ProbDomain()[level] != Domain())
      amrex::Abort("initData: problem domains do not match");
    
    int idX = -1;
    for (int i = 0; i < plotnames.size(); ++i)
      if (plotnames[i] == "x_velocity") idX = i;

    if (idX == -1)
      amrex::Abort("Could not find velocity fields in supplied velocity_plotfile");

    MultiFab tmp(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
      amrData.FillVar(tmp, level, plotnames[idX+i], 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(tmp,true); mfi.isValid(); ++mfi) {
        S_new[mfi].plus(tmp[mfi], mfi.tilebox(), 0, Xvel+i, 1);
      }
      amrData.FlushGrids(idX+i);
    }

    if (verbose)
      amrex::Print() << "initData: finished init from velocity_plotfile" << '\n';
  }
#endif /*AMREX_USE_VELOCITY*/

  make_rho_prev_time();
  make_rho_curr_time();
  //
  // Initialize other types.
  //
  initDataOtherTypes();

  //
  // Load typical values for each state component
  //
  set_typical_values(false);
  
  //
  // Initialize mixture fraction data.
  //
  init_mixture_fraction();

  //
  // Initialize active control
  //
  initActiveControl();

  //
  // Initialize divU
  //
  if (have_divu)
  {
    const Real dt       = 1.0;
    const Real dtin     = -1.0; // Dummy value denotes initialization.
    const Real tnp1 = state[Divu_Type].curTime();
    MultiFab&  Divu_new = get_new_data(Divu_Type);

    state[State_Type].setTimeLevel(tnp1,dt,dt);

    calcDiffusivity(tnp1);

    calc_divu(tnp1,dtin,Divu_new);
  }

  is_first_step_after_regrid = false;
  old_intersect_new          = grids;

#ifdef AMREX_PARTICLES
  NavierStokesBase::initParticleData();
#endif
}

void
PeleLM::initDataOtherTypes ()
{
  // Fill RhoH component using EOS function explicitly
  const Real tnp1  = state[State_Type].curTime();
  compute_rhohmix(tnp1,get_new_data(State_Type),RhoH);

  // Set initial omegadot = 0
  get_new_data(RhoYdot_Type).setVal(0);

  // Put something reasonable into the FuncCount variable
  get_new_data(FuncCount_Type).setVal(1);

  // Fill the thermodynamic pressure into the state
  setThermoPress(state[State_Type].curTime());
}

void
PeleLM::compute_instantaneous_reaction_rates (MultiFab&       R,
                                              const MultiFab& S,
                                              Real            time,
                                              int             nGrow,
                                              HowToFillGrow   how)
{
   if (hack_nochem)
   {
      R.setVal(0.0,0,NUM_SPECIES,R.nGrow());
      return;
   }

   const Real strt_time = ParallelDescriptor::second();

   AMREX_ASSERT((nGrow==0)  ||  (how == HT_ZERO_GROW_CELLS) || (how == HT_EXTRAP_GROW_CELLS));

   if ((nGrow>0) && (how == HT_ZERO_GROW_CELLS))
   {
       R.setBndry(0,0,NUM_SPECIES);
   }

// TODO: the mask is not used right now.
//       we need something that's zero in covered and 1 otherwise to avoid an if test
//       in the kernel.
   amrex::FabArray<amrex::BaseFab<int>> maskMF;
   maskMF.define(grids, dmap, 1, 0);
   maskMF.setVal(1.0);
#ifdef AMREX_USE_EB
   maskMF.copy(ebmask);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhoY    = S.array(mfi,first_spec);
      auto const& rhoH    = S.array(mfi,RhoH);
      auto const& T       = S.array(mfi,Temp);
      auto const& mask    = maskMF.array(mfi);
      auto const& rhoYdot = R.array(mfi);

      amrex::ParallelFor(bx, [rhoY, rhoH, T, mask, rhoYdot]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         reactionRateRhoY( i, j, k, rhoY, rhoH, T, mask,
                           rhoYdot );
      });
   }

   if ((nGrow>0) && (how == HT_EXTRAP_GROW_CELLS))
   {
      R.FillBoundary(0,NUM_SPECIES, geom.periodicity());
      AMREX_ASSERT(R.nGrow() == 1);
      Extrapolater::FirstOrderExtrap(R, geom, 0, NUM_SPECIES);
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::compute_instantaneous_reaction_rates(): lev: " << level
                     << ", time: " << run_time << '\n';
   }
} 

void
PeleLM::init (AmrLevel& old)
{
   NavierStokesBase::init(old);

   PeleLM* oldht    = (PeleLM*) &old;
   const Real    tnp1 = oldht->state[State_Type].curTime();

   //
   // Get good version of rhoYdot and Function count
   //
   MultiFab& Ydot = get_new_data(RhoYdot_Type);
   MultiFab& FuncCount = get_new_data(FuncCount_Type);
   RhoH_to_Temp(get_new_data(State_Type));

   FillPatchIterator Ydotfpi(*oldht,Ydot,Ydot.nGrow(),tnp1,RhoYdot_Type,0,NUM_SPECIES);
   const MultiFab& Ydot_old = Ydotfpi.get_mf();

   FillPatchIterator FctCntfpi(*oldht,FuncCount,FuncCount.nGrow(),tnp1,FuncCount_Type,0,1);
   const MultiFab& FuncCount_old = FctCntfpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
   for (MFIter mfi(Ydot_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx         = mfi.tilebox();
      auto const& rhoYdot_n = Ydot.array(mfi); 
      auto const& rhoYdot_o = Ydot_old.array(mfi); 
      auto const& FctCnt_n  = FuncCount.array(mfi); 
      auto const& FctCnt_o  = FuncCount_old.array(mfi); 
      amrex::ParallelFor(bx, [rhoYdot_n, rhoYdot_o, FctCnt_n, FctCnt_o]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         for (int n = 0; n < NUM_SPECIES; n++) {
            rhoYdot_n(i,j,k,n) = rhoYdot_o(i,j,k,n);
         }
         FctCnt_n(i,j,k) = FctCnt_o(i,j,k);
      });
   }
}

//
// Inits the data on a new level that did not exist before regridding.
//
void
PeleLM::init ()
{
   NavierStokesBase::init();

   PeleLM& coarser  = getLevel(level-1);
   const Real    tnp1 = coarser.state[State_Type].curTime();
   //
   // Get best ydot data.
   //
   FillCoarsePatch(get_new_data(RhoYdot_Type),0,tnp1,RhoYdot_Type,0,NUM_SPECIES);

   RhoH_to_Temp(get_new_data(State_Type));
   get_new_data(State_Type).setBndry(1.e30);
   showMF("sdc",get_new_data(State_Type),"sdc_new_state",level);

   FillCoarsePatch(get_new_data(FuncCount_Type),0,tnp1,FuncCount_Type,0,1);
}

void
PeleLM::post_timestep (int crse_iteration)
{
  BL_PROFILE("PLM::post_timestep()");

  NavierStokesBase::post_timestep(crse_iteration);

  if (plot_reactions && level == 0)
  {
    const int Ndiag = auxDiag["REACTIONS"]->nComp();
    //
    // Multiply by the inverse of the coarse timestep.
    //
    const Real factor = 1.0 / crse_dt;

    for (int i = parent->finestLevel(); i >= 0; --i)
      getLevel(i).auxDiag["REACTIONS"]->mult(factor);

    for (int i = parent->finestLevel(); i > 0; --i)
    {
      PeleLM& clev = getLevel(i-1);
      PeleLM& flev = getLevel(i);

      MultiFab& Ydot_crse = *(clev.auxDiag["REACTIONS"]);
      MultiFab& Ydot_fine = *(flev.auxDiag["REACTIONS"]);

      clev.average_down(Ydot_fine, Ydot_crse, 0, Ndiag);
    }
  }
}

void
PeleLM::post_restart ()
{
   NavierStokesBase::post_restart();

   make_rho_prev_time();
   make_rho_curr_time();

   Real dummy  = 0;
   int MyProc  = ParallelDescriptor::MyProc();
   int step    = parent->levelSteps(0);
   int is_restart = 1;

   activeControl(step,is_restart,0.0,crse_dt);
}

void
PeleLM::postCoarseTimeStep (Real cumtime)
{
  //
  // postCoarseTimeStep() is only called by level 0.
  //
  AMREX_ASSERT(level == 0);
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
PeleLM::post_regrid (int lbase,
                     int new_finest)
{
   NavierStokesBase::post_regrid(lbase, new_finest);
   //
   // FIXME: This may be necessary regardless, unless the interpolation
   //        to fine from coarse data preserves rho=sum(rho.Y)
   //
   if (!do_set_rho_to_species_sum) return;

   int       minzero = 0;           // Flag to clip the species mass density to zero before summing it up to build rho.
   if (parent->levelSteps(0)>0 && level>lbase) {
      MultiFab& Snew = get_new_data(State_Type);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
         for (MFIter mfi(Snew,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& rho     = Snew.array(mfi,Density);
            auto const& rhoY    = Snew.array(mfi,first_spec);
            if (minzero) {
               amrex::ParallelFor(bx, [rhoY]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  fabMinMax( i, j, k, NUM_SPECIES, 0.0, Real_MAX, rhoY);
               });
            }
            amrex::ParallelFor(bx, [rho, rhoY]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               rho(i,j,k) = 0.0;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  rho(i,j,k) += rhoY(i,j,k,n);
               }
            });
         }
      }
   }
   make_rho_curr_time();
}

void
PeleLM::checkPoint (const std::string& dir,
                    std::ostream&      os,
                    VisMF::How         how,
                    bool               dump_old)
{
  NavierStokesBase::checkPoint(dir,os,how,dump_old);

  if (level == 0)
  {
    if (ParallelDescriptor::IOProcessor())
    {
      const std::string tvfile = dir + "/" + typical_values_filename;
      std::ofstream tvos;
      tvos.open(tvfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
      if (!tvos.good())
        amrex::FileOpenFailed(tvfile);
      Box tvbox(IntVect(),(NUM_STATE-1)*amrex::BASISV(0));
      int nComp = typical_values.size();
      FArrayBox tvfab(tvbox,nComp);
      for (int i=0; i<nComp; ++i) {
        tvfab.dataPtr()[i] = typical_values[i];
      }
      tvfab.writeOn(tvos);
    }
  }

  if (closed_chamber) {

      VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

      if (ParallelDescriptor::IOProcessor()) {

          std::ofstream PAMBFile;
          PAMBFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
          std::string PAMBFileName(dir + "/PAMB");
          PAMBFile.open(PAMBFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);

          if( !PAMBFile.good()) {
              amrex::FileOpenFailed(PAMBFileName);
          }
          
          PAMBFile.precision(17);
          
          int step = parent->levelSteps(0);

          // write out title line
          PAMBFile << "Writing p_amb to checkpoint\n";
          if (step == 0) {
              PAMBFile << p_amb_old << "\n";
          }
          else {
              PAMBFile << p_amb_new << "\n";
          }
      }
  }
}

void
PeleLM::post_init (Real stop_time)
{
  if (level > 0)
    //
    // Nothing to sync up at level > 0.
    //
    return;

  const Real strt_time    = ParallelDescriptor::second();
  const Real tnp1     = state[State_Type].curTime();
  const int  finest_level = parent->finestLevel();
  Real        dt_init     = 0.0;
  Real Sbar;

  Vector<Real> dt_save(finest_level+1);
  Vector<int>  nc_save(finest_level+1);
  Real        dt_init2 = 0.0;
  Vector<Real> dt_save2(finest_level+1);
  Vector<int>  nc_save2(finest_level+1);

  // ensure system is solvable by creating deltaS = S - Sbar
  if (closed_chamber == 1)
  {

    // ensure divu is average down so computing deltaS = S - Sbar works for multilevel
    for (int k = finest_level-1; k >= 0; k--)
    {
      PeleLM&   fine_lev = getLevel(k+1);
      PeleLM&   crse_lev = getLevel(k);

      MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
      MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

      crse_lev.average_down(Divu_fine, Divu_crse, 0, 1);
    }

    // compute number of cells
    Real num_cells = grids.numPts();

    // compute Sbar and subtract from S
    for (int lev = 0; lev <= finest_level; lev++)
    {
      // pointer to S
      MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
      if (lev == 0)
      {
        Sbar = divu_lev.sum() / num_cells;
      }
      divu_lev.plus(-Sbar,0,1);
    }
  }

  //
  // Ensure state is consistent, i.e. velocity field satisfies initial
  // estimate of constraint, coarse levels are fine level averages, pressure
  // is zero.
  //
  post_init_state();

  if (closed_chamber == 1)
  {
    // restore S
    for (int lev = 0; lev <= finest_level; lev++)
    {
      MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
      divu_lev.plus(Sbar,0,1);
    }
  }

  //
  // Estimate the initial timestepping.
  //
  post_init_estDT(dt_init,nc_save,dt_save,stop_time);
  //
  // Better estimate needs dt to estimate divu
  //
  const bool do_iter        = do_init_proj && projector;
  const int  init_divu_iter = do_iter ? num_divu_iters : 0;

  if (verbose) {
    amrex::Print() << "doing num_divu_iters = " << num_divu_iters << '\n';
  }
  for (int iter = 0; iter < init_divu_iter; ++iter)
  {
    //
    // Update species destruction rates in each level but not state.
    //
    if (NUM_SPECIES > 0)
    {

      for (int k = 0; k <= finest_level; k++)
      {
        
        MultiFab& S_new = getLevel(k).get_new_data(State_Type);

        //
        // Don't update S_new in this strang_chem() call ...
        //
        MultiFab S_tmp(S_new.boxArray(),S_new.DistributionMap(),S_new.nComp(),0,MFInfo(),getLevel(k).Factory());
        MultiFab Forcing_tmp(S_new.boxArray(),S_new.DistributionMap(),NUM_SPECIES+1,0,MFInfo(),getLevel(k).Factory());
        Forcing_tmp.setVal(0);

        getLevel(k).advance_chemistry(S_new,S_tmp,dt_save[k]/2.0,Forcing_tmp,0);
      }
    }
    //
    // Recompute the velocity to obey constraint with chemistry and
    // divqrad and then average that down.
    //
    if (NUM_SPECIES > 0)
    {
      for (int k = 0; k <= finest_level; k++)
      {
        MultiFab&  Divu_new = getLevel(k).get_new_data(Divu_Type);
        getLevel(k).calc_divu(tnp1,dt_save[k],Divu_new);
      }


      if (!hack_noavgdivu)
      {
        for (int k = finest_level-1; k >= 0; k--)
        {
          PeleLM&   fine_lev = getLevel(k+1);
          PeleLM&   crse_lev = getLevel(k);

          MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
          MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

	  crse_lev.average_down(Divu_fine, Divu_crse, 0, 1);
        }
      }
      //
      // Recompute the initial velocity field based on this new constraint
      //
      const Real divu_time = state[Divu_Type].curTime();

      int havedivu = 1;

      // ensure system is solvable by creating deltaS = S - Sbar
      if (closed_chamber == 1)
      {
        // compute number of cells
        Real num_cells = grids.numPts();

        // compute Sbar and subtract from S
        for (int lev = 0; lev <= finest_level; lev++)
        {
          // pointer to S
          MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
          if (lev == 0)
          {
            Sbar = divu_lev.sum() / num_cells;
          }
          divu_lev.plus(-Sbar,0,1);
        }
      }

      projector->initialVelocityProject(0,divu_time,havedivu,1);

      if (closed_chamber == 1)
      {
        // restore S
        for (int lev = 0; lev <= finest_level; lev++)
        {
          MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
          divu_lev.plus(Sbar,0,1);
        }
      }

      //
      // Average down the new velocity
      //
      for (int k = finest_level-1; k >= 0; k--)
      {
        PeleLM&   fine_lev = getLevel(k+1);
        PeleLM&   crse_lev = getLevel(k);
        MultiFab&       S_fine  = getLevel(k+1).get_new_data(State_Type);
        MultiFab&       S_crse  = getLevel(k).get_new_data(State_Type);

	crse_lev.average_down(S_fine, S_crse, Xvel, AMREX_SPACEDIM);
      }
    }
    //
    // Estimate the initial timestepping again, using new velocity
    // (necessary?) Need to pass space to save dt, nc, but these are
    // hacked, just pass something.
    //
    // reset Ncycle to nref...
    //
    parent->setNCycle(nc_save);
    post_init_estDT(dt_init2, nc_save2, dt_save2, stop_time);
    //
    // Compute dt_init,dt_save as the minimum of the values computed
    // in the calls to post_init_estDT
    // Then setTimeLevel and dt_level to these values.
    //
    dt_init = std::min(dt_init,dt_init2);
    Vector<Real> dt_level(finest_level+1,dt_init);

    parent->setDtLevel(dt_level);
    for (int k = 0; k <= finest_level; k++)
    {
      dt_save[k] = std::min(dt_save[k],dt_save2[k]);
      getLevel(k).setTimeLevel(tnp1,dt_init,dt_init);
    }
  } // end divu_iters

  //
  // Initialize the pressure by iterating the initial timestep.
  //
  post_init_press(dt_init, nc_save, dt_save);
  //
  // Compute the initial estimate of conservation.
  //
  if (sum_interval > 0)
    sum_integrated_quantities();

  if (verbose)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    amrex::Print() << "PeleLM::post_init(): lev: " << level << ", time: " << run_time << '\n';
  }
}

void
PeleLM::sum_integrated_quantities ()
{
  const int finest_level = parent->finestLevel();
  const Real tnp1        = state[State_Type].curTime();

  if (verbose)
  {
    Real mass = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
       mass += getLevel(lev).volWgtSum("density",tnp1);

    Print() << "TIME= " << tnp1 << " MASS= " << mass;
  }
  
  if (getSpeciesIdx(fuelName) >= 0)
  {
    Real fuelmass = 0.0;
    std::string fuel = "rho.Y(" + fuelName + ")";
    for (int lev = 0; lev <= finest_level; lev++) {
       fuelmass += getLevel(lev).volWgtSum(fuel,tnp1);
    }
    if (verbose) amrex::Print() << " FUELMASS= " << fuelmass;
  }

  if (verbose)
  {
    if (getSpeciesIdx(productName) >= 0)
    {
      Real productmass = 0.0;
      std::string product = "rho.Y(" + productName + ")";
      for (int lev = 0; lev <= finest_level; lev++) {
        productmass += getLevel(lev).volWgtSum(product,tnp1);
      }	  
      Print() << " PRODUCTMASS= " << productmass;
    }
    Print() << '\n';
  }

  int step    = parent->levelSteps(0);
  int is_restart = 0;
  activeControl(step,is_restart,tnp1,crse_dt);

  {
     Real rho_h    = 0.0;
     Real rho_temp = 0.0;
     for (int lev = 0; lev <= finest_level; lev++)
     {
       rho_temp += getLevel(lev).volWgtSum("rho_temp",tnp1);
       rho_h     = getLevel(lev).volWgtSum("rhoh",tnp1);
       if (verbose) Print() << "TIME= " << tnp1 << " LEV= " << lev << " RHOH= " << rho_h << '\n';
     }
     if (verbose) Print() << "TIME= " << tnp1 << " RHO*TEMP= " << rho_temp << '\n';
  }
  
  if (verbose)
  {
    int old_prec = std::cout.precision(15);
    
    Vector<const MultiFab*> Slevs(parent->finestLevel()+1);  
    for (int lev = 0; lev <= parent->finestLevel(); lev++) {
      Slevs[lev] = &(getLevel(lev).get_new_data(State_Type));
    }
    auto Smin = VectorMin(Slevs,FabArrayBase::mfiter_tile_size,0,NUM_STATE,0);
    auto Smax = VectorMax(Slevs,FabArrayBase::mfiter_tile_size,0,NUM_STATE,0);

    Print() << "TIME= " << tnp1 << " min,max temp = " << Smin[Temp] << ", " << Smax[Temp] << '\n';
    for (int n=0; n<AMREX_SPACEDIM; ++n) {
      std::string str = (n==0 ? "xvel" : (n==1 ? "yvel" : "zvel") );
      Print() << "TIME= " << tnp1  << " min,max "<< str << "  = " << Smin[Xvel+n] << ", " << Smax[Xvel+n] << '\n';
    }

    if (NUM_SPECIES > 0)
    {
      Real min_sum, max_sum;
      for (int lev = 0; lev <= finest_level; lev++) {
        auto mf = getLevel(lev).derive("rhominsumrhoY",tnp1,0);
        Real this_min = mf->min(0);
        Real this_max = mf->max(0);
        if (lev==0) {
          min_sum = this_min;
          max_sum = this_max;
        } else {
          min_sum = std::min(this_min,min_sum);
          max_sum = std::max(this_max,max_sum);
        }
      }
            
      Print() << "TIME= " << tnp1 
              << " min,max rho-sum rho Y_l = "
              << min_sum << ", " << max_sum << '\n';
      
      for (int lev = 0; lev <= finest_level; lev++)
      {
        auto mf = getLevel(lev).derive("sumRhoYdot",tnp1,0);
        Real this_min = mf->min(0);
        Real this_max = mf->max(0);
        if (lev==0) {
          min_sum = this_min;
          max_sum = this_max;
        } else {
          min_sum = std::min(this_min,min_sum);
          max_sum = std::max(this_max,max_sum);
        }
      }
            
      Print() << "TIME= " << tnp1 
              << " min,max sum RhoYdot = "
              << min_sum << ", " << max_sum << '\n';
    }
    std::cout.precision(old_prec);
  }
}
    
void
PeleLM::post_init_press (Real&        dt_init,
                         Vector<int>&  nc_save,
                         Vector<Real>& dt_save)
{
  const int  nState          = desc_lst[State_Type].nComp();
  const int  nGrow           = 0;
  const Real tnp1        = state[State_Type].curTime();
  const int  finest_level    = parent->finestLevel();
  NavierStokesBase::initial_iter = true;
  Real Sbar_old, Sbar_new;
  //
  // Make space to save a copy of the initial State_Type state data
  //
  Vector<std::unique_ptr<MultiFab> > saved_state(finest_level+1);

  if (init_iter > 0) {
      for (int k = 0; k <= finest_level; k++) {
          saved_state[k].reset(new MultiFab(getLevel(k).grids,
                                            getLevel(k).dmap,
                                            nState,nGrow));
      }
  }

  //
  // Iterate over the advance function.
  //
  for (int iter = 0; iter < init_iter; iter++)
  {
    //
    // Squirrel away copy of pre-advance State_Type state
    //
    for (int k = 0; k <= finest_level; k++)
      MultiFab::Copy(*saved_state[k],
                     getLevel(k).get_new_data(State_Type),
                     0,
                     0,
                     nState,
                     nGrow);

    for (int k = 0; k <= finest_level; k++ )
    {
      getLevel(k).advance(tnp1,dt_init,1,1);
    }
    //
    // This constructs a guess at P, also sets p_old == p_new.
    //
    Vector<MultiFab*> sig(finest_level+1,nullptr);
    for (int k = 0; k <= finest_level; k++)
    {
      sig[k] = &(getLevel(k).get_rho_half_time());
    }

    // ensure system is solvable by creating deltaS = S - Sbar
    if (closed_chamber == 1)
    {

      // compute number of cells
      Real num_cells = grids.numPts();

      // ensure divu_new is average down so computing deltaS = S - Sbar works for multilevel
      for (int k = finest_level-1; k >= 0; k--)
      {
        PeleLM&   fine_lev = getLevel(k+1);
        PeleLM&   crse_lev = getLevel(k);
	    
        MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

	crse_lev.average_down(Divu_fine, Divu_crse, 0, 1);
      }

      // compute Sbar and subtract from S
      for (int lev = 0; lev <= finest_level; lev++)
      {
        // pointer to S
        MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
        if (lev == 0)
        {
          Sbar_new = divu_lev.sum() / num_cells;
        }
        divu_lev.plus(-Sbar_new,0,1);
      }
	  
      // ensure divu_old is average down so computing deltaS = S - Sbar works for multilevel
      for (int k = finest_level-1; k >= 0; k--)
      {
        PeleLM&   fine_lev = getLevel(k+1);
        PeleLM&   crse_lev = getLevel(k);
		  
        MultiFab& Divu_crse = crse_lev.get_old_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_old_data(Divu_Type);

	crse_lev.average_down(Divu_fine, Divu_crse, 0, 1);
      }

      // compute Sbar and subtract from S
      for (int lev = 0; lev <= finest_level; lev++)
      {
        // pointer to S
        MultiFab& divu_lev = getLevel(lev).get_old_data(Divu_Type);
        if (lev == 0)
        {
          Sbar_old = divu_lev.sum() / num_cells;
        }
        divu_lev.plus(-Sbar_old,0,1);
      }
    }

    if (projector)
    {
      int havedivu = 1;
      projector->initialSyncProject(0,sig,parent->dtLevel(0),tnp1,
                                    havedivu);
    }
	
    if (closed_chamber == 1)
    {
      // restore S_new
      for (int lev = 0; lev <= finest_level; lev++)
      {
        MultiFab& divu_lev = getLevel(lev).get_new_data(Divu_Type);
        divu_lev.plus(Sbar_new,0,1);
      }

      // restore S_old
      for (int lev = 0; lev <= finest_level; lev++)
      {
        MultiFab& divu_lev = getLevel(lev).get_old_data(Divu_Type);
        divu_lev.plus(Sbar_old,0,1);
      }

    }

    for (int k = finest_level-1; k>= 0; k--)
    {
      getLevel(k).avgDown();
    }
    for (int k = 0; k <= finest_level; k++)
    {
      //
      // Reset state variables to initial time, but
      // do not pressure variable, only pressure time.
      //
      getLevel(k).resetState(tnp1, dt_init, dt_init);
    }
    //
    // For State_Type state, restore state we saved above
    //
    for (int k = 0; k <= finest_level; k++) {
      MultiFab::Copy(getLevel(k).get_new_data(State_Type),
                     *saved_state[k],
                     0,
                     0,
                     nState,
                     nGrow);
    }

    NavierStokesBase::initial_iter = false;
  }

  if (init_iter <= 0)
    NavierStokesBase::initial_iter = false; // Just being compulsive -- rbp.

  NavierStokesBase::initial_step = false;
  //
  // Re-instate timestep.
  //
  for (int k = 0; k <= finest_level; k++)
  {
    getLevel(k).setTimeLevel(tnp1,dt_save[k],dt_save[k]);
    if (getLevel(k).state[Press_Type].descriptor()->timeType() ==
        StateDescriptor::Point)
      getLevel(k).state[Press_Type].setNewTimeLevel(tnp1+.5*dt_init);
  }
  parent->setDtLevel(dt_save);
  parent->setNCycle(nc_save);
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure
// iteration section.
//
void
PeleLM::resetState (Real time,
                    Real dt_old,
                    Real dt_new)
{
   NavierStokesBase::resetState(time,dt_old,dt_new);

   state[RhoYdot_Type].reset();
   state[RhoYdot_Type].setTimeLevel(time,dt_old,dt_new);

   state[FuncCount_Type].reset();
   state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

void
PeleLM::avgDown ()
{
   if (level == parent->finestLevel()) return;

   const Real      strt_time = ParallelDescriptor::second();

   //
   // Average down the State and Pressure at the new time.
   //
   avgDown_StatePress();

   //
   // Reset the temperature
   //
   MultiFab&       S_crse    = get_new_data(State_Type);
   RhoH_to_Temp(S_crse);

   //
   // Next average down divu at new time.
   //
   PeleLM&         fine_lev  = getLevel(level+1);
   if (hack_noavgdivu) 
   {
     //
     // Now that state averaged down, recompute divu (don't avgDown,
     // since that will give a very different value, and screw up mac)
     //
     StateData& divuSD = get_state_data(Divu_Type);// should be const...
     const Real time   = divuSD.curTime();
     const Real dt     = time - divuSD.prevTime();
     calc_divu(time,dt,divuSD.newData());
   }
   else
   {
     MultiFab& Divu_crse = get_new_data(Divu_Type);
     MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

     average_down(Divu_fine, Divu_crse, 0, Divu_crse.nComp());
   }

   if (verbose > 1)
   {
     const int IOProc   = ParallelDescriptor::IOProcessorNumber();
     Real      run_time = ParallelDescriptor::second() - strt_time;

     ParallelDescriptor::ReduceRealMax(run_time,IOProc);

     amrex::Print() << "PeleLM::avgDown(): lev: " << level << ", time: " << run_time << '\n';
   }
}

static
Vector<const MultiFab *>
GetVecOfPtrs(const MultiFab* const* a, int scomp, int ncomp)
{
    Vector<const MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        r.push_back(new MultiFab(*a[d], amrex::make_alias, scomp, ncomp));
    }
    return r;
}

static
Vector<MultiFab *>
GetVecOfPtrs(MultiFab* const* a, int scomp, int ncomp)
{
    Vector<MultiFab*> r;
    r.reserve(AMREX_SPACEDIM);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        r.push_back(new MultiFab(*a[d], amrex::make_alias, scomp, ncomp));
    }
    return r;
}

void
PeleLM::diffusionFJDriver(ForkJoin&                   fj,
                  Real                        prev_time,
                  Real                        curr_time,
                  Real                        be_cn_theta,
                  int                         rho_flag,
                  const Vector<Real>&         visc_coef,
                  int                         visc_coef_comp,
                  const IntVect&              cratio,
                  const BCRec&                bc,
                  const Geometry&             in_geom,
                  bool                        add_hoop_stress,
                  const Diffusion::SolveMode& solve_mode,
                  bool                        add_old_time_divFlux,
                  const amrex::Vector<int>&   is_diffusive,
                  bool                        has_coarse_data,
                  bool                        has_delta_rhs,
                  bool                        has_alpha_in,
                  bool                        has_betan,
                  bool                        has_betanp1)
{
  int S_comp=0, Rho_comp=0, fluxComp=0, rhsComp=0, alpha_in_comp=0, betaComp=0;

  int nlev = (has_coarse_data ? 2 : 1);

  Vector<MultiFab*> S_old(nlev,0), S_new(nlev,0), Rho_old(nlev,0), Rho_new(nlev,0);

  S_old[0] = &(fj.get_mf("S_old_fine"));
  S_new[0] = &(fj.get_mf("S_new_fine"));
  if (nlev>1) {
    S_old[1] = &(fj.get_mf("S_old_crse"));
    S_new[1] = &(fj.get_mf("S_new_crse"));
  }
  int num_comp = S_old[0]->nComp();

  if (rho_flag == 2) {
    Rho_old[0] = &(fj.get_mf("Rho_old_fine"));
    Rho_new[0] = &(fj.get_mf("Rho_new_fine"));
    if (nlev>1) {
      Rho_old[1] = &(fj.get_mf("Rho_old_crse"));
      Rho_new[1] = &(fj.get_mf("Rho_new_crse"));
    }
  }

  const ForkJoin::ComponentSet compSet = fj.ComponentBounds("S_old_fine");
  int vStart = visc_coef_comp + compSet.lo;
  Vector<Real> visc_coef_shifted(&visc_coef[vStart],&visc_coef[vStart+num_comp]);

  MultiFab *rho_half_mf, *delta_rhs=0, *alpha_in=0;
  if (rho_flag == 1) {
    rho_half_mf = &(fj.get_mf("rho_half"));
  }

  Vector<MultiFab *> fluxn   = fj.get_mf_vec("fluxn");
  Vector<MultiFab *> fluxnp1 = fj.get_mf_vec("fluxnp1");

  if (has_delta_rhs) {
    delta_rhs = &(fj.get_mf("delta_rhs"));
  }

  if (has_alpha_in) {
    alpha_in = &(fj.get_mf("alpha_in"));
  }

  Vector<MultiFab *> betan(AMREX_SPACEDIM,0), betanp1(AMREX_SPACEDIM,0);
  if (has_betan) {
    betan   = fj.get_mf_vec("betan");  
  }
  if (has_betanp1) {
    betanp1 = fj.get_mf_vec("betanp1");
  }

  Vector<MultiFab *> area_mf = fj.get_mf_vec("area");

  diffusion->diffuse_scalar (S_old,Rho_old,S_new,Rho_new,S_comp,num_comp,Rho_comp,
                             prev_time,curr_time,be_cn_theta,*rho_half_mf,rho_flag,
                             &(fluxn[0]),&(fluxnp1[0]),fluxComp,delta_rhs,rhsComp,
                             alpha_in,alpha_in_comp,&(betan[0]),&(betanp1[0]),betaComp,
                             cratio,bc,in_geom,
                             solve_mode,add_old_time_divFlux,is_diffusive);
}

void
PeleLM::diffuse_scalar_fj  (const Vector<MultiFab*>&  S_old,
                            const Vector<MultiFab*>&  Rho_old,
                            Vector<MultiFab*>&        S_new,
                            const Vector<MultiFab*>&  Rho_new,
                            int                       S_comp,
                            int                       num_comp,
                            int                       Rho_comp,
                            Real                      prev_time,
                            Real                      curr_time,
                            Real                      be_cn_theta,
                            const MultiFab&           rho_mid,
                            int                       rho_flag,
                            MultiFab* const*          fluxn,
                            MultiFab* const*          fluxnp1,
                            int                       fluxComp,
                            MultiFab*                 delta_rhs, 
                            int                       rhsComp,
                            const MultiFab*           alpha_in, 
                            int                       alpha_in_comp,
                            const MultiFab* const*    betan, 
                            const MultiFab* const*    betanp1,
                            int                       betaComp,
                            const Vector<Real>&       visc_coef,
                            int                       visc_coef_comp,
                            const MultiFab&           theVolume,
                            const MultiFab* const*    theArea,
                            const IntVect&            cratio,
                            const BCRec&              bc,
                            const Geometry&           theGeom,
                            bool                      add_hoop_stress,
                            const Diffusion::SolveMode& solve_mode,
                            bool                      add_old_time_divFlux,
                            const amrex::Vector<int>& diffuse_this_comp)
{
  int n_procs = ParallelDescriptor::NProcs();
  int n_tasks_suggest = std::min(num_comp,n_procs);
  int n_tasks = std::min(std::max(1,num_forkjoin_tasks),n_tasks_suggest);

  if (n_tasks == 1)
  {
    diffusion->diffuse_scalar(S_old,Rho_old,S_new,Rho_new,S_comp,num_comp,Rho_comp,
                                  prev_time,curr_time,be_cn_theta,rho_mid,rho_flag,
                                  fluxn,fluxnp1,fluxComp,delta_rhs,rhsComp,
                                  alpha_in,alpha_in_comp,betan,betanp1,betaComp,
                                  cratio,bc,theGeom,
                                  solve_mode,add_old_time_divFlux,diffuse_this_comp);
  }
  else
  {

    Print() << "Diffusion: using " << n_tasks << " fork-join tasks for "
            << num_comp << " diffusion calls (on a total of " << n_procs << " ranks)" << std::endl;

    ForkJoin fj(n_tasks);
    fj.SetVerbose(forkjoin_verbose);

    MultiFab S_old_fine(*S_old[0], amrex::make_alias, S_comp, num_comp);
    MultiFab S_new_fine(*S_new[0], amrex::make_alias, S_comp, num_comp);
    MultiFab Rho_old_fine(*Rho_old[0], amrex::make_alias, Rho_comp, 1);
    MultiFab Rho_new_fine(*Rho_new[0], amrex::make_alias, Rho_comp, 1);
    MultiFab *S_old_crse, *S_new_crse, *Rho_old_crse, *Rho_new_crse, *Alpha, *Rhs;

    const int ng = 1;
    AMREX_ALWAYS_ASSERT(S_old_fine.nGrow() >= 1 && S_new_fine.nGrow() >= 1);
    fj.reg_mf(S_old_fine,  "S_old_fine",  ForkJoin::Strategy::split,    ForkJoin::Intent::in,    ng);
    fj.reg_mf(S_new_fine,  "S_new_fine",  ForkJoin::Strategy::split,    ForkJoin::Intent::inout, ng);
    if (rho_flag == 2) {
      AMREX_ALWAYS_ASSERT(Rho_old_fine.nGrow() >= 1 && Rho_new_fine.nGrow() >= 1);
      fj.reg_mf(Rho_old_fine,"Rho_old_fine",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in, ng);
      fj.reg_mf(Rho_new_fine,"Rho_new_fine",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in, ng);
    }

    bool has_coarse_data = S_old.size() > 1;
    if (has_coarse_data) {
      S_old_crse = new MultiFab(*S_old[1], amrex::make_alias, S_comp, num_comp);
      S_new_crse = new MultiFab(*S_new[1], amrex::make_alias, S_comp, num_comp);
      AMREX_ALWAYS_ASSERT(S_old_crse->nGrow() >= 1 && S_new_crse->nGrow() >= 1);
      fj.reg_mf(*S_old_crse,  "S_old_crse",  ForkJoin::Strategy::split,    ForkJoin::Intent::in, ng);
      fj.reg_mf(*S_new_crse,  "S_new_crse",  ForkJoin::Strategy::split,    ForkJoin::Intent::in, ng);
      if (rho_flag == 2) {
        Rho_old_crse = new MultiFab(*Rho_old[1], amrex::make_alias, Rho_comp, 1);
        Rho_new_crse = new MultiFab(*Rho_new[1], amrex::make_alias, Rho_comp, 1);
        AMREX_ALWAYS_ASSERT(Rho_old_crse->nGrow() >= 1 && Rho_new_crse->nGrow() >= 1);
        fj.reg_mf(*Rho_old_crse,"Rho_old_crse",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in, ng);
        fj.reg_mf(*Rho_new_crse,"Rho_new_crse",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in, ng);
      }
    }

    if (rho_flag == 1) {
      fj.reg_mf(rho_mid,"rho_half",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in);
    }

    fj.reg_mf_vec(GetVecOfPtrs(fluxn  ,fluxComp,num_comp),"fluxn",  ForkJoin::Strategy::split,ForkJoin::Intent::inout);
    fj.reg_mf_vec(GetVecOfPtrs(fluxnp1,fluxComp,num_comp),"fluxnp1",ForkJoin::Strategy::split,ForkJoin::Intent::inout);

    bool has_delta_rhs = false;
    if (delta_rhs != 0) {
      AMREX_ALWAYS_ASSERT(delta_rhs->nComp() >= rhsComp + num_comp);
      Rhs = new MultiFab(*delta_rhs, amrex::make_alias, rhsComp, num_comp);
      fj.reg_mf(*Rhs,"delta_rhs",ForkJoin::Strategy::split,ForkJoin::Intent::in);
      has_delta_rhs = true;
    }

    bool has_alpha_in = false;
    if (alpha_in != 0) {
      AMREX_ALWAYS_ASSERT(alpha_in->nComp() >= alpha_in_comp + num_comp);
      Alpha = new MultiFab(*alpha_in, amrex::make_alias, alpha_in_comp, num_comp);
      fj.reg_mf(*Alpha,"alpha_in",ForkJoin::Strategy::split,ForkJoin::Intent::in);
      has_alpha_in = true;
    }

    int allnull, allthere;
    Diffusion::checkBeta(betan,   allthere, allnull);
    bool has_betan = false;
    if (allthere) {
      fj.reg_mf_vec(GetVecOfPtrs(betan,  betaComp,num_comp),"betan",  ForkJoin::Strategy::split,ForkJoin::Intent::in);
      has_betan = true;
    }

    Diffusion::checkBeta(betanp1, allthere, allnull);
    bool has_betanp1 = false;
    if (allthere) {
      fj.reg_mf_vec(GetVecOfPtrs(betanp1,betaComp,num_comp),"betanp1",ForkJoin::Strategy::split,ForkJoin::Intent::in);
      has_betanp1 = true;
    }

    fj.reg_mf(theVolume,"volume",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in);
    fj.reg_mf_vec(GetVecOfPtrs(theArea,0,1),"area",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in);

    fj.fork_join(
      [=,&bc,&theGeom,&visc_coef] (ForkJoin &f)
      {
        diffusionFJDriver(f,
                          prev_time,
                          curr_time,
                          be_cn_theta,
                          rho_flag,
                          visc_coef,
                          visc_coef_comp,
                          cratio,
                          bc,
                          theGeom,
                          add_hoop_stress,
                          solve_mode,
                          add_old_time_divFlux,
                          diffuse_this_comp,
                          has_coarse_data, has_delta_rhs, has_alpha_in, has_betan, has_betanp1);
      }
      );
  }
}

void
PeleLM::differential_diffusion_update (MultiFab& Force,
                                       MultiFab& Dwbar,
                                       int       FComp,
                                       MultiFab& Dnew,
                                       int       DComp,
                                       MultiFab& DDnew)
{
   BL_PROFILE_REGION_START("R::PLM::differential_diffusion_update()");
   BL_PROFILE("PLM::differential_diffusion_update()");

   // Recompute the D[Nspec+1] = Div(Fi.Hi) using
   // Fi.Hi based on solution here.
   //
   AMREX_ASSERT(Force.boxArray() == grids);
   AMREX_ASSERT(FComp+Force.nComp()>=NUM_SPECIES+1);
   AMREX_ASSERT(Dnew.boxArray() == grids);
   AMREX_ASSERT(DDnew.boxArray() == grids);
   AMREX_ASSERT(DComp+Dnew.nComp()>=NUM_SPECIES+2);

   const Real strt_time = ParallelDescriptor::second();

   if (hack_nospecdiff)
   {
      amrex::Error("differential_diffusion_update: hack_nospecdiff not implemented");
   }

   MultiFab Rh; // allocated memeory not needed for this, since rho_flag=2 for Y solves

   int nGrow=1;
   const Real prev_time = state[State_Type].prevTime();
   const Real curr_time = state[State_Type].curTime();
   const Real dt = curr_time - prev_time;

   int sComp = std::min((int)Density, std::min((int)first_spec,(int)Temp) );
   int eComp = std::max((int)Density, std::max((int)last_spec,(int)Temp) );
   int nComp = eComp - sComp + 1;

   // Set new data to old on valid, but FillPatch old and new to get Dirichlet boundary data for each
   MultiFab::Copy(get_new_data(State_Type),get_old_data(State_Type),first_spec,first_spec,NUM_SPECIES,0);

   FillPatch(*this,get_old_data(State_Type),nGrow,prev_time,State_Type,sComp,nComp,sComp);
   FillPatch(*this,get_new_data(State_Type),nGrow,curr_time,State_Type,sComp,nComp,sComp);

   auto Snc = std::unique_ptr<MultiFab>(new MultiFab());
   auto Snp1c = std::unique_ptr<MultiFab>(new MultiFab());

   if (level > 0) {
      auto& crselev = getLevel(level-1);
      Snc->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, nGrow);  Snc->setVal(0,0,NUM_STATE,nGrow);
      FillPatch(crselev,*Snc  ,nGrow,prev_time,State_Type,sComp,nComp,sComp);
      Snp1c->define(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, nGrow);
      FillPatch(crselev,*Snp1c,nGrow,curr_time,State_Type,sComp,nComp,sComp);
   }

   const int nlev = (level ==0 ? 1 : 2);
   Vector<MultiFab*> Sn(nlev,0), Snp1(nlev,0);
   Sn[0]   = &(get_old_data(State_Type));
   Snp1[0] = &(get_new_data(State_Type));

   if (nlev>1) {
      Sn[1]   =  Snc.get();
      Snp1[1] =  Snp1c.get();
   }

   const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();

   MultiFab *delta_rhs = &Force;
   const int rhsComp = FComp;

   const MultiFab *alpha = 0;
   const int alphaComp = 0, fluxComp = 0;

   FluxBoxes fb_dnp1;
   MultiFab **betan = 0; // not needed
   MultiFab **betanp1 = fb_dnp1.define(this,NUM_SPECIES+1);
   getDiffusivity(betanp1, curr_time, first_spec, 0, NUM_SPECIES); // species (rhoD)
   getDiffusivity(betanp1, curr_time, Temp, NUM_SPECIES, 1); // temperature (lambda)

   Vector<int> diffuse_comp(NUM_SPECIES+1);
   for (int icomp=0; icomp<NUM_SPECIES+1; ++icomp) {
     diffuse_comp[icomp] = is_diffusive[first_spec + icomp];
   }

   const int rho_flag = Diffusion::set_rho_flag(diffusionType[first_spec]);
   const BCRec& bc = theBCs[first_spec];
   for (int icomp=1; icomp<NUM_SPECIES; ++icomp) {
     AMREX_ALWAYS_ASSERT(rho_flag == Diffusion::set_rho_flag(diffusionType[first_spec+icomp]));
     AMREX_ALWAYS_ASSERT(bc == theBCs[first_spec+icomp]);
   }

   const bool add_hoop_stress = false; // Only true if sigma == Xvel && Geometry::IsRZ())
   const Diffusion::SolveMode& solve_mode = Diffusion::ONEPASS;
   const bool add_old_time_divFlux = false; // rhs contains the time-explicit diff terms already
   const Real be_cn_theta_SDC = 1;

   const int betaComp = 0;
   const int visc_coef_comp = first_spec;
   const int Rho_comp = Density;
   const int bc_comp = first_spec;

   const MultiFab *a[AMREX_SPACEDIM];
   for (int d=0; d<AMREX_SPACEDIM; ++d) {
       a[d] = &(area[d]);
   }

   // Diffuse all the species
   diffuse_scalar_fj(Sn, Sn, Snp1, Snp1, first_spec, NUM_SPECIES, Rho_comp,
                     prev_time,curr_time,be_cn_theta_SDC,Rh,rho_flag,
                     SpecDiffusionFluxn,SpecDiffusionFluxnp1,fluxComp,
                     delta_rhs,rhsComp,alpha,alphaComp,
                     betan,betanp1,betaComp,visc_coef,visc_coef_comp,
                     volume,a,crse_ratio,theBCs[bc_comp],geom,
                     add_hoop_stress,solve_mode,add_old_time_divFlux,diffuse_comp);

// Here we apply to covered cells the saved reference state data 
// in order to avoid non-physical values such as Ys=0 for all species
#ifdef AMREX_USE_EB
   set_body_state(*Snp1[0]);
   set_body_state(*Sn[0]);
#endif

   if (use_wbar) {
      // add lagged grad Wbar fluxes (SpecDiffusionFluxWbar) to time-advanced 
      // species diffusion fluxes (SpecDiffusionFluxnp1)(kp1) at this point.
      // These fluxes are similar to the one we put into the RHS in advance() since they are not
      // solved implicitly.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*Sn[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         for (int dir=0; dir<AMREX_SPACEDIM; ++dir)
         {
            const Box& ebx = mfi.nodaltilebox(dir);
            auto const& flux_spec = SpecDiffusionFluxnp1[dir]->array(mfi);
            auto const& flux_wbar = SpecDiffusionFluxWbar[dir]->array(mfi);
            amrex::ParallelFor(ebx, NUM_SPECIES, [ flux_spec, flux_wbar ]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                flux_spec(i,j,k,n) += flux_wbar(i,j,k,n);
            });
         }
      }
   }

   //
   // Modify/update new-time fluxes to ensure sum of species fluxes = 0, then compute Dnew[m] = -Div(flux[m])
   //
   const BCRec& Tbc = AmrLevel::desc_lst[State_Type].getBCs()[Temp];

   adjust_spec_diffusion_fluxes(SpecDiffusionFluxnp1,get_new_data(State_Type),
                                Tbc,curr_time);

   // Compute Dnp1kp1 (here called Dnew)
   flux_divergence(Dnew,DComp,SpecDiffusionFluxnp1,0,NUM_SPECIES,-1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   // Update species after diffusion solve using input Force 
   for (MFIter mfi(*Snp1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& rhoY_old = Sn[0]->array(mfi,first_spec);
      auto const& rhoY_new = Snp1[0]->array(mfi,first_spec);
      auto const& force    = Force.array(mfi);
      auto const& dnp1kp1  = Dnew.array(mfi);
      int use_wbar_lcl     = use_wbar;
      auto const& dwbar    = (use_wbar) ? Dwbar.array(mfi) : Dnew.array(mfi); 
      amrex::ParallelFor(bx, NUM_SPECIES, [ rhoY_old, rhoY_new, force, dnp1kp1, dwbar, dt, use_wbar_lcl]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
          rhoY_new(i,j,k,n) = rhoY_old(i,j,k,n) + dt * ( force(i,j,k,n) + dnp1kp1(i,j,k,n) );
          if ( use_wbar_lcl ) {  // We remove the dwbar term since it was included in both force (in PeleLM::advance) and in the fluxes (above)
             rhoY_new(i,j,k,n) -= dt * dwbar(i,j,k,n);
          }
      });
   }

   //
   // Do iterative enthalpy/temperature solve
   //

   // build energy fluxes based on species fluxes, Gamma_m, and cell-centered states
   // 1. flux[NUM_SPECIES+1] = sum_m (H_m Gamma_m)
   // 2. compute flux[NUM_SPECIES+2] = - lambda grad T
   //
   compute_enthalpy_fluxes(SpecDiffusionFluxnp1,betanp1,curr_time);

   // Divergence of energy fluxes:
   // 1. Dnew[N+1] = -Div(flux[N+2])
   // 2. DD = -Sum{ Div(H_m Gamma_m) }
   flux_divergence(Dnew,DComp+NUM_SPECIES+1,SpecDiffusionFluxnp1,NUM_SPECIES+2,1,-1);

   flux_divergence(DDnew,0,SpecDiffusionFluxnp1,NUM_SPECIES+1,1,-1);

   if (deltaT_verbose) {
      Print() << "Iterative solve for deltaT: " << std::endl;
   }
   MultiFab Trhs(grids,dmap,1,0,MFInfo(),Factory());
   MultiFab Told(grids,dmap,1,1,MFInfo(),Factory());
   MultiFab RhoCp(grids,dmap,1,0,MFInfo(),Factory());
   MultiFab RhT(get_new_data(State_Type), amrex::make_alias, Density, 1);

   Real deltaT_iter_norm = 0;
   for (int L=0; L<num_deltaT_iters_MAX && (L==0 || deltaT_iter_norm >= deltaT_norm_max); ++L)
   {
      // Bundle in a single kernel assembling the Trhs and compute rhoCpmix
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*Snp1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         // Things needed for the Trhs
         auto const& rhs      = Trhs.array(mfi); 
         auto const& rhoH_o   = Sn[0]->array(mfi,RhoH);
         auto const& rhoH_n   = Snp1[0]->array(mfi,RhoH);
         auto const& force    = Force.array(mfi,NUM_SPECIES); 
         auto const& dnp1kp1  = Dnew.array(mfi,DComp+NUM_SPECIES+1);
         auto const& ddnp1kp1 = DDnew.array(mfi);
         Real        dtinv    = 1.0/dt;

         // Things needed to compute RhoCpmix
         auto const& rho     = Snp1[0]->array(mfi,Density);
         auto const& rhoY    = Snp1[0]->array(mfi,first_spec);
         auto const& T       = Snp1[0]->array(mfi,Temp);
         auto const& RhoCpm  = RhoCp.array(mfi);

         // Things needed to store T
         auto const& T_old   = Told.array(mfi);

         amrex::ParallelFor(bx, [rhs, rhoH_o, rhoH_n, force, dnp1kp1, ddnp1kp1, dtinv,
                                 rho, rhoY, T, RhoCpm, T_old]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            // Trhs computation
            rhs(i,j,k) =  ( rhoH_o(i,j,k) - rhoH_n(i,j,k) ) * dtinv +
                        + force(i,j,k) + dnp1kp1(i,j,k) + ddnp1kp1(i,j,k);
            // rhoCpmix computation
            getCpmixGivenRYT( i, j, k, rho, rhoY, T, RhoCpm );
            RhoCpm(i,j,k) *= rho(i,j,k);

            // Store Told
            T_old(i,j,k) = T(i,j,k);
         });
      }

      Snp1[0]->setVal(0.0,Temp,1,1);
      if (nlev>1) {
         Snp1[1]->setVal(0.0,Temp,1,1);
      }

      int rho_flagT = 0; // Do not do rho-based hacking of the diffusion problem
      const Vector<int> diffuse_this_comp = {1};
      diffusion->diffuse_scalar(Sn, Sn, Snp1, Snp1, Temp, 1, Rho_comp,
                                prev_time,curr_time,be_cn_theta_SDC,RhT,rho_flagT,
                                SpecDiffusionFluxn,SpecDiffusionFluxnp1,NUM_SPECIES+2,
                                &Trhs,0,&RhoCp,0,
                                betan,betanp1,NUM_SPECIES,
                                crse_ratio,theBCs[Temp],geom,
                                solve_mode,add_old_time_divFlux,diffuse_this_comp);

#ifdef AMREX_USE_EB
      EB_set_covered_faces({AMREX_D_DECL(SpecDiffusionFluxnp1[0],SpecDiffusionFluxnp1[1],SpecDiffusionFluxnp1[2])},0.);
      EB_set_covered_faces({AMREX_D_DECL(SpecDiffusionFluxn[0],SpecDiffusionFluxn[1],SpecDiffusionFluxn[2])},0.);
#endif
          
      deltaT_iter_norm = Snp1[0]->norm0(Temp);
      if (deltaT_verbose) {
        Print() << "   DeltaT solve norm [" << L << "] = " << deltaT_iter_norm << std::endl;
      }

      // T <= T + deltaT
      MultiFab::Add(get_new_data(State_Type),Told,0,Temp,1,0);

      // Update energy fluxes, divergences
      compute_enthalpy_fluxes(SpecDiffusionFluxnp1,betanp1,curr_time);
      flux_divergence(Dnew,DComp+NUM_SPECIES+1,SpecDiffusionFluxnp1,NUM_SPECIES+2,1,-1);
      flux_divergence(DDnew,0,SpecDiffusionFluxnp1,NUM_SPECIES+1,1,-1);
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      // Update (RhoH)^{k+1,L+1}
      for (MFIter mfi(*Snp1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = Snp1[0]->array(mfi,Density);  
         auto const& rhoY    = Snp1[0]->array(mfi,first_spec);
         auto const& T       = Snp1[0]->array(mfi,Temp);
         auto const& rhoHm   = Snp1[0]->array(mfi,RhoH);

         amrex::ParallelFor(bx, [rho, rhoY, T, rhoHm]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
         });
      }

      if ( L==(num_deltaT_iters_MAX-1) 
           && deltaT_iter_norm >= deltaT_norm_max
           && deltaT_crashOnConvFail ) {
         Abort("deltaT_iters not converged");
      }
   } // end deltaT iters

#ifdef AMREX_USE_EB
   set_body_state(*Snp1[0]);
   EB_set_covered_faces({AMREX_D_DECL(SpecDiffusionFluxnp1[0],SpecDiffusionFluxnp1[1],SpecDiffusionFluxnp1[2])},0.);
#endif

   showMF("mysdc",Dnew,"sdc_Dnew_afterDeltaTiter_inDDupdate",level,1,parent->levelSteps(level));

   // We have just performed the correction diffusion solve for Y_m and h
   // If updateFluxReg=T, we update VISCOUS flux registers:
   //   ADD Gamma_{m,AD}^(k+1),           (in flux[0:NUM_SPECIES-1])
   //   ADD -lambda^(k) grad T_AD^(k+1)   (in flux[NUM_SPECIES+2])
   // 
   // And update the ADVECTIVE flux registers:
   //   ADD h_m . Gamma_{m,AD}            (in flux[NUM_SPECIES+1])
   //
   if (do_reflux && updateFluxReg)
   {
     showMF("DBGSync",*SpecDiffusionFluxnp1[0],"DBGSync_DiffFluxX_Dhat",level,parent->levelSteps(level));
     showMF("DBGSync",*SpecDiffusionFluxnp1[1],"DBGSync_DiffFluxY_Dhat",level,parent->levelSteps(level));
     for (int d = 0; d < AMREX_SPACEDIM; d++)
     {
       if (level > 0)
       {
         auto& vfr = getViscFluxReg();
         auto& afr = getAdvFluxReg();
         vfr.FineAdd(*SpecDiffusionFluxnp1[d],d,0,first_spec,NUM_SPECIES,dt);
         afr.FineAdd(*SpecDiffusionFluxnp1[d],d,NUM_SPECIES+1,RhoH,1,dt);
         vfr.FineAdd(*SpecDiffusionFluxnp1[d],d,NUM_SPECIES+2,RhoH,1,dt);
       }
       if (level < parent->finestLevel())
       {
         auto& vfr = getLevel(level+1).getViscFluxReg();
         auto& afr = getLevel(level+1).getAdvFluxReg();
         vfr.CrseInit((*SpecDiffusionFluxnp1[d]),d,0,first_spec,NUM_SPECIES,-dt,FluxRegister::ADD);
         afr.CrseInit((*SpecDiffusionFluxnp1[d]),d,NUM_SPECIES+1,RhoH,1,-dt,FluxRegister::ADD);
         vfr.CrseInit((*SpecDiffusionFluxnp1[d]),d,NUM_SPECIES+2,RhoH,1,-dt,FluxRegister::ADD);
       }
     }
   }

   //
   // Ensure consistent grow cells
   //
   if (Dnew.nGrow() > 0)
   {
      Dnew.FillBoundary(DComp, NUM_SPECIES+1, geom.periodicity());
      Extrapolater::FirstOrderExtrap(Dnew, geom, DComp, NUM_SPECIES+1);
   }

   if (DDnew.nGrow() > 0)
   {
      DDnew.FillBoundary(0, 1, geom.periodicity());
      Extrapolater::FirstOrderExtrap(DDnew, geom, 0, 1);
   }

   if (verbose)
   {
     const int IOProc   = ParallelDescriptor::IOProcessorNumber();
     Real      run_time = ParallelDescriptor::second() - strt_time;

     ParallelDescriptor::ReduceRealMax(run_time,IOProc);

     amrex::Print() << "PeleLM::differential_diffusion_update(): lev: " << level
                    << ", time: " << run_time << '\n';
   }

   BL_PROFILE_REGION_STOP("R::PLM::differential_diffusion_update()");
}

void
PeleLM::adjust_spec_diffusion_fluxes (MultiFab* const * flux,
                                      const MultiFab&   S,
                                      const BCRec&      bc,
                                      Real              time)
{
   BL_PROFILE("PeleLM:::adjust_spec_diffusion_fluxes()");
   //
   // Adjust the species diffusion fluxes so that their sum is zero.
   //
   const Real strt_time = ParallelDescriptor::second();
   const Box& domain = geom.Domain();

   int ngrow = 3;
   MultiFab TT(grids,dmap,NUM_SPECIES,ngrow,MFInfo(),Factory());
   FillPatch(*this,TT,ngrow,time,State_Type,first_spec,NUM_SPECIES,0);

#ifdef AMREX_USE_EB
   Vector<BCRec> math_bc(NUM_SPECIES);
   math_bc = fetchBCArray(State_Type,first_spec,NUM_SPECIES);

   MultiFab edgstate[AMREX_SPACEDIM];
   int nghost(4);         // Use 4 for now

   for (int i(0); i < AMREX_SPACEDIM; i++)
   {
     const BoxArray& ba = getEdgeBoxArray(i);
     edgstate[i].define(ba, dmap, NUM_SPECIES, nghost, MFInfo(), Factory());
   }
   EB_interp_CellCentroid_to_FaceCentroid(TT, AMREX_D_DECL(edgstate[0],edgstate[1],edgstate[2]), 0, 0, NUM_SPECIES, geom, math_bc);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {    
      for (int dir =0; dir < AMREX_SPACEDIM; ++dir)
      {
         const Box& ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         auto const& rhoY     = TT.array(mfi);
         auto const& flux_dir = flux[dir]->array(mfi);

#ifdef AMREX_USE_EB
         const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(S[mfi]);
         const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();

         if (flags.getType(amrex::grow(ebx,0)) != FabType::covered )
         {
            // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
            if (flags.getType(amrex::grow(ebx,nghost)) == FabType::regular )
            {
               amrex::ParallelFor(ebx, [dir, rhoY, flux_dir, edomain, bc]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  int idx[3] = {i,j,k};
                  bool on_lo = ( ( bc.lo(dir) == EXT_DIR ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
                  bool on_hi = ( ( bc.hi(dir) == EXT_DIR ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
                  repair_flux( i, j, k, dir, on_lo, on_hi, rhoY, flux_dir );
               });
            }
            else
            {         
               auto const& rhoYed_d   = edgstate[dir].array(mfi);
               auto const& areafrac_d = areafrac[dir]->array(mfi);
               amrex::ParallelFor(ebx, [dir, rhoY, flux_dir, rhoYed_d, areafrac_d, edomain, bc]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  int idx[3] = {i,j,k};
                  bool on_lo = ( ( bc.lo(dir) == EXT_DIR ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
                  bool on_hi = ( ( bc.hi(dir) == EXT_DIR ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
                  repair_flux_eb( i, j, k, dir, on_lo, on_hi, rhoY, rhoYed_d, areafrac_d, flux_dir );
               });
            }
         }
#else
         amrex::ParallelFor(ebx, [dir, rhoY, flux_dir, edomain, bc]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc.lo(dir) == EXT_DIR ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc.hi(dir) == EXT_DIR ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            repair_flux( i, j, k, dir, on_lo, on_hi, rhoY, flux_dir );
         });

#endif
      }
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::adjust_spec_diffusion_fluxes(): lev: " << level 
                     << ", time: " << run_time << '\n';
   }
}

void
PeleLM::compute_enthalpy_fluxes (MultiFab* const*       flux,
                                 const MultiFab* const* beta,
                                 Real                   time)
{
  BL_PROFILE("PeleLM:::compute_enthalpy_fluxes()");
  /*
    Build heat fluxes based on species fluxes, Gamma_m, and fill-patched cell-centered states
    Set:
         flux[NUM_SPECIES+1] = sum_m (H_m Gamma_m)
         flux[NUM_SPECIES+2] = - lambda grad T
  */

  AMREX_ASSERT(beta && beta[0]->nComp() == NUM_SPECIES+1);

  const Real strt_time = ParallelDescriptor::second();

  //
  // First step, we create an operator to get (EB-aware) fluxes from, it will provide flux[NUM_SPECIES+2]
  //
  
  LPInfo info;
  info.setAgglomeration(1);
  info.setConsolidation(1);
  info.setMetricTerm(false);
  info.setMaxCoarseningLevel(0);
  
#ifdef AMREX_USE_EB
  const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>((parent->getLevel(level)).Factory());
  MLEBABecLap op({geom}, {grids}, {dmap}, info, {ebf});
#else
  MLABecLaplacian op({geom}, {grids}, {dmap}, info);
#endif

  op.setMaxOrder(diffusion->maxOrder());
  MLMG mg(op);

  {
    const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();
    const BCRec& bc = theBCs[Temp];

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    Diffusion::setDomainBC(mlmg_lobc, mlmg_hibc, bc); // Same for all comps, by assumption
    op.setDomainBC(mlmg_lobc, mlmg_hibc);
  }

  MultiFab TTc;
  if (level > 0)
  {
    PeleLM& clev = getLevel(level-1);
    TTc.define(clev.grids,clev.dmap,1,0,MFInfo(),clev.Factory());
    FillPatch(clev,TTc,0,time,State_Type,Temp,1,0);
    op.setCoarseFineBC(&TTc, crse_ratio[0]);
  }
  int ngrow = 1;
  MultiFab TT(grids,dmap,1,ngrow,MFInfo(),Factory());
  FillPatch(*this,TT,ngrow,time,State_Type,Temp,1,0);
  op.setLevelBC(0, &TT);

  // Creating alpha and beta coefficients.
  Real      a               = 0;
  Real      b               = 1.;
  op.setScalars(a, b);


  // Here it is NUM_SPECIES because lambda is stored after the last species (first starts at 0)
  Diffusion::setBeta(op,beta,NUM_SPECIES);
  
  AMREX_D_TERM( flux[0]->setVal(0., NUM_SPECIES+2, 1);,
                flux[1]->setVal(0., NUM_SPECIES+2, 1);,
                flux[2]->setVal(0., NUM_SPECIES+2, 1););

  // Here it is NUM_SPECIES+2 because this is the heat flux (NUM_SPECIES+3 in enth_diff_terms in fortran)
  // No multiplication by dt here
  Diffusion::computeExtensiveFluxes(mg, TT, flux, NUM_SPECIES+2, 1, &geom, b);

  //
  // Now we have flux[NUM_SPECIES+2]
  // Second step, let's compute flux[NUM_SPECIES+1] = sum_m (H_m Gamma_m)
  

   // First we want to gather h_i from T and then interpolate it to face centroids
   MultiFab Enth(grids,dmap,NUM_SPECIES,ngrow,MFInfo(),Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(TT,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx  = mfi.tilebox();
         const Box gbx  = grow(bx,1);
         auto const& T  = TT.array(mfi);
         auto const& Hi = Enth.array(mfi,0);

         amrex::ParallelFor(gbx, [T, Hi]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getHGivenT( i, j, k, T, Hi );
         });
      }
   }

  Vector<BCRec> math_bc(NUM_SPECIES);
  math_bc = fetchBCArray(State_Type,first_spec,NUM_SPECIES);

  MultiFab enth_edgstate[AMREX_SPACEDIM];

  for (int i(0); i < AMREX_SPACEDIM; i++)
  {
    const BoxArray& ba = getEdgeBoxArray(i);
    enth_edgstate[i].define(ba, dmap, NUM_SPECIES, 1, MFInfo(), Factory());
  }

#ifdef AMREX_USE_EB
   EB_interp_CellCentroid_to_FaceCentroid(Enth, AMREX_D_DECL(enth_edgstate[0],enth_edgstate[1],enth_edgstate[2]), 0, 0, NUM_SPECIES, geom, math_bc);

#else
   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(Enth,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      const Box& bx  = mfi.tilebox();
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         //const Box ebx = surroundingNodes(bx,dir);
         const Box ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         const auto& enth_c  = Enth.array(mfi,0);
         const auto& enth_ed = enth_edgstate[dir].array(mfi,0);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, use_harmonic_avg, enth_c, enth_ed, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, NUM_SPECIES, use_harmonic_avg, on_lo, on_hi, enth_c, enth_ed);
         });
      }
   }
#endif

  //
  // Now we construct the actual fluxes: sum[ (species flux).(species enthalpy) ]
  //
  
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
  FArrayBox fab_tmp;
  for (MFIter mfi(TT,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();

#if AMREX_USE_EB
    // this is to check efficiently if this tile contains any eb stuff
    const EBFArrayBox& in_fab = static_cast<EBFArrayBox const&>(TT[mfi]);
    const EBCellFlagFab& flags = in_fab.getEBCellFlagFab();

    if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
    {
      // If tile is completely covered by EB geometry, set 
      // value to some very large number so we know if
      // we accidentaly use these covered vals later in calculations
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      {
        const Box&  ebox_d          = mfi.nodaltilebox(dir);
        auto const& flux_enthaply_d = (*flux[dir]).array(mfi, NUM_SPECIES+1);

        amrex::ParallelFor(ebox_d, [=]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              flux_enthaply_d(i,j,k) = 1.2345e30;
          });
      }
    }
    else
    {
#endif
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      {
        const Box&  ebox_d          = mfi.nodaltilebox(dir);
        auto const& flux_species_d  = (*flux[dir]).array(mfi, 0);
        auto const& flux_enthaply_d = (*flux[dir]).array(mfi, NUM_SPECIES+1);
        auto const& enth_d          = enth_edgstate[dir].array(mfi);

        amrex::ParallelFor(ebox_d, [=]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              flux_enthaply_d(i,j,k) = 0.0;
              for (int n = 0; n < NUM_SPECIES; n++) {
                  flux_enthaply_d(i,j,k) += flux_species_d(i,j,k,n)*enth_d(i,j,k,n);
              }
          });
      }
#if AMREX_USE_EB
    }
#endif
  }
}

  if (verbose > 1)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(run_time,IOProc);
    amrex::Print() << "PeleLM::compute_enthalpy_fluxes(): lev: " << level 
                   << ", time: " << run_time << '\n';
  }
}
    
void
PeleLM::velocity_diffusion_update (Real dt)
{
   BL_PROFILE("PeleLM::velocity_diffusion_update()");
   //
   // Do implicit c-n solve for velocity
   // compute the viscous forcing
   // do following except at initial iteration--rbp, per jbb
   //
   if (is_diffusive[Xvel])
   {
      const Real strt_time = ParallelDescriptor::second();

      const Real time = state[State_Type].prevTime();

      int rho_flag;
      if (do_mom_diff == 0)
      {
         rho_flag = 1;
      }
      else
      {
         rho_flag = 3;
      }

      MultiFab *delta_rhs = 0;

      FluxBoxes fb_betan(this,1,0);
      FluxBoxes fb_betanp1(this,1,0);
      MultiFab** betan = fb_betan.get();
      MultiFab** betanp1 = fb_betanp1.get();
      getViscosity(betan, time);
      getViscosity(betanp1, time+dt);

      int rhsComp  = 0;
      int betaComp = 0;
      diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                  delta_rhs,rhsComp,betan,viscn_cc,
                                  betanp1,viscnp1_cc,betaComp);

      delete delta_rhs;

      if (verbose > 1)
      {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "PeleLM::velocity_diffusion_update(): lev: " << level
                       << ", time: " << run_time << '\n';
      }
  }
}

void
PeleLM::getViscTerms (MultiFab& visc_terms,
                      int       src_comp,
                      int       num_comp,
                      Real      time)
{
   BL_PROFILE("PeleLM::getViscTerms()");
   const Real strt_time = ParallelDescriptor::second();
   //
   // Load "viscous" terms, starting from component = 0.
   //
   // JFG: for species, this procedure returns the *negative* of the divergence of
   // of the diffusive fluxes.  specifically, in the mixture averaged case, the
   // diffusive flux vector for species k is
   //
   //       j_k = - rho D_k,mix grad Y_k
   //
   // so the divergence of the flux, div dot j_k, has a negative in it.  instead
   // this procedure returns - div dot j_k to remove the negative.
   //
   // note the fluxes used in the code are extensive, that is, scaled by the areas
   // of the cell edges.  the calculation of the divergence is the sum of un-divided
   // differences of the extensive fluxes, all divided by volume.  so the effect is
   // to give the true divided difference approximation to the divergence of the
   // intensive flux.

   const int  nGrow     = visc_terms.nGrow();
   //
   // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
   //
   visc_terms.setBndry(1.e30);
   if (src_comp < AMREX_SPACEDIM)
   {
      if (src_comp != Xvel || num_comp < AMREX_SPACEDIM)
         amrex::Error("tensor v -> getViscTerms needs all v-components at once");

      FluxBoxes fb(this);
      MultiFab** vel_visc = fb.get();
      getViscosity(vel_visc, time);

      for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
         showMF("velVT",*(vel_visc[dir]),amrex::Concatenate("velVT_viscn_",dir,1),level);
      }
      int viscComp = 0;
      const TimeLevel whichTime = which_time(State_Type,time);
      AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
      auto vel_visc_cc = (whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
      diffusion->getTensorViscTerms(visc_terms,time,vel_visc,vel_visc_cc,viscComp);

      showMF("velVT",visc_terms,"velVT_visc_terms_1",level);
   }
   else
   {
     amrex::Abort("Should only call getViscTerms for velocity");
   }

   //
   // Ensure consistent grow cells
   //
   if (nGrow > 0)
   {
      visc_terms.FillBoundary(0,num_comp, geom.periodicity());
      AMREX_ASSERT(visc_terms.nGrow() == 1);
      Extrapolater::FirstOrderExtrap(visc_terms, geom, 0, num_comp);
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::getViscTerms(): lev: " << level
                     << ", time: " << run_time << '\n';
   }
}

void dumpProfileFab(const FArrayBox& fab,
                    const std::string file)
{
  const Box& box = fab.box();
  int imid = (int)(0.5*(box.smallEnd()[0] + box.bigEnd()[0]));
  IntVect iv1 = box.smallEnd(); iv1[0] = imid;
  IntVect iv2 = box.bigEnd(); iv2[0] = imid;
  iv1[1] = 115; iv2[1]=125;
  Box pb(iv1,iv2);
  amrex::Print() << "dumping over " << pb << '\n';

  std::ofstream osf(file.c_str());
  for (IntVect iv = pb.smallEnd(); iv <= pb.bigEnd(); pb.next(iv))
  {
    osf << iv[1] << " " << (iv[1]+0.5)*3.5/256 << " ";
    for (int n=0; n<fab.nComp(); ++n) 
      osf << fab(iv,n) << " ";
    osf << '\n';
  }
  amrex::Abort();
}

void dumpProfile(const MultiFab& mf,
                 const std::string file)
{
  const FArrayBox& fab = mf[0];
  dumpProfileFab(fab,file);
}

void
PeleLM::compute_differential_diffusion_fluxes (const MultiFab& S,
                                               const MultiFab* Scrse,
                                               MultiFab* const * flux,
                                               const MultiFab* const * beta,
                                               Real dt,
                                               Real time,
                                               bool include_Wbar_fluxes)
{
   BL_PROFILE("PeleLM:::compute_differential_diffusion_fluxes_msd()");
   const Real strt_time = ParallelDescriptor::second();

   if (hack_nospecdiff)
   {
      amrex::Error("compute_differential_diffusion_fluxes: hack_nospecdiff not implemented");
   }

   MultiFab  rh;      // Never need alpha for RhoY
   MultiFab* alpha_in        = 0;
   int       alpha_in_comp   = 0;
   int       fluxComp        = 0;
   int       betaComp        = 0;
   Real      a               = 0;
   Real      b               = 1.;
   bool      add_hoop_stress = false; // Only true if sigma == Xvel && Geometry::IsRZ())
//   bool      add_hoop_stress = true; // Only true if sigma == Xvel && Geometry::IsRZ())

// Some stupid check below to ensure we don't do stupid things
#if (AMREX_SPACEDIM == 3)
   // Here we ensure that R-Z related routines cannot be called in 3D
   if (add_hoop_stress){
      amrex::Abort("in diffuse_scalar: add_hoop_stress for R-Z geometry called in 3D !");
   }
#endif

#ifdef AMREX_USE_EB
   // Here we ensure that R-Z cannot work with EB (for now)
   if (add_hoop_stress){
      amrex::Abort("in diffuse_scalar: add_hoop_stress for R-Z geometry not yet working with EB support !");
   }
#endif

   const int nGrow = 1;
   AMREX_ASSERT(S.nGrow()>=nGrow);
   bool has_coarse_data = bool(Scrse);

   const DistributionMapping* dmc = (has_coarse_data ? &(Scrse->DistributionMap()) : 0);
   const BoxArray* bac = (has_coarse_data ? &(Scrse->boxArray()) : 0);

   MultiFab Soln(grids,dmap,1,nGrow,MFInfo(),Factory());
   auto Solnc = std::unique_ptr<MultiFab>(new MultiFab());
   if (has_coarse_data) {
      Solnc->define(*bac, *dmc, 1, nGrow,MFInfo(),getLevel(level-1).Factory());
   }

   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   info.setMaxCoarseningLevel(0);
#ifdef AMREX_USE_EB
   const auto& ebf = &(dynamic_cast<EBFArrayBoxFactory const&>(Factory()));
   MLEBABecLap op({geom}, {grids}, {dmap}, info, {ebf});
#else
   MLABecLaplacian op({geom}, {grids}, {dmap}, info);
#endif

   op.setMaxOrder(diffusion->maxOrder());
   MLMG mg(op);

   const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();
   const int rho_flag = 2;
   const BCRec& bc = theBCs[first_spec];
   for (int icomp=1; icomp<NUM_SPECIES; ++icomp) {
      AMREX_ALWAYS_ASSERT(bc == theBCs[first_spec+icomp]);
   }

   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
   Diffusion::setDomainBC(mlmg_lobc, mlmg_hibc, bc); // Same for all comps, by assumption
   op.setDomainBC(mlmg_lobc, mlmg_hibc);

   for (int icomp = 0; icomp < NUM_SPECIES+1; ++icomp)
   {
      const int sigma = first_spec + icomp;
      {
         if (has_coarse_data) {
            MultiFab::Copy(*Solnc,*Scrse,sigma,0,1,nGrow);
            if (rho_flag == 2) {
              MultiFab::Divide(*Solnc,*Scrse,Density,0,1,nGrow);
            }
            op.setCoarseFineBC(Solnc.get(), crse_ratio[0]);
         }
         MultiFab::Copy(Soln,S,sigma,0,1,nGrow);
         if (rho_flag == 2) {
            MultiFab::Divide(Soln,S,Density,0,1,nGrow);
         }
         op.setLevelBC(0, &Soln);
      }
      {
         Real* rhsscale = 0;
         // above sets rho_flag = 2;
         // this is potentially dangerous if rho_flag==1, because rh is an undefined MF
         const MultiFab& rho = (rho_flag == 1) ? rh : S;
         const int Rho_comp = (rho_flag ==1) ? 0 : Density;
         op.setScalars(a,b);
      }

      // Set beta as \rho D_m
      Diffusion::setBeta(op,beta,betaComp+icomp);

      // No multiplication by dt here.
      Diffusion::computeExtensiveFluxes(mg, Soln, flux, fluxComp+icomp, 1, &geom, b);
   }

   Soln.clear();

   if (use_wbar && include_Wbar_fluxes) {
      // Total flux : - \rho D_m \nabla Y_m (above) - \rho D_m Y_m / \overline{W} \nabla \overline{W} (below)
      compute_Wbar_fluxes(time,0,1.0);
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
         MultiFab::Add(*flux[d],*SpecDiffusionFluxWbar[d],0,fluxComp,NUM_SPECIES,0);
      }
   }

   //
   // Modify update/fluxes to preserve flux sum = 0 (conservatively correct Gamma_m)
   adjust_spec_diffusion_fluxes(flux, S, bc, time);

   // build heat fluxes based on species fluxes, Gamma_m, and cell-centered states
   // compute flux[NUM_SPECIES+1] = sum_m (H_m Gamma_m)
   // compute flux[NUM_SPECIES+2] = - lambda grad T
   compute_enthalpy_fluxes(flux,beta,time);

#ifdef AMREX_USE_EB
   // Get rid of flux on covered faces
   EB_set_covered_faces({AMREX_D_DECL(flux[0],flux[1],flux[2])},0.);
#endif

   // We have just computed "DD" and heat fluxes given an input state.
   //
   // Update DIFFUSIVE flux registers as follows.
   // If we are in the predictor and sdc_iterMAX>1:
   //   ADD -(1/2)*lambda^n GradT^n
   // If we are in the predictor and sdc_iterMAX=1:
   //   ADD -1.0*lambda^n GradT^n
   // If updateFluxReg=T (we are in the final corrector):
   //   ADD -(1/2)*lambda^(k) GradT^(k)
   //
   // Update ADVECTIVE flux registers as follows.
   // If we are in the predictor and sdc_iterMAX>1:
   //   ADD (1/2)*h_m^n Gamma_m^n
   // If we are in the predictor and sdc_iterMAX=1:
   //   ADD 1.0*h_m^n Gamma_m^n
   // If updateFluxReg=T (we are in the final corrector):
   //   ADD (1/2)*h_m^(k) Gamma_m^(k)

   if ( do_reflux && ( is_predictor || updateFluxReg ) )
   {
     if (is_predictor) {
        showMF("DBGSync",*flux[0],"DBGSync_DiffFluxX_Dn",level,parent->levelSteps(level));
        showMF("DBGSync",*flux[1],"DBGSync_DiffFluxY_Dn",level,parent->levelSteps(level));
     }
     if (updateFluxReg) {
        showMF("DBGSync",*flux[0],"DBGSync_DiffFluxX_Dnp1",level,parent->levelSteps(level));
        showMF("DBGSync",*flux[1],"DBGSync_DiffFluxY_Dnp1",level,parent->levelSteps(level));
     }
     for (int d = 0; d < AMREX_SPACEDIM; d++)
     {
       if (level > 0)
       {
         if (is_predictor)
         {
           const Real fac = (sdc_iterMAX==1) ? dt : 0.5*dt;
           getViscFluxReg().FineAdd(*flux[d],d,0,first_spec,NUM_SPECIES,fac);
           getAdvFluxReg().FineAdd(*flux[d],d,NUM_SPECIES+1,RhoH,1,fac);
           getViscFluxReg().FineAdd(*flux[d],d,NUM_SPECIES+2,RhoH,1,fac);
         }
         if (updateFluxReg)
         {
           getViscFluxReg().FineAdd(*flux[d],d,0,first_spec,NUM_SPECIES,-0.5*dt);
           getAdvFluxReg().FineAdd(*flux[d],d,NUM_SPECIES+1,RhoH,1,-0.5*dt);
           getViscFluxReg().FineAdd(*flux[d],d,NUM_SPECIES+2,RhoH,1,-0.5*dt);
         }
       }

       if (level < parent->finestLevel())
       {
         if (is_predictor)
         {
           const Real fac = (sdc_iterMAX==1) ? dt : 0.5*dt;
           getViscFluxReg(level+1).CrseInit((*flux[d]),d,0,first_spec,NUM_SPECIES,-fac,FluxRegister::ADD);
           getAdvFluxReg(level+1).CrseInit((*flux[d]),d,NUM_SPECIES+1,RhoH,1,-fac,FluxRegister::ADD);
           getViscFluxReg(level+1).CrseInit((*flux[d]),d,NUM_SPECIES+2,RhoH,1,-fac,FluxRegister::ADD);
         }
         if (updateFluxReg)
         {
           getViscFluxReg(level+1).CrseInit((*flux[d]),d,0,first_spec,NUM_SPECIES,0.5*dt,FluxRegister::ADD);
           getAdvFluxReg(level+1).CrseInit((*flux[d]),d,NUM_SPECIES+1,RhoH,1,0.5*dt,FluxRegister::ADD);
           getViscFluxReg(level+1).CrseInit((*flux[d]),d,NUM_SPECIES+2,RhoH,1,0.5*dt,FluxRegister::ADD);
         }
       }
     }
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::compute_differential_diffusion_fluxes(): lev: " << level
                     << ", time: " << run_time << '\n';
   }
}

void
PeleLM::scalar_advection_update (Real dt,
                                 int  first_scalar,
                                 int  last_scalar)
{
   //
   // Careful: If here, the sign of aofs is flipped (wrt the usual NS treatment).
   //
   MultiFab&       S_new = get_new_data(State_Type);
   const MultiFab& S_old = get_old_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
   for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx   = mfi.tilebox();
      auto const& snew   = S_new.array(mfi);
      auto const& sold   = S_old.array(mfi);
      auto const& adv    = aofs->array(mfi);

      amrex::ParallelFor(bx, [snew, sold, adv, dt, first_scalar, last_scalar]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         for (int n = first_scalar; n <= last_scalar; n++) {
            snew(i,j,k,n) = sold(i,j,k,n) + dt * adv(i,j,k,n);
         }
      });
   }

}

void
PeleLM::flux_divergence (MultiFab&        fdiv,
                         int              fdivComp,
                         const MultiFab* const* f,
                         int              fluxComp,
                         int              nComp,
                         Real             scale) const
{
   BL_PROFILE("PeleLM::flux_divergence()");
   AMREX_ASSERT(fdiv.nComp() >= fdivComp+nComp);

//////////////////////////////////////////////////////
// Version using AMReX computeDivergence function
// Will be updated later once more options are available in AMReX
//////////////////////////////////////////////////////
//   // Need aliases for the fluxes and divergence
//   Array<MultiFab const*, AMREX_SPACEDIM> flux_alias;
//   AMREX_D_TERM(flux_alias[0] = new MultiFab(*f[0], amrex::make_alias, fluxComp, nComp);,
//                flux_alias[1] = new MultiFab(*f[1], amrex::make_alias, fluxComp, nComp);,
//                flux_alias[2] = new MultiFab(*f[2], amrex::make_alias, fluxComp, nComp););
//
//   MultiFab div_alias(fdiv, amrex::make_alias, fdivComp, nComp);
//
//#ifdef AMREX_USE_EB
//   EB_computeDivergence(div_alias, flux_alias, geom, true);
//#else
//   computeDivergence(div_alias, flux_alias, geom);
//#endif
//   div_alias.mult(scale,0,nComp,0);


//////////////////////////////////////////////////////
//  PeleLM divergence function
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(fdiv,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      AMREX_D_TERM(auto const& fluxX = f[0]->array(mfi,fluxComp);,
                   auto const& fluxY = f[1]->array(mfi,fluxComp);,
                   auto const& fluxZ = f[2]->array(mfi,fluxComp););
      auto const& divergence   = fdiv.array(mfi,fdivComp);
      auto const& vol          = volume.const_array(mfi);

#ifdef AMREX_USE_EB
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
      auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
      auto const& flag    = flagfab.const_array();
#endif

#ifdef AMREX_USE_EB
      if (flagfab.getType(bx) == FabType::covered) {              // Covered boxes
         amrex::ParallelFor(bx, nComp, [divergence]
         AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
         {
            divergence(i,j,k,n) = 0.0;
         });
      } else if (flagfab.getType(bx) != FabType::regular ) {     // EB containing boxes 
         auto vfrac = ebfactory.getVolFrac().const_array(mfi);
         amrex::ParallelFor(bx, [nComp, flag, vfrac, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( flag(i,j,k).isCovered() ) {
               for (int n = 0; n < nComp; n++) {
                  divergence(i,j,k,n) = 0.0;
               }
            } else if ( flag(i,j,k).isRegular() ) {
               fluxDivergence( i, j, k, nComp,
                               AMREX_D_DECL(fluxX, fluxY, fluxZ),
                               vol, scale, divergence);
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               fluxDivergence( i, j, k, nComp,
                               AMREX_D_DECL(fluxX, fluxY, fluxZ),
                               vol, scale, divergence);
               for (int n = 0; n < nComp; n++) {
                  divergence(i,j,k,n) *= vfracinv;
               }
            }
         });
      } else {                                                   // Regular boxes
#endif
         amrex::ParallelFor(bx, [nComp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fluxDivergence( i, j, k, nComp,
                            AMREX_D_DECL(fluxX, fluxY, fluxZ),
                            vol, scale, divergence);
         });
#ifdef AMREX_USE_EB
      }
#endif
   }

#ifdef AMREX_USE_EB
   {
      MultiFab fdiv_SrcGhostCell(grids,dmap,nComp,fdiv.nGrow()+2,MFInfo(),Factory());
      fdiv_SrcGhostCell.setVal(0.);
      fdiv_SrcGhostCell.copy(fdiv, fdivComp, 0, nComp);
      amrex::single_level_weighted_redistribute( {fdiv_SrcGhostCell}, {fdiv}, *volfrac, fdivComp, nComp, {geom} );
   }
   EB_set_covered(fdiv,0.);
#endif
}

void
PeleLM::compute_differential_diffusion_terms (MultiFab& D,
                                              MultiFab& DD,
                                              Real      time,
                                              Real      dt,
                                              bool      include_Wbar_terms)
{
   BL_PROFILE("PeleLM:::compute_differential_diffusion_terms()");
   //
   // Sets vt for species, RhoH and Temp together
   // Uses state at time to explicitly compute fluxes, and resets internal
   //  data for fluxes, etc
   //
   AMREX_ASSERT(D.boxArray() == grids);
   AMREX_ASSERT(D.nComp() >= NUM_SPECIES+2); // room for spec+RhoH+Temp

   if (hack_nospecdiff)
   {
     amrex::Error("compute_differential_diffusion_terms: hack_nospecdiff not implemented");
   }

   const TimeLevel whichTime = which_time(State_Type,time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
   MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

   //
   // Compute/adjust species fluxes/heat flux/conduction, save in class data
   int nGrow = 1;
   int sComp = std::min((int)Density, std::min((int)first_spec,(int)Temp) );
   int eComp = std::max((int)Density, std::max((int)last_spec,(int)Temp) );
   int nComp = eComp - sComp + 1;
   FillPatch(*this,get_new_data(State_Type),nGrow,time,State_Type,sComp,nComp,sComp);
   std::unique_ptr<MultiFab> Scrse;
   if (level > 0) {
     auto& crselev = getLevel(level-1);
     Scrse.reset(new MultiFab(crselev.boxArray(), crselev.DistributionMap(), NUM_STATE, nGrow));
     FillPatch(crselev,*Scrse,nGrow,time,State_Type,Density,NUM_SPECIES+2,Density);
   }

   FluxBoxes fb_diff;
   MultiFab **beta = fb_diff.define(this,NUM_SPECIES+1);    // Local transport coeff face-centroid container
   getDiffusivity(beta, time, first_spec, 0, NUM_SPECIES);  // species (rhoD)
   getDiffusivity(beta, time, Temp, NUM_SPECIES, 1);        // temperature (lambda)

   compute_differential_diffusion_fluxes(get_new_data(State_Type),Scrse.get(),flux,beta,dt,time,include_Wbar_terms);

   D.setVal(0.0);
   DD.setVal(0.0);

   // Compute "D":
   // D[0:NUM_SPECIES-1] = -Div( Fk )
   // D[ NUM_SPECIES+1 ] = Div( lambda Grad(T) )
   flux_divergence(D,0,flux,0,NUM_SPECIES,-1.0);
   flux_divergence(D,NUM_SPECIES+1,flux,NUM_SPECIES+2,1,-1.0);

   // Compute "DD":
   // DD = -Sum{ Div( hk . Fk ) } a.k.a. the "diffdiff" terms
   flux_divergence(DD,0,flux,NUM_SPECIES+1,1,-1.0);

   if (D.nGrow() > 0 && DD.nGrow() > 0)
   {
     const int nc = NUM_SPECIES+2;

     D.FillBoundary(0,nc, geom.periodicity());
     DD.FillBoundary(0,1, geom.periodicity());

     AMREX_ASSERT(D.nGrow() == 1);
     AMREX_ASSERT(DD.nGrow() == 1);

     Extrapolater::FirstOrderExtrap(D, geom, 0, nc);
     Extrapolater::FirstOrderExtrap(DD, geom, 0, 1);
   }
}

void
PeleLM::state_stats (MultiFab& S)
{
   if (verbose) {
      //
      // Calculate some minimums and maximums.
      //
      auto scaleMin = VectorMin({&S},FabArrayBase::mfiter_tile_size,Density,NUM_STATE-AMREX_SPACEDIM,0);
      auto scaleMax = VectorMax({&S},FabArrayBase::mfiter_tile_size,Density,NUM_STATE-AMREX_SPACEDIM,0);

      bool aNegY = false;
      for (int i = 0; i < NUM_SPECIES && !aNegY; ++i) {
         if (scaleMin[first_spec+i-AMREX_SPACEDIM] < 0) aNegY = true;
      }

      amrex::Print() << "  Min,max temp = " << scaleMin[Temp   - AMREX_SPACEDIM]
                     << ", "                << scaleMax[Temp   - AMREX_SPACEDIM] << '\n';
      amrex::Print() << "  Min,max rho  = " << scaleMin[Density- AMREX_SPACEDIM]
                     << ", "                << scaleMax[Density- AMREX_SPACEDIM] << '\n';
      amrex::Print() << "  Min,max rhoh = " << scaleMin[RhoH   - AMREX_SPACEDIM]
                     << ", "                << scaleMax[RhoH   - AMREX_SPACEDIM] << '\n';

      if (aNegY){
         Vector<std::string> names;
         EOS::speciesNames(names);
         amrex::Print() << "  Species w/min < 0: ";
         for (int i = 0; i < NUM_SPECIES; ++i) {
            int idx = first_spec + i - AMREX_SPACEDIM;
            if ( scaleMin[idx] < 0) {
               amrex::Print() << "Y(" << names[i] << ") [" << scaleMin[idx] << "]  ";
            }
         }
         amrex::Print() << '\n';
      }
   }
}

void
PeleLM::compute_rhoRT (const MultiFab& S,
                             MultiFab& Press,
                             int       pComp)
{
   AMREX_ASSERT(pComp<Press.nComp());

   const Real strt_time = ParallelDescriptor::second();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = S.array(mfi,Density);
         auto const& rhoY    = S.array(mfi,first_spec);
         auto const& T       = S.array(mfi,Temp);
         auto const& P       = Press.array(mfi, pComp);

         amrex::ParallelFor(bx, [rho, rhoY, T, P]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getPGivenRTY( i, j, k, rho, rhoY, T, P );
         });
      }
   }
   
   if (verbose > 1)
   {
     const int IOProc   = ParallelDescriptor::IOProcessorNumber();
     Real      run_time = ParallelDescriptor::second() - strt_time;

     ParallelDescriptor::ReduceRealMax(run_time,IOProc);

     amrex::Print() << "PeleLM::compute_rhoRT(): lev: " << level
                    << ", time: " << run_time << '\n';
   }
}

//
// Setup for the advance function.
//

#ifndef NDEBUG
#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e30;
#endif
#endif

void
PeleLM::set_htt_hmixTYP ()
{
   const int finest_level = parent->finestLevel();

   // set typical value for hmix, needed for TfromHY solves if not provided explicitly
   if (typical_values[RhoH]==typical_RhoH_value_default)
   {
      htt_hmixTYP = 0;
      std::vector< std::pair<int,Box> > isects;
      for (int lvl = 0; lvl <= finest_level; lvl++)
      {
         AmrLevel&       ht = getLevel(lvl);
         const MultiFab& S  = ht.get_new_data(State_Type);
         const BoxArray& ba = ht.boxArray();
         const DistributionMapping& dm = ht.DistributionMap();
         MultiFab hmix(ba,dm,1,0,MFInfo(),Factory());
         MultiFab::Copy(hmix,S,RhoH,0,1,0);
         MultiFab::Divide(hmix,S,Density,0,1,0);
         if (lvl != finest_level)
         {
            AmrLevel& htf = getLevel(lvl+1);
            BoxArray  baf = htf.boxArray();
            baf.coarsen(parent->refRatio(lvl));
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(hmix,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               auto const& h_mix = hmix.array(mfi);
               baf.intersections(ba[mfi.index()],isects);
               for (int is = 0; is < isects.size(); is++) {
                  amrex::ParallelFor(isects[is].second, [h_mix]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     h_mix(i,j,k) = 0.0;
                  });
               }
            }
         }
         htt_hmixTYP = std::max(htt_hmixTYP,hmix.norm0(0));
      }
      ParallelDescriptor::ReduceRealMax(htt_hmixTYP);
      if (verbose)
      amrex::Print() << "setting htt_hmixTYP(via domain scan) = " << htt_hmixTYP << '\n';
   }
   else
   {
      htt_hmixTYP = typical_values[RhoH];
      if (verbose)
         amrex::Print() << "setting htt_hmixTYP(from user input) = " << htt_hmixTYP << '\n';
   }
}

void
PeleLM::advance_setup (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
{
   NavierStokesBase::advance_setup(time, dt, iteration, ncycle);

   for (int k = 0; k < num_state_type; k++)
   {
      MultiFab& nstate = get_new_data(k);
      MultiFab& ostate = get_old_data(k);
      MultiFab::Copy(nstate,ostate,0,0,nstate.nComp(),nstate.nGrow());
   }
   if (level == 0)
      set_htt_hmixTYP();

   make_rho_curr_time();

   RhoH_to_Temp(get_old_data(State_Type));

   if (plot_reactions && level == 0)
   {
      for (int i = parent->finestLevel(); i >= 0; --i)
         getLevel(i).auxDiag["REACTIONS"]->setVal(0);
   }
}

void
PeleLM::setThermoPress(Real time)
{
  const TimeLevel whichTime = which_time(State_Type,time);
    
  AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
  MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
    
  compute_rhoRT (S,S,RhoRT);
}

void
PeleLM::set_reasonable_grow_cells_for_R (Real time)
{
   //
   // Ensure reasonable grow cells for R.
   //
   MultiFab& React = get_data(RhoYdot_Type, time);
   React.FillBoundary(0,NUM_SPECIES, geom.periodicity());
   AMREX_ASSERT(React.nGrow() == 1);
   Extrapolater::FirstOrderExtrap(React, geom, 0, NUM_SPECIES); //FIXME: Is this in the wrong order?
}

Real
PeleLM::advance (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
{
  BL_PROFILE("PLM::advance()");

  //====================================
  BL_PROFILE_VAR("PeleLM::advance::mac", PLM_MAC);
  if (closed_chamber == 1 && level == 0)
  {
    // set new-time ambient pressure to be a copy of old-time ambient pressure
    p_amb_new = p_amb_old;
  }
  BL_PROFILE_VAR_STOP(PLM_MAC);
  //====================================

  BL_PROFILE_REGION_START("PeleLM::advance()");
  BL_PROFILE("PeleLM::advance()[src_sdc]");
  is_predictor = true;
  updateFluxReg = false;

  if (level == 0)
  {
    crse_dt = dt;
  }

  if (verbose)
  {
    amrex::Print() << "PeleLM::advance(): at start of time step\n"
                   << "SDC Advancing level " << level
                   << " : starting time = " << time
                   << " with dt = "         << dt << '\n';
  }

  //====================================
  BL_PROFILE_VAR("PeleLM::advance::setup", PLM_SETUP);
  // swaps old and new states for all state types
  // then copies each of the old state types into the new state types
  advance_setup(time,dt,iteration,ncycle);

  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& S_old = get_old_data(State_Type);

  const Real prev_time = state[State_Type].prevTime();
  const Real tnp1  = state[State_Type].curTime();

  //
  // Calculate the time N viscosity and diffusivity
  //   Note: The viscosity and diffusivity at time N+1 are
  //         initialized here to the time N values just to
  //         have something reasonable.
  //
  const int num_diff = NUM_STATE-AMREX_SPACEDIM-1;
  calcViscosity(prev_time,dt,iteration,ncycle);
  calcDiffusivity(prev_time);

  MultiFab::Copy(*viscnp1_cc, *viscn_cc, 0, 0, 1, viscn_cc->nGrow());
  MultiFab::Copy(*diffnp1_cc, *diffn_cc, 0, 0, num_diff, diffn_cc->nGrow());
  
  if (level==0 && reset_typical_vals_int>0)
  {
    int L0_steps = parent->levelSteps(0);
    if (L0_steps>0 && L0_steps%reset_typical_vals_int==0)
    {
      reset_typical_values(get_old_data(State_Type));
    }
  }

  if (do_check_divudt)
  {
    checkTimeStep(dt);
  }

  Real dt_test = 0.0;
  BL_PROFILE_VAR_STOP(PLM_SETUP);
  //====================================

  //====================================
  BL_PROFILE_VAR("PeleLM::advance::diffusion", PLM_DIFF);
  if (floor_species == 1)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
     {
        const Box& bx = mfi.tilebox();            
        auto const& rhoY    = S_old.array(mfi,first_spec);  
        amrex::ParallelFor(bx, [rhoY]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           fabMinMax( i, j, k, NUM_SPECIES, 0.0, Real_MAX, rhoY);
        });
     }
  }

  // Build a copy of old-time rho with grow cells for use in the diffusion solves
  make_rho_prev_time();
  BL_PROFILE_VAR_STOP(PLM_DIFF);
  //====================================

  //====================================
  BL_PROFILE_VAR_START(PLM_MAC);
  // compute old-time thermodynamic pressure
  setThermoPress(prev_time);  
  BL_PROFILE_VAR_STOP(PLM_MAC);
  //====================================

  //====================================
  BL_PROFILE_VAR_START(PLM_DIFF);
  // Compute Dn and DDn (based on state at tn)
  // (Note that coeffs at tn and tnp1 were intialized in above)
  MultiFab Dn(grids,dmap,NUM_SPECIES+2,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab DDn(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab DWbar;
  if (use_wbar) {
     DWbar.define(grids,dmap,NUM_SPECIES,nGrowAdvForcing,MFInfo(),Factory());
  }

  if (verbose) amrex::Print() << "Computing Dn & DDn \n";

  bool include_Wbar_terms = true;
  compute_differential_diffusion_terms(Dn,DDn,prev_time,dt,include_Wbar_terms);
  BL_PROFILE_VAR_STOP(PLM_DIFF);
  //====================================

  //====================================
  BL_PROFILE_VAR("PeleLM::advance::reactions", PLM_REAC);
  /*
    You could compute instantaneous I_R here but for now it's using either the
    previous step or divu_iter's version of I_R.  Either way, we have to make 
    sure that the nGrowAdvForcing grow cells have something reasonable in them
  */
  set_reasonable_grow_cells_for_R(tnp1);
  BL_PROFILE_VAR_STOP(PLM_REAC);
  //====================================

  // copy old state into new state for Dn and DDn.
  // Note: this was already done above for transport coefficients
  // and in advance_setup for scalar & divu

  MultiFab Dnp1(grids,dmap,NUM_SPECIES+2,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab DDnp1(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab Dhat(grids,dmap,NUM_SPECIES+2,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab DDhat(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab chi(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab chi_increment(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());
  MultiFab mac_divu(grids,dmap,1,nGrowAdvForcing,MFInfo(),Factory());

  //====================================
  BL_PROFILE_VAR_START(PLM_DIFF);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Dn,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
     const Box& gbx = mfi.growntilebox();
     auto const& dn       = Dn.array(mfi);
     auto const& dnp1k    = Dnp1.array(mfi);
     auto const& dnp1kp1  = Dhat.array(mfi);
     auto const& ddn      = DDn.array(mfi);
     auto const& ddnp1k   = DDnp1.array(mfi);
     auto const& ddnp1kp1 = DDhat.array(mfi);
     auto const& chi_ar      = chi.array(mfi);
     amrex::ParallelFor(gbx, [dn, dnp1k, dnp1kp1, ddn, ddnp1k, ddnp1kp1, chi_ar]
     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
     {
        for (int n = 0; n < NUM_SPECIES+2; n++) {
           dnp1k(i,j,k,n) = dn(i,j,k,n);
           dnp1kp1(i,j,k,n) = 0.0;
        }
        ddnp1k(i,j,k) = ddn(i,j,k);
        ddnp1kp1(i,j,k) = 0.0;
        chi_ar(i,j,k) = 0.0;
     });
  }
  BL_PROFILE_VAR_STOP(PLM_DIFF);
  //====================================

  is_predictor = false;

  BL_PROFILE_VAR_NS("PeleLM::advance::velocity_adv", PLM_VEL);

  for (int sdc_iter=1; sdc_iter<=sdc_iterMAX; ++sdc_iter)
  {

    if (sdc_iter == sdc_iterMAX)
      updateFluxReg = true;

    if (sdc_iter > 1)
    {
      //====================================
      // compute new-time transport coefficients
      BL_PROFILE_VAR_START(PLM_DIFF);
      calcDiffusivity(tnp1);

      // compute Dnp1 and DDnp1 iteratively lagged
      bool include_Wbar_terms_np1 = true;
      compute_differential_diffusion_terms(Dnp1,DDnp1,tnp1,dt,include_Wbar_terms_np1);
      BL_PROFILE_VAR_STOP(PLM_DIFF);
      //====================================

      //====================================
      // compute new-time DivU with instantaneous reaction rates
      BL_PROFILE_VAR_START(PLM_MAC);
      calc_divu(tnp1, dt, get_new_data(Divu_Type));
      BL_PROFILE_VAR_STOP(PLM_MAC);
      //====================================
    }

    //====================================
    // compute U^{ADV,*}
    BL_PROFILE_VAR_START(PLM_VEL);
    dt_test = predict_velocity(dt);
    BL_PROFILE_VAR_STOP(PLM_VEL);
    //====================================

    //====================================
    BL_PROFILE_VAR_START(PLM_MAC);
    // create S^{n+1/2} by averaging old and new

    MultiFab Forcing(grids,dmap,NUM_SPECIES+1,nGrowAdvForcing,MFInfo(),Factory());
    Forcing.setBndry(1.e30);
    FillPatch(*this,mac_divu,nGrowAdvForcing,time+0.5*dt,Divu_Type,0,1,0);

    // compute new-time thermodynamic pressure and chi_increment
    setThermoPress(tnp1);

    chi_increment.setVal(0.0,nGrowAdvForcing);
    calc_dpdt(tnp1,dt,chi_increment,u_mac);
    
#ifdef AMREX_USE_EB
    {
      MultiFab chi_tmp(grids,dmap,1,chi.nGrow()+2,MFInfo(),Factory());
      chi_tmp.setVal(0.);
      chi_tmp.copy(chi_increment);
      amrex::single_level_weighted_redistribute(  {chi_tmp}, {chi_increment}, *volfrac, 0, 1, {geom} );
      EB_set_covered(chi_increment,0.0);
    }
#endif

    // Add chi_increment to chi and add chi to time-centered mac_divu
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(chi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& gbx = mfi.growntilebox();
       auto const& chi_ar      = chi.array(mfi);
       auto const& chi_inc_ar  = chi_increment.array(mfi);
       auto const& divu_ar     = mac_divu.array(mfi);
       amrex::ParallelFor(gbx, [chi_ar, chi_inc_ar, divu_ar]
       AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
          chi_ar(i,j,k) += chi_inc_ar(i,j,k);
          divu_ar(i,j,k) += chi_ar(i,j,k);
       });
    }

    Real Sbar = 0;
    if (closed_chamber == 1 && level == 0)
    {
      Sbar = adjust_p_and_divu_for_closed_chamber(mac_divu);
    }

    // MAC-project... and overwrite U^{ADV,*}
    mac_project(time,dt,S_old,&mac_divu,nGrowAdvForcing,updateFluxReg);

    if (closed_chamber == 1 && level == 0 && Sbar != 0)
    {
      mac_divu.plus(Sbar,0,1); // add Sbar back to mac_divu
    }
    BL_PROFILE_VAR_STOP(PLM_MAC);
    //====================================

    //====================================
    // Scalar advection TPROF
    BL_PROFILE_VAR("PeleLM::advance::scal_advection", PLM_ADV);
    if (verbose) amrex::Print() << "A (SDC iter " << sdc_iter << ")\n";
    // Compute A (advection terms) with F = Dn + R
    //
    //  F[Temp] = ( Div(lambda.Grad(T)) - Sum{ hk.( R_k + Div( Fk ) ) } )/(rho.Cp)    NOTE: DD added below
    //          = ( D[N+1] + Sum{ hk.( R_k + D[k] ) } ) / (rho.Cp)
    //
    Forcing.setVal(0.0);

    int sComp = std::min(RhoH, std::min((int)Density, std::min((int)first_spec,(int)Temp) ) );
    int eComp = std::max(RhoH, std::max((int)Density, std::max((int)last_spec, (int)Temp) ) );
    int nComp = eComp - sComp + 1;

    FillPatchIterator S_fpi(*this,get_old_data(State_Type),nGrowAdvForcing,prev_time,State_Type,sComp,nComp);
    MultiFab& Smf=S_fpi.get_mf();

    int Rcomp = Density - sComp;
    int RYcomp = first_spec - sComp;
    int Tcomp = Temp - sComp;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Smf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& gbx = mfi.growntilebox();
       auto const& rho     = Smf.array(mfi,Rcomp);
       auto const& rhoY    = Smf.array(mfi,RYcomp);
       auto const& T       = Smf.array(mfi,Tcomp);
       auto const& dn      = Dn.array(mfi,0);
       auto const& ddn     = DDn.array(mfi);
       auto const& r       = get_new_data(RhoYdot_Type).array(mfi,0);
       auto const& fY      = Forcing.array(mfi,0);
       auto const& fT      = Forcing.array(mfi,NUM_SPECIES);
       Real        dp0dt_d = dp0dt;
       int     closed_ch_d = closed_chamber;

       amrex::ParallelFor(gbx, [rho, rhoY, T, dn, ddn, r, fY, fT, dp0dt_d, closed_ch_d]
       AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
          buildAdvectionForcing( i, j, k, rho, rhoY, T, dn, ddn, r, dp0dt_d, closed_ch_d, fY, fT );
       });
    }

    Forcing.FillBoundary(0,NUM_SPECIES+1,geom.periodicity());

    aofs->setVal(1.e30,aofs->nGrow());

    compute_scalar_advection_fluxes_and_divergence(Forcing,mac_divu,dt);
    showMF("DBGSync",u_mac[0],"DBGSync_umacX",level,sdc_iter,parent->levelSteps(level));
    showMF("DBGSync",u_mac[1],"DBGSync_umacY",level,sdc_iter,parent->levelSteps(level));
    
    // update rho since not affected by diffusion
    scalar_advection_update(dt, Density, Density);
    make_rho_curr_time();

    BL_PROFILE_VAR_STOP(PLM_ADV);
    // Close Scalar advection TPROF
    //==================================== 

    //==================================== 
    // Open Scalar diffusion TPROF
    BL_PROFILE_VAR_START(PLM_DIFF);
    if (verbose) amrex::Print() << "D (SDC iter " << sdc_iter << ")\n";
    // Compute Dhat, diffuse with F
    //                 = A + R + 0.5(Dn + Dnp1) - Dnp1 + Dhat + 0.5(DDn + DDnp1) - DDnp1 + DDHat
    //                 = A + R + 0.5(Dn - Dnp1) + Dhat + 0.5(DDn - DDnp1) + DDhat
    // 

    Forcing.setVal(0.0);

    // Get the Wbar term is required
    if (use_wbar) {
       // Compute Wbar fluxes from state np1k (lagged)
       // Compute DWbar term: - \nabla \cdot \Gamma_{\overline{W}_m}
       // They will be added to the Fnp1kp1 directly in diff_diff_update()
       const Real  cur_time  = state[State_Type].curTime();
       compute_Wbar_fluxes(cur_time,0,1.0);
       flux_divergence(DWbar,0,SpecDiffusionFluxWbar,0,NUM_SPECIES,-1);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      for (MFIter mfi(Forcing,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& bx = mfi.tilebox();
        auto const& dn      = Dn.array(mfi,0);
        auto const& ddn     = DDn.array(mfi);
        auto const& dnp1k   = Dnp1.array(mfi,0);
        auto const& ddnp1k  = DDnp1.array(mfi);
        auto const& r       = get_new_data(RhoYdot_Type).array(mfi,0);
        auto const& a       = aofs->array(mfi,first_spec);
        auto const& fY      = Forcing.array(mfi,0);
        auto const& fT      = Forcing.array(mfi,NUM_SPECIES);
        int use_wbar_lcl    = use_wbar;
        auto const& dwbar   = (use_wbar) ? DWbar.array(mfi) : Dn.array(mfi,0);

        Real        dp0dt_d = dp0dt;
        int     closed_ch_d = closed_chamber;

        amrex::ParallelFor(bx, [dn, ddn, dnp1k, ddnp1k, dwbar,
                                r, a, fY, fT, dp0dt_d, closed_ch_d, use_wbar_lcl]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           buildDiffusionForcing( i, j, k, dn, ddn, dnp1k, ddnp1k, r, a, dp0dt_d, closed_ch_d, fY, fT );
           if (use_wbar_lcl) {
              for (int n = 0; n < NUM_SPECIES; n++) {
                 fY(i,j,k,n) += dwbar(i,j,k,n);
              }
           }
        });
      }
    }

#ifdef AMREX_USE_EB    
    EB_set_covered(Forcing,0.);
#endif

    differential_diffusion_update(Forcing,DWbar,0,Dhat,0,DDhat);
    BL_PROFILE_VAR_STOP(PLM_DIFF);
    // Close Scalar diffusion TPROF
    //==================================== 

    //==================================== 
    BL_PROFILE_VAR_START(PLM_REAC);
    if (verbose) amrex::Print() << "R (SDC iter " << sdc_iter << ")\n";
    /////////////////////
    // Compute R (F = A + 0.5(Dn - Dnp1 + DDn - DDnp1) + Dhat + DDhat )
    // hack: for advance_chemistry, use same Forcing used for species eqn (just subtract S_old and the omegaDot term)
    /////////////////////
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
       for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          const Box& bx = mfi.tilebox();
          auto const& sold    = S_old.array(mfi,first_spec);
          auto const& snew    = S_new.array(mfi,first_spec);
          auto const& r       = get_new_data(RhoYdot_Type).array(mfi);
          auto const& force   = Forcing.array(mfi);
          amrex::Real dtinv   = 1.0/dt;
          amrex::ParallelFor(bx, [sold, snew, r, force, dtinv]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
             for (int n = 0; n < NUM_SPECIES; n++) {
                force(i,j,k,n) = ( snew(i,j,k,n) - sold(i,j,k,n) ) * dtinv - r(i,j,k,n);
             }
             force(i,j,k,NUM_SPECIES) = ( snew(i,j,k,NUM_SPECIES) - sold(i,j,k,NUM_SPECIES) ) * dtinv;
          });
       }
    }

    showMF("mysdc",Forcing,"sdc_F_befAdvChem",level,sdc_iter,parent->levelSteps(level));
    showMF("mysdc",S_old,"sdc_Sold_befAdvChem",level,sdc_iter,parent->levelSteps(level));
    showMF("mysdc",S_new,"sdc_Snew_befAdvChem",level,sdc_iter,parent->levelSteps(level));

    advance_chemistry(S_old,S_new,dt,Forcing,0);

#ifdef AMREX_USE_EB
    set_body_state(S_new);
#endif
    RhoH_to_Temp(S_new);
    BL_PROFILE_VAR_STOP(PLM_REAC);
    // Close species chemistry TPROF
    //====================================

    //====================================
    BL_PROFILE_VAR_START(PLM_DIFF);
    if (floor_species == 1)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          const Box& bx = mfi.tilebox();            
          auto const& rhoY    = S_new.array(mfi,first_spec);  
          amrex::ParallelFor(bx, [rhoY]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
             fabMinMax( i, j, k, NUM_SPECIES, 0.0, Real_MAX, rhoY);
          });
       }
    }
    BL_PROFILE_VAR_STOP(PLM_DIFF);
    //====================================

    showMF("mysdc",S_new,"sdc_Snew_end_sdc",level,sdc_iter,parent->levelSteps(level));
    showMF("mysdc",S_old,"sdc_Sold_end_sdc",level,sdc_iter,parent->levelSteps(level));

    state_stats(S_new);
    if (verbose) amrex::Print() << "DONE WITH (SDC corrector " << sdc_iter << ")\n";

    //====================================
    BL_PROFILE_VAR_START(PLM_MAC);
    setThermoPress(tnp1);
    BL_PROFILE_VAR_STOP(PLM_MAC);
    //====================================

    showMF("DBGSync",S_new,"DBGSync_Snew_end_sdc",level,sdc_iter,parent->levelSteps(level));
  }

   Dn.clear();
   DDn.clear();
   Dnp1.clear();
   DDnp1.clear();
   Dhat.clear();
   DDhat.clear();
   chi_increment.clear();

   if (verbose) amrex::Print() << " SDC iterations complete \n";

   if (plot_consumption)
   {
     for (int j=0; j<consumptionName.size(); ++j)
     {
       int consumptionComp = getSpeciesIdx(consumptionName[j]);
       MultiFab::Copy((*auxDiag["CONSUMPTION"]),get_new_data(RhoYdot_Type),consumptionComp,j,1,0);
       auxDiag["CONSUMPTION"]->mult(-1,j,1); // Convert production to consumption
     }
   }

   if (plot_heat_release)
   {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
         FArrayBox EnthFab;
         for (MFIter mfi((*auxDiag["HEATRELEASE"]),TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            EnthFab.resize(bx,NUM_SPECIES); 
            Elixir  Enthi   = EnthFab.elixir();
            auto const& T   = get_new_data(State_Type).array(mfi,Temp);
            auto const& Hi  = EnthFab.array();
            auto const& r   = get_new_data(RhoYdot_Type).array(mfi);
            auto const& HRR = (*auxDiag["HEATRELEASE"]).array(mfi);

            amrex::ParallelFor(bx, [T, Hi, HRR, r]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               getHGivenT( i, j, k, T, Hi );
               HRR(i,j,k) = 0.0;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  HRR(i,j,k) -= Hi(i,j,k,n) * r(i,j,k,n);
               } 
            });
         }
      }
   }

   //====================================
   BL_PROFILE_VAR_START(PLM_DIFF);
   // Update transport coefficient with post-SDC state
   calcDiffusivity(tnp1);
   calcViscosity(tnp1,dt,iteration,ncycle);
   BL_PROFILE_VAR_STOP(PLM_DIFF);
   //====================================

   //
   // Set the dependent value of RhoRT to be the thermodynamic pressure.  By keeping this in
   // the state, we can use the average down stuff to be sure that RhoRT_avg is avg(RhoRT),
   // not ave(Rho)avg(R)avg(T), which seems to give the p-relax stuff in the mac Rhs troubles.
   //
   //====================================
   BL_PROFILE_VAR_START(PLM_MAC);
   setThermoPress(tnp1);
   BL_PROFILE_VAR_STOP(PLM_MAC);
   //====================================

   //====================================
   BL_PROFILE_VAR("PeleLM::advance::project", PLM_PROJ);
   calc_divu(time+dt, dt, get_new_data(Divu_Type));

   if (!NavierStokesBase::initial_step && level != parent->finestLevel())
   {
      //
      // Set new divu to old divu where covered by fine.
      //
      BoxArray crsndgrids = getLevel(level+1).grids;
      crsndgrids.coarsen(fine_ratio);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif            
      for (MFIter mfi(get_new_data(Divu_Type),TilingIfNotGPU()); mfi.isValid();++mfi)
      {
         auto isects = crsndgrids.intersections(mfi.tilebox());
         for (int it = 0, N = isects.size(); it < N; it++)
         {
            const Box& ovlp = isects[it].second;
            auto const& divu_n  = get_new_data(Divu_Type).array(mfi);
            auto const& divu_o  = get_old_data(Divu_Type).array(mfi);
            amrex::ParallelFor(ovlp, [divu_n, divu_o]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               divu_n(i,j,k) = divu_o(i,j,k);
            });
         }
      }
   }
   BL_PROFILE_VAR_STOP(PLM_PROJ);
   //====================================

   //====================================
   BL_PROFILE_VAR_START(PLM_VEL);
   //
   // Add the advective and other terms to get velocity (or momentum) at t^{n+1}.
   //
   
   if (do_mom_diff == 0) {
      velocity_advection(dt);
   }
   
   velocity_update(dt);
   BL_PROFILE_VAR_STOP(PLM_VEL);
   //====================================

   //====================================
   BL_PROFILE_VAR_START(PLM_PROJ);
   // compute chi correction
   // place to take dpdt stuff out of nodal project
   //    calc_dpdt(tnp1,dt,chi_increment,u_mac);
   //    MultiFab::Add(get_new_data(Divu_Type),chi_increment,0,0,1,0);

   // subtract mean from divu
   Real Sbar_old = 0;
   Real Sbar_new = 0;
   if (closed_chamber == 1 && level == 0)
   {
      MultiFab& divu_old = get_old_data(Divu_Type);
      MultiFab& divu_new = get_new_data(Divu_Type);

      // compute number of cells
      Real num_cells = grids.numPts();

      // compute average of S at old and new times
      Sbar_old = divu_old.sum() / num_cells;
      Sbar_new = divu_new.sum() / num_cells;

      // subtract mean from divu
      divu_old.plus(-Sbar_old,0,1);
      divu_new.plus(-Sbar_new,0,1);
   }

   //
   // Increment rho average.
   //
   if (!initial_step)
   {
      if (level > 0)
      {
         Real alpha = 1.0/Real(ncycle);
         if (iteration == ncycle)
           alpha = 0.5/Real(ncycle);
         incrRhoAvg(alpha);
      }

      //
      // Do a level project to update the pressure and velocity fields.
      //

      level_projector(dt,time,iteration);

      // restore Sbar
      if (closed_chamber == 1 && level == 0)
      {
         MultiFab& divu_old = get_old_data(Divu_Type);
         MultiFab& divu_new = get_new_data(Divu_Type);

         divu_old.plus(Sbar_old,0,1);
         divu_new.plus(Sbar_new,0,1);
      }

      if (level > 0 && iteration == 1) p_avg.setVal(0);
   }
   BL_PROFILE_VAR_STOP(PLM_PROJ);
   //====================================

#ifdef AMREX_PARTICLES
   if (theNSPC() != 0 && !initial_step)
   {
      theNSPC()->AdvectWithUmac(u_mac, level, dt);
   }
#endif

   //====================================
   BL_PROFILE_VAR("PeleLM::advance::cleanup", PLM_CLEANUP);
   advance_cleanup(iteration,ncycle);
   BL_PROFILE_VAR_STOP(PLM_CLEANUP);
   //====================================

   //
   // Update estimate for allowable time step.
   //
   if (fixed_dt > 0)
   {
      dt_test = estTimeStep();
   }
   else
   {
      dt_test = std::min(dt_test, estTimeStep());
   }

   if (verbose) amrex::Print() << "PeleLM::advance(): at end of time step\n";

   state_stats(S_new);
   
   //====================================
   BL_PROFILE_VAR_START(PLM_MAC);
   // during initialization, reset time 0 ambient pressure
   if (closed_chamber == 1 && level == 0 && !initial_step)
   {
      p_amb_old = p_amb_new;
   }
   BL_PROFILE_VAR_STOP(PLM_MAC);
   //====================================

   BL_PROFILE_REGION_STOP("PeleLM::advance()");

   return dt_test;
}

Real
PeleLM::adjust_p_and_divu_for_closed_chamber(MultiFab& mac_divu)
{
   MultiFab& S_new = get_new_data(State_Type);
   MultiFab& S_old = get_old_data(State_Type);

   const Real prev_time = state[State_Type].prevTime();
   const Real cur_time  = state[State_Type].curTime();
   const Real dt = cur_time - prev_time;

   // used for closed chamber algorithm
   MultiFab theta_halft(grids,dmap,1,nGrowAdvForcing);

   // compute old, new, and time-centered theta = 1 / (gamma P)
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();            
      auto const& rhoY_o  = S_old.array(mfi,first_spec);
      auto const& rhoY_n  = S_new.array(mfi,first_spec);
      auto const& T_o     = S_old.array(mfi,Temp);
      auto const& T_n     = S_new.array(mfi,Temp);
      auto const& theta   = theta_halft.array(mfi);
      amrex::Real pamb_o  = p_amb_old;
      amrex::Real pamb_n  = p_amb_new;

      amrex::ParallelFor(bx, [rhoY_o, rhoY_n, T_o, T_n, theta, pamb_o, pamb_n]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         amrex::Real gammaInv_o = getGammaInv( i, j, k, rhoY_o, T_o );
         amrex::Real gammaInv_n = getGammaInv( i, j, k, rhoY_n, T_n );
         theta(i,j,k) = 0.5 * ( gammaInv_o/pamb_o + gammaInv_n/pamb_n );
      });
   }

   // compute number of cells
   Real num_cells = grids.numPts();

   // compute the average of mac_divu theta
   Real Sbar = mac_divu.sum() / num_cells;
   thetabar = theta_halft.sum() / num_cells;

   // subtract mean from mac_divu and theta_nph
   mac_divu.plus(-Sbar,0,1);
   theta_halft.plus(-thetabar,0,1);

   p_amb_new = p_amb_old + dt*(Sbar/thetabar);
   dp0dt = Sbar/thetabar;

   // update mac rhs by adding delta_theta * (Sbar / thetabar)
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(mac_divu,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& m_divu  = mac_divu.array(mfi);
      auto const& theta   = theta_halft.array(mfi);
      amrex::Real scaling = Sbar/thetabar;

      amrex::ParallelFor(bx, [m_divu, theta, scaling]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         m_divu(i,j,k) -= theta(i,j,k) * scaling;
      });
   }

   amrex::Print() << "level 0: prev_time, p_amb_old, p_amb_new, delta = " 
                  << prev_time << " " << p_amb_old << " " << p_amb_new << " "
                  << p_amb_new-p_amb_old << std::endl;

   return Sbar;
}

DistributionMapping
PeleLM::getFuncCountDM (const BoxArray& bxba, int ngrow)
{
  //
  // Sometimes "mf" is the valid region of the State.
  // Sometimes it's the region covered by AuxBoundaryData.
  // When ngrow>0 were doing AuxBoundaryData with nGrow()==ngrow.
  //
  DistributionMapping rr;
  rr.RoundRobinProcessorMap(bxba.size(),ParallelDescriptor::NProcs());

  MultiFab fctmpnew(bxba, rr, 1, 0);
  fctmpnew.setVal(1);

  const MultiFab& FC = get_new_data(FuncCount_Type);

  // Return original DM if FC is zero (propably Null reaction dir)
  if ( FC.norm0(0) <= 0.0 ) {
     DistributionMapping new_dmap{bxba};
     return new_dmap;
  }

  fctmpnew.copy(FC,0,0,1,std::min(ngrow,FC.nGrow()),0);

  int count = 0;
  Vector<long> vwrk(bxba.size());

// FIXME need to be tiled  
  for (MFIter mfi(fctmpnew); mfi.isValid(); ++mfi)
     vwrk[count++] = static_cast<long>(fctmpnew[mfi].sum<RunOn::Host>(0));

  fctmpnew.clear();

#if AMREX_USE_MPI
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  Vector<int> nmtags(ParallelDescriptor::NProcs(),0);
  Vector<int> offset(ParallelDescriptor::NProcs(),0);

  const Vector<int>& procmap = rr.ProcessorMap();

  for (int i = 0; i < vwrk.size(); i++)
    nmtags[procmap[i]]++;

  AMREX_ASSERT(nmtags[ParallelDescriptor::MyProc()] == count);

  for (int i = 1; i < offset.size(); i++)
    offset[i] = offset[i-1] + nmtags[i-1];

  Vector<long> vwrktmp;

  if (ParallelDescriptor::IOProcessor()) vwrktmp = vwrk;

  MPI_Gatherv(vwrk.dataPtr(),
              count,
              ParallelDescriptor::Mpi_typemap<long>::type(),
              ParallelDescriptor::IOProcessor() ? vwrktmp.dataPtr() : 0,
              nmtags.dataPtr(),
              offset.dataPtr(),
              ParallelDescriptor::Mpi_typemap<long>::type(),
              IOProc,
              ParallelDescriptor::Communicator());

  if (ParallelDescriptor::IOProcessor())
  {
    //
    // We must now assemble vwrk in the proper order.
    //
    std::vector< std::vector<int> > table(ParallelDescriptor::NProcs());

    for (int i = 0; i < vwrk.size(); i++)
      table[procmap[i]].push_back(i);

    int idx = 0;
    for (int i = 0; i < table.size(); i++)
    {
      std::vector<int>& tbl = table[i];
      for (int j = 0; j < tbl.size(); j++)
        vwrk[tbl[j]] = vwrktmp[idx++];
    }
  }

  vwrktmp.clear();
  //
  // Send the properly-ordered vwrk to all processors.
  //
  ParallelDescriptor::Bcast(vwrk.dataPtr(), vwrk.size(), IOProc);
#endif

  DistributionMapping res;

  res.KnapSackProcessorMap(vwrk,ParallelDescriptor::NProcs());

  return res;
}

void
PeleLM::advance_chemistry (MultiFab&       mf_old,
                           MultiFab&       mf_new,
                           Real            dt,
                           const MultiFab& Force,
                           int             nCompF,
                           bool            use_stiff_solver)
{
  BL_PROFILE("PLM::advance_chemistry()");
  
  const Real strt_time = ParallelDescriptor::second();

  const bool do_avg_down_chem = avg_down_chem
    && level < parent->finestLevel()
    && getLevel(level+1).state[RhoYdot_Type].hasOldData();

  if (hack_nochem)
  {
    MultiFab::Copy(mf_new,mf_old,first_spec,first_spec,NUM_SPECIES+2,0);
    MultiFab::Saxpy(mf_new,dt,Force,0,first_spec,NUM_SPECIES+1,0);
    get_new_data(RhoYdot_Type).setVal(0);
    get_new_data(FuncCount_Type).setVal(0);
  }
  else
  {
    BoxArray cf_grids;

    if (do_avg_down_chem)
    {
      cf_grids = getLevel(level+1).boxArray(); cf_grids.coarsen(fine_ratio);
    }

    MultiFab&  React_new = get_new_data(RhoYdot_Type);
    const int  ngrow     = std::min(std::min(React_new.nGrow(),mf_old.nGrow()),mf_new.nGrow()); 

    BoxArray  ba           = mf_new.boxArray();
    DistributionMapping dm = mf_new.DistributionMap();
#ifndef AMREX_USE_GPU
    //
    // Chop the grids to level out the chemistry work when on the CPU.
    // We want enough grids so that KNAPSACK works well,
    // but not too many to make unweildy BoxArrays.
    //
    const int Threshold = chem_box_chop_threshold * ParallelDescriptor::NProcs();
    bool      done      = (ba.size() >= Threshold);

    for (int cnt = 1; !done; cnt *= 2)
    {
      const IntVect ChunkSize = parent->maxGridSize(level)/cnt;

      if ( AMREX_D_TERM(ChunkSize[0] < 16, || ChunkSize[1] < 16, || ChunkSize[2] < 16) )
        //
        // Don't let grids get too small.
        //
        break;

      IntVect chunk(ChunkSize);

      for (int j = AMREX_SPACEDIM-1; j >=0 && ba.size() < Threshold; j--)
      {
        chunk[j] /= 2;
        ba.maxSize(chunk);
        if (ba.size() >= Threshold) done = true;
      }
    }

    dm = getFuncCountDM(ba,ngrow);
#endif

    MultiFab STemp(ba, dm, NUM_SPECIES+3, 0);
    MultiFab FTemp(ba, dm, Force.nComp(), 0);

    STemp.copy(mf_old,first_spec,0,NUM_SPECIES+3); // Parallel copy.
    FTemp.copy(Force);                          // Parallel copy.

    MultiFab fcnCntTemp(ba, dm, 1, 0);
    MultiFab diagTemp;

    // Setup a mask for chemistry. Right used only for EB.
    // TODO: find a way to set that up properly when running boxes on GPU.
    amrex::FabArray<amrex::BaseFab<int>>  react_mask;
    react_mask.define(ba, dm,  1, 0);
#ifdef AMREX_USE_EB     
    react_mask.copy(ebmask);
#else
    react_mask.setVal(1);
#endif
    
    const bool do_diag = plot_reactions && amrex::intersect(ba,auxDiag["REACTIONS"]->boxArray()).size() != 0;

    if (do_diag)
    {
        diagTemp.define(ba, dm, auxDiag["REACTIONS"]->nComp(), 0);
        diagTemp.copy(*auxDiag["REACTIONS"]); // Parallel copy
    }

    if (verbose) 
      amrex::Print() << "*** advance_chemistry: FABs in tmp MF: " << STemp.size() << '\n';

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter Smfi(STemp,amrex::TilingIfNotGPU()); Smfi.isValid(); ++Smfi)
    {
        const Box& bx       = Smfi.tilebox();
        auto const& rhoY     = STemp.array(Smfi);
        auto const& rhoH     = STemp.array(Smfi,NUM_SPECIES);
        auto const& temp     = STemp.array(Smfi,NUM_SPECIES+1);
        auto const& fcl      = fcnCntTemp.array(Smfi);
        auto const& frc_rhoY = FTemp.array(Smfi);
        auto const& frc_rhoH = FTemp.array(Smfi, NUM_SPECIES);
        auto const& mask     = react_mask.array(Smfi);
        int ncells           = bx.numPts();

        // Convert MKS -> CGS
        ParallelFor(bx, [rhoY,rhoH, frc_rhoY, frc_rhoH]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           for (int n = 0; n < NUM_SPECIES; n++) {
              rhoY(i,j,k,n) *= 1.0e-3;
              frc_rhoY(i,j,k,n) *= 1.0e-3;
           }
           rhoH(i,j,k) *= 10.0;
           frc_rhoH(i,j,k) *= 10.0;
        });

#ifdef AMREX_USE_GPU
        const auto ec = Gpu::ExecutionConfig(ncells);
#endif

        BL_PROFILE_VAR("React()", ReactInLoop);
        Real dt_incr     = dt;
        Real time_chem   = 0;
#ifndef AMREX_USE_GPU
        /* Solve */
        int tmp_fctCn = react(bx, rhoY, frc_rhoY, temp, rhoH, frc_rhoH, fcl, mask, dt_incr, time_chem);
        dt_incr   = dt;
        time_chem = 0;
#else
        int reactor_type = 2;
        int tmp_fctCn = react(bx, rhoY, frc_rhoY, temp, rhoH, frc_rhoH, fcl, mask, 
                              dt_incr, time_chem, reactor_type, amrex::Gpu::gpuStream());
        dt_incr = dt;
        time_chem = 0;
#endif
        BL_PROFILE_VAR_STOP(ReactInLoop);

        // Convert CGS -> MKS
        ParallelFor(bx, [rhoY,rhoH]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
        {
           for (int n = 0; n < NUM_SPECIES; n++) {
              rhoY(i,j,k,n) *= 1.0e3;
           }
           rhoH(i,j,k) *= 0.1;
        });

#ifdef AMREX_USE_GPU
        Gpu::Device::streamSynchronize();
#endif

    }

    FTemp.clear();

    mf_new.copy(STemp,0,first_spec,NUM_SPECIES+3); // Parallel copy.

    STemp.clear();

    //
    // Set React_new (I_R).
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_old,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx         = mfi.tilebox();
        auto const& rhoY_o    = mf_old.array(mfi,first_spec);
        auto const& rhoY_n    = mf_new.array(mfi,first_spec);
        auto const& frc_rhoY  = Force.array(mfi);
        auto const& rhoYdot   = React_new.array(mfi);
        Real dt_inv = 1.0/dt;
        ParallelFor(bx, NUM_SPECIES, [rhoY_o, rhoY_n, frc_rhoY, rhoYdot, dt_inv]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           rhoYdot(i,j,k,n) = - ( rhoY_o(i,j,k,n) - rhoY_n(i,j,k,n) ) * dt_inv - frc_rhoY(i,j,k,n);
        });
    }

    if (do_diag)
    {
      auxDiag["REACTIONS"]->copy(diagTemp); // Parallel copy
      diagTemp.clear();
    }

    MultiFab& FC = get_new_data(FuncCount_Type);
    FC.copy(fcnCntTemp,0,0,1,0,std::min(ngrow,FC.nGrow()));
    fcnCntTemp.clear();
    //
    // Approximate covered crse chemistry (I_R) with averaged down fine I_R from previous time step.
    //
    if (do_avg_down_chem)
    {
      MultiFab& fine_React = getLevel(level+1).get_old_data(RhoYdot_Type);
      average_down(fine_React, React_new, 0, NUM_SPECIES);
    }
    //
    // Ensure consistent grow cells.
    //
    if (ngrow > 0)
    {
      AMREX_ASSERT(React_new.nGrow() == 1);
      React_new.FillBoundary(0,NUM_SPECIES, geom.periodicity());
      Extrapolater::FirstOrderExtrap(React_new, geom, 0, NUM_SPECIES);
    }
  }

  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real mx = ParallelDescriptor::second() - strt_time, mn = mx;

    ParallelDescriptor::ReduceRealMin(mn,IOProc);
    ParallelDescriptor::ReduceRealMax(mx,IOProc);

    amrex::Print() << "PeleLM::advance_chemistry(): lev: " << level << ", time: ["
                   << mn << " ... " << mx << "]\n";
  }
}

void
PeleLM::compute_scalar_advection_fluxes_and_divergence (const MultiFab& Force,
                                                        const MultiFab& DivU,
                                                        Real            dt)
{

  BL_PROFILE("PLM::comp_sc_adv_fluxes_and_div()");
  //
  // Compute -Div(advective fluxes)  [ which is -aofs in NS, BTW ... careful...
  //
  if (verbose) amrex::Print() << "... computing advection terms\n";

  const Real  strt_time = ParallelDescriptor::second();
  const Real* dx        = geom.CellSize();
  const Real  prev_time = state[State_Type].prevTime();  

  int ng = godunov_hyp_grow;
  int sComp = std::min((int)Density, std::min((int)first_spec,(int)Temp) );
  int eComp = std::max((int)Density, std::max((int)last_spec,(int)Temp) );
  int nComp = eComp - sComp + 1;

  FillPatchIterator S_fpi(*this,get_old_data(State_Type),ng,prev_time,State_Type,sComp,nComp);
  MultiFab& Smf=S_fpi.get_mf();
  
  int rhoYcomp = first_spec - sComp;
  int Rcomp = Density - sComp;
  int Tcomp = Temp - sComp;

  // Floor small values of states to be extrapolated
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Smf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
     const Box& gbx = mfi.growntilebox(godunov_hyp_grow);
     auto fab = Smf.array(mfi);
     AMREX_HOST_DEVICE_FOR_4D ( gbx, NUM_SPECIES, i, j, k, n,
     {
        auto val = fab(i,j,k,n+Rcomp);
        val = std::abs(val) > 1.e-20 ? val : 0.0;
     });
  }

#ifdef AMREX_USE_EB
  //////////////////////////////////////
  //
  // HERE IS THE EB PROCEDURE
  //
  //////////////////////////////////////

  // Initialize accumulation for rho = Sum(rho.Y)
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    EdgeState[d]->setVal(0);
    EdgeFlux[d]->setVal(0);
  }

  // Advect RhoY  
  {
     Vector<BCRec> math_bcs(NUM_SPECIES);
     math_bcs = fetchBCArray(State_Type, first_spec,NUM_SPECIES);
     BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[first_spec-Density];

     MOL::ComputeAofs( *aofs, first_spec, NUM_SPECIES, Smf, rhoYcomp,
                       AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
                       AMREX_D_DECL(*EdgeState[0],*EdgeState[1],*EdgeState[2]), first_spec, false,
                       AMREX_D_DECL(*EdgeFlux[0],*EdgeFlux[1],*EdgeFlux[2]), first_spec,
                       math_bcs, d_bcrec_ptr, geom );
     EB_set_covered(*aofs, 0.);
  }

   // Set flux, flux divergence, and face values for rho as sums of the corresponding RhoY quantities
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
   {
      const Box bx = S_mfi.tilebox();

      // this is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox& in_fab = static_cast<EBFArrayBox const&>(Smf[S_mfi]);
      const EBCellFlagFab& flags = in_fab.getEBCellFlagFab();

      if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
      {
         auto const& adv_rho   = aofs->array(S_mfi,Density);  
         amrex::ParallelFor(bx, [ adv_rho ]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            adv_rho(i,j,k) = 0.0; 
         });
         for (int dir=0; dir<AMREX_SPACEDIM; ++dir)
         {
            const Box& ebx = S_mfi.nodaltilebox(dir);
            auto const& state_ed = EdgeState[dir]->array(S_mfi,0);
            auto const& fluxes   = EdgeFlux[dir]->array(S_mfi,0);
            amrex::ParallelFor(ebx, NUM_STATE, [state_ed,fluxes]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               state_ed(i,j,k,n) = 0.0; 
               fluxes(i,j,k,n) = 0.0; 
            });
         }
      }
      else
      {
         auto const& adv_rho   = aofs->array(S_mfi,Density);  
         auto const& adv_rhoY  = aofs->array(S_mfi,first_spec);  
         amrex::ParallelFor(bx, [ adv_rho, adv_rhoY ]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            adv_rho(i,j,k) = 0.0; 
            for (int n = 0; n < NUM_SPECIES; n++) {
               adv_rho(i,j,k) += adv_rhoY(i,j,k,n);
            }
         });

         for (int dir=0; dir<AMREX_SPACEDIM; dir++)
         {
            const Box& ebx = S_mfi.grownnodaltilebox(dir,EdgeState[dir]->nGrow());
            auto const& rho    = EdgeState[dir]->array(S_mfi,Density);
            auto const& rho_F  = EdgeFlux[dir]->array(S_mfi,Density);
            auto const& rhoY   = EdgeState[dir]->array(S_mfi,first_spec);
            auto const& rhoY_F = EdgeFlux[dir]->array(S_mfi,first_spec);
            amrex::ParallelFor(ebx, [rho, rho_F, rhoY, rhoY_F]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               rho(i,j,k) = 0.0; 
               rho_F(i,j,k) = 0.0; 
               for (int n = 0; n < NUM_SPECIES; n++) {
                  rho(i,j,k) += rhoY(i,j,k,n);
                  rho_F(i,j,k) += rhoY_F(i,j,k,n);
               }
            });
         }
      }
   }

  //  Set covered values of density not to zero in roder to use fab.invert
  //  Get typical values for Rho
  {
    Vector<Real> typvals;
    typvals.resize(NUM_STATE);
    typvals[Density] = typical_values[Density];
    typvals[Temp] = typical_values[Temp];
    typvals[RhoH] = typical_values[RhoH];
    for (int k = 0; k < NUM_SPECIES; ++k) {
       typvals[first_spec+k] = typical_values[first_spec+k]*typical_values[Density];
    }
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},first_spec,NUM_SPECIES,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Temp,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Density,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},RhoH,1,typvals);
  }
  EB_set_covered_faces({AMREX_D_DECL(EdgeFlux[0],EdgeFlux[1],EdgeFlux[2])},0.);

  // Extrapolate Temp, then compute flux divergence and value for RhoH from face values of T,Y,Rho
  {
    Vector<BCRec> math_bcs(1);
    math_bcs = fetchBCArray(State_Type, Temp, 1);
    BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[Temp-Density];

    MOL::ComputeAofs( *aofs, Temp, 1, Smf, Tcomp,
		      AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
		      AMREX_D_DECL(*EdgeState[0],*EdgeState[1],*EdgeState[2]), Temp, false,
		      AMREX_D_DECL(*EdgeFlux[0],*EdgeFlux[1],*EdgeFlux[2]), Temp,
		      math_bcs, d_bcrec_ptr, geom );
    EB_set_covered(*aofs, 0.);
  }

  //  Set covered values of density not to zero in roder to use fab.invert
  //  Get typical values for Rho
  {
    Vector<Real> typvals;
    typvals.resize(NUM_STATE);
    typvals[Density] = typical_values[Density];
    typvals[Temp] = typical_values[Temp];
    typvals[RhoH] = typical_values[RhoH];
    for (int k = 0; k < NUM_SPECIES; ++k) {
       typvals[first_spec+k] = typical_values[first_spec+k]*typical_values[Density];
    }
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},first_spec,NUM_SPECIES,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Temp,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Density,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},RhoH,1,typvals);
  }

  // Compute RhoH on faces, store in NUM_SPECIES+1 component of edgestate[d]
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
     for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
     {
        Box bx = S_mfi.tilebox();
        // this is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox& in_fab = static_cast<EBFArrayBox const&>(Smf[S_mfi]);
        const EBCellFlagFab& flags = in_fab.getEBCellFlagFab();

        if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
        {
           for (int d=0; d<AMREX_SPACEDIM; ++d)
           {
              const Box& ebox = S_mfi.nodaltilebox(d);
              auto const& rhoHm   = EdgeState[d]->array(S_mfi,RhoH);
              auto const& rhoHm_F = EdgeFlux[d]->array(S_mfi,RhoH);
              amrex::ParallelFor(ebox, [rhoHm,rhoHm_F]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                 rhoHm(i,j,k) = 0.0;
                 rhoHm_F(i,j,k) = 0.0;
              });
           }
        }
        else
        {
           for (int d=0; d<AMREX_SPACEDIM; ++d)
           {
              const Box& ebox = S_mfi.grownnodaltilebox(d,EdgeState[d]->nGrow());
              auto const& rho     = EdgeState[d]->array(S_mfi,Density);  
              auto const& rhoY    = EdgeState[d]->array(S_mfi,first_spec);
              auto const& T       = EdgeState[d]->array(S_mfi,Temp);
              auto const& rhoHm   = EdgeState[d]->array(S_mfi,RhoH);

              amrex::ParallelFor(ebox, [rho, rhoY, T, rhoHm]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                 getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
              });
           }
        }
     }
  }

  //  Set covered values of density not to zero in roder to use fab.invert
  //  Get typical values for Rho
  {
    Vector<Real> typvals;
    typvals.resize(NUM_STATE);
    typvals[Density] = typical_values[Density];
    typvals[Temp] = typical_values[Temp];
    typvals[RhoH] = typical_values[RhoH];
    for (int k = 0; k < NUM_SPECIES; ++k) {
       typvals[first_spec+k] = typical_values[first_spec+k]*typical_values[Density];
    }
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},first_spec,NUM_SPECIES,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Temp,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Density,1,typvals);
    EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},RhoH,1,typvals);
  }

  for (int d=0; d<AMREX_SPACEDIM; ++d)
  {
    EdgeState[d]->FillBoundary(geom.periodicity());
    EdgeFlux[d]->FillBoundary(geom.periodicity());
  }


  // Compute -Div(flux.Area) for RhoH, return Area-scaled (extensive) fluxes
  {
     Vector<BCRec> math_bcs(1);
     math_bcs = fetchBCArray(State_Type, RhoH, 1);
     BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[RhoH-Density];

     MOL::ComputeAofs( *aofs, RhoH, 1, Smf, NUM_SPECIES+1,
		       AMREX_D_DECL(u_mac[0],u_mac[1],u_mac[2]),
		       AMREX_D_DECL(*EdgeState[0],*EdgeState[1],*EdgeState[2]), RhoH, true,
		       AMREX_D_DECL(*EdgeFlux[0],*EdgeFlux[1],*EdgeFlux[2]), RhoH,
		       math_bcs, d_bcrec_ptr, geom ); 
     EB_set_covered(*aofs, 0.);
  }

  for (int d=0; d<AMREX_SPACEDIM; ++d)
  {
     EdgeState[d]->FillBoundary(geom.periodicity());
     EdgeFlux[d]->FillBoundary(geom.periodicity());
  }

  //  Set covered values of density not to zero in roder to use fab.invert
  //  Get typical values for Rho
  {
     Vector<Real> typvals;
     typvals.resize(NUM_STATE);
     typvals[Density] = typical_values[Density];
     typvals[Temp] = typical_values[Temp];
     typvals[RhoH] = typical_values[RhoH];
     for (int k = 0; k < NUM_SPECIES; ++k) {
        typvals[first_spec+k] = typical_values[first_spec+k]*typical_values[Density];
     }
     EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},first_spec,NUM_SPECIES,typvals);
     EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Temp,1,typvals);
     EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},Density,1,typvals);
     EB_set_covered_faces({AMREX_D_DECL(EdgeState[0],EdgeState[1],EdgeState[2])},RhoH,1,typvals);
  }
  EB_set_covered_faces({AMREX_D_DECL(EdgeFlux[0],EdgeFlux[1],EdgeFlux[2])},0.);
 
#else
  //////////////////////////////////////
  //
  // HERE IS THE NON-EB PROCEDURE
  //
  //////////////////////////////////////

  // Initialize accumulation for rho = Sum(rho.Y)
  for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
    EdgeState[dir]->setVal(0.0);
    EdgeFlux[dir]->setVal(0.0);
  }

  // Advect RhoY 
{
  BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[first_spec-Density];
    
  amrex::Gpu::DeviceVector<int> iconserv;
  iconserv.resize(NUM_SPECIES, 0);
  for (int comp = 0; comp < NUM_SPECIES; ++comp)
  {
     iconserv[comp] = (advectionType[first_spec+comp] == Conservative) ? 1 : 0;
  }

  Godunov::ComputeAofs(*aofs, first_spec, NUM_SPECIES,
                        Smf, rhoYcomp,
                        AMREX_D_DECL( u_mac[0], u_mac[1], u_mac[2] ),
                        AMREX_D_DECL(*EdgeState[0], *EdgeState[1], *EdgeState[2]), first_spec, false,
                        AMREX_D_DECL(*EdgeFlux[0], *EdgeFlux[1], *EdgeFlux[2]), first_spec,
                        Force, 0, DivU, d_bcrec_ptr, geom, iconserv,
                        dt, godunov_use_ppm, godunov_use_forces_in_trans, false );
}

  // Set flux, flux divergence, and face values for rho as sums of the corresponding RhoY quantities
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
  {
    const Box bx = S_mfi.tilebox();

    auto const& adv_rho   = aofs->array(S_mfi,Density);
    auto const& adv_rhoY  = aofs->array(S_mfi,first_spec);
    amrex::ParallelFor(bx, [ adv_rho, adv_rhoY ]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      adv_rho(i,j,k) = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        adv_rho(i,j,k) += adv_rhoY(i,j,k,n);
      }
    });

    for (int dir=0; dir<AMREX_SPACEDIM; dir++)
    {
      const Box& ebx =amrex::surroundingNodes(bx,dir);
      auto const& rho    = EdgeState[dir]->array(S_mfi,Density);
      auto const& rho_F  = EdgeFlux[dir]->array(S_mfi,Density);
      auto const& rhoY   = EdgeState[dir]->array(S_mfi,first_spec);
      auto const& rhoY_F = EdgeFlux[dir]->array(S_mfi,first_spec);
      amrex::ParallelFor(ebx, [rho, rho_F, rhoY, rhoY_F]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        rho(i,j,k) = 0.0;
        rho_F(i,j,k) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          rho(i,j,k) += rhoY(i,j,k,n);
          rho_F(i,j,k) += rhoY_F(i,j,k,n);
        }
      });
    }
  }

  // Extrapolate Temp, then compute flux divergence and value for RhoH from face values of T,Y,Rho
{
  BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[Temp-Density];
  
  amrex::Gpu::DeviceVector<int> iconserv;
  iconserv.resize(1, 0);
  iconserv[0] = (advectionType[Temp] == Conservative) ? 1 : 0;


  Godunov::ComputeAofs(*aofs, Temp, 1,
                       Smf, Tcomp,
                       AMREX_D_DECL( u_mac[0], u_mac[1], u_mac[2] ),
                       AMREX_D_DECL(*EdgeState[0], *EdgeState[1], *EdgeState[2]), Temp, false,
                       AMREX_D_DECL(*EdgeFlux[0], *EdgeFlux[1], *EdgeFlux[2]), Temp,
                       Force, NUM_SPECIES, DivU, d_bcrec_ptr, geom, iconserv,
                       dt, godunov_use_ppm, godunov_use_forces_in_trans, false );

}

  // Compute RhoH on faces, store in NUM_SPECIES+1 component of edgestate[d]
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{  
  for (MFIter S_mfi(Smf,TilingIfNotGPU()); S_mfi.isValid(); ++S_mfi)
  {  
    Box bx = S_mfi.tilebox();
        
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
      const Box& ebox = amrex::surroundingNodes(bx,d);
      auto const& rho     = EdgeState[d]->array(S_mfi,Density);
      auto const& rhoY    = EdgeState[d]->array(S_mfi,first_spec);
      auto const& T       = EdgeState[d]->array(S_mfi,Temp);
      auto const& rhoHm   = EdgeState[d]->array(S_mfi,RhoH);

      amrex::ParallelFor(ebox, [rho, rhoY, T, rhoHm]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
      });
    }
  }
}

{
  BCRec  const* d_bcrec_ptr = &(m_bcrec_scalars_d.dataPtr())[RhoH-Density];

  amrex::Gpu::DeviceVector<int> iconserv;
  iconserv.resize(1, 0);
  iconserv[0] = (advectionType[RhoH] == Conservative) ? 1 : 0;

  Godunov::ComputeAofs(*aofs, RhoH, 1,
                       Smf, NUM_SPECIES+1,
                       AMREX_D_DECL( u_mac[0], u_mac[1], u_mac[2] ),
                       AMREX_D_DECL(*EdgeState[0], *EdgeState[1], *EdgeState[2]), RhoH, true,
                       AMREX_D_DECL(*EdgeFlux[0], *EdgeFlux[1], *EdgeFlux[2]), RhoH,
                       Force, NUM_SPECIES+1, DivU, d_bcrec_ptr, geom, iconserv,
                       dt, godunov_use_ppm, godunov_use_forces_in_trans, false );
}


#endif

   AMREX_D_TERM(showMF("sdc",*EdgeState[0],"sdc_ESTATE_x",level,parent->levelSteps(level));,
                showMF("sdc",*EdgeState[1],"sdc_ESTATE_y",level,parent->levelSteps(level));,
                showMF("sdc",*EdgeState[2],"sdc_ESTATE_z",level,parent->levelSteps(level)););
   AMREX_D_TERM(showMF("sdc",*EdgeFlux[0],"sdc_FLUX_x",level,parent->levelSteps(level));,
                showMF("sdc",*EdgeFlux[1],"sdc_FLUX_y",level,parent->levelSteps(level));,
                showMF("sdc",*EdgeFlux[2],"sdc_FLUX_z",level,parent->levelSteps(level)););
   showMF("sdc",*aofs,"sdc_aofs",level,parent->levelSteps(level));

// NOTE: Change sense of aofs here so that d/dt ~ aofs...be sure we use our own update function!
   aofs->mult(-1.0,Density,NUM_SCALARS);

   // AJN FLUXREG
   // We have just computed the advective fluxes.  If updateFluxReg=T,
   // ADD advective fluxes to ADVECTIVE flux register
   //
   //   FIXME: Since the fluxes and states are class data, perhaps these flux register calls should be
   //          managed in the advance function directly rather than hidden/encoded inside here??
   if (do_reflux && updateFluxReg)
   {
      if (level > 0)
      {
         for (int d = 0; d < AMREX_SPACEDIM; d++)
         {
           advflux_reg->FineAdd((*EdgeFlux[d]),d,Density,Density,NUM_SCALARS,dt);
         }
      }
      if (level < parent->finestLevel())
      {
         for (int d = 0; d < AMREX_SPACEDIM; d++)
         {
           getAdvFluxReg(level+1).CrseInit((*EdgeFlux[d]),d,Density,Density,NUM_SCALARS,-dt,FluxRegister::ADD);
         }
      }
   }

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::compute_scalar_advection_fluxes_and_divergence(): lev: "
                     << level << ", time: " << run_time << '\n';
   }
}

//
// An enum to clean up mac_sync...questionable usefullness
//
enum SYNC_SCHEME {ReAdvect, UseEdgeState, Other};

void
PeleLM::mac_sync ()
{
   BL_PROFILE("PeleLM::mac_sync()");
   if (!do_reflux) return;
   if (verbose) amrex::Print() << "... mac_sync\n";

   const Real strt_time = ParallelDescriptor::second();

   const int  finest_level   = parent->finestLevel();
   const int  ngrids         = grids.size();
   const Real prev_time      = state[State_Type].prevTime();
   const Real curr_time      = state[State_Type].curTime();
   const Real prev_pres_time = state[Press_Type].prevTime();
   const Real dt             = parent->dtLevel(level);
   MultiFab&  rh             = get_rho_half_time();

   MultiFab& S_new = get_new_data(State_Type);

   ////////////////////////
   // save states that we need to reset with each mac sync iteration
   ////////////////////////
   const int numscal = NUM_STATE - AMREX_SPACEDIM;

   MultiFab chi_sync(grids,dmap,1,0,MFInfo(),Factory());
   chi_sync.setVal(0.0);

   Vector<std::unique_ptr<MultiFab> > S_new_sav(finest_level+1);

   // Save new pre-sync S^{n+1,p} state for my level and all the finer ones
   for (int lev=level; lev<=finest_level; lev++)
   {
      const MultiFab& S_new_lev = getLevel(lev).get_new_data(State_Type);
      S_new_sav[lev].reset(new MultiFab(S_new_lev.boxArray(),
                                        S_new_lev.DistributionMap(),
                                        NUM_STATE,1,MFInfo(),getLevel(lev).Factory()));
      MultiFab::Copy(*S_new_sav[lev],S_new_lev,0,0,NUM_STATE,1);
      showMF("DBGSync",*S_new_sav[lev],"sdc_Snew_BeginSync",lev,0,parent->levelSteps(level));
   }

   Vector<std::unique_ptr<MultiFab> > Ssync_sav(finest_level);
   Vector<std::unique_ptr<MultiFab> > Vsync_sav(finest_level);

   // Save Sync RHS & Vsync for my level and all the finer ones (but the finest)
   for (int lev=level; lev<=finest_level-1; lev++)
   {
      const MultiFab& Ssync_lev = getLevel(lev).Ssync;
      Ssync_sav[lev].reset(new MultiFab(Ssync_lev.boxArray(),
                                        Ssync_lev.DistributionMap(),
                                        numscal,1,MFInfo(),getLevel(lev).Factory()));
      MultiFab::Copy(*Ssync_sav[lev],Ssync_lev,0,0,numscal,1);
      showMF("DBGSync",*Ssync_sav[lev],"sdc_Ssync_BeginSync",level,0,parent->levelSteps(level));

      const MultiFab& Vsync_lev = getLevel(lev).Vsync;
      Vsync_sav[lev].reset(new MultiFab(Vsync_lev.boxArray(),
                                        Vsync_lev.DistributionMap(),
                                        AMREX_SPACEDIM,1,MFInfo(),getLevel(lev).Factory()));
      MultiFab::Copy(*Vsync_sav[lev],Vsync_lev,0,0,AMREX_SPACEDIM,1);
      showMF("DBGSync",*Vsync_sav[lev],"sdc_Vsync_BeginSync",level,0,parent->levelSteps(level));
   }

   if (use_wbar) {
      // compute (overwrite) beta grad Wbar terms using the pre-sync state
      compute_Wbar_fluxes(curr_time,0,1.0);
   }

   ////////////////////////
   // begin mac_sync_iter loop here
   // The loop allows an update of chi to ensure that the sync correction remains
   // on the constrain.
   ////////////////////////

   MultiFab DeltaYsync(grids,dmap,NUM_SPECIES,0,MFInfo(),Factory());
   MultiFab chi_sync_increment(grids,dmap,1,0,MFInfo(),Factory());

   // save pressure
   Real p_amb_new_temp = p_amb_new;

   BL_PROFILE_VAR("PeleLM::mac_sync::Ssync", PLM_SSYNC);
   for (int mac_sync_iter=0; mac_sync_iter < num_mac_sync_iter; mac_sync_iter++)
   {
      bool last_mac_sync_iter = (mac_sync_iter == num_mac_sync_iter-1);

      if ( mac_sync_iter != 0 ) {
         // Restore saved copy of S_new^{n+1,p} state into S_new
         for (int lev=level; lev<=finest_level; lev++)
         {
           MultiFab& S_new_lev = getLevel(lev).get_new_data(State_Type);
           MultiFab::Copy(S_new_lev,*S_new_sav[lev],0,0,NUM_STATE,1);
         }

         // Restore stored Ssync/Vsync from the one saved before sync_iter
         for (int lev=level; lev<=finest_level-1; lev++)
         {
           MultiFab& Ssync_lev = getLevel(lev).Ssync;
           MultiFab::Copy(Ssync_lev,*Ssync_sav[lev],0,0,numscal,1);

           MultiFab& Vsync_lev = getLevel(lev).Vsync;
           MultiFab::Copy(Vsync_lev,*Vsync_sav[lev],0,0,AMREX_SPACEDIM,1);
         }
      }

      // Update back rho^{n+1,p} and rho^{n+1/2,p}
      make_rho_curr_time();
      get_rho_half_time();

      //
      // Compute the corrective pressure, mac_sync_phi, used to 
      // compute U^{ADV,corr} in mac_sync_compute
      //
      //TODO: offset/subtract_avg is not used ... ?
//      bool subtract_avg = (closed_chamber && level == 0);
      Real offset = 0.0;

      BL_PROFILE_VAR("PeleLM::mac_sync::ucorr", PLM_UCORR);
      Array<MultiFab*,AMREX_SPACEDIM> Ucorr;
#ifdef AMREX_USE_EB
      const int ng = 4; // For redistribution ... We may not need 4 but for now we play safe
#else
      const int ng = 0;
#endif
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
         const BoxArray& edgeba = getEdgeBoxArray(idim);
         // fixme? unsure how many ghost cells...
         Ucorr[idim]= new MultiFab(edgeba,dmap,1,ng,MFInfo(),Factory());
      }
      Vector<BCRec> rho_math_bc = fetchBCArray(State_Type,Density,1);
      mac_projector->mac_sync_solve(level,dt,rh,rho_math_bc[0],fine_ratio,Ucorr,&chi_sync);
      showMF("DBGSync",chi_sync,"sdc_chi_sync_inSync",level,mac_sync_iter,parent->levelSteps(level));
      showMF("DBGSync",*Ucorr[0],"sdc_UcorrX_inSync",level,mac_sync_iter,parent->levelSteps(level));
      showMF("DBGSync",*Ucorr[1],"sdc_UcorrY_inSync",level,mac_sync_iter,parent->levelSteps(level));
      BL_PROFILE_VAR_STOP(PLM_UCORR);


      if (closed_chamber && level == 0)
      {
        Print() << "0-1 MAC SYNC OFFSET " << offset << std::endl;
        Print() << "0-1 MAC SYNC PAMB ADJUSTMENT " << -dt*(offset/thetabar) << std::endl;
        p_amb_new = p_amb_new_temp - dt*(offset/thetabar);
        p_amb_old = p_amb_new;
      }

      Vector<SYNC_SCHEME> sync_scheme(NUM_STATE,ReAdvect);

      if (do_mom_diff == 1) {
         for (int i=0; i<AMREX_SPACEDIM; ++i) {
            sync_scheme[i] = UseEdgeState;
         }
      }

      for (int i=AMREX_SPACEDIM; i<NUM_STATE; ++i) {
         sync_scheme[i] = UseEdgeState;
      }

      Vector<int> incr_sync(NUM_STATE,0);
      for (int i=0; i<sync_scheme.size(); ++i) {
         if (sync_scheme[i] == ReAdvect) {
            incr_sync[i] = 1;
         }
      }
 
      // After solving for mac_sync_phi in mac_sync_solve(), we
      // can now do the sync advect step in mac_sync_compute().
      // This consists of two steps
      //
      // 1. compute U^{ADV,corr} as the gradient of mac_sync_phi : done above now !
      // 2. add -D^MAC ( U^{ADV,corr} * rho * q)^{n+1/2} ) to flux registers,
      //    which already contain the delta F adv/diff flux mismatches

      BL_PROFILE_VAR("PeleLM::mac_sync::Vsync", PLM_VSYNC);
      if (do_mom_diff == 0) 
      {
         mac_projector->mac_sync_compute(level,Ucorr,u_mac,Vsync,Ssync,rho_half,
                                         (level > 0) ? &getAdvFluxReg(level) : 0,
                                         advectionType,prev_time,
                                         prev_pres_time,dt,NUM_STATE,
                                         be_cn_theta,
                                         modify_reflux_normal_vel,
                                         do_mom_diff,
                                         incr_sync,
                                         last_mac_sync_iter);
      }
      else
      {
         for (int comp=0; comp<AMREX_SPACEDIM; ++comp)
         {
            if (sync_scheme[comp]==UseEdgeState)
            {
               mac_projector->mac_sync_compute(level,Ucorr,Vsync,comp,
                                               comp,EdgeState, comp,rho_half,
                                               (level > 0 ? &getAdvFluxReg(level):0),
                                               advectionType,modify_reflux_normal_vel,dt,
                                               last_mac_sync_iter);
            }
         }
      }
      showMF("DBGSync",Vsync,"sdc_Vsync_AfterMACSync",level,mac_sync_iter,parent->levelSteps(level));
      BL_PROFILE_VAR_STOP(PLM_VSYNC);

      //
      // Scalars.
      //
      showMF("DBGSync",Ssync,"sdc_Ssync_BefMinusUcorr",level,mac_sync_iter,parent->levelSteps(level));
      for (int comp=AMREX_SPACEDIM; comp<NUM_STATE; ++comp)
      {
         if (sync_scheme[comp]==UseEdgeState)
         {
            int s_ind = comp - AMREX_SPACEDIM;
            //
            // Ssync contains the adv/diff coarse-fine flux mismatch divergence
            // This routine does a sync advect step for a single scalar component,
            // i.e., subtracts the D(Ucorr rho q) term from Ssync
            // The half-time edge states are passed in.
            // This routine is useful when the edge states are computed
            // in a physics-class-specific manner. (For example, as they are
            // in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
            // Note: the density component now contains (delta rho)^sync since there
            // is no diffusion for this term
            //
            mac_projector->mac_sync_compute(level,Ucorr,Ssync,comp,s_ind,
                                            EdgeState,comp,rho_half,
                                            (level > 0 ? &getAdvFluxReg(level):0),
                                            advectionType,modify_reflux_normal_vel,dt,
                                            last_mac_sync_iter);
         }
      }
      showMF("DBGSync",Ssync,"sdc_Ssync_MinusUcorr",level,mac_sync_iter,parent->levelSteps(level));

      //
      // Delete Ucorr; we're done with it.
      //
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
         delete Ucorr[idim];

      Ssync.mult(dt); // Turn this into an increment over dt

      MultiFab DdWbar;
      if (use_wbar) {
         // Substract pre-sync beta grad Wbar fluxes
         // If mac_sync_iter = 0 -> dWbar fluxes == 0, otherwise post-sync - pre-sync
         compute_Wbar_fluxes(curr_time,1,-1.0);

         // take divergence of beta grad delta Wbar
         DdWbar.define(grids,dmap,NUM_SPECIES,nGrowAdvForcing,MFInfo(),Factory());
         MultiFab* const * fluxWbar = SpecDiffusionFluxWbar;
         flux_divergence(DdWbar,0,fluxWbar,0,NUM_SPECIES,-1);
         DdWbar.mult(dt);
      }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho      = S_new.array(mfi,Density);
         auto const& rhoY     = S_new.array(mfi,first_spec);
         auto const& dYsync   = DeltaYsync.array(mfi);
         auto const& drhosync = Ssync.array(mfi,Density-AMREX_SPACEDIM);
         auto const& ssync    = Ssync.array(mfi,first_spec-AMREX_SPACEDIM);
         int use_wbar_lcl     = use_wbar;
         auto const& dwbar    = (use_wbar) ? DdWbar.array(mfi) : DeltaYsync.array(mfi);

         amrex::ParallelFor(bx, [rho, rhoY, dYsync, drhosync, dwbar,
                                 ssync, use_wbar_lcl]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            Real rhoinv = 1.0/rho(i,j,k);
            for (int n = 0; n < NUM_SPECIES; n++) {
               dYsync(i,j,k,n) = rhoY(i,j,k,n) * rhoinv * drhosync(i,j,k);   // DeltaYSsync = Y^{n+1,p} * (delta rho)^sync
               ssync(i,j,k,n) -= dYsync(i,j,k,n);                            // Ssync = Ssync - DeltaYSync
               if (use_wbar_lcl) {
                  ssync(i,j,k,n) += dwbar(i,j,k,n);                          // Ssync = Ssync + DdWbar
               }
            }
         });
      }

      // TODO: the folloxing should be merged with the above kernel
      // but I need to find a way to EB_set_covered.
      // 
      // Increment density, rho^{n+1} = rho^{n+1,p} + (delta_rho)^sync
      //
#ifdef AMREX_USE_EB
      EB_set_covered(Ssync,0.);
#endif
      MultiFab::Add(S_new,Ssync,Density-AMREX_SPACEDIM,Density,1,0);
      make_rho_curr_time();
      BL_PROFILE_VAR_STOP(PLM_SSYNC);

      //
      // If mom_diff, scale Vsync by rho so we can diffuse with the same call below
      //
      BL_PROFILE_VAR_START(PLM_VSYNC);
      if (do_mom_diff == 1)
      {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(rho_ctime, TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& rho_c    = rho_ctime.array(mfi);
            auto const& vsync    = Vsync.array(mfi,Xvel);
            amrex::ParallelFor(bx, [rho_c, vsync]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               Real rhoinv = 1.0/rho_c(i,j,k);
               for (int n = 0; n < AMREX_SPACEDIM; n++) {
                  vsync(i,j,k,n) *= rhoinv;
               }
            });
         }
      }
      BL_PROFILE_VAR_STOP(PLM_VSYNC);

      if (do_diffuse_sync)
      {
         FluxBoxes fb_beta(this);
         MultiFab** beta = fb_beta.get();

         BL_PROFILE_VAR_STOP(PLM_VSYNC);
         if (is_diffusive[Xvel])
         {
           int rho_flag = (do_mom_diff == 0) ? 1 : 3;
           getViscosity(beta, curr_time);
           diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,rho_half,rho_flag,beta,0,
                                    last_mac_sync_iter);
         }
         BL_PROFILE_VAR_STOP(PLM_VSYNC);

         // species and temperature diffusion

         BL_PROFILE_VAR_START(PLM_SSYNC);
         // FIXME: Really wish there was a way to avoid this temporary....
         FluxBoxes fb_GammaKp1(this, NUM_SPECIES+3, 0);
         MultiFab** GammaKp1 = fb_GammaKp1.get();
         for (int d=0; d<AMREX_SPACEDIM; ++d) {
           MultiFab::Copy(*GammaKp1[d],*SpecDiffusionFluxnp1[d],0,0,NUM_SPECIES+3,0); // get Gamma^{presync}
         }

         FluxBoxes fb_betanp1(this, NUM_SPECIES+1, 0);
         MultiFab **betanp1 = fb_betanp1.get();
         getDiffusivity(betanp1, curr_time, first_spec, 0, NUM_SPECIES);  // species
         getDiffusivity(betanp1, curr_time, Temp, NUM_SPECIES, 1);        // temperature (lambda)
         compute_enthalpy_fluxes(GammaKp1,betanp1,curr_time);             // Compute F[N+1] = sum_m (H_m Gamma_m), 
                                                                          //         F[N+2] = - lambda grad T

         // Note: DT (comp 1) and DD (comp 0) are in the same multifab
         MultiFab DiffTerms_pre(grids,dmap,2,0,MFInfo(),Factory());
         flux_divergence(DiffTerms_pre,0,GammaKp1,NUM_SPECIES+1,2,-1);

         // Diffuse species sync to get dGamma^{sync}, then form Gamma^{postsync} = Gamma^{presync} + dGamma^{sync}
         // Before this call Ssync = \nabla \cdot \delta F_{rhoY} - \delta_ADV - Y^{n+1,p}\delta rho^{sync}
         bool include_deltaWbar = true;
         differential_spec_diffuse_sync(dt, include_deltaWbar, last_mac_sync_iter);
         // Ssync for species now contains rho^{n+1}(delta_Y)^sync

         for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Add(*SpecDiffusionFluxnp1[d],*GammaKp1[d],0,0,NUM_SPECIES,0);
         }

         //
         // For all species increment sync by (sync_for_rho)*Y_presync.
         // Before this, Ssync holds rho^{n+1} (delta Y)^sync
         // DeltaYsync holds Y^{n+1,p} * (delta rho)^sync
         //
         MultiFab::Add(Ssync,DeltaYsync,0,first_spec-AMREX_SPACEDIM,NUM_SPECIES,0);
         MultiFab::Add(S_new,Ssync,first_spec-AMREX_SPACEDIM,first_spec,NUM_SPECIES,0);

         if (use_wbar) {
            // compute (overwrite) beta grad Wbar terms using the post-sync state
            compute_Wbar_fluxes(curr_time,0,1.0);
         }

         // Trying to solve for:
         // \rho^{n+1} * Cp{n+1,\eta} * ∆T^{\eta+1} - dt / 2 \nabla \cdot \lambda^{n+1,p} \nabla ∆T^{\eta+1} = 
         // \rho^{n+1,p}*h^{n+1,p} - \rho^{n+1}*h^{n+1,\eta} + dt*Ssync + dt/2*(DT^{n+1,\eta} - DT^{n+1,p} + H^{n+1,\eta} - H^{n+1,p})

         // Here Ssync contains refluxed enthalpy fluxes (from -lambda.Grad(T) and hm.Gamma_m)
         MultiFab Trhs(grids,dmap,1,0,MFInfo(),Factory());
         MultiFab Told(grids,dmap,1,0,MFInfo(),Factory());
         MultiFab RhoCp_post(grids,dmap,1,0,MFInfo(),Factory());
         MultiFab DiffTerms_post(grids,dmap,2,0,MFInfo(),Factory());

         if (deltaT_verbose) Print() << "Starting deltaT iters in mac_sync... " << std::endl;

         Real deltaT_iter_norm = 0;
         for (int L=0; L<num_deltaT_iters_MAX && (L==0 || deltaT_iter_norm >= deltaT_norm_max); ++L)
         {
            // Bundle in a single kernel assembling the Trhs and compute rhoCpmix
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& bx = mfi.tilebox();
               // Things needed for the Trhs
               auto const& rhs      = Trhs.array(mfi);
               auto const& dT_pre   = DiffTerms_pre.array(mfi,1);
               auto const& dd_pre   = DiffTerms_pre.array(mfi,0);
               auto const& dT_post  = DiffTerms_post.array(mfi,1);
               auto const& dd_post  = DiffTerms_post.array(mfi,0);
               auto const& rhoH_eta = S_new.array(mfi,RhoH);
               auto const& rhoH_0   = S_new_sav[level]->array(mfi,RhoH);
               auto const& ssync    = Ssync.array(mfi,RhoH-AMREX_SPACEDIM);

               // Things needed to compute RhoCpmix
               auto const& rho     = S_new.array(mfi,Density);
               auto const& rhoY    = S_new.array(mfi,first_spec);
               auto const& T       = S_new.array(mfi,Temp);
               auto const& RhoCpm  = RhoCp_post.array(mfi);

               // Things needed to store T
               auto const& T_old   = Told.array(mfi);

               amrex::ParallelFor(bx, [rhs, dT_pre, dd_pre, dT_post, dd_post, rhoH_eta, rhoH_0, ssync, dt, L,
                                       rho, rhoY, T, RhoCpm,
                                       T_old]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  // Trhs computation
                  rhs(i,j,k) =  rhoH_0(i,j,k) - rhoH_eta(i,j,k) + ssync(i,j,k);
                  if ( L > 0 ) {   // On first pass, pre & post are equals. Do that to avoid a copy.
                     rhs(i,j,k) += dt * ( dT_post(i,j,k) + dd_post(i,j,k) - dT_pre(i,j,k) - dd_pre(i,j,k) );
                  }

                  // rhoCpmix computation
                  getCpmixGivenRYT( i, j, k, rho, rhoY, T, RhoCpm );
                  RhoCpm(i,j,k) *= rho(i,j,k);

                  // Store Told
                  T_old(i,j,k) = T(i,j,k);
               });
            }

            showMF("DBGSync",Trhs,"sdc_Trhs_indeltaTiter",level,(mac_sync_iter+1)*1000+L,parent->levelSteps(level));
            showMF("DBGSync",RhoCp_post,"sdc_RhoCp_post_indeltaTiter",level,(mac_sync_iter+1)*1000+L,parent->levelSteps(level));

            // Let's deltaT be an alias to T in S_new
            get_new_data(State_Type).setVal(0.0,Temp,1,1);        // Can't do that in the kernel above, I need 1 GC here.
            MultiFab deltaT(get_new_data(State_Type), amrex::make_alias, Temp, 1);

            const Real be_cn_theta_SDC = 1.0;

            // Set-up our own MG solver for the deltaT
            LPInfo info;
            info.setAgglomeration(1);
            info.setConsolidation(1);
            info.setMetricTerm(false);
#ifdef AMREX_USE_EB
            const auto& ebf = &dynamic_cast<EBFArrayBoxFactory const&>((parent->getLevel(level)).Factory());
            MLEBABecLap deltaTSyncOp({geom}, {grids}, {dmap}, info, {ebf});
#else
            MLABecLaplacian deltaTSyncOp({geom}, {grids}, {dmap}, info);
#endif
            deltaTSyncOp.setMaxOrder(diffusion->maxOrder());
            deltaTSyncOp.setScalars(1.0, dt*be_cn_theta_SDC);
            deltaTSyncOp.setACoeffs(0.0, RhoCp_post);

            AMREX_D_TERM(MultiFab lambdax(*betanp1[0], amrex::make_alias, NUM_SPECIES, 1);,
                         MultiFab lambday(*betanp1[1], amrex::make_alias, NUM_SPECIES, 1);,
                         MultiFab lambdaz(*betanp1[2], amrex::make_alias, NUM_SPECIES, 1););
            std::array<MultiFab*,AMREX_SPACEDIM> bcoeffs{AMREX_D_DECL(&lambdax,&lambday,&lambdaz)};
            deltaTSyncOp.setBCoeffs(0,amrex::GetArrOfConstPtrs(bcoeffs));

            std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
            std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
            const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();
            const BCRec& bc = theBCs[Temp];
            Diffusion::setDomainBC(mlmg_lobc, mlmg_hibc, bc); // Same for all comps, by assumption
            deltaTSyncOp.setDomainBC(mlmg_lobc, mlmg_hibc);
            if (level > 0) {
               deltaTSyncOp.setCoarseFineBC(nullptr, crse_ratio[0]);
            }
            deltaTSyncOp.setLevelBC(0, &deltaT);

            MLMG deltaTSyncSolve(deltaTSyncOp);
            deltaTSyncSolve.setVerbose(0);

            // Ignore EB covered cells when computing Trhs.norm. In certain cases, including
            // covered cells can result in a very large S_tol_abs that tricks MLMG into
            // thinking the system is sufficiently solved when it's not
            const Real S_tol     = visc_tol;
            const Real S_tol_abs = visc_tol * Trhs.norm0(0,0,false,true);

            deltaTSyncSolve.solve({&deltaT}, {&Trhs}, S_tol, S_tol_abs); 

            deltaT_iter_norm = deltaT.norm0(0);
            if (deltaT_verbose) {
               Print() << "DeltaTsync solve norm [" << L << "] = " << deltaT_iter_norm << std::endl;
            }

            showMF("DBGSync",deltaT,"sdc_deltaT_indeltaTiter",level,(mac_sync_iter+1)*1000+L,parent->levelSteps(level));

            MultiFab::Add(get_new_data(State_Type),Told,0,Temp,1,0);
            compute_enthalpy_fluxes(SpecDiffusionFluxnp1,betanp1,curr_time); // Compute F[N+1], F[N+2]
            flux_divergence(DiffTerms_post,0,SpecDiffusionFluxnp1,NUM_SPECIES+1,2,-1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            // Update (RhoH)^{postsync,L}
            for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& bx = mfi.tilebox();
               auto const& rho     = S_new.array(mfi,Density);  
               auto const& rhoY    = S_new.array(mfi,first_spec);
               auto const& T       = S_new.array(mfi,Temp);
               auto const& rhoHm   = S_new.array(mfi,RhoH);

               amrex::ParallelFor(bx, [rho, rhoY, T, rhoHm]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
               });
            }

            if ( L==(num_deltaT_iters_MAX-1) 
                 && deltaT_iter_norm >= deltaT_norm_max
                 && deltaT_crashOnConvFail ) {
               Abort("deltaT_iters not converged in mac_sync");
            }
         } // deltaT_iters

         // Construct the \delta(\rho h) from (post - pre) deltaT iters
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& rhoHm   = S_new.array(mfi,RhoH);
            auto const& rhoHm_0 = S_new_sav[level]->array(mfi,RhoH);
            auto const& ssync   = Ssync.array(mfi,RhoH-AMREX_SPACEDIM);
            amrex::ParallelFor(bx, [ rhoHm, rhoHm_0, ssync ]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               ssync(i,j,k) = rhoHm(i,j,k) - rhoHm_0(i,j,k);
            });
         }

         // Need to compute increment of enthalpy fluxes if last_sync_iter
         // Enthalpy fluxes for pre-sync state stored in GammaKp1 above 
         // Enthalpy fluxes for post-sync computed in the last deltaT iter and stored in SpecDiffusionFluxnp1
         // To avoid creating yet another container, subtract SpecDiffusionFluxnp1 from GammaKp1 for the two enthalpy pieces
         // NOTE: FluxReg incremented with -dt scale whereas it should be dt within a SDC context, but what we have
         // in GammaKp1 is currently of the opposite sign.
         if (do_reflux && level > 0 && last_mac_sync_iter) {
            for (int d=0; d<AMREX_SPACEDIM; ++d)
            {
               MultiFab::Subtract(*GammaKp1[d],*SpecDiffusionFluxnp1[d],NUM_SPECIES+1,NUM_SPECIES+1,2,0);
               getViscFluxReg().FineAdd(*GammaKp1[d],d,NUM_SPECIES+2,RhoH,1,-dt);
               getAdvFluxReg().FineAdd(*GammaKp1[d],d,NUM_SPECIES+1,RhoH,1,-dt);
            }
         }
         BL_PROFILE_VAR_STOP(PLM_SSYNC);
      }
      else
      {
         Abort("FIXME: Properly deal with do_diffuse_sync=0");
      }

      BL_PROFILE_VAR_START(PLM_SSYNC);
      // Update coarse post-sync temp from rhoH and rhoY
      RhoH_to_Temp(S_new);
      setThermoPress(curr_time);

      showMF("DBGSync",Ssync,"sdc_SsyncToInterp_inSyncIter",level,mac_sync_iter,parent->levelSteps(level));
      //
      // Interpolate the sync correction to the finer levels.
      // Interpolate rhoY and rhoH Syncs and recompute rho and Temp on fine
      //
      Real mult = 1.0;
      Vector<int*>         sync_bc(grids.size());
      Vector< Vector<int> > sync_bc_array(grids.size());
      for (int i = 0; i < ngrids; i++)
      {
        sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
        sync_bc[i]       = sync_bc_array[i].dataPtr();
      }
      IntVect ratio = IntVect::TheUnitVector();

      for (int lev = level+1; lev <= finest_level; lev++)
      {
         ratio               *= parent->refRatio(lev-1);
         PeleLM& fine_level  = getLevel(lev);
         MultiFab& S_new_lev = fine_level.get_new_data(State_Type);
         showMF("DBGSync",S_new_lev,"sdc_SnewBefIncr_inSync",level,mac_sync_iter,parent->levelSteps(level));
         //
         // New way of interpolating syncs to make sure mass is conserved
         // and to ensure freestream preservation for species & temperature.
         //
         const BoxArray& fine_grids            = S_new_lev.boxArray();
         const DistributionMapping& fine_dmap  = S_new_lev.DistributionMap();
         const int nghost                      = S_new_lev.nGrow();

         MultiFab increment(fine_grids, fine_dmap, numscal, nghost,MFInfo(),getLevel(lev).Factory());
         increment.setVal(0.0,nghost);

         SyncInterp(Ssync, level, increment, lev, ratio, 
                    first_spec-AMREX_SPACEDIM, first_spec-AMREX_SPACEDIM, NUM_SPECIES+1, 1, mult,
                    sync_bc.dataPtr());

         // Set rhoIncr = Sum rhoYincr if required and update S_new += incr
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(increment,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& rhoincr     = increment.array(mfi,Density-AMREX_SPACEDIM);
            auto const& rhoYincr    = increment.array(mfi,first_spec-AMREX_SPACEDIM);
            auto const& S_new_incr  = S_new_lev.array(mfi,Density);
            int  rhoSumrhoY_flag    = do_set_rho_to_species_sum;

            amrex::ParallelFor(bx, [rhoincr, rhoYincr, S_new_incr, rhoSumrhoY_flag, numscal]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               if ( rhoSumrhoY_flag ) {
                  rhoincr(i,j,k) = 0.0;   
                  for (int n = 0; n < NUM_SPECIES; n++) { 
                     rhoincr(i,j,k) += rhoYincr(i,j,k,n);
                  }
               }   
               for (int n = 0; n < numscal; n++) { 
                  S_new_incr(i,j,k,n) += rhoincr(i,j,k,n);
               }   
            });
         }
         showMF("DBGSync",increment,"sdc_increment_inSync",level,mac_sync_iter,parent->levelSteps(level));

         if (last_mac_sync_iter)
         {
           fine_level.make_rho_curr_time();
           fine_level.incrRhoAvg(increment,Density-AMREX_SPACEDIM,1.0);
         }
         //
         // Recompute temperature and rho R T after interpolation of the mac_sync correction
         //   of the individual quantities rho, Y, T.
         //
         RhoH_to_Temp(S_new_lev);
         fine_level.setThermoPress(curr_time);
      }

      //
      // Average down rho R T after interpolation of the mac_sync correction
      //   of the individual quantities rho, Y, T.
      //
      for (int lev = finest_level-1; lev >= level; lev--)
      {
        PeleLM& fine_level = getLevel(lev+1);
        PeleLM& crse_level = getLevel(lev);

        MultiFab& S_crse_loc = crse_level.get_new_data(State_Type);
        MultiFab& S_fine_loc = fine_level.get_new_data(State_Type);

	crse_level.average_down(S_fine_loc, S_crse_loc, RhoRT, 1);
      }
      PeleLM& fine_level = getLevel(level+1);
      showMF("DBGSync",fine_level.get_new_data(State_Type),"sdc_SnewFine_EndSyncIter",level+1,mac_sync_iter,parent->levelSteps(level));
      showMF("DBGSync",get_new_data(State_Type),"sdc_SnewCoarse_EndSyncIter",level,mac_sync_iter,parent->levelSteps(level));

      // Compute chi increment and update chi_sync
      chi_sync_increment.setVal(0,0);
      calc_dpdt(curr_time,dt,chi_sync_increment,u_mac);
      MultiFab::Subtract(chi_sync,chi_sync_increment,0,0,1,0);
      BL_PROFILE_VAR_STOP(PLM_SSYNC);
   } // end loop over mac_sync_iters

   BL_PROFILE_VAR_START(PLM_SSYNC);
   DeltaYsync.clear();
   chi_sync.clear();
   S_new_sav.clear();
   Ssync_sav.clear();
   Vsync_sav.clear();

   //
   // Diffuse all other state quantities
   //

   // Start by getting a list of components to do
   // Do not diffuse {Density, RhoH, Spec, Temp}
   //         or anything that is not diffusive
   Vector<int> comps_to_diffuse;
   for (int sigma = 0; sigma < numscal; sigma++)
   {
      const int state_ind = AMREX_SPACEDIM + sigma;
      const bool is_spec = state_ind<=last_spec && state_ind>=first_spec;
      int do_it
            =  state_ind!=Density 
            && state_ind!=Temp
            && state_ind!=RhoH
            && !is_spec
            && is_diffusive[state_ind];

      if (do_it)
      {
         comps_to_diffuse.push_back(sigma);
      }
   }

   if (comps_to_diffuse.size() > 0)
   {
      amrex::Print() << " Doing some comps_to_diffuse \n"; 
      FluxBoxes fb_flux(this);
      MultiFab **flux = fb_flux.get();

      FluxBoxes fb_beta(this);
      MultiFab** beta = fb_beta.get();

      MultiFab DeltaSsync(grids,dmap,1,0);

      for (int n=0; n<comps_to_diffuse.size(); ++n)
      {
         int sigma = comps_to_diffuse[n];
         const int state_ind = AMREX_SPACEDIM + sigma;

         //
         // For variables of the form rhoQ, increment RHS
         // by -DeltaSsync = - (sync_for_rho)*Q_presync.
         //
         if (advectionType[state_ind] == Conservative)
         {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& bx = mfi.tilebox();
               auto const& rho      = S_new.array(mfi,Density);
               auto const& scalar   = S_new.array(mfi,state_ind);
               auto const& dssync   = DeltaSsync.array(mfi);
               auto const& drhosync = Ssync.array(mfi,Density-AMREX_SPACEDIM);
               auto const& scalsync = Ssync.array(mfi,sigma);
               amrex::ParallelFor(bx, [rho, scalar, dssync, drhosync, scalsync ]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  dssync(i,j,k) = scalar(i,j,k) / rho(i,j,k) * drhosync(i,j,k);
                  scalsync(i,j,k) -= dssync(i,j,k);
               });
            }
         }   

         MultiFab* alpha = 0;
         getDiffusivity(beta, curr_time, state_ind, 0, 1);
         int rho_flag = 0;

         // on entry, Ssync = RHS for sync_for_Q diffusive solve
         // on exit,  Ssync = rho^{n+1} * sync_for_Q
         // on exit,  flux = - coeff * grad Q
         diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,rho_half,
                                  rho_flag,flux,0,beta,0,alpha,0);

         // rho Scal{n+1} = rho Scal{n+1,p} + rho^{n+1} * delta{scal}^{sync} + delta{rho}^{sync} * scal{n+1,p}
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(S_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& scalar   = S_new.array(mfi,state_ind);
            auto const& dssync   = DeltaSsync.array(mfi);
            auto const& scalsync = Ssync.array(mfi,sigma);
            amrex::ParallelFor(bx, [scalar, dssync, scalsync ]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               scalar(i,j,k) += scalsync(i,j,k) + dssync(i,j,k);
            });
         }

#ifdef AMREX_USE_EB
         EB_set_covered_faces({AMREX_D_DECL(flux[0],flux[1],flux[2])},0.);
#endif

         if (level > 0)
         {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
               getViscFluxReg().FineAdd(*flux[d],d,0,state_ind,1,dt);
            }
         }

         // Interpolate sync to finer levels
         Vector<int*>         sync_bc(grids.size());
         Vector< Vector<int> > sync_bc_array(grids.size());
         for (int i = 0; i < ngrids; i++)
         {
            sync_bc_array[i] = getBCArray(State_Type,i,state_ind,1);
            sync_bc[i]       = sync_bc_array[i].dataPtr();
         }

         IntVect ratio = IntVect::TheUnitVector();
         for (int lev = level+1; lev <= finest_level; lev++)
         {
            ratio                                *= parent->refRatio(lev-1);
            PeleLM& fine_level                   =  getLevel(lev);
            MultiFab& S_new_lev                  =  fine_level.get_new_data(State_Type);
            const BoxArray& fine_grids           =  S_new_lev.boxArray();
            const DistributionMapping& fine_dmap =  S_new_lev.DistributionMap();
            const int nghost                     =  S_new_lev.nGrow();

            MultiFab increment(fine_grids, fine_dmap, numscal, nghost);
            increment.setVal(0,nghost);

            Real mult = 1.0;
            SyncInterp(Ssync, level, increment, lev, ratio, 
                       sigma, sigma, 1, 1, mult, sync_bc.dataPtr());

            MultiFab::Add(S_new_lev,increment,0,state_ind,1,nghost);
         }
      }
   }

   //
   // Average down rho R T after interpolation of the mac_sync correction
   //   of the individual quantities rho, Y, T.
   //
   for (int lev = finest_level-1; lev >= level; lev--)
   {
     PeleLM& fine_level = getLevel(lev+1);
     PeleLM& crse_level = getLevel(lev);

     MultiFab& S_crse_loc = crse_level.get_new_data(State_Type);
     MultiFab& S_fine_loc = fine_level.get_new_data(State_Type);

     crse_level.average_down(S_fine_loc, S_crse_loc, RhoRT, 1);
   }
   BL_PROFILE_VAR_STOP(PLM_SSYNC);

   if (verbose)
   {
     const int IOProc   = ParallelDescriptor::IOProcessorNumber();
     Real      run_time = ParallelDescriptor::second() - strt_time;

     ParallelDescriptor::ReduceRealMax(run_time,IOProc);

     amrex::Print() << "PeleLM::mac_sync(): lev: " << level << ", time: " << run_time << '\n';
   }

}

void
PeleLM::compute_Wbar_fluxes(Real time,
                            int increment_flag,
                            Real increment_coeff)
{
   BL_PROFILE("PLM::compute_Wbar_fluxes()");

   // average species transport coefficients to edged
   FluxBoxes fb_beta(this, NUM_SPECIES, 0);
   MultiFab **beta = fb_beta.get();
   getDiffusivity(beta, time, first_spec, 0, NUM_SPECIES); // species

   int nGrowOp = 1;

   // Define Wbar
   MultiFab Wbar;
   Wbar.define(grids,dmap,1,nGrowOp,MFInfo(),Factory());

   // Get fillpatched rho and rhoY
   FillPatchIterator fpi(*this,Wbar,nGrowOp,time,State_Type,Density,NUM_SPECIES+1);
   MultiFab& mf= fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& gbx = mfi.growntilebox(); 
      auto const& rho     = mf.const_array(mfi,0);
      auto const& rhoY    = mf.const_array(mfi,1);
      auto const& Wbar_ar = Wbar.array(mfi); 
      amrex::ParallelFor(gbx, [rho, rhoY, Wbar_ar]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         getMwmixGivenRY(i, j, k, rho, rhoY, Wbar_ar);
      });
   }

   //
   // Use a LinOp + MLMG to get the gradient of Wbar
   //
   LPInfo info;
   info.setAgglomeration(0);
   info.setConsolidation(0);
   info.setMetricTerm(false);
   info.setMaxCoarseningLevel(0);

#ifdef AMREX_USE_EB
   const auto& ebf = &(dynamic_cast<EBFArrayBoxFactory const&>(Factory()));
   MLEBABecLap visc_op({geom}, {grids}, {dmap}, info, {ebf});
#else
   MLABecLaplacian visc_op({geom}, {grids}, {dmap}, info);
#endif

   visc_op.setMaxOrder(diffusion->maxOrder());

   // Set the domain BC as one of the species
   {
     const Vector<BCRec>& theBCs = AmrLevel::desc_lst[State_Type].getBCs();
     const BCRec& bc = theBCs[first_spec];

     std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
     std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
     Diffusion::setDomainBC(mlmg_lobc, mlmg_hibc, bc); // Same for all comps, by assumption
     visc_op.setDomainBC(mlmg_lobc, mlmg_hibc);
   }

   const int nGrowCrse = 1;

   MultiFab Wbar_crse;
   if (level > 0)
   {
      auto& crselev = getLevel(level-1);

      // Define Wbar coarse
      Wbar_crse.define(crselev.boxArray(), crselev.DistributionMap(), 1, nGrowCrse);

      // Get fillpatched coarse rho and rhoY
      FillPatchIterator fpic(crselev,Wbar_crse,nGrowCrse,time,State_Type,Density,NUM_SPECIES+1);
      MultiFab& mfc = fpic.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mfc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx = mfi.growntilebox(); 
         auto const& rho     = mfc.array(mfi,0);
         auto const& rhoY    = mfc.array(mfi,1);
         auto const& Wbar_ar = Wbar_crse.array(mfi); 
         amrex::ParallelFor(gbx, [rho, rhoY, Wbar_ar]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getMwmixGivenRY(i, j, k, rho, rhoY, Wbar_ar);
         });
      }
      visc_op.setCoarseFineBC(&Wbar_crse, crse_ratio[0]);
   }
   visc_op.setLevelBC(0, &Wbar);

   visc_op.setScalars(0.0,1.0);
   {
      FluxBoxes fb_Coeff(this,1,0);
      MultiFab** bcoeffs = fb_Coeff.get();
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
         bcoeffs[d]->setVal(1.0);
      }
      std::array<MultiFab*,AMREX_SPACEDIM> fp{AMREX_D_DECL(bcoeffs[0],bcoeffs[1],bcoeffs[2])};
      visc_op.setBCoeffs(0, amrex::GetArrOfConstPtrs(fp));
   }

   MLMG mg(visc_op);   

   FluxBoxes fb_flux(this,1,0);
   MultiFab** gradWbar = fb_flux.get();

   Diffusion::computeExtensiveFluxes(mg,Wbar,gradWbar,0, 1, &geom, -1.0);

   Vector<BCRec> math_bc(1);
   math_bc = fetchBCArray(State_Type,first_spec,1);

   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif   
   {
   FArrayBox rhoY_ed;
   for (MFIter mfi(Wbar,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      for (int d=0; d<AMREX_SPACEDIM; ++d) 
      {
         // Get edge centered rhoYs
         const Box ebx = mfi.nodaltilebox(d);
         rhoY_ed.resize(ebx,NUM_SPECIES); 
         Elixir rhoY_el = rhoY_ed.elixir();

         const Box& edomain = amrex::surroundingNodes(domain,d);
         const auto& rhoY_arr   = mf.array(mfi,1);     // mf from FPI on Density+rhoYs
         const auto& rhoYed_arr = rhoY_ed.array(0);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(d));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(d));
         amrex::ParallelFor(ebx, [d, bc_lo, bc_hi, use_harmonic_avg, rhoY_arr, rhoYed_arr, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[d] <= edomain.smallEnd(d) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[d] >= edomain.bigEnd(d) ) );
            cen2edg_cpp( i, j, k, d, NUM_SPECIES, use_harmonic_avg, on_lo, on_hi, rhoY_arr, rhoYed_arr);
         });

         auto const& rhoY        = rhoY_ed.array(0);
         auto const& gradWbar_ar = gradWbar[d]->array(mfi);
         auto const& beta_ar     = beta[d]->array(mfi);
         auto const& wbarFlux    = SpecDiffusionFluxWbar[d]->array(mfi);

         // Wbar flux is : - \rho Y_m / \overline{W} * D_m * \nabla \overline{W}
         // with beta_m = \rho * D_m below
         if ( increment_flag == 0 ) {                 // Overwrite wbar fluxes
            amrex::ParallelFor(ebx, [gradWbar_ar, beta_ar, rhoY, wbarFlux]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               // Get Wbar from rhoYs
               amrex::Real rho = 0.0;
               amrex::Real rho_inv;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  rho += rhoY(i,j,k,n);
               }
               rho_inv = 1.0 / rho;
               amrex::Real y[NUM_SPECIES];
               for (int n = 0; n < NUM_SPECIES; n++) {
                  y[n] = rhoY(i,j,k,n) * rho_inv;
               }
               amrex::Real Wbar;
               EOS::Y2WBAR(y, Wbar);
               Wbar *= 0.001;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  wbarFlux(i,j,k,n) = - y[n] / Wbar * beta_ar(i,j,k,n) * gradWbar_ar(i,j,k);
               }   
            });
         } else {                                     // Increment wbar fluxes
            amrex::ParallelFor(ebx, [gradWbar_ar, beta_ar, rhoY, wbarFlux, increment_coeff]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               // Get Wbar from rhoYs
               amrex::Real rho = 0.0;
               amrex::Real rho_inv;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  rho += rhoY(i,j,k,n);
               }
               rho_inv = 1.0 / rho;
               amrex::Real y[NUM_SPECIES];
               for (int n = 0; n < NUM_SPECIES; n++) {
                  y[n] = rhoY(i,j,k,n) * rho_inv;
               }
               amrex::Real Wbar;
               EOS::Y2WBAR(y, Wbar);
               Wbar *= 0.001;
               for (int n = 0; n < NUM_SPECIES; n++) {
                  wbarFlux(i,j,k,n) -= increment_coeff * y[n] / Wbar * beta_ar(i,j,k,n) * gradWbar_ar(i,j,k);
               }   
            });
         }
      }
   }
   }
}

void
PeleLM::differential_spec_diffuse_sync (Real dt,
                                        bool Wbar_corrector,
                                        bool last_mac_sync_iter)
{
   BL_PROFILE_REGION_START("R::PeleLM::differential_spec_diffuse_sync()");
   BL_PROFILE("PeleLM::differential_spec_diffuse_sync()");

   // Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
   // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
   // Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync

   if (hack_nospecdiff)
   {
     amrex::Error("differential_spec_diffuse_sync: hack_nospecdiff not implemented");
   }

   const Real strt_time = ParallelDescriptor::second();

   if (verbose) amrex::Print() << "Doing differential sync diffusion ..." << '\n';
   //
   // Do implicit c-n solve for each scalar...but dont reflux.
   // Save the fluxes, coeffs and source term, we need 'em for later
   // Actually, since Ssync multiplied by dt in mac_sync before
   // a call to this routine (I hope...), Ssync is in units of s,
   // not the "usual" ds/dt...lets convert it (divide by dt) so
   // we can use a generic flux adjustment function
   //
   const Real tnp1 = state[State_Type].curTime();
   FluxBoxes fb_betanp1(this, NUM_SPECIES, 0);
   MultiFab **betanp1 = fb_betanp1.get();
   getDiffusivity(betanp1, tnp1, first_spec, 0, NUM_SPECIES); // species

   MultiFab Rhs(grids,dmap,NUM_SPECIES,0,MFInfo(),Factory());
   const int spec_Ssync_sComp = first_spec - AMREX_SPACEDIM;

   //
   // Rhs and Ssync contain the RHS of DayBell:2000 Eq (18)
   // with the additional -Y_m^{n+1,p} * (delta rho)^sync term
   // Copy this into Rhs; we will need this later since we overwrite SSync
   // in the solves.
   //
   MultiFab::Copy(Rhs,Ssync,spec_Ssync_sComp,0,NUM_SPECIES,0);
   //
   // Some standard settings
   //
   const Vector<int> rho_flag(NUM_SPECIES,2);
   const MultiFab* alpha = 0;

   FluxBoxes fb_fluxSC(this);
   MultiFab** fluxSC = fb_fluxSC.get();

   const MultiFab& RhoHalftime = get_rho_half_time();

   int sdc_theta = 1.0;

   for (int sigma = 0; sigma < NUM_SPECIES; ++sigma)
   {
      //
      // Here, we use Ssync as a source in units of s, as expected by diffuse_Ssync
      // (i.e., ds/dt ~ d(Ssync)/dt, vs. ds/dt ~ Rhs in diffuse_scalar).  This was
      // apparently done to mimic diffuse_Vsync, which does the same, because the
      // diffused result is an acceleration, not a velocity, req'd by the projection.
      //
      const int ssync_ind = first_spec + sigma - Density;

      // on entry, Ssync = RHS for (delta Ytilde)^sync diffusive solve
      // on exit, Ssync = rho^{n+1} * (delta Ytilde)^sync
      // on exit, fluxSC = rhoD grad (delta Ytilde)^sync
      diffusion->diffuse_Ssync(Ssync,ssync_ind,dt,sdc_theta,
                               RhoHalftime,rho_flag[sigma],fluxSC,0,
                               betanp1,sigma,alpha,0);
      //
      // Pull fluxes into flux array
      // this is the rhoD grad (delta Ytilde)^sync terms in DayBell:2000 Eq (18)
      //
      for (int d=0; d<AMREX_SPACEDIM; ++d)
      {
        MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSC[d],0,sigma,1,0);
      }
   }
   fb_betanp1.clear();
   fb_fluxSC.clear();

   // if in the Wbar corrector, add the grad delta Wbar fluxes
   if (use_wbar && Wbar_corrector)
   {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         for (int dir=0; dir<AMREX_SPACEDIM; ++dir)
         {
            const Box& ebx = mfi.nodaltilebox(dir);
            auto const& flux_spec = SpecDiffusionFluxnp1[dir]->array(mfi);
            auto const& flux_wbar = SpecDiffusionFluxWbar[dir]->array(mfi);
            amrex::ParallelFor(ebx, NUM_SPECIES, [ flux_spec, flux_wbar ]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                flux_spec(i,j,k,n) += flux_wbar(i,j,k,n);
            });
         }
      }
   }

   //
   // Modify update/fluxes to preserve flux sum = 0
   // (Be sure to pass the "normal" looking Rhs to this generic function)
   //
   // need to correct SpecDiffusionFluxnp1 to contain rhoD grad (delta Y)^sync
   int ng = 1;
   FillPatch(*this,get_new_data(State_Type),ng,tnp1,State_Type,Density,NUM_SPECIES+2,Density);

   adjust_spec_diffusion_fluxes(SpecDiffusionFluxnp1, get_new_data(State_Type),
                                AmrLevel::desc_lst[State_Type].getBCs()[Temp],tnp1);

   //
   // Need to correct Ssync to contain rho^{n+1} * (delta Y)^sync.
   // Do this by setting
   // Ssync = "RHS from diffusion solve" + (dt*sdc_theta)*div(delta Gamma)
   //
   // Recompute update with adjusted diffusion fluxes
   //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(Ssync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx          = mfi.tilebox();
      AMREX_D_TERM(auto const& fluxX = SpecDiffusionFluxnp1[0]->array(mfi,0);,
                   auto const& fluxY = SpecDiffusionFluxnp1[1]->array(mfi,0);,
                   auto const& fluxZ = SpecDiffusionFluxnp1[2]->array(mfi,0););
      auto const& ssync      = Ssync.array(mfi,first_spec-AMREX_SPACEDIM);
      auto const& vol        = volume.const_array(mfi);
      auto const& rhs_Dsolve = Rhs.const_array(mfi);

#ifdef AMREX_USE_EB
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
      auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
      auto const& flag    = flagfab.const_array();
#endif

      Real scale = -dt * sdc_theta;

#ifdef AMREX_USE_EB
      if (flagfab.getType(bx) == FabType::covered) {              // Covered boxes
         amrex::ParallelFor(bx, NUM_SPECIES, [ssync]
         AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
         {
            ssync(i,j,k,n) = 0.0;
         });
      } else if (flagfab.getType(bx) != FabType::regular ) {     // EB containing boxes 
         auto vfrac = ebfactory.getVolFrac().const_array(mfi);
         amrex::ParallelFor(bx, [flag, vfrac, ssync, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, rhs_Dsolve, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( flag(i,j,k).isCovered() ) {
               for (int n = 0; n < NUM_SPECIES; n++) {
                  ssync(i,j,k,n) = 0.0;
               }
            } else if ( flag(i,j,k).isRegular() ) {
               fluxDivergence( i, j, k, NUM_SPECIES,
                               AMREX_D_DECL(fluxX, fluxY, fluxZ),
                               vol, scale, ssync);
               for (int n = 0; n < NUM_SPECIES; n++) {
                  ssync(i,j,k,n) += rhs_Dsolve(i,j,k,n);
               }
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               fluxDivergence( i, j, k, NUM_SPECIES,
                               AMREX_D_DECL(fluxX, fluxY, fluxZ),
                               vol, scale, ssync);
               for (int n = 0; n < NUM_SPECIES; n++) {
                  ssync(i,j,k,n) *= vfracinv;
                  ssync(i,j,k,n) += rhs_Dsolve(i,j,k,n);
               }
            }
         });
      } else {                                                   // Regular boxes
#endif
         amrex::ParallelFor(bx, [ssync, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, rhs_Dsolve, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fluxDivergence( i, j, k, NUM_SPECIES,
                            AMREX_D_DECL(fluxX, fluxY, fluxZ),
                            vol, scale, ssync);
            for (int n = 0; n < NUM_SPECIES; n++) {
               ssync(i,j,k,n) += rhs_Dsolve(i,j,k,n);
            }
         });
#ifdef AMREX_USE_EB
      }
#endif
   }

   Rhs.clear();

#ifdef AMREX_USE_EB
   EB_set_covered_faces({AMREX_D_DECL(SpecDiffusionFluxnp1[0],SpecDiffusionFluxnp1[1],SpecDiffusionFluxnp1[2])},0.);
#endif

   //
   // Do refluxing AFTER flux adjustment
   //
   if (do_reflux && level > 0 && last_mac_sync_iter)
   {
      for (int d=0; d<AMREX_SPACEDIM; ++d)
      {
         getViscFluxReg().FineAdd(*SpecDiffusionFluxnp1[d],d,0,first_spec,NUM_SPECIES,dt);
      }
   }

   if (verbose)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);
      amrex::Print() << "PeleLM::differential_spec_diffuse_sync(): lev: " << level
                     << ", time: " << run_time << '\n';
   }

   BL_PROFILE_REGION_STOP("R::PeleLM::differential_spec_diffuse_sync()");
}

void
PeleLM::reflux ()
{
   BL_PROFILE("PeleLM::reflux()");

   // no need to reflux if this is the finest level
   if (level == parent->finestLevel()) return;

   const Real strt_time = ParallelDescriptor::second();

   AMREX_ASSERT(do_reflux);

   //
   // First do refluxing step.
   //
   FluxRegister& fr_adv  = getAdvFluxReg(level+1);
   FluxRegister& fr_visc = getViscFluxReg(level+1);
   const Real    dt_crse = parent->dtLevel(level);
   const Real    scale   = 1.0/dt_crse;
   //
   // It is important, for do_mom_diff == 0, to do the viscous
   //   refluxing first, since this will be divided by rho_half
   //   before the advective refluxing is added.  In the case of
   //   do_mom_diff == 1, both components of the refluxing will
   //   be divided by rho^(n+1) in NavierStokesBase::level_sync.
   //
   // take divergence of diffusive flux registers into cell-centered RHS
   fr_visc.Reflux(Vsync,volume,scale,0,0,AMREX_SPACEDIM,geom);
   if (do_reflux_visc)
     fr_visc.Reflux(Ssync,volume,scale,AMREX_SPACEDIM,0,NUM_STATE-AMREX_SPACEDIM,geom);

   showMF("DBGSync",Ssync,"sdc_Ssync_inReflux_AftVisc",level,parent->levelSteps(level));

   const MultiFab& RhoHalftime = get_rho_half_time();

   if (do_mom_diff == 0) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& rhohalf = RhoHalftime.array(mfi);

         amrex::ParallelFor(bx, AMREX_SPACEDIM, [vsync, rhohalf]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            vsync(i,j,k,n) /= rhohalf(i,j,k);
         });
      }
   }

   // for any variables that used non-conservative advective differencing,
   // divide the sync by rhohalf
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(Ssync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhohalf = RhoHalftime.array(mfi);
      auto const& ssync   = Ssync.array(mfi);
      for (int istate = AMREX_SPACEDIM; istate < NUM_STATE; istate++)
      {
         if (advectionType[istate] == NonConservative)
         {
            const int sigma = istate - AMREX_SPACEDIM;
            amrex::ParallelFor(bx, [ssync, rhohalf, sigma]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               ssync(i,j,k,sigma) /= rhohalf(i,j,k);
            });
         }
      }
   }

   // take divergence of advective flux registers into cell-centered RHS
   fr_adv.Reflux(Vsync,volume,scale,0,0,AMREX_SPACEDIM,geom);
   fr_adv.Reflux(Ssync,volume,scale,AMREX_SPACEDIM,0,NUM_STATE-AMREX_SPACEDIM,geom);
   showMF("DBGSync",Ssync,"sdc_Ssync_inReflux_AftAdvReflux",level,parent->levelSteps(level));

   //
   // This is necessary in order to zero out the contribution to any
   // coarse grid cells which underlie fine grid cells.
   //
   BoxArray baf = getLevel(level+1).boxArray();
   baf.coarsen(fine_ratio);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      std::vector< std::pair<int,Box> > isects;
      for (MFIter mfi(Vsync,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.growntilebox();
         auto const& vsync   = Vsync.array(mfi);
         auto const& ssync   = Ssync.array(mfi);
         int nstate          = NUM_STATE;

         baf.intersections(bx,isects);

         for (int it = 0, N = isects.size(); it < N; it++) {
            amrex::ParallelFor(isects[it].second, [vsync, ssync, nstate]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               for (int n = 0; n < AMREX_SPACEDIM; n++) {
                  vsync(i,j,k,n) = 0.0;
               }
               for (int n = 0; n < nstate-AMREX_SPACEDIM; n++) {
                  ssync(i,j,k,n) = 0.0;
               }
            });
         }
      }
   }
   showMF("DBGSync",Ssync,"sdc_Ssync_inReflux_AftZerosFine",level,parent->levelSteps(level));

   if (verbose > 1)
   {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "PeleLM::Reflux(): lev: " << level << ", time: " << run_time << '\n';
   }
}

void
PeleLM::calcViscosity (const Real time,
                       const Real dt,
                       const int  iteration,
                       const int  ncycle)
{
   BL_PROFILE("PLM::calcViscosity()");

   const TimeLevel whichTime = which_time(State_Type, time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab& visc = (whichTime == AmrOldTime ? *viscn_cc : *viscnp1_cc);
   MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
  
   // Index management
   int sComp = amrex::min((int)first_spec,(int)Temp);
   int eComp = amrex::max((int)last_spec, (int)Temp);
   int nComp = eComp - sComp + 1;
   int nGrow = 1;
   int Tcomp =  Temp       - sComp;
   int RYcomp = first_spec - sComp;

   // Fillpatch the state   
   FillPatchIterator fpi(*this,S,nGrow,time,State_Type,sComp,nComp);
   MultiFab& S_cc = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& gbx = mfi.growntilebox();
      auto const& rhoY    = S_cc.array(mfi,RYcomp);
      auto const& T       = S_cc.array(mfi,Tcomp);
      auto const& mu      = visc.array(mfi);

      amrex::ParallelFor(gbx, [rhoY, T, mu]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         getVelViscosity( i, j, k, rhoY, T, mu);
      });
   }


}

void
PeleLM::calcDiffusivity (const Real time)
{
   BL_PROFILE("PLM::calcDiffusivity()");

   const TimeLevel whichTime = which_time(State_Type, time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab& diff = (whichTime == AmrOldTime) ? *diffn_cc : *diffnp1_cc;
   MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);

   // for open chambers, ambient pressure is constant in time
   Real p_amb = p_amb_old;

   if (closed_chamber == 1)
   {
     // for closed chambers, use piecewise-linear interpolation
     // we need level 0 prev and tnp1 for closed chamber algorithm
     AmrLevel& amr_lev = parent->getLevel(0);
     StateData& state_data = amr_lev.get_state_data(0);
     const Real lev_0_prevtime = state_data.prevTime();
     const Real lev_0_curtime = state_data.curTime();

     // linearly interpolate from level 0 ambient pressure
     p_amb = (lev_0_curtime - time )/(lev_0_curtime-lev_0_prevtime) * p_amb_old +
             (time - lev_0_prevtime)/(lev_0_curtime-lev_0_prevtime) * p_amb_new;
   }

   // Index management
   int sComp = amrex::min((int)first_spec, (int)Temp);
   int eComp = amrex::max((int)last_spec,  (int)Temp);
   int nComp = eComp - sComp + 1;
   int nGrow = 1;
   int Tcomp  = Temp       - sComp;
   int RYcomp = first_spec - sComp;

   // Fillpatch the state   
   FillPatchIterator fpi(*this,S,nGrow,time,State_Type,sComp,nComp);
   MultiFab& S_cc = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& gbx = mfi.growntilebox();
      auto const& rhoY    = S_cc.array(mfi,RYcomp);
      auto const& T       = S_cc.array(mfi,Tcomp);
      auto const& rhoD    = diff.array(mfi,0);
      auto const& lambda  = diff.array(mfi,NUM_SPECIES);
      auto const& mu      = diff.array(mfi,NUM_SPECIES+1);

      if ( unity_Le ) {
         amrex::Real ScInv = 1.0/schmidt;
         amrex::Real PrInv = 1.0/prandtl;
         amrex::ParallelFor(gbx, [rhoY, T, rhoD, lambda, mu, ScInv, PrInv] 
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTransportCoeffUnityLe( i, j, k, ScInv, PrInv, rhoY, T, rhoD, lambda, mu);
         });
      } else {
         amrex::ParallelFor(gbx, [rhoY, T, rhoD, lambda, mu]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTransportCoeff( i, j, k, rhoY, T, rhoD, lambda, mu);
         });
      }
   }
}

void
PeleLM::getViscosity (MultiFab* viscosity[AMREX_SPACEDIM],
                      const Real time)
{
   //
   // Select time level to work with (N or N+1)
   //
   const TimeLevel whichTime = which_time(State_Type,time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab* visc      = (whichTime == AmrOldTime) ? viscn_cc : viscnp1_cc;

#ifdef AMREX_USE_EB
   // EB : use EB CCentroid -> FCentroid
   auto math_bc = fetchBCArray(State_Type,0,1);
   EB_interp_CellCentroid_to_FaceCentroid(*visc, AMREX_D_DECL(*viscosity[0],*viscosity[1],*viscosity[2]), 
                                          0, 0, 1, geom, math_bc);
   EB_set_covered_faces({AMREX_D_DECL(viscosity[0],viscosity[1],viscosity[2])},0.0);

#else
   // NON-EB : simply use center_to_edge_fancy
   auto math_bc = fetchBCArray(State_Type,Temp,1);
   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(*visc,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         const Box ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         const auto& visc_c  = visc->array(mfi,0);
         const auto& visc_ed = viscosity[dir]->array(mfi,0);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, use_harmonic_avg, visc_c, visc_ed, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, 1, use_harmonic_avg, on_lo, on_hi, visc_c, visc_ed);
         });
      }
   }
#endif

// WARNING: maybe something specific to EB has to be done here

   if (do_LES){
      FluxBoxes mu_LES(this,1,0);
      MultiFab** mu_LES_mf = mu_LES.get();
      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
         mu_LES_mf[dir]->setVal(0., 0, mu_LES_mf[dir]->nComp(), mu_LES_mf[dir]->nGrow());
      }

      NavierStokesBase::calc_mut_LES(mu_LES_mf,time);

      for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
         MultiFab::Add(*viscosity[dir], *mu_LES_mf[dir], 0, 0, 1, 0);
     }
   }

}

void
PeleLM::getDiffusivity (MultiFab* diffusivity[AMREX_SPACEDIM],
                        const Real time,
                        const int state_comp,
                        const int dst_comp,
                        const int ncomp)
{
   BL_PROFILE("PLM::getDiffusivity()");
   AMREX_ASSERT(state_comp > Density);
   AMREX_ASSERT(diffusivity[0]->nComp() >= dst_comp+ncomp);
   AMREX_ASSERT( ( state_comp == first_spec && ncomp == NUM_SPECIES ) ||
              ( state_comp == Temp && ncomp == 1 ) );
   
   //
   // Select time level to work with (N or N+1)
   //
   const TimeLevel whichTime = which_time(State_Type,time);
   AMREX_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

   MultiFab* diff      = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;

   const int offset    = AMREX_SPACEDIM + 1;          // No diffusion coeff for vels or rho
   int       diff_comp = state_comp - offset;
   if (state_comp == Temp) diff_comp -= 1;            // Because RhoH is squeezed in between.

#ifdef AMREX_USE_EB
   // EB : use EB CCentroid -> FCentroid
   auto math_bc = fetchBCArray(State_Type,state_comp,ncomp);
   EB_interp_CellCentroid_to_FaceCentroid(*diff, AMREX_D_DECL(*diffusivity[0],*diffusivity[1],*diffusivity[2]), 
                                          diff_comp, dst_comp, ncomp, geom, math_bc);
   EB_set_covered_faces({AMREX_D_DECL(diffusivity[0],diffusivity[1],diffusivity[2])},0.0);

#else
   // NON-EB : simply use center_to_edge_fancy
   auto math_bc = fetchBCArray(State_Type,Temp,1);
   const Box& domain = geom.Domain();
   bool use_harmonic_avg = def_harm_avg_cen2edge ? true : false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(*diff,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
      {
         const Box ebx = mfi.nodaltilebox(dir);
         const Box& edomain = amrex::surroundingNodes(domain,dir);
         const auto& diff_c  = diff->array(mfi,diff_comp);
         const auto& diff_ed = diffusivity[dir]->array(mfi,dst_comp);
         const auto bc_lo = fpi_phys_loc(math_bc[0].lo(dir));
         const auto bc_hi = fpi_phys_loc(math_bc[0].hi(dir));
         amrex::ParallelFor(ebx, [dir, bc_lo, bc_hi, ncomp, use_harmonic_avg, diff_c, diff_ed, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == HT_Edge ) && ( idx[dir] <= edomain.smallEnd(dir) ) );
            bool on_hi = ( ( bc_hi == HT_Edge ) && ( idx[dir] >= edomain.bigEnd(dir) ) );
            cen2edg_cpp( i, j, k, dir, ncomp, use_harmonic_avg, on_lo, on_hi, diff_c, diff_ed);
         });
      }
   }
#endif

   if (zeroBndryVisc > 0) {
      zeroBoundaryVisc(diffusivity,time,state_comp,dst_comp,ncomp);
   }
}

void
PeleLM::zeroBoundaryVisc (MultiFab*  beta[AMREX_SPACEDIM],
                          const Real time,
                          const int  state_comp,
                          const int  dst_comp,
                          const int  ncomp) const
{
  AMREX_ASSERT(state_comp > Density);

  const auto geomdata = geom.data();

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      const Box edom = amrex::surroundingNodes(geom.Domain(),dir);
      for (MFIter mfi(*(beta[dir]),TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        const Box& ebox     = amrex::surroundingNodes(mfi.growntilebox(),dir);
        auto const& beta_arr = beta[dir]->array(mfi,dst_comp);
        amrex::ParallelFor(ebox, [beta_arr,dir,geomdata,edom,state_comp,ncomp]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          zero_visc(i, j, k, beta_arr, geomdata, edom, dir, state_comp, ncomp);
        });
      }
    }
  }
}

void
PeleLM::calc_divu (Real      time,
                   Real      dt,
                   MultiFab& divu)
{
   BL_PROFILE("PLM::calc_divu()");

   const int nGrow = 0;
   int vtCompT = NUM_SPECIES + 1;
   int vtCompY = 0;
   MultiFab mcViscTerms(grids,dmap,NUM_SPECIES+2,nGrow,MFInfo(),Factory());

   // we don't want to update flux registers due to fluxes in divu computation
   bool do_reflux_hold = do_reflux;
   do_reflux = false;

   // DD is computed and stored in divu, and passed in to initialize the calc of divu
   bool include_Wbar_terms = true;
   compute_differential_diffusion_terms(mcViscTerms,divu,time,dt,include_Wbar_terms);
 
   do_reflux = do_reflux_hold;

   // if we are in the initial projection (time=0, dt=-1), set RhoYdot=0
   // if we are in a divu_iter (time=0, dt>0), use I_R
   // if we are in an init_iter or regular time step (time=dt, dt>0), use instananeous
   MultiFab   RhoYdotTmp;
   MultiFab&  S       = get_data(State_Type,time);
   const bool use_IR  = (time == 0 && dt > 0);

   MultiFab& RhoYdot = (use_IR) ? get_new_data(RhoYdot_Type) : RhoYdotTmp;

   if (!use_IR)
   {
     if (time == 0)
     {
       // initial projection, set omegadot to zero
       RhoYdot.define(grids,dmap,NUM_SPECIES,0,MFInfo(),Factory());
       RhoYdot.setVal(0.0);
     }
     else if (dt > 0)
     {
       // init_iter or regular time step, use instantaneous omegadot
       RhoYdot.define(grids,dmap,NUM_SPECIES,0,MFInfo(),Factory());
       compute_instantaneous_reaction_rates(RhoYdot,S,time,nGrow);
     }
     else
     {
       amrex::Abort("bad divu_logic - shouldn't be here");
     }
   }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhoY    = S.array(mfi,first_spec);
      auto const& T       = S.array(mfi,Temp);
      auto const& vtT     = mcViscTerms.array(mfi,vtCompT);
      auto const& vtY     = mcViscTerms.array(mfi,vtCompY);
      auto const& rhoYdot = RhoYdot.array(mfi);
      auto const& du      = divu.array(mfi);

      amrex::ParallelFor(bx, [ rhoY, T, vtT, vtY, rhoYdot, du]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         compute_divu( i, j, k, rhoY, T, vtT, vtY, rhoYdot, du );
      });
   }
  
#ifdef AMREX_USE_EB    
   EB_set_covered(divu,0.);
#endif
}

//
// Compute the Eulerian Dp/Dt for use in pressure relaxation.
//
void
PeleLM::calc_dpdt (Real      time,
                   Real      dt,
                   MultiFab& dpdt,
                   MultiFab* umac)
{
   BL_PROFILE("PeleLM::calc_dpdt()");

   // for open chambers, ambient pressure is constant in time
   Real p_amb = p_amb_old;

   if (closed_chamber == 1)
   {
     // for closed chambers, use piecewise-linear interpolation
     // we need level 0 prev and tnp1 for closed chamber algorithm
     AmrLevel& amr_lev = parent->getLevel(0);
     StateData& state_data = amr_lev.get_state_data(0);
     const Real lev_0_prevtime = state_data.prevTime();
     const Real lev_0_curtime = state_data.curTime();

     // linearly interpolate from level 0 ambient pressure
     p_amb = (lev_0_curtime - time )/(lev_0_curtime-lev_0_prevtime) * p_amb_old +
             (time - lev_0_prevtime)/(lev_0_curtime-lev_0_prevtime) * p_amb_new;
   }
  
   if (dt <= 0.0 || dpdt_factor <= 0)
   {
     dpdt.setVal(0);
     return;
   }

   int nGrow = dpdt.nGrow();
   FillPatchIterator S_fpi(*this,get_new_data(State_Type),nGrow,time,State_Type,RhoRT,1);
   MultiFab& Peos=S_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(dpdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& dPdt  = dpdt.array(mfi);
      auto const& P     = Peos.array(mfi);
      Real dpdt_fac     = dpdt_factor;

      amrex::ParallelFor(bx, [dPdt, P, p_amb, dt, dpdt_fac]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         dPdt(i,j,k) = (P(i,j,k) - p_amb) / ( dt * P(i,j,k) ) * dpdt_fac;
      });
   }

   Peos.clear();

   if (nGrow > 0) {
     dpdt.FillBoundary(0,1, geom.periodicity());
     AMREX_ASSERT(dpdt.nGrow() == 1);
     Extrapolater::FirstOrderExtrap(dpdt, geom, 0, 1);
   }
}

void
PeleLM::RhoH_to_Temp (MultiFab& S,
                      int       nGrow,
                      int       dominmax)
{
  BL_PROFILE("PLM::RhoH_to_Temp()");
  const Real strt_time = ParallelDescriptor::second();

  // TODO: simplified version of that function for now: no iters, no tols, ... PPhys need to be fixed

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
     const Box& bx    = mfi.tilebox();
     auto const& T    = S.array(mfi,Temp);
     auto const& rho  = S.array(mfi,Density);
     auto const& rhoY = S.array(mfi,first_spec);
     auto const& rhoH = S.array(mfi,RhoH);
     amrex::ParallelFor(bx, [T,rho,rhoY,rhoH]
     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
     {
        getTfromHY(i,j,k, rho, rhoY, rhoH, T);
     });
  }

  /*
  //
  // If this hasn't been set yet, we cannnot do it correct here (scan multilevel),
  // just wing it.
  //
  const Real htt_hmixTYP_SAVE = htt_hmixTYP; 
  if (htt_hmixTYP <= 0)
  {
    if (typical_values[RhoH]==typical_RhoH_value_default)
    {
      htt_hmixTYP = S.norm0(RhoH);
    }
    else
    {
      htt_hmixTYP = typical_values[RhoH];
    }        
    if (verbose) amrex::Print() << "setting htt_hmixTYP = " << htt_hmixTYP << '\n';
  }

  int max_iters = 0;
  AMREX_ALWAYS_ASSERT(nGrow <= S.nGrow());

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(max:max_iters)
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& box = mfi.growntilebox(nGrow);
    max_iters = std::max(max_iters, RhoH_to_Temp(S[mfi],box,dominmax));
  }
  */

   if (dominmax) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(S,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& T = S.array(mfi,Temp);
         amrex::Real tempmin = htt_tempmin;
         amrex::Real tempmax = htt_tempmax;

         amrex::ParallelFor(bx, [tempmin, tempmax, T]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fabMinMax( i, j, k, 1, tempmin, tempmax, T);
         });
      }
   }

  if (verbose > 1)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

   // ParallelDescriptor::ReduceIntMax(max_iters,IOProc);
    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (verbose) amrex::Print() << "PeleLM::RhoH_to_Temp: time: " << run_time << '\n';
  }
  //
  // Reset it back.
  //
  //htt_hmixTYP = htt_hmixTYP_SAVE;
}

void
PeleLM::getForce(FArrayBox&       force,
                 const Box&       bx,
                 int              scomp,
                 int              ncomp,
                 const Real       time,
                 const FArrayBox& Vel,
                 const FArrayBox& Scal,
                 int              scalScomp)
{
   AMREX_ASSERT(force.box().contains(bx));

   const auto& velocity = Vel.array();
   const auto& scalars  = Scal.array(scalScomp);
   const auto& f        = force.array(scomp);

   const auto  dx       = geom.CellSizeArray();
   RealBox gridloc      = RealBox(bx,geom.CellSize(),geom.ProbLo());
   const Real  grav     = gravity;
   const int   nscal    = NUM_SCALARS;

   // Get non-static info for the pseudo gravity forcing
   int pseudo_gravity    = ctrl_pseudoGravity;
   const Real dV_control = ctrl_dV;

   amrex::ParallelFor(bx, [f, scalars, velocity, time, grav, pseudo_gravity, dV_control, dx, scomp, ncomp]
   AMREX_GPU_DEVICE(int i, int j, int k) noexcept
   {
      makeForce(i,j,k, scomp, ncomp, pseudo_gravity, time, grav, dV_control, dx, velocity, scalars, f);
   });
}

void
PeleLM::setPlotVariables ()
{
  AmrLevel::setPlotVariables();

  Vector<std::string> names;
  EOS::speciesNames(names);

  
// Here we specify to not plot all the rho.Y from the state variables
// because it takes a lot of space in memory and disk usage
// To plot rho.Y, we pass into a derive function (see PeleLM_setup.cpp)
// and it can be activated with "amr.derive_plot_vars=rhoY" in the input file
  for (int i = 0; i < names.size(); i++)
  {
    const std::string name = "rho.Y("+names[i]+")";
    parent->deleteStatePlotVar(name);
  }

  if (verbose)
  {
    amrex::Print() << "\nState Plot Vars: ";

    std::list<std::string>::const_iterator li = 
      parent->statePlotVars().begin(), end = parent->statePlotVars().end();

    for ( ; li != end; ++li)
      amrex::Print() << *li << ' ';
    amrex::Print() << '\n';

    amrex::Print() << "\nDerive Plot Vars: ";

    li  = parent->derivePlotVars().begin();
    end = parent->derivePlotVars().end();

    for ( ; li != end; ++li)
      amrex::Print() << *li << ' ';
    amrex::Print() << '\n';
  }
}

void
PeleLM::writePlotFile (const std::string& dir,
                       std::ostream&  os,
                       VisMF::How     how)
{
  if ( ! Amr::Plot_Files_Output() ) return;

  ParmParse pp("ht");

  pp.query("plot_rhoydot",plot_rhoydot);
  //
  // Note that this is really the same as its NavierStokes counterpart,
  // but in order to add diagnostic MultiFabs into the plotfile, code had
  // to be interspersed within this function.
  //

  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int,int> > plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++)
  {
    for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
    {
      const std::string& name = desc_lst[typ].name(comp);

      if (parent->isStatePlotVar(name) && desc_lst[typ].getType() == IndexType::TheCellType())
      {
        //
        // When running SDC we get things of this form in the State.
        // I want a simple way not to write'm out to plotfiles.
        //
        if (name.find("I_R[") != std::string::npos)
        {
          if (plot_rhoydot)
            plot_var_map.push_back(std::pair<int,int>(typ,comp));
        }
        else
        {
          plot_var_map.push_back(std::pair<int,int>(typ,comp));
        }
      }
    }
  }

  int num_derive = 0;
  std::list<std::string> derive_names;
  const std::list<DeriveRec>& dlist = derive_lst.dlist();

  for (std::list<DeriveRec>::const_iterator it = dlist.begin(), end = dlist.end();
       it != end;
       ++it)
  {
    if (parent->isDerivePlotVar(it->name()))
    {
      derive_names.push_back(it->name());
      num_derive += it->numDerive();
    }
  }

  int num_auxDiag = 0;
  for (const auto& kv : auxDiag)
  {
    num_auxDiag += kv.second->nComp();
  }

  int n_data_items = plot_var_map.size() + num_derive + num_auxDiag;

#ifdef AMREX_USE_EB
    // add in vol frac
    n_data_items++;
#endif

  Real tnp1 = state[State_Type].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor())
  {
    //
    // The first thing we write out is the plotfile type.
    //
    os << thePlotFileType() << '\n';

    if (n_data_items == 0)
      amrex::Error("Must specify at least one valid data item to plot");

    os << n_data_items << '\n';

    //
    // Names of variables -- first state, then derived
    //
    for (int i =0; i < plot_var_map.size(); i++)
    {
      int typ  = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      os << desc_lst[typ].name(comp) << '\n';
    }

    for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
         it != end;
         ++it)
    {
      const DeriveRec* rec = derive_lst.get(*it);
      for (int i = 0; i < rec->numDerive(); i++)
        os << rec->variableName(i) << '\n';
    }

    //
    // Hack in additional diagnostics.
    //
    for (std::map<std::string,Vector<std::string> >::const_iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
         it != end;
         ++it)
    {
      for (int i=0; i<it->second.size(); ++i)
        os << it->second[i] << '\n';
    }

#ifdef AMREX_USE_EB
	//add in vol frac
	os << "volFrac\n";
#endif

    os << AMREX_SPACEDIM << '\n';
    os << parent->cumTime() << '\n';
    int f_lev = parent->finestLevel();
    os << f_lev << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; i++)
      os << Geom().ProbLo(i) << ' ';
    os << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; i++)
      os << Geom().ProbHi(i) << ' ';
    os << '\n';
    for (int i = 0; i < f_lev; i++)
      os << parent->refRatio(i)[0] << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
      os << parent->Geom(i).Domain() << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
      os << parent->levelSteps(i) << ' ';
    os << '\n';
    for (int i = 0; i <= f_lev; i++)
    {
      for (int k = 0; k < AMREX_SPACEDIM; k++)
        os << parent->Geom(i).CellSize()[k] << ' ';
      os << '\n';
    }
    os << (int) Geom().Coord() << '\n';
    os << "0\n"; // Write bndry data.

    // job_info file with details about the run
    std::ofstream jobInfoFile;
    std::string FullPathJobInfoFile = dir;
    std::string PrettyLine = "===============================================================================\n";

    FullPathJobInfoFile += "/job_info";
    jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

    // job information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Job Information\n";
    jobInfoFile << PrettyLine;
	
    jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
    jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
    jobInfoFile << "\n\n";

    // plotfile information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Plotfile Information\n";
    jobInfoFile << PrettyLine;

    time_t now = time(0);

    // Convert now to tm struct for local timezone
    tm* localtm = localtime(&now);
    jobInfoFile   << "output data / time: " << asctime(localtm);

    char currentDir[FILENAME_MAX];
    if (getcwd(currentDir, FILENAME_MAX)) {
      jobInfoFile << "output dir:         " << currentDir << "\n";
    }

    jobInfoFile << "\n\n";


    // build information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Build Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
    jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
    jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
    jobInfoFile << "BoxLib dir:    " << buildInfoGetAMReXDir() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
    jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
    jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
    jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

    jobInfoFile << "\n";

    for (int nn = 1; nn <= buildInfoGetNumModules(); nn++) {
      jobInfoFile << buildInfoGetModuleName(nn) << ": " << buildInfoGetModuleVal(nn) << "\n";
    }

    jobInfoFile << "\n";

    const char* githash1 = buildInfoGetGitHash(1);
    const char* githash2 = buildInfoGetGitHash(2);
    const char* githash3 = buildInfoGetGitHash(3);
    const char* githash4 = buildInfoGetGitHash(4);
    if (strlen(githash1) > 0) {
      jobInfoFile << "PeleLM git hash: " << githash1 << "\n";
    }
    if (strlen(githash2) > 0) {
      jobInfoFile << "AMReX git hash: " << githash2 << "\n";
    }
    if (strlen(githash3) > 0) {
      jobInfoFile << "IAMR   git hash: " << githash3 << "\n";
    }
    if (strlen(githash4) > 0) {
      jobInfoFile << "PelePhysics git hash: " << githash4 << "\n";
    }

    jobInfoFile << "\n\n";


    // runtime parameters
    jobInfoFile << PrettyLine;
    jobInfoFile << " Inputs File Parameters\n";
    jobInfoFile << PrettyLine;
	
    ParmParse::dumpTable(jobInfoFile, true);

    jobInfoFile.close();	

  }
  // Build the directory to hold the MultiFab at this level.
  // The name is relative to the directory containing the Header file.
  //
  static const std::string BaseName = "/Cell";

  std::string LevelStr = amrex::Concatenate("Level_", level, 1);
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    FullPath += '/';
  FullPath += LevelStr;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (ParallelDescriptor::IOProcessor())
    if (!amrex::UtilCreateDirectory(FullPath, 0755))
      amrex::CreateDirectoryFailed(FullPath);
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();

  if (ParallelDescriptor::IOProcessor())
  {
    os << level << ' ' << grids.size() << ' ' << tnp1 << '\n';
    os << parent->levelSteps(level) << '\n';

    for (int i = 0; i < grids.size(); ++i)
    {
      RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
      for (int n = 0; n < AMREX_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0)
    {
      std::string PathNameInHeader = LevelStr;
      PathNameInHeader += BaseName;
      os << PathNameInHeader << '\n';
    }
    
#ifdef AMREX_USE_EB
	// volfrac threshhold for amrvis
	// fixme? pulled directly from CNS, might need adjustment
        if (level == parent->finestLevel()) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                os << "1.0e-6\n";
            }
        }
#endif
    
    
  }
  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  // NOTE: we are assuming that each state variable has one component,
  // but a derived variable is allowed to have multiple components.
  int       cnt   = 0;
  int       ncomp = 1;
  const int nGrow = 0;
  //MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
  MultiFab plotMF(grids,dmap,n_data_items,nGrow,MFInfo(),Factory());
  MultiFab* this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (int i = 0; i < plot_var_map.size(); i++)
  {
    int typ  = plot_var_map[i].first;
    int comp = plot_var_map[i].second;
    this_dat = &state[typ].newData();
    MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
    cnt+= ncomp;
  }
  //
  // Cull data from derived variables.
  // 
  Real plot_time;

  if (derive_names.size() > 0)
  {
    for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
         it != end;
         ++it) 
    {
      if (*it == "avg_pressure" ||
          *it == "gradpx"       ||
          *it == "gradpy"       ||
          *it == "gradpz") 
      {
        if (state[Press_Type].descriptor()->timeType() == 
            StateDescriptor::Interval) 
        {
          plot_time = tnp1;
        }
        else
        {
          int f_lev = parent->finestLevel();
          plot_time = getLevel(f_lev).state[Press_Type].curTime();
        }
      }
      else
      {
        plot_time = tnp1;
      } 
      const DeriveRec* rec = derive_lst.get(*it);
      ncomp = rec->numDerive();
      auto derive_dat = derive(*it,plot_time,nGrow);
      MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
      cnt += ncomp;
    }
  }
  
 
  //
  // Cull data from diagnostic multifabs.
  //
  for (const auto& kv : auxDiag)
  {
      int nComp = kv.second->nComp();
      MultiFab::Copy(plotMF,*kv.second,0,cnt,nComp,nGrow);
      cnt += nComp;
  }

#ifdef AMREX_USE_EB
    // add volume fraction to plotfile
    plotMF.setVal(0.0, cnt, 1, nGrow);
    MultiFab::Copy(plotMF,*volfrac,0,cnt,1,nGrow);

    // set covered values for ease of viewing
    EB_set_covered(plotMF, 0.0);
#endif

  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;

//#ifdef AMREX_USE_EB
//  bool plot_centroid_data = true;
//  pp.query("plot_centroid_data",plot_centroid_data);
//  if (plot_centroid_data) {
//    MultiFab TT(grids,dmap,plotMF.nComp(),plotMF.nGrow()+1,MFInfo(),Factory());
//    MultiFab::Copy(TT,plotMF,0,0,plotMF.nComp(),plotMF.nGrow());
//    TT.FillBoundary(geom.periodicity());
//    EB_interp_CC_to_Centroid(plotMF,TT,0,0,plotMF.nComp(),geom);
//  }
//#endif

  VisMF::Write(plotMF,TheFullPath,how);
}

std::unique_ptr<MultiFab>
PeleLM::derive (const std::string& name,
                Real               time,
                int                ngrow)
{        
  AMREX_ASSERT(ngrow >= 0);
  
  std::unique_ptr<MultiFab> mf;
  const DeriveRec* rec = derive_lst.get(name);
  if (rec)
  {
    mf.reset(new MultiFab(grids, dmap, rec->numDerive(), ngrow));
    int dcomp = 0;
    derive(name,time,*mf,dcomp);
  }
  else
  {
    mf = std::move(AmrLevel::derive(name,time,ngrow));
  }

  if (mf==nullptr) {
    std::string msg("PeleLM::derive(): unknown variable: ");
    msg += name;
    amrex::Error(msg.c_str());
  }
  return mf;
}
 
void
PeleLM::derive (const std::string& name,
                Real               time,
                MultiFab&          mf,
                int                dcomp)
{
#ifdef AMREX_PARTICLES
   ParticleDerive(name,time,mf,dcomp);
#else
   AmrLevel::derive(name,time,mf,dcomp);
#endif
}

void
PeleLM::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{
  BL_PROFILE("PLM::errorEst()");

  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

#ifdef AMREX_USE_EB
  if (refine_cutcells) {
        const MultiFab& S_new = get_new_data(State_Type);
        amrex::TagCutCells(tags, S_new);
  }
#endif

  for (int j=0; j<errtags.size(); ++j) {
    std::unique_ptr<MultiFab> mf;
    if (errtags[j].Field() != std::string()) {
      mf = derive(errtags[j].Field(), time, errtags[j].NGrow());
    }
    errtags[j](tags,mf.get(),clearval,tagval,time,level,geom);
  }
}

void
PeleLM::activeControl(const int  step,
                      const int restart,
                      const Real time,
                      const Real dt)
{

   if (!ctrl_active) return;

   const int finest_level = parent->finestLevel();

   // Get the current target state (either amount of fuel or T-iso position)
   Real coft = 0.0;

   if ( !ctrl_use_temp ) {
      // Get fuel mass
      Real fuelmass = 0.0;
      std::string fuel = "rho.Y(" + fuelName + ")";
      for (int lev = 0; lev <= finest_level; lev++) {
        fuelmass += getLevel(lev).volWgtSum(fuel,time);
      }
      coft = fuelmass;
   } else {
      // Get the low T position
      coft = 1.e37;
      for (int lev = 0; lev <= finest_level; lev++) {

         // Level data
         MultiFab&   mf     = getLevel(lev).get_new_data(State_Type);
         const auto geomdata = getLevel(lev).geom.data();

         // Static -> local AC data
         int  AC_FlameDir   = ctrl_flameDir;
         Real AC_Tcross     = ctrl_temperature;

         // FPI on temp
         FillPatchIterator Tfpi(getLevel(lev),mf,1,time,State_Type,Temp,1);
         MultiFab& Tmf=Tfpi.get_mf();

         // If not finest, set covered to 0.0 : is it really necessary ?
         if ( lev < finest_level ) {
            auto bac = getLevel(lev).boxArray();
            auto baf = getLevel(lev+1).boxArray();
            baf.coarsen(fine_ratio);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
               std::vector< std::pair<int,Box> > isects;
               for (MFIter mfi(Tmf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
               {
                  auto const& T_arr = Tmf.array(mfi);
                  baf.intersections(bac[mfi.index()],isects);
                  for (int is = 0; is < isects.size(); is++) {
                     amrex::ParallelFor(isects[is].second, [T_arr]
                     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                     {
                        T_arr(i,j,k) = 0.0;
                     });
                  }
               }
            }
         }

         Real lowT = amrex::ReduceMin(Tmf, 0, [geomdata,AC_Tcross,AC_FlameDir]
                     AMREX_GPU_HOST_DEVICE(Box const& bx, Array4<Real const> const& T_arr ) noexcept -> Real
                     {
                        using namespace amrex::literals;
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);
                        amrex::Real tmp_pos = 1.e37_rt;
                        const amrex::Real* prob_lo = geomdata.ProbLo();
                        const amrex::Real* dx = geomdata.CellSize();
                        Real coor[3] = {0.0};
                        for       (int k = lo.z; k <= hi.z; ++k) {
                            for    (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    Real lcl_pos = 1.e37_rt;
                                    if ( T_arr(i,j,k) > AC_Tcross ) {
                                        int idx[3] = {i,j,k};
                                        idx[AC_FlameDir] -= 1;
                                        if ( T_arr(idx[0],idx[1],idx[2]) < AC_Tcross ) {
                                            AMREX_D_TERM(coor[0] = prob_lo[0] + (i+0.5)*dx[0];,
                                                         coor[1] = prob_lo[1] + (j+0.5)*dx[1];,
                                                         coor[2] = prob_lo[2] + (k+0.5)*dx[2];);
                                            Real slope = ((T_arr(i,j,k) ) - T_arr(idx[0],idx[1],idx[2]))/dx[AC_FlameDir];
                                            lcl_pos = coor[AC_FlameDir] - dx[AC_FlameDir] + ( AC_Tcross - T_arr(idx[0],idx[1],idx[2]) ) / slope;
                                        }
                                    }
                                    tmp_pos = amrex::min(tmp_pos,lcl_pos);
                                }
                            }
                        }
                        return tmp_pos;
                     });
          coft = amrex::min(coft,lowT);
      }
      ParallelDescriptor::ReduceRealMin(coft);
   }

   // If first time, set old target state
   if ( ctrl_coftOld < 0.0 ) ctrl_coftOld = coft;

   if ( step == 0 ) ctrl_nfilled = ctrl_NavgPts+1;

   // Manage AC history I/O
   if ( restart ) {
      ctrl_nfilled = ctrl_NavgPts+1;
      struct stat buffer;
      bool have_history = (stat (ctrl_AChistory.c_str(), &buffer) == 0);
      if ( have_history ) {
         if ( ctrl_verbose ) Print() << " Setting AC from history from " << ctrl_AChistory << "\n";
         int step_io;
         Real time_io, vel_io, slocal_io, dV_io, s_est_io, coft_old_io;
         std::fstream ACfile (ctrl_AChistory.c_str(), std::fstream::in);
         if (ACfile.is_open()) {
            while ( ACfile.good() ) {
               ACfile >> step_io
                      >> time_io
                      >> vel_io
                      >> slocal_io
                      >> dV_io
                      >> s_est_io
                      >> coft_old_io;
               if ( ( ( step - step_io ) >= 0 ) &&
                    ( ( step - step_io ) <= ctrl_NavgPts ) ) {   // Fill the previous step data
                  int pos = ( step - step_io );
                  ctrl_nfilled -= 1;
                  ctrl_time_pts[pos] = time_io;
                  ctrl_velo_pts[pos] = vel_io;
                  ctrl_cntl_pts[pos] = coft_old_io;
               }
               if ( step_io == step ) {                          // Found the right step
                  ctrl_V_in = vel_io;
                  ctrl_tBase = time_io;
                  ctrl_dV = dV_io;
                  ctrl_sest = s_est_io;
                  ctrl_coftOld = coft_old_io;
               }
            }
            ACfile.close();
         }
         if ( ctrl_verbose ) {
            Print() << " AC history arrays: \n";
            for (int n = 0; n < ctrl_time_pts.size(); n++) {
               Print() << "  ["<<n<<"] time: "<< ctrl_time_pts[n]
                                <<", velo: "<< ctrl_velo_pts[n]
                                <<", coft: "<< ctrl_cntl_pts[n] << "\n";
            }
         }
      } else {
         if ( ctrl_verbose ) Print() << " AC history file " << ctrl_AChistory << " not found, restarting from scratch \n";
      }
   }

   //Real zbase_control += ctrl_V_in * dt + ctrl_dV * dt * dt;
   ctrl_V_in_old = ctrl_V_in;
   if ( dt > 0.0 ) {
      ctrl_V_in += dt * ctrl_dV;
   }

   Real slocal = 0.5 * (ctrl_V_in_old + ctrl_V_in) - (coft - ctrl_coftOld) / ( dt * ctrl_scale );

   // -------------------------------------------
   // Get s_est averaged from previous N steps
   // -------------------------------------------
   // Update ctrl_* Vectors if not restarting
   if (!restart) {
      ctrl_time_pts.insert(ctrl_time_pts.begin(),time);
      ctrl_velo_pts.insert(ctrl_velo_pts.begin(),ctrl_V_in);
      ctrl_cntl_pts.insert(ctrl_cntl_pts.begin(),coft);
      if ( ctrl_time_pts.size() > ctrl_NavgPts ) { // Pop_back only if it's filled
         ctrl_time_pts.pop_back();
         ctrl_velo_pts.pop_back();
         ctrl_cntl_pts.pop_back();
      }
   }

   if ( ctrl_nfilled <= 0 ) {
      Real velIntegral = 0.0;
      for (int n = 1; n <= ctrl_NavgPts; n++ ) { // Piecewise constant velocity over NavgPts last steps
         velIntegral += 0.5 * ( ctrl_velo_pts[n-1] + ctrl_velo_pts[n] )
                            * ( ctrl_time_pts[n-1] - ctrl_time_pts[n] );
      }
      ctrl_sest = ( velIntegral - ( ctrl_cntl_pts[0] - ctrl_cntl_pts[ctrl_NavgPts]) / ctrl_scale )
                 / ( ctrl_time_pts[0] - ctrl_time_pts[ctrl_NavgPts] );

   } else {
      ctrl_nfilled -= 1;
      if ( step != 0 ) {
         ctrl_sest = (1.0 - ctrl_corr ) * ctrl_sest + ctrl_corr * slocal;
      }
   }

   // Compute Vnew
   int ctrl_method = 3;
   Real Vnew = 0.0;

   if (ctrl_method == 1) {         // linear
      Real vslope = 2.0 * ( (ctrl_cfix - coft)/(ctrl_scale * ctrl_tauControl) + ctrl_sest - ctrl_V_in ) / ctrl_tauControl;
      Vnew = ctrl_V_in + dt * vslope;
   } else if (ctrl_method == 2) { // Quadratic 1
      Real vslope = 3.0 * ( (ctrl_cfix - coft)/(ctrl_scale * ctrl_tauControl) + ctrl_sest - ctrl_V_in ) / ctrl_tauControl;
      Vnew = ctrl_V_in + ( dt - 0.5 * dt*dt / ctrl_tauControl ) * vslope;
   } else if (ctrl_method == 3) { // Quadratic 2
      Real rhs2 = 2.0 * ( (ctrl_cfix - coft)/(ctrl_scale * ctrl_tauControl) + ctrl_sest - ctrl_V_in ) / ctrl_tauControl;
      Real rhs1 = ( ctrl_sest - ctrl_V_in ) / ctrl_tauControl;
      Real vt_tay = 3.0 * rhs2 - 2.0 * rhs1;
      Real vtt_tay = 6.0 * ( rhs1 - rhs2 ) / ctrl_tauControl;
      Vnew = ctrl_V_in + dt * vt_tay + 0.5 * dt * dt * vtt_tay;
   }

   // Limit Vnew
   Real dVmax = ctrl_changeMax * 1.0;
   Real dVmin = ctrl_changeMax * std::max(1.0,ctrl_V_in);
   Vnew = std::max(Vnew,0.0);
   Vnew = std::min(std::max(Vnew,ctrl_V_in-dVmin),ctrl_V_in+dVmax);
   Vnew = std::min(Vnew,ctrl_velMax);

   if (!restart && step > 0) {
      ctrl_tBase = time;

      // Compute
      ctrl_dV = ( Vnew - ctrl_V_in ) / dt;
   }

   ctrl_coftOld = coft;

   // Pass dV and ctrl_V_in to GPU compliant problem routines
   ac_parm->ctrl_V_in = ctrl_V_in;
   ac_parm->ctrl_dV = ctrl_dV;
   ac_parm->ctrl_tBase = ctrl_tBase;

   // Verbose
   if ( ctrl_verbose && !restart) {
      Print() << "\n------------------------AC CONTROL------------------------ \n";
      Print() << " |     Time: " << time << " -     dt: " << dt << "\n";
      Print() << " | Position: " << coft/ctrl_scale << " - Target: " << ctrl_cfix/ctrl_scale << "\n";
      Print() << " |    V_new: " << Vnew << " -   V_in: " << ctrl_V_in << " - dV: " << ctrl_dV << "\n";
      Print() << "---------------------------------------------------------- \n";
   }

   // Append to (or create) AC history file
   if (!restart) {
      std::ofstream ACfile(ctrl_AChistory.c_str(), std::ofstream::out | std::ofstream::app );
      Print(ACfile).SetPrecision(15) << step << "  "
                                     << time << "  "
                                     << ctrl_V_in << "  "
                                     << slocal << "  "
                                     << ctrl_dV << "  "
                                     << ctrl_sest << "  "
                                     << ctrl_coftOld << std::endl;
   }
}

void
PeleLM::initActiveControl()
{
   if (level!=0) return;

   if (verbose) Print() <<  "PeleLM::initActiveControl()\n";

   // Parse input file
   ParmParse ppAC("active_control");
   ppAC.query("on",ctrl_active);
   ppAC.query("use_temp",ctrl_use_temp);
   ppAC.query("temperature",ctrl_temperature);
   ppAC.query("tau",ctrl_tauControl);
   ppAC.query("v",ctrl_verbose);
   ppAC.query("height",ctrl_h);
   ppAC.query("changeMax",ctrl_changeMax);
   ppAC.query("velMax",ctrl_velMax);
   ppAC.query("pseudo_gravity",ctrl_pseudoGravity);
   ppAC.query("flame_direction",ctrl_flameDir);
   ppAC.query("AC_history",ctrl_AChistory);
   ppAC.query("npoints_average",ctrl_NavgPts);

   // Active control checks
   if ( ctrl_use_temp && (ctrl_temperature <= 0.0) )
      amrex::Error("active_control.temperature MUST be set with active_control.use_temp = 1");

   if ( ctrl_active && (ctrl_tauControl <= 0.0) )
      amrex::Error("active_control.tau MUST be set when using active_control");

   if ( ctrl_active && (ctrl_h <= 0.0) )
      amrex::Error("active_control.height MUST be set when using active_control");

   if ( ctrl_active && ( ctrl_flameDir > 2 ) )
      amrex::Error("active_control.flame_direction MUST be 0, 1 or 2 for X, Y and Z resp.");

   // Exit here if AC not activated
   if ( !ctrl_active ) return;

   // Activate AC in problem specific / GPU compliant sections
   ac_parm->ctrl_active = ctrl_active;

   // Resize vector for temporal average
   ctrl_time_pts.resize(ctrl_NavgPts+1,-1.0);
   ctrl_velo_pts.resize(ctrl_NavgPts+1,-1.0);
   ctrl_cntl_pts.resize(ctrl_NavgPts+1,-1.0);

   // Compute some active control parameters
   Real area_tot = 1.0;
   const Real* problo = geom.ProbLo();
   const Real* probhi = geom.ProbHi();
   for (int dim = 0; dim < AMREX_SPACEDIM; dim++) {
      if (dim != ctrl_flameDir) area_tot *= (probhi[dim] - problo[dim]);
   }

   // Extract data from bc: assumes flow comes in from lo side of ctrl_flameDir
   ProbParm const* lprobparm = prob_parm.get();
   ACParm const* lacparm = ac_parm.get(); 
   PmfData const* lpmfdata = pmf_data_g;
   amrex::Gpu::DeviceVector<amrex::Real> s_ext_v(DEF_NUM_STATE);
   amrex::Real* s_ext_d = s_ext_v.data();
   amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(problo[0],problo[1],problo[2])};
   x[ctrl_flameDir] -= 1.0;
   const int ctrl_flameDir_l = ctrl_flameDir;
   const amrex::Real time_l = -1.0;
   const auto geomdata = geom.data();
   Box dumbx({AMREX_D_DECL(0,0,0)},{AMREX_D_DECL(0,0,0)});
   amrex::ParallelFor(dumbx, [x,s_ext_d,ctrl_flameDir_l,time_l,geomdata,lprobparm,lacparm, lpmfdata]
   AMREX_GPU_DEVICE(int i, int j, int k) noexcept
   {
      bcnormal(x, s_ext_d, ctrl_flameDir_l, 1, time_l, geomdata, *lprobparm, *lacparm, lpmfdata);
   });
   amrex::Real s_ext[DEF_NUM_STATE];
#if AMREX_USE_GPU
   amrex::Gpu::dtoh_memcpy(s_ext,s_ext_d,sizeof(amrex::Real)*DEF_NUM_STATE);
#else
   std::memcpy(s_ext,s_ext_d,sizeof(amrex::Real)*DEF_NUM_STATE);
#endif

   if ( !ctrl_use_temp ) {
      // Get the fuel rhoY
      Vector<std::string> specNames;
      EOS::speciesNames(specNames);
      int fuelidx = -1;
      for (int k = 0; k < NUM_SPECIES; k++) {
         if ( !specNames[k].compare(fuelName) ) fuelidx = k;
      }
      ctrl_scale = area_tot * s_ext[first_spec+fuelidx];
      ctrl_cfix  = ctrl_h * ctrl_scale;
   } else {
      // If using temp, scale is 1.0
      ctrl_scale = 1.0;
      ctrl_cfix = ctrl_h;
   }

   // Extract initial ctrl_V_in from bc
   ctrl_V_in  = s_ext[ctrl_flameDir];
   ctrl_V_in_old = ctrl_V_in;

   // Pass V_in to bc
   ac_parm->ctrl_V_in = ctrl_V_in;

   if ( ctrl_verbose && ctrl_active ) {
      if ( ctrl_use_temp ) {
         Print() << " Active control based on temperature iso-level activated."
                 << " Maintaining the flame at " << ctrl_h << " in " << ctrl_flameDir << " direction. \n";
      } else {
         Print() << " Active control based on fuel mass activated."
                 << " Maintaining the flame at " << ctrl_h << " in " << ctrl_flameDir << " direction. \n";
      }
   }
}

void
PeleLM::parseComposition(Vector<std::string> compositionIn,
                         std::string         compositionType,
                         Real               *massFrac)
{
   Real compoIn[NUM_SPECIES] = {0.0};

   // Get species names
   Vector<std::string> specNames;
   EOS::speciesNames(specNames);

   // For each entry in the user-provided composition, parse name and value
   std::string delimiter = ":";
   int specCountIn = compositionIn.size();
   for (int i = 0; i < specCountIn; i++ ) {
      int sep = compositionIn[i].find(delimiter);
      if ( sep == std::string::npos ) {
         Abort("Error parsing '"+compositionIn[i]+"' --> unable to find delimiter :");
      }
      std::string specNameIn = compositionIn[i].substr(0, sep);
      Real value = std::stod(compositionIn[i].substr(sep+1,compositionIn[i].length()));
      int foundIt = 0;
      for (int k = 0; k < NUM_SPECIES; k++ ) {
         if ( specNameIn == specNames[k] ) {
            compoIn[k] = value;
            foundIt = 1;
         }
      }
      if ( !foundIt ) {
         Abort("Error parsing '"+compositionIn[i]+"' --> unable to match to any species name");
      }
   }

   // Ensure that it sums to 1.0:
   Real sum = 0.0;
   for (int k = 0; k < NUM_SPECIES; k++ ) {
      sum += compoIn[k];
   }
   for (int k = 0; k < NUM_SPECIES; k++ ) {
      compoIn[k] /= sum;
   }

   // Fill the massFrac array, convert from mole fraction if necessary
   if ( compositionType == "mass" ) {                // mass
      for (int i = 0; i < NUM_SPECIES; i++ ) {
         massFrac[i] = compoIn[i];
      }
   } else if ( compositionType == "mole" ) {         // mole
      EOS::X2Y(compoIn,massFrac);
   } else {
      Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
   }
}
