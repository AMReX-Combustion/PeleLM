//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
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

#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ErrorList.H>
#include <PeleLM.H>
#include <PeleLM_F.H>
#include <Prob_F.H>
#include <DIFFUSION_F.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_SPACE.H>
#include <AMReX_Interpolater.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Utility.H>
#include <NS_util.H>

#if defined(BL_USE_NEWMECH) || defined(BL_USE_VELOCITY)
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#ifdef USE_SUNDIALS_PP
#include <actual_Creactor.h>
#else
#include <actual_reactor.H> 
#endif

#include <Prob_F.H>
#include <NAVIERSTOKES_F.H>
#include <DERIVE_F.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

static Box stripBox; // used for debugging

#ifdef BL_USE_FLOAT
const Real Real_MIN = FLT_MIN;
const Real Real_MAX = FLT_MAX;
#else
const Real Real_MIN = DBL_MIN;
const Real Real_MAX = DBL_MAX;
#endif

const int  GEOM_GROW   = 1;
const int  PRESS_GROW  = 1;
const int  DIVU_GROW   = 1;
const int  DSDT_GROW   = 1;
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
  bool                  do_active_control;
  bool                  do_active_control_temp;
  Real                  temp_control;
  Real                  crse_dt;
  int                   chem_box_chop_threshold;
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
int  PeleLM::have_trac;
int  PeleLM::have_rhort;
int  PeleLM::Trac;
int  PeleLM::RhoRT;
int  PeleLM::first_spec;
int  PeleLM::last_spec;
int  PeleLM::nspecies;
int  PeleLM::nreactions;
Vector<std::string>  PeleLM::spec_names;
int  PeleLM::floor_species;
int  PeleLM::do_set_rho_to_species_sum;
Real PeleLM::rgas;
Real PeleLM::prandtl;
Real PeleLM::schmidt;
Real PeleLM::constant_thick_val;
int  PeleLM::unity_Le;
Real PeleLM::htt_tempmin;
Real PeleLM::htt_tempmax;
Real PeleLM::htt_hmixTYP;
int  PeleLM::zeroBndryVisc;
int  PeleLM::do_add_nonunityLe_corr_to_rhoh_adv_flux;
int  PeleLM::do_check_divudt;
int  PeleLM::hack_nochem;
int  PeleLM::hack_nospecdiff;
int  PeleLM::hack_noavgdivu;
int  PeleLM::use_tranlib;
Real PeleLM::trac_diff_coef;
Real PeleLM::P1atm_MKS;
bool PeleLM::plot_reactions;
bool PeleLM::plot_consumption;
bool PeleLM::plot_heat_release;
int  PeleLM::cvode_iE;
int  PeleLM::cvode_ncells;
static bool plot_rhoydot;
bool PeleLM::flag_active_control;
Real PeleLM::new_T_threshold;
int  PeleLM::nGrowAdvForcing=1;
bool PeleLM::avg_down_chem;
int  PeleLM::reset_typical_vals_int=-1;
Real PeleLM::typical_Y_val_min=1.e-10;
std::map<std::string,Real> PeleLM::typical_values_FileVals;

std::string                                PeleLM::turbFile;
std::map<std::string, Vector<std::string> > PeleLM::auxDiag_names;

Vector<Real> PeleLM::typical_values;

int PeleLM::sdc_iterMAX;
int PeleLM::num_mac_sync_iter;

int PeleLM::mHtoTiterMAX;
Vector<amrex::Real> PeleLM::mTmpData;

std::string  PeleLM::probin_file = "probin";

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
                         MultiFab& rhohmix)
{
  const Real strt_time = ParallelDescriptor::second();

  const int ngrow  = 0; // We only do this on the valid region
  const int sComp  = std::min(std::min((int)Density,(int)Temp),first_spec);
  const int eComp  = std::max(std::max((int)Density,(int)Temp),first_spec+nspecies-1);
  const int nComp  = eComp - sComp + 1;
  const int sCompR = Density - sComp;
  const int sCompT = Temp - sComp;
  const int sCompY = first_spec - sComp;
  const int sCompH = 0;

  FillPatchIterator fpi(*this,rhohmix,ngrow,time,State_Type,sComp,nComp);
  MultiFab& smf = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox tmp;

    for (MFIter mfi(rhohmix, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      FArrayBox& sfab  = smf[mfi];

      //
      // Convert rho*Y to Y for this operation
      //
      tmp.resize(bx,1);
      tmp.copy(sfab,sCompR,0,1);
      tmp.invert(1.0);
	    
      for (int k = 0; k < nspecies; k++) {
        sfab.mult(tmp,0,sCompY+k,1);
      }
	    
      getHmixGivenTY_pphys(rhohmix[mfi],sfab,sfab,bx,
                                    sCompT,sCompY,sCompH);
      //
      // Convert hmix to rho*hmix
      //
      rhohmix[mfi].mult(sfab,bx,sCompR,sCompH,1);
    }
  }

  if (verbose)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    amrex::Print() << "PeleLM::compute_rhohmix(): lev: " << level << ", time: " << run_time << '\n';
  }
}

void
PeleLM::init_network ()
{
	pphys_network_init(); 
}

void
PeleLM::init_transport (int iEG)
{
	pphys_transport_init(&iEG); 
}

void
PeleLM::init_extern ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  amrex::Print() << "reading extern runtime parameters ... \n" << std::endl;

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
  {
    probin_file_name[i] = probin_file[i];
  }

  plm_extern_init(probin_file_name.dataPtr(),&probin_file_length);
}

void
PeleLM::reactionRateRhoY_pphys(FArrayBox&       RhoYdot,
                             const FArrayBox& RhoY,
                             const FArrayBox& RhoH,
                             const FArrayBox& T,
                             const Box&       box,
                             int              sCompRhoY,
                             int              sCompRhoH,
                             int              sCompT,
                             int              sCompRhoYdot) const
{
    BL_ASSERT(RhoYdot.nComp() >= sCompRhoYdot + nspecies);
    BL_ASSERT(RhoY.nComp() >= sCompRhoY + nspecies);
    BL_ASSERT(RhoH.nComp() > sCompRhoH);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = RhoY.box();
    const Box& mbbx = RhoH.box();
    const Box& mcbx = T.box();
    const Box& mobx = RhoYdot.box();
    
    const Box& ovlp = box & mabx & mbbx & mcbx & mobx;
    if( ! ovlp.ok() ) return;
    
    pphys_RRATERHOY(BL_TO_FORTRAN_BOX(ovlp),
                    BL_TO_FORTRAN_N_ANYD(RhoY,sCompRhoY),
                    BL_TO_FORTRAN_N_ANYD(RhoH,sCompRhoH),
                    BL_TO_FORTRAN_N_ANYD(T,sCompT),
                    BL_TO_FORTRAN_N_ANYD(RhoYdot,sCompRhoYdot) );
}

void
PeleLM::getPGivenRTY_pphys(FArrayBox&       p,
			    const FArrayBox& Rho,
			    const FArrayBox& T,
			    const FArrayBox& Y,
			    const Box&       box,
			    int              sCompR,
			    int              sCompT,
			    int              sCompY,
			    int              sCompP) const
{
    BL_ASSERT(p.nComp() > sCompP);
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+nspecies);
    
    BL_ASSERT(p.box().contains(box));
    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    pphys_PfromRTY(BL_TO_FORTRAN_BOX(box),
		             BL_TO_FORTRAN_N_ANYD(p,sCompP),
		             BL_TO_FORTRAN_N_ANYD(Rho,sCompR),
		             BL_TO_FORTRAN_N_ANYD(T,sCompT),
		             BL_TO_FORTRAN_N_ANYD(Y,sCompY));
}

void
PeleLM::getHmixGivenTY_pphys(FArrayBox&       hmix,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompT,
			      int              sCompY,
			      int              sCompH) const
{
    BL_ASSERT(hmix.nComp() > sCompH);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+nspecies);
    
    BL_ASSERT(hmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    pphys_HMIXfromTY(BL_TO_FORTRAN_BOX(box),
		               BL_TO_FORTRAN_N_ANYD(hmix,sCompH),
		               BL_TO_FORTRAN_N_ANYD(T,sCompT),
		               BL_TO_FORTRAN_N_ANYD(Y,sCompY));
}

void
PeleLM::getHGivenT_pphys(FArrayBox&       h,
			  const FArrayBox& T,
			  const Box&       box,
			  int              sCompT,
			  int              sCompH) const
{
    BL_ASSERT(h.nComp() >= sCompH + nspecies);
    BL_ASSERT(T.nComp() > sCompT);

    BL_ASSERT(h.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    
    pphys_HfromT(BL_TO_FORTRAN_BOX(box),
		           BL_TO_FORTRAN_N_ANYD(h,sCompH),
		           BL_TO_FORTRAN_N_ANYD(T,sCompT));
}

void
PeleLM::getMwmixGivenY_pphys(FArrayBox&       mwmix,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompY,
			      int              sCompMw) const
{
    BL_ASSERT(mwmix.nComp() > sCompMw);
    BL_ASSERT(Y.nComp() >= sCompY+nspecies);

    BL_ASSERT(mwmix.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    pphys_MWMIXfromY(BL_TO_FORTRAN_BOX(box),
		               BL_TO_FORTRAN_N_ANYD(mwmix,sCompMw),
		               BL_TO_FORTRAN_N_ANYD(Y,sCompY));
}

void
PeleLM::getCpmixGivenTY_pphys(FArrayBox&       cpmix,
			       const FArrayBox& T,
			       const FArrayBox& Y,
			       const Box&       box,
			       int              sCompT,
			       int              sCompY,
			       int              sCompCp) const
{
    BL_ASSERT(cpmix.nComp() > sCompCp);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+nspecies);

    BL_ASSERT(cpmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    pphys_CPMIXfromTY(BL_TO_FORTRAN_BOX(box),
                      BL_TO_FORTRAN_N_ANYD(cpmix,sCompCp),
                      BL_TO_FORTRAN_N_ANYD(T,sCompT),
                      BL_TO_FORTRAN_N_ANYD(Y,sCompY));
}


int 
PeleLM::getTGivenHY_pphys(FArrayBox&     T,
                        const FArrayBox& H,
                        const FArrayBox& Y,
                        const Box&       box,
                        int              sCompH,
                        int              sCompY,
                        int              sCompT,
	                const Real&      errMAX)
{
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(H.nComp() > sCompH);
    BL_ASSERT(Y.nComp() >= sCompY+nspecies);
    
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(H.box().contains(box));
    BL_ASSERT(Y.box().contains(box));

    Real solveTOL = (errMAX<0 ? 1.0e-8 : errMAX);
    
    return pphys_TfromHY(BL_TO_FORTRAN_BOX(box),
			                BL_TO_FORTRAN_N_ANYD(T,sCompT),
                         BL_TO_FORTRAN_N_ANYD(H,sCompH),
			                BL_TO_FORTRAN_N_ANYD(Y,sCompY),
			                &solveTOL, &mHtoTiterMAX, mTmpData.dataPtr());
}


int 
PeleLM::getSpeciesIdx(const std::string& spName)
{
    for (int i=0; i<nspecies; i++) {
        if (spName == spec_names[i]) {
            return i;
        }
    }
    return -1;
}

void
PeleLM::getSpeciesNames(Vector<std::string>& spn)
{
    for (int i = 0; i < nspecies; i++) {
        int len = 20;
        Vector<int> int_spec_name(len);
        pphys_get_spec_name(int_spec_name.dataPtr(),&i,&len);
        char* spec_name = new char[len+1];
        for (int j = 0; j < len; j++) {
            spec_name[j] = int_spec_name[j];
        }
        spec_name[len] = '\0';
        //spec_names.push_back(spec_name);
        spn.push_back(spec_name);
        delete [] spec_name;
    }
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
  do_active_control       = false;
  do_active_control_temp  = false;
  temp_control            = -1;
  crse_dt                 = -1;
  chem_box_chop_threshold = -1;
  
  // These two flags below belong to the NavierStokesBase class in IAMR,
  // we set them to 1 by default even if later we do a constant mu and/or lambda
  // Note that constant mu and lambda will be imposed from Constant module in PelePhysics
  variable_vel_visc = 1;
  variable_scal_diff = 1;

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
  PeleLM::have_trac                 = 0;
  PeleLM::have_rhort                = 0;
  PeleLM::Trac                      = -1;
  PeleLM::RhoRT                     = -1;
  PeleLM::first_spec                = -1;
  PeleLM::last_spec                 = -2;
  PeleLM::nspecies                  = 0;
  PeleLM::floor_species             = 1;
  PeleLM::do_set_rho_to_species_sum = 1;
  PeleLM::rgas                      = -1.0;
  PeleLM::prandtl                   = .7;
  PeleLM::schmidt                   = .7;
  PeleLM::constant_thick_val        = -1;
  PeleLM::unity_Le                  = 1;
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
  PeleLM::P1atm_MKS                 = -1.0;
  PeleLM::turbFile                  = "";
  PeleLM::plot_reactions            = false;
  PeleLM::plot_consumption          = true;
  PeleLM::plot_heat_release         = true;
  plot_rhoydot                            = false;
  PeleLM::new_T_threshold           = -1;  // On new AMR level, max change in lower bound for T, not used if <=0
  PeleLM::avg_down_chem             = false;
  PeleLM::reset_typical_vals_int    = -1;
  PeleLM::typical_values_FileVals.clear();

  PeleLM::do_add_nonunityLe_corr_to_rhoh_adv_flux = 1;

  PeleLM::sdc_iterMAX               = 1;
  PeleLM::num_mac_sync_iter         = 1;
  PeleLM::mHtoTiterMAX              = 20;
  PeleLM::cvode_iE                  = 2;
  PeleLM::cvode_ncells              = 1;

  ParmParse pp("ns");

  pp.query("do_diffuse_sync",do_diffuse_sync);
  BL_ASSERT(do_diffuse_sync == 0 || do_diffuse_sync == 1);
  pp.query("do_reflux_visc",do_reflux_visc);
  BL_ASSERT(do_reflux_visc == 0 || do_reflux_visc == 1);
  pp.query("do_active_control",do_active_control);
  pp.query("do_active_control_temp",do_active_control_temp);
  pp.query("temp_control",temp_control);

  if (do_active_control_temp && temp_control <= 0)
    amrex::Error("temp_control MUST be set with do_active_control_temp");

  PeleLM::flag_active_control = do_active_control;

  verbose = pp.contains("v");

  pp.query("divu_ceiling",divu_ceiling);
  BL_ASSERT(divu_ceiling >= 0 && divu_ceiling <= 3);
  pp.query("divu_dt_factor",divu_dt_factor);
  BL_ASSERT(divu_dt_factor>0 && divu_dt_factor <= 1.0);
  pp.query("min_rho_divu_ceiling",min_rho_divu_ceiling);
  if (divu_ceiling) BL_ASSERT(min_rho_divu_ceiling >= 0.0);

  pp.query("htt_tempmin",htt_tempmin);
  pp.query("htt_tempmax",htt_tempmax);

  pp.query("floor_species",floor_species);
  BL_ASSERT(floor_species == 0 || floor_species == 1);

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

  pp.query("sdc_iterMAX",sdc_iterMAX);
  pp.query("num_mac_sync_iter",num_mac_sync_iter);

  pp.query("thickening_factor",constant_thick_val);
  if (constant_thick_val != -1)
  {
    if (verbose)
      amrex::Print() << "PeleLM::read_params: using a constant thickening factor = " 
		     << constant_thick_val << '\n';
  }

  pp.query("do_add_nonunityLe_corr_to_rhoh_adv_flux", do_add_nonunityLe_corr_to_rhoh_adv_flux);
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
  // Check that we are setting correctly default variability/visc for velocities.
  //
  if (variable_vel_visc != 1)
    amrex::Error("PeleLM::read_params() -- must use variable viscosity");
  //
  // Check that we are setting correctly default variability/visc for diffusivity
  //
  if (variable_scal_diff != 1)
    amrex::Error("PeleLM::read_params() -- must use variable scalar diffusivity");
  //
  // Read in scalar value and use it as tracer.
  //
  pp.query("scal_diff_coefs",trac_diff_coef);

  for (int i = 0; i < visc_coef.size(); i++)
    visc_coef[i] = bogus_value;

  // Get some useful amr inputs
  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);

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
    ParmParse pp("peleLM");
    
    // Get boundary conditions
    Vector<std::string> lo_bc_char(BL_SPACEDIM);
    Vector<std::string> hi_bc_char(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc_char,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc_char,0,BL_SPACEDIM);

    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    bool flag_closed_chamber = false;
    for (int dir = 0; dir<BL_SPACEDIM; dir++){
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
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (DefaultGeometry().isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "NavierStokesBase::variableSetUp:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Abort("NavierStokesBase::Initialize()");
                }
            } 
        }
    }

    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (!DefaultGeometry().isPeriodic(dir))
            {
              if (lo_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
              if (hi_bc[dir] == Interior)
              {
                  std::cerr << "NavierStokesBase::variableSetUp:Interior bc in direction "
                            << dir
                            << " but not defined as periodic\n";
                  amrex::Abort("NavierStokesBase::Initialize()");
              }
            }
        }
    } 

    PeleLM::closed_chamber            = 1;
    if (flag_closed_chamber == true){
      PeleLM::closed_chamber            = 0;
    }

    PeleLM::dpdt_factor = 1.0;
    pp.query("dpdt_factor",dpdt_factor);

}

void
PeleLM::Finalize ()
{
  initialized = false;
}

static
void
FabMinMax (FArrayBox& fab,
           const Box& box,
           Real       fmin,
           Real       fmax,
           int        sComp,
           int        nComp)
{
  BL_ASSERT(fab.box().contains(box));
  BL_ASSERT(sComp + nComp <= fab.nComp());

  fab_minmax(BL_TO_FORTRAN_BOX(box),
             BL_TO_FORTRAN_ANYD(fab),
             &fmin, &fmax, &nComp);
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

    mf.copyTo(sub,0,0,mf.nComp(),mf.nGrow());

    if (ShowMF_Check_Nans)
    {
      BL_ASSERT(!sub.contains_nan(box,0,mf.nComp()));
    }
    std::ofstream os;
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
        //                BL_ASSERT(!mf[mfi].contains_nan(mfi.validbox(),0,mf.nComp()));
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
    : lmef(0), lmef_box(_lmef_box), box(_box), min_time(_min_time), max_time(_max_time), max_level(_max_level)
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
    BL_ASSERT(lmef);

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
    BL_ASSERT(lmef_box);

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
PeleLM::center_to_edge_fancy (const FArrayBox& cfab,
                              FArrayBox&       efab,
                              const Box&       ccBox,
                              const Box&       ebox,
                              int              sComp,
                              int              dComp,
                              int              nComp,
                              const Box&       domain,
                              const FPLoc&     bc_lo,
                              const FPLoc&     bc_hi)
{
  //const Box&      ebox = efab.box();
  const IndexType ixt  = ebox.ixType();

  BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));

  int dir = -1;
  for (int d = 0; d < BL_SPACEDIM; d++)
    if (ixt.test(d))
      dir = d;

  BL_ASSERT(amrex::grow(ccBox,-amrex::BASISV(dir)).contains(amrex::enclosedCells(ebox)));
  BL_ASSERT(sComp+nComp <= cfab.nComp() && dComp+nComp <= efab.nComp());
  //
  // Exclude unnecessary cc->ec calcs
  //
  Box ccVBox = ccBox;
  if (bc_lo!=HT_Center)
    ccVBox.setSmall(dir,std::max(domain.smallEnd(dir),ccVBox.smallEnd(dir)));
  if (bc_hi!=HT_Center)
    ccVBox.setBig(dir,std::min(domain.bigEnd(dir),ccVBox.bigEnd(dir)));
  //
  // Shift cell-centered data to edges
  //
  const int isharm = def_harm_avg_cen2edge?1:0;
  cen2edg(ebox.loVect(),ebox.hiVect(),
               ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),cfab.dataPtr(sComp),
               ARLIM(efab.loVect()),ARLIM(efab.hiVect()),efab.dataPtr(dComp),
               &nComp, &dir, &isharm);
  //
  // Fix boundary...i.e. fill-patched data in cfab REALLY lives on edges
  //
  if ( !(domain.contains(ccBox)) )
  {
    if (bc_lo==HT_Edge)
    {
      BoxList gCells = amrex::boxDiff(ccBox,domain);
      if (gCells.ok())
      {
        const int inc = +1;
        FArrayBox ovlpFab;
        for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
             bli != end;
             ++bli)
        {
          if (bc_lo == HT_Edge)
          {
            ovlpFab.resize(*bli,nComp);
            ovlpFab.copy(cfab,sComp,0,nComp);
            ovlpFab.shiftHalf(dir,inc);
            efab.copy(ovlpFab,0,dComp,nComp);
          }
        }
      }
    }
    if (bc_hi==HT_Edge)
    {
      BoxList gCells = amrex::boxDiff(ccBox,domain);
      if (gCells.ok())
      {
        const int inc = -1;
        FArrayBox ovlpFab;
        for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
             bli != end;
             ++bli)
        {
          if (bc_hi == HT_Edge)
          {
            ovlpFab.resize(*bli,nComp);
            ovlpFab.copy(cfab,sComp,0,nComp);
            ovlpFab.shiftHalf(dir,inc);
            efab.copy(ovlpFab,0,dComp,nComp);
          }
        }
      }
    }
  }
}    

void
PeleLM::variableCleanUp ()
{
  NavierStokesBase::variableCleanUp();


  ShowMF_Sets.clear();
  auxDiag_names.clear();
  typical_values.clear();
}

PeleLM::PeleLM ()
{
  if (!init_once_done)
    init_once();

  if (!do_temp)
    amrex::Abort("do_temp MUST be true");

  if (!have_divu)
    amrex::Abort("have_divu MUST be true");

  if (!have_dsdt)
    amrex::Abort("have_dsdt MUST be true");

  // p_amb_old and p_amb_new contain the old-time and new-time
  // pressure at level 0.  For closed chamber problems they change over time.
  // set p_amb_old and new if they haven't been set yet
  // to the value in mod_Fvar_def.F90 set in PROB_F.F90
  // only the coarse level advance and the level 0-1 mac_sync 
  // can modify these later
  if (p_amb_old == -1.0)
  {
    get_pamb(&p_amb_old);
  }
  if (p_amb_new == -1.0)
  {
    get_pamb(&p_amb_new);
  }

  updateFluxReg = false;

  EdgeState              = 0;
  EdgeFlux               = 0;
  SpecDiffusionFluxn     = 0;
  SpecDiffusionFluxnp1   = 0;
#ifdef USE_WBAR
  SpecDiffusionFluxWbar  = 0;
#endif
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

  if (!have_dsdt)
    amrex::Abort("have_dsdt MUST be true");

  // p_amb_old and p_amb_new contain the old-time and new-time
  // pressure at level 0.  For closed chamber problems they change over time.
  // set p_amb_old and new if they haven't been set yet
  // to the value in mod_Fvar_def.F90 set in PROB_F.F90
  // only the coarse level advance and the level 0-1 mac_sync 
  // can modify these later
  if (p_amb_old == -1.0)
  {
    get_pamb(&p_amb_old);
  }
  if (p_amb_new == -1.0)
  {
    get_pamb(&p_amb_new);
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
  const int nEdgeStates = desc_lst[State_Type].nComp();

  mTmpData.resize(mHtoTiterMAX);

  raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nEdgeStates, nGrow)});
  EdgeState = raii_fbs.back()->get();

  raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nEdgeStates, nGrow)});
  EdgeFlux  = raii_fbs.back()->get();
    
  if (nspecies>0 && !unity_Le)
  {
    raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nspecies+3, nGrow)});
    SpecDiffusionFluxn   = raii_fbs.back()->get();

    raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nspecies+3, nGrow)});
    SpecDiffusionFluxnp1 = raii_fbs.back()->get();

#ifdef USE_WBAR
    raii_fbs.push_back(std::unique_ptr<FluxBoxes>{new FluxBoxes(this, nspecies, nGrow)});
    SpecDiffusionFluxWbar = raii_fbs.back()->get();
#endif
    sumSpecFluxDotGradHn.define(grids,dmap,1,0);
    sumSpecFluxDotGradHnp1.define(grids,dmap,1,0);
  }

  for (const auto& kv : auxDiag_names)
  {
      auxDiag[kv.first] = std::unique_ptr<MultiFab>(new MultiFab(grids,dmap,kv.second.size(),0));
      auxDiag[kv.first]->setVal(0);
  }

  // HACK for debugging
  if (level==0)
    stripBox = getStrip(geom);

#ifdef USE_WBAR
  // this will hold the transport coefficients for Wbar
  diffWbar_cc.define(grids,dmap,nspecies,1);
#endif
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

  int tTrac = -1;
  have_trac = isStateVariable("tracer", dummy_State_Type, tTrac);
  AMREX_ALWAYS_ASSERT(tTrac == Trac);
  have_trac = have_trac && State_Type == dummy_State_Type;

  int tRhoRT = -1;
  have_rhort = isStateVariable("RhoRT", dummy_State_Type, tRhoRT);
  AMREX_ALWAYS_ASSERT(tRhoRT == RhoRT);
  have_rhort = have_rhort && State_Type == dummy_State_Type;

  if (!have_temp)
    amrex::Abort("PeleLM::init_once(): RhoH & Temp must both be the state");
    
  if (!have_rhort && verbose)
    amrex::Warning("PeleLM::init_once(): RhoRT being stored in the Tracer slot");
    
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
  BL_ASSERT(Temp > RhoH && RhoH > Density);
  //
  // Here we want to count relative to Density instead of relative
  // to RhoH, so we can put RhoH after the species instead of before.  
  // This logic should work in both cases.
  //
  // Got to make sure nspecies is initialized though
  last_spec  = first_spec + nspecies - 1;
    
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

  set_scal_numb(&density, &Temp, &Trac, &RhoH, &first_spec, &last_spec);
  //
  // Load constants from Fortran module to thickenig factor, etc.
  //
    
  set_ht_adim_common( &constant_thick_val,
                      &prandtl,  &schmidt, &unity_Le);

  //
  // make space for typical values
  //
  typical_values.resize(NUM_STATE,-1); // -ve means don't use for anything
  typical_values[RhoH] = typical_RhoH_value_default;

  Vector<std::string> speciesNames;
  PeleLM::getSpeciesNames(speciesNames);
  ParmParse pp("ht");
  for (int i=0; i<nspecies; ++i) {
    const std::string ppStr = std::string("typValY_") + speciesNames[i];
    if (pp.countval(ppStr.c_str())>0) {
      pp.get(ppStr.c_str(),typical_values_FileVals[speciesNames[i]]);
    }
  }
  std::string otherKeys[4] = {"Temp", "RhoH", "Vel", "Trac"};
  for (int i=0; i<4; ++i) {
    const std::string ppStr(std::string("typVal_")+otherKeys[i]);
    if (pp.countval(ppStr.c_str())>0) {
      pp.get(ppStr.c_str(),typical_values_FileVals[otherKeys[i]]);
    }
  }
  //
  // Get universal gas constant from Fortran.
  //
  rgas = pphys_getRuniversal();
  P1atm_MKS = pphys_getP1atm_MKS();

  if (rgas <= 0.0)
  {
    std::cerr << "PeleLM::init_once(): bad rgas: " << rgas << '\n';
    amrex::Abort();
  }
  if (P1atm_MKS <= 0.0)
  {
    std::cerr << "PeleLM::init_once(): bad P1atm_MKS: " << P1atm_MKS << '\n';
    amrex::Abort();
  }
  //
  // Chemistry.
  //
  int ydot_good = RhoYdot_Type >= 0 && RhoYdot_Type <  desc_lst.size()
                                    && RhoYdot_Type != Divu_Type
                                    && RhoYdot_Type != Dsdt_Type
                                    && RhoYdot_Type != State_Type;


  //
  // This is the minimum number of boxes per MPI proc that I want
  // when chopping up chemistry work.
  //
  pp.query("chem_box_chop_threshold",chem_box_chop_threshold);

  if (chem_box_chop_threshold <= 0)
  {
#ifdef BL_USE_OMP
    chem_box_chop_threshold = 8;
#else
    chem_box_chop_threshold = 4;
#endif
  }
    
  if (!ydot_good)
    amrex::Error("PeleLM::init_once(): need RhoYdot_Type if do_chemistry");

  if (desc_lst[RhoYdot_Type].nComp() < nspecies)
    amrex::Error("PeleLM::init_once(): RhoYdot_Type needs nspecies components");
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
    auxDiag_names["REACTIONS"].resize(nreactions);
    amrex::Print() << "nreactions "<< nreactions << '\n';
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
    
  pp.query("new_T_threshold",new_T_threshold);

#ifdef BL_COMM_PROFILING
  auxDiag_names["COMMPROF"].resize(3);
  auxDiag_names["COMMPROF"][0] = "mpiRank";
  auxDiag_names["COMMPROF"][1] = "proximityRank";
  auxDiag_names["COMMPROF"][2] = "proximityOrder";
#endif

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
#pragma omp parallel
#endif
    for (MFIter mfi(rYdot_old,true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      FArrayBox& ry1 = rYdot_old[mfi];
      FArrayBox& ry2 = rYdot_new[mfi];
      const FArrayBox& S1 = S_old[mfi];
      const FArrayBox& S2 = S_new[mfi];
      for (int i=0; i<nspecies; ++i)
      {
        ry1.mult(S1,bx,Density,i,1);
        ry2.mult(S2,bx,Density,i,1);
      }
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
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      // read in title line
      std::getline(is, line);

      is >> p_amb_old;
      p_amb_new = p_amb_old;
  }
}

void
PeleLM::set_typical_values(bool is_restart)
{
  if (level==0)
  {
    const int nComp = typical_values.size();

    if (is_restart) 
    {
      BL_ASSERT(nComp==NUM_STATE);

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

      auto scaleMax = VectorMax(Slevs,FabArrayBase::mfiter_tile_size,Density,NUM_STATE-BL_SPACEDIM,0);
      auto scaleMin = VectorMin(Slevs,FabArrayBase::mfiter_tile_size,Density,NUM_STATE-BL_SPACEDIM,0);      
      auto velMaxV = VectorMaxAbs(Slevs,FabArrayBase::mfiter_tile_size,0,BL_SPACEDIM,0);

      auto velMax = *max_element(std::begin(velMaxV), std::end(velMaxV));

      for (int i=0; i<NUM_STATE - BL_SPACEDIM; ++i) {
        typical_values[i + BL_SPACEDIM] = std::abs(scaleMax[i] - scaleMin[i]);
		  if ( typical_values[i + BL_SPACEDIM] <= 0.0)
			  typical_values[i + BL_SPACEDIM] = 0.5*std::abs(scaleMax[i] + scaleMin[i]);
      }

      for (int i=0; i<BL_SPACEDIM; ++i) {
        typical_values[i] = velMax;
      }

      AMREX_ALWAYS_ASSERT(typical_values[Density] > 0);
      typical_values[RhoH] = typical_values[RhoH] / typical_values[Density];
      for (int i=0; i<nspecies; ++i) {
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
        else if (it->first == "Trac")
          typical_values[Trac] = it->second;
        else if (it->first == "Vel")
        {
          for (int d=0; d<BL_SPACEDIM; ++d)
            typical_values[d] = std::max(it->second,typical_Y_val_min);
        }
      }
    }

    amrex::Print() << "Typical vals: " << '\n';
    amrex::Print() << "\tVelocity: ";
    for (int i=0; i<BL_SPACEDIM; ++i) {
      amrex::Print() << typical_values[i] << ' ';
    }
    amrex::Print() << '\n';
    amrex::Print() << "\tDensity: " << typical_values[Density] << '\n';
    amrex::Print() << "\tTemp:    " << typical_values[Temp]    << '\n';
    amrex::Print() << "\tRhoH:    " << typical_values[RhoH]    << '\n';
    for (int i=0; i<nspecies; ++i)
      {
        amrex::Print() << "\tY_" << spec_names[i] << ": " << typical_values[first_spec+i] <<'\n';
      }
  }
}

void
PeleLM::reset_typical_values (const MultiFab& S)
{
  //
  // NOTE: Assumes that this level has valid data everywhere.
  //
  const int nComp = typical_values.size();

  BL_ASSERT(nComp == S.nComp());

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
        else if (it->first == "Trac")
          typical_values[Trac] = it->second;
        else if (it->first == "Vel")
        {
          for (int d=0; d<BL_SPACEDIM; ++d)
            typical_values[d] = std::max(it->second,typical_Y_val_min);
        }
      }
    }
  }

  amrex::Print() << "New typical vals: " << '\n';
  amrex::Print() << "\tVelocity: ";
  for (int i=0; i<BL_SPACEDIM; ++i) {
    amrex::Print() << typical_values[i] << ' ';
  }
  amrex::Print() << '\n';
  amrex::Print() << "\tDensity:  " << typical_values[Density] << '\n';
  amrex::Print() << "\tTemp:     " << typical_values[Temp]    << '\n';
  amrex::Print() << "\tRhoH:     " << typical_values[RhoH]    << '\n';
  for (int i=0; i<nspecies; ++i)
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

  const int   n_grow   = 1;
  const Real  cur_time = state[State_Type].curTime();
  const Real* dx       = geom.CellSize();
  MultiFab*   dsdt     = getDsdt(0,cur_time);
  MultiFab*   divu     = getDivCond(0,cur_time);

  FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
  MultiFab& Umf=U_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(min:divu_dt)
#endif
{
  for (MFIter mfi(Umf,true); mfi.isValid();++mfi)
  {
    Real dt;
    const Box& bx = mfi.tilebox();
    const FArrayBox& dsdt_mf   =  (*dsdt)[mfi];
    const FArrayBox& divu_mf   =  (*divu)[mfi];
    const FArrayBox& U   = Umf[mfi];
    const FArrayBox& Rho = rho_ctime[mfi];
    const FArrayBox& vol = volume[mfi];
    const FArrayBox& areax = area[0][mfi];
    const FArrayBox& areay = area[1][mfi];
#if (AMREX_SPACEDIM==3) 
    const FArrayBox& areaz = area[2][mfi];
#endif

    est_divu_dt(divu_ceiling, &divu_dt_factor, dx,
                BL_TO_FORTRAN_ANYD(divu_mf),
                BL_TO_FORTRAN_ANYD(dsdt_mf),
                BL_TO_FORTRAN_ANYD(Rho),
                BL_TO_FORTRAN_ANYD(U),
                BL_TO_FORTRAN_ANYD(vol),
                BL_TO_FORTRAN_ANYD(areax),
                BL_TO_FORTRAN_ANYD(areay),
#if (AMREX_SPACEDIM==3) 
                BL_TO_FORTRAN_ANYD(areaz), 
#endif 
                BL_TO_FORTRAN_BOX(bx), &dt, &min_rho_divu_ceiling);

    divu_dt = std::min(divu_dt,dt);

  }
}

  delete divu;
  delete dsdt;

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

  const int   n_grow    = 1;
  const Real  cur_time  = state[State_Type].curTime();
  const Real* dx        = geom.CellSize();
  MultiFab*   dsdt      = getDsdt(0,cur_time);
  MultiFab*   divu      = getDivCond(0,cur_time);

  FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
  MultiFab& Umf=U_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Umf,true); mfi.isValid();++mfi)
  {
    const Box&    grdbx  = mfi.tilebox();
    const int        i   = mfi.index();
    FArrayBox&       U   = Umf[mfi];
    const FArrayBox& Rho = rho_ctime[mfi];

    check_divu_dt( divu_ceiling, &divu_dt_factor, dx,
                   BL_TO_FORTRAN_ANYD((*divu)[mfi]),
                   BL_TO_FORTRAN_ANYD((*dsdt)[mfi]),
                   BL_TO_FORTRAN_ANYD(Rho),
                   BL_TO_FORTRAN_ANYD(U),
                   BL_TO_FORTRAN_ANYD(volume[i]),
                   BL_TO_FORTRAN_ANYD(area[0][i]),
                   BL_TO_FORTRAN_ANYD(area[1][i]),
#if (AMREX_SPACEDIM==3)
                   BL_TO_FORTRAN_ANYD(area[2][i]),
#endif
                   BL_TO_FORTRAN_BOX(grdbx),
                   &dt, &min_rho_divu_ceiling);
  }

  delete dsdt;
  delete divu;
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
  int         ns       = NUM_STATE - BL_SPACEDIM;
  const Real* dx       = geom.CellSize();
  MultiFab&   S_new    = get_new_data(State_Type);
  MultiFab&   P_new    = get_new_data(Press_Type);
  const Real  cur_time = state[State_Type].curTime();

#ifdef BL_USE_NEWMECH
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
  int nspecies;
  pphys_get_num_spec(&nspecies);
  Vector<std::string> names;
  PeleLM::getSpeciesNames(names);
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
  int idSpec = idX + BL_SPACEDIM;

  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    amrData.FillVar(S_new, level, plotnames[idX+i], Xvel+i);
    amrData.FlushGrids(idX+i);
  }
  amrData.FillVar(S_new, level, plotnames[idT], Temp);
  amrData.FlushGrids(idT);

  for (int i = 0; i < nspecies; i++)
  {
    amrData.FillVar(S_new, level, plotnames[idSpec+i], first_spec+i);
    amrData.FlushGrids(idSpec+i);
  }

  if (verbose) amrex::Print() << "initData: finished init from pltfile" << '\n';
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter snewmfi(S_new,true); snewmfi.isValid(); ++snewmfi)
  {
    //BL_ASSERT(grids[snewmfi.index()] == snewmfi.validbox());

    P_new[snewmfi].setVal(0);

    const Box& vbx = snewmfi.tilebox();
    RealBox    gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());

#ifdef BL_USE_NEWMECH
    init_data_new_mech (&level, &cur_time,
                        BL_TO_FORTRAN_BOX(vbx), &ns,
                        S_new[snewmfi].dataPtr(Xvel),
                        BL_TO_FORTRAN_N_ANYD(S_new[snewmfi],AMREX_SPACEDIM),
                        BL_TO_FORTRAN_ANYD(P_new[snewmfi]),
                        dx, AMREX_ZFILL(gridloc.lo()), AMREX_ZFILL(gridloc.hi()) );
#else
    init_data (&level, &cur_time,
               BL_TO_FORTRAN_BOX(vbx), &ns,
               S_new[snewmfi].dataPtr(Xvel),
               BL_TO_FORTRAN_N_ANYD(S_new[snewmfi],AMREX_SPACEDIM),
               BL_TO_FORTRAN_ANYD(P_new[snewmfi]),
               dx, AMREX_ZFILL(gridloc.lo()), AMREX_ZFILL(gridloc.hi()) );
#endif
  }

  showMFsub("1D",S_new,stripBox,"1D_S",level);

#ifdef BL_USE_VELOCITY
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
    for (int i = 0; i < BL_SPACEDIM; i++)
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
#endif /*BL_USE_VELOCITY*/

  make_rho_prev_time();
  make_rho_curr_time();
  //
  // Initialize other types.
  //
  initDataOtherTypes();
  //
  // Initialize divU and dSdt.
  //
  if (have_divu)
  {
    const Real dt       = 1.0;
    const Real dtin     = -1.0; // Dummy value denotes initialization.
    const Real tnp1 = state[Divu_Type].curTime();
    MultiFab&  Divu_new = get_new_data(Divu_Type);

    state[State_Type].setTimeLevel(tnp1,dt,dt);

    calcDiffusivity(tnp1);
#ifdef USE_WBAR
    calcDiffusivity_Wbar(tnp1);
#endif

    calc_divu(tnp1,dtin,Divu_new);

    if (have_dsdt)
      get_new_data(Dsdt_Type).setVal(0);
  }

  if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
  {
    get_new_data(Dpdt_Type).setVal(0);
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
  // to be consistent with Strang code, compute enthalpy with the EOS
  const Real tnp1  = state[State_Type].curTime();
  {
    MultiFab rhoh(grids,dmap,1,0);
    compute_rhohmix(tnp1,rhoh);
    get_new_data(State_Type).copy(rhoh,0,RhoH,1);
  }
  //
  // Assume that by now, S_new has "good" data
  //
  MultiFab& R = get_new_data(RhoYdot_Type);

  // It is our current belief that at the beginning of a simulation we want omegadot=0
  // This is how the Strang code does it
  //    compute_instantaneous_reaction_rates(R,get_new_data(State_Type),nGrow);
  get_new_data(RhoYdot_Type).setVal(0);

  showMFsub("1D",R,stripBox,"1D_dd_idot_R",level);

  get_new_data(FuncCount_Type).setVal(1);

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
    R.setVal(0);
    R.setBndry(0,0,nspecies);
    return;
  }

  const Real strt_time = ParallelDescriptor::second();

  BL_ASSERT((nGrow==0)  ||  (how == HT_ZERO_GROW_CELLS) || (how == HT_EXTRAP_GROW_CELLS));

  if ((nGrow>0) && (how == HT_ZERO_GROW_CELLS))
  {
    R.setBndry(0,0,nspecies);
  }
    
  const int sCompRhoH    = RhoH;
  const int sCompRhoY    = first_spec;
  const int sCompT       = Temp;
  const int sCompRhoYdot = 0;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const FArrayBox& rhoY = S[mfi];
    const FArrayBox& rhoH = S[mfi];
    const FArrayBox& T    = S[mfi];
    const Box& box = mfi.tilebox();
    FArrayBox& rhoYdot = R[mfi];

    reactionRateRhoY_pphys(rhoYdot,rhoY,rhoH,T,box,sCompRhoY,sCompRhoH,sCompT,sCompRhoYdot);

  }
    
  if ((nGrow>0) && (how == HT_EXTRAP_GROW_CELLS))
  {
    R.FillBoundary(0,nspecies, geom.periodicity());
    BL_ASSERT(R.nGrow() == 1);
    Extrapolater::FirstOrderExtrap(R, geom, 0, nspecies);
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
  // Get best ydot data.
  //
  MultiFab& Ydot = get_new_data(RhoYdot_Type);
  {
    FillPatchIterator fpi(*oldht,Ydot,Ydot.nGrow(),tnp1,RhoYdot_Type,0,nspecies);
    const MultiFab& mf_fpi = fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif  
    for (MFIter mfi(mf_fpi,true); mfi.isValid(); ++mfi)
    {
      const Box& vbx  = mfi.tilebox();
      const FArrayBox& pfab = mf_fpi[mfi];
       Ydot[mfi].copy(pfab,vbx,0,vbx,0,nspecies);
    }
  }
  
  RhoH_to_Temp(get_new_data(State_Type));

  MultiFab& FuncCount = get_new_data(FuncCount_Type);
  {
    FillPatchIterator fpi(*oldht,FuncCount,FuncCount.nGrow(),tnp1,FuncCount_Type,0,1);
    const MultiFab& mf_fpi = fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel
#endif  
      for (MFIter mfi(mf_fpi,true); mfi.isValid(); ++mfi)
    {
      const Box& vbx  = mfi.tilebox();
      const FArrayBox& pfab = mf_fpi[mfi];
      FuncCount[mfi].copy(pfab,vbx,0,vbx,0,1);
    }
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
  FillCoarsePatch(get_new_data(RhoYdot_Type),0,tnp1,RhoYdot_Type,0,nspecies);

  RhoH_to_Temp(get_new_data(State_Type));
  get_new_data(State_Type).setBndry(1.e30);
  showMF("sdc",get_new_data(State_Type),"sdc_new_state",level);

  if (new_T_threshold>0)
  {
    MultiFab& crse = coarser.get_new_data(State_Type);
    MultiFab& fine = get_new_data(State_Type);

    RhoH_to_Temp(crse,0); // Make sure T is current
    Real min_T_crse = crse.min(Temp);
    Real min_T_fine = min_T_crse * std::min(1.0, new_T_threshold);

    const int* ratio = crse_ratio.getVect();
    int Tcomp = (int)Temp;
    int Rcomp = (int)Density;
    int n_tmp = std::max( (int)Density, std::max( (int)Temp, std::max( last_spec, RhoH) ) );
    Vector<Real> tmp(n_tmp);
    int num_cells_hacked = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine,true); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = fine[mfi];
      const Box& box = mfi.tilebox();
      num_cells_hacked += conservative_T_floor(BL_TO_FORTRAN_BOX(box),
                                               BL_TO_FORTRAN_ANYD(fab),
                                               &min_T_fine, &Tcomp, &Rcomp,
                                               &first_spec, &last_spec, &RhoH,
                                               ratio, tmp.dataPtr(), &n_tmp);
    }

    ParallelDescriptor::ReduceIntSum(num_cells_hacked);

    if (num_cells_hacked > 0) {

      Real old_min = fine.min(Temp);
      RhoH_to_Temp(get_new_data(State_Type));
      Real new_min = fine.min(Temp);

      if (verbose) {
          amrex::Print() << "...level data adjusted to reduce new extrema (" << num_cells_hacked
                         << " cells affected), new min = " << new_min 
                         << " (old min = " << old_min << ")" << '\n';
      }
    }
  }
  showMF("sdc",get_new_data(State_Type),"sdc_new_state_floored",level);

  FillCoarsePatch(get_new_data(FuncCount_Type),0,tnp1,FuncCount_Type,0,1);
}

void
PeleLM::post_timestep (int crse_iteration)
{
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

      amrex::average_down(Ydot_fine, Ydot_crse, flev.geom, clev.geom,
                           0, Ndiag, parent->refRatio(i-1));
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

  if (do_active_control)
  {
    int usetemp = 0;
    active_control(&dummy,&dummy,&crse_dt,&MyProc,&step,&is_restart,&usetemp);
  }
  else if (do_active_control_temp)
  {
    int usetemp = 1;
    active_control(&dummy,&dummy,&crse_dt,&MyProc,&step,&is_restart,&usetemp);
  }
}

void
PeleLM::postCoarseTimeStep (Real cumtime)
{
  //
  // postCoarseTimeStep() is only called by level 0.
  //
  BL_ASSERT(level == 0);
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
  if (do_set_rho_to_species_sum)
  {
    const int nGrow = 0;
    if (parent->levelSteps(0)>0 && level>lbase)
      set_rho_to_species_sum(get_new_data(State_Type),0,nGrow,0);
  }
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
  //
  // Load typical values for each state component
  //
  set_typical_values(false);

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

      amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, crse_lev.geom,
                           0, 1, crse_lev.fine_ratio);
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

  amrex::Print() << "doing num_divu_iters = " << num_divu_iters << '\n';

  for (int iter = 0; iter < init_divu_iter; ++iter)
  {
    //
    // Update species destruction rates in each level but not state.
    //
    if (nspecies > 0)
    {
      for (int k = 0; k <= finest_level; k++)
      {
        MultiFab& S_new = getLevel(k).get_new_data(State_Type);
        //
        // Don't update S_new in this strang_chem() call ...
        //
        MultiFab S_tmp(S_new.boxArray(),S_new.DistributionMap(),S_new.nComp(),0);

        MultiFab Forcing_tmp(S_new.boxArray(),S_new.DistributionMap(),nspecies+1,0);
        Forcing_tmp.setVal(0);

        getLevel(k).advance_chemistry(S_new,S_tmp,dt_save[k]/2.0,Forcing_tmp,0);
      }
    }
    //
    // Recompute the velocity to obey constraint with chemistry and
    // divqrad and then average that down.
    //
    if (nspecies > 0)
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

          amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, crse_lev.geom,
                               0, 1, crse_lev.fine_ratio);
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

      projector->initialVelocityProject(0,divu_time,havedivu);

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

        amrex::average_down(S_fine, S_crse, fine_lev.geom, crse_lev.geom,
                             Xvel, BL_SPACEDIM, getLevel(k).fine_ratio);
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
    int MyProc  = ParallelDescriptor::MyProc();
    int step    = parent->levelSteps(0);
    int is_restart = 0;

    if (do_active_control)
    {
      Real fuelmass = 0.0;
      std::string fuel = "rho.Y(" + fuelName + ")";
      for (int lev = 0; lev <= finest_level; lev++)
        fuelmass += getLevel(lev).volWgtSum(fuel,tnp1);
        
      if (verbose) amrex::Print() << " FUELMASS= " << fuelmass;
        
      int usetemp = 0;
      active_control(&fuelmass,&tnp1,&crse_dt,&MyProc,&step,&is_restart,&usetemp);      
    }
    else if (do_active_control_temp)
    {
      const int   DM     = BL_SPACEDIM-1;
      const Real* dx     = geom.CellSize();
      const Real* problo = geom.ProbLo();
      Real        hival  = geom.ProbHi(DM);
      MultiFab&   mf     = get_new_data(State_Type);

      FillPatchIterator Tfpi(*this,mf,1,tnp1,State_Type,Temp,1);
      MultiFab& Tmf=Tfpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif  
      for (MFIter mfi(Tmf,true); mfi.isValid();++mfi)
      {
        const FArrayBox& fab  = Tmf[mfi];
        const Box&       vbox = mfi.tilebox();

        for (IntVect iv = vbox.smallEnd(); iv <= vbox.bigEnd(); vbox.next(iv))
        {
          const Real T_hi = fab(iv);

          if (T_hi > temp_control)
          {
            const Real hi = problo[DM] + (iv[DM] + .5) * dx[DM];

            if (hi < hival)
            {
              hival = hi;

              IntVect lo_iv = iv;
              
              lo_iv[DM] -= 1;

              const Real T_lo = fab(lo_iv);

              if (T_lo < temp_control)
              {
                const Real lo    = problo[DM] + (lo_iv[DM] + .5) * dx[DM];
                const Real slope = (T_hi - T_lo) / (hi - lo);
                  
                hival = (temp_control - T_lo) / slope + lo;
              }
            }
          }
        }
      }

      ParallelDescriptor::ReduceRealMin(hival);
      
      if (verbose) Print() << " HIVAL= " << hival;
        
      int usetemp = 1;      
      active_control(&hival,&tnp1,&crse_dt,&MyProc,&step,&is_restart,&usetemp);      
    }
    else
    {
      Real fuelmass = 0.0;
      std::string fuel = "rho.Y(" + fuelName + ")";
      for (int lev = 0; lev <= finest_level; lev++) {
        fuelmass += getLevel(lev).volWgtSum(fuel,tnp1);
      }
      if (verbose) amrex::Print() << " FUELMASS= " << fuelmass;
    }
  }

  if (getSpeciesIdx(productName) >= 0)
  {
    Real productmass = 0.0;
    std::string product = "rho.Y(" + productName + ")";
    for (int lev = 0; lev <= finest_level; lev++) {
      productmass += getLevel(lev).volWgtSum(product,tnp1);
    }	  
    if (verbose) Print() << " PRODUCTMASS= " << productmass;
  }
  if (verbose) Print() << '\n';

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
    for (int n=0; n<BL_SPACEDIM; ++n) {
      std::string str = (n==0 ? "xvel" : (n==1 ? "yvel" : "zvel") );
      Print() << "TIME= " << tnp1  << " min,max "<< str << "  = " << Smin[Xvel+n] << ", " << Smax[Xvel+n] << '\n';
    }

    if (nspecies > 0)
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

        amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, crse_lev.geom,
                             0, 1, crse_lev.fine_ratio);
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

        amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, crse_lev.geom,
                             0, 1, crse_lev.fine_ratio);
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
  PeleLM&   fine_lev  = getLevel(level+1);
  MultiFab&       S_crse    = get_new_data(State_Type);
  MultiFab&       S_fine    = fine_lev.get_new_data(State_Type);

  amrex::average_down(S_fine, S_crse, fine_lev.geom, geom,
                       0, S_crse.nComp(), fine_ratio);
  //
  // Fill rho_ctime at the current and finer levels with the correct data.
  //
  for (int lev = level; lev <= parent->finestLevel(); lev++)
  {
    getLevel(lev).make_rho_curr_time();
  }
  //
  // Reset the temperature
  //
  RhoH_to_Temp(S_crse);
  //
  // Now average down pressure over time n-(n+1) interval.
  //
  MultiFab&       P_crse      = get_new_data(Press_Type);
  MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
  MultiFab&       P_fine_avg  = fine_lev.p_avg;
  MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
  const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();
  const DistributionMapping& P_fdmap = fine_lev.state[Press_Type].DistributionMap();

  BoxArray crse_P_fine_BA = P_fgrids;

  crse_P_fine_BA.coarsen(fine_ratio);

  MultiFab crse_P_fine(crse_P_fine_BA,P_fdmap,1,0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(crse_P_fine,true); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();

    injectDown(bx,crse_P_fine[mfi],P_fine[mfi],fine_ratio);
  }

  P_crse.copy(crse_P_fine, parent->Geom(level).periodicity());  // Parallel copy

  crse_P_fine.clear();
  //
  // Next average down divu and dSdT at new time.
  //
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
            
    amrex::average_down(Divu_fine, Divu_crse, fine_lev.geom, geom,
                         0, Divu_crse.nComp(), fine_ratio);
  }

  if (hack_noavgdivu)
  {
    //
    //  Now that have good divu, recompute dsdt (don't avgDown,
    //   since that will give a very different value, and screw up mac)
    //
    StateData& dsdtSD = get_state_data(Dsdt_Type);// should be const...
    MultiFab&  dsdt   = dsdtSD.newData();

    if (get_state_data(Divu_Type).hasOldData())
    {
      const Real time = dsdtSD.curTime();
      const Real dt   = time - dsdtSD.prevTime();
      calc_dsdt(time,dt,dsdt);
    }
    else
    {
      dsdt.setVal(0.0);
    }
  }
  else
  {
    MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
    MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
            
    amrex::average_down(Dsdt_fine, Dsdt_crse, fine_lev.geom, geom,
                         0, Dsdt_crse.nComp(), fine_ratio);
  }

  if (verbose > 1)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    amrex::Print() << "PeleLM::avgDown(): lev: " << level << ", time: " << run_time << '\n';
  }
}

void
PeleLM::differential_diffusion_update (MultiFab& Force,
                                       int       FComp,
                                       MultiFab& Dnew,
                                       int       DComp,
                                       MultiFab& DDnew)
{
  BL_PROFILE_REGION_START("R::HT::differential_diffusion_update()");
  BL_PROFILE("HT::differential_diffusion_update()");
  //
  // Recompute the D[Nspec+1] = Div(Fi.Hi) using
  // Fi.Hi based on solution here.
  //
  BL_ASSERT(Force.boxArray() == grids);
  BL_ASSERT(FComp+Force.nComp()>=nspecies+1);
  BL_ASSERT(Dnew.boxArray() == grids);
  BL_ASSERT(DDnew.boxArray() == grids);
  BL_ASSERT(DComp+Dnew.nComp()>=nspecies+1);

  const Real strt_time = ParallelDescriptor::second();

  if (hack_nospecdiff)
  {
    amrex::Error("differential_diffusion_update: hack_nospecdiff not implemented");
  }
  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);
  Real prev_time = state[State_Type].prevTime();
  Real curr_time = state[State_Type].curTime();
  Real dt = curr_time - prev_time;
  BL_ASSERT(dt>0);
  const int nGrow = 0;

  MultiFab::Copy(S_new,S_old,first_spec,first_spec,nspecies+1,0);
  int rho_flag = 2;
  MultiFab rhonph;
  Diffusion::SolveMode solve_mode = Diffusion::ONEPASS;
  MultiFab **betanp1, **betan = 0; // Will not need betan since time-explicit pieces computed above
  FluxBoxes fb_betap1(this, nspecies+2, nGrow);
  betanp1 = fb_betap1.get();
  getDiffusivity(betanp1, curr_time, first_spec, 0, nspecies+1); // species (rhoD) and RhoH (lambda/cp)
  getDiffusivity(betanp1, curr_time, Temp, nspecies+1, 1); // temperature (lambda)
  //
  // Diffuse RhoY (and RhoH, unless theta_enthalpy>0 -- if so, will modify the source term for
  //   RhoH with the results of the RhoY diffusion)
  //
  MultiFab* alpha = 0; // Never need alpha for RhoY, RhoH
  int alphaComp = 0;
  int nComp = nspecies+1;
  for (int sigma = 0; sigma < nComp; ++sigma)
  {
    int betaComp = sigma;
    const int state_ind = first_spec + sigma;
    bool add_old_time_divFlux = false; // indicate that the rhs contains the time-explicit diff terms already
    diffusion->diffuse_scalar(dt,state_ind,1.0,rhonph,rho_flag,
			      SpecDiffusionFluxn,SpecDiffusionFluxnp1,sigma,&Force,sigma,alpha,
			      alphaComp,betan,betanp1,betaComp,solve_mode,add_old_time_divFlux);
  }

#ifdef USE_WBAR
  // add lagged grad Wbar fluxes (SpecDiffusionFluxWbar) to time-advanced 
  // species diffusion fluxes (SpecDiffusionFluxnp1)
  for (int d=0; d<BL_SPACEDIM; ++d)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*SpecDiffusionFluxWbar[d],true); mfi.isValid(); ++mfi)
    {
      const Box& ebox = mfi.tilebox();
      (*SpecDiffusionFluxnp1[d])[mfi].plus((*SpecDiffusionFluxWbar[d])[mfi],ebox,0,0,nspecies);
    }
  }
#endif

  //
  // Modify/update new-time fluxes to ensure sum of species fluxes = 0
  //
  adjust_spec_diffusion_fluxes(curr_time);
  //
  // AJN FLUXREG
  // We have just performed the correction diffusion solve for Y_m and h
  // If updateFluxReg=T, we update VISCOUS flux registers:
  //   ADD Gamma_{m,AD}^(k+1),
  //   ADD lambda^(k)/cp^(k) grad h_AD^(k+1).
  //
  if (do_reflux && updateFluxReg)
  {
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
      if (level > 0)
      {
	getViscFluxReg().FineAdd(*SpecDiffusionFluxnp1[d],d,
				 0,first_spec,nspecies+1,dt);
      }
      if (level < parent->finestLevel())
      {
	getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,
						    0,first_spec,nspecies+1,-dt,
						    FluxRegister::ADD);
      }
    }
  }

  flux_divergence(Dnew,DComp,SpecDiffusionFluxnp1,0,nComp,-1);
    
  if (verbose)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    amrex::Print() << "PeleLM::differential_diffusion_update(): lev: " << level 
		   << ", time: " << run_time << '\n';
  }

  BL_PROFILE_REGION_STOP("R::HT::differential_diffusion_update()");
}

void
PeleLM::adjust_spec_diffusion_fluxes (Real time)
{
  //
  // In this function we explicitly adjust the species diffusion fluxes so that their sum
  // is zero.  The fluxes are class member data, either SpecDiffusionFluxn or 
  // SpecDiffusionFluxnp1, depending on time
  //
  const Real      strt_time = ParallelDescriptor::second();
  const TimeLevel whichTime = which_time(State_Type,time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    

  MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

  MultiFab& S = get_data(State_Type,time);
  //
  // Fill grow cells in the state for RhoY, Temp
  // For Dirichlet physical boundary grow cells, this data will live on the cell face, otherwise
  // it will live at cell centers.
  //
  const int nGrow = 1;
  BL_ASSERT(S.nGrow()>=nGrow);
  FillPatchIterator Tfpi(*this,S,nGrow,time,State_Type,Temp,1),
         Yfpi(*this,S,nGrow,time,State_Type,first_spec,nspecies);

  MultiFab& Tmf = Tfpi.get_mf();
  MultiFab& Ymf = Yfpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Tmf,true); mfi.isValid();++mfi)
  {
    const Box& vbox = mfi.tilebox();
    FArrayBox& fab = S[mfi];
    const FArrayBox& Tpfab = Tmf[mfi];
    const FArrayBox& Ypfab = Ymf[mfi];
    BoxList gcells = amrex::boxDiff(Box(vbox).grow(nGrow),vbox);
    for (BoxList::const_iterator it = gcells.begin(), end = gcells.end(); it != end; ++it)
    {
      const Box& gbox = *it;
      fab.copy(Tpfab,gbox,0,gbox,Temp,1);
      fab.copy(Ypfab,gbox,0,gbox,first_spec,nspecies);
    }
  }
  //
  // Get boundary info for Y (assume all Ys have the same boundary type.
  //
  const BCRec& Ybc = get_desc_lst()[State_Type].getBC(first_spec);
  //
  // The following REPAIR_FLUX routine modifies the fluxes of all the species
  // to ensure that they sum to zero.  It requires the RhoY on valid + 1 grow (and, at least
  // as of this writing actually used the values of RhoY on edges that it gets by arithmetic
  // averaging.
  //
  const Box& domain = geom.Domain();
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {

    FArrayBox& Y = S[mfi];
    int sCompY=first_spec;

    for (int d =0; d < BL_SPACEDIM; ++d)
    {
      FArrayBox& f = (*flux[d])[mfi];
      const Box& ebox   = mfi.nodaltilebox(d);
      const Box& edomain = amrex::surroundingNodes(domain,d);
      repair_flux(BL_TO_FORTRAN_BOX(ebox),
                  BL_TO_FORTRAN_BOX(edomain),
                  BL_TO_FORTRAN_ANYD(f), 
                  BL_TO_FORTRAN_N_ANYD(Y,sCompY),
                  &d, Ybc.vect());
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
PeleLM::compute_enthalpy_fluxes (Real                   time,
                                 const MultiFab* const* beta)
{
  BL_ASSERT(beta && beta[0]->nComp() == nspecies+2);

  const Real      strt_time = ParallelDescriptor::second();
  const TimeLevel whichTime = which_time(State_Type,time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    

  MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

  MultiFab& S = get_data(State_Type,time);

  MultiFab& sumSpecFluxDotGradH = (whichTime == AmrOldTime) ? sumSpecFluxDotGradHn : sumSpecFluxDotGradHnp1;
  //
  //  Compute species enthalpy on the edges, and increment the heat flux with the 
  //  flux of enthalpy due to species diffusion.  While here, we also compute the Fi.Grad(Hi) term
  //  required for the temperature equation.  Both the fluxes and the Fi.Grad(Hi) terms are stored
  //  in the class data.
  //
  const Box&   domain = geom.Domain();
  const BCRec& Tbc    = get_desc_lst()[State_Type].getBC(Temp);
  //
  // Fill ghost cells for rhoY and Temp.
  //
  FillPatchIterator rYfpi(*this,S,1,time,State_Type,first_spec,nspecies),
         Tfpi(*this,S,1,time,State_Type,Temp,1);
         
  MultiFab& Tmf = Tfpi.get_mf();
  MultiFab& rYmf = rYfpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(Tmf,true); mfi.isValid();++mfi)
  {
    const Box& vbx  = mfi.tilebox();
    const FArrayBox& Tfab = Tmf[mfi];
    const FArrayBox& rYfab = rYmf[mfi];

    S[mfi].copy(rYfab,vbx,0,vbx,first_spec,nspecies);
    S[mfi].copy(Tfab,vbx,0,vbx,Temp,1);
  }
}

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox ftmp[BL_SPACEDIM];

    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      int              FComp    = 0;
      int              TComp    = Temp;
      int              RhoYComp = first_spec;
      int              dComp    = 0;
      FArrayBox&       fh       = sumSpecFluxDotGradH[mfi];            
      FArrayBox&       T        = S[mfi];
      FArrayBox&       RhoY     = S[mfi];
      const Real*      dx       = geom.CellSize();

      for (int d=0; d<BL_SPACEDIM; ++d) {
        Box ebox = surroundingNodes(box,d);
        ftmp[d].resize(ebox,nspecies+3);
        ftmp[d].copy((*flux[d])[mfi],ebox,0,ebox,0,nspecies+1);
      }

      enth_diff_terms(BL_TO_FORTRAN_BOX(box),
                      BL_TO_FORTRAN_BOX(domain),
                      dx,
                      BL_TO_FORTRAN_N_ANYD(T,TComp),
                      BL_TO_FORTRAN_N_ANYD(RhoY,RhoYComp),

                      BL_TO_FORTRAN_N_ANYD((*beta[0])[mfi],dComp),
                      BL_TO_FORTRAN_N_ANYD(ftmp[0],FComp),
                      BL_TO_FORTRAN_ANYD(area[0][mfi]),

                      BL_TO_FORTRAN_N_ANYD((*beta[1])[mfi],dComp),
                      BL_TO_FORTRAN_N_ANYD(ftmp[1],FComp),
                      BL_TO_FORTRAN_ANYD(area[1][mfi]),
#if AMREX_SPACEDIM == 3
                      BL_TO_FORTRAN_N_ANYD((*beta[2])[mfi],dComp),
                      BL_TO_FORTRAN_N_ANYD(ftmp[2],FComp),
                      BL_TO_FORTRAN_ANYD(area[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(fh),
                      Tbc.vect() );
      
      for (int d=0; d<BL_SPACEDIM; ++d) {
        Box etbox = mfi.nodaltilebox(d);
        (*flux[d])[mfi].copy(ftmp[d],etbox,nspecies+1,etbox,nspecies+1,2);
      }

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

  //
  // Do implicit c-n solve for velocity
  // compute the viscous forcing
  // do following except at initial iteration--rbp, per jbb
  //
  if (is_diffusive[Xvel])
  {
    const Real strt_time = ParallelDescriptor::second();

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
    FluxBoxes fb_betan, fb_betanp1;

    diffuse_velocity_setup(dt, delta_rhs, fb_betan, fb_betanp1);

    int rhsComp  = 0;
    int betaComp = 0;
    diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                delta_rhs,rhsComp,fb_betan.get(),fb_betanp1.get(),betaComp);

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
PeleLM::diffuse_velocity_setup (Real        dt,
                                MultiFab*&  delta_rhs,
                                FluxBoxes&  fb_betan,
                                FluxBoxes&  fb_betanp1)
{
  //
  // Do setup for implicit c-n solve for velocity.
  //
  BL_ASSERT(delta_rhs==0);
  const Real time = state[State_Type].prevTime();
  //
  // Assume always variable viscosity.
  //
  MultiFab** betan = fb_betan.define(this);
  MultiFab** betanp1 = fb_betanp1.define(this);
    
  getViscosity(betan, time);
  getViscosity(betanp1, time+dt);
  //
  // Do the following only if it is a non-trivial operation.
  //
  if (S_in_vel_diffusion)
  {
    //
    // Include div mu S*I terms in rhs
    //  (i.e. make nonzero delta_rhs to add into RHS):
    //
    // The scalar and tensor solvers incorporate the relevant pieces of
    //  of Div(tau), provided the flow is divergence-free.  Howeever, if
    //  Div(U) =/= 0, there is an additional piece not accounted for,
    //  which is of the form A.Div(U).  For constant viscosity, Div(tau)_i
    //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.  For mu not constant,
    //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
    //
    // As a convenience, we treat this additional term as a "source" in
    // the diffusive solve, computing Div(U) in the "normal" way we
    // always do--via a call to calc_divu.  This routine computes delta_rhs
    // if necessary, and stores it as an auxilliary rhs to the viscous solves.
    // This is a little strange, but probably not bad.
    //
    delta_rhs = new MultiFab(grids,dmap,BL_SPACEDIM,0);
    delta_rhs->setVal(0);
      
    MultiFab divmusi(grids,dmap,BL_SPACEDIM,0);
    //
    // Assume always variable viscosity.
    //
    diffusion->compute_divmusi(time,betan,divmusi);
    divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
    (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    //
    // Assume always variable viscosity.
    //
    diffusion->compute_divmusi(time+dt,betanp1,divmusi);
    divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
    (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
  }
}

void
PeleLM::getViscTerms (MultiFab& visc_terms,
                      int       src_comp, 
                      int       num_comp,
                      Real      time)
{
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
  //
  MultiFab** vel_visc  = 0;        // Potentially reused, raise scope
  FluxBoxes fb;
  const int  nGrow     = visc_terms.nGrow();
  //
  // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
  //
  visc_terms.setBndry(1.e30);
  if (src_comp < BL_SPACEDIM)
  {
    if (src_comp != Xvel || num_comp < BL_SPACEDIM)
      amrex::Error("tensor v -> getViscTerms needs all v-components at once");

    vel_visc = fb.define(this);
    getViscosity(vel_visc, time);

    showMF("velVT",*viscn_cc,"velVT_viscn_cc",level);
    for (int dir=0; dir<BL_SPACEDIM; ++dir) {
      showMF("velVT",*(vel_visc[dir]),amrex::Concatenate("velVT_viscn_",dir,1),level);
    }

    int viscComp = 0;
    diffusion->getTensorViscTerms(visc_terms,time,vel_visc,viscComp);
    showMF("velVT",visc_terms,"velVT_visc_terms_1",level);
  }
  else
  {
    amrex::Abort("Should only call getViscTerms for velocity");
  }
  //
  // Add Div(u) term if desired, if this is velocity, and if Div(u) is nonzero
  // If const-visc, term is mu.Div(u)/3, else it's -Div(mu.Div(u).I)*2/3
  //
  if (src_comp < BL_SPACEDIM && S_in_vel_diffusion)
  {
    if (num_comp < BL_SPACEDIM)
      amrex::Error("getViscTerms() need all v-components at once");
    
    MultiFab divmusi(grids,dmap,BL_SPACEDIM,0);
    showMF("velVT",get_old_data(Divu_Type),"velVT_divu",level);
    //
    // Assume always using variable viscosity.
    //
    diffusion->compute_divmusi(time,vel_visc,divmusi); // pre-computed visc above
    divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
    showMF("velVT",divmusi,"velVT_divmusi",level);
    visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
    showMF("velVT",visc_terms,"velVT_visc_terms_3",level);
  }
  //
  // Ensure consistent grow cells
  //
  if (nGrow > 0)
  {
    visc_terms.FillBoundary(0,num_comp, geom.periodicity());
    BL_ASSERT(visc_terms.nGrow() == 1);
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

Real MFnorm(const MultiFab& mf,
            int             sComp,
            int             nComp)
{
  const int p = 0; // max norm
  Real normtot = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    normtot = std::max(normtot,mf[mfi].norm(mfi.tilebox(),p,sComp,nComp));
  ParallelDescriptor::ReduceRealMax(normtot);
  return normtot;
}

void
PeleLM::compute_differential_diffusion_fluxes (const Real& time,
                                               const Real& dt)
{
  BL_PROFILE("HT:::compute_differential_diffusion_fluxes()");
  // explicit computation of species and enthalpy (heat) diffusion fluxes
  // save fluxes in class data
  const Real      strt_time = ParallelDescriptor::second();
  const TimeLevel whichTime = which_time(State_Type,time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    

  MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

  if (hack_nospecdiff)
  {
    amrex::Error("compute_differential_diffusion_fluxes: hack_nospecdiff not implemented");
  }

  const int  rho_flag  = 2;
  MultiFab   rhonph;      // Never need alpha for RhoY
  MultiFab*  alpha     = 0;
  int        alphaComp = 0;
  const Real a         = 0;
  const Real b         = 1;
  MultiFab&  S         = get_data(State_Type,time);

  // allocate edge-beta for species, RhoH, and Temp
  FluxBoxes fb(this, nspecies+2, 0);
  MultiFab** beta = fb.get();

  // average transport coefficients for species, RhoH, and Temp to edges
  getDiffusivity(beta, time, first_spec, 0, nspecies+1);
  getDiffusivity(beta, time, Temp, nspecies+1, 1);

  showMF("dd",S,"dd_preFP",level);
  //
  // Fill grow cells in the state for (Rho,RhoY,RhoH,T).  For Dirichlet physical
  // boundary grow cells, this data will live on the cell face, otherwise
  // it will live at cell centers.
  //
  const int nGrow = 1;

  BL_ASSERT(S.nGrow()>=nGrow);

  FillPatchIterator Yfpi(*this,S,nGrow,time,State_Type,Density,nspecies+2),
         Tfpi(*this,S,nGrow,time,State_Type,Temp,1);
         
  MultiFab& Tmf = Tfpi.get_mf();
  MultiFab& Ymf = Yfpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(Tmf,true); mfi.isValid();++mfi)
  {
    const Box& vbox   = mfi.tilebox();
    FArrayBox& fab    = S[mfi];
    const FArrayBox& Tfab = Tmf[mfi];
    const FArrayBox& Yfab = Ymf[mfi];

    BoxList    gcells = amrex::boxDiff(Box(vbox).grow(nGrow),vbox);
    for (BoxList::const_iterator it = gcells.begin(), end = gcells.end(); it != end; ++it)
    {
      const Box& gbox = *it;
      fab.copy(Yfab,gbox,0,gbox,Density,nspecies+2);
      fab.copy(Tfab,gbox,0,gbox,Temp,1);
    }
  }
}
  showMF("dd",S,"dd_rsT_fp",level);
  //
  // Create and fill (Rho,RhoY,RhoH,T) at coarser level
  //    
  MultiFab S_crse;
  MultiFab Phi_crse;
  if (level > 0) 
  {
    PeleLM* coarser = (PeleLM*) &(parent->getLevel(level-1));

    Phi_crse.define(coarser->grids,coarser->dmap,        1, 2);
    S_crse.define  (coarser->grids,coarser->dmap,S.nComp(), 0);

    FillPatchIterator S_fpi(*coarser,S_crse,0,time,State_Type,0,S.nComp());
    MultiFab& Smf=S_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif  
  for (MFIter mfi(Smf,true); mfi.isValid();++mfi)
    {
      const Box& vbx  = mfi.tilebox();
      const FArrayBox& sfab = Smf[mfi];

      S_crse[mfi].copy(sfab,vbx,0,vbx,0,S.nComp());
    }
  }

  MultiFab Phi(grids,dmap,1,nGrow);
  const int phiComp = 0;

  for (int sigma = 0; sigma < nspecies+1; ++sigma)
  {
    const int state_ind = first_spec + sigma;
    ViscBndry visc_bndry;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.growntilebox();
      Phi[mfi].copy(S[mfi],box,state_ind,box,0,1);
      Phi[mfi].divide(S[mfi],box,Density,0,1);
    }

    if (level > 0) 
    {
      MultiFab::Copy(Phi_crse,S_crse,state_ind,0,1,0);
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_crse,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.growntilebox();
        Phi_crse[mfi].divide(S_crse[mfi],box,Density,0,1);
      }
      Phi_crse.FillBoundary(getLevel(level-1).geom.periodicity());
    }

    diffusion->getBndryDataGivenS(visc_bndry,Phi,Phi_crse,state_ind,0,1);

    const bool     do_applyBC           = true;
    const bool     bndry_already_filled = true;            
    const int      betaComp             = sigma;
    Real*          rhsscale             = 0;
    ABecLaplacian* visc_op = diffusion->getViscOp(state_ind,a,b,time,visc_bndry,
                                                  rhonph,rho_flag,rhsscale,
                                                  beta,betaComp,alpha,alphaComp,bndry_already_filled);
    visc_op->maxOrder(diffusion->maxOrder());
    visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Phi,do_applyBC,LinOp::Inhomogeneous_BC,phiComp,sigma);

    delete visc_op;
  }

  Phi.clear();

  //
  // Remove scaling left in fluxes from solve.
  //
  for (int d=0; d < BL_SPACEDIM; ++d)
    flux[d]->mult(b/geom.CellSize()[d],0,nspecies+1);

#ifdef USE_WBAR
  compute_Wbar_fluxes(time,0);

  // add grad Wbar fluxes (SpecDiffusionFluxWbar) to 
  // species diffusion fluxes (flux)
  for (int d=0; d<BL_SPACEDIM; ++d)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*SpecDiffusionFluxWbar[d],true); mfi.isValid(); ++mfi)
    {
      const Box& ebox = mfi.tilebox();
      (*flux[d])[mfi].plus((*SpecDiffusionFluxWbar[d])[mfi],ebox,0,0,nspecies);
    }
  }
#endif

  //
  // Modify update/fluxes to preserve flux sum = 0 (conservatively correct Gamma_m)
  adjust_spec_diffusion_fluxes(time);
  //
  // AJN FLUXREG
  // We have just explicitly computed "D" for Y_m and h given an input state.
  // Update VISCOUS flux registers as follows.
  // If this is the predictor and sdc_iterMAX>1:
  //   COPY (1/2)*Gamma_m^n and (1/2)*lambda^n/cp^n grad h^n to flux registers
  // If this is the predictor AND sdc_iterMAX=1:
  //   DO NOTHING - we never call compute_differential_diffusion_terms again to
  //   subtract (1/2)*Gamma_m^{(k)} and (1/2)*lambda/cp grad h^{(k)}, which in this
  //   case is equal to the t^n version
  // If updateFluxReg=T (we are in the final corrector):
  //   SUBTRACT (1/2)*Gamma_m^(k) and (1/2)*lambda^(k)/cp^(k) grad h^(k)
  //
  if ( do_reflux && ( (is_predictor && sdc_iterMAX>1) || updateFluxReg ) )
  {
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
      if (level > 0)
      {
        if (is_predictor && sdc_iterMAX>1)
        {
          getViscFluxReg().FineAdd(*SpecDiffusionFluxn[d],d,0,
                                   first_spec,nspecies+1,0.5*dt);
        }
        if (updateFluxReg)
        {
          getViscFluxReg().FineAdd(*SpecDiffusionFluxnp1[d],d,0,
                                   first_spec,nspecies+1,-0.5*dt);
        }
      }

      if (level < parent->finestLevel())
      {
        if (is_predictor && sdc_iterMAX>1)
        {
          getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxn[d]),d,0,
                                                      first_spec,nspecies+1,-0.5*dt,
                                                      FluxRegister::COPY);
        }
        if (updateFluxReg)
        {
          getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,0,
                                                      first_spec,nspecies+1,0.5*dt,
                                                      FluxRegister::ADD);
        }
      }
    }
  }
  // build heat flux and temperature source terms
  // compute lambda grad T (for temperature and divu)
  // compute sum_m (Gamma_m + lambda/cp grad Y) (for enthalpy)
  // compute sum_m Gamma_m dot grad h_m (for divu)
  //
  compute_enthalpy_fluxes(time,beta);

  //
  // AJN FLUXREG
  // We have just computed "DD" given an input state.
  // Update ADVECTIVE flux registers as follows.
  // If we are in the predictor and sdc_iterMAX>1:
  //   ADD (1/2)*h_m^n(Gamma_m^n-lambda^n/cp^n grad Y_m^n)
  // If we are in the predictor and sdc_iterMAX=1:
  //   ADD 1.0*h_m^n(Gamma_m^n-lambda^n/cp^n grad Y_m^n)
  // If updateFluxReg=T (we are in the final corrector):
  //   ADD (1/2)*h_m^(k)(Gamma_m^(k)-lambda^(k)/cp^(k) grad Y_m^(k))
  //
  if ( do_reflux && ( is_predictor || updateFluxReg ) )
  {
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
      if (level > 0)
      {
        if (is_predictor)
        {
          const Real fac = (sdc_iterMAX==1) ? dt : 0.5*dt;
          advflux_reg->FineAdd(*SpecDiffusionFluxn[d],d,
                               nspecies+1,RhoH,1,fac);
        }
        if (updateFluxReg)
        {
          advflux_reg->FineAdd(*SpecDiffusionFluxnp1[d],d,
                               nspecies+1,RhoH,1,0.5*dt);
        }
      }

      if (level < parent->finestLevel())
      {
        if (is_predictor)
        {
          const Real fac = (sdc_iterMAX==1) ? dt : 0.5*dt;
          getAdvFluxReg(level+1).CrseInit((*SpecDiffusionFluxn[d]),d,nspecies+1,RhoH,
                                          1,-fac,FluxRegister::ADD);
        }
        if (updateFluxReg)
        {
          getAdvFluxReg(level+1).CrseInit((*SpecDiffusionFluxnp1[d]),d,nspecies+1,RhoH,
                                          1,-0.5*dt,FluxRegister::ADD);
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
#pragma omp parallel
#endif  
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
  {
    const Box& box   = mfi.tilebox();
    const int nc = last_scalar - first_scalar + 1; 
    FArrayBox& snew = S_new[mfi];
      
    snew.copy( (*aofs)[mfi],box,first_scalar,box,first_scalar,nc);
    snew.mult(dt,box,first_scalar,nc);
    snew.plus(S_old[mfi],box,first_scalar,first_scalar,nc);            
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
  BL_ASSERT(fdiv.nComp() >= fdivComp+nComp);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(fdiv,true); mfi.isValid(); ++mfi)
  {
    const Box& box = mfi.tilebox();
    FArrayBox& fab = fdiv[mfi];

    flux_div( BL_TO_FORTRAN_BOX(box),
              BL_TO_FORTRAN_N_ANYD(fab,fdivComp),
              BL_TO_FORTRAN_N_ANYD((*f[0])[mfi],fluxComp),
              BL_TO_FORTRAN_N_ANYD((*f[1])[mfi],fluxComp),
#if ( AMREX_SPACEDIM == 3 )
              BL_TO_FORTRAN_N_ANYD((*f[2])[mfi],fluxComp),
#endif
              BL_TO_FORTRAN_ANYD(volume[mfi]),
              &nComp, &scale);
  }
}

void
PeleLM::compute_differential_diffusion_terms (MultiFab& D,
                                              MultiFab& DD,
#ifdef USE_WBAR
                                              MultiFab& DWbar,
#endif
                                              Real      time,
                                              Real      dt)
{
  BL_PROFILE("HT:::compute_differential_diffusion_terms()");
  // 
  // Sets vt for species, RhoH and Temp together
  // Uses state at time to explicitly compute fluxes, and resets internal
  //  data for fluxes, etc
  //
  BL_ASSERT(D.boxArray() == grids);
  BL_ASSERT(D.nComp() >= nspecies+2); // room for spec+RhoH+Temp

  if (hack_nospecdiff)
  {
    amrex::Error("compute_differential_diffusion_terms: hack_nospecdiff not implemented");
  }

  const TimeLevel whichTime = which_time(State_Type,time);
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
  MultiFab& sumSpecFluxDotGradH = (whichTime == AmrOldTime) ? sumSpecFluxDotGradHn : sumSpecFluxDotGradHnp1;
  MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;
#ifdef USE_WBAR
  MultiFab* const * fluxWbar = SpecDiffusionFluxWbar;
#endif
  //
  // Compute/adjust species fluxes/heat flux/conduction, save in class data
  compute_differential_diffusion_fluxes(time,dt);


  // Then compute:
  // -Div(Fi),
  // -Div((lambda/cp).Grad(h) + Fi.(Lei-1)) + heating,
  // -Div(lambda.Grad(T)) + heating + sum(Fi.Grad(Hi)) 
  //

  // compute div Gamma_m for species AND 
  // div lambda/cp grad h for enthalpy
  flux_divergence(D,0,flux,0,nspecies+1,-1);

  // compute div sum_m h_m (Gamma_m + lambda/cp grad Y_m), a.k.a. the "diffdiff" terms
  flux_divergence(DD,0,flux,nspecies+1,1,-1);

  // compute div lambda grad T for temperature
  flux_divergence(D,nspecies+1,flux,nspecies+2,1,-1);

#ifdef USE_WBAR
  // compute div beta_for_Wbar grad Wbar
  flux_divergence(DWbar,0,fluxWbar,0,nspecies,-1);
#endif

  // add sum_m Gamma_m dot grad h_m to D for temperature
  MultiFab::Add(D,sumSpecFluxDotGradH,0,nspecies+1,1,0);

  if (D.nGrow() > 0 && DD.nGrow() > 0)
  {
    const int nc = nspecies+2;

    D.FillBoundary(0,nc, geom.periodicity());
    DD.FillBoundary(0,1, geom.periodicity());

    BL_ASSERT(D.nGrow() == 1);
    BL_ASSERT(DD.nGrow() == 1);
	
    Extrapolater::FirstOrderExtrap(D, geom, 0, nc);
    Extrapolater::FirstOrderExtrap(DD, geom, 0, 1);
  }
}

void
PeleLM::set_rho_to_species_sum (MultiFab& S,
                                int       strtcomp, 
                                int       nghost_in,
                                int       minzero)
{
  set_rho_to_species_sum(S, strtcomp, S, strtcomp, nghost_in, minzero);
}

//
// This function
//       sets the density component in S_out to the sum of 
//       the species components in S_in
// if minzero = 1, the species components are first "max-ed" w/ zero
//  thus changing S_in values    
//
// s_in_strt is the state component corresponding to the
// 0-th component of S_in. It is otherwise assumed that 
// that the components of S_in "align" with those in S_out.
//
void
PeleLM::set_rho_to_species_sum (MultiFab& S_in, 
                                int       s_in_strt,
                                MultiFab& S_out,
                                int       s_out_strt,  
                                int       nghost_in,
                                int       minzero)

{
  const BoxArray& sgrids = S_in.boxArray();

  BL_ASSERT(sgrids == S_out.boxArray());

  const int s_first_spec = first_spec - s_in_strt;
  const int s_last_spec  = last_spec  - s_in_strt;
  const int s_num_spec   = last_spec - first_spec + 1;
  const int s_density    = Density-s_out_strt;
  const int nghost       = std::min(S_in.nGrow(), std::min(S_out.nGrow(), nghost_in));

  if (minzero)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_in,true); mfi.isValid(); ++mfi)
    {
      const Box&  box = mfi.tilebox();
      FabMinMax(S_in[mfi], box, 0.0, Real_MAX, s_first_spec, s_num_spec);
    }
  }

  BL_ASSERT(s_density >= 0);

  S_out.setVal(0, s_density, 1, nghost);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_in,true); mfi.isValid(); ++mfi)
  {
    const Box&       box = mfi.tilebox();
    for (int spec = s_first_spec; spec <= s_last_spec; spec++)
    {
      S_out[mfi].plus(S_in[mfi],box,spec,s_density,1);
    }
  }

  make_rho_curr_time();
}

void
PeleLM::temperature_stats (MultiFab& S)
{
  if (verbose)
  {
    //
    // Calculate some minimums and maximums.
    //
    auto scaleMin = VectorMin({&S},FabArrayBase::mfiter_tile_size,Density,NUM_STATE-BL_SPACEDIM,0);
    auto scaleMax = VectorMax({&S},FabArrayBase::mfiter_tile_size,Density,NUM_STATE-BL_SPACEDIM,0);

    bool aNegY = false;
    for (int i = 0; i < nspecies && !aNegY; ++i) {
      if (scaleMin[first_spec+i-BL_SPACEDIM] < 0) aNegY = true;
    }

    amrex::Print() << "  Min,max temp = " << scaleMin[Temp   - BL_SPACEDIM]
                   << ", "                << scaleMax[Temp   - BL_SPACEDIM] << '\n';
    amrex::Print() << "  Min,max rho  = " << scaleMin[Density- BL_SPACEDIM]
                   << ", "                << scaleMax[Density- BL_SPACEDIM] << '\n';
    amrex::Print() << "  Min,max rhoh = " << scaleMin[RhoH   - BL_SPACEDIM]
                   << ", "                << scaleMax[RhoH   - BL_SPACEDIM] << '\n';

    if (aNegY){
        Vector<std::string> names;
        PeleLM::getSpeciesNames(names);
        amrex::Print() << "  Species w/min < 0: ";
        for (int i = 0; i < nspecies; ++i) {
        int idx = first_spec + i - BL_SPACEDIM;
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
                       MultiFab&       p, 
                       int             pComp,
                       const MultiFab* temp)
{
  BL_ASSERT(pComp<p.nComp());

  const Real strt_time = ParallelDescriptor::second();

  int nCompY = last_spec - first_spec + 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox sfab, tmp;
  
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
      const int  i      = mfi.index();
      //const Box& box    = sgrids[i];
      const Box& box = mfi.tilebox();
      const int  sCompR = 0;
      const int  sCompT = 1;
      const int  sCompY = 2;
    
      sfab.resize(box,nCompY+2);
      BL_ASSERT(S[mfi].box().contains(box));
      sfab.copy(S[mfi],box,Density,box,sCompR,1);
  
      if (temp)
      {
        BL_ASSERT(temp->boxArray()[i].contains(box));
        sfab.copy((*temp)[mfi],box,0,box,sCompT,1);
      }
      else
      {
        sfab.copy(S[mfi],box,Temp,box,sCompT,1);
      }
      sfab.copy(S[mfi],box,first_spec,box,sCompY,nCompY);
  
      tmp.resize(box,1);
      tmp.copy(sfab,box,sCompR,box,0,1);
      tmp.invert(1);
  
      for (int k = 0; k < nCompY; k++)
        sfab.mult(tmp,box,0,sCompY+k,1);
  
      getPGivenRTY_pphys(p[mfi],sfab,sfab,sfab,box,sCompR,sCompT,sCompY,pComp);
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
    for (int k = 0; k <= finest_level; k++)
    {
      AmrLevel&       ht = getLevel(k);
      const MultiFab& S  = ht.get_new_data(State_Type);
      const BoxArray& ba = ht.boxArray();
      const DistributionMapping& dm = ht.DistributionMap();
      MultiFab hmix(ba,dm,1,0);
      MultiFab::Copy(hmix,S,RhoH,0,1,0);
      MultiFab::Divide(hmix,S,Density,0,1,0);
      if (k != finest_level) 
      {
        AmrLevel& htf = getLevel(k+1);
        BoxArray  baf = htf.boxArray();
        baf.coarsen(parent->refRatio(k));
#ifdef _OPENMP
#pragma omp parallel
#endif  
        for (MFIter mfi(hmix,true); mfi.isValid(); ++mfi)
        {
          baf.intersections(ba[mfi.index()],isects);
          for (int i = 0; i < isects.size(); i++)
            hmix[mfi].setVal(0,isects[i].second,0,1);
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

#ifdef USE_WBAR
  calcDiffusivity_Wbar(time);
#endif

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
    
  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
  MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
    
  const int pComp = (have_rhort ? RhoRT : Trac);

  compute_rhoRT (S,S,pComp);
}

Real
PeleLM::predict_velocity (Real  dt)
{
  if (verbose) Print() << "... predict edge velocities\n";
  //
  // Get simulation parameters.
  //
  const int   nComp          = BL_SPACEDIM;
  const Real* dx             = geom.CellSize();
  const Real  prev_time      = state[State_Type].prevTime();
  const Real  prev_pres_time = state[Press_Type].prevTime();
  const Real  strt_time      = ParallelDescriptor::second();
  //
  // Compute viscous terms at level n.
  // Ensure reasonable values in 1 grow cell.  Here, do extrap for
  // c-f/phys boundary, since we have no interpolator fn, also,
  // preserve extrap for corners at periodic/non-periodic intersections.
  //
  MultiFab visc_terms(grids,dmap,nComp,1);

  if (be_cn_theta != 1.0)
  {
    getViscTerms(visc_terms,Xvel,nComp,prev_time);
  }
  else
  {
    visc_terms.setVal(0);
  }

  MultiFab Gp(grids,dmap,BL_SPACEDIM,1);
  getGradP(Gp, prev_pres_time);
    
  FillPatchIterator U_fpi(*this,visc_terms,Godunov::hypgrow(),prev_time,State_Type,Xvel,BL_SPACEDIM);
  MultiFab& Umf=U_fpi.get_mf();

  // Floor small values of states to be extrapolated
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Umf,true); mfi.isValid(); ++mfi)
  {
    Box gbx=mfi.growntilebox(Godunov::hypgrow());
    auto fab = Umf.array(mfi);
    AMREX_HOST_DEVICE_FOR_4D ( gbx, BL_SPACEDIM, i, j, k, n,
    {
      auto& val = fab(i,j,k,n);
      val = std::abs(val) > 1.e-20 ? val : 0;
    });
  }

  FillPatchIterator S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);
  MultiFab& Smf=S_fpi.get_mf();

  //
  // Compute "grid cfl number" based on cell-centered time-n velocities
  //
  auto umax = VectorMaxAbs({&Umf},FabArrayBase::mfiter_tile_size,0,BL_SPACEDIM,Umf.nGrow());
  Real cflmax = dt*umax[0]/dx[0];
  for (int d=1; d<BL_SPACEDIM; ++d) {
    cflmax = std::max(cflmax,dt*umax[d]/dx[d]);
  }
  Real tempdt = std::min(change_max,cfl/cflmax);

#ifdef _OPENMP
#pragma omp parallel
#endif      
  {
    FArrayBox tforces, Uface[BL_SPACEDIM];
    Vector<int> bndry[BL_SPACEDIM];

    for (MFIter U_mfi(Umf,true); U_mfi.isValid(); ++U_mfi)
    {
      Box bx=U_mfi.tilebox();
      FArrayBox& Ufab = Umf[U_mfi];
      for (int d=0; d<BL_SPACEDIM; ++d) {
        Uface[d].resize(surroundingNodes(bx,d),1);
      }

      if (getForceVerbose)
        amrex::Print() << "---\nA - Predict velocity:\n Calling getForce..." << '\n';
      getForce(tforces,bx,1,Xvel,BL_SPACEDIM,prev_time,Ufab,Smf[U_mfi],0);

      //
      // Compute the total forcing.
      //
      godunov->Sum_tf_gp_visc(tforces,0,visc_terms[U_mfi],0,Gp[U_mfi],0,rho_ptime[U_mfi],0);

      D_TERM(bndry[0] = fetchBCArray(State_Type,bx,0,1);,
             bndry[1] = fetchBCArray(State_Type,bx,1,1);,
             bndry[2] = fetchBCArray(State_Type,bx,2,1););

      godunov->ExtrapVelToFaces(bx, dx, dt,
                                D_DECL(Uface[0], Uface[1], Uface[2]),
                                D_DECL(bndry[0], bndry[1], bndry[2]),
                                Ufab, tforces);


      for (int d=0; d<BL_SPACEDIM; ++d) {
        const Box& ebx = U_mfi.nodaltilebox(d);
        u_mac[d][U_mfi].copy(Uface[d],ebx,0,ebx,0,1);
      }
    }
  }

  if (verbose > 1)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    Print() << "PeleLM::predict_velocity(): lev: " << level 
            << ", time: " << run_time << '\n';
  }

  return dt*tempdt;
}

void
PeleLM::set_reasonable_grow_cells_for_R (Real time)
{
  //
  // Ensure reasonable grow cells for R.
  //
  MultiFab& React = get_data(RhoYdot_Type, time);
  React.FillBoundary(0,nspecies, geom.periodicity());
  BL_ASSERT(React.nGrow() == 1);
  Extrapolater::FirstOrderExtrap(React, geom, 0, nspecies);
}

Real
PeleLM::advance (Real time,
                 Real dt,
                 int  iteration,
                 int  ncycle)
{

  BL_PROFILE_VAR("HT::advance::mac", HTMAC);
  if (closed_chamber == 1 && level == 0)
  {
    // set new-time ambient pressure to be a copy of old-time ambient pressure
    p_amb_new = p_amb_old;
  }

  BL_PROFILE_VAR_STOP(HTMAC);

  BL_PROFILE_REGION_START("R::HT::advance()[src_sdc]");
  BL_PROFILE("HT::advance()[src_sdc]");
  is_predictor = true;
  updateFluxReg = false;

  if (level == 0)
  {
    crse_dt = dt;
    int thisLevelStep = parent->levelSteps(0);
    set_common(&time,&thisLevelStep);
  }

  if (verbose)
  {
    amrex::Print() << "PeleLM::advance(): at start of time step\n"
		   << "SDC Advancing level " << level
		   << " : starting time = " << time
		   << " with dt = "         << dt << '\n';
  }

  // swaps old and new states for all state types
  // then copies each of the old state types into the new state types
  // computes old transport coefficients and copies them into new
  BL_PROFILE_VAR("HT::advance::setup", HTSETUP);
  advance_setup(time,dt,iteration,ncycle);
  BL_PROFILE_VAR_STOP(HTSETUP);

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
    
  MultiFab& S_new = get_new_data(State_Type);
  MultiFab& S_old = get_old_data(State_Type);

  const Real prev_time = state[State_Type].prevTime();
  const Real tnp1  = state[State_Type].curTime();


  Real dt_test = 0.0;    

  BL_PROFILE_VAR("HT::advance::diffusion", HTDIFF);
  if (floor_species == 1)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();            
      FArrayBox& species = S_old[mfi];
      floor_spec(BL_TO_FORTRAN_BOX(box), 
                 BL_TO_FORTRAN_N_ANYD(species,first_spec));
    }
  }
  BL_PROFILE_VAR_STOP(HTDIFF);

  // Build a copy of old-time rho with grow cells for use in the diffusion solves
  BL_PROFILE_VAR_START(HTDIFF);
  make_rho_prev_time();
  BL_PROFILE_VAR_STOP(HTDIFF);

  // compute old-time thermodynamic pressure
  BL_PROFILE_VAR_START(HTMAC);
  setThermoPress(prev_time);  
  BL_PROFILE_VAR_STOP(HTMAC);

  MultiFab Dn(grids,dmap,nspecies+2,nGrowAdvForcing);
  MultiFab DDn(grids,dmap,1,nGrowAdvForcing);
#ifdef USE_WBAR
  MultiFab DWbar(grids,dmap,nspecies,nGrowAdvForcing);
#endif

  // Compute Dn and DDn (based on state at tn)
  //  (Note that coeffs at tn and tnp1 were intialized in _setup)
  if (verbose) amrex::Print() << "Computing Dn, DDn, and DWbar \n";

  BL_PROFILE_VAR_START(HTDIFF);
#ifdef USE_WBAR
  compute_differential_diffusion_terms(Dn,DDn,DWbar,prev_time,dt);
#else
  compute_differential_diffusion_terms(Dn,DDn,prev_time,dt);
#endif
  BL_PROFILE_VAR_STOP(HTDIFF);

  showMF("sdc",Dn,"sdc_Dn",level,parent->levelSteps(level));
  showMF("sdc",DDn,"sdc_DDn",level,parent->levelSteps(level));
#ifdef USE_WBAR
  showMF("sdc",DWbar,"sdc_DWbar",level,parent->levelSteps(level));
#endif

  /*
    You could compute instantaneous I_R here but for now it's using either the
    previous step or divu_iter's version of I_R.  Either way, we have to make 
    sure that the nGrowAdvForcing grow cells have something reasonable in them
  */
  BL_PROFILE_VAR("HT::advance::reactions", HTREAC);
  set_reasonable_grow_cells_for_R(tnp1);
  BL_PROFILE_VAR_STOP(HTREAC);

  // copy old state into new state for Dn and DDn.
  // Note: this was already done for scalars, transport coefficients,
  // and divu in advance_setup
  MultiFab Dnp1(grids,dmap,nspecies+2,nGrowAdvForcing);
  MultiFab DDnp1(grids,dmap,1,nGrowAdvForcing);
  MultiFab chi(grids,dmap,1,nGrowAdvForcing);
  MultiFab chi_increment(grids,dmap,1,nGrowAdvForcing);
  MultiFab mac_divu(grids,dmap,1,nGrowAdvForcing);

  // used for closed chamber algorithm
  MultiFab theta_old(grids,dmap,1,nGrowAdvForcing);
  MultiFab theta_nph(grids,dmap,1,nGrowAdvForcing);
  Real Sbar;
  Real Sbar_old, Sbar_new;

  BL_PROFILE_VAR_START(HTDIFF);
  MultiFab::Copy(Dnp1,Dn,0,0,nspecies+2,nGrowAdvForcing);
  MultiFab::Copy(DDnp1,DDn,0,0,1,nGrowAdvForcing);
  BL_PROFILE_VAR_STOP(HTDIFF);

  BL_PROFILE_VAR_START(HTMAC);
  chi.setVal(0,nGrowAdvForcing);
  BL_PROFILE_VAR_STOP(HTMAC);

  is_predictor = false;

  BL_PROFILE_VAR_NS("HT::advance::velocity_adv", HTVEL);
  for (int sdc_iter=1; sdc_iter<=sdc_iterMAX; ++sdc_iter)
  {

    if (sdc_iter == sdc_iterMAX)
    {
      updateFluxReg = true;
    }

    if (sdc_iter > 1)
    {
      // compute new-time transport coefficients
      BL_PROFILE_VAR_START(HTDIFF);
      calcDiffusivity(tnp1);
#ifdef USE_WBAR
      calcDiffusivity_Wbar(tnp1);
#endif
      BL_PROFILE_VAR_STOP(HTDIFF);

      // compute Dnp1 and DDnp1
      // iteratively lagged
      if (verbose) 
	amrex::Print() << "Computing Dnp1 and DDnp1 (SDC iteration " << sdc_iter << ")\n";
      BL_PROFILE_VAR_START(HTDIFF);
#ifdef USE_WBAR
      compute_differential_diffusion_terms(Dnp1,DDnp1,DWbar,tnp1,dt);
#else
      compute_differential_diffusion_terms(Dnp1,DDnp1,tnp1,dt);
#endif
      BL_PROFILE_VAR_STOP(HTDIFF);

      // compute new-time DivU with instantaneous reaction rates
      BL_PROFILE_VAR_START(HTMAC);
      calc_divu(tnp1, dt, get_new_data(Divu_Type));
      BL_PROFILE_VAR_STOP(HTMAC);
    }

    // compute U^{ADV,*}
    BL_PROFILE_VAR_START(HTVEL);
    dt_test = predict_velocity(dt);  
    BL_PROFILE_VAR_STOP(HTVEL);

    // create S^{n+1/2} by averaging old and new
    BL_PROFILE_VAR_START(HTMAC);
    MultiFab Forcing(grids,dmap,nspecies+1,nGrowAdvForcing);
    Forcing.setBndry(1.e30);
    FillPatch(*this,mac_divu,nGrowAdvForcing,time+0.5*dt,Divu_Type,0,1,0);
    BL_PROFILE_VAR_STOP(HTMAC);
    showMF("sdc",Forcing,"sdc_mac_rhs1",level,sdc_iter,parent->levelSteps(level));
      
    // compute new-time thermodynamic pressure
    BL_PROFILE_VAR_START(HTMAC);
    setThermoPress(tnp1);

    // compute chi_increment
    chi_increment.setVal(0.0,nGrowAdvForcing);
    calc_dpdt(tnp1,dt,chi_increment,u_mac);
    BL_PROFILE_VAR_STOP(HTMAC);

    showMF("sdc",chi_increment,"sdc_chi_increment",level,sdc_iter,parent->levelSteps(level));

    // add chi_increment to chi
    BL_PROFILE_VAR_START(HTMAC);
    MultiFab::Add(chi,chi_increment,0,0,1,nGrowAdvForcing);
    BL_PROFILE_VAR_STOP(HTMAC);
    showMF("sdc",chi,"sdc_chi",level,sdc_iter,parent->levelSteps(level));

    // add chi to time-centered mac_divu
    BL_PROFILE_VAR_START(HTMAC);
    MultiFab::Add(mac_divu,chi,0,0,1,nGrowAdvForcing);
    BL_PROFILE_VAR_STOP(HTMAC);

    if (closed_chamber == 1 && level == 0)
    {	

      BL_PROFILE_VAR_START(HTMAC);

      // compute old, new, and time-centered theta = 1 / (gamma P)
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();            
        FArrayBox& thetafab = theta_old[mfi];
        const FArrayBox& rhoY = S_old[mfi];
        const FArrayBox& T = S_old[mfi];
        calc_gamma_pinv(BL_TO_FORTRAN_BOX(box),
                        BL_TO_FORTRAN_ANYD(thetafab),
                        BL_TO_FORTRAN_N_ANYD(rhoY,first_spec),
                        BL_TO_FORTRAN_N_ANYD(T,Temp),
                        &p_amb_old);
      }
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();            
        FArrayBox& thetafab = theta_nph[mfi];
        const FArrayBox& rhoY = S_new[mfi];
        const FArrayBox& T = S_new[mfi];
        calc_gamma_pinv(BL_TO_FORTRAN_BOX(box),
                        BL_TO_FORTRAN_ANYD(thetafab),
                        BL_TO_FORTRAN_N_ANYD(rhoY,first_spec),
                        BL_TO_FORTRAN_N_ANYD(T,Temp),
                        &p_amb_new);
      }
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(theta_nph,true); mfi.isValid(); ++mfi)
      {
        FArrayBox& th_nph = theta_nph[mfi];
        const FArrayBox& th_old = theta_old[mfi];
        const Box& box = mfi.tilebox();
        th_nph.plus(th_old,box,box,0,0,1);
        th_nph.mult(0.5,box,0,1);
      }

      // compute number of cells
      Real num_cells = grids.numPts();

      // compute the average of mac_divu theta
      Sbar = mac_divu.sum() / num_cells;
      thetabar = theta_nph.sum() / num_cells;

      // subtract mean from mac_divu and theta_nph
      mac_divu.plus(-Sbar,0,1);
      theta_nph.plus(-thetabar,0,1);

      p_amb_new = p_amb_old + dt*(Sbar/thetabar);
      dp0dt = Sbar/thetabar;

      // update mac rhs by adding delta_theta * (Sbar / thetabar)
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mac_divu,true); mfi.isValid(); ++mfi)
      {
        FArrayBox& m_du = mac_divu[mfi];
        FArrayBox& th_nph = theta_nph[mfi];
        const Box& box = mfi.tilebox();
        th_nph.mult(Sbar/thetabar,box,0,1);
        m_du.minus(th_nph,box,box,0,0,1);
      }
      BL_PROFILE_VAR_STOP(HTMAC);

      amrex::Print() << "level 0: prev_time, p_amb_old, p_amb_new, delta = " 
		     << prev_time << " " << p_amb_old << " " << p_amb_new << " "
                     << p_amb_new-p_amb_old << std::endl;
    }

    // MAC-project... and overwrite U^{ADV,*}
    showMF("sdc",Forcing,"sdc_Forcing_for_mac",level,sdc_iter,parent->levelSteps(level));
    BL_PROFILE_VAR_START(HTMAC);
    mac_project(time,dt,S_old,&mac_divu,nGrowAdvForcing,updateFluxReg);

    if (closed_chamber == 1 && level == 0)
    {
      // add Sbar back to mac_divu
      mac_divu.plus(Sbar,0,1);
    }
    BL_PROFILE_VAR_STOP(HTMAC);

    //
    // Compute A (advection terms) with F = Dn + R
    //
    showMF("sdc",get_new_data(RhoYdot_Type),"sdc_R_for_A_forcing",level,sdc_iter,parent->levelSteps(level));
    BL_PROFILE_VAR("HT::advance::advection", HTADV);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Forcing,true); mfi.isValid(); ++mfi)
    {
      FArrayBox& f = Forcing[mfi];
      const FArrayBox& dn = Dn[mfi];
      const FArrayBox& ddn = DDn[mfi];
      const FArrayBox& r = get_new_data(RhoYdot_Type)[mfi];
      const Box& gbox = mfi.growntilebox();
      f.copy(dn,gbox,0,gbox,0,nspecies+1); // add Dn to RhoY and RhoH
      f.plus(ddn,gbox,gbox,0,nspecies,1); // add DDn to RhoH forcing
      if (closed_chamber == 1)
      {
        f.plus(dp0dt,gbox,nspecies,1); // add dp0/dt to enthalpy forcing
      }
      f.plus(r,gbox,gbox,0,0,nspecies); // add R to RhoY, no contribution for RhoH
    }
    BL_PROFILE_VAR_STOP(HTADV);

    showMF("sdc",Forcing,"sdc_Forcing_for_A",level,sdc_iter,parent->levelSteps(level));
    BL_PROFILE_VAR_START(HTADV);
    Forcing.FillBoundary(0,nspecies+1,geom.periodicity());
    BL_PROFILE_VAR_STOP(HTADV);

    if (verbose) amrex::Print() << "A (SDC iter " << sdc_iter << ")\n";

    // compute A
    BL_PROFILE_VAR_START(HTADV);
    aofs->setVal(1.e30,aofs->nGrow());
    compute_scalar_advection_fluxes_and_divergence(Forcing,mac_divu,dt);
    BL_PROFILE_VAR_STOP(HTADV);
    showMF("sdc",*aofs,"sdc_A_pred",level,sdc_iter,parent->levelSteps(level));

    // update rho, rho*Y, and rho*h
    BL_PROFILE_VAR_START(HTADV);
    scalar_advection_update(dt, Density, RhoH);

    // update pointer to new-time density
    make_rho_curr_time();
    BL_PROFILE_VAR_STOP(HTADV);

    // 
    // Compute Dhat, diffuse with F 
    //                 = A + R + 0.5(Dn + Dnp1) - Dnp1 + Dhat + 0.5(DDn + DDnp1)
    //                 = A + R + 0.5(Dn - Dnp1) + Dhat + 0.5(DDn + DDnp1)
    // NOTE: Here we use 0.5*DDnp1 from the previous iteration
    // 
    if (verbose) amrex::Print() << "Dhat (SDC corrector " << sdc_iter << ")\n";

    showMF("sdc",*aofs,"sdc_A_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",get_new_data(RhoYdot_Type),"sdc_R_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",Dn,"sdc_Dn_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",DDn,"sdc_DDn_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",Dnp1,"sdc_Dnp1_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",DDnp1,"sdc_DDnp1_before_Dhat",level,sdc_iter,parent->levelSteps(level));

    BL_PROFILE_VAR_START(HTDIFF);
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    for (MFIter mfi(Forcing,true); mfi.isValid(); ++mfi) 
    {
      const Box& box = mfi.tilebox();
      FArrayBox& f = Forcing[mfi];
      const FArrayBox& a = (*aofs)[mfi];
      const FArrayBox& r = get_new_data(RhoYdot_Type)[mfi];
      const FArrayBox& dn = Dn[mfi];
      const FArrayBox& ddn = DDn[mfi];
      const FArrayBox& dnp1 = Dnp1[mfi];
      const FArrayBox& ddnp1 = DDnp1[mfi];

      f.copy(dn,box,0,box,0,nspecies+1); // copy Dn into RhoY and RhoH
      f.minus(dnp1,box,box,0,0,nspecies+1); // subtract Dnp1 from RhoY and RhoH
      f.plus(ddn  ,box,box,0,nspecies,1); // add DDn to RhoH, no contribution for RhoY
      f.plus(ddnp1,box,box,0,nspecies,1); // add DDnp1 to RhoH, no contribution for RhoY
      f.mult(0.5,box,0,nspecies+1);
#ifdef USE_WBAR
      const FArrayBox& dwbar = DWbar[mfi];
      f.plus(dwbar,box,box,0,0,nspecies); // add DWbar to RhoY
#endif
      if (closed_chamber == 1)
      {
        f.plus(dp0dt,box,nspecies,1); // add dp0/dt to enthalpy forcing
      }
      f.plus(a,box,box,first_spec,0,nspecies+1); // add A into RhoY and RhoH
      f.plus(r,box,box,0,0,nspecies); // no reactions for RhoH
    }
  }
    BL_PROFILE_VAR_STOP(HTDIFF);

    MultiFab Dhat(grids,dmap,nspecies+2,nGrowAdvForcing);

    // advection-diffusion solve
    showMF("sdc",Forcing,"sdc_Forcing_before_Dhat",level,sdc_iter,parent->levelSteps(level));
    BL_PROFILE_VAR_START(HTDIFF);
    differential_diffusion_update(Forcing,0,Dhat,0,DDnp1);
    BL_PROFILE_VAR_STOP(HTDIFF);

    showMF("sdc",Dn,"sdc_Dn_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",Dnp1,"sdc_Dnp1_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",DDn,"sdc_DDn_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",DDnp1,"sdc_DDnp1_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",Dhat,"sdc_Dhat_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",*aofs,"sdc_A_before_R",level,sdc_iter,parent->levelSteps(level));

    // 
    // Compute R (F = A + 0.5(Dn - Dnp1 + DDn + DDnp1) + Dhat )
    // 
    BL_PROFILE_VAR_START(HTREAC);
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    for (MFIter mfi(Forcing,true); mfi.isValid(); ++mfi) 
    {
      const Box& box = mfi.tilebox();
      FArrayBox& f = Forcing[mfi];
      const FArrayBox& a = (*aofs)[mfi];
      const FArrayBox& dn = Dn[mfi];
      const FArrayBox& dnp1 = Dnp1[mfi];
      const FArrayBox& dhat = Dhat[mfi];
      const FArrayBox& ddn = DDn[mfi];
      const FArrayBox& ddnp1 = DDnp1[mfi];
      f.copy(dn,box,0,box,0,nspecies+1);    // copy Dn into RhoY and RhoH
      f.minus(dnp1,box,box,0,0,nspecies+1); // subtract Dnp1 from RhoY and RhoH
      f.plus(ddn  ,box,box,0,nspecies,1);   // add DDn to RhoH, no contribution for RhoY
      f.plus(ddnp1,box,box,0,nspecies,1);   // add DDnp1 to RhoH, no contribution for RhoY
      f.mult(0.5,box,0,nspecies+1);
      if (closed_chamber == 1)
      {
        f.plus(dp0dt,box,nspecies,1); // add dp0/dt to enthalpy forcing
      }
      f.plus(dhat,box,box,0,0,nspecies+1);       // add Dhat to RhoY and RHoH
      f.plus(a,box,box,first_spec,0,nspecies+1); // add A to RhoY and RhoH
    }
  }
  
    BL_PROFILE_VAR_STOP(HTREAC);
    Dhat.clear();

    if (verbose) amrex::Print() << "R (SDC corrector " << sdc_iter << ")\n";

    showMF("sdc",S_old,"sdc_Sold_before_R",level,sdc_iter,parent->levelSteps(level));
    showMF("sdc",Forcing,"sdc_Forcing_before_R",level,sdc_iter,parent->levelSteps(level));
    BL_PROFILE_VAR_START(HTREAC);
    advance_chemistry(S_old,S_new,dt,Forcing,0);
    BL_PROFILE_VAR_STOP(HTREAC);
    showMF("sdc",S_new,"sdc_Snew_after_R",level,sdc_iter,parent->levelSteps(level));

    BL_PROFILE_VAR_START(HTDIFF);
    if (floor_species == 1)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();            
        FArrayBox& species = S_new[mfi];
        floor_spec(BL_TO_FORTRAN_BOX(box), 
                   BL_TO_FORTRAN_N_ANYD(species,first_spec));
      }
    }
    BL_PROFILE_VAR_STOP(HTDIFF);

    if (verbose) amrex::Print() << "DONE WITH R (SDC corrector " << sdc_iter << ")\n";

    BL_PROFILE_VAR_START(HTMAC);
    setThermoPress(tnp1);
    BL_PROFILE_VAR_STOP(HTMAC);
  }

  Dn.clear();
  DDn.clear();
  Dnp1.clear();
  DDnp1.clear();
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
    const MultiFab& R = get_new_data(RhoYdot_Type);
    const MultiFab& dat = get_new_data(State_Type);
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      FArrayBox enthi, T;

      for (MFIter mfi(R,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();
        T.resize(box,1);
        T.setVal(298.15,box);
        T.copy(dat[mfi],Temp,0,1);

        enthi.resize(box,R.nComp());
        getHGivenT_pphys(enthi,T,box,0,0);
        enthi.mult(R[mfi],box,0,0,R.nComp());
        
        // Form heat release
        (*auxDiag["HEATRELEASE"])[mfi].setVal(0,box);
        for (int j=0; j<R.nComp(); ++j)
        {
          (*auxDiag["HEATRELEASE"])[mfi].minus(enthi,box,j,0,1);
        }
      }
    }
  }
#ifdef BL_COMM_PROFILING
  for (MFIter mfi(*auxDiag["COMMPROF"]); mfi.isValid(); ++mfi)
  {
    int rank(ParallelDescriptor::MyProc());
    (*auxDiag["COMMPROF"])[mfi].setVal(rank, 0);
    (*auxDiag["COMMPROF"])[mfi].setVal(DistributionMapping::ProximityMap(rank),   1);
    (*auxDiag["COMMPROF"])[mfi].setVal(DistributionMapping::ProximityOrder(rank), 2);
  }
#endif

  BL_PROFILE_VAR_START(HTDIFF);
  calcDiffusivity(tnp1);
#ifdef USE_WBAR
  calcDiffusivity_Wbar(tnp1);
#endif

  calcViscosity(tnp1,dt,iteration,ncycle);
  BL_PROFILE_VAR_STOP(HTDIFF);
  //
  // Set the dependent value of RhoRT to be the thermodynamic pressure.  By keeping this in
  // the state, we can use the average down stuff to be sure that RhoRT_avg is avg(RhoRT),
  // not ave(Rho)avg(R)avg(T), which seems to give the p-relax stuff in the mac Rhs troubles.
  //
  BL_PROFILE_VAR_START(HTMAC);
  setThermoPress(tnp1);
  BL_PROFILE_VAR_STOP(HTMAC);

  BL_PROFILE_VAR("HT::advance::project", HTPROJ);
  calc_divu(time+dt, dt, get_new_data(Divu_Type));
  BL_PROFILE_VAR_STOP(HTPROJ);
  showMF("sdc",get_new_data(Divu_Type),"sdc_Divu",level,parent->levelSteps(level));

  BL_PROFILE_VAR_START(HTPROJ);
  if (!NavierStokesBase::initial_step && level != parent->finestLevel())
  {
    //
    // Set new divu to old div + dt*dsdt_old where covered by fine.
    //
    BoxArray crsndgrids = getLevel(level+1).grids;

    crsndgrids.coarsen(fine_ratio);
            
    MultiFab& divu_new = get_new_data(Divu_Type);
    MultiFab& divu_old = get_old_data(Divu_Type);
    MultiFab& dsdt_old = get_old_data(Dsdt_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif            
    for (MFIter mfi(divu_new,true); mfi.isValid();++mfi)
    {
      auto isects = crsndgrids.intersections(mfi.tilebox());

      for (int i = 0, N = isects.size(); i < N; i++)
      {
        const Box& ovlp = isects[i].second;
        divu_new[mfi].copy(dsdt_old[mfi],ovlp,0,ovlp,0,1);
        divu_new[mfi].mult(dt,ovlp,0,1);
        divu_new[mfi].plus(divu_old[mfi],ovlp,0,0,1);
      }
    }
  }
  BL_PROFILE_VAR_STOP(HTPROJ);
  showMF("sdc",get_new_data(Divu_Type),"sdc_DivUnew",level,parent->levelSteps(level));
        
  BL_PROFILE_VAR_START(HTPROJ);
  calc_dsdt(time, dt, get_new_data(Dsdt_Type));
  BL_PROFILE_VAR_STOP(HTPROJ);
  showMF("sdc",get_new_data(Dsdt_Type),"sdc_Dsdtnew",level,parent->levelSteps(level));

  BL_PROFILE_VAR_START(HTPROJ);
  if (NavierStokesBase::initial_step)
    MultiFab::Copy(get_old_data(Dsdt_Type),get_new_data(Dsdt_Type),0,0,1,0);
  BL_PROFILE_VAR_STOP(HTPROJ);
  //
  // Add the advective and other terms to get velocity (or momentum) at t^{n+1}.
  //
  BL_PROFILE_VAR_START(HTVEL);
  if (do_mom_diff == 0) {
    velocity_advection(dt);
  }
  velocity_update(dt);
  BL_PROFILE_VAR_STOP(HTVEL);

  showMF("sdc",get_new_data(State_Type),"sdc_Snew_preProj",level,parent->levelSteps(level));

  // compute chi correction
  // place to take dpdt stuff out of nodal project
  //    calc_dpdt(tnp1,dt,chi_increment,u_mac);
  //    MultiFab::Add(get_new_data(Divu_Type),chi_increment,0,0,1,0);

  // subtract mean from divu
  BL_PROFILE_VAR_START(HTPROJ);
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
  BL_PROFILE_VAR_STOP(HTPROJ);

  showMF("sdc",get_new_data(State_Type),"sdc_Snew_postProj",level,parent->levelSteps(level));

#ifdef AMREX_PARTICLES
  if (theNSPC() != 0)
  {
    theNSPC()->AdvectWithUmac(u_mac, level, dt);
  }
#endif

  BL_PROFILE_VAR("HT::advance::cleanup", HTCLEANUP);
  advance_cleanup(iteration,ncycle);
  BL_PROFILE_VAR_STOP(HTCLEANUP);
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

  temperature_stats(S_new);

  BL_PROFILE_VAR_START(HTMAC);

  // during initialization, reset time 0 ambient pressure
  if (closed_chamber == 1 && level == 0 && !initial_step)
  {
    p_amb_old = p_amb_new;
  }
  BL_PROFILE_VAR_STOP(HTMAC);

  BL_PROFILE_REGION_STOP("R::HT::advance()[src_sdc]");

  return dt_test;
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
  fctmpnew.copy(FC,0,0,1,std::min(ngrow,FC.nGrow()),0);

  int count = 0;
  Vector<long> vwrk(bxba.size());

// FIXME need to be tiled  
  for (MFIter mfi(fctmpnew); mfi.isValid(); ++mfi)
    vwrk[count++] = static_cast<long>(fctmpnew[mfi].sum(0));

  fctmpnew.clear();

#if BL_USE_MPI
  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  Vector<int> nmtags(ParallelDescriptor::NProcs(),0);
  Vector<int> offset(ParallelDescriptor::NProcs(),0);

  const Vector<int>& procmap = rr.ProcessorMap();

  for (int i = 0; i < vwrk.size(); i++)
    nmtags[procmap[i]]++;

  BL_ASSERT(nmtags[ParallelDescriptor::MyProc()] == count);

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
  BL_PROFILE("HT:::advance_chemistry()");

  const Real strt_time = ParallelDescriptor::second();

  const bool do_avg_down_chem = avg_down_chem
    && level < parent->finestLevel()
    && getLevel(level+1).state[RhoYdot_Type].hasOldData();

  if (hack_nochem)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox tmp;
        for (MFIter mfi(mf_old,true); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            tmp.resize(box,nspecies+1);
            const FArrayBox& f = Force[mfi];
            tmp.copy(f,box,0,box,0,nspecies+1);
            tmp.mult(dt,box,0,nspecies+1);
            FArrayBox& Sold = mf_old[mfi];
            FArrayBox& Snew = mf_new[mfi];
            Snew.copy(Sold,box,first_spec,box,first_spec,nspecies+1);
            Snew.plus(tmp,box,box,0,first_spec,nspecies+1);
            Snew.copy(Sold,box,Temp,box,Temp,1);
        }
    }
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
    //
    // Chop the grids to level out the chemistry work.
    // We want enough grids so that KNAPSACK works well,
    // but not too many to make unweildy BoxArrays.
    //
    const int Threshold = chem_box_chop_threshold * ParallelDescriptor::NProcs();
    BoxArray  ba        = mf_new.boxArray();
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

      for (int j = BL_SPACEDIM-1; j >=0 && ba.size() < Threshold; j--)
      {
        chunk[j] /= 2;
        ba.maxSize(chunk);
        if (ba.size() >= Threshold) done = true;
      }
    }

    DistributionMapping dm = getFuncCountDM(ba,ngrow);

    MultiFab diagTemp;
    MultiFab STemp(ba, dm, nspecies+3, 0);
    MultiFab fcnCntTemp(ba, dm, 1, 0);
    MultiFab FTemp(ba, dm, Force.nComp(), 0);

    const bool do_diag = plot_reactions && amrex::intersect(ba,auxDiag["REACTIONS"]->boxArray()).size() != 0;

    if (do_diag)
    {
        diagTemp.define(ba, dm, auxDiag["REACTIONS"]->nComp(), 0);
        diagTemp.copy(*auxDiag["REACTIONS"]); // Parallel copy
    }

    if (verbose) 
      amrex::Print() << "*** advance_chemistry: FABs in tmp MF: " << STemp.size() << '\n';

    STemp.copy(mf_old,first_spec,0,nspecies+3); // Parallel copy.
    FTemp.copy(Force);                          // Parallel copy.

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter Smfi(STemp,true); Smfi.isValid(); ++Smfi)
    {
      const auto&      rhoY     = STemp.array(Smfi);
      const Box&       bx       = Smfi.tilebox();
      const auto&      fcl      = fcnCntTemp.array(Smfi);
      const auto&      frcing   = FTemp.array(Smfi);
      FArrayBox*       chemDiag = (do_diag ? &(diagTemp[Smfi]) : 0);

      Real dt_incr = dt;
      Real time_init = 0;
      int reInit = 1;
      double pressure = 1.0; // dummy FIXME

      const auto len = amrex::length(bx);
      const auto lo  = amrex::lbound(bx);

      double tmp_vect[(nspecies+1)];
      double tmp_src_vect[nspecies];
      double tmp_vect_energy[1];
      double tmp_src_vect_energy[1];

      For(bx,
      [&] AMREX_GPU_DEVICE (int i, int j, int k)
      {
         for (int sp=0;sp<nspecies; sp++){
             tmp_vect[sp]       = rhoY(i,j,k,sp) * 1.e-3;
             tmp_src_vect[sp]   = frcing(i,j,k,sp) * 1.e-3;
         }
         tmp_vect[nspecies]     = rhoY(i,j,k,nspecies+2);
         tmp_vect_energy[0]     = rhoY(i,j,k,nspecies) * 10.0;
         tmp_src_vect_energy[0] = frcing(i,j,k,nspecies) * 10.0;

         fcl(i,j,k) = react(tmp_vect, tmp_src_vect,
                            tmp_vect_energy, tmp_src_vect_energy,
#ifndef USE_SUNDIALS_PP
                            &pressure,
#endif
                            &dt_incr, &time_init);

         dt_incr = dt;
         for (int sp=0;sp<nspecies; sp++){
            rhoY(i,j,k,sp)      = tmp_vect[sp] * 1.e+3;
            if (rhoY(i,j,k,sp) != rhoY(i,j,k,sp)) {
               amrex::Abort("NaNs !! ");
            }
         }
		     rhoY(i,j,k,nspecies+2)  = tmp_vect[nspecies];
         if (rhoY(i,j,k,nspecies+2) != rhoY(i,j,k,nspecies+2)) {
            amrex::Abort("NaNs !! ");
         }
         rhoY(i,j,k,nspecies) = tmp_vect_energy[0] * 1.e-01;
         if (rhoY(i,j,k,nspecies) != rhoY(i,j,k,nspecies)) {
            amrex::Abort("NaNs !! ");
         }

      });

    }

    FTemp.clear();

    mf_new.copy(STemp,0,first_spec,nspecies+3); // Parallel copy.

    STemp.clear();
    //
    // Set React_new (I_R).
    //
    MultiFab::Copy(React_new, mf_old, first_spec, 0, nspecies, 0);

    MultiFab::Subtract(React_new, mf_new, first_spec, 0, nspecies, 0);

    React_new.mult(-1/dt);

    MultiFab::Subtract(React_new, Force, 0, 0, nspecies, 0);

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
      amrex::average_down(fine_React, React_new, getLevel(level+1).geom, geom,
                           0, nspecies, fine_ratio);
    }
    //
    // Ensure consistent grow cells.
    //
    if (ngrow > 0)
    {
      BL_ASSERT(React_new.nGrow() == 1);
      React_new.FillBoundary(0,nspecies, geom.periodicity());
      Extrapolater::FirstOrderExtrap(React_new, geom, 0, nspecies);
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
  BL_PROFILE("HT::comp_sc_adv_fluxes_and_div()");
  //
  // Compute -Div(advective fluxes)  [ which is -aofs in NS, BTW ... careful...
  //
  if (verbose) amrex::Print() << "... computing advection terms\n";

  const Real  strt_time = ParallelDescriptor::second();
  const Real* dx        = geom.CellSize();
  const Real  prev_time = state[State_Type].prevTime();  

  MultiFab dummy(grids,dmap,1,0,MFInfo().SetAlloc(false));
  
  // We advect rho.Y and rho.h
  FillPatchIterator S_fpi(*this,dummy,Godunov::hypgrow(),prev_time,State_Type,first_spec,nspecies+1);
  MultiFab& Smf=S_fpi.get_mf();

  // Floor small values of states to be extrapolated
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Smf,true); mfi.isValid(); ++mfi)
  {
    Box gbx=mfi.growntilebox(Godunov::hypgrow());
    auto fab = Smf.array(mfi);
    AMREX_HOST_DEVICE_FOR_4D ( gbx, nspecies+1, i, j, k, n,
    {
      auto& val = fab(i,j,k,n);
      val = std::abs(val) > 1.e-20 ? val : 0;
    });
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  FArrayBox cflux[BL_SPACEDIM];
  FArrayBox edgstate[BL_SPACEDIM];
  Vector<int> state_bc;
 
  for (MFIter S_mfi(Smf,true); S_mfi.isValid(); ++S_mfi)
  {
    const Box& bx = S_mfi.tilebox();
    const FArrayBox& divu = DivU[S_mfi];
    const FArrayBox& force = Force[S_mfi];

    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      const Box& ebx = amrex::surroundingNodes(bx,d);
      cflux[d].resize(ebx,nspecies+1);
      edgstate[d].resize(ebx,nspecies+1);
    }
     
    (*aofs)[S_mfi].setVal(0,bx,Density,NUM_SCALARS);
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      const Box& ebx = S_mfi.nodaltilebox(d);
      (*EdgeState[d])[S_mfi].setVal(0,ebx,Density,NUM_SCALARS);
      (*EdgeFlux[d])[S_mfi].setVal(0,ebx,Density,NUM_SCALARS);
    } 
        
    state_bc = fetchBCArray(State_Type,bx,first_spec,nspecies+1);

    // Note that the FPU argument is no longer used in IAMR->Godunov.cpp because FPU is now default
    godunov->AdvectScalars(bx, dx, dt, 
                           D_DECL(  area[0][S_mfi],  area[1][S_mfi],  area[2][S_mfi]),
                           D_DECL( u_mac[0][S_mfi], u_mac[1][S_mfi], u_mac[2][S_mfi]),
                           D_DECL(cflux[0],cflux[1],cflux[2]),
                           D_DECL(edgstate[0],edgstate[1],edgstate[2]),
                           Smf[S_mfi], 0, nspecies+1 , force, 0, divu, 0,
                           (*aofs)[S_mfi], first_spec, advectionType, state_bc, FPU, volume[S_mfi]);
       
    // Accumulate rho flux divergence, rho on edges, and rho flux on edges
    
     for (int d=0; d<BL_SPACEDIM; ++d)
     {
       const Box& ebx = S_mfi.nodaltilebox(d);
       (*EdgeFlux[d])[S_mfi].copy(cflux[d],ebx,0,ebx,first_spec,nspecies+1);
       (*EdgeState[d])[S_mfi].copy(edgstate[d],ebx,0,ebx,first_spec,nspecies+1);
     }
          
     for (int comp = 0 ; comp < nspecies+1 ; comp++)
     {      
       int state_ind = first_spec + comp;
       if (state_ind >= first_spec && state_ind <= last_spec)
       {
         (*aofs)[S_mfi].plus((*aofs)[S_mfi],bx,bx,state_ind,Density,1);
         for (int d=0; d<BL_SPACEDIM; d++)
         {
           const Box& ebx = S_mfi.nodaltilebox(d);
           (*EdgeState[d])[S_mfi].plus(edgstate[d],ebx,ebx,comp,Density,1);
           (*EdgeFlux[d])[S_mfi].plus(cflux[d],ebx,ebx,comp,Density,1);
         }
       }
     }

  }    
}

  showMF("sdc",*EdgeState[0],"sdc_ESTATE_x",level,parent->levelSteps(level));
  showMF("sdc",*EdgeState[1],"sdc_ESTATE_y",level,parent->levelSteps(level));
#if BL_SPACEDIM==3
  showMF("sdc",*EdgeState[2],"sdc_ESTATE_z",level,parent->levelSteps(level));
#endif

  // NOTE: Change sense of aofs here so that d/dt ~ aofs...be sure we use our own update function!
  aofs->mult(-1,Density,NUM_SCALARS);

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
      for (int d = 0; d < BL_SPACEDIM; d++)
      {
        advflux_reg->FineAdd((*EdgeFlux[d]),d,Density,Density,NUM_SCALARS,dt);
      }
    }
    if (level < parent->finestLevel())
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
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

  BL_PROFILE("HT::mac_sync()");
  if (verbose) amrex::Print() << "... mac_sync\n";

  const Real strt_time = ParallelDescriptor::second();

  int        sigma;
  const int  finest_level   = parent->finestLevel();
  const int  ngrids         = grids.size();
  const Real prev_time      = state[State_Type].prevTime();
  const Real tnp1       = state[State_Type].curTime();
  const Real prev_pres_time = state[Press_Type].prevTime();
  const Real dt             = parent->dtLevel(level);
  //
  // DeltaSsync Will hold q^{n+1,p} * (delta rho)^sync for conserved quantities
  // as defined before Eq (18) in DayBell:2000.  Note that in the paper, 
  // Eq (18) is missing Y_m^{n+1,p} * (delta rho)^sync in the RHS
  // and Eq (19) is missing the h^{n+1,p} * (delta rho)^sync in the RHS.
  //
  MultiFab*  DeltaSsync = 0;

  MultiFab& S_new = get_new_data(State_Type);

  ////////////////////////
  // save states that we need to reset with each mac sync iteration
  ////////////////////////
  const int numscal = NUM_STATE - BL_SPACEDIM;

  MultiFab chi_sync(grids,dmap,1,0);
  chi_sync.setVal(0);

  Vector<std::unique_ptr<MultiFab> > S_new_sav(finest_level+1);

  for (int lev=level; lev<=finest_level; lev++)
  {
    const MultiFab& S_new_lev = getLevel(lev).get_new_data(State_Type);
    S_new_sav[lev].reset(new MultiFab(S_new_lev.boxArray(),
                                      S_new_lev.DistributionMap(),
                                      NUM_STATE,1));
    MultiFab::Copy(*S_new_sav[lev],S_new_lev,0,0,NUM_STATE,1);
  }

  Vector<std::unique_ptr<MultiFab> > Ssync_sav(finest_level);
  Vector<std::unique_ptr<MultiFab> > Vsync_sav(finest_level);

  for (int lev=level; lev<=finest_level-1; lev++)
  {
    const MultiFab& Ssync_lev = getLevel(lev).Ssync;
    Ssync_sav[lev].reset(new MultiFab(Ssync_lev.boxArray(),
                                      Ssync_lev.DistributionMap(),
                                      numscal,1));
    MultiFab::Copy(*Ssync_sav[lev],Ssync_lev,0,0,numscal,1);

    const MultiFab& Vsync_lev = getLevel(lev).Vsync;
    Vsync_sav[lev].reset(new MultiFab(Vsync_lev.boxArray(),
                                      Vsync_lev.DistributionMap(),
                                      BL_SPACEDIM,1));
    MultiFab::Copy(*Vsync_sav[lev],Vsync_lev,0,0,BL_SPACEDIM,1);
  }

  ////////////////////////
  // begin mac_sync_iter loop here
  ////////////////////////

  // save pressure
  Real p_amb_new_temp = p_amb_new;

  for (int mac_sync_iter=0; mac_sync_iter < num_mac_sync_iter; mac_sync_iter++)
  {
    bool last_mac_sync_iter;
    if (mac_sync_iter == num_mac_sync_iter-1)
    {
      last_mac_sync_iter=true;
    }
    else
    {
      last_mac_sync_iter=false;
    }

    MultiFab chi_sync_increment(grids,dmap,1,0);

    for (int lev=level; lev<=finest_level; lev++)
    {
      MultiFab& S_new_lev = getLevel(lev).get_new_data(State_Type);
      MultiFab::Copy(S_new_lev,*S_new_sav[lev],0,0,NUM_STATE,1);
    }

    for (int lev=level; lev<=finest_level-1; lev++)
    {
      MultiFab& Ssync_lev = getLevel(lev).Ssync;
      MultiFab::Copy(Ssync_lev,*Ssync_sav[lev],0,0,numscal,1);

      MultiFab& Vsync_lev = getLevel(lev).Vsync;
      MultiFab::Copy(Vsync_lev,*Vsync_sav[lev],0,0,BL_SPACEDIM,1);
    }

    make_rho_curr_time();
    get_rho_half_time();

    //
    // Compute the corrective pressure, mac_sync_phi, used to 
    // compute U^{ADV,corr} in mac_sync_compute
    //
    bool subtract_avg = false;
    if (closed_chamber && level == 0)
    {
      subtract_avg = true;
    }
    Real offset = 0.0;

    BL_PROFILE_VAR("HT::mac_sync::ucorr", HTUCORR);
    mac_projector->mac_sync_solve(level,dt,rho_half,fine_ratio,
                                  &chi_sync,subtract_avg,offset);
    BL_PROFILE_VAR_STOP(HTUCORR);

    if (closed_chamber && level == 0)
    {
      Print() << "0-1 MAC SYNC OFFSET " << offset << std::endl;
      Print() << "0-1 MAC SYNC PAMB ADJUSTMENT " << -dt*(offset/thetabar) << std::endl;
      p_amb_new = p_amb_new_temp - dt*(offset/thetabar);
      p_amb_old = p_amb_new;
    }

    if (!do_reflux) return;

    Vector<SYNC_SCHEME> sync_scheme(NUM_STATE,ReAdvect);

    if (do_mom_diff == 1)
      for (int i=0; i<BL_SPACEDIM; ++i)
        sync_scheme[i] = UseEdgeState;

    for (int i=BL_SPACEDIM; i<NUM_STATE; ++i)
      sync_scheme[i] = UseEdgeState;
        
    Vector<int> incr_sync(NUM_STATE,0);
    for (int i=0; i<sync_scheme.size(); ++i)
      if (sync_scheme[i] == ReAdvect)
        incr_sync[i] = 1;
    //
    // After solving for mac_sync_phi in mac_sync_solve(), we
    // can now do the sync advect step in mac_sync_compute().
    // This consists of two steps
    //
    // 1. compute U^{ADV,corr} as the gradient of mac_sync_phi
    // 2. add -D^MAC ( U^{ADV,corr} * rho * q)^{n+1/2} ) to flux registers,
    //    which already contain the delta F adv/diff flux mismatches
    //

    //
    // velocities
    //
    BL_PROFILE_VAR("HT::mac_sync::Vsync", HTVSYNC);
    if (do_mom_diff == 0) 
    {
      mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,rho_half,
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
      for (int comp=0; comp<BL_SPACEDIM; ++comp)
      {
        if (sync_scheme[comp]==UseEdgeState)
        {
          mac_projector->mac_sync_compute(level,Vsync,comp,
                                          comp,EdgeState, comp,rho_half,
                                          (level > 0 ? &getAdvFluxReg(level):0),
                                          advectionType,modify_reflux_normal_vel,dt,
                                          last_mac_sync_iter);
        }
      }
    }
    BL_PROFILE_VAR_STOP(HTVSYNC);

    showMF("sdcSync",Ssync,"sdc_Ssync_no_Ucorr",level,parent->levelSteps(level));
    showMF("sdcSync",*EdgeState[0],"sdc_ESTATE_x_no_Ucorr",level,parent->levelSteps(level));
    showMF("sdcSync",*EdgeState[1],"sdc_ESTATE_y_no_Ucorr",level,parent->levelSteps(level));
    //
    // Scalars.
    //
    BL_PROFILE_VAR("HT::mac_sync::Ssync", HTSSYNC);
    for (int comp=BL_SPACEDIM; comp<NUM_STATE; ++comp)
    {
      if (sync_scheme[comp]==UseEdgeState)
      {
        int s_ind = comp - BL_SPACEDIM;
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
        mac_projector->mac_sync_compute(level,Ssync,comp,s_ind,
                                        EdgeState,comp,rho_half,
                                        (level > 0 ? &getAdvFluxReg(level):0),
                                        advectionType,modify_reflux_normal_vel,dt,
                                        last_mac_sync_iter);
      }
    }
    BL_PROFILE_VAR_STOP(HTSSYNC);
        
    showMF("sdcSync",Ssync,"sdc_Ssync_after_Ucorr",level,parent->levelSteps(level));
    showMF("sdcSync",*EdgeState[0],"sdc_ESTATE_x_no_Ucorr",level,parent->levelSteps(level));
    showMF("sdcSync",*EdgeState[1],"sdc_ESTATE_y_no_Ucorr",level,parent->levelSteps(level));

    BL_PROFILE_VAR_START(HTSSYNC);
    Ssync.mult(dt,Ssync.nGrow());

    sync_setup(DeltaSsync);
    //
    // For all conservative variables Q (other than density)
    // set DeltaSsync = q^{n+1,p} * (delta rho)^sync,
    // then subtract DeltaSsync from Ssync
    // (these are the terms that were accidentally omitted in (18) and (19)
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      FArrayBox delta_ssync;

      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
      {
        const Box& grd = mfi.tilebox();

        int iconserved = -1;

        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
          if (istate != Density && advectionType[istate] == Conservative)
          {
            iconserved++;
            delta_ssync.resize(grd,1);
            delta_ssync.copy(S_new[mfi],grd,istate,grd,0,1); // delta_ssync = (rho*q)^{n+1,p}
            delta_ssync.divide(S_new[mfi],grd,Density,0,1); // delta_ssync = q^{n+1,p}
            FArrayBox& s_sync = Ssync[mfi]; // Ssync = RHS of Eq (18), (19) without the q^{n+1,p} * (delta rho)^sync terms
            // note that the density component contains (delta rho)^sync
            delta_ssync.mult(s_sync,grd,Density-BL_SPACEDIM,0,1); // delta_ssync = q^{n+1,p} * (delta rho)^sync
            (*DeltaSsync)[mfi].copy(delta_ssync,grd,0,grd,iconserved,1); // DeltaSsync = q^{n+1,p} * (delta rho)^sync
            s_sync.minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1); // Ssync = Ssync - q^{n+1,p} * (delta rho)^sync
          }
        }
      }

      delta_ssync.clear();
    }

    //
    // Now, increment density.
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      S_new[mfi].plus(Ssync[mfi],box,Density-BL_SPACEDIM,Density,1);
    }

    make_rho_curr_time();
    BL_PROFILE_VAR_STOP(HTSSYNC);
    //
    // Set do_diffuse_sync to 0 for debugging reasons only.
    //
    BL_PROFILE_VAR_START(HTVSYNC);
    if (do_mom_diff == 1)
    {
      for (int d=0; d<BL_SPACEDIM; ++d) {
        MultiFab::Divide(Vsync,rho_ctime,0,Xvel+d,1,0); 
      }
    }
    BL_PROFILE_VAR_STOP(HTVSYNC);
      
    BL_PROFILE_VAR_START(HTSSYNC);
    if (do_diffuse_sync)
    {
      FluxBoxes fb_beta(this);
      MultiFab** beta = fb_beta.get();
      if (is_diffusive[Xvel])
      {
        int rho_flag = (do_mom_diff == 0) ? 1 : 3;
        getViscosity(beta, tnp1);
        diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,rho_half,rho_flag,beta,0,
                                 last_mac_sync_iter);
      }
	    
      if (!unity_Le 
          && nspecies>0 
          && do_add_nonunityLe_corr_to_rhoh_adv_flux) 
      {
        //
        // Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
        // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
        // Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync
        //
        differential_spec_diffuse_sync(dt,false,last_mac_sync_iter);

#ifdef USE_WBAR
        for (int dir=0; dir<BL_SPACEDIM; ++dir)
        {
          (*SpecDiffusionFluxWbar[dir]).setVal(0.);
        }

        // compute beta grad Wbar terms using the n+1,p state
        // we want this to have a negative sign since we will add
        // the beta grad Wbar terms using the updated state later to create the delta
        // store in SpecDiffusionFluxWbar
        compute_Wbar_fluxes(tnp1,-1);

        // create an updated (but stil temporary) new state species
        // by adding Ssync plus DeltaSsync to the new state
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        for (MFIter mfi(Ssync,true); mfi.isValid(); ++mfi)
        {
          const Box& bx = mfi.tilebox();
		
          int iconserved = -1;
		
          for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
          {
            if (istate != Density && advectionType[istate] == Conservative)
            {
              iconserved++;
	
              // only actually update the species
              // this is the Y_m^{n+1,p} * delta rho^sync piece
              if (istate >= first_spec && istate <= last_spec)
              {
                Ssync[mfi].plus((*DeltaSsync)[mfi],bx,iconserved,istate-BL_SPACEDIM,1);
              }
            }
          }
        }
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
          const Box& bx = mfi.tilebox();
		
          for (int sigma = 0; sigma < numscal; sigma++)
          {
            // only actually update the species
            // this is the rho^{n+1} * delta Y^sync piece
            if (BL_SPACEDIM+sigma >= first_spec && BL_SPACEDIM+sigma <= last_spec)
            {
              S_new[mfi].plus(Ssync[mfi],bx,sigma,BL_SPACEDIM+sigma,1);
            }
          }
        }
}
        // compute beta grad Wbar terms using the temporary new state
        // add these to SpecDiffusionFluxWbar so it contains the delta
        compute_Wbar_fluxes(tnp1,1);
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        // subtract Ssync and DeltaSsync from the new state to restore the n+1,p state
        for (MFIter mfi(Ssync,true); mfi.isValid(); ++mfi)
        {
          const Box& bx = mfi.tilebox();
		
          int iconserved = -1;
		
          for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
          {
            if (istate != Density && advectionType[istate] == Conservative)
            {
              iconserved++;
	
              // only actually update the species
              // this is the Y_m^{n+1,p} * delta rho^sync piece
              if (istate >= first_spec && istate <= last_spec)
              {
                Ssync[mfi].minus((*DeltaSsync)[mfi],bx,iconserved,istate-BL_SPACEDIM,1);
              }
            }
          }
        }
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
          const Box& bx = mfi.tilebox();
		
          for (int sigma = 0; sigma < numscal; sigma++)
          {
            // only actually update the species
            // this is the rho^{n+1} * delta Y^sync piece
            if (BL_SPACEDIM+sigma >= first_spec && BL_SPACEDIM+sigma <= last_spec)
            {
              S_new[mfi].minus(Ssync[mfi],bx,sigma,BL_SPACEDIM+sigma,1);
            }
          }
        }
}

        // take divergence of beta grad delta Wbar and multiply divergence by dt/2
        MultiFab DWbar(grids,dmap,nspecies,nGrowAdvForcing);
        MultiFab* const * fluxWbar = SpecDiffusionFluxWbar;
        flux_divergence(DWbar,0,fluxWbar,0,nspecies,-1);
        DWbar.mult(dt/2.0);
        
        
        // reset Ssync to be the same RHS as above, but with the (dt/2) div beta grad delta Wbar term
        // use the code above, but add on the grad delta Wbar term
#ifdef _OPENMP
#pragma omp parallel
#endif
{
        FArrayBox delta_ssync;
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
          const int  i   = mfi.index();
          const Box& grd = grids[i];
          const FArrayBox& DWbarFab = DWbar[mfi];
          int iconserved = -1;
		
          for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
          {
            if (istate != Density && advectionType[istate] == Conservative)
            {
              iconserved++;
              if (istate >= first_spec && istate <= last_spec)
              {
                delta_ssync.resize(grd,1);
                delta_ssync.copy(S_new[mfi],grd,istate,grd,0,1); // delta_ssync = (rho*q)^{n+1,p}
                delta_ssync.divide(S_new[mfi],grd,Density,0,1); // delta_ssync = q^{n+1,p}
                FArrayBox& s_sync = Ssync[mfi]; // Ssync = RHS of Eq (18), (19) without the q^{n+1,p} * (delta rho)^sync terms
                delta_ssync.mult(s_sync,grd,Density-BL_SPACEDIM,0,1); // delta_ssync = q^{n+1,p} * (delta rho)^sync
                (*DeltaSsync)[mfi].copy(delta_ssync,grd,0,grd,iconserved,1); // DeltaSsync = q^{n+1,p} * (delta rho)^sync
                s_sync.minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1); // Ssync = Ssync - q^{n+1,p} * (delta rho)^sync
                s_sync.plus(DWbarFab,grd,istate-first_spec,istate-BL_SPACEDIM,1); // add grad delta Wbar terms
              }
            }
          }
        }
        delta_ssync.clear();
}
        

        // call differential_spec_diffuse_sync again, but this time the conservative
        // correction needs to be the sum of the delta Y_m and SpecDiffusionFluxWbar terms
        // FIXME
        differential_spec_diffuse_sync(dt, true, last_mac_sync_iter);

#endif

        const Real a = 1.0;     // Passed around, but not used
        Real rhsscale;          //  -ditto-
        const int rho_flag = 2; // FIXME: Messy assumption
        MultiFab *alpha=0;      //  -ditto-
        FluxBoxes fb_SC  (this, 1, 0);
        FluxBoxes fb_NULN(this, nspecies, 0);
        FluxBoxes fb_visc(this, 1, 0);
        MultiFab **fluxSC    =   fb_SC.get();
        MultiFab **fluxNULN  = fb_NULN.get();
        MultiFab **rhoh_visc = fb_visc.get();
	  
        const int nGrow    = 1; // Size to grow fil-patched fab for T below
        const int dataComp = 0; // coeffs loaded into 0-comp for all species
                  
        getDiffusivity(rhoh_visc, tnp1, RhoH, 0, 1); // RhoH (lambda/cp)
	  
        MultiFab Soln(grids,dmap,1,1);

        // compute lambda/cp grad (delta Y_m^sync)
        for (int comp = 0; comp < nspecies; ++comp)
        {
          const Real b     = be_cn_theta;
          const int  sig = first_spec + comp;
          //
          //  start by getting lambda/cp.Grad(delta Y^sync)
          //   (note: neg of usual diff flux)
          //
          ViscBndry      visc_bndry;
          ABecLaplacian* visc_op;

          visc_op = diffusion->getViscOp(sig,a,b,tnp1,
                                         visc_bndry,rho_half,
                                         rho_flag,&rhsscale,rhoh_visc,dataComp,
                                         alpha,dataComp);

          visc_op->maxOrder(diffusion->maxOrder());

          // copy rho^{n+1} * (delta Y)^sync into Soln
          MultiFab::Copy(Soln,Ssync,sig-BL_SPACEDIM,0,1,0);

          // divide Soln by rho
#ifdef _OPENMP
#pragma omp parallel
#endif
          for (MFIter Smfi(Soln,true); Smfi.isValid(); ++Smfi)
          {
            Soln[Smfi].divide(S_new[Smfi],Smfi.tilebox(),Density,0,1);
          }

          // compute lambda/cp.Grad(delta Y^sync) and weight
          visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);

          for (int d = 0; d < BL_SPACEDIM; ++d)
            fluxSC[d]->mult(-b/geom.CellSize()[d]);
          //
          // Here, get fluxNULN = (lambda/cp - rho.D)Grad(delta Ysync)
          //                    = lambda/cp.Grad(delta Ysync) + SpecSyncDiffFlux
          //
          for (int d = 0; d < BL_SPACEDIM; ++d)
          {
#ifdef _OPENMP
#pragma omp parallel
#endif	      
            for (MFIter SDF_mfi(*SpecDiffusionFluxnp1[d],true) ; SDF_mfi.isValid(); ++SDF_mfi)
            {
              FArrayBox& fluxSC_fab   = (*fluxSC[d])[SDF_mfi];
              FArrayBox& fluxNULN_fab = (*fluxNULN[d])[SDF_mfi];
              FArrayBox& SDF_fab = (*SpecDiffusionFluxnp1[d])[SDF_mfi];
              const Box& ebox    = SDF_mfi.tilebox();
              // copy in (delta Gamma)^sync
              fluxNULN_fab.copy(SDF_fab,ebox,comp,ebox,comp,1);
              // add in (lambda/cp) grad (delta Y^sync)
              fluxNULN_fab.plus(fluxSC_fab,ebox,0,comp,1);
            }
          }
          delete visc_op;
        }

        Soln.clear();

        fb_SC.clear();
        fb_visc.clear();
        //
        // Multiply fluxi by h_i (let FLXDIV routine below sum up the fluxes)
        //
        //FArrayBox eTemp, h;

        FillPatchIterator Tnew_fpi(*this,S_new,nGrow,tnp1,State_Type,Temp,1);
        const MultiFab& Tnew_mf = Tnew_fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
{
        FArrayBox eTemp, h;

    for (MFIter mfi(Tnew_mf,true); mfi.isValid(); ++mfi)
        {
          const Box& box = mfi.tilebox();

          for (int d = 0; d < BL_SPACEDIM; ++d)
          {
            const Box& ebox = mfi.nodaltilebox(d);
            eTemp.resize(ebox,1);

            FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
            FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
            center_to_edge_fancy(Tnew_mf[mfi],eTemp,amrex::grow(box,amrex::BASISV(d)),ebox,
                                 0,0,1,geom.Domain(),bc_lo,bc_hi);

            h.resize(ebox,nspecies);
            getHGivenT_pphys(h,eTemp,ebox,0,0);

            // multiply fluxNULN by h_m
            (*fluxNULN[d])[mfi].mult(h,ebox,0,0,nspecies);
          }
        }
}
        //h.clear();
        //eTemp.clear();

        //
        // Add the NULN fluxes to the RHS of the (delta h)^sync diffusion solve
        // afterwards, the entire RHS should be ready.
        //
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter Ssync_mfi(Ssync,true); Ssync_mfi.isValid(); ++Ssync_mfi)
        {
          FArrayBox& syncn = Ssync[Ssync_mfi];
          const FArrayBox& synco = Ssync[Ssync_mfi];
          const Box& box = Ssync_mfi.tilebox();

          //
          // Multiply by dt*dt, one to make it extensive, and one because
          // Ssync multiplied above by dt, need same units here.
          //
          const Real mult = dt*dt;
          const int sigmaRhoH = RhoH - BL_SPACEDIM; // RhoH comp in Ssync

          incrwext_flx_div( BL_TO_FORTRAN_BOX(box),
                            BL_TO_FORTRAN_ANYD((*fluxNULN[0])[Ssync_mfi]),
                            BL_TO_FORTRAN_ANYD((*fluxNULN[1])[Ssync_mfi]),
#if ( AMREX_SPACEDIM == 3 )
                            BL_TO_FORTRAN_ANYD((*fluxNULN[2])[Ssync_mfi]),
#endif
                            BL_TO_FORTRAN_N_ANYD(synco,sigmaRhoH),
                            BL_TO_FORTRAN_N_ANYD(syncn,sigmaRhoH),
                            BL_TO_FORTRAN_ANYD(volume[Ssync_mfi]),
                            &nspecies, &mult);
        }
      }

      FluxBoxes fb_flux(this);
      MultiFab **flux = fb_flux.get();

      showMF("sdcSync",Ssync,"sdc_Sync_preDiff",level,parent->levelSteps(level));

      for (sigma = 0; sigma < numscal; sigma++)
      {
        int rho_flag = 0;
        int do_viscsyncflux = do_reflux;
        const int state_ind = BL_SPACEDIM + sigma;
        //
        // To diffuse, or not?
        // (1) Density, no
        // (2) RhoH...if diffusive
        // (3) Trac...if diffusive
        // (4) Spec:
        //    (a) if Le==1, and spec diffusive
        //    (b) if Le!=1, do differential diffusion instead, done above
        // (5) Temp, no (set instead by RhoH to Temp)
        //
        const bool is_spec = state_ind<=last_spec && state_ind>=first_spec;
        int do_it
          =  state_ind!=Density 
          && state_ind!=Temp
          && is_diffusive[state_ind]
          && !(is_spec && !unity_Le);

        if (do_it && (is_spec || state_ind==RhoH))
          rho_flag = 2;

        if (do_it)
        {
          MultiFab* alpha = 0;
          getDiffusivity(beta, tnp1, state_ind, 0, 1);
                    
          // on entry, Ssync = RHS for (delta h)^sync diffusive solve
          // on exit, Ssync = rho^{n+1} * (delta h)^sync
          // on exit, flux = coeff * grad phi
          diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,rho_half,
                                   rho_flag,flux,0,beta,0,alpha,0);
          if (do_viscsyncflux && level > 0 && last_mac_sync_iter)
          {
            for (int d=0; d<BL_SPACEDIM; ++d)
              getViscFluxReg().FineAdd(*flux[d],d,0,state_ind,1,dt);
          }
        }
      }
    }
    BL_PROFILE_VAR_STOP(HTSSYNC);
    showMF("sdcSync",Ssync,"sdc_Sync_postDiff",level,parent->levelSteps(level));

    //
    // For all conservative variables Q (other than density)
    // increment sync by (sync_for_rho)*q_presync.
    // Before this loop, Ssync holds rho^{n+1} (delta phi)^sync
    // DeltaSsync holds (delta rho)^sync phi^p
    //
    BL_PROFILE_VAR_START(HTSSYNC);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Ssync,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
	  
      int iconserved = -1;
	  
      for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
      {
        if (istate != Density && advectionType[istate] == Conservative)
        {
          iconserved++;
          Ssync[mfi].plus((*DeltaSsync)[mfi],box,iconserved,istate-BL_SPACEDIM,1);
        }
      }
    }
    sync_cleanup(DeltaSsync);
    BL_PROFILE_VAR_STOP(HTSSYNC);
    showMF("sdcSync",Ssync,"sdc_Sync_post_add_QdRho",level,parent->levelSteps(level));
    //
    // Increment the state (for all but rho, since that was done above)
    //
    BL_PROFILE_VAR_START(HTSSYNC);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();

      for (int sig = 0; sig < numscal; sig++)
      {
        if (!(BL_SPACEDIM+sig == Density))
        {
          S_new[mfi].plus(Ssync[mfi],box,sig,BL_SPACEDIM+sig,1);
        }
      }
    }

    //
    // Recompute temperature and rho R T after the mac_sync.
    //
    RhoH_to_Temp(S_new);
    setThermoPress(tnp1);
    BL_PROFILE_VAR_STOP(HTSSYNC);
      
    showMF("sdcSync",S_new,"sdc_Snew_postSync",level,parent->levelSteps(level));
    //
    // Get boundary conditions.
    //
    BL_PROFILE_VAR_START(HTSSYNC);
    Real mult = 1.0;
    Vector<int*>         sync_bc(grids.size());
    Vector< Vector<int> > sync_bc_array(grids.size());
    for (int i = 0; i < ngrids; i++)
    {
      sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
      sync_bc[i]       = sync_bc_array[i].dataPtr();
    }
    //
    // Interpolate the sync correction to the finer levels.
    //
    IntVect ratio = IntVect::TheUnitVector();
      
    for (int lev = level+1; lev <= finest_level; lev++)
    {
      ratio                   *= parent->refRatio(lev-1);
      PeleLM& fine_level = getLevel(lev);
      MultiFab& S_new_lev      = fine_level.get_new_data(State_Type);
      //
      // New way of interpolating syncs to make sure mass is conserved
      // and to ensure freestream preservation for species & temperature.
      //
      const BoxArray& fine_grids = S_new_lev.boxArray();
      const DistributionMapping& fine_dmap  = S_new_lev.DistributionMap();
      const int nghost           = S_new_lev.nGrow();
      MultiFab increment(fine_grids, fine_dmap, numscal, nghost);
      increment.setVal(0,nghost);
      //
      // Note: we use the lincc_interp (which_interp==3) for density,
      // rho*h and rho*Y, cell_cons_interp for everything else. Doing
      // so is needed for freestream preservation of Y and T.  The setting
      // which_interp=5 calls the lincc_interp, but then follows with a
      // "protection" step to be sure that all the components but the
      // first and the last have their sync adjusted to try to preserve
      // positivity after the sync is applied.  The density sync is then
      // adjusted to be the sum of the species syncs.
      //
      // HACK note: Presently, the species mass syncs are redistributed 
      //            without consequence to the enthalpy sync.  Clearly
      //            the species carry enthalphy, so the enthalph sync should
      //            be adjusted as well.  Note yet sure how to do this correctly.
      //            Punt for now...
      //
      const SyncInterpType which_interp = CellConsProt_T;

      const int nComp = 2+nspecies;

      SyncInterp(Ssync, level, increment, lev, ratio, 
                 Density-BL_SPACEDIM, Density-BL_SPACEDIM, nComp, 1, mult, 
                 sync_bc.dataPtr(), which_interp, Density);

      if (have_trac)
        SyncInterp(Ssync, level, increment, lev, ratio, 
                   Trac-BL_SPACEDIM, Trac-BL_SPACEDIM, 1, 1, mult, 
                   sync_bc.dataPtr());
	
      if (have_rhort)
        SyncInterp(Ssync, level, increment, lev, ratio, 
                   RhoRT-BL_SPACEDIM, RhoRT-BL_SPACEDIM, 1, 1, mult, 
                   sync_bc.dataPtr());
	
      SyncInterp(Ssync, level, increment, lev, ratio, 
                 Temp-BL_SPACEDIM, Temp-BL_SPACEDIM, 1, 1, mult, 
                 sync_bc.dataPtr());
      if (do_set_rho_to_species_sum)
      {
        increment.setVal(0,Density-BL_SPACEDIM,1,0);

        for (int istate = first_spec; istate <= last_spec; istate++)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
          for (MFIter mfi(increment,true); mfi.isValid(); ++mfi)
          {
            const Box& box = mfi.tilebox();
            increment[mfi].plus(increment[mfi],box,
                                istate-BL_SPACEDIM,Density-BL_SPACEDIM,1);
          }
        }
      }
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(increment,true); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.tilebox();
        S_new_lev[mfi].plus(increment[mfi],box,0,Density,numscal);
      }

      if (last_mac_sync_iter)
      {
        fine_level.make_rho_curr_time();
        fine_level.incrRhoAvg(increment,Density-BL_SPACEDIM,1.0);
      }
      //
      // Recompute temperature and rho R T after interpolation of the mac_sync correction
      //   of the individual quantities rho, Y, T.
      //
      RhoH_to_Temp(S_new_lev);
      fine_level.setThermoPress(tnp1);
    }

    //
    // Average down Trac = rho R T after interpolation of the mac_sync correction
    //   of the individual quantities rho, Y, T.
    //
    for (int lev = finest_level-1; lev >= level; lev--)
    {
      PeleLM& fine_level = getLevel(lev+1);
      PeleLM& crse_level = getLevel(lev);

      MultiFab& S_crse_loc = crse_level.get_new_data(State_Type);
      MultiFab& S_fine_loc = fine_level.get_new_data(State_Type);

      const int pComp = (have_rhort ? RhoRT : Trac);
      amrex::average_down(S_fine_loc, S_crse_loc, fine_level.geom, crse_level.geom,
                           pComp, 1, crse_level.fine_ratio);
    }
    BL_PROFILE_VAR_STOP(HTSSYNC);
    showMF("sdcSync",S_new,"sdc_Snew_postInterpAvgSync",level,parent->levelSteps(level));

    chi_sync_increment.setVal(0,0);
    calc_dpdt(tnp1,dt,chi_sync_increment,u_mac);

    MultiFab::Subtract(chi_sync,chi_sync_increment,0,0,1,0);

  }  // end loop over mac_sync_iters

  if (verbose)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;
      
    ParallelDescriptor::ReduceRealMax(run_time,IOProc);
      
    amrex::Print() << "PeleLM::mac_sync(): lev: " << level << ", time: " << run_time << '\n';
  }

}

#ifdef USE_WBAR
void
PeleLM::compute_Wbar_fluxes(Real time,
                            Real inc)
{
  BL_PROFILE("HT::compute_Wbar_fluxes()");

  // allocate edge-beta for Wbar
  FluxBoxes fb_betaWbar(this, nspecies, 0);
  MultiFab** betaWbar =fb_betaWbar.get();

  // average transport coefficients for Wbar to edges
  getDiffusivity_Wbar(betaWbar,time);

  int nGrowOp = 1;

  MultiFab rho_and_species(grids,dmap,nspecies+1,nGrowOp);

  FillPatchIterator fpi(*this,rho_and_species,nGrowOp,time,State_Type,Density,nspecies+1);
  MultiFab& mf= fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  FArrayBox tmp;
  for (MFIter mfi(mf,true); mfi.isValid();++mfi)
  {
    const Box& box = mfi.growntilebox(); 
    FArrayBox& rho_and_spec = rho_and_species[mfi];
    FArrayBox& Umf = mf[mfi];
    rho_and_spec.copy(Umf,box,0,box,0,nspecies+1);

    tmp.resize(box,1);
    tmp.copy(rho_and_spec,0,0,1);
    tmp.invert(1,box);

    for (int comp = 0; comp < nspecies; ++comp) 
      rho_and_spec.mult(tmp,box,0,comp+1,1);
  }
}

  // add in grad wbar term
  MultiFab Wbar;

  Wbar.define(grids,dmap,1,nGrowOp);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rho_and_species,true); mfi.isValid(); ++mfi)
  {
    const Box& gbox = mfi.growntilebox(); // amrex::grow(mfi.tilebox(),nGrowOp);
    getMwmixGivenY_pphys(Wbar[mfi],rho_and_species[mfi],gbox,1,0);
  }

  //
  // Here, we'll use the same LinOp as Y for filling grow cells.
  //
  const Real* dx    = geom.CellSize();
  const BCRec& bc = get_desc_lst()[State_Type].getBC(first_spec);
  ViscBndry*     bndry   = new ViscBndry(grids,dmap,1,geom);
  ABecLaplacian* visc_op = new ABecLaplacian(bndry,dx);

  visc_op->maxOrder(diffusion->maxOrder());

  MultiFab rho_and_species_crse;
  const int nGrowCrse = InterpBndryData::IBD_max_order_DEF - 1;

  if (level == 0)
  {
    bndry->setBndryValues(Wbar,0,0,1,bc);
  }
  else
  {
    PeleLM& coarser = *(PeleLM*) &(parent->getLevel(level-1));

    rho_and_species_crse.define(coarser.grids,coarser.dmap,nspecies+1,nGrowCrse);

    FillPatchIterator fpi(coarser,rho_and_species_crse,nGrowCrse,time,State_Type,Density,nspecies+1);
    MultiFab& mf= fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  FArrayBox tmp;
  for (MFIter mfi(mf,true); mfi.isValid();++mfi)
    {
      const Box& box = mfi.growntilebox(); 
      FArrayBox& Umf = mf[mfi];
      FArrayBox& rho_and_spec = rho_and_species_crse[mfi];
      rho_and_spec.copy(Umf,box,0,box,0,nspecies+1);
      tmp.resize(box,1);
      tmp.copy(rho_and_spec,box,0,box,0,1);
      tmp.invert(1,box);
	
      for (int comp = 0; comp < nspecies; ++comp) 
        rho_and_spec.mult(tmp,box,0,comp+1,1);
    }
}
    BoxArray cgrids = grids;
    cgrids.coarsen(crse_ratio);
    BndryRegister crse_br(cgrids,dmap,0,1,nGrowCrse,1);
    crse_br.setVal(1.e200);
    MultiFab Wbar_crse(rho_and_species_crse.boxArray(),
                       rho_and_species_crse.DistributionMap(),
                       1,nGrowCrse);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho_and_species_crse,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.growntilebox();
      getMwmixGivenY_pphys(Wbar_crse[mfi],rho_and_species_crse[mfi],box,1,0);
    }	  
    crse_br.copyFrom(Wbar_crse,nGrowCrse,0,0,1);
    bndry->setBndryValues(crse_br,0,Wbar,0,0,1,crse_ratio,bc);
  }

  visc_op->setScalars(0,1);
  visc_op->bCoefficients(1);
  visc_op->applyBC(Wbar,0,1,0,LinOp::Inhomogeneous_BC);

  delete visc_op;

#ifdef _OPENMP
#pragma omp parallel
#endif   
  for (MFIter mfi(Wbar,true); mfi.isValid(); ++mfi)
  {
    
    const FArrayBox& wbar = Wbar[mfi];
    const Real       mult = -1.0;

    for (int d=0; d<AMREX_SPACEDIM; ++d) 
    {
      const Box&        vbox = mfi.nodaltilebox(d);
      const FArrayBox& rhoDe = (*betaWbar[d])[mfi];
      FArrayBox&    fluxfab  = (*SpecDiffusionFluxWbar[d])[mfi];

      for (int ispec=0; ispec<nspecies; ++ispec)
      {
        grad_wbar(BL_TO_FORTRAN_BOX(vbox),
                  BL_TO_FORTRAN_ANYD(wbar),
                  BL_TO_FORTRAN_ANYD(rhoDe),
                  BL_TO_FORTRAN_N_ANYD(fluxfab,ispec),
                  BL_TO_FORTRAN_ANYD(area[d][mfi]),
                  &dx[d], &d, &mult, &inc);
      }
    }
  }
}

#endif

void
PeleLM::differential_spec_diffuse_sync (Real dt,
                                        bool Wbar_corrector,
                                        bool last_mac_sync_iter)
{
  BL_PROFILE_REGION_START("R::HT::differential_spec_diffuse_sync()");
  BL_PROFILE("HT::differential_spec_diffuse_sync()");
  
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
  FluxBoxes fb_betanp1(this, nspecies, 0);
  MultiFab **betanp1 = fb_betanp1.get();
  getDiffusivity(betanp1, tnp1, first_spec, 0, nspecies); // species

  MultiFab Rhs(grids,dmap,nspecies,0);
  const int spec_Ssync_sComp = first_spec - BL_SPACEDIM;

  //
  // Rhs and Ssync contain the RHS of DayBell:2000 Eq (18)
  // with the additional -Y_m^{n+1,p} * (delta rho)^sync term
  // Copy this into Rhs; we will need this later since we overwrite SSync
  // in the solves.
  //
  MultiFab::Copy(Rhs,Ssync,spec_Ssync_sComp,0,nspecies,0);
  //
  // Some standard settings
  //
  const Vector<int> rho_flag(nspecies,2);
  const MultiFab* alpha = 0;

  FluxBoxes fb_fluxSC(this);
  MultiFab** fluxSC = fb_fluxSC.get();

  const MultiFab& RhoHalftime = get_rho_half_time();

  for (int sigma = 0; sigma < nspecies; ++sigma)
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
    diffusion->diffuse_Ssync(Ssync,ssync_ind,dt,be_cn_theta,
                             RhoHalftime,rho_flag[sigma],fluxSC,0,
                             betanp1,sigma,alpha,0);
    //
    // Pull fluxes into flux array
    // this is the rhoD grad (delta Ytilde)^sync terms in DayBell:2000 Eq (18)
    //
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSC[d],0,sigma,1,0);
    }
  }
  fb_betanp1.clear();
  fb_fluxSC.clear();
  //
  // Modify update/fluxes to preserve flux sum = 0
  // (Be sure to pass the "normal" looking Rhs to this generic function)
  //

#ifdef USE_WBAR

  // if in the Wbar corrector, add the grad delta Wbar fluxes
  if (Wbar_corrector)
  {
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(*SpecDiffusionFluxWbar[d],true); mfi.isValid(); ++mfi)
      {
        const Box& ebox =mfi.tilebox();
        (*SpecDiffusionFluxnp1[d])[mfi].plus((*SpecDiffusionFluxWbar[d])[mfi],ebox,0,0,nspecies);
      }
    }
  }

#endif

  // need to correct SpecDiffusionFluxnp1 to contain rhoD grad (delta Y)^sync
  adjust_spec_diffusion_fluxes(tnp1);
  //
  // Need to correct Ssync to contain rho^{n+1} * (delta Y)^sync.
  // Do this by setting
  // Ssync = "RHS from diffusion solve" + (dt/2)*div(delta Gamma)
  //
  // Recompute update with adjusted diffusion fluxes
  //
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  FArrayBox update, efab[BL_SPACEDIM];

  for (MFIter mfi(Ssync,true); mfi.isValid(); ++mfi)
  {
    int        iGrid = mfi.index();
    const Box& box   = mfi.tilebox();

    // copy corrected (delta gamma) on edges into efab
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
      const Box& ebox = amrex::surroundingNodes(box,d);
      efab[d].resize(ebox,nspecies);
            
      efab[d].copy((*SpecDiffusionFluxnp1[d])[mfi],ebox,0,ebox,0,nspecies);
      efab[d].mult(be_cn_theta,ebox,0,nspecies);
    }

    update.resize(box,nspecies);
    update.setVal(0);

    // is this right? - I'm not sure what the scaling on SpecDiffusionFluxnp1 is
    Real scale = -dt;

    // take divergence of (delta gamma) and put it in update
    // update = scale * div(flux) / vol
    // we want update to contain (dt/2) div (delta gamma)
    flux_div ( BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_ANYD(update),
               BL_TO_FORTRAN_ANYD(efab[0]),
               BL_TO_FORTRAN_ANYD(efab[1]),
#if ( AMREX_SPACEDIM == 3 )
               BL_TO_FORTRAN_ANYD(efab[2]),
#endif
               BL_TO_FORTRAN_ANYD(volume[mfi]),
               &nspecies,&scale);

    // add RHS from diffusion solve
    update.plus(Rhs[iGrid],box,0,0,nspecies);

    // Ssync = "RHS from diffusion solve" + (dt/2) * div (delta gamma)
    Ssync[mfi].copy(update,box,0,box,first_spec-BL_SPACEDIM,nspecies);
  }
}
  Rhs.clear();
  //
  // Do refluxing AFTER flux adjustment
  //
  if (do_reflux && level > 0 && last_mac_sync_iter)
  {
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      getViscFluxReg().FineAdd(*SpecDiffusionFluxnp1[d],d,0,first_spec,nspecies,dt);
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

  BL_PROFILE_REGION_STOP("R::HT::differential_spec_diffuse_sync()");
}

void
PeleLM::reflux ()
{
  BL_PROFILE("HT::reflux()");
  // no need to reflux if this is the finest level
  if (level == parent->finestLevel()) return;

  const Real strt_time = ParallelDescriptor::second();

  BL_ASSERT(do_reflux);

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
  fr_visc.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
  if (do_reflux_visc)
    fr_visc.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

  showMF("sdcSync",Ssync,"sdc_Ssync_after_viscReflux",level);

  const MultiFab& RhoHalftime = get_rho_half_time();

  if (do_mom_diff == 0) 
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Vsync,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();

      D_TERM(Vsync[mfi].divide(RhoHalftime[mfi],box,0,Xvel,1);,
             Vsync[mfi].divide(RhoHalftime[mfi],box,0,Yvel,1);,
             Vsync[mfi].divide(RhoHalftime[mfi],box,0,Zvel,1););
    }
  }
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    FArrayBox tmp;
  
    // for any variables that used non-conservative advective differencing,
    // divide the sync by rhohalf
    for (MFIter mfi(Ssync,true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      tmp.resize(bx,1);
      tmp.copy(RhoHalftime[mfi],0,0,1);
      tmp.invert(1);
  
      for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
      {
        if (advectionType[istate] == NonConservative)
        {
          const int sigma = istate -  BL_SPACEDIM;
  
          Ssync[mfi].mult(tmp,bx,0,sigma,1);
        }
      }
    }
  
    tmp.clear();
  }
  // take divergence of advective flux registers into cell-centered RHS
  fr_adv.Reflux(Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
  fr_adv.Reflux(Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);
  showMF("sdcSync",Ssync,"sdc_Ssync_after_advReflux",level);

  BoxArray baf = getLevel(level+1).boxArray();

  baf.coarsen(fine_ratio);
  //
  // This is necessary in order to zero out the contribution to any
  // coarse grid cells which underlie fine grid cells.
  //
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector< std::pair<int,Box> > isects;
  
    for (MFIter mfi(Vsync,true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.growntilebox(); 
      baf.intersections(box,isects);

      for (int i = 0, N = isects.size(); i < N; i++)
      {
        Vsync[mfi].setVal(0,isects[i].second,0,BL_SPACEDIM);
        Ssync[mfi].setVal(0,isects[i].second,0,NUM_STATE-BL_SPACEDIM);
      }
    }
  }
  showMF("sdcSync",Ssync,"sdc_Ssync_after_zero",level);

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
  const TimeLevel whichTime = which_time(State_Type, time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  compute_vel_visc(time, whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
}

void
PeleLM::calcDiffusivity (const Real time)
{
  BL_PROFILE("HT::calcDiffusivity()");

  const TimeLevel whichTime = which_time(State_Type, time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  const int  nGrow           = 1;
  const int  offset          = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
  MultiFab&  diff            = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
  

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

  FillPatchIterator Rho_and_spec_fpi(*this,diff,nGrow,time,State_Type,Density,nspecies+1),
         Temp_fpi(*this,diff,nGrow,time,State_Type,Temp,1);
  
  MultiFab& Rho_and_spec_mf = Rho_and_spec_fpi.get_mf();
  MultiFab& Temp_mf = Temp_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  FArrayBox tmp, bcen;
  for (MFIter mfi(Rho_and_spec_mf,true); mfi.isValid();++mfi)
  {
    FArrayBox& Tfab = Temp_mf[mfi];
    FArrayBox& RYfab = Rho_and_spec_mf[mfi];
    const Box& gbox = mfi.growntilebox();

    const int  vflag   = false;
    // rhoD + lambda + mu
    const int nc_bcen = nspecies+2; 
    int       dotemp  = 1;
    bcen.resize(gbox,nc_bcen);
    
    spec_temp_visc(BL_TO_FORTRAN_BOX(gbox),
                   BL_TO_FORTRAN_N_ANYD(Temp_mf[mfi],0),
                   BL_TO_FORTRAN_N_ANYD(Rho_and_spec_mf[mfi],1),
                   BL_TO_FORTRAN_N_ANYD(bcen,0),
                   &nc_bcen, &P1atm_MKS, &dotemp, &vflag, &p_amb);    
        
    FArrayBox& Dfab = diff[mfi];

    // beta for Y's
    Dfab.copy(bcen,gbox,0,gbox,first_spec-offset,nspecies);
    // lambda in Temp slot
    Dfab.copy(bcen,gbox,nspecies,gbox,Temp-offset,1);

    //
    // Convert from tmp=RhoY_l to Y_l
    //
    tmp.resize(gbox,1);
    tmp.copy(RYfab,gbox,0,gbox,0,1);
    tmp.invert(1,gbox);
    for (int n = 1; n < nspecies+1; n++)
    {
      RYfab.mult(tmp,gbox,0,n,1);
    }

    for (int icomp = RhoH; icomp <= NUM_STATE; icomp++)
    {
      if (icomp == RhoH)
      {
        // lambda/cp in RhoH slot
        const int sCompT = 0, sCompY = 1, sCompCp = RhoH-offset;
        getCpmixGivenTY_pphys(Dfab,Tfab,RYfab,gbox,sCompT,sCompY,sCompCp);
        Dfab.invert(1,gbox,sCompCp,1);
        Dfab.mult(Dfab,gbox,Temp-offset,RhoH-offset,1);
      }
      else if (icomp == Trac || icomp == RhoRT)
      {
        // fill Trac and RhoRT slot so trac_diff_coef (typically zero)
        Dfab.setVal(trac_diff_coef, gbox, icomp-offset, 1);
      }
    }
  }
}
  showMFsub("1D",diff,stripBox,"1D_calcD_visc",level);
}

#ifdef USE_WBAR
void
PeleLM::calcDiffusivity_Wbar (const Real time)
{
  BL_PROFILE("HT::calcDiffusivity_Wbar()");
  // diffn_cc or diffnp1_cc contains cell-centered transport coefficients from Y's
  //

  const TimeLevel whichTime = which_time(State_Type, time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab& diff       = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
  const int nGrow      = diff.nGrow();

  BL_ASSERT(diffWbar_cc.nGrow() >= nGrow);

  FillPatchIterator Rho_and_spec_fpi(*this,diff,nGrow,time,State_Type,Density,nspecies+1);
  MultiFab& mf=Rho_and_spec_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif  
  for (MFIter mfi(mf,true); mfi.isValid();++mfi)
  {
    const FArrayBox& RD = diff[mfi];
    const FArrayBox& RYfab = mf[mfi];
    FArrayBox& Dfab_Wbar = diffWbar_cc[mfi];
    const Box& gbox = mfi.growntilebox();
        
    beta_wbar(BL_TO_FORTRAN_BOX(gbox),
              BL_TO_FORTRAN_ANYD(RD),
              BL_TO_FORTRAN_ANYD(Dfab_Wbar),
              BL_TO_FORTRAN_N_ANYD(RYfab,1));
  }
}
#endif

void
PeleLM::getViscosity (MultiFab*  beta[BL_SPACEDIM],
                      const Real time)
{
  BL_PROFILE("HT::getViscosity()");
  const TimeLevel whichTime = which_time(State_Type, time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* visc = (whichTime == AmrOldTime) ? viscn_cc : viscnp1_cc;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter viscMfi(*visc,true); viscMfi.isValid(); ++viscMfi)
  {
    const Box& box = viscMfi.tilebox();

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).lo(dir));
      FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).hi(dir));
      const Box& ebox = viscMfi.nodaltilebox(dir);
      center_to_edge_fancy((*visc)[viscMfi],(*beta[dir])[viscMfi],
                           amrex::grow(box,amrex::BASISV(dir)), ebox, 0, 0, 1,
                           geom.Domain(), bc_lo, bc_hi);
    }
  }
}

void
PeleLM::getDiffusivity (MultiFab*  beta[BL_SPACEDIM],
                        const Real time,
                        const int  state_comp,
                        const int  dst_comp,
                        const int  ncomp)
{
  BL_PROFILE("HT::getDiffusivity()");
  BL_ASSERT(state_comp > Density);

  const TimeLevel whichTime = which_time(State_Type, time);

  BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

  MultiFab* diff      = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;
  const int offset    = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
  int       diff_comp = state_comp - offset;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter diffMfi(*diff,true); diffMfi.isValid(); ++diffMfi)
  {
    const Box& box = diffMfi.tilebox();

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).lo(dir));
      FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).hi(dir));
      const Box& ebox = diffMfi.nodaltilebox(dir);
      center_to_edge_fancy((*diff)[diffMfi],(*beta[dir])[diffMfi],
                           amrex::grow(box,amrex::BASISV(dir)), ebox, diff_comp, 
                           dst_comp, ncomp, geom.Domain(), bc_lo, bc_hi);
    }
  }

  if (zeroBndryVisc > 0)
    zeroBoundaryVisc(beta,time,state_comp,dst_comp,ncomp);
}

#ifdef USE_WBAR
void
PeleLM::getDiffusivity_Wbar (MultiFab*  beta[BL_SPACEDIM],
                             const Real time)	   
{
  BL_PROFILE("HT::getDiffusivity_Wbar()");
  MultiFab& diff = diffWbar_cc;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter diffMfi(diff,true); diffMfi.isValid(); ++diffMfi)
  {
    const Box& box = diffMfi.tilebox();

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(first_spec).lo(dir));
      FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(first_spec).hi(dir));
      const Box& ebox = diffMfi.nodaltilebox(dir);
      center_to_edge_fancy(diff[diffMfi],(*beta[dir])[diffMfi],
                           amrex::grow(box,amrex::BASISV(dir)), ebox, 0,
                           0, nspecies, geom.Domain(), bc_lo, bc_hi);
    }
  }

  if (zeroBndryVisc > 0)
    zeroBoundaryVisc(beta,time,BL_SPACEDIM+1,0,nspecies);
}
#endif

void
PeleLM::zeroBoundaryVisc (MultiFab*  beta[BL_SPACEDIM],
                          const Real time,
                          const int  state_comp,
                          const int  dst_comp,
                          const int  ncomp) const
{
  BL_ASSERT(state_comp > Density);

  const int isrz = (int) geom.IsRZ();
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      Box edom = amrex::surroundingNodes(geom.Domain(),dir);
  
      for (MFIter mfi(*(beta[dir]),true); mfi.isValid(); ++mfi)
      {
        FArrayBox& beta_fab = (*(beta[dir]))[mfi];
        const Box& ebox     = amrex::surroundingNodes(mfi.growntilebox(),dir);
        zero_visc(BL_TO_FORTRAN_N_ANYD(beta_fab,dst_comp),
                  BL_TO_FORTRAN_BOX(ebox),
                  BL_TO_FORTRAN_BOX(edom),
                  geom.CellSize(), geom.ProbLo(), phys_bc.vect(),
                  &dir, &isrz, &state_comp, &ncomp);
      }
    }
  }
}

void
PeleLM::compute_vel_visc (Real      time,
                          MultiFab* beta)
{
  const int nGrow = beta->nGrow();

  BL_ASSERT(nGrow == 1);

  MultiFab dummy(grids,dmap,1,0,MFInfo().SetAlloc(false));

  FillPatchIterator Temp_fpi(*this,dummy,nGrow,time,State_Type,Temp,1),
                    Spec_fpi(*this,dummy,nGrow,time,State_Type,first_spec,nspecies);
  MultiFab& Temp_mf = Temp_fpi.get_mf();
  MultiFab& Spec_mf = Spec_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    for (MFIter mfi(Temp_mf,true); mfi.isValid();++mfi)
    {
      const Box& box  = mfi.growntilebox();
      FArrayBox& temp = Temp_mf[mfi];
      FArrayBox& spec = Spec_mf[mfi];
      
      vel_visc(BL_TO_FORTRAN_BOX(box),
               BL_TO_FORTRAN_N_ANYD(temp,0),
               BL_TO_FORTRAN_N_ANYD(spec,0),
               BL_TO_FORTRAN_N_ANYD((*beta)[mfi],0));
    }
  }
}

void
PeleLM::calc_divu (Real      time,
                   Real      dt,
                   MultiFab& divu)
{
  BL_PROFILE("HT::calc_divu()");
  const int nGrow = 0;

  int       vtCompY=0, vtCompT=0;
  MultiFab  mcViscTerms;

#ifdef USE_WBAR
  MultiFab DWbar_temp(grids,dmap,nspecies,nGrowAdvForcing);
#endif

  vtCompT = nspecies + 1;
  vtCompY = 0;
  mcViscTerms.define(grids,dmap,nspecies+2,nGrow);

  // we don't want to update flux registers due to fluxes in divu computation
  bool do_reflux_hold = do_reflux;
  do_reflux = false;
    
  // DD is computed and stored in divu, but we don't need it and overwrite
  // divu in CALCDIVU.
#ifdef USE_WBAR
  compute_differential_diffusion_terms(mcViscTerms,divu,DWbar_temp,time,dt);
#else
  compute_differential_diffusion_terms(mcViscTerms,divu,time,dt);
#endif

  do_reflux = do_reflux_hold;

  //
  // if we are in the initial projection (time=0, dt=-1), set RhoYdot=0
  // if we are in a divu_iter (time=0, dt>0), use I_R
  // if we are in an init_iter or regular time step (time=dt, dt>0), use instananeous
  MultiFab   RhoYdotTmp;
  MultiFab&  S       = get_data(State_Type,time);
  const bool use_IR  = (time == 0 && dt > 0);

  MultiFab&  RhoYdot = (use_IR) ? get_new_data(RhoYdot_Type) : RhoYdotTmp;

  if (!use_IR)
  {
    if (time == 0)
    {
      // initial projection, set omegadot to zero
      RhoYdot.define(grids,dmap,nspecies,0);
      RhoYdot.setVal(0);
    }
    else if (dt > 0)
    {
      // init_iter or regular time step, use instantaneous omegadot
      RhoYdot.define(grids,dmap,nspecies,0);
      compute_instantaneous_reaction_rates(RhoYdot,S,time,nGrow);
    }
    else
    {
      amrex::Abort("bad divu_logic - shouldn't be here");
    }
  }
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& box = mfi.tilebox();            
    FArrayBox& du = divu[mfi];
    const FArrayBox& vtY = mcViscTerms[mfi];
    const FArrayBox& vtT = mcViscTerms[mfi];
    const FArrayBox& rhoY = S[mfi];
    const FArrayBox& rYdot = RhoYdot[mfi];
    const FArrayBox& T = S[mfi];
    calc_divu_fortran(BL_TO_FORTRAN_BOX(box),
                      BL_TO_FORTRAN_ANYD(du),
                      BL_TO_FORTRAN_ANYD(rYdot),
                      BL_TO_FORTRAN_N_ANYD(vtY, vtCompY),
                      BL_TO_FORTRAN_N_ANYD(vtT, vtCompT),
                      BL_TO_FORTRAN_N_ANYD(rhoY, first_spec),
                      BL_TO_FORTRAN_N_ANYD(T,Temp));
  }
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
  BL_PROFILE("HT::calc_dpdt()");

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

  const int pComp = (have_rhort ? RhoRT : Trac);
  int nGrow = dpdt.nGrow();
  MultiFab Peos(grids,dmap,1,nGrow);
  FillPatchIterator S_fpi(*this,Peos,nGrow,time,State_Type,pComp,1);
  MultiFab& Smf=S_fpi.get_mf();
  
#ifdef _OPENMP
#pragma omp parallel
#endif  
  for (MFIter mfi(Smf,true); mfi.isValid();++mfi)
  {
    const Box& bx = mfi.tilebox();
    FArrayBox&       S   = Smf[mfi];
    Peos[mfi].copy(S,bx,0,bx,0,1);
  }
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(dpdt,true); mfi.isValid(); ++mfi)
  {
    const Box& vbox = mfi.tilebox();

    dpdt[mfi].copy(Peos[mfi],vbox,0,vbox,0,1);
    dpdt[mfi].plus(-p_amb,vbox);
    dpdt[mfi].mult(1.0/dt,vbox,0,1);
    dpdt[mfi].divide(Peos[mfi],vbox,0,0,1);
    dpdt[mfi].mult(dpdt_factor,vbox,0,1);
  }

  Peos.clear();

  if (nGrow > 0) {
    dpdt.FillBoundary(0,1, geom.periodicity());
    BL_ASSERT(dpdt.nGrow() == 1);
    Extrapolater::FirstOrderExtrap(dpdt, geom, 0, 1);
  }
}

//
// Function to use if Divu_Type and Dsdt_Type are in the state.
//

void
PeleLM::calc_dsdt (Real      time,
                   Real      dt,
                   MultiFab& dsdt)
{
  dsdt.setVal(0);
}


void
PeleLM::RhoH_to_Temp (MultiFab& S,
                      int       nGrow,
                      int       dominmax)
{
  BL_PROFILE("HT::RhoH_to_Temp()");
  const Real strt_time = ParallelDescriptor::second();
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

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(max:max_iters)
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
  {
    const Box& box = mfi.growntilebox(nGrow);
    max_iters = std::max(max_iters, RhoH_to_Temp(S[mfi],box,dominmax));
  }

  if (verbose > 1)
  {
    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_time = ParallelDescriptor::second() - strt_time;

    ParallelDescriptor::ReduceIntMax(max_iters,IOProc);
    ParallelDescriptor::ReduceRealMax(run_time,IOProc);

    if (verbose) amrex::Print() << "PeleLM::RhoH_to_Temp: max_iters = " << max_iters << ", time: " << run_time << '\n';
  }
  //
  // Reset it back.
  //
  htt_hmixTYP = htt_hmixTYP_SAVE;
}

int PeleLM::RhoH_to_Temp_DoIt(FArrayBox&       Tfab,
                      const FArrayBox& Hfab,
                      const FArrayBox& Yfab,
                      const Box&       box,
                      int              sCompH,
                      int              sCompY,
                      int              dCompT,
                      const Real&      htt_hmixTYP)
{
	// FIXME
  const Real eps = 1.0e-8; //getHtoTerrMAX_pphys();
  Real errMAX = eps*htt_hmixTYP;

  int iters = getTGivenHY_pphys(Tfab,Hfab,Yfab,box,sCompH,sCompY,dCompT,errMAX);

  if (iters < 0)
    amrex::Error("PeleLM::RhoH_to_Temp(fab): error in H->T");

  return iters;
}


int
PeleLM::RhoH_to_Temp (FArrayBox& S,
                      const Box& box,
                      int        dominmax)
{
  BL_ASSERT(S.box().contains(box));

  // temporary storage for 1/rho, H, and Y
  FArrayBox rhoInv, H, Y;
  rhoInv.resize(box,1);
  H.resize(box,1);
  Y.resize(box,nspecies);

  // copy in rho, rhoH, and rhoY
  rhoInv.copy(S,box,Density,   box,0,1);
  H.copy     (S,box,RhoH,      box,0,1);
  Y.copy     (S,box,first_spec,box,0,nspecies);

  // invert rho, then multiply 1/rho by rhoH and rhoY to get H and Y
  rhoInv.invert(1,box,0,1);
  H.mult(rhoInv,box,0,0,1);
  for (int i=0; i<nspecies; i++)
  {
    Y.mult(rhoInv,box,0,i,1);
  }    

  // we index into Temperature component in S.  H and Y begin with the 0th component
  int iters = RhoH_to_Temp_DoIt(S,H,Y,box,0,0,Temp,htt_hmixTYP);

  if (dominmax)
    FabMinMax(S, box, htt_tempmin, htt_tempmax, Temp, 1);    

  return iters;
}

void
PeleLM::setPlotVariables ()
{
  AmrLevel::setPlotVariables();

  Vector<std::string> names;
  PeleLM::getSpeciesNames(names);

  
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

    for (int n = 1; n <= buildInfoGetNumModules(); n++) {
      jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
    }

    jobInfoFile << "\n";

    const char* githash1 = buildInfoGetGitHash(1);
    const char* githash2 = buildInfoGetGitHash(2);
    const char* githash3 = buildInfoGetGitHash(3);
    if (strlen(githash1) > 0) {
      jobInfoFile << "PeleLM git hash: " << githash1 << "\n";
    }
    if (strlen(githash2) > 0) {
      jobInfoFile << "AMReX git hash: " << githash2 << "\n";
    }
    if (strlen(githash3) > 0) {
      jobInfoFile << "IAMR   git hash: " << githash3 << "\n";
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

  std::string Level = amrex::Concatenate("Level_", level, 1);
  //
  // Now for the full pathname of that directory.
  //
  std::string FullPath = dir;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    FullPath += '/';
  FullPath += Level;
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
      for (int n = 0; n < BL_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0)
    {
      std::string PathNameInHeader = Level;
      PathNameInHeader += BaseName;
      os << PathNameInHeader << '\n';
    }
  }
  //
  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  // NOTE: we are assuming that each state variable has one component,
  // but a derived variable is allowed to have multiple components.
  int       cnt   = 0;
  int       ncomp = 1;
  const int nGrow = 0;
  MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
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
  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  VisMF::Write(plotMF,TheFullPath,how);
}

std::unique_ptr<MultiFab>
PeleLM::derive (const std::string& name,
                Real               time,
                int                ngrow)
{        
  BL_ASSERT(ngrow >= 0);
  
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
  if (name == "mean_progress_curvature")
  {
    // 
    // Smooth the temperature, then use smoothed T to compute mean progress curvature
    //

    // Assert because we do not know how to un-convert the destination
    //   and also, implicitly assume the convert that was used to build mf BA is trivial
    BL_ASSERT(mf.boxArray()[0].ixType()==IndexType::TheCellType());

    const DeriveRec* rec = derive_lst.get(name);

    Box srcB = mf.boxArray()[0];
    Box dstB = rec->boxMap()(srcB);

    // Find nGrowSRC
    int nGrowSRC = 0;
    for ( ; !srcB.contains(dstB); ++nGrowSRC)
    {
      srcB.grow(1);
    }
    BL_ASSERT(nGrowSRC);  // Need grow cells for this to work!

    FillPatchIterator fpi(*this,mf,nGrowSRC,time,State_Type,Temp,1);
    MultiFab& tmf = fpi.get_mf();

    int num_smooth_pre = 3;
        
    for (int i=0; i<num_smooth_pre; ++i)
    {
      // Fix up fine-fine and periodic
      if (i!=0) 
        tmf.FillBoundary(0,1,geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(tmf,true); mfi.isValid(); ++mfi)
      {
        // 
        // Use result MultiFab for temporary container to hold smooth T field
        // 
        const Box& box = mfi.tilebox();
        smooth( BL_TO_FORTRAN_BOX(box),
                BL_TO_FORTRAN_ANYD(tmf[mfi]),
                BL_TO_FORTRAN_N_ANYD(mf[mfi],dcomp));
      }
      
      MultiFab::Copy(tmf,mf,0,0,1,0);
    }

    tmf.FillBoundary(0,1,geom.periodicity());

    const Real* dx = geom.CellSize();
        
    FArrayBox nWork;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
      const FArrayBox& Tg = tmf[mfi];
      FArrayBox& MC = mf[mfi];
      const Box& box = mfi.tilebox();
      const Box& nodebox = amrex::surroundingNodes(box);
      nWork.resize(nodebox,BL_SPACEDIM);
            
      mcurve( BL_TO_FORTRAN_BOX(box),
              BL_TO_FORTRAN_ANYD(Tg),
              BL_TO_FORTRAN_N_ANYD(MC,dcomp),
              BL_TO_FORTRAN_ANYD(nWork),
              dx);
    }
  } 
  else
  {
#ifdef AMREX_PARTICLES
    ParticleDerive(name,time,mf,dcomp);
#else
    AmrLevel::derive(name,time,mf,dcomp);
#endif
  }
}

void
PeleLM::errorEst (TagBoxArray& tags,
                  int         clearval,
                  int         tagval,
                  Real        time,
                  int         n_error_buf, 
                  int         ngrow)
{
  BL_PROFILE("HT::errorEst()");
  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  for (int j = 0; j < err_list.size(); j++)
  {
    const ErrorRec::ErrorFunc& efunc = err_list[j].errFunc();
    const LM_Error_Value* lmfunc = dynamic_cast<const LM_Error_Value*>(&efunc);
    bool box_tag = lmfunc && lmfunc->BoxTag();

    auto mf = box_tag ? 0 : derive(err_list[j].name(), time, err_list[j].nGrow());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grids,dmap,true); mfi.isValid(); ++mfi)
    {
      const Box&  vbx     = mfi.tilebox();
      RealBox     gridloc = RealBox(vbx,geom.CellSize(),geom.ProbLo());
      Vector<int> itags   = tags[mfi].tags();
      int*        tptr    = itags.dataPtr();
      const int*  tlo     = tags[mfi].box().loVect();
      const int*  thi     = tags[mfi].box().hiVect();
      const Real* xlo     = gridloc.lo();

      if (box_tag)
      {
        lmfunc->tagCells1(tptr, tlo, thi,
                          &tagval, &clearval,
                          BL_TO_FORTRAN_BOX(vbx),
                          BL_TO_FORTRAN_BOX(geom.Domain()),
                          dx, xlo, prob_lo, &time, &level);
      }
      else
      {
        FArrayBox&  fab     = (*mf)[mfi];
        Real*       dat     = fab.dataPtr();
        const int*  dlo     = fab.box().loVect();
        const int*  dhi     = fab.box().hiVect();
        const int   ncomp   = fab.nComp();
        
        if (lmfunc==0) 
        {
          err_list[j].errFunc()(tptr, tlo, thi,
                                &tagval, &clearval,
                                BL_TO_FORTRAN_ANYD(fab),
                                BL_TO_FORTRAN_BOX(vbx), &ncomp,
                                BL_TO_FORTRAN_BOX(geom.Domain()),
                                dx, xlo, prob_lo, &time, &level);
        }
        else
        {
          lmfunc->tagCells(tptr, tlo, thi,
                           &tagval, &clearval,
                           dat, dlo, dhi, 
                           BL_TO_FORTRAN_BOX(vbx), &ncomp,
                           BL_TO_FORTRAN_BOX(geom.Domain()),
                           dx, xlo, prob_lo, &time, &level);
        }
      }
                      
      tags[mfi].tags(itags);

    }
  }
}
