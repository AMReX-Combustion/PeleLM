
#ifndef _PeleLM_H_
#define _PeleLM_H_

#include <AMReX_AuxBoundaryData.H>
#include <NavierStokesBase.H>

#include <AMReX_ForkJoin.H>

#include <mechanism.h>
#include <pelelm_prob_parm.H>
#include <PeleLM_parm.H>

#include <list>
#include <map>
#include <utility>

//
// Note: define TEMPERATURE if you want a variable T in the 
//       State_Type part of the state
//       whose component index is Temp and 
//       whose evoution equation is 
//       \pd (rho T)/\pd t + diver (\rho U T) = \diver k/c_p grad T
//       define RADIATION only if TEMPERATURE is also defined and the evolution equation
//       for T is
//       \pd (rho T)/\pd t + diver (\rho U T) = \diver k/c_p grad T - 1/c_p diver q_rad
//
//       Note that component Temp is T, not rho T. This was done so
//       that we could use the existing diffusion operator and 
//       multigrid code.
//

extern "C"
{
  typedef void (*LMEF)(BL_FORT_IFAB_ARG_ANYD(tag),
                       const int* tagval, const int* clearval,
                       const BL_FORT_FAB_ARG_ANYD(data),
                       const int* lo, const int* hi, const int* nvar,
                       const int* domain_lo, const int* domain_hi,
                       const amrex::Real* dx, const amrex::Real* xlo,
                       const amrex::Real* prob_lo, const amrex::Real* time,
                       const int* level, const amrex::Real* value);

  typedef void (*LMEF_BOX) (BL_FORT_IFAB_ARG_ANYD(tag),
                            const int* tagval, const int* clearval,
                            const amrex::Real* blo, const amrex::Real* bhi,
                            const int* lo, const int* hi,
                            const int* domain_lo, const int* domain_hi,
                            const amrex::Real* dx, const amrex::Real* xlo,
                            const amrex::Real* prob_lo, const amrex::Real* time,
                            const int* level);
}

class LM_Error_Value
  :
  public amrex::ErrorRec::ErrorFunc
{
public:
  LM_Error_Value()
    :
    lmef(0), lmef_box(0), value(), min_time(), max_time(), max_level() {}
  
  LM_Error_Value (amrex::Real min_time, amrex::Real max_time, int max_level);
  
  LM_Error_Value (LMEF lmef, amrex::Real value, amrex::Real min_time,
                  amrex::Real max_time, int max_level);
  
  LM_Error_Value (LMEF_BOX lmef_box, const amrex::RealBox& box, amrex::Real min_time,
                  amrex::Real max_time, int max_level);
  
  virtual ~LM_Error_Value () {}
  
  virtual amrex::ErrorRec::ErrorFunc* clone () const {
    if (BoxTag()) {
      return new LM_Error_Value(lmef_box,box,min_time,max_time,max_level);
    }
    return new LM_Error_Value(lmef,value,min_time,max_time,max_level);
  }
  
  void tagCells(int* tag, const int* tlo, const int* thi,
                const int* tagval, const int* clearval,
                const amrex::Real* data, const int* dlo, const int* dhi,
                const int* lo, const int* hi, const int* nvar,
                const int* domain_lo, const int* domain_hi,
                const amrex::Real* dx, const amrex::Real* xlo,
                const amrex::Real* prob_lo, const amrex::Real* time,
                const int* level) const;

  void tagCells1(int* tag, const int* tlo, const int* thi,
                 const int* tagval, const int* clearval,
                 const int* lo, const int* hi,
                 const int* domain_lo, const int* domain_hi,
                 const amrex::Real* dx, const amrex::Real* xlo,
                 const amrex::Real* prob_lo, const amrex::Real* time,
                 const int* level) const;
  
  int MaxLevel() const {return max_level;}
  amrex::Real MinTime() const {return min_time;}
  amrex::Real MaxTime() const {return max_time;}
  amrex::Real Value() const {return value;}
  bool BoxTag() const {return lmef==0 && lmef_box!=0;}
  
protected:
  LMEF lmef;
  LMEF_BOX lmef_box;
  amrex::Real value, min_time, max_time;
  amrex::RealBox box;
  int max_level;
};

class PeleLM
  :
  public NavierStokesBase
{

public:

  PeleLM ();
  PeleLM (amrex::Amr&            papa,
          int                    lev,
          const amrex::Geometry& level_geom,
          const amrex::BoxArray& bl,
          const amrex::DistributionMapping& dm,
          amrex::Real            time);

  virtual ~PeleLM ();

  ////////////////////////////////////////////////////////////////////////////
  //    AmrLevel virtual functions                                          //
  ////////////////////////////////////////////////////////////////////////////

  //
  // Advance grids at this level in time.
  //
  virtual amrex::Real advance (amrex::Real time,
                               amrex::Real dt,
                               int  iteration,
                               int  ncycle) override;

  virtual void checkPoint (const std::string& dir,
                           std::ostream&      os,
                           amrex::VisMF::How  how = amrex::VisMF::OneFilePerCPU,
                           bool               dump_old = true) override;
  //
  // Returns a MultiFab containing the derived data for this level.
  // The user is responsible for deleting this pointer when done
  // with it.  If ngrow>0 the MultiFab is built on the appropriately
  // grown amrex::BoxArray.
  //
  virtual std::unique_ptr<amrex::MultiFab> derive (const std::string& name,
                                                   amrex::Real        time,
                                                   int                ngrow) override;
  //
  // This version of derive() fills the dcomp'th component of mf with the derived quantity.
  //
  virtual void derive (const std::string& name,
                       amrex::Real        time,
                       amrex::MultiFab&   mf,
                       int                dcomp) override;
  //
  // Init data on this level after regridding if old level
  // did not exist previously.
  //
  virtual void init () override;
  //
  // Init data on this level from another NavierStokes (during regrid).
  //
  virtual void init (amrex::AmrLevel& old) override;

  virtual void initData () override;
  //
  // Contains operations to be done only after a full coarse timestep.
  //
  virtual void postCoarseTimeStep (amrex::Real cumtime) override;

  virtual void post_init (amrex::Real stop_time) override;

  virtual void post_regrid (int lbase, int new_finest) override;

  virtual void post_restart () override;

  virtual void post_timestep (int iteration) override;

  virtual void restart (amrex::Amr&   papa,
                        std::istream& is,
                        bool          bReadSpecial = false) override;
  //
  // Set the variables that are put in the plotfile....
  //
  virtual void setPlotVariables () override;
  //
  // Set time levels of state data.
  //
  virtual void setTimeLevel (amrex::Real time,
                             amrex::Real dt_old,
                             amrex::Real dt_new) override;
  //
  // Write plot file stuff to specified directory.
  //
  virtual void writePlotFile (const std::string& dir,
                              std::ostream&      os,
                              amrex::VisMF::How  how) override;

  ////////////////////////////////////////////////////////////////////////////
  //    PeleLM public static functions                                //
  ////////////////////////////////////////////////////////////////////////////

  //
  // Define data descriptors.
  //
  static void variableSetUp ();

  static void rhoydotSetUp ();
  //
  // Cleanup data descriptors at end of run.
  //
  static void variableCleanUp ();

  static void getSpeciesNames(amrex::Vector<std::string>& spn);

  static int getSpeciesIdx(const std::string& spName);

  void advance_chemistry (amrex::MultiFab&       mf_old,
                          amrex::MultiFab&       mf_new,
                          amrex::Real            dt,
                          const amrex::MultiFab& Force,
                          int             nCompF,
                          bool            use_stiff_solver = true);

  void compute_scalar_advection_fluxes_and_divergence (const amrex::MultiFab& Force,
                                                       const amrex::MultiFab& divu,
                                                       amrex::Real            dt);

  // Mixture fraction machinery
  static amrex::Array<amrex::Real, 4> Beta_mix;
  static amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;
  static amrex::Real Zfu;
  static amrex::Real Zox;
  static bool        mixture_fraction_ready;

  ////////////////////////////////////////////////////////////////////////////
  //    Overriding Virtual Functions in NavierStokesBase                    //
  ////////////////////////////////////////////////////////////////////////////

  //
  // Setup for a level timestep.
  //
  virtual void advance_setup (amrex::Real time,
                              amrex::Real dt,
                              int  iteration,
                              int  ncycle) override;

  virtual void avgDown () override; // Average down for all the state types.
  //
  // Note: these two functions must be supplied.
  //       If divu is not included in the state, then
  //       they can be no-op functions
  //
  virtual void calc_divu (amrex::Real      time,
                          amrex::Real      dt,
                          amrex::MultiFab& divu) override;

  virtual void calcViscosity (const amrex::Real time,
                              const amrex::Real dt,
                              const int  iteration,
                              const int  ncycle) override;

  virtual void calcDiffusivity (const amrex::Real time) override;

  virtual amrex::Real estTimeStep () override;

  virtual void getViscosity (amrex::MultiFab*  viscosity[BL_SPACEDIM],
                             const amrex::Real time) override;

  virtual void getViscTerms (amrex::MultiFab& visc_terms,
                             int       src_comp, 
                             int       num_comp,
                             amrex::Real      time) override;

  virtual void getForce (amrex::FArrayBox&       force,
                         const amrex::Box&       bx,
                         int                     strt_comp,
                         int                     num_comp,
                         const amrex::Real       time,
                         const amrex::FArrayBox& Vel,
                         const amrex::FArrayBox& Scal,
                         int                     scalScomp) override;


  virtual void mac_sync () override;
  //
  // Crse/fine fixup functions.
  //
  virtual void reflux () override;
  //
  // Reset time levels for the initial iterations.
  //
  virtual void resetState (amrex::Real time,
                           amrex::Real dt_old,
                           amrex::Real dt_new) override;

  virtual void sum_integrated_quantities () override;

  virtual void scalar_advection_update (amrex::Real dt,
                                        int  first_scalar,
                                        int  last_scalar) override;

  virtual void velocity_diffusion_update (amrex::Real dt) override;

  virtual void errorEst (amrex::TagBoxArray& tags,
                         int                 clearval,
                         int                 tagval,
                         amrex::Real         time,
                         int                 n_error_buf, 
                         int                 ngrow) override;

  ////////////////////////////////////////////////////////////////////////////
  //    PeleLM protected static functions                             //
  ////////////////////////////////////////////////////////////////////////////

  static void Initialize ();
  static void Initialize_specific ();
  static void init_mixture_fraction ();
  static void Finalize ();
    
#ifdef AMREX_PARTICLES
  virtual int timestamp_num_extras () override;
  virtual void timestamp_add_extras (int lev, amrex::Real time, amrex::MultiFab& mf) override;
  static void read_particle_params ();
#endif


  // enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };
  // enum Solver_Status {HT_InProgress, HT_Stalled, HT_Solved};

  void adjust_spec_diffusion_fluxes (amrex::MultiFab* const * flux,
                                     const amrex::MultiFab&   S,
                                     const amrex::BCRec&      bc,
                                     amrex::Real              time);

  void calcDiffusivity_Wbar (const amrex::Real time);

  void calc_dpdt (amrex::Real      time,
                  amrex::Real      dt,
                  amrex::MultiFab& dpdt,
                  amrex::MultiFab* u_mac);

  void checkTimeStep (amrex::Real dt);

  void compute_differential_diffusion_fluxes (const amrex::MultiFab& S,
                                              const amrex::MultiFab* Scrse,
                                              amrex::MultiFab* const * flux,
                                              const amrex::MultiFab* const * beta,
                                              amrex::Real dt,
					      amrex::Real time,
                                              bool  include_Wbar_fluxes);

  void compute_differential_diffusion_terms (amrex::MultiFab& D,
                                             amrex::MultiFab& DD,
                                             amrex::Real      time,
                                             amrex::Real      dt,
                                             bool             include_Wbar_terms);

  void compute_enthalpy_fluxes (amrex::MultiFab* const*       flux,
                                const amrex::MultiFab* const* beta,
				amrex::Real		      time);

  enum HowToFillGrow {HT_ZERO_GROW_CELLS, HT_EXTRAP_GROW_CELLS, HT_NUM_GROW_OPTIONS};

  void compute_instantaneous_reaction_rates (amrex::MultiFab&       R,
                                             const amrex::MultiFab& S,
                                             amrex::Real            time,
                                             int             nGrow = 0,
                                             HowToFillGrow   how = HT_ZERO_GROW_CELLS);

  amrex::Real adjust_p_and_divu_for_closed_chamber(amrex::MultiFab& mac_divu);

  void compute_rhohmix (amrex::Real      time,
                        amrex::MultiFab& rhohmix,
                        int              dComp);

  void compute_rhoRT (const amrex::MultiFab& S,
                            amrex::MultiFab& P,
                                   int       pComp);

  void compute_Wbar_fluxes(amrex::Real time,
                           int         inc_flag,
                           amrex::Real inc_coeff);

  void define_data ();

  void differential_diffusion_update (amrex::MultiFab& Force,
                                      amrex::MultiFab& DWbar,
                                      int       FComp,
                                      amrex::MultiFab& D,
                                      int       DComp,
                                      amrex::MultiFab& DD);

  void differential_spec_diffuse_sync (amrex::Real dt,
                                       bool Wbar_corrector,
                                       bool last_mac_sync_iter);

  void diffuse_velocity_setup (amrex::Real        dt,
                               amrex::MultiFab*&  delta_rhs,
                               FluxBoxes&  fb_betan, 
                               FluxBoxes&  fb_betanp1);

  void flux_divergence (amrex::MultiFab&        fdiv,
                        int              fdivComp,
                        const amrex::MultiFab* const* extensive_fluxes,
                        int              fluxComp,
                        int              nComp,
                        amrex::Real             scale) const;

  void getDiffusivity (amrex::MultiFab* diffusivity[BL_SPACEDIM],
                       const amrex::Real time,
                       const int state_comp,
                       const int dst_comp,
                       const int num_comp);

  void getDiffusivity_Wbar (amrex::MultiFab* diffusivity[BL_SPACEDIM],
                            const amrex::Real time);

  amrex::DistributionMapping getFuncCountDM (const amrex::BoxArray& bxba, int ngrow);

  PeleLM& getLevel (int lev)
    {
      return *(PeleLM*) &parent->getLevel(lev);
    }

  void initDataOtherTypes ();

  void post_init_press (amrex::Real&        dt_init,
                        amrex::Vector<int>&  nc_save,
                        amrex::Vector<amrex::Real>& dt_save);

  void set_htt_hmixTYP ();

  void set_reasonable_grow_cells_for_R(amrex::Real time);

  void setThermoPress(amrex::Real time);

  void set_typical_values(bool restart);
  void reset_typical_values(const amrex::MultiFab& S);

  void activeControl(const int         step,
                     const int         restart,
                     const amrex::Real time,
                     const amrex::Real dt);

  void initActiveControl();

  static void parseComposition(amrex::Vector<std::string> compositionIn,
                               std::string                compositionType,
                               amrex::Real               *massFrac);

  void state_stats (amrex::MultiFab& S);

  void zeroBoundaryVisc (amrex::MultiFab*  beta[AMREX_SPACEDIM],
                         const amrex::Real time,
                         const int  state_comp,
                         const int  dst_comp,
                         const int  ncomp) const;

  void diffuse_scalar_fj  (const amrex::Vector<amrex::MultiFab*>& S_old,
                           const amrex::Vector<amrex::MultiFab*>& Rho_old,
                           amrex::Vector<amrex::MultiFab*>&       S_new,
                           const amrex::Vector<amrex::MultiFab*>& Rho_new,
                           int                                    S_comp,
                           int                                    num_comp,
                           int                                    Rho_comp,
                           amrex::Real                            prev_time,
                           amrex::Real                            curr_time,
                           amrex::Real                            be_cn_theta,
                           const amrex::MultiFab&                 rho_half,
                           int                                    rho_flag,
                           amrex::MultiFab* const*                fluxn,
                           amrex::MultiFab* const*                fluxnp1,
                           int                                    fluxComp,
                           amrex::MultiFab*                       delta_rhs, 
                           int                                    rhsComp,
                           const amrex::MultiFab*                 alpha, 
                           int                                    alphaComp,
                           const amrex::MultiFab* const*          betan, 
                           const amrex::MultiFab* const*          betanp1,
                           int                                    betaComp,
                           const amrex::Vector<amrex::Real>&      visc_coef,
                           int                                    visc_coef_comp,
                           const amrex::MultiFab&                 volume,
                           const amrex::MultiFab* const*          area,
                           const amrex::IntVect&                  cratio,
                           const amrex::BCRec&                    bc,
                           const amrex::Geometry&                 geom,
                           bool                                   add_hoop_stress,
                           const Diffusion::SolveMode&            solve_mode = Diffusion::ONEPASS,
                           bool                                   add_old_time_divFlux = true,
                           const amrex::Vector<int>&              is_diffusive = amrex::Vector<int>());

  void diffusionFJDriver(amrex::ForkJoin&                   fj,
                  amrex::Real                        prev_time,
                  amrex::Real                        curr_time,
                  amrex::Real                        be_cn_theta,
                  int                         rho_flag,
                  const amrex::Vector<amrex::Real>&         visc_coef,
                  int                         visc_coef_comp,
                  const amrex::IntVect&              cratio,
                  const amrex::BCRec&                bc,
                  const amrex::Geometry&             geom,
                  bool                        add_hoop_stress,
                  const Diffusion::SolveMode& solve_mode,
                  bool                        add_old_time_divFlux,
                  const amrex::Vector<int>&   is_diffusive,
                  bool                        has_coarse_data,
                  bool                        has_delta_rhs,
                  bool                        has_alpha_in,
                  bool                        has_betan,
                  bool                        has_betanp1);
                           
  //
  // Functions for interpolating from cell centers to cell edges
  //
  enum FPLoc { HT_Edge = 0, HT_Center };
  //
  static FPLoc fpi_phys_loc (int p_bc);
  //
  static void center_to_edge_fancy (const amrex::FArrayBox& cfab,
                                    amrex::FArrayBox&       efab,
                                    const amrex::Box&       ccBox,
                                    const amrex::Box&       eBox,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp,
                                    const amrex::Box&       domain,
                                    const FPLoc&     bc_lo,
                                    const FPLoc&     bc_hi);

  static void init_once ();

  static void init_transport (int iEG);

  void reactionRateRhoY_pphys(amrex::FArrayBox&       RhoYdot,
                          const amrex::FArrayBox& RhoY,
                          const amrex::FArrayBox& RhoH,
                          const amrex::FArrayBox& T,
                          const amrex::BaseFab<int>& mask,
                          const amrex::Box&       box,
                          int              sCompRhoY,
                          int              sCompRhoH,
                          int              sCompT,
                          int              sCompRhoYdot) const;

  void getHmixGivenTY_pphys (amrex::FArrayBox&       hmix,
                         const amrex::FArrayBox& T,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

  static int getTGivenHY_pphys (amrex::FArrayBox&       T,
                     const amrex::FArrayBox& H,
                     const amrex::FArrayBox& Y,
                     const amrex::Box&       box,
                     int              sCompH,
                     int              sCompY,
                     int              sCompT,
                     const amrex::Real&      errMAX = -1);

  void HfromT_pphys(const int* lo, const int* hi,
		     amrex::Real*       H, ARLIM_P(Hlo), ARLIM_P(Hhi),
		     const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

  static void RhoH_to_Temp (amrex::MultiFab& S,
                            int       nGrow = 0,
                            int       dominmax = false);

  ////////////////////////////////////////////////////////////////////////////
  //    Private Data                                                        //
  ////////////////////////////////////////////////////////////////////////////

  static amrex::Vector<amrex::AMRErrorTag> errtags;

  amrex::Vector<std::unique_ptr<FluxBoxes> > raii_fbs;
    
  amrex::MultiFab**  EdgeState;
  amrex::MultiFab**  EdgeFlux;
  amrex::MultiFab**  SpecDiffusionFluxn;
  amrex::MultiFab**  SpecDiffusionFluxnp1;

  amrex::MultiFab**  SpecDiffusionFluxWbar;

  static bool                     plot_reactions;
  static bool                     plot_consumption;
  static bool                     plot_heat_release;
  std::map<std::string,std::unique_ptr<amrex::MultiFab> > auxDiag;
  static std::map<std::string,amrex::Vector<std::string> > auxDiag_names;

  bool updateFluxReg;
  bool is_predictor;

  // these refer to the old and new-time ambient pressure for level 0
  static amrex::Real p_amb_old;
  static amrex::Real p_amb_new;
  static amrex::Real dp0dt;
  static amrex::Real thetabar;

  static int closed_chamber;
  static amrex::Real dpdt_factor;

  // Active control data
  static bool        ctrl_active;
  static bool        ctrl_use_temp;
  static amrex::Real ctrl_tauControl;
  static amrex::Real ctrl_cfix;
  static amrex::Real ctrl_coftOld;
  static amrex::Real ctrl_sest;
  static amrex::Real ctrl_corr;
  static amrex::Real ctrl_V_in;
  static amrex::Real ctrl_V_in_old;
  static amrex::Real ctrl_changeMax;
  static amrex::Real ctrl_tBase;
  static amrex::Real ctrl_dV;
  static amrex::Real ctrl_scale;
  static amrex::Real ctrl_zBase;
  static amrex::Real ctrl_h;
  static amrex::Real ctrl_velMax;
  static amrex::Real ctrl_temperature;
  static int         ctrl_verbose;
  static int         ctrl_NavgPts;
  static int         ctrl_nfilled;
  static int         ctrl_flameDir;
  static int         ctrl_pseudoGravity;
  static std::string ctrl_AChistory;
  static amrex::Vector<amrex::Real> ctrl_time_pts;
  static amrex::Vector<amrex::Real> ctrl_velo_pts;
  static amrex::Vector<amrex::Real> ctrl_cntl_pts;

  static std::unique_ptr<ProbParm> prob_parm;
  static std::unique_ptr<ACParm> ac_parm;
  //
  // Protected static data.
  //
  static int         num_divu_iters;
  static int         init_once_done;
  static int         RhoH;
  static int         do_OT_radiation;
  static int         do_heat_sink;
  static int         unity_Le;
  static int         use_wbar;
  static int         FuncCount_Type;
  static int         RhoYdot_Type;
  static int         divu_ceiling;
  static amrex::Real min_rho_divu_ceiling;
  static amrex::Real divu_dt_factor;
  static amrex::Real htt_tempmin;
  static amrex::Real htt_tempmax;
  static amrex::Real htt_hmixTYP;
  static amrex::Real typical_Y_val_min;
  static int         RhoRT;
  static int         have_rhort;
  static int         first_spec;
  static int         last_spec;
  static amrex::Vector<std::string> spec_names;
  static int         floor_species;
  static amrex::Real prandtl;
  static amrex::Real schmidt;
  static amrex::Real constant_thick_val;
  static int         do_diffuse_sync;
  static int         do_reflux_visc;
  static int         do_set_rho_to_species_sum;
  static int         zeroBndryVisc;
  static int         do_check_divudt;
  static int         hack_nochem;
  static int         hack_nospecdiff;
  static int         hack_noavgdivu;
  static amrex::Real trac_diff_coef;
  static int         use_tranlib;
  static std::string turbFile;
  static std::string fuelName;
  static std::string productName;
  static amrex::Vector<std::string> consumptionName;
  static int         calcDiffusivity_count;
  static amrex::Vector<amrex::Real> typical_values;
  static bool        avg_down_chem;
  static int         reset_typical_vals_int;
  static std::map<std::string,amrex::Real> typical_values_FileVals;
  static int         sdc_iterMAX;
  static int         num_mac_sync_iter;
  static int         deltaT_verbose;
  static int         deltaT_crashOnConvFail;
  static int         iter_debug;
  static int         mHtoTiterMAX;
  static amrex::Vector<amrex::Real> mTmpData;
  

  static int  ncells_chem;
  static bool use_typ_vals_chem; 
  static amrex::Real relative_tol_chem;
  static amrex::Real absolute_tol_chem;

  static int nGrowAdvForcing;

};

void pelelm_cc_ext_fill (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);

void pelelm_press_fill (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);

void pelelm_face_fill  (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);

void pelelm_dummy_fill  (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);

#endif /*_PeleLM_H_*/
