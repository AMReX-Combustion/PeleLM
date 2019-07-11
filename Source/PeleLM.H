
#ifndef _PeleLM_H_
#define _PeleLM_H_

#include <AMReX_AuxBoundaryData.H>
#include <NavierStokesBase.H>

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
  typedef void (*LMEF)(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                       const int* tagval, const int* clearval,
                       const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                       const int* lo, const int* hi, const int* nvar,
                       const int* domain_lo, const int* domain_hi,
                       const amrex::Real* dx, const amrex::Real* xlo,
                       const amrex::Real* prob_lo, const amrex::Real* time,
                       const int* level, const amrex::Real* value);

  typedef void (*LMEF_BOX) (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
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
  
  void tagCells(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                const int* tagval, const int* clearval,
                const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                const int* lo, const int* hi, const int* nvar,
                const int* domain_lo, const int* domain_hi,
                const amrex::Real* dx, const amrex::Real* xlo,
                const amrex::Real* prob_lo, const amrex::Real* time,
                const int* level) const;

  void tagCells1(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
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
    
protected:

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
  // Note: these two functions must be supplied in NS_setup.C
  //       If divu (dsdt) is not included in the state, then
  //       they can be no-op functions
  //       See example in the _tube_flow problem dependent directory
  //
  virtual void calc_divu (amrex::Real      time,
                          amrex::Real      dt,
                          amrex::MultiFab& divu) override;

  virtual void calc_dsdt (amrex::Real      time,
                          amrex::Real      dt,
                          amrex::MultiFab& dsdt) override;

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
  static void Finalize ();
    
#ifdef AMREX_PARTICLES
  virtual int timestamp_num_extras () override;
  virtual void timestamp_add_extras (int lev, amrex::Real time, amrex::MultiFab& mf) override;
  static void read_particle_params ();
#endif

private:

  // enum YdotAction { HT_EstimateYdotNew, HT_ImproveYdotOld, HT_LeaveYdotAlone };
  // enum Solver_Status {HT_InProgress, HT_Stalled, HT_Solved};

  void adjust_spec_diffusion_fluxes (amrex::Real time);

  void advance_chemistry (amrex::MultiFab&       mf_old,
                          amrex::MultiFab&       mf_new,
                          amrex::Real            dt,
                          const amrex::MultiFab& Force,
                          int             nCompF,
                          bool            use_stiff_solver = true);

#ifdef USE_WBAR
  void calcDiffusivity_Wbar (const amrex::Real time);
#endif

  void calc_dpdt (amrex::Real      time,
                  amrex::Real      dt,
                  amrex::MultiFab& dpdt,
                  amrex::MultiFab* u_mac);

  void checkTimeStep (amrex::Real dt);

  void compute_differential_diffusion_fluxes (const amrex::Real& time,
                                              const amrex::Real& dt);

  void compute_differential_diffusion_terms (amrex::MultiFab& D,
                                             amrex::MultiFab& DD,
#ifdef USE_WBAR
                                             amrex::MultiFab& DWbar,
#endif
                                             amrex::Real      time,
                                             amrex::Real      dt);

  void compute_enthalpy_fluxes (amrex::Real                   time,
                                const amrex::MultiFab* const* beta);

  enum HowToFillGrow {HT_ZERO_GROW_CELLS, HT_EXTRAP_GROW_CELLS, HT_NUM_GROW_OPTIONS};

  void compute_instantaneous_reaction_rates (amrex::MultiFab&       R,
                                             const amrex::MultiFab& S,
                                             amrex::Real            time,
                                             int             nGrow = 0,
                                             HowToFillGrow   how = HT_ZERO_GROW_CELLS);

  void compute_rhohmix (amrex::Real       time,
                        amrex::MultiFab& rhohmix);

  void compute_rhoRT (const amrex::MultiFab& S,
                      amrex::MultiFab&       P,
                      int             pComp,
                      const amrex::MultiFab* T = 0);

  void compute_scalar_advection_fluxes_and_divergence (const amrex::MultiFab& Force,
                                                       const amrex::MultiFab& divu,
                                                       amrex::Real            dt);

  void compute_vel_visc (amrex::Real      time,
                         amrex::MultiFab* beta);

#ifdef USE_WBAR
  void compute_Wbar_fluxes(amrex::Real time,
                           amrex::Real inc);
#endif

  void define_data ();

  void differential_diffusion_update (amrex::MultiFab& Force,
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

#ifdef USE_WBAR
  void getDiffusivity_Wbar (amrex::MultiFab* diffusivity[BL_SPACEDIM],
                            const amrex::Real time);
#endif

  amrex::DistributionMapping getFuncCountDM (const amrex::BoxArray& bxba, int ngrow);

  PeleLM& getLevel (int lev)
    {
      return *(PeleLM*) &parent->getLevel(lev);
    }

  void initDataOtherTypes ();

  void post_init_press (amrex::Real&        dt_init,
                        amrex::Vector<int>&  nc_save,
                        amrex::Vector<amrex::Real>& dt_save);

  amrex::Real predict_velocity (amrex::Real dt);

  void set_htt_hmixTYP ();

  void set_reasonable_grow_cells_for_R(amrex::Real time);

  void set_rho_to_species_sum (amrex::MultiFab& S,
                               int       strtcomp, 
                               int       nghost_in = 0,
                               int       minzero = 0);

  void set_rho_to_species_sum (amrex::MultiFab& S_in,
                               int       s_in_start,
                               amrex::MultiFab& S_out,
                               int       s_out_start, 
                               int       nghost_in = 0,
                               int       minzero = 0);

  void setThermoPress(amrex::Real time);

  void set_typical_values(bool restart);
  void reset_typical_values(const amrex::MultiFab& S);

  void temperature_stats (amrex::MultiFab& S);

  void zeroBoundaryVisc (amrex::MultiFab*  beta[BL_SPACEDIM],
                         const amrex::Real time,
                         const int  state_comp,
                         const int  dst_comp,
                         const int  ncomp) const;

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

  static void init_extern ();

  static void init_network ();

  static void init_transport (int iEG);

  //bool solveChemistry_sdc(amrex::FArrayBox&        rhoYnew,
  //                          amrex::FArrayBox&        rhoHnew,
  //                          amrex::FArrayBox&        Tnew,
  //                          const amrex::FArrayBox&  rhoYold,
  //                          const amrex::FArrayBox&  rhoHold,
  //                          const amrex::FArrayBox&  Told,
  //                          const amrex::FArrayBox&  const_src,
  //                          amrex::FArrayBox&        FuncCount,
  //                          const amrex::Box&        box,
  //                          int               sComprhoY,
  //                          int               sComprhoH,
  //                          int               sCompT,
  //                          amrex::Real              dt,
  //                          amrex::FArrayBox*        chemDiag,
  //                          bool              use_stiff_solver = true) const;

  void reactionRateRhoY_pphys(amrex::FArrayBox&       RhoYdot,
                          const amrex::FArrayBox& RhoY,
                          const amrex::FArrayBox& RhoH,
                          const amrex::FArrayBox& T,
                          const amrex::Box&       box,
                          int              sCompRhoY,
                          int              sCompRhoH,
                          int              sCompT,
                          int              sCompRhoYdot) const;

  void getPGivenRTY_pphys(amrex::FArrayBox&       p,
                         const amrex::FArrayBox& Rho,
                         const amrex::FArrayBox& T,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompR,
                         int              sCompT,
                         int              sCompY,
                         int              sCompP) const;

  void getMwmixGivenY_pphys(amrex::FArrayBox&       mwmix,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompY,
                         int              sCompMw) const;

  void getHmixGivenTY_pphys (amrex::FArrayBox&       hmix,
                         const amrex::FArrayBox& T,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

  void getHGivenT_pphys (amrex::FArrayBox&       h,
                     const amrex::FArrayBox& T,
                     const amrex::Box&       box,
                     int              sCompT,
                     int              sCompH) const;

  void getCpmixGivenTY_pphys (amrex::FArrayBox&       cpmix,
                          const amrex::FArrayBox& T,
                          const amrex::FArrayBox& Y,
                          const amrex::Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCp) const;

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

  static int RhoH_to_Temp (amrex::FArrayBox& S,
                           const amrex::Box& box,
                           int        dominmax = false);

  static int RhoH_to_Temp_DoIt(amrex::FArrayBox&       T,
		  const amrex::FArrayBox& H,
		  const amrex::FArrayBox& Y,
		  const amrex::Box&       box,
		  int              sCompH, 
		  int              sCompY,
		  int              sCompT,  
		  const amrex::Real&      htt_hmixTYP );


  ////////////////////////////////////////////////////////////////////////////
  //    Private Data                                                        //
  ////////////////////////////////////////////////////////////////////////////

  amrex::Vector<std::unique_ptr<FluxBoxes> > raii_fbs;
    
  amrex::MultiFab**  EdgeState;
  amrex::MultiFab**  EdgeFlux;
  amrex::MultiFab**  SpecDiffusionFluxn;
  amrex::MultiFab**  SpecDiffusionFluxnp1;
#ifdef USE_WBAR
  amrex::MultiFab**  SpecDiffusionFluxWbar;
#endif

  amrex::MultiFab    sumSpecFluxDotGradHn;
  amrex::MultiFab    sumSpecFluxDotGradHnp1;

#ifdef USE_WBAR
  amrex::MultiFab    diffWbar_cc;
#endif

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

  //
  // Protected static data.
  //
  static int         num_divu_iters;
  static int         init_once_done;
  static int         RhoH;
  static int         do_OT_radiation;
  static int         do_heat_sink;
  static int         unity_Le;
  static int         FuncCount_Type;
  static int         RhoYdot_Type;
  static int         divu_ceiling;
  static amrex::Real min_rho_divu_ceiling;
  static amrex::Real divu_dt_factor;
  static amrex::Real htt_tempmin;
  static amrex::Real htt_tempmax;
  static amrex::Real htt_hmixTYP;
  static amrex::Real typical_Y_val_min;
  static int         Trac;
  static int         have_trac;
  static int         RhoRT;
  static int         have_rhort;
  static int         first_spec;
  static int         last_spec;
  static int         nspecies;
  static int         nreactions;
  static amrex::Vector<std::string> spec_names;
  static int         floor_species;
  static amrex::Real rgas;
  static amrex::Real prandtl;
  static amrex::Real schmidt;
  static amrex::Real constant_thick_val;
  static int         do_diffuse_sync;
  static int         do_reflux_visc;
  static int         dpdt_option;
  static int         do_set_rho_to_species_sum;
  static int         zeroBndryVisc;
  static int         do_add_nonunityLe_corr_to_rhoh_adv_flux;
  static int         do_check_divudt;
  static int         hack_nochem;
  static int         hack_nospecdiff;
  static int         hack_noavgdivu;
  static amrex::Real trac_diff_coef;
  static int         use_tranlib;
  static amrex::Real P1atm_MKS;
  static std::string turbFile;
  static std::string fuelName;
  static std::string productName;
  static amrex::Vector<std::string> consumptionName;
  static int         calcDiffusivity_count;
  static amrex::Real new_T_threshold;
  static amrex::Vector<amrex::Real> typical_values;
  static bool        avg_down_chem;
  static int         reset_typical_vals_int;
  static std::map<std::string,amrex::Real> typical_values_FileVals;
  static int         sdc_iterMAX;
  static int         num_mac_sync_iter;
  static int         iter_debug;
  static int         mHtoTiterMAX;
  static amrex::Vector<amrex::Real> mTmpData;
  
  static bool         flag_active_control;

  static int  cvode_iE;
  static int  cvode_ncells;

  static int nGrowAdvForcing;

  static std::string probin_file;
};

#endif /*_PeleLM_H_*/
