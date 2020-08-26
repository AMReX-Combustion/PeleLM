#ifndef _Efield_H_
#define _Efield_H_

   static int  ef_verbose;
   static int  ef_debug;
   static int  ef_substep;

// Additional scalar indices
   static int         nE;        // electron density
   static int         iE_sp;     // electron species in list of specs
   static int         PhiV;      // electric potential

   static amrex::BCRec   phiV_bc;   // Separate set of BC for PhiV

// Extra transport properties for electrons and charged species
   FluxBoxes diff_e;
   amrex::MultiFab**          De_ec;

// Scaling of non-linear system
   amrex::Real                 nE_scale;
   amrex::Real                 phiV_scale;
   amrex::Real                 FnE_scale;
   amrex::Real                 FphiV_scale;

// Algo misc.
   amrex::Real                 dtsub;              // sub time step size
   amrex::Real                 curtime;            // current time within the sub-iteration
   amrex::Real                 nl_stateNorm;
   amrex::Real                 nl_residNorm;

// Class data
   amrex::MultiFab            ef_state_old;
   amrex::MultiFab            ef_state_refGhostCell;
   amrex::MultiFab            bg_charge;
   amrex::MultiFab            nl_state;
   amrex::MultiFab            nl_resid;

// PhiV poisson solve options
   static amrex::Real ef_PoissonTol;
   static int         ef_PoissonMaxIter;
   static int         ef_PoissonVerbose;
   static int         ef_PoissonMaxOrder;

// JFNK options
   static int         ef_use_PETSC_direct;      // Flag to trigger using PETSC direct solve
   static amrex::Real ef_lambda_jfnk;           // Perturbation mag. for JFNK differencing
   static int         ef_diffT_jfnk;            // Differencing type: 1 - one-sided first order, 2 - centered second order
   static int         ef_maxNewtonIter;         // Maximum number of newton iterations
   static amrex::Real ef_newtonTol;             // Tolerance of newton iteration (on inf. norm.)
   static int         ef_GMRES_size;            // Size of the GMRES Krylov base
   static amrex::Real ef_GMRES_reltol;          // Relarive tolerance of the GMRES solve
   static int         ef_GMRES_maxRst;          // Maximum number of GMRES restart

   static void ef_init();        // Init EF params

// SDC predictor : solve for PhiV with 'old' state and compute ef transport properties
   void ef_advance_setup(const amrex::Real &time);

// Poisson solve for PhiV with current charge distribution
   void ef_solve_phiv(const amrex::Real &time);

// Define some multifabs for transport properties
   void ef_define_data();

// Compute EF related transport properties: Kp_sp, Kp_e, Diff_e
   void ef_calc_transport(const amrex::Real &time);

// Solve PNP function: return updated PhiV and nE
   void ef_solve_PNP(             int       misdc,
                     const amrex::Real      &dt,
                     const amrex::Real      &time,
                     const amrex::MultiFab  &Dn,
                     const amrex::MultiFab  &Dnp1,
                     const amrex::MultiFab  &Dhat,
                           amrex::MultiFab  &ForcingnE);

// non-linear residual function
   void ef_nlResidual(const amrex::Real        &dt,
                      const amrex::MultiFab    &a_nl_state,
                            amrex::MultiFab    &a_nl_resid,
                            int                update_res_scaling = 0,
                            int                update_precond = 0);

// Laplacian of phiV
   void compPhiVLap(amrex::MultiFab    &a_phi,
                    amrex::MultiFab    &phiLap);

// Electron diffusion term
   void compElecDiffusion(amrex::MultiFab    &a_nE,
                          amrex::MultiFab    &elecDiff);

// Provisional charge distribution
   void ef_bg_chrg(const amrex::Real      &dt,
                   const amrex::MultiFab  &Dn,
                   const amrex::MultiFab  &Dnp1,
                   const amrex::MultiFab  &Dhat);

// Newton exit test
   int testExitNewton(const amrex::MultiFab  &res,
                            int              newtonIter);

// MultiFab norm function
   void ef_normMF(const amrex::MultiFab  &a_vec,
                        amrex::Real  &norm);

   void jtimesv(const amrex::MultiFab  &x, amrex::MultiFab  &Ax);

// Setup PhiV BC for linear Poisson solve
   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc);

// Setup PhiV BC for linear Poisson solve
   void ef_set_neBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                    std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc);

// Species additional properties (should be in PelePhysics probably)
   static amrex::GpuArray<amrex::Real, NUM_SPECIES> zk;

#endif
