#ifndef _Efield_H_
#define _Efield_H_

   static int  ef_verbose;

// Additional scalar indices
   static int         nE;        // electron density
   static int         iE_sp;     // electron species in list of specs
   static int         PhiV;      // electric potential

   static amrex::BCRec   phiV_bc;   // Separate set of BC for PhiV

// Extra transport properties for electrons and charged species
   amrex::MultiFab            Ke_cc;
   amrex::MultiFab            De_cc;

// PhiV poisson solve options
   static amrex::Real ef_PoissonTol;
   static int         ef_PoissonMaxIter;
   static int         ef_PoissonVerbose;
   static int         ef_PoissonMaxOrder;

   static void ef_init();        // Init EF params

// SDC predictor : solve for PhiV with 'old' state and compute ef transport properties
   void ef_advance_setup(const amrex::Real &time);

// Poisson solve for PhiV with current charge distribution
   void ef_solve_phiv(const amrex::Real &time);

// Define some multifabs for transport properties
   void ef_define_data();

// Compute EF related transport properties: Kp_sp, Kp_e, Diff_e
   void ef_calc_transport(const amrex::Real &time);

// Setup PhiV BC for linear Poisson solve
   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc);

// Species additional properties (should be in PelePhysics probably)
   static amrex::GpuArray<amrex::Real, NUM_SPECIES> zk;

#endif
