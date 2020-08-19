#ifndef _Efield_H_
#define _Efield_H_


// Additional scalar indices
   static int         nE;        // electron density
   static int         iE_sp;     // electron species in list of specs
   static int         PhiV;      // electric potential

   static amrex::BCRec   phiV_bc;   // Separate set of BC for PhiV

// PhiV poisson solve options
   static amrex::Real ef_PoissonTol;
   static int         ef_PoissonMaxIter;
   static int         ef_PoissonVerbose;
   static int         ef_PoissonMaxOrder;

   static void ef_init();        // Init EF params

// Poisson solve for PhiV with current charge distribution
   void ef_solve_phiv(const amrex::Real &time);

// Setup PhiV BC for linear Poisson solve
   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc);

// Species additional properties (should be in PelePhysics probably)
   static amrex::GpuArray<amrex::Real, NUM_SPECIES> zk;

#endif
