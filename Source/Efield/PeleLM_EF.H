#ifndef _Efield_H_
#define _Efield_H_

// Additional scalar indices
   static int         nE;        // electron density 
   static int         iE_sp;     // electron species in list of specs
   static int         PhiV;      // electric potential

   static amrex::BCRec   phiV_bc;   // Separate set of BC for PhiV
#endif
