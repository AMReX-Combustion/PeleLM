#ifndef _PROB_F_H_
#define _PROB_F_H_

#ifndef BL_LANG_FORT


#ifdef __cplusplus
  extern "C" {
#endif

#if 0
      void FORT_HACK(const int* lo, const int* hi,
                     const REAL* dat, ARLIM_P(lo), ARLIM_P(hi),
                     const REAL* datc, ARLIM_P(loc), ARLIM_P(hic),
                     const int* nc, const int* mapOldToNew,
                     const int* lenmap, const int* ratio);
#endif

      void init_data_new_mech (const int* level, const amrex::Real* time,
                               const int* lo, const int* hi,
                               const int* nscal,
                               amrex::Real* vel, 
                               BL_FORT_FAB_ARG_ANYD(scal),
                               BL_FORT_FAB_ARG_ANYD(press),
                               const amrex::Real* dx,
                               const amrex::Real* xlo, const amrex::Real* xhi);

      void init_data (const int* level, const amrex::Real* time,
                      const int* lo, const int* hi,
                      const int* nscal,
                      amrex::Real* vel, 
                      BL_FORT_FAB_ARG_ANYD(scal),
                      BL_FORT_FAB_ARG_ANYD(press),
                      const amrex::Real* dx,
                      const amrex::Real* xlo, const amrex::Real* xhi);

     void FORT_DERNULL (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                        const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const amrex::Real* delta, const amrex::Real* xlo,
                        const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                        const int* level, const int* grid_no) ;

     void zero_visc (BL_FORT_FAB_ARG_ANYD(visc),
                     const int* lo,   const int* hi,
                     const int* domlo, const int* domhi,
                     const amrex::Real* dx, const amrex::Real* problo,
                     const int* bc, const int* idir, const int* isrz,
                     const int* state_comp, const int* ncomp);

     void FORT_DERFORCING (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                           const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEX (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                          const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* delta, const amrex::Real* xlo,
                          const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                          const int* level, const int* grid_no) ;

     void FORT_DERFORCEY (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                          const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* delta, const amrex::Real* xlo,
                          const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                          const int* level, const int* grid_no) ;

     void FORT_DERFORCEZ (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                          const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* delta, const amrex::Real* xlo,
                          const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                          const int* level, const int* grid_no) ;

     void derUserDefined (BL_FORT_FAB_ARG_ANYD(der),const int* nvar,
                          const BL_FORT_FAB_ARG_ANYD(data),const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* delta, const amrex::Real* xlo,
                          const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
                          const int* level, const int* grid_no) ;

#ifdef __cplusplus
  }
#endif
#endif

#endif
