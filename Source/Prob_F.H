#ifndef _PROB_F_H_
#define _PROB_F_H_

#include <PROB_NS_F.H>

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
                                 amrex::Real* vel, amrex::Real* scal, 
                                 ARLIM_P(s_lo), ARLIM_P(s_hi),
                                 amrex::Real* press,
                                 ARLIM_P(p_lo), ARLIM_P(p_hi),
                                 const amrex::Real* dx,
                                 const amrex::Real* xlo, const amrex::Real* xhi);

      void init_data (const int* level, const amrex::Real* time,
                                 const int* lo, const int* hi,
                                 const int* nscal,
                                 amrex::Real* vel, amrex::Real* scal, 
                                 ARLIM_P(s_lo), ARLIM_P(s_hi),
                                 amrex::Real* press,
                                 ARLIM_P(p_lo), ARLIM_P(p_hi),
                                 const amrex::Real* dx,
                                 const amrex::Real* xlo, const amrex::Real* xhi);
                                 

     void drhomry
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void dsrhoydot
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void derRhoY
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;
                          
     void dermassfrac
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void dermolefrac
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void derconcentration
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void dertransportcoeff
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;

     void dermolweight
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
                           ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const amrex::Real* delta, const amrex::Real* xlo,
                           const amrex::Real* time,  const amrex::Real* dt, const int* bcrec,
                           const int* level, const int* grid_no) ;
     
     void FORT_DERNULL
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			                     ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			                     const int* lo, const int* hi,
			                     const int* domain_lo, const int* domain_hi,
			                     const amrex::Real* delta, const amrex::Real* xlo,
			                     const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void rhoh_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void vel_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void xvel_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);

     void yvel_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void zvel_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void den_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void temp_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void adv_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void press_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);     
                        
     void chem_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc,
			                   const int* stateID);

     void all_chem_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                            const int* dom_lo, const int* dom_hi,
                            const amrex::Real* dx, const amrex::Real* grd_lo,
                            const amrex::Real* time, const int* bc);

     void drhort  (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
			  const int* nvar, const amrex::Real* compdat,
			  ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			  const int* lo, const int* hi,
			  const int* domain_lo, const int* domain_hi,
			  const amrex::Real* delta, const amrex::Real* xlo,
			  const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
			  const int* level, const int* grid_no) ;
     
     void zero_visc (amrex::Real* data, ARLIM_P(blo), ARLIM_P(bhi),
                         const int* lo,   const int* hi,
                         const int* domlo, const int* domhi,
                         const amrex::Real* dx, const amrex::Real* problo,
                         const int* bc, const int* idir, const int* isrz,
                         const int* state_comp, const int* ncomp);

     void FORT_DERFORCING
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEX
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEY
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEZ
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#ifdef __cplusplus
  }
#endif
#endif

#endif
