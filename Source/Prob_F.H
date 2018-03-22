#ifndef _PROB_F_H_
#define _PROB_F_H_

#include <PROB_NS_F.H>

#ifdef BL_LANG_FORT

#else
#  if defined(BL_FORT_USE_UPPERCASE)

#  elif defined(BL_FORT_USE_LOWERCASE)

#  elif defined(BL_FORT_USE_UNDERSCORE)

#  endif

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
                                 
      void FORT_SET_PROB_SPEC(const int* fuelID,
			      const int* oxidID,
			      const int* prodID,
			      const int* numspec);

     void FORT_CENTERLINE_TEMP_STATS(ARLIM_P(templo), ARLIM_P(temphi),
				     const amrex::Real* tempdat, const amrex::Real* zcen,
				     const int& first_step, const amrex::Real* hy,
				     const amrex::Real* max_temp,
				     const amrex::Real* zmax_temp);
     
     void FORT_FLAMETRACERROR (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real* xlo,
                          const amrex::Real* prob_lo, const amrex::Real* time,
                          const int* level);

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

     void FORT_DERNULL
                          (amrex::Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const amrex::Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const amrex::Real* delta, const amrex::Real* xlo,
			   const amrex::Real* time,  const amrex::Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void temp_error(int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
                     const int* tagval, const int* clearval,
                     amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                     const int* lo, const int* hi,
                     const int* ncomp,
                     const int* domlo, const int* domhi,
                     const amrex::Real* dx, const amrex::Real* grd_lo, const amrex::Real* problo,
                     const amrex::Real* time, const int* level);
    
     void mv_error(int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
                     const int* tagval, const int* clearval,
                     amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                     const int* lo, const int* hi,
                     const int* ncomp,
                     const int* domlo, const int* domhi,
                     const amrex::Real* dx, const amrex::Real* grd_lo, const amrex::Real* problo,
                     const amrex::Real* time, const int* level);
     
     void adv_error(int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
                     const int* tagval, const int* clearval,
                     amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                     const int* lo, const int* hi,
                     const int* ncomp,
                     const int* domlo, const int* domhi,
                     const amrex::Real* dx, const amrex::Real* grd_lo, const amrex::Real* problo,
                     const amrex::Real* time, const int* level);
    
     void FORT_RHOHFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

     void vel_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);
     
     void FORT_XVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const amrex::Real* dx, const amrex::Real* grd_lo,
                        const amrex::Real* time, const int* bc);

     void FORT_YVELFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
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
                        
     void FORT_CHEMFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc,
			 const int* stateID);

     void FORT_ALLCHEMFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
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
     
     void FORT_RADLOSS(const int* lo, const int* hi,
                       amrex::Real*       f, ARLIM_P(flo), ARLIM_P(fhi),
                       const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                       const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                       const amrex::Real* dx, const amrex::Real* Patm, const amrex::Real* time);
     
     void FORT_ZEROVISC (amrex::Real* data, ARLIM_P(blo), ARLIM_P(bhi),
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
