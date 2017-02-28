#ifndef _PROB_F_H_
#define _PROB_F_H_

#include <PROB_NS_F.H>

#ifdef BL_LANG_FORT
#    define FORT_INITDATANEWMECH        initdatanewmech
#    define FORT_SET_PROB_SPEC          setpspec
#    define FORT_CENTERLINE_TEMP_STATS  cltstats
#    define FORT_FLAMETRACERROR         fterror
#    define FORT_DERNULL                dernull
#    define FORT_DERRHOMINUSSUMRHOY     drhomry
#ifdef LMC_SDC
#    define FORT_DERSUMRHOYDOT          dsrhoydot
#else
#    define FORT_DERSUMYDOT             dsydot
#endif
#    define FORT_RHOHFILL               rhohfill
#    define FORT_VELFILL                velfill
#    define FORT_CHEMFILL               chemfill
#    define FORT_ALLCHEMFILL            allchemfill
#    define FORT_DERRHORT               drhort
#    define FORT_HACK                   hack
#    define FORT_DERMOLEFRAC            dermolefrac
#    define FORT_RADLOSS                radloss
#    define FORT_ZEROVISC               zerovisc
#    define FORT_DERCONCENTRATION       derconcentration
#    define FORT_DERFORCING             drforcing
#    define FORT_DERFORCEX              drforcex
#    define FORT_DERFORCEY              drforcey
#    define FORT_DERFORCEZ              drforcez
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_INITDATANEWMECH        INITDATANEWMECH
#    define FORT_SET_PROB_SPEC          SETPSPEC
#    define FORT_CENTERLINE_TEMP_STATS  CLTSTATS
#    define FORT_FLAMETRACERROR         FTERROR
#    define FORT_DERNULL                DERNULL
#    define FORT_DERRHOMINUSSUMRHOY     DRHOMRY
#ifdef LMC_SDC
#    define FORT_DERSUMRHOYDOT          DSRHOYDOT
#else
#    define FORT_DERSUMYDOT             DSYDOT
#endif
#    define FORT_RHOHFILL               RHOHFILL
#    define FORT_VELFILL                VELFILL
#    define FORT_CHEMFILL               CHEMFILL
#    define FORT_ALLCHEMFILL            ALLCHEMFILL
#    define FORT_DERRHORT               DRHORT
#    define FORT_HACK                   HACK
#    define FORT_DERMOLEFRAC            DERMOLEFRAC
#    define FORT_RADLOSS                RADLOSS
#    define FORT_ZEROVISC               ZEROVISC
#    define FORT_DERCONCENTRATION       DERCONCENTRATION
#    define FORT_DERFORCING             DRFORCING
#    define FORT_DERFORCEX              DRFORCEX
#    define FORT_DERFORCEY              DRFORCEY
#    define FORT_DERFORCEZ              DRFORCEZ
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_INITDATANEWMECH        initdatanewmech
#    define FORT_SET_PROB_SPEC          setpspec
#    define FORT_CENTERLINE_TEMP_STATS  cltstats
#    define FORT_FLAMETRACERROR         fterror
#    define FORT_DERNULL                dernull
#    define FORT_DERRHOMINUSSUMRHOY     drhomry
#ifdef LMC_SDC
#    define FORT_DERSUMRHOYDOT          dsrhoydot
#else
#    define FORT_DERSUMYDOT             dsydot
#endif
#    define FORT_RHOHFILL               rhohfill
#    define FORT_VELFILL                velfill
#    define FORT_CHEMFILL               chemfill
#    define FORT_ALLCHEMFILL            allchemfill
#    define FORT_DERRHORT               drhort
#    define FORT_DERMOLEFRAC            dermolefrac
#    define FORT_RADLOSS                radloss
#    define FORT_ZEROVISC               zerovisc
#    define FORT_DERCONCENTRATION       derconcentration
#    define FORT_DERFORCING             drforcing
#    define FORT_DERFORCEX              drforcex
#    define FORT_DERFORCEY              drforcey
#    define FORT_DERFORCEZ              drforcez
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_INITDATANEWMECH        initdatanewmech_
#    define FORT_SET_PROB_SPEC          setpspec_
#    define FORT_CENTERLINE_TEMP_STATS  cltstats_
#    define FORT_FLAMETRACERROR         fterror_
#    define FORT_DERNULL                dernull_
#    define FORT_DERRHOMINUSSUMRHOY     drhomry_
#ifdef LMC_SDC
#    define FORT_DERSUMRHOYDOT          dsrhoydot_
#else
#    define FORT_DERSUMYDOT             dsydot_
#endif
#    define FORT_RHOHFILL               rhohfill_
#    define FORT_VELFILL                velfill_
#    define FORT_CHEMFILL               chemfill_
#    define FORT_ALLCHEMFILL            allchemfill_
#    define FORT_DERRHORT               drhort_
#    define FORT_HACK                   hack_
#    define FORT_DERMOLEFRAC            dermolefrac_
#    define FORT_RADLOSS                radloss_
#    define FORT_ZEROVISC               zerovisc_
#    define FORT_DERCONCENTRATION       derconcentration_
#    define FORT_DERFORCING             drforcing_
#    define FORT_DERFORCEX              drforcex_
#    define FORT_DERFORCEY              drforcey_
#    define FORT_DERFORCEZ              drforcez_
#  endif

  extern "C" {
#if 0
      void FORT_HACK(const int* lo, const int* hi,
                     const REAL* dat, ARLIM_P(lo), ARLIM_P(hi),
                     const REAL* datc, ARLIM_P(loc), ARLIM_P(hic),
                     const int* nc, const int* mapOldToNew,
                     const int* lenmap, const int* ratio);
#endif

      void FORT_INITDATANEWMECH (const int* level, const Real* time,
                                 const int* lo, const int* hi,
                                 const int* nscal,
                                 Real* vel, Real* scal, 
                                 ARLIM_P(s_lo), ARLIM_P(s_hi),
                                 Real* press,
                                 ARLIM_P(p_lo), ARLIM_P(p_hi),
                                 const Real* dx,
                                 const Real* xlo, const Real* xhi);

      void FORT_SET_PROB_SPEC(const int* fuelID,
			      const int* oxidID,
			      const int* prodID,
			      const int* numspec);

     void FORT_CENTERLINE_TEMP_STATS(ARLIM_P(templo), ARLIM_P(temphi),
				     const Real* tempdat, const Real* zcen,
				     const int& first_step, const Real* hy,
				     const Real* max_temp,
				     const Real* zmax_temp);
     
     void FORT_FLAMETRACERROR (int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                          const int* tagval, const int* clearval,
                          Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                          const int* lo, const int* hi, const int* nvar,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* prob_lo, const Real* time,
                          const int* level);

     void FORT_DERNULL
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERRHOMINUSSUMRHOY
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

#ifdef LMC_SDC
     void FORT_DERSUMRHOYDOT
#else
     void FORT_DERSUMYDOT
#endif
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERMOLEFRAC
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERCONCENTRATION
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_RHOHFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc);

     void FORT_VELFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                        const int* dom_lo, const int* dom_hi,
                        const Real* dx, const Real* grd_lo,
                        const Real* time, const int* bc);

     void FORT_CHEMFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const Real* dx, const Real* grd_lo,
                         const Real* time, const int* bc,
			 const int* stateID);

     void FORT_ALLCHEMFILL (Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                            const int* dom_lo, const int* dom_hi,
                            const Real* dx, const Real* grd_lo,
                            const Real* time, const int* bc);

     void FORT_DERRHORT  (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
			  const int* nvar, const Real* compdat,
			  ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			  const int* lo, const int* hi,
			  const int* domain_lo, const int* domain_hi,
			  const Real* delta, const Real* xlo,
			  const Real* time, const Real* dt, const int* bcrec, 
			  const int* level, const int* grid_no) ;
     
     void FORT_RADLOSS(const int* lo, const int* hi,
                       Real*       f, ARLIM_P(flo), ARLIM_P(fhi),
                       const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                       const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                       const Real* dx, const Real* Patm, const Real* time);
     
     void FORT_ZEROVISC (Real* data, ARLIM_P(blo), ARLIM_P(bhi),
                         const int* lo,   const int* hi,
                         const int* domlo, const int* domhi,
                         const Real* dx, const Real* problo,
                         const int* bc, const int* idir, const int* isrz,
                         const int* state_comp, const int* ncomp);

     void FORT_DERFORCING
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEX
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEY
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

     void FORT_DERFORCEZ
                          (Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                           const int* nvar, const Real* compdat,
			   ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time,  const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no) ;

  }
#endif

#endif
