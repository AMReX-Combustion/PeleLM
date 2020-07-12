#ifndef _PeleLM_F_H_
#define _PeleLM_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

//PelePhysics stuffs

//  int pphys_CONPsolv_SDC(const int* lo, const int* hi,
//                 amrex::Real* rhoYnew, ARLIM_P(rhoYnlo), ARLIM_P(rhoYnhi),
//                 amrex::Real* rhoHnew, ARLIM_P(rhoHnlo), ARLIM_P(rhoHnhi),
//                 amrex::Real* Tnew, ARLIM_P(Tnlo), ARLIM_P(Tnhi),
//                 const amrex::Real* rhoYold, ARLIM_P(rhoYolo), ARLIM_P(rhoYohi),
//                 const amrex::Real* rhoHold, ARLIM_P(rhoHolo), ARLIM_P(rhoHohi),
//                 const amrex::Real* Told, ARLIM_P(Tolo), ARLIM_P(Tohi),
//                 const amrex::Real* const_src, ARLIM_P(cslo), ARLIM_P(cshi),
//                 amrex::Real* FuncCount, ARLIM_P(FuncCountlo), ARLIM_P(FuncCounthi),
//                 const amrex::Real* dt, amrex::Real* diag, 
//                         const int* do_diag, const int* do_stiff);

   void pphys_network_init();

   void pphys_network_close();

  //void pphys_reactor_close();

   void pphys_get_spec_index(int* spec,
                             int* len,
                             int* idx);

   void pphys_set_verbose_vode ();

   void active_control(const amrex::Real* fuelmass,
                       const amrex::Real* time,
                       const amrex::Real* crsedt,
                       const int*  myproc,
                       const int*  step,
                       const int*  restart,
                       const int*  usetemp);

   void get_pamb(amrex::Real* pamb);

   void set_ht_adim_common(const amrex::Real* thickFacVal,
                           const amrex::Real* Pr, const amrex::Real* Sc, const int* unityLe);

   void part_cnt_err(BL_FORT_IFAB_ARG_ANYD(tag),
                     const int* tagval, const int* clearval,
                     BL_FORT_FAB_ARG_ANYD(data),
                     const int* lo, const int* hi, const int* nvar,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real* xlo,
                     const amrex::Real* prob_lo, const amrex::Real* time,
                     const int* level);

   void valgt_error(BL_FORT_IFAB_ARG_ANYD(tag),
                    const int* tagval, const int* clearval,
                    const BL_FORT_FAB_ARG_ANYD(dat),
                    const int* lo, const int* hi, const int* nvar,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* dx, const amrex::Real* xlo,
                    const amrex::Real* prob_lo, const amrex::Real* time,
                    const int* level, const amrex::Real* value);

   void vallt_error(BL_FORT_IFAB_ARG_ANYD(tag),
                    const int* tagval, const int* clearval,
                    const BL_FORT_FAB_ARG_ANYD(dat),
                    const int* lo, const int* hi, const int* nvar,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* dx, const amrex::Real* xlo,
                    const amrex::Real* prob_lo, const amrex::Real* time,
                    const int* level, const amrex::Real* value);

   void magvort_error(BL_FORT_IFAB_ARG_ANYD(tag),
                      const int* tagval, const int* clearval,
                      const BL_FORT_FAB_ARG_ANYD(dat),
                      const int* lo, const int* hi, const int* nvar,
                      const int* domain_lo, const int* domain_hi,
                      const amrex::Real* dx, const amrex::Real* xlo,
                      const amrex::Real* prob_lo, const amrex::Real* time,
                      const int* level, const amrex::Real* value);

   void diffgt_error(BL_FORT_IFAB_ARG_ANYD(tag),
                     const int* tagval, const int* clearval,
                     const BL_FORT_FAB_ARG_ANYD(dat),
                     const int* lo, const int* hi, const int* nvar,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real* xlo,
                     const amrex::Real* prob_lo, const amrex::Real* time,
                     const int* level, const amrex::Real* value);

   void box_error(BL_FORT_IFAB_ARG_ANYD(tag),
                  const int* tagval, const int* clearval,
                  const amrex::Real* blo, const amrex::Real* bhi,
                  const int* lo, const int* hi,
                  const int* domain_lo, const int* domain_hi,
                  const amrex::Real* dx, const amrex::Real* xlo,
                  const amrex::Real* prob_lo, const amrex::Real* time,
                  const int* level);

   void set_prob_spec (const int& dm,
                       const amrex::Real* problo, const amrex::Real* probhi,
                       const int* bathID,  const int* fuelID,
                       const int* oxidID,
                       const int* prodID,
                       const int* numspec,
                       const int& flag_active_control);

   void dcma_error(BL_FORT_IFAB_ARG_ANYD(tag),
                   const int* tagval, const int* clearval,
                   const BL_FORT_FAB_ARG_ANYD(dat),
                   const int* lo, const int* hi, const int* nvar,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* dx, const amrex::Real* xlo,
                   const amrex::Real* prob_lo, const amrex::Real* time,
                   const int* level, const amrex::Real* value);

#ifdef __cplusplus
}
#endif
#endif

#endif /*_PeleLM_F_H_*/
