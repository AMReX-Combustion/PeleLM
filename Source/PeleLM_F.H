#ifndef _PeleLM_F_H_
#define _PeleLM_F_H_ 

#ifndef BL_LANG_FORT

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

//PelePhysics stuffs

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
