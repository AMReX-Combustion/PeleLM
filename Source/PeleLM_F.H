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

   void pphys_RRATERHOY(const int* lo, const int* hi,
                        const BL_FORT_FAB_ARG_ANYD(RhoY),
                        const BL_FORT_FAB_ARG_ANYD(RhoH),
                        const BL_FORT_FAB_ARG_ANYD(T),
                              BL_FORT_FAB_ARG_ANYD(RhoYdot));

   void pphys_PfromRTY(const int* lo, const int* hi,
                             BL_FORT_FAB_ARG_ANYD(P),
                       const BL_FORT_FAB_ARG_ANYD(R),
                       const BL_FORT_FAB_ARG_ANYD(T),
                       const BL_FORT_FAB_ARG_ANYD(Y));

   void pphys_HMIXfromTY(const int* lo, const int* hi,
                               BL_FORT_FAB_ARG_ANYD(Hmix),
                         const BL_FORT_FAB_ARG_ANYD(T),
                         const BL_FORT_FAB_ARG_ANYD(Y));

   void pphys_HfromT(const int* lo, const int* hi,
                           BL_FORT_FAB_ARG_ANYD(H),
                     const BL_FORT_FAB_ARG_ANYD(T));

   void pphys_MWMIXfromY(const int* lo, const int* hi,
                               BL_FORT_FAB_ARG_ANYD(MWMIX),
                         const BL_FORT_FAB_ARG_ANYD(Y));

   void pphys_CPMIXfromTY(const int* lo, const int* hi,
                                BL_FORT_FAB_ARG_ANYD(CPMIX),
                          const BL_FORT_FAB_ARG_ANYD(T),
                          const BL_FORT_FAB_ARG_ANYD(Y));

   int pphys_TfromHY(const int* lo, const int* hi,
                           BL_FORT_FAB_ARG_ANYD(T), 
                     const BL_FORT_FAB_ARG_ANYD(H),
                     const BL_FORT_FAB_ARG_ANYD(Y),
                     const amrex::Real* errMAX, const int* NiterMAX,
                     const amrex::Real* tmp);

  //void pphys_TfromHYpt(amrex::Real* T, const amrex::Real* Hin, const amrex::Real* Yin,
  //                      const amrex::Real* errMAX, const int* NiterMAX, const amrex::Real* res,
  //                      int* Niter);

   void pphys_network_init();

   void pphys_network_close();

  //void pphys_reactor_close();

   void pphys_transport_init(int* ieg);

   void pphys_transport_close();

   void plm_extern_init(const int* name, const int* namlen);

   void pphys_get_num_spec(int* nspec);

   void pphys_get_spec_name(int* spec_names,
                            int* ispec,
                            int* len);

   void pphys_get_spec_index(int* spec,
                             int* len,
                             int* idx);

   amrex::Real pphys_getRuniversal();

   amrex::Real pphys_getP1atm_MKS();

   int pphys_numReactions();

   void pphys_set_verbose_vode ();

   void pphys_calc_src_sdc(int* n, amrex::Real* time,  amrex::Real* temp,
                           amrex::Real* Z, amrex::Real* ZP);

   void active_control(const amrex::Real* fuelmass,
                       const amrex::Real* time,
                       const amrex::Real* crsedt,
                       const int*  myproc,
                       const int*  step,
                       const int*  restart,
                       const int*  usetemp);

   void set_scal_numb(const int* Density,
                      const int* Temp,
                      const int* Trac,
                      const int* RhoH,
                      const int* FirstSpec,
                      const int* LastSpec);

   void get_pamb(amrex::Real* pamb);

   void divu_fill (BL_FORT_FAB_ARG_ANYD(data),
                   const int* dom_lo, const int* dom_hi,
                   const amrex::Real* dx, const amrex::Real* grd_lo,
                   const amrex::Real* time, const int* bc);

   void dsdt_fill (BL_FORT_FAB_ARG_ANYD(data),
                   const int* dom_lo, const int* dom_hi,
                   const amrex::Real* dx, const amrex::Real* grd_lo,
                   const amrex::Real* time, const int* bc);

   void dqrad_fill (BL_FORT_FAB_ARG_ANYD(data),
                    const int* dom_lo, const int* dom_hi,
                    const amrex::Real* dx, const amrex::Real* grd_lo,
                    const amrex::Real* time, const int* bc);

   void ydot_fill (BL_FORT_FAB_ARG_ANYD(data),
                   const int* dom_lo, const int* dom_hi,
                   const amrex::Real* dx, const amrex::Real* grd_lo,
                   const amrex::Real* time, const int* bc);

   void rhoYdot_fill (BL_FORT_FAB_ARG_ANYD(data),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

   void est_divu_dt(const int& flag,  const amrex::Real* divu_dt_factor,
                    const amrex::Real* dx, 
                    const BL_FORT_FAB_ARG_ANYD(sdat),
                    const BL_FORT_FAB_ARG_ANYD(dsdtdat),
                    const BL_FORT_FAB_ARG_ANYD(rhodat),
                    const BL_FORT_FAB_ARG_ANYD(vel),
                    const BL_FORT_FAB_ARG_ANYD(volume),
                    const BL_FORT_FAB_ARG_ANYD(areax),
                    const BL_FORT_FAB_ARG_ANYD(areay),
#if (AMREX_SPACEDIM == 3) 
                    const BL_FORT_FAB_ARG_ANYD(areaz),
#endif
                    const int* lo,const int* hi,
                    amrex::Real* dt, const amrex::Real* min_rho_divu_ceiling);

   void check_divu_dt(const int& flag, const amrex::Real* divu_dt_factor,
                    const amrex::Real* dx, 
                    const BL_FORT_FAB_ARG_ANYD(sdat),
                    const BL_FORT_FAB_ARG_ANYD(dsdtdat),
                    const BL_FORT_FAB_ARG_ANYD(rhodat),
                    const BL_FORT_FAB_ARG_ANYD(vel),
                    const BL_FORT_FAB_ARG_ANYD(volume),
                    const BL_FORT_FAB_ARG_ANYD(areax),
                    const BL_FORT_FAB_ARG_ANYD(areay),
#if (AMREX_SPACEDIM == 3) 
                    const BL_FORT_FAB_ARG_ANYD(areaz),
#endif
                    const int* lo,const int* hi,
                    const amrex::Real* dt, const amrex::Real* min_rho_divu_ceiling);

   void set_ht_adim_common(const amrex::Real* thickFacVal,
                           const amrex::Real* Pr, const amrex::Real* Sc, const int* unityLe);

   void vel_visc(const int* lo, const int* hi,
                 const BL_FORT_FAB_ARG_ANYD(T),
                 const BL_FORT_FAB_ARG_ANYD(Y),
                       BL_FORT_FAB_ARG_ANYD(mu));

   void beta_wbar(const int* lo, const int* hi,
                  const BL_FORT_FAB_ARG_ANYD(RD),
                        BL_FORT_FAB_ARG_ANYD(RD_Wbar),
                  const BL_FORT_FAB_ARG_ANYD(Y));

   void spec_temp_visc(const int* lo, const int* hi,
                       const BL_FORT_FAB_ARG_ANYD(T),
                       const BL_FORT_FAB_ARG_ANYD(Y),
                             BL_FORT_FAB_ARG_ANYD(beta),
                       const int* nCompD, const amrex::Real* P1ATM_MKS,
                       const int* do_temp, const int* do_VelVisc, 
                       const amrex::Real* Pamb_in);

   void fab_minmax (const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_ANYD(fabdat),
                    const amrex::Real* min, const amrex::Real* max, const int* nc);

   void repair_flux(const int* lo,  const int* hi,
                    const int* dlo, const int* dhi,
                          BL_FORT_FAB_ARG_ANYD(flux),
                    const BL_FORT_FAB_ARG_ANYD(RY),
                    const int* dir, const int* bc);

   void enth_diff_terms(const int* lo, const int* hi, 
                        const int* dlo, const int* dhi, const amrex::Real* dx,
                        const BL_FORT_FAB_ARG_ANYD(T),

                        const BL_FORT_FAB_ARG_ANYD(rDx),
                        const BL_FORT_FAB_ARG_ANYD(Fix),
                        const BL_FORT_FAB_ARG_ANYD(ax),
                        const BL_FORT_FAB_ARG_ANYD(rDy),
                        const BL_FORT_FAB_ARG_ANYD(Fiy),
                        const BL_FORT_FAB_ARG_ANYD(ay),
#if AMREX_SPACEDIM == 3
                        const BL_FORT_FAB_ARG_ANYD(rDz),
                        const BL_FORT_FAB_ARG_ANYD(Fiz),
                        const BL_FORT_FAB_ARG_ANYD(az),
#endif
                        const int* Tbc);

   void incrwext_flx_div (const int* lo, const int* hi,
                          const BL_FORT_FAB_ARG_ANYD(xflux),
                          const BL_FORT_FAB_ARG_ANYD(yflux),
#if AMREX_SPACEDIM == 3
                          const BL_FORT_FAB_ARG_ANYD(zflux),
#endif
                          const BL_FORT_FAB_ARG_ANYD(sto),
                                BL_FORT_FAB_ARG_ANYD(stn),
                          const BL_FORT_FAB_ARG_ANYD(vol),
                          const int* nc, const amrex::Real* dt);

   void flux_div (const int* lo, const int* hi,
                        BL_FORT_FAB_ARG_ANYD(update), 
                  const BL_FORT_FAB_ARG_ANYD(xflux),
                  const BL_FORT_FAB_ARG_ANYD(yflux),
#if AMREX_SPACEDIM == 3
                  const BL_FORT_FAB_ARG_ANYD(zflux),
#endif
                  const BL_FORT_FAB_ARG_ANYD(vol),
                  const int* nc, const amrex::Real* scale);

   void calc_divu_fortran(const int* lo,  const int* hi,
                                BL_FORT_FAB_ARG_ANYD(divu),
                          const BL_FORT_FAB_ARG_ANYD(rYdot),  
                          const BL_FORT_FAB_ARG_ANYD(vtY),  
                          const BL_FORT_FAB_ARG_ANYD(vtT),  
                          const BL_FORT_FAB_ARG_ANYD(rhoY),  
                          const BL_FORT_FAB_ARG_ANYD(T));

   void calc_gamma_pinv(const int* lo,  const int* hi,
                              BL_FORT_FAB_ARG_ANYD(theta),
                        const BL_FORT_FAB_ARG_ANYD(rhoY),
                        const BL_FORT_FAB_ARG_ANYD(T),
                        const amrex::Real* Pamb_in); 

   void floor_spec(const int* lo,  const int* hi, 
                   BL_FORT_FAB_ARG_ANYD(spec));

   void set_common (const amrex::Real* time, const int* iteration);

   void compute_ugradp (const BL_FORT_FAB_ARG_ANYD(p),
                              BL_FORT_FAB_ARG_ANYD(ugradp),
                        const BL_FORT_FAB_ARG_ANYD(umac),
                        const BL_FORT_FAB_ARG_ANYD(vmac),
#if AMREX_SPACEDIM == 3
                        const BL_FORT_FAB_ARG_ANYD(wmac),
#endif
                        const int* lo, const int* hi, const amrex::Real* dx);

   int conservative_T_floor(const int* lo, const int* hi,
                            BL_FORT_FAB_ARG_ANYD(fabdat),
                            const amrex::Real* min_T_fine, const int* Tcomp, const int* Rcomp,
                            const int* first_comp, const int* last_spec, const int* RhoH,
                            const int* ratio, amrex::Real* tmp, const int* nt);

   void part_cnt_err(BL_FORT_IFAB_ARG_ANYD(tag),
                     const int* tagval, const int* clearval,
                     BL_FORT_FAB_ARG_ANYD(data),
                     const int* lo, const int* hi, const int* nvar,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real* xlo,
                     const amrex::Real* prob_lo, const amrex::Real* time,
                     const int* level);

   void recomp_update(const int* lo, const int* hi,
                            BL_FORT_FAB_ARG_ANYD(update),
                      const BL_FORT_FAB_ARG_ANYD(xflux),
                      const BL_FORT_FAB_ARG_ANYD(yflux),
#if AMREX_SPACEDIM == 3
                      const BL_FORT_FAB_ARG_ANYD(zflux),
#endif
                      const BL_FORT_FAB_ARG_ANYD(volume),
                      const int* nc);

   void mcurve(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_ANYD(T),
                     BL_FORT_FAB_ARG_ANYD(K),
                     BL_FORT_FAB_ARG_ANYD(nW),
               const amrex::Real* delta);

   void smooth(const int* lo, const int* hi,
                     BL_FORT_FAB_ARG_ANYD(i),
               const BL_FORT_FAB_ARG_ANYD(o));

   void grad_wbar(const int* lo, const int* hi,
                  const BL_FORT_FAB_ARG_ANYD(Wbar),
                  const BL_FORT_FAB_ARG_ANYD(rDe),
                        BL_FORT_FAB_ARG_ANYD(flx),
                  const BL_FORT_FAB_ARG_ANYD(a),
                  const amrex::Real* delta, const int* dir,
                  const amrex::Real* mult, const amrex::Real* inc);

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
