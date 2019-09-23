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
//	                 amrex::Real* rhoYnew, ARLIM_P(rhoYnlo), ARLIM_P(rhoYnhi),
//	                 amrex::Real* rhoHnew, ARLIM_P(rhoHnlo), ARLIM_P(rhoHnhi),
//	                 amrex::Real* Tnew, ARLIM_P(Tnlo), ARLIM_P(Tnhi),
//	                 const amrex::Real* rhoYold, ARLIM_P(rhoYolo), ARLIM_P(rhoYohi),
//	                 const amrex::Real* rhoHold, ARLIM_P(rhoHolo), ARLIM_P(rhoHohi),
//	                 const amrex::Real* Told, ARLIM_P(Tolo), ARLIM_P(Tohi),
//	                 const amrex::Real* const_src, ARLIM_P(cslo), ARLIM_P(cshi),
//	                 amrex::Real* FuncCount, ARLIM_P(FuncCountlo), ARLIM_P(FuncCounthi),
//	                 const amrex::Real* dt, amrex::Real* diag, 
//                         const int* do_diag, const int* do_stiff);

  void pphys_RRATERHOY(const int* lo, const int* hi,
                        const amrex::Real* RhoY, ARLIM_P(RhoYlo),    ARLIM_P(RhoYhi),
                        const amrex::Real* RhoH, ARLIM_P(RhoHlo),    ARLIM_P(RhoHhi),
                        const amrex::Real* T,    ARLIM_P(Tlo),       ARLIM_P(Thi),
                        amrex::Real*    RhoYdot, ARLIM_P(RhoYdotlo), ARLIM_P(RhoYdothi) );

  void pphys_PfromRTY(const int* lo, const int* hi,
		       amrex::Real*       P, ARLIM_P(Plo), ARLIM_P(Phi),
		       const amrex::Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
		       const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
		       const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));

  void pphys_HMIXfromTY(const int* lo, const int* hi,
			 amrex::Real*    HMIX, ARLIM_P(Hlo), ARLIM_P(Hhi),
			 const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));

  void pphys_HfromT(const int* lo, const int* hi,
		     amrex::Real*       H, ARLIM_P(Hlo), ARLIM_P(Hhi),
		     const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

  void pphys_MWMIXfromY(const int* lo, const int* hi,
			 amrex::Real*    MWMIX,ARLIM_P(MWlo),ARLIM_P(MWhi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));

  void pphys_CPMIXfromTY(const int* lo, const int* hi,
                         BL_FORT_FAB_ARG_3D(CPMIX),
                         const BL_FORT_FAB_ARG_3D(T),
                         const BL_FORT_FAB_ARG_3D(Y));

  int pphys_TfromHY(const int* lo, const int* hi,
		     amrex::Real* T,       ARLIM_P(Tlo), ARLIM_P(Thi),
		     const amrex::Real* h, ARLIM_P(hlo), ARLIM_P(hhi),
		     const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
		     const amrex::Real* errMAX, const int* NiterMAX,const amrex::Real* tmp);

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
                          const int* RhoH,
                          const int* FirstSpec,
                          const int* LastSpec);
    
  void get_pamb(amrex::Real* pamb);

  void compute_rho_dgrad_hdot_grad_Y(const amrex::Real* dx,
                                      const int* lo, const int* hi,
                                      ARLIM_P(speclo),ARLIM_P(spechi), const amrex::Real* specdat,
                                      ARLIM_P(hlo),ARLIM_P(hhi), const amrex::Real* hdat,
                                      ARLIM_P(betaxlo),ARLIM_P(betaxhi), const amrex::Real* betax,
                                      ARLIM_P(betaylo),ARLIM_P(betayhi), const amrex::Real* betay,
#if (BL_SPACEDIM == 3) 
                                      ARLIM_P(betazlo),ARLIM_P(betazhi), const amrex::Real* betaz,
#endif            
                                      ARLIM_P(ddrhslo),ARLIM_P(ddrhshi), amrex::Real* ddrhsdat);

  void divu_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void dsdt_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void dqrad_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                       const int* dom_lo, const int* dom_hi,
                       const amrex::Real* dx, const amrex::Real* grd_lo,
                       const amrex::Real* time, const int* bc);

  void ydot_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void rhoYdot_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

  void est_divu_dt(const int& flag,  const amrex::Real* divu_dt_factor,
                        const amrex::Real* dx, 
                        const amrex::Real* sdat,ARLIM_P(slo),ARLIM_P(shi),
                        const amrex::Real* dsdtdat,ARLIM_P(dslo),ARLIM_P(dshi),
                        const amrex::Real* rhodat,ARLIM_P(rholo),ARLIM_P(rhohi),
                        const amrex::Real* vel,ARLIM_P(ulo),ARLIM_P(uhi),
                        const amrex::Real* volume, ARLIM_P (v_lo), ARLIM_P (v_hi), 
                        const amrex::Real* areax,
                        ARLIM_P(areax_lo),ARLIM_P(areax_hi),
                        const amrex::Real* areay,
                        ARLIM_P(areay_lo),ARLIM_P(areay_hi),
#if (BL_SPACEDIM == 3) 
                        const amrex::Real* areaz,
                        ARLIM_P(areaz_lo),ARLIM_P(areaz_hi),
#endif            
                        const int* lo,const int* hi,
                        amrex::Real* dt, const amrex::Real* min_rho_divu_ceiling);

  void check_divu_dt(const int& flag, const amrex::Real* divu_dt_factor,
                          const amrex::Real* dx, 
                          const amrex::Real* sdat, ARLIM_P(slo), ARLIM_P(shi),
                          const amrex::Real* dsdtdat,
                          const amrex::Real* rhodat,ARLIM_P(rholo),ARLIM_P(rhohi),
                          const amrex::Real* vel,ARLIM_P(ulo),ARLIM_P(uhi),
                          const amrex::Real* volume, ARLIM_P (v_lo), ARLIM_P (v_hi), 
                          const amrex::Real* areax,
                          ARLIM_P(areax_lo),ARLIM_P(areax_hi),
                          const amrex::Real* areay,
                          ARLIM_P(areay_lo),ARLIM_P(areay_hi),
#if (BL_SPACEDIM == 3) 
                          const amrex::Real* areaz,
                          ARLIM_P(areaz_lo),ARLIM_P(areaz_hi),
#endif            
                          const int* lo,const int* hi,
                          const amrex::Real* dt, const amrex::Real* min_rho_divu_ceiling);

  
  void set_ht_adim_common(const amrex::Real* thickFacVal,
                               const amrex::Real* Pr, const amrex::Real* Sc, const int* unityLe);

  void vel_visc(const int* lo, const int* hi,
                const BL_FORT_FAB_ARG_3D(T),
                const BL_FORT_FAB_ARG_3D(Y),
                BL_FORT_FAB_ARG_3D(mu));

  void beta_wbar(const int* lo, const int* hi,
                 const BL_FORT_FAB_ARG_3D(RD),
                 BL_FORT_FAB_ARG_3D(RD_Wbar),
                 const BL_FORT_FAB_ARG_3D(Y));

  void spec_temp_visc(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_3D(T),
                      const BL_FORT_FAB_ARG_3D(Y),
                      BL_FORT_FAB_ARG_3D(beta),
                      const int* nCompD, const amrex::Real* P1ATM_MKS,
                      const int* do_temp, const int* do_VelVisc, 
                      const amrex::Real* Pamb_in);

  void fab_minmax (const int* lo, const int* hi,
                      amrex::Real* fabdat, ARLIM_P(flo), ARLIM_P(fhi),
                      const amrex::Real* min, const amrex::Real* max, const int* nc);

  void repair_flux(const int* lo,  const int* hi,
                        const int* dlo, const int* dhi,
                        amrex::Real* flux,    ARLIM_P(flo),  ARLIM_P(fhi),
                        const amrex::Real* RY,ARLIM_P(RYlo), ARLIM_P(RYhi),
                        const int* dir, const int* bc);

  void enth_diff_terms(const int* lo, const int* hi, 
                       const int* dlo, const int* dhi, const amrex::Real* dx,
                       const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),

                       const amrex::Real* rDx,ARLIM_P(rDxlo),ARLIM_P(rDxhi),
                       const amrex::Real* Fix,ARLIM_P(Fixlo),ARLIM_P(Fixhi),
                       const amrex::Real* ax, ARLIM_P(axlo), ARLIM_P(axhi),

                       const amrex::Real* rDy,ARLIM_P(rDylo),ARLIM_P(rDyhi),
                       const amrex::Real* Fiy,ARLIM_P(Fiylo),ARLIM_P(Fiyhi),
                       const amrex::Real* ay, ARLIM_P(aylo), ARLIM_P(ayhi),
#if BL_SPACEDIM == 3
                       const amrex::Real* rDz,ARLIM_P(rDzlo),ARLIM_P(rDzhi),
                       const amrex::Real* Fiz,ARLIM_P(Fizlo),ARLIM_P(Fizhi),
                       const amrex::Real* az, ARLIM_P(azlo), ARLIM_P(azhi),
#endif
                       const int* Tbc);

  void incrwext_flx_div (const int* lo, const int* hi,
                           const amrex::Real* xflux, ARLIM_P(fxlo),  ARLIM_P(fxhi),
                           const amrex::Real* yflux, ARLIM_P(fylo),  ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                           const amrex::Real* zflux, ARLIM_P(fzlo),  ARLIM_P(fzhi),
#endif
                           const amrex::Real* sto,     ARLIM_P(stolo),    ARLIM_P(stohi),
                           amrex::Real*       stn,     ARLIM_P(stnlo),    ARLIM_P(stnhi),
                           const amrex::Real* vol,     ARLIM_P(vollo),    ARLIM_P(volhi),
                           const int* nc, const amrex::Real* dt);

  void flux_div(const int* lo, const int* hi,
                    const amrex::Real* updat, ARLIM_P(uplo), ARLIM_P(uphi),
                    const amrex::Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
                    const amrex::Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                    const amrex::Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
                    const amrex::Real* vol,   ARLIM_P(vlo),  ARLIM_P(vhi),
                    const int* nc, const amrex::Real* scale);

  void calc_divu_fortran(const int* lo,  const int* hi,
                     amrex::Real* divu,         ARLIM_P(dulo),    ARLIM_P(duhi),
                     const amrex::Real* rYdot,  ARLIM_P(rYdotlo), ARLIM_P(rYdothi),
                     const amrex::Real* vtY,    ARLIM_P(vtYlo),   ARLIM_P(vtYhi),
                     const amrex::Real* vtT,    ARLIM_P(vtTlo),   ARLIM_P(vtThi),
                     const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                     const amrex::Real* T,      ARLIM_P(Tlo),     ARLIM_P(Thi));

  void calc_gamma_pinv(const int* lo,  const int* hi,
                          amrex::Real* theta,        ARLIM_P(thlo),    ARLIM_P(thhi),
                          const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                          const amrex::Real* T,      ARLIM_P(Tlo),     ARLIM_P(Thi),
                          const amrex::Real* Pamb_in);
    
  void floor_spec(const int* lo,  const int* hi,
                      const amrex::Real* spec,   ARLIM_P(speclo),    ARLIM_P(spechi));

  void set_common (const amrex::Real* time, const int* iteration);

  void compute_ugradp (const amrex::Real* p, ARLIM_P(p_lo),  ARLIM_P(p_hi),
                           const amrex::Real* ugradp, ARLIM_P(up_lo),  ARLIM_P(up_hi),
                           const amrex::Real* umac, ARLIM_P(umac_lo),  ARLIM_P(umac_hi),
                           const amrex::Real* vmac, ARLIM_P(vmac_lo),  ARLIM_P(vmac_hi),
#if BL_SPACEDIM == 3
                           const amrex::Real* wmac, ARLIM_P(wmac_lo),  ARLIM_P(wmac_hi),
#endif
                           const int* lo, const int* hi, const amrex::Real* dx);

  int conservative_T_floor(const int* lo, const int* hi,
                                amrex::Real* fabdat, ARLIM_P(flo), ARLIM_P(fhi),
                                const amrex::Real* min_T_fine, const int* Tcomp, const int* Rcomp,
                                const int* first_comp, const int* last_spec, const int* RhoH,
                                const int* ratio, amrex::Real* tmp, const int* nt);
    
  void part_cnt_err(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                         const int* tagval, const int* clearval,
                         amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* lo, const int* hi, const int* nvar,
                         const int* domain_lo, const int* domain_hi,
                         const amrex::Real* dx, const amrex::Real* xlo,
                         const amrex::Real* prob_lo, const amrex::Real* time,
                         const int* level);

  void mcurve(const int* lo, const int* hi,
                   const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
                   amrex::Real*       K, ARLIM_P(Klo),  ARLIM_P(Khi),
                   amrex::Real*      nW, ARLIM_P(nwlo), ARLIM_P(nwhi),
                   const amrex::Real* delta);

  void smooth(const int* lo, const int* hi,
                   amrex::Real*       i, ARLIM_P(ilo), ARLIM_P(ihi),
                   const amrex::Real* o, ARLIM_P(olo), ARLIM_P(ohi));

  void grad_wbar(const int* lo, const int* hi,
                     const amrex::Real* Wbar, ARLIM_P(Wlo),    ARLIM_P(Whi),
                     const amrex::Real* rDe,  ARLIM_P(rDelo),  ARLIM_P(rDehi),
                     amrex::Real* flx,  ARLIM_P(flxlo),  ARLIM_P(flxhi),
                     const amrex::Real* a,    ARLIM_P(alo),    ARLIM_P(ahi),
                     const amrex::Real* delta, const int* dir,
                     const amrex::Real* mult, const amrex::Real* inc);
    
  void valgt_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                   const int* tagval, const int* clearval,
                   const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                   const int* lo, const int* hi, const int* nvar,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* dx, const amrex::Real* xlo,
                   const amrex::Real* prob_lo, const amrex::Real* time,
                   const int* level, const amrex::Real* value);

  void vallt_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                   const int* tagval, const int* clearval,
                   const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                   const int* lo, const int* hi, const int* nvar,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* dx, const amrex::Real* xlo,
                   const amrex::Real* prob_lo, const amrex::Real* time,
                   const int* level, const amrex::Real* value);

  void magvort_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                   const int* tagval, const int* clearval,
                   const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                   const int* lo, const int* hi, const int* nvar,
                   const int* domain_lo, const int* domain_hi,
                   const amrex::Real* dx, const amrex::Real* xlo,
                   const amrex::Real* prob_lo, const amrex::Real* time,
                   const int* level, const amrex::Real* value);

  void diffgt_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                    const int* tagval, const int* clearval,
                    const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                    const int* lo, const int* hi, const int* nvar,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* dx, const amrex::Real* xlo,
                    const amrex::Real* prob_lo, const amrex::Real* time,
                    const int* level, const amrex::Real* value);
  
  void box_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
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
   void dcma_error(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                    const int* tagval, const int* clearval,
                    const amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
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
