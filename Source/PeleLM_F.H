#ifndef _PeleLM_F_H_
#define _PeleLM_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_ACTIVECONTROL             activecontrol
#    define FORT_SET_SCAL_NUMB             setscnum
#    define FORT_GETPAMB                   getpamb
#    define FORT_GETCLOSEDCHAMBER          getclosedchamber
#    define FORT_GETDPDT                   getdpdt
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady
#    define FORT_DIVUFILL                  divufill
#    define FORT_DSDTFILL                  dsdtfill
#    define FORT_DQRADFILL                 dqradfill
#    define FORT_YDOTFILL                  ydotfill
#    define FORT_RHOYDOTFILL               rhoydotfill
#    define FORT_EST_DIVU_DT               estdivudt
#    define FORT_CHECK_DIVU_DT             checkdivudt
#    define FORT_INIT_VISC_DATA            initviscdata
#    define FORT_SET_HT_VISC_COMMON        setvisccm
#    define FORT_INIT_TYPVALS_COMMON       inittvcm
#    define FORT_VELVISC                   velvisc
#    define FORT_SPECTEMPVISC              spectempvisc
#    define FORT_FABMINMAX                 fabminmax
#    define FORT_REPAIR_FLUX               rparflux
#    define FORT_RECOMP_UPDATE             recompupdt
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv
#    define FORT_FLUXDIV                   flxdiv
#    define FORT_CALCDIVU                  calcdivu
#    define FORT_CALCGAMMAPINV             calcgammapinv
#    define FORT_FLOORSPEC                 floorspec
#    define FORT_ENTH_DIFF_TERMS           enthdiff
#    define FORT_SET_COMMON                setcommon
#    define FORT_COMPUTE_UGRADP            computeugradp
#    define FORT_HTDDRELAX                 htddrelax
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent
#    define FORT_RK_STEP_SELECTION         rkstepselection
#    define FORT_CONSERVATIVE_T_FLOOR      constfloor
#    define FORT_GETTYPICALVALS            gettypvals
#    define FORT_SETTYPICALVALS            settypvals
#    define FORT_PART_CNT_ERR              particlecnterr
#    define FORT_MCURVE                    mcurve
#    define FORT_SMOOTH                    smooth
#    define FORT_GRADWBAR                  gradwbar
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_ACTIVECONTROL             ACTIVECONTROL
#    define FORT_SET_SCAL_NUMB             SETSCNUM
#    define FORT_GETPAMB                   GETPAMB
#    define FORT_GETCLOSEDCHAMBER          GETCLOSEDCHAMBER
#    define FORT_GETDPDT                   GETDPDT
#    define FORT_COMPUTE_RHODGRADHDOTGRADY COMPUTERHODGRADHDOTGRADY
#    define FORT_DIVUFILL                  DIVUFILL
#    define FORT_DSDTFILL                  DSDTFILL
#    define FORT_DQRADFILL                 DQRADFILL
#    define FORT_YDOTFILL                  YDOTFILL
#    define FORT_RHOYDOTFILL               RHOYDOTFILL
#    define FORT_EST_DIVU_DT               ESTDIVUDT
#    define FORT_CHECK_DIVU_DT             CHECKDIVUDT
#    define FORT_INIT_VISC_DATA            INITVISCDATA
#    define FORT_SET_HT_VISC_COMMON        SETVISCCM
#    define FORT_INIT_TYPVALS_COMMON       INITTVCM
#    define FORT_VELVISC                   VELVISC
#    define FORT_SPECTEMPVISC              SPECTEMPVISC
#    define FORT_FABMINMAX                 FABMINMAX
#    define FORT_REPAIR_FLUX               RPARFLUX
#    define FORT_RECOMP_UPDATE             RECOMPUPDT
#    define FORT_INCRWEXTFLXDIV            INCRWEXTFLXDIV
#    define FORT_FLUXDIV                   FLXDIV
#    define FORT_CALCDIVU                  CALCDIVU
#    define FORT_CALCGAMMAPINV             CALCGAMMAPINV
#    define FORT_FLOORSPEC                 FLOORSPEC
#    define FORT_ENTH_DIFF_TERMS           ENTHDIFF
#    define FORT_SET_COMMON                SETCOMMON
#    define FORT_COMPUTE_UGRADP            COMPUTEUGRADP
#    define FORT_HTDDRELAX                 HTDDRELAX
#    define FORT_RK_MIXTURE_AVERAGED       RKMIXTUREAVERAGED
#    define FORT_RK_MULTICOMPONENT         RKMULTICOMPONENT
#    define FORT_RK_STEP_SELECTION         RKSTEPSELECTION
#    define FORT_CONSERVATIVE_T_FLOOR      CONSTFLOOR
#    define FORT_GETTYPICALVALS            GETTYPVALS
#    define FORT_SETTYPICALVALS            SETTYPVALS
#    define FORT_PART_CNT_ERR              PARTICLECNTERR
#    define FORT_MCURVE                    MCURVE
#    define FORT_SMOOTH                    SMOOTH
#    define FORT_GRADWBAR                  GRADWBAR
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_ACTIVECONTROL             activecontrol
#    define FORT_SET_SCAL_NUMB             setscnum
#    define FORT_GETPAMB                   getpamb
#    define FORT_GETCLOSEDCHAMBER          getclosedchamber
#    define FORT_GETDPDT                   getdpdt
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady
#    define FORT_DIVUFILL                  divufill
#    define FORT_DSDTFILL                  dsdtfill
#    define FORT_DQRADFILL                 dqradfill
#    define FORT_YDOTFILL                  ydotfill
#    define FORT_RHOYDOTFILL               rhoydotfill
#    define FORT_EST_DIVU_DT               estdivudt
#    define FORT_CHECK_DIVU_DT             checkdivudt
#    define FORT_INIT_VISC_DATA            initviscdata
#    define FORT_SET_HT_VISC_COMMON        setvisccm
#    define FORT_INIT_TYPVALS_COMMON       inittvcm
#    define FORT_VELVISC                   velvisc
#    define FORT_SPECTEMPVISC              spectempvisc
#    define FORT_FABMINMAX                 fabminmax
#    define FORT_REPAIR_FLUX               rparflux
#    define FORT_RECOMP_UPDATE             recompupdt
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv
#    define FORT_FLUXDIV                   flxdiv
#    define FORT_CALCDIVU                  calcdivu
#    define FORT_CALCGAMMAPINV             calcgammapinv
#    define FORT_FLOORSPEC                 floorspec
#    define FORT_ENTH_DIFF_TERMS           enthdiff
#    define FORT_SET_COMMON                setcommon
#    define FORT_COMPUTE_UGRADP            computeugradp
#    define FORT_HTDDRELAX                 htddrelax
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent
#    define FORT_RK_STEP_SELECTION         rkstepselection
#    define FORT_CONSERVATIVE_T_FLOOR      constfloor
#    define FORT_GETTYPICALVALS            gettypvals
#    define FORT_SETTYPICALVALS            settypvals
#    define FORT_PART_CNT_ERR              particlecnterr
#    define FORT_MCURVE                    mcurve
#    define FORT_SMOOTH                    smooth
#    define FORT_GRADWBAR                  gradwbar
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_ACTIVECONTROL             activecontrol_
#    define FORT_SET_SCAL_NUMB             setscnum_
#    define FORT_GETPAMB                   getpamb_
#    define FORT_GETCLOSEDCHAMBER          getclosedchamber_
#    define FORT_GETDPDT                   getdpdt_
#    define FORT_COMPUTE_RHODGRADHDOTGRADY computerhodgradhdotgrady_
#    define FORT_DIVUFILL                  divufill_
#    define FORT_DSDTFILL                  dsdtfill_
#    define FORT_DQRADFILL                 dqradfill_
#    define FORT_YDOTFILL                  ydotfill_
#    define FORT_RHOYDOTFILL               rhoydotfill_
#    define FORT_EST_DIVU_DT               estdivudt_
#    define FORT_CHECK_DIVU_DT             checkdivudt_
#    define FORT_INIT_VISC_DATA            initviscdata_
#    define FORT_SET_HT_VISC_COMMON        setvisccm_
#    define FORT_INIT_TYPVALS_COMMON       inittvcm_
#    define FORT_VELVISC                   velvisc_
#    define FORT_SPECTEMPVISC              spectempvisc_
#    define FORT_FABMINMAX                 fabminmax_
#    define FORT_REPAIR_FLUX               rparflux_
#    define FORT_RECOMP_UPDATE             recompupdt_
#    define FORT_INCRWEXTFLXDIV            incrwextflxdiv_
#    define FORT_FLUXDIV                   flxdiv_
#    define FORT_CALCDIVU                  calcdivu_
#    define FORT_CALCGAMMAPINV             calcgammapinv_
#    define FORT_FLOORSPEC                 floorspec_
#    define FORT_ENTH_DIFF_TERMS           enthdiff_
#    define FORT_SET_COMMON                setcommon_
#    define FORT_COMPUTE_UGRADP            computeugradp_
#    define FORT_HTDDRELAX                 htddrelax_
#    define FORT_RK_MIXTURE_AVERAGED       rkmixtureaveraged_
#    define FORT_RK_MULTICOMPONENT         rkmulticomponent_
#    define FORT_RK_STEP_SELECTION         rkstepselection_
#    define FORT_CONSERVATIVE_T_FLOOR      constfloor_
#    define FORT_GETTYPICALVALS            gettypvals_
#    define FORT_SETTYPICALVALS            settypvals_
#    define FORT_PART_CNT_ERR              particlecnterr_
#    define FORT_MCURVE                    mcurve_
#    define FORT_SMOOTH                    smooth_
#    define FORT_GRADWBAR                  gradwbar_
#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

  void FORT_ACTIVECONTROL(const amrex::Real* fuelmass,
                          const amrex::Real* time,
                          const amrex::Real* crsedt,
                          const int*  myproc,
                          const int*  step,
                          const int*  restart,
                          const int*  usetemp);

  void FORT_SET_SCAL_NUMB(const int* Density,
                          const int* Temp,
                          const int* Trac,
                          const int* RhoH,
                          const int* FirstSpec,
                          const int* LastSpec);
    
  void FORT_GETPAMB(amrex::Real* pamb);

  void FORT_GETCLOSEDCHAMBER(int* closed_chamber);

  void FORT_GETDPDT(amrex::Real* dpdt_factor);

  void FORT_COMPUTE_RHODGRADHDOTGRADY(const amrex::Real* dx,
                                      const int* lo, const int* hi,
                                      ARLIM_P(speclo),ARLIM_P(spechi), const amrex::Real* specdat,
                                      ARLIM_P(hlo),ARLIM_P(hhi), const amrex::Real* hdat,
                                      ARLIM_P(betaxlo),ARLIM_P(betaxhi), const amrex::Real* betax,
                                      ARLIM_P(betaylo),ARLIM_P(betayhi), const amrex::Real* betay,
#if (BL_SPACEDIM == 3) 
                                      ARLIM_P(betazlo),ARLIM_P(betazhi), const amrex::Real* betaz,
#endif            
                                      ARLIM_P(ddrhslo),ARLIM_P(ddrhshi), amrex::Real* ddrhsdat);

  void FORT_DIVUFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void FORT_DSDTFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void FORT_DQRADFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                       const int* dom_lo, const int* dom_hi,
                       const amrex::Real* dx, const amrex::Real* grd_lo,
                       const amrex::Real* time, const int* bc);

  void FORT_YDOTFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                      const int* dom_lo, const int* dom_hi,
                      const amrex::Real* dx, const amrex::Real* grd_lo,
                      const amrex::Real* time, const int* bc);

  void FORT_RHOYDOTFILL (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* dom_lo, const int* dom_hi,
                         const amrex::Real* dx, const amrex::Real* grd_lo,
                         const amrex::Real* time, const int* bc);

  void FORT_EST_DIVU_DT(const int& flag,  const amrex::Real* divu_dt_factor,
                        const amrex::Real* dx, 
                        const amrex::Real* sdat,ARLIM_P(slo),ARLIM_P(shi),
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
                        amrex::Real* dt, const amrex::Real* min_rho_divu_ceiling);

  void FORT_CHECK_DIVU_DT(const int& flag, const amrex::Real* divu_dt_factor,
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

  void FORT_INIT_VISC_DATA();
  void FORT_SET_HT_VISC_COMMON(const int* muIsVar,     const amrex::Real* muVal,
                               const int* lambdaIsVar, const amrex::Real* lambdaVal,
                               const int* rhoDIsVar,   const amrex::Real* rhoDVal,
                               const amrex::Real* Pr, const amrex::Real* Sc, const int* unityLe);

  void FORT_INIT_TYPVALS_COMMON();
      
  void FORT_VELVISC(const int* lo, const int* hi,
                    ARLIM_P(tlo), ARLIM_P(thi), const amrex::Real* tempdat,
                    ARLIM_P(slo), ARLIM_P(shi), const amrex::Real* specdat,
                    ARLIM_P(blo), ARLIM_P(bhi), amrex::Real* betadat);

  void FORT_SPECTEMPVISC(const int* lo, const int* hi,
                         ARLIM_P(tlo), ARLIM_P(thi), const amrex::Real* tempdat,
                         ARLIM_P(slo), ARLIM_P(shi), const amrex::Real* specdat,
                         ARLIM_P(blo), ARLIM_P(bhi), amrex::Real* betadat,
                         const int* nCompD, const amrex::Real* P1ATM_MKS,
                         const int* do_temp, const int* do_VelVisc, 
                         const amrex::Real* Pamb_in);

  void FORT_FABMINMAX(const int* lo, const int* hi,
                      amrex::Real* fabdat, ARLIM_P(flo), ARLIM_P(fhi),
                      const amrex::Real* min, const amrex::Real* max, const int* nc);

  void FORT_REPAIR_FLUX(const int* lo,  const int* hi,
                        const int* dlo, const int* dhi,
                        amrex::Real* flux,    ARLIM_P(flo),  ARLIM_P(fhi),
                        const amrex::Real* RY,ARLIM_P(RYlo), ARLIM_P(RYhi),
                        const int* dir, const int* bc);

  void FORT_ENTH_DIFF_TERMS(const int* lo, const int* hi, 
                            const int* dlo, const int* dhi, const amrex::Real* dx,
                            const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
                            const amrex::Real* RhoY, ARLIM_P(rYlo),  ARLIM_P(rYhi),

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
                            amrex::Real* fh,      ARLIM_P(fhlo), ARLIM_P(fhhi),
                            const int* Tbc);

  void FORT_INCRWEXTFLXDIV(const int* lo, const int* hi,
                           const amrex::Real* xflux, ARLIM_P(fxlo),  ARLIM_P(fxhi),
                           const amrex::Real* yflux, ARLIM_P(fylo),  ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                           const amrex::Real* zflux, ARLIM_P(fzlo),  ARLIM_P(fzhi),
#endif
                           const amrex::Real* sto,     ARLIM_P(stolo),    ARLIM_P(stohi),
                           amrex::Real*       stn,     ARLIM_P(stnlo),    ARLIM_P(stnhi),
                           const amrex::Real* vol,     ARLIM_P(vollo),    ARLIM_P(volhi),
                           const int* nc, const amrex::Real* dt);

  void FORT_FLUXDIV(const int* lo, const int* hi,
                    const amrex::Real* updat, ARLIM_P(uplo), ARLIM_P(uphi),
                    const amrex::Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
                    const amrex::Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                    const amrex::Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
                    const amrex::Real* vol,   ARLIM_P(vlo),  ARLIM_P(vhi),
                    const int* nc, const amrex::Real* scale);

  void FORT_CALCDIVU(const int* lo,  const int* hi,
                     amrex::Real* divu,         ARLIM_P(dulo),    ARLIM_P(duhi),
                     const amrex::Real* rYdot,  ARLIM_P(rYdotlo), ARLIM_P(rYdothi),
                     const amrex::Real* vtY,    ARLIM_P(vtYlo),   ARLIM_P(vtYhi),
                     const amrex::Real* vtT,    ARLIM_P(vtTlo),   ARLIM_P(vtThi),
                     const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                     const amrex::Real* T,      ARLIM_P(Tlo),     ARLIM_P(Thi));

  void FORT_CALCGAMMAPINV(const int* lo,  const int* hi,
                          amrex::Real* theta,        ARLIM_P(thlo),    ARLIM_P(thhi),
                          const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                          const amrex::Real* T,      ARLIM_P(Tlo),     ARLIM_P(Thi),
                          const amrex::Real* Pamb_in);
    
  void FORT_FLOORSPEC(const int* lo,  const int* hi,
                      const amrex::Real* spec,   ARLIM_P(speclo),    ARLIM_P(spechi));

  void FORT_SET_COMMON (const amrex::Real* time, const int* iteration);

  void FORT_COMPUTE_UGRADP(const amrex::Real* p, ARLIM_P(p_lo),  ARLIM_P(p_hi),
                           const amrex::Real* ugradp, ARLIM_P(up_lo),  ARLIM_P(up_hi),
                           const amrex::Real* umac, ARLIM_P(umac_lo),  ARLIM_P(umac_hi),
                           const amrex::Real* vmac, ARLIM_P(vmac_lo),  ARLIM_P(vmac_hi),
#if BL_SPACEDIM == 3
                           const amrex::Real* wmac, ARLIM_P(wmac_lo),  ARLIM_P(wmac_hi),
#endif
                           const int* lo, const int* hi, const amrex::Real* dx);

  void FORT_HTDDRELAX(const int* lo, const int* hi,
                      amrex::Real* s, ARLIM_P(slo),  ARLIM_P(shi),
                      const int& sCompY, const int& sCompT,
                      const int& sCompH, const int& sCompR,
                      const amrex::Real* L, ARLIM_P(Llo),  ARLIM_P(Lhi),
                      const amrex::Real* a, ARLIM_P(alo),  ARLIM_P(ahi),
                      const amrex::Real* r, ARLIM_P(rlo),  ARLIM_P(rhi),
                      const amrex::Real& thetaDt, const amrex::Real* fac, amrex::Real* maxRes, amrex::Real* maxCor,
                      const int& for_T0_H1, const int& res_only, const amrex::Real& mult);

  void FORT_GETTYPICALVALS(amrex::Real* vals, const int* nVals);

  void FORT_SETTYPICALVALS(const amrex::Real* vals, const int* nVals);

  void FORT_RK_MIXTURE_AVERAGED (
    const int* domain_lo, const int* domain_hi,
    const int* lo, const int* hi,
    amrex::Real* areax, ARLIM_P(areax_lo), ARLIM_P(areax_hi),
    amrex::Real* areay, ARLIM_P(areay_lo), ARLIM_P(areay_hi),
    const int* bc,
    const amrex::Real* dx,
    int* index_of_firstY,
    int* index_of_lastY,
    int* index_of_rho,
    int* index_of_rhoH,
    int* index_of_T,
    amrex::Real* maximum_error,
    int* maximum_iterations,
    int* ncomps,
    int* nspcecies,
    const amrex::Real* scale_factor,
    amrex::Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi),
    amrex::Real* volume, ARLIM_P(volume_lo), ARLIM_P(volume_hi),
    amrex::Real* div_of_flux_for_H, ARLIM_P(div_of_flux_for_H_lo), ARLIM_P(div_of_flux_for_H_hi),
    amrex::Real* div_of_flux_for_Y, ARLIM_P(div_of_flux_for_Y_lo), ARLIM_P(div_of_flux_for_Y_hi),
    amrex::Real* xflux_for_H, ARLIM_P(xflux_for_H_lo), ARLIM_P(xflux_for_H_hi),
    amrex::Real* xflux_for_Y, ARLIM_P(xflux_for_Y_lo), ARLIM_P(xflux_for_Y_hi),
    amrex::Real* yflux_for_Y, ARLIM_P(yflux_for_Y_lo), ARLIM_P(yflux_for_Y_hi),
    amrex::Real* yflux_for_H, ARLIM_P(yflux_for_H_lo), ARLIM_P(yflux_for_H_hi)
    );

  void FORT_RK_MULTICOMPONENT (
    const int* domain_lo, const int* domain_hi,
    const int* lo, const int* hi,
    amrex::Real* areax, ARLIM_P(areax_lo), ARLIM_P(areax_hi),
    amrex::Real* areay, ARLIM_P(areay_lo), ARLIM_P(areay_hi),
    const int* bc,
    const amrex::Real* dx,
    int* index_of_firstY,
    int* index_of_lastY,
    int* index_of_rho,
    int* index_of_rhoH,
    int* index_of_T,
    amrex::Real* maximum_error,
    int* maximum_iterations,
    int* ncomps,
    int* nspcecies,
    const amrex::Real* scale_factor,
    amrex::Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi),
    amrex::Real* volume, ARLIM_P(volume_lo), ARLIM_P(volume_hi),
    amrex::Real* div_of_flux_for_H, ARLIM_P(div_of_flux_for_H_lo), ARLIM_P(div_of_flux_for_H_hi),
    amrex::Real* div_of_flux_for_Y, ARLIM_P(div_of_flux_for_Y_lo), ARLIM_P(div_of_flux_for_Y_hi),
    amrex::Real* xflux_for_H, ARLIM_P(xflux_for_H_lo), ARLIM_P(xflux_for_H_hi),
    amrex::Real* xflux_for_Y, ARLIM_P(xflux_for_Y_lo), ARLIM_P(xflux_for_Y_hi),
    amrex::Real* yflux_for_Y, ARLIM_P(yflux_for_Y_lo), ARLIM_P(yflux_for_Y_hi),
    amrex::Real* yflux_for_H, ARLIM_P(yflux_for_H_lo), ARLIM_P(yflux_for_H_hi)
    );

  void FORT_RK_STEP_SELECTION (
    const int* lo, const int* hi,
    const amrex::Real* dx,
    int* index_of_firstY,
    int* index_of_lastY,
    int* index_of_rho,
    int* index_of_rhoH,
    int* index_of_T,
    amrex::Real* maximum_error,
    int* maximum_iterations,
    int* ncomps,
    int* nspcecies,
    amrex::Real* smallest_time_step,
    amrex::Real* state, ARLIM_P(state_lo), ARLIM_P(state_hi)
    );

  int FORT_CONSERVATIVE_T_FLOOR(const int* lo, const int* hi,
                                amrex::Real* fabdat, ARLIM_P(flo), ARLIM_P(fhi),
                                const amrex::Real* min_T_fine, const int* Tcomp, const int* Rcomp,
                                const int* first_comp, const int* last_spec, const int* RhoH,
                                const int* ratio, amrex::Real* tmp, const int* nt);
    
  void FORT_PART_CNT_ERR(int* tag, ARLIM_P(tlo), ARLIM_P(thi),
                         const int* tagval, const int* clearval,
                         amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                         const int* lo, const int* hi, const int* nvar,
                         const int* domain_lo, const int* domain_hi,
                         const amrex::Real* dx, const amrex::Real* xlo,
                         const amrex::Real* prob_lo, const amrex::Real* time,
                         const int* level);

  void FORT_RECOMP_UPDATE(const int* lo, const int* hi,
                          const amrex::Real* updat, ARLIM_P(uplo), ARLIM_P(uphi),
                          const amrex::Real* xflux, ARLIM_P(fxlo), ARLIM_P(fxhi),
                          const amrex::Real* yflux, ARLIM_P(fylo), ARLIM_P(fyhi),
#if BL_SPACEDIM == 3
                          const amrex::Real* zflux, ARLIM_P(fzlo), ARLIM_P(fzhi),
#endif
                          const amrex::Real* vol,   ARLIM_P(vlo),  ARLIM_P(vhi),
                          const int* nc);

  void FORT_MCURVE(const int* lo, const int* hi,
                   const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
                   amrex::Real*       K, ARLIM_P(Klo),  ARLIM_P(Khi),
                   amrex::Real*      nW, ARLIM_P(nwlo), ARLIM_P(nwhi),
                   const amrex::Real* delta);

  void FORT_SMOOTH(const int* lo, const int* hi,
                   amrex::Real*       i, ARLIM_P(ilo), ARLIM_P(ihi),
                   const amrex::Real* o, ARLIM_P(olo), ARLIM_P(ohi));

  void FORT_GRADWBAR(const int* lo, const int* hi,
                     const amrex::Real* Wbar, ARLIM_P(Wlo),    ARLIM_P(Whi),
                     const amrex::Real* rDe,  ARLIM_P(rDelo),  ARLIM_P(rDehi),
                     amrex::Real* flx,  ARLIM_P(flxlo),  ARLIM_P(flxhi),
                     const amrex::Real* a,    ARLIM_P(alo),    ARLIM_P(ahi),
                     const amrex::Real* delta, const int* dir,
                     const amrex::Real* mult, const amrex::Real* inc);
    
#ifdef __cplusplus
}
#endif
#endif

#endif /*_PeleLM_F_H_*/
