#ifndef _ChemDriver_F_H_
#define _ChemDriver_F_H_

#ifdef BL_LANG_FORT
#    define CD_XTY_PT            cdxtypt
#    define CD_YTX_PT            cdytxpt
#    define CD_MWT               cdmwt
#    define FORT_USINGEG         dusingeg
#    define FORT_USINGMC         dusingmc
#    define FORT_MCLINKNAME      dmclink
#    define FORT_SETVERBOSEVODE  dverbose
#    define FORT_SETVODETOLS     dvodetols
#    define FORT_SETVODESUBCYC   dmxsubcy
#    define FORT_SETSPECSCALY    dsetscal
#    define FORT_INITCHEM        dinitchem
#    define FORT_FINALIZECHEM    dfinalchem
#    define FORT_GETCKMAXNAMELEN dckmxnml
#    define FORT_GETCKNUMSPEC    dcknspec
#    define FORT_GETCKNUMELT     dcknelt
#    define FORT_GETCKSPECNAME   dckspnam
#    define FORT_GETCKELTNAME    dckeltnam
#    define FORT_GETCKMWT        dckmwt
#    define FORT_GETCKAWT        dckawt
#    define FORT_GETCKDIMPARAMS  dckdims
#    define FORT_MASSTOMOLE      dmstomol
#    define FORT_FRrateXTP       dfrrtxtp
#    define FORT_HTRLS           dhtrls
#    define FORT_RRATEY          drratey
#ifdef LMC_SDC
#    define FORT_RRATERHOY       drraterhoy
#endif
#    define FORT_RRATEC          drratec
#    define FORT_MOLETOMASS      dmoltoms
#    define FORT_MASSTP_TO_CONC  dmstptocon
#    define FORT_MASSR_TO_CONC   dmsrtocon
#    define FORT_CONC_TO_MOLE    dcontomol
#    define FORT_CONPSOLV        dcpsolv
#    define FORT_GET_REACTION_MAP dgrm
#ifdef LMC_SDC
#    define FORT_CONPSOLV_SDC    dcpsolvsdc
#endif
#    define FORT_BETA_WBAR       dbetawbar
#    define FORT_MIXAVG_RHODIFF_TEMP dmavgdiftemp
#    define FORT_MIX_SHEAR_VISC  dmixshvis
#    define FORT_RHOfromPTY      drhopty
#    define FORT_RHOfromPvTY     drhopvty
#    define FORT_PfromRTY        dprhoty
#    define FORT_TfromPRY        dtprhoy
#    define FORT_CPMIXfromTY     dcpmty
#    define FORT_CVMIXfromTY     dcvmty
#    define FORT_HMIXfromTY      dhmty
#    define FORT_MWMIXfromY      dmwmy
#    define FORT_CPfromT         dcpt
#    define FORT_HfromT          dht
#    define FORT_TfromHY         dthy
#    define FORT_TfromHYpt       dthypt
#    define FORT_GETCKNUMREAC    dcknreac
#    define FORT_MOLPROD         dmolprod
#    define FORT_GETELTMOLES     dgetemol
#    define FORT_RUNIV           dgetr
#    define FORT_P1ATMMKS        dgetp1atm
#    define FORT_OTrad_TDF       dgetrad
#    define FORT_FINDLHS         dfindlhs
#    define FORT_FINDRHS         dfindrhs
#    define FORT_CKINU           dckinu
#    define FORT_CKSYMR          dcksymr
#    define FORT_CKELTXINSPY     deltxinspy
#    define FORT_NORMMASS        dnormmass
#    define FORT_SETTMINTRANS    dsetttrans
#    define FORT_SETNU           dsetnu
#    define FORT_GET_REACTION_MAP dgrm
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define CD_XTY_PT            CDXTYPT
#    define CD_YTX_PT            CDYTXPT
#    define CD_MWT               CDMWT
#    define FORT_USINGEG         DUSINGEG
#    define FORT_USINGMC         DUSINGMC
#    define FORT_MCLINKNAME      DMCLINK
#    define FORT_SETVERBOSEVODE  DVERBOSE
#    define FORT_SETVVODETOLS    DVODETOLS
#    define FORT_SETVODESUBCYC   DMXSUBCY
#    define FORT_SETSPECSCALY    DSETSCAL
#    define FORT_INITCHEM        DINITCHEM
#    define FORT_FINALIZECHEM    DFINALCHEM
#    define FORT_GETCKMAXNAMELEN DCKMXNML
#    define FORT_GETCKNUMSPEC    DCKNSPEC
#    define FORT_GETCKNUMELT     DCKNELT
#    define FORT_GETCKSPECNAME   DCKSPNAM
#    define FORT_GETCKELTNAME    DCKELTNAM
#    define FORT_GETCKAWT        DCKAWT
#    define FORT_GETCKMWT        DCKMWT
#    define FORT_GETCKDIMPARAMS  DCKDIMS
#    define FORT_MASSTOMOLE      DMSTOMOL
#    define FORT_FRrateXTP       DFRRTXTP
#    define FORT_HTRLS           DHTRLS
#    define FORT_RRATEY          DRRATEY
#ifdef LMC_SDC
#    define FORT_RRATERHOY       DRRATERHOY
#endif
#    define FORT_RRATEC          DRRATEC
#    define FORT_MOLETOMASS      DMOLTOMS
#    define FORT_MASSTP_TO_CONC  DMSTPTOCON
#    define FORT_MASSR_TO_CONC   DMSRTOCON
#    define FORT_CONC_TO_MOLE    DCONTOMOL
#    define FORT_CONPSOLV        DCPSOLV
#ifdef LMC_SDC
#    define FORT_CONPSOLV_SDC    DCPSOLVSDC
#endif
#    define FORT_BETA_WBAR       DBETAWBAR
#    define FORT_MIXAVG_RHODIFF_TEMP DMAVGDIFTEMP
#    define FORT_MIX_SHEAR_VISC  DMIXSHVIS
#    define FORT_RHOfromPTY      DRHOPTY
#    define FORT_RHOfromPvTY     DRHOPVTY
#    define FORT_PfromRTY        DPRHOTY
#    define FORT_TfromPRY        DTPRHOY
#    define FORT_CPMIXfromTY     DCPMTY
#    define FORT_CVMIXfromTY     DCVMTY
#    define FORT_HMIXfromTY      DHMTY
#    define FORT_MWMIXfromY      DMWMY
#    define FORT_CPfromT         DCPT
#    define FORT_HfromT          DHT
#    define FORT_TfromHY         DTHY
#    define FORT_TfromHYpt       DTHYPT
#    define FORT_GETCKNUMREAC    DCKNREAC
#    define FORT_MOLPROD         DMOLPROD
#    define FORT_GETELTMOLES     DGETEMOL
#    define FORT_RUNIV           DGETR
#    define FORT_P1ATMMKS        DGETP1ATM
#    define FORT_OTrad_TDF       DGETRAD
#    define FORT_FINDLHS         DFINDLHS
#    define FORT_FINDRHS         DFINDRHS
#    define FORT_CKINU           DCKINU
#    define FORT_CKSYMR          DCKSYMR
#    define FORT_CKELTXINSPY     DELTXINSPY
#    define FORT_NORMMASS        DNORMMASS
#    define FORT_SETTMINTRANS    DSETTTRANS
#    define FORT_SETNU           DSETNU
#    define FORT_GET_REACTION_MAP DGRM
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define CD_XTY_PT            cdxtypt
#    define CD_YTX_PT            cdytxpt
#    define CD_MWT               cdmwt
#    define FORT_USINGEG         dusingeg
#    define FORT_USINGMC         dusingmc
#    define FORT_MCLINKNAME      dmclink
#    define FORT_SETVERBOSEVODE  dverbose
#    define FORT_SETVODETOLS     dvodetols
#    define FORT_SETVODESUBCYC   dmxsubcy
#    define FORT_SETSPECSCALY    dsetscal
#    define FORT_INITCHEM        dinitchem
#    define FORT_FINALIZECHEM    dfinalchem
#    define FORT_GETCKMAXNAMELEN dckmxnml
#    define FORT_GETCKNUMSPEC    dcknspec
#    define FORT_GETCKNUMELT     dcknelt
#    define FORT_GETCKSPECNAME   dckspnam
#    define FORT_GETCKELTNAME    dckeltnam
#    define FORT_GETCKMWT        dckmwt
#    define FORT_GETCKAWT        dckawt
#    define FORT_GETCKDIMPARAMS  dckdims
#    define FORT_MASSTOMOLE      dmstomol
#    define FORT_MASSTP_TO_CONC  dmstptocon
#    define FORT_MASSR_TO_CONC   dmsrtocon
#    define FORT_CONC_TO_MOLE    dcontomol
#    define FORT_FRrateXTP       dfrrtxtp
#    define FORT_HTRLS           dhtrls
#    define FORT_RRATEY          drratey
#ifdef LMC_SDC
#    define FORT_RRATERHOY       drraterhoy
#endif
#    define FORT_RRATEC          drratec
#    define FORT_MOLETOMASS      dmoltoms
#    define FORT_CONPSOLV        dcpsolv
#ifdef LMC_SDC
#    define FORT_CONPSOLV_SDC    dcpsolvsdc
#endif
#    define FORT_BETA_WBAR       dbetawbar
#    define FORT_MIXAVG_RHODIFF_TEMP dmavgdiftemp
#    define FORT_MIX_SHEAR_VISC  dmixshvis
#    define FORT_RHOfromPTY      drhopty
#    define FORT_RHOfromPvTY     drhopvty
#    define FORT_PfromRTY        dprhoty
#    define FORT_TfromPRY        dtprhoy
#    define FORT_CPMIXfromTY     dcpmty
#    define FORT_CVMIXfromTY     dcvmty
#    define FORT_HMIXfromTY      dhmty
#    define FORT_MWMIXfromY      dmwmy
#    define FORT_CPfromT         dcpt
#    define FORT_HfromT          dht
#    define FORT_TfromHY         dthy
#    define FORT_TfromHYpt       dthypt
#    define FORT_GETCKNUMREAC    dcknreac
#    define FORT_MOLPROD         dmolprod
#    define FORT_GETELTMOLES     dgetemol
#    define FORT_RUNIV           dgetr
#    define FORT_P1ATMMKS        dgetp1atm
#    define FORT_OTrad_TDF       dgetrad
#    define FORT_FINDLHS         dfindlhs
#    define FORT_FINDRHS         dfindrhs
#    define FORT_CKINU           dckinu
#    define FORT_CKSYMR          dcksymr
#    define FORT_CKELTXINSPY     deltxinspy
#    define FORT_NORMMASS        dnormmass
#    define FORT_SETTMINTRANS    dsetttrans
#    define FORT_SETNU           dsetnu
#    define FORT_GET_REACTION_MAP dgrm
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define CD_XTY_PT            cdxtypt_
#    define CD_YTX_PT            cdytxpt_
#    define CD_MWT               cdmwt_
#    define FORT_USINGEG         dusingeg_
#    define FORT_USINGMC         dusingmc_
#    define FORT_MCLINKNAME      dmclink_
#    define FORT_SETVERBOSEVODE  dverbose_
#    define FORT_SETVODETOLS     dvodetols_
#    define FORT_SETVODESUBCYC   dmxsubcy_
#    define FORT_SETSPECSCALY    dsetscal_
#    define FORT_INITCHEM        dinitchem_
#    define FORT_FINALIZECHEM    dfinalchem_
#    define FORT_GETCKMAXNAMELEN dckmxnml_
#    define FORT_GETCKNUMSPEC    dcknspec_
#    define FORT_GETCKNUMELT     dcknelt_
#    define FORT_GETCKSPECNAME   dckspnam_
#    define FORT_GETCKELTNAME    dckeltnam_
#    define FORT_GETCKMWT        dckmwt_
#    define FORT_GETCKAWT        dckawt_
#    define FORT_GETCKDIMPARAMS  dckdims_
#    define FORT_MASSTOMOLE      dmstomol_
#    define FORT_MASSTP_TO_CONC  dmstptocon_
#    define FORT_MASSR_TO_CONC   dmsrtocon_
#    define FORT_CONC_TO_MOLE    dcontomol_
#    define FORT_FRrateXTP       dfrrtxtp_
#    define FORT_HTRLS           dhtrls_
#    define FORT_RRATEY          drratey_
#ifdef LMC_SDC
#    define FORT_RRATERHOY       drraterhoy_
#endif
#    define FORT_RRATEC          drratec_
#    define FORT_MOLETOMASS      dmoltoms_
#    define FORT_CONPSOLV        dcpsolv_
#ifdef LMC_SDC
#    define FORT_CONPSOLV_SDC    dcpsolvsdc_
#endif
#    define FORT_BETA_WBAR       dbetawbar_
#    define FORT_MIXAVG_RHODIFF_TEMP dmavgdiftemp_
#    define FORT_MIX_SHEAR_VISC  dmixshvis_
#    define FORT_RHOfromPTY      drhopty_
#    define FORT_RHOfromPvTY     drhopvty_
#    define FORT_PfromRTY        dprhoty_
#    define FORT_TfromPRY        dtprhoy_
#    define FORT_CPMIXfromTY     dcpmty_
#    define FORT_CVMIXfromTY     dcvmty_
#    define FORT_HMIXfromTY      dhmty_
#    define FORT_MWMIXfromY      dmwmy_
#    define FORT_CPfromT         dcpt_
#    define FORT_HfromT          dht_
#    define FORT_TfromHY         dthy_
#    define FORT_TfromHYpt       dthypt_
#    define FORT_GETCKNUMREAC    dcknreac_
#    define FORT_MOLPROD         dmolprod_
#    define FORT_GETELTMOLES     dgetemol_
#    define FORT_RUNIV           dgetr_
#    define FORT_P1ATMMKS        dgetp1atm_
#    define FORT_SETTRLIB        dsettrlb_
#    define FORT_GET_TR_ID       dgettrid_
#    define FORT_GET_EG_ID       dgetegid_
#    define FORT_OTrad_TDF       dgetrad_
#    define FORT_FINDLHS         dfindlhs_
#    define FORT_FINDRHS         dfindrhs_
#    define FORT_CKINU           dckinu_
#    define FORT_CKSYMR          dcksymr_
#    define FORT_CKELTXINSPY     deltxinspy_
#    define FORT_NORMMASS        dnormmass_
#    define FORT_SETTMINTRANS    dsetttrans_
#    define FORT_SETNU           dsetnu_
#    define FORT_GET_REACTION_MAP dgrm_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_SETTMINTRANS(const Real* Tmin);
    void FORT_SETVERBOSEVODE();
    void FORT_SETVODETOLS(const Real* rtol, const Real* atol, const int* itol);
    void FORT_SETVODESUBCYC(const int* maxcyc);
    void FORT_SETSPECSCALY(const int* name, const int* length);
    void FORT_INITCHEM();
    void FORT_FINALIZECHEM();
    int FORT_GETCKMAXNAMELEN();
    int FORT_GETCKNUMSPEC();
    int FORT_GETCKNUMELT();
    int FORT_GETCKSPECNAME(const int* str_len, const int* coded);
    int FORT_GETCKELTNAME(const int* str_len, const int* coded);
    int FORT_CKELTXINSPY(const int* eltID, const int* spID);
    void FORT_GETCKMWT(Real* mwt);
    void FORT_GETCKAWT(Real* awt);
    void FORT_FINDLHS(int* reactions, int* listLen, const int* idx);
    void FORT_FINDRHS(int* reactions, int* listLen, const int* idx);
    void FORT_CKINU(int* Nids,int* KI,const int* lenKI,int* NU,const int* lenNU,
                    const int* rxnID, const int* nuAll);
    int FORT_CKSYMR(const int* fortReacIdx, const int* coded);

    void FORT_GETCKDIMPARAMS(const int* maxreac,  const int* maxspec,
                             const int* maxelts,  const int* maxord,
                             const int* maxthrdb, const int* maxtp,
                             const int* maxsp,    const int* maxspnml);

    void FORT_FRrateXTP(const int* lo, const int* hi,
                        const Real* X, ARLIM_P(Xlo),    ARLIM_P(Xhi),
                        const Real* T, ARLIM_P(Tlo),    ARLIM_P(Thi),
                        Real*    FwdK, ARLIM_P(FwdKlo), ARLIM_P(FwdKhi),
                        Real*    RevK, ARLIM_P(RevKlo), ARLIM_P(RevKhi),
                        const Real* Patm, const int* rxns, const int* Nrxns);

    void FORT_HTRLS(const int* lo, const int* hi,
                    const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                    const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                    Real*       Q, ARLIM_P(Qlo), ARLIM_P(Qhi), const Real* Patm);

    void FORT_RRATEY(const int* lo, const int* hi,
                     const Real* Y, ARLIM_P(Ylo),    ARLIM_P(Yhi),
                     const Real* T, ARLIM_P(Tlo),    ARLIM_P(Thi),
                     Real*    Ydot, ARLIM_P(Ydotlo), ARLIM_P(Ydothi),
                     Real* Patm);    

#ifdef LMC_SDC
    void FORT_RRATERHOY(const int* lo, const int* hi,
                        const Real* RhoY, ARLIM_P(RhoYlo),    ARLIM_P(RhoYhi),
                        const Real* RhoH, ARLIM_P(RhoHlo),    ARLIM_P(RhoHhi),
                        const Real* T,    ARLIM_P(Tlo),       ARLIM_P(Thi),
                        Real*    RhoYdot, ARLIM_P(RhoYdotlo), ARLIM_P(RhoYdothi) );
#endif

    void FORT_RRATEC(const int* lo, const int* hi,
                     const Real* C, ARLIM_P(Clo),    ARLIM_P(Chi),
                     const Real* T, ARLIM_P(Tlo),    ARLIM_P(Thi),
                     Real*    Cdot, ARLIM_P(Cdotlo), ARLIM_P(Cdothi),
                     Real* Patm);

    void FORT_MASSTOMOLE(const int* lo, const int* hi,
			 const Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 Real*       mole, ARLIM_P(molelo), ARLIM_P(molehi));

    void FORT_MOLETOMASS(const int* lo, const int* hi,
			 const Real* mole, ARLIM_P(molelo), ARLIM_P(molehi),
			 Real*       mass, ARLIM_P(masslo), ARLIM_P(masshi));

    void FORT_MASSTP_TO_CONC(const int* lo, const int* hi, const Real* Patm,
			 const Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 const Real* temp, ARLIM_P(templo), ARLIM_P(temphi),
			 Real*       conc, ARLIM_P(conclo), ARLIM_P(conchi));

    void FORT_MASSR_TO_CONC(const int* lo, const int* hi, 
			 const Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 const Real* temp, ARLIM_P(templo), ARLIM_P(temphi),
			 const Real* rho,  ARLIM_P(rholo),  ARLIM_P(rhohi),
			 Real*       conc, ARLIM_P(conclo), ARLIM_P(conchi));

    void FORT_CONC_TO_MOLE(const int* lo, const int* hi, 
			 const Real* conc, ARLIM_P(conclo), ARLIM_P(conchi),
			 Real*       mole, ARLIM_P(molelo), ARLIM_P(molehi));

    void FORT_NORMMASS(const int* lo, const int* hi, const int* id,
                       const Real* Y,  ARLIM_P(Ylo),  ARLIM_P(Yhi),
                       const Real* Yn, ARLIM_P(Ynlo), ARLIM_P(Ynhi));

    int FORT_CONPSOLV(const int* lo, const int* hi,
		       Real* Ynew, ARLIM_P(Ynlo), ARLIM_P(Ynhi),
		       Real* Tnew, ARLIM_P(Tnlo), ARLIM_P(Tnhi),
		       const Real* Yold, ARLIM_P(Yolo), ARLIM_P(Yohi),
		       const Real* Told, ARLIM_P(Tolo), ARLIM_P(Tohi),
		       Real* FuncCount, ARLIM_P(FuncCountlo), ARLIM_P(FuncCounthi),
		       const Real* p, const Real* dt, Real* diag, 
                       const int* do_diag, const int* do_stiff);

#ifdef LMC_SDC
    int FORT_CONPSOLV_SDC(const int* lo, const int* hi,
			   Real* rhoYnew, ARLIM_P(rhoYnlo), ARLIM_P(rhoYnhi),
			   Real* rhoHnew, ARLIM_P(rhoHnlo), ARLIM_P(rhoHnhi),
			   Real* Tnew, ARLIM_P(Tnlo), ARLIM_P(Tnhi),
			   const Real* rhoYold, ARLIM_P(rhoYolo), ARLIM_P(rhoYohi),
			   const Real* rhoHold, ARLIM_P(rhoHolo), ARLIM_P(rhoHohi),
			   const Real* Told, ARLIM_P(Tolo), ARLIM_P(Tohi),
			   const Real* const_src, ARLIM_P(cslo), ARLIM_P(cshi),
			   Real* FuncCount, ARLIM_P(FuncCountlo), ARLIM_P(FuncCounthi),
			   const Real* dt, Real* diag, 
                           const int* do_diag, const int* do_stiff);
#endif

    void FORT_BETA_WBAR(const int* lo, const int* hi,
			const Real* RD, ARLIM_P(RDlo),      ARLIM_P(RDhi),
			Real* RD_Wbar,  ARLIM_P(RDlo_Wbar), ARLIM_P(RDhi_Wbar),
			const Real* Y,  ARLIM_P(Ylo),       ARLIM_P(Yhi));

    void FORT_MIXAVG_RHODIFF_TEMP(const int* lo, const int* hi,
                                  Real*       RD, ARLIM_P(RDlo), ARLIM_P(RDhi),
                                  const Real* T,  ARLIM_P(Tlo),  ARLIM_P(Thi),
                                  const Real* Y,  ARLIM_P(Ylo),  ARLIM_P(Yhi),
                                  const Real* Patm, const int* do_temp,
                                  const int* do_VelVisc);
    
    void FORT_MIX_SHEAR_VISC(const int* lo, const int* hi,
                             Real*       D, ARLIM_P(Dlo), ARLIM_P(Dhi),
                             const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                             const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
    
    void FORT_RHOfromPTY(const int* lo, const int* hi,
			 Real*       R, ARLIM_P(Rlo), ARLIM_P(Rhi),
			 const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
			 const Real* Patm);
      
    void FORT_RHOfromPvTY(const int* lo, const int* hi,
			 Real*       R, ARLIM_P(Rlo), ARLIM_P(Rhi),
			 const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
			 const Real* P, ARLIM_P(Plo), ARLIM_P(Phi));
      
    void FORT_PfromRTY(const int* lo, const int* hi,
		       Real*       P, ARLIM_P(Plo), ARLIM_P(Phi),
		       const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
		       const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
		       const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_TfromPRY(const int* lo, const int* hi,
		       const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
		       const Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
		       const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                       const Real* Patm);
      
    void FORT_CPMIXfromTY(const int* lo, const int* hi,
			  Real*   CPMIX, ARLIM_P(CPlo), ARLIM_P(CPhi),
			  const Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
			  const Real* Y, ARLIM_P(Ylo),  ARLIM_P(Yhi));
    
    void FORT_CVMIXfromTY(const int* lo, const int* hi,
			  Real*   CVMIX, ARLIM_P(CVlo), ARLIM_P(CVhi),
			  const Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
			  const Real* Y, ARLIM_P(Ylo),  ARLIM_P(Yhi));
    
    void FORT_HMIXfromTY(const int* lo, const int* hi,
			 Real*    HMIX, ARLIM_P(Hlo), ARLIM_P(Hhi),
			 const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_MWMIXfromY(const int* lo, const int* hi,
			 Real*    MWMIX,ARLIM_P(MWlo),ARLIM_P(MWhi),
			 const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_CPfromT(const int* lo, const int* hi,
		      Real*      CP, ARLIM_P(CPlo), ARLIM_P(CPhi),
		      const Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi));
    
    void FORT_HfromT(const int* lo, const int* hi,
		     Real*       H, ARLIM_P(Hlo), ARLIM_P(Hhi),
		     const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

    int FORT_TfromHY(const int* lo, const int* hi,
		     Real* T,       ARLIM_P(Tlo), ARLIM_P(Thi),
		     const Real* h, ARLIM_P(hlo), ARLIM_P(hhi),
		     const Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
		     const Real* errMAX, const int* NiterMAX,const Real* tmp);

    void FORT_TfromHYpt(Real* T, const Real* Hin, const Real* Yin,
                        const Real* errMAX, const int* NiterMAX, const Real* res,
                        int* Niter);

    int FORT_GETCKNUMREAC();


    void FORT_MOLPROD(const int* lo, const int* hi, const int* id, 
		      Real*       Q, ARLIM_P(Qlo), ARLIM_P(Qhi),
		      const Real* C, ARLIM_P(Clo), ARLIM_P(Chi),
		      const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

    void FORT_GETELTMOLES(const int* name, const int* nameinlen,
                          const int* lo, const int* hi,
                          Real*       C_elt, ARLIM_P(Celtlo), ARLIM_P(Celthi),
                          const Real* C,     ARLIM_P(Clo),    ARLIM_P(Chi));

    void FORT_OTrad_TDF(const int* lo, const int* hi,
                        Real* Qloss,   ARLIM_P(Qlo), ARLIM_P(Qhi),
                        const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                        const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                        const Real* Patm, const Real* T_bg);
    Real FORT_RUNIV();
    Real FORT_P1ATMMKS();

    int FORT_GET_EG_ID();
    int FORT_GET_TR_ID();
    void FORT_SETTRLIB(const int* lib_idx);

    void FORT_SETNU(int* nu, const int& N);

    void FORT_GET_REACTION_MAP(int* rmap);
}

#endif
#endif

