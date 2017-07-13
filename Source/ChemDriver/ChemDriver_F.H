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
#    define FORT_RRATERHOY       drraterhoy
#    define FORT_RRATEC          drratec
#    define FORT_MOLETOMASS      dmoltoms
#    define FORT_MASSTP_TO_CONC  dmstptocon
#    define FORT_MASSR_TO_CONC   dmsrtocon
#    define FORT_CONC_TO_MOLE    dcontomol
#    define FORT_GET_REACTION_MAP dgrm
#    define FORT_CONPSOLV_SDC    dcpsolvsdc
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
#    define FORT_RRATERHOY       DRRATERHOY
#    define FORT_RRATEC          DRRATEC
#    define FORT_MOLETOMASS      DMOLTOMS
#    define FORT_MASSTP_TO_CONC  DMSTPTOCON
#    define FORT_MASSR_TO_CONC   DMSRTOCON
#    define FORT_CONC_TO_MOLE    DCONTOMOL
#    define FORT_CONPSOLV_SDC    DCPSOLVSDC
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
#    define FORT_RRATERHOY       drraterhoy
#    define FORT_RRATEC          drratec
#    define FORT_MOLETOMASS      dmoltoms
#    define FORT_CONPSOLV_SDC    dcpsolvsdc
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
#    define FORT_RRATERHOY       drraterhoy_
#    define FORT_RRATEC          drratec_
#    define FORT_MOLETOMASS      dmoltoms_
#    define FORT_CONPSOLV_SDC    dcpsolvsdc_
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

#include <AMReX_ArrayLim.H>

extern "C" {

    void FORT_SETTMINTRANS(const amrex::Real* Tmin);
    void FORT_SETVERBOSEVODE();
    void FORT_SETVODETOLS(const amrex::Real* rtol, const amrex::Real* atol, const int* itol);
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
    void FORT_GETCKMWT(amrex::Real* mwt);
    void FORT_GETCKAWT(amrex::Real* awt);
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
                        const amrex::Real* X, ARLIM_P(Xlo),    ARLIM_P(Xhi),
                        const amrex::Real* T, ARLIM_P(Tlo),    ARLIM_P(Thi),
                        amrex::Real*    FwdK, ARLIM_P(FwdKlo), ARLIM_P(FwdKhi),
                        amrex::Real*    RevK, ARLIM_P(RevKlo), ARLIM_P(RevKhi),
                        const amrex::Real* Patm, const int* rxns, const int* Nrxns);

    void FORT_HTRLS(const int* lo, const int* hi,
                    const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                    const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                    amrex::Real*       Q, ARLIM_P(Qlo), ARLIM_P(Qhi), const amrex::Real* Patm);

    void FORT_RRATERHOY(const int* lo, const int* hi,
                        const amrex::Real* RhoY, ARLIM_P(RhoYlo),    ARLIM_P(RhoYhi),
                        const amrex::Real* RhoH, ARLIM_P(RhoHlo),    ARLIM_P(RhoHhi),
                        const amrex::Real* T,    ARLIM_P(Tlo),       ARLIM_P(Thi),
                        amrex::Real*    RhoYdot, ARLIM_P(RhoYdotlo), ARLIM_P(RhoYdothi) );

    void FORT_RRATEC(const int* lo, const int* hi,
                     const amrex::Real* C, ARLIM_P(Clo),    ARLIM_P(Chi),
                     const amrex::Real* T, ARLIM_P(Tlo),    ARLIM_P(Thi),
                     amrex::Real*    Cdot, ARLIM_P(Cdotlo), ARLIM_P(Cdothi),
                     amrex::Real* Patm);

    void FORT_MASSTOMOLE(const int* lo, const int* hi,
			 const amrex::Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 amrex::Real*       mole, ARLIM_P(molelo), ARLIM_P(molehi));

    void FORT_MOLETOMASS(const int* lo, const int* hi,
			 const amrex::Real* mole, ARLIM_P(molelo), ARLIM_P(molehi),
			 amrex::Real*       mass, ARLIM_P(masslo), ARLIM_P(masshi));

    void FORT_MASSTP_TO_CONC(const int* lo, const int* hi, const amrex::Real* Patm,
			 const amrex::Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 const amrex::Real* temp, ARLIM_P(templo), ARLIM_P(temphi),
			 amrex::Real*       conc, ARLIM_P(conclo), ARLIM_P(conchi));

    void FORT_MASSR_TO_CONC(const int* lo, const int* hi, 
			 const amrex::Real* mass, ARLIM_P(masslo), ARLIM_P(masshi),
			 const amrex::Real* temp, ARLIM_P(templo), ARLIM_P(temphi),
			 const amrex::Real* rho,  ARLIM_P(rholo),  ARLIM_P(rhohi),
			 amrex::Real*       conc, ARLIM_P(conclo), ARLIM_P(conchi));

    void FORT_CONC_TO_MOLE(const int* lo, const int* hi, 
			 const amrex::Real* conc, ARLIM_P(conclo), ARLIM_P(conchi),
			 amrex::Real*       mole, ARLIM_P(molelo), ARLIM_P(molehi));

    void FORT_NORMMASS(const int* lo, const int* hi, const int* id,
                       const amrex::Real* Y,  ARLIM_P(Ylo),  ARLIM_P(Yhi),
                       const amrex::Real* Yn, ARLIM_P(Ynlo), ARLIM_P(Ynhi));

    int FORT_CONPSOLV_SDC(const int* lo, const int* hi,
			   amrex::Real* rhoYnew, ARLIM_P(rhoYnlo), ARLIM_P(rhoYnhi),
			   amrex::Real* rhoHnew, ARLIM_P(rhoHnlo), ARLIM_P(rhoHnhi),
			   amrex::Real* Tnew, ARLIM_P(Tnlo), ARLIM_P(Tnhi),
			   const amrex::Real* rhoYold, ARLIM_P(rhoYolo), ARLIM_P(rhoYohi),
			   const amrex::Real* rhoHold, ARLIM_P(rhoHolo), ARLIM_P(rhoHohi),
			   const amrex::Real* Told, ARLIM_P(Tolo), ARLIM_P(Tohi),
			   const amrex::Real* const_src, ARLIM_P(cslo), ARLIM_P(cshi),
			   amrex::Real* FuncCount, ARLIM_P(FuncCountlo), ARLIM_P(FuncCounthi),
			   const amrex::Real* dt, amrex::Real* diag, 
                           const int* do_diag, const int* do_stiff);

    void FORT_BETA_WBAR(const int* lo, const int* hi,
			const amrex::Real* RD, ARLIM_P(RDlo),      ARLIM_P(RDhi),
			amrex::Real* RD_Wbar,  ARLIM_P(RDlo_Wbar), ARLIM_P(RDhi_Wbar),
			const amrex::Real* Y,  ARLIM_P(Ylo),       ARLIM_P(Yhi));

    void FORT_MIXAVG_RHODIFF_TEMP(const int* lo, const int* hi,
                                  amrex::Real*       RD, ARLIM_P(RDlo), ARLIM_P(RDhi),
                                  const amrex::Real* T,  ARLIM_P(Tlo),  ARLIM_P(Thi),
                                  const amrex::Real* Y,  ARLIM_P(Ylo),  ARLIM_P(Yhi),
                                  const amrex::Real* Patm, const int* do_temp,
                                  const int* do_VelVisc);
    
    void FORT_MIX_SHEAR_VISC(const int* lo, const int* hi,
                             amrex::Real*       D, ARLIM_P(Dlo), ARLIM_P(Dhi),
                             const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                             const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
    
    void FORT_RHOfromPTY(const int* lo, const int* hi,
			 amrex::Real*       R, ARLIM_P(Rlo), ARLIM_P(Rhi),
			 const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
			 const amrex::Real* Patm);
      
    void FORT_RHOfromPvTY(const int* lo, const int* hi,
			 amrex::Real*       R, ARLIM_P(Rlo), ARLIM_P(Rhi),
			 const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
			 const amrex::Real* P, ARLIM_P(Plo), ARLIM_P(Phi));
      
    void FORT_PfromRTY(const int* lo, const int* hi,
		       amrex::Real*       P, ARLIM_P(Plo), ARLIM_P(Phi),
		       const amrex::Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
		       const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
		       const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_TfromPRY(const int* lo, const int* hi,
		       const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
		       const amrex::Real* R, ARLIM_P(Rlo), ARLIM_P(Rhi),
		       const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
                       const amrex::Real* Patm);
      
    void FORT_CPMIXfromTY(const int* lo, const int* hi,
			  amrex::Real*   CPMIX, ARLIM_P(CPlo), ARLIM_P(CPhi),
			  const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
			  const amrex::Real* Y, ARLIM_P(Ylo),  ARLIM_P(Yhi));
    
    void FORT_CVMIXfromTY(const int* lo, const int* hi,
			  amrex::Real*   CVMIX, ARLIM_P(CVlo), ARLIM_P(CVhi),
			  const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi),
			  const amrex::Real* Y, ARLIM_P(Ylo),  ARLIM_P(Yhi));
    
    void FORT_HMIXfromTY(const int* lo, const int* hi,
			 amrex::Real*    HMIX, ARLIM_P(Hlo), ARLIM_P(Hhi),
			 const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_MWMIXfromY(const int* lo, const int* hi,
			 amrex::Real*    MWMIX,ARLIM_P(MWlo),ARLIM_P(MWhi),
			 const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi));
      
    void FORT_CPfromT(const int* lo, const int* hi,
		      amrex::Real*      CP, ARLIM_P(CPlo), ARLIM_P(CPhi),
		      const amrex::Real* T, ARLIM_P(Tlo),  ARLIM_P(Thi));
    
    void FORT_HfromT(const int* lo, const int* hi,
		     amrex::Real*       H, ARLIM_P(Hlo), ARLIM_P(Hhi),
		     const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

    int FORT_TfromHY(const int* lo, const int* hi,
		     amrex::Real* T,       ARLIM_P(Tlo), ARLIM_P(Thi),
		     const amrex::Real* h, ARLIM_P(hlo), ARLIM_P(hhi),
		     const amrex::Real* Y, ARLIM_P(Ylo), ARLIM_P(Yhi),
		     const amrex::Real* errMAX, const int* NiterMAX,const amrex::Real* tmp);

    void FORT_TfromHYpt(amrex::Real* T, const amrex::Real* Hin, const amrex::Real* Yin,
                        const amrex::Real* errMAX, const int* NiterMAX, const amrex::Real* res,
                        int* Niter);

    int FORT_GETCKNUMREAC();


    void FORT_MOLPROD(const int* lo, const int* hi, const int* id, 
		      amrex::Real*       Q, ARLIM_P(Qlo), ARLIM_P(Qhi),
		      const amrex::Real* C, ARLIM_P(Clo), ARLIM_P(Chi),
		      const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

    void FORT_GETELTMOLES(const int* name, const int* nameinlen,
                          const int* lo, const int* hi,
                          amrex::Real*       C_elt, ARLIM_P(Celtlo), ARLIM_P(Celthi),
                          const amrex::Real* C,     ARLIM_P(Clo),    ARLIM_P(Chi));

    void FORT_OTrad_TDF(const int* lo, const int* hi,
                        amrex::Real* Qloss,   ARLIM_P(Qlo), ARLIM_P(Qhi),
                        const amrex::Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                        const amrex::Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                        const amrex::Real* Patm, const amrex::Real* T_bg);
    amrex::Real FORT_RUNIV();
    amrex::Real FORT_P1ATMMKS();

    int FORT_GET_EG_ID();
    int FORT_GET_TR_ID();
    void FORT_SETTRLIB(const int* lib_idx);

    void FORT_SETNU(int* nu, const int& N);

    void FORT_GET_REACTION_MAP(int* rmap);
}

#endif
#endif

