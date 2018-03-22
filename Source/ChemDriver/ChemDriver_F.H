#ifndef _ChemDriver_F_H_
#define _ChemDriver_F_H_

#ifdef BL_LANG_FORT

#else
#  if defined(BL_FORT_USE_UPPERCASE)

#  elif defined(BL_FORT_USE_LOWERCASE)

#  elif defined(BL_FORT_USE_UNDERSCORE)

#  endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif
#endif

