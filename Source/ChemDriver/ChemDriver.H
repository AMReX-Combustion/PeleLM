#ifndef _ChemDriver_H_
#define _ChemDriver_H_

#include <string>
#include <map>

#include <Array.H>
#include <FArrayBox.H>
#include <ChemDriver_F.H>

extern "C" {
#include <ReactionData.H>
}

// Selected pointwise conversions w/static linkage callable from Fortran/C, internally verifies
//   that the chemistry has been initialized prior.  BEWARE: argument array lengths cannot be checked!!
extern "C" {
  void CD_XTY_PT(const Real* X, Real* Y); // Mole fraction, X, to mass fraction, Y
  void CD_YTX_PT(const Real* X, Real* Y); // Mass fraction, Y, to mole fraction, X
  void CD_MWT(Real* mwt); // Species molecular weight, mwt
};

class ChemDriver
{
public:

  enum TRANSPORT
  {
    CD_EG, 
    CD_TRANLIB
  };

  class Parameter
  {
  public:
    Parameter(int                reaction_id,
              REACTION_PARAMETER param_id,
              int                species_id = -1);
    std::ostream& operator<<(std::ostream& os) const;
    void operator=(Real new_value);
    const Real& Value() const;
    Real& Value();
    Real DefaultValue() const;
    void ResetToDefault();
    int SpeciesID() const {return species_id;}
    int RxnID() const {return reaction_id;}
    std::string GetParamString() const;
  private:
    int reaction_id;
    REACTION_PARAMETER param_id;
    int species_id;
    double *rp, *rdefp;
  };
  void ResetAllParamsToDefault();


    ChemDriver ();	  

    ~ChemDriver ();

    static void SetTransport(const ChemDriver::TRANSPORT& tran_in);
    static ChemDriver::TRANSPORT Transport ();

    bool solveTransient (FArrayBox&        Ynew,
                         FArrayBox&        Tnew,
                         const FArrayBox&  Yold,
                         const FArrayBox&  Told,
                         FArrayBox&        FuncCount,
                         const Box&        box,
                         int               sCompY,
                         int               sCompT,
                         Real              dt,
                         Real              Patm,
                         FArrayBox*        chemDiag=0,
                         bool              use_stiff_solver = true) const;

#ifdef LMC_SDC
    bool solveTransient_sdc(FArrayBox&        rhoYnew,
                            FArrayBox&        rhoHnew,
                            FArrayBox&        Tnew,
                            const FArrayBox&  rhoYold,
                            const FArrayBox&  rhoHold,
                            const FArrayBox&  Told,
                            const FArrayBox&  const_src,
                            FArrayBox&        FuncCount,
                            const Box&        box,
                            int               sComprhoY,
                            int               sComprhoH,
                            int               sCompT,
                            Real              dt,
                            FArrayBox*        chemDiag,
                            bool              use_stiff_solver = true) const;
#endif

    void set_verbose_vode ();
    void set_max_vode_subcycles (int max_cyc);
    void set_species_Yscales (const std::string& scalesFile);
    //
    // Species info.
    //
    int numSpecies () const;
    int numElements () const;
    int numReactions () const;
    const Array<std::string>& speciesNames () const;
    const Array<std::string>& elementNames () const;
    Array<Real> speciesMolecWt () const;
    Array<Real> elementAtomicWt () const;
    int index (const std::string& specName) const;
    int indexElt (const std::string& eltName) const;
    Array<int> reactionsWithXonL (const std::string& lSpecName) const;
    Array<int> reactionsWithXonR (const std::string& rSpecName) const;
    Array<std::pair<std::string,int> > specCoeffsInReactions (int reacIdx) const;
    std::string reactionString (int reacIdx) const;
    std::string reactionStringBuild (int reacIdx) const;
    void printReactions() const;
    const Array<int>& reactionMap() const;
    const Array<int>& reactionReverseMap() const;

    int numberOfElementXinSpeciesY (const std::string& eltX,
                                    const std::string& spcY) const;
    //
    // Thermo info.
    //
    void getRhoGivenPTY (FArrayBox&       Rho,
                         Real             Patm,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompR) const;

    void getRhoGivenPvTY (FArrayBox&       Rho,
                          const FArrayBox& P,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompP,
                          int              sCompT,
                          int              sCompY,
                          int              sCompR) const;

    void getPGivenRTY (FArrayBox&       p,
                       const FArrayBox& Rho,
                       const FArrayBox& T,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompT,
                       int              sCompY,
                       int              sCompP) const;

    void getTGivenPRY (FArrayBox&       T,
                       Real             p,
                       const FArrayBox& Rho,
                       const FArrayBox& Y,
                       const Box&       box,
                       int              sCompR,
                       int              sCompY,
                       int              sCompT) const;

    void getCpmixGivenTY (FArrayBox&       cpmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCp) const;

    void getCvmixGivenTY (FArrayBox&       cvmix,
                          const FArrayBox& T,
                          const FArrayBox& Y,
                          const Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCv) const;

    void getHmixGivenTY (FArrayBox&       hmix,
                         const FArrayBox& T,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

    void getMwmixGivenY (FArrayBox&       mwmix,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompY,
                         int              sCompMw) const;

    void getCpGivenT (FArrayBox&       cp,
                      const FArrayBox& T,
                      const Box&       box,
                      int              sCompT,
                      int              sCompCp) const;

    void getHGivenT (FArrayBox&       h,
                     const FArrayBox& T,
                     const Box&       box,
                     int              sCompT,
                     int              sCompH) const;

    Real getRuniversal () const;
    Real getP1atm_MKS () const;
    //
    // Compute T that satisfies hmix=sum(h(T)), returns max Newton iterations
    // on any point over grid, return -1 if T jumped out of bounds during solve,
    // and -2 if solve failed anywhere.  Note that a temporary is not used, and
    // the solver kicks out when/if it fails, so it may exit after converting
    // only part of the temperature array.  Save a copy of T and check the return
    // code if you want to be extra careful...
    //
    int getTGivenHY (FArrayBox&       T,
                     const FArrayBox& H,
                     const FArrayBox& Y,
                     const Box&       box,
                     int              sCompH,
                     int              sCompY,
                     int              sCompT,
                     const Real&      errMAX = -1) const;

    Array<Real> massFracToMoleFrac (const Array<Real>& Y) const;
    Array<Real> moleFracToMassFrac (const Array<Real>& X) const;

    void molarProduction (FArrayBox&       Q,
                          const std::string&   specName,
                          const FArrayBox& C,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompC,
                          int              sCompT,
                          int              sCompQ) const;
    //
    // Compute heat release (J/(m^3.s)) based on the temp, press and mass fractions
    //
    void heatRelease (FArrayBox&       Q,
                      const FArrayBox& Y,
                      const FArrayBox& T,
                      Real             Patm,
                      const Box&       box,
                      int              sCompY,
                      int              sCompT,
                      int              sCompQ) const;
    //
    // Compute dY/dt based on the input temp, press and mass fractions.
    //
    void reactionRateY (FArrayBox&       Ydot,
                        const FArrayBox& Y,
                        const FArrayBox& T,
                        Real             Patm,
                        const Box&       box,
                        int              sCompY,
                        int              sCompT,
                        int              sCompYdot) const;

#ifdef LMC_SDC
    //
    // Compute dRhoY/dt based on the input temp, press and mass densities.
    //
    void reactionRateRhoY(FArrayBox&       RhoYdot,
                          const FArrayBox& RhoY,
                          const FArrayBox& RhoH,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompRhoY,
                          int              sCompRhoH,
                          int              sCompT,
                          int              sCompRhoYdot) const;
#endif

    void fwdRevReacRatesGivenXTP (FArrayBox&        FwdK,
                                  FArrayBox&        RevK,
                                  const Array<int>& rxnIDs,
                                  const FArrayBox&  X,
                                  const FArrayBox&  T,
                                  Real              Patm,
                                  const Box&        box,
                                  int               sCompX,
                                  int               sCompT,
                                  int               sCompFwdK,
                                  int               sCompRevK) const;
    
    void massFracToMoleFrac (FArrayBox&       X,
                             const FArrayBox& Y,
                             const Box&       box,
                             int              sCompY,
                             int              sCompX) const;

    void moleFracToMassFrac (FArrayBox&       Y,
                             const FArrayBox& X,
                             const Box&       box,
                             int              sCompX,
                             int              sCompY) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              Real             Patm,
                              const Box&       box,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void massFracToMolarConc (FArrayBox&       C,
                              const FArrayBox& Y,
                              const FArrayBox& T,
                              const FArrayBox& Rho,
                              const Box&       box,
                              int              sCompR,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void molarConcToMoleFrac (FArrayBox&       X,
                              const FArrayBox& C,
                              const Box&       box,
                              int              sCompC,
                              int              sCompX) const;
    //
    // Normalize mass fractions to prevent negatives.
    //
    void normalizeMassFrac (FArrayBox&       Ynorm,
                            const FArrayBox& Y,
                            const std::string&   excessSpecies,
                            const Box&       box,
                            int              sCompY,
                            int              sCompYnorm) const;
    //
    // Chemical Diffusivities.
    //
    void getMixAveragedRhoDiff (FArrayBox&       rhoD,
                                const FArrayBox& Y,
                                const FArrayBox& T,
                                Real             Patm,
                                const Box&       box,
                                int              sCompY,
                                int              sCompT,
                                int              sCompRD) const;
    //
    // Viscosity.
    //
    void getMixShearVisc (FArrayBox&       eta,
                          const FArrayBox& Y,
                          const FArrayBox& T,
                          const Box&       box,
                          int              sCompY,
                          int              sCompT,
                          int              sCompE) const;
    
    void getElementMoles (FArrayBox&       C_elt,
                          const std::string&   name,
                          const FArrayBox& C,
                          const Box&       box,
                          int              sCompC,
                          int              sCompC_elt) const;
    //
    // Optically thin radiation model.
    //
    void getOTradLoss_TDF (FArrayBox&       Qloss,
                           const FArrayBox& T,
                           const FArrayBox& X,
                           const Real       Patm,
                           const Real       T_bg,
                           const Box&       box,
                           int              sCompX,
                           int              sCompT,
                           int              sCompQ) const;
    //
    // H - to - T solve parameter access.
    //
    Real getHtoTerrMAX () const;
    void setHtoTerrMAX (Real err);
    int getHtoTiterMAX () const;
    void setHtoTiterMAX (int err);    
    //
    // Handy functions.
    //
    static Array<int> encodeStringForFortran(const std::string& astr);
    static std::string decodeStringFromFortran(const int* coded, int length);

    //
    // Data structure for generating chemical path diagrams,
    //   an Edge represents the transfer of an atom from one chemical
    //   species, sp1, to another, sp2, which occurs through a list
    //   of reactions.  For each reaction, n of these atoms are
    //   transferred (RWL is a list of pairs of rxn ID to atom transfer cnt)
    //
    class Edge
    {
    public:
        friend std::ostream& operator<< (std::ostream& os, const Edge& e);

        Edge () {}

        Edge (const std::string& n1,
	      const std::string& n2,
	      const Array<std::pair<int,Real> > rwl,
	      const ChemDriver* CD);
        
        Edge (const std::string& n1,
	      const std::string& n2,
	      int reac,
	      Real weight,
	      const ChemDriver* CD);
        
        int equivSign (const Edge& rhs) const;

        void combine (const Edge& rhs, int sgn);

        bool touchesSp(const std::string& rhs) const;
        
        void reverse();

        const Array<std::pair<int,Real> >& rwl () const;

        const std::string left() const;

        const std::string right() const;

        bool operator< (const Edge& rhs) const;

        const Array<std::pair<int,Real> > RateWeightList() const {return RWL;}

        const ChemDriver* CD() const {return cd;}

    private:
        std::string sp1, sp2;
        Array<std::pair<int,Real> > RWL; // RateWeightList, each pair is (reac,coeff)
        const ChemDriver* cd;
    };

    // 
    // Helper class for building chem edges.  A group is a list of constituent 
    // atoms, and this class allows a coupleof useful math operations on these
    // groups.  
    //
    class Group
    {
    public:
        friend std::ostream& operator<< (std::ostream& os, const Group& g);
        
        Group () {}

        Group (const std::map<std::string,int>& eltCnts);

        Group (const Group& rhs);

        Group operator- (const Group& rhs) const;

        Group operator* (int rhs) const;

        int operator[](const std::string& id) const;
            
        bool sameSign() const;

        bool contains(const std::string& id) const;
    
        Real awt(); // non-const because embedded lazy eval

        int size() const;
    
    private:
        void FillAtomicWeights ();
        std::map<std::string,int> mEltCnts;
        static std::map<std::string,Real> AtomicWeight;
    };

    // 
    // Compute edges from chem mechanism
    //
    std::list<Edge> getEdges (const std::string& trElt,
                              int PrintVerbose=0,
                              int HackSplitting=1) const;
    
protected:

    void getSpeciesNames ();
    void getElementNames ();
    void getStoichCoeffs ();

private:

    void initOnce ();

    Array<std::string> mSpeciesNames;
    Array<std::string> mElementNames;
    Real               mHtoTerrMAX;
    int                mHtoTiterMAX;
    Array<Real>        mTmpData;
    int mMaxreac, mMaxspec, mMaxelts, mMaxord, mMaxthrdb, mMaxtp, mMaxsp, mMaxspnml;
    Array<int>         mNu;
    Array<int> reaction_map, reaction_rev_map;
};

std::ostream& operator<< (std::ostream&  os, const ChemDriver::Parameter& param);

ChemDriver::Group operator* (int n, const ChemDriver::Group& g);

inline
int
ChemDriver::numSpecies () const
{
    return mSpeciesNames.size();
}

inline
int
ChemDriver::numElements () const
{
    return mElementNames.size();
}

inline
const Array<std::string>&
ChemDriver::speciesNames () const
{
    return mSpeciesNames;
}

inline
const Array<std::string>&
ChemDriver::elementNames () const
{
    return mElementNames;
}

inline
int
ChemDriver::index (const std::string& specName) const
{
    for (int i=0; i<mSpeciesNames.size(); i++)
    {
        if (specName == mSpeciesNames[i])
            return i;
    }
    return -1;
}

inline
int
ChemDriver::indexElt (const std::string& eltName) const
{
    for (int i=0; i<mElementNames.size(); i++)
    {
        if (eltName == mElementNames[i])
            return i;
    }
    return -1;
}

inline
Real
ChemDriver::getHtoTerrMAX () const
{
    return mHtoTerrMAX;
}

inline
void
ChemDriver::setHtoTerrMAX (Real err)
{
    mHtoTerrMAX = err;
}

inline
int
ChemDriver::getHtoTiterMAX () const
{
    return mHtoTiterMAX;
}

inline
void
ChemDriver::setHtoTiterMAX (int iter)
{
    if (iter != mHtoTiterMAX)
    {
        mHtoTiterMAX = iter;
        mTmpData.resize(mHtoTiterMAX);
    }
}

inline
const Array<int>&
ChemDriver::reactionMap() const
{
  return reaction_map;
}

inline
const Array<int>&
ChemDriver::reactionReverseMap() const
{
  return reaction_rev_map;
}

#endif /*_ChemDriver_H_*/
