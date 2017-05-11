#ifndef _ChemDriver_H_
#define _ChemDriver_H_

#include <string>
#include <map>

#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <ChemDriver_F.H>

extern "C" {
#include <ReactionData.H>
}

// Selected pointwise conversions w/static linkage callable from Fortran/C, internally verifies
//   that the chemistry has been initialized prior.  BEWARE: argument array lengths cannot be checked!!
extern "C" {
  void CD_XTY_PT(const amrex::Real* X, amrex::Real* Y); // Mole fraction, X, to mass fraction, Y
  void CD_YTX_PT(const amrex::Real* X, amrex::Real* Y); // Mass fraction, Y, to mole fraction, X
  void CD_MWT(amrex::Real* mwt); // Species molecular weight, mwt
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
    void operator=(amrex::Real new_value);
    const amrex::Real& Value() const;
    amrex::Real& Value();
    amrex::Real DefaultValue() const;
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

    bool solveTransient_sdc(amrex::FArrayBox&        rhoYnew,
                            amrex::FArrayBox&        rhoHnew,
                            amrex::FArrayBox&        Tnew,
                            const amrex::FArrayBox&  rhoYold,
                            const amrex::FArrayBox&  rhoHold,
                            const amrex::FArrayBox&  Told,
                            const amrex::FArrayBox&  const_src,
                            amrex::FArrayBox&        FuncCount,
                            const amrex::Box&        box,
                            int               sComprhoY,
                            int               sComprhoH,
                            int               sCompT,
                            amrex::Real              dt,
                            amrex::FArrayBox*        chemDiag,
                            bool              use_stiff_solver = true) const;

    void set_verbose_vode ();
    void set_max_vode_subcycles (int max_cyc);
    void set_species_Yscales (const std::string& scalesFile);
    //
    // Species info.
    //
    int numSpecies () const;
    int numElements () const;
    int numReactions () const;
    const amrex::Array<std::string>& speciesNames () const;
    const amrex::Array<std::string>& elementNames () const;
    amrex::Array<amrex::Real> speciesMolecWt () const;
    amrex::Array<amrex::Real> elementAtomicWt () const;
    int index (const std::string& specName) const;
    int indexElt (const std::string& eltName) const;
    amrex::Array<int> reactionsWithXonL (const std::string& lSpecName) const;
    amrex::Array<int> reactionsWithXonR (const std::string& rSpecName) const;
    amrex::Array<std::pair<std::string,int> > specCoeffsInReactions (int reacIdx) const;
    std::string reactionString (int reacIdx) const;
    std::string reactionStringBuild (int reacIdx) const;
    void printReactions() const;
    const amrex::Array<int>& reactionMap() const;
    const amrex::Array<int>& reactionReverseMap() const;

    int numberOfElementXinSpeciesY (const std::string& eltX,
                                    const std::string& spcY) const;
    //
    // Thermo info.
    //
    void getRhoGivenPTY (amrex::FArrayBox&       Rho,
                         amrex::Real             Patm,
                         const amrex::FArrayBox& T,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompR) const;

    void getRhoGivenPvTY (amrex::FArrayBox&       Rho,
                          const amrex::FArrayBox& P,
                          const amrex::FArrayBox& T,
                          const amrex::FArrayBox& Y,
                          const amrex::Box&       box,
                          int              sCompP,
                          int              sCompT,
                          int              sCompY,
                          int              sCompR) const;

    void getPGivenRTY (amrex::FArrayBox&       p,
                       const amrex::FArrayBox& Rho,
                       const amrex::FArrayBox& T,
                       const amrex::FArrayBox& Y,
                       const amrex::Box&       box,
                       int              sCompR,
                       int              sCompT,
                       int              sCompY,
                       int              sCompP) const;

    void getTGivenPRY (amrex::FArrayBox&       T,
                       amrex::Real             p,
                       const amrex::FArrayBox& Rho,
                       const amrex::FArrayBox& Y,
                       const amrex::Box&       box,
                       int              sCompR,
                       int              sCompY,
                       int              sCompT) const;

    void getCpmixGivenTY (amrex::FArrayBox&       cpmix,
                          const amrex::FArrayBox& T,
                          const amrex::FArrayBox& Y,
                          const amrex::Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCp) const;

    void getCvmixGivenTY (amrex::FArrayBox&       cvmix,
                          const amrex::FArrayBox& T,
                          const amrex::FArrayBox& Y,
                          const amrex::Box&       box,
                          int              sCompT,
                          int              sCompY,
                          int              sCompCv) const;

    void getHmixGivenTY (amrex::FArrayBox&       hmix,
                         const amrex::FArrayBox& T,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompT,
                         int              sCompY,
                         int              sCompH) const;

    void getMwmixGivenY (amrex::FArrayBox&       mwmix,
                         const amrex::FArrayBox& Y,
                         const amrex::Box&       box,
                         int              sCompY,
                         int              sCompMw) const;

    void getCpGivenT (amrex::FArrayBox&       cp,
                      const amrex::FArrayBox& T,
                      const amrex::Box&       box,
                      int              sCompT,
                      int              sCompCp) const;

    void getHGivenT (amrex::FArrayBox&       h,
                     const amrex::FArrayBox& T,
                     const amrex::Box&       box,
                     int              sCompT,
                     int              sCompH) const;

    amrex::Real getRuniversal () const;
    amrex::Real getP1atm_MKS () const;
    //
    // Compute T that satisfies hmix=sum(h(T)), returns max Newton iterations
    // on any point over grid, return -1 if T jumped out of bounds during solve,
    // and -2 if solve failed anywhere.  Note that a temporary is not used, and
    // the solver kicks out when/if it fails, so it may exit after converting
    // only part of the temperature array.  Save a copy of T and check the return
    // code if you want to be extra careful...
    //
    int getTGivenHY (amrex::FArrayBox&       T,
                     const amrex::FArrayBox& H,
                     const amrex::FArrayBox& Y,
                     const amrex::Box&       box,
                     int              sCompH,
                     int              sCompY,
                     int              sCompT,
                     const amrex::Real&      errMAX = -1) const;

    amrex::Array<amrex::Real> massFracToMoleFrac (const amrex::Array<amrex::Real>& Y) const;
    amrex::Array<amrex::Real> moleFracToMassFrac (const amrex::Array<amrex::Real>& X) const;

    void molarProduction (amrex::FArrayBox&       Q,
                          const std::string&   specName,
                          const amrex::FArrayBox& C,
                          const amrex::FArrayBox& T,
                          const amrex::Box&       box,
                          int              sCompC,
                          int              sCompT,
                          int              sCompQ) const;
    //
    // Compute heat release (J/(m^3.s)) based on the temp, press and mass fractions
    //
    void heatRelease (amrex::FArrayBox&       Q,
                      const amrex::FArrayBox& Y,
                      const amrex::FArrayBox& T,
                      amrex::Real             Patm,
                      const amrex::Box&       box,
                      int              sCompY,
                      int              sCompT,
                      int              sCompQ) const;

    //
    // Compute dRhoY/dt based on the input temp, press and mass densities.
    //
    void reactionRateRhoY(amrex::FArrayBox&       RhoYdot,
                          const amrex::FArrayBox& RhoY,
                          const amrex::FArrayBox& RhoH,
                          const amrex::FArrayBox& T,
                          const amrex::Box&       box,
                          int              sCompRhoY,
                          int              sCompRhoH,
                          int              sCompT,
                          int              sCompRhoYdot) const;

    void fwdRevReacRatesGivenXTP (amrex::FArrayBox&        FwdK,
                                  amrex::FArrayBox&        RevK,
                                  const amrex::Array<int>& rxnIDs,
                                  const amrex::FArrayBox&  X,
                                  const amrex::FArrayBox&  T,
                                  amrex::Real              Patm,
                                  const amrex::Box&        box,
                                  int               sCompX,
                                  int               sCompT,
                                  int               sCompFwdK,
                                  int               sCompRevK) const;
    
    void massFracToMoleFrac (amrex::FArrayBox&       X,
                             const amrex::FArrayBox& Y,
                             const amrex::Box&       box,
                             int              sCompY,
                             int              sCompX) const;

    void moleFracToMassFrac (amrex::FArrayBox&       Y,
                             const amrex::FArrayBox& X,
                             const amrex::Box&       box,
                             int              sCompX,
                             int              sCompY) const;

    void massFracToMolarConc (amrex::FArrayBox&       C,
                              const amrex::FArrayBox& Y,
                              const amrex::FArrayBox& T,
                              amrex::Real             Patm,
                              const amrex::Box&       box,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void massFracToMolarConc (amrex::FArrayBox&       C,
                              const amrex::FArrayBox& Y,
                              const amrex::FArrayBox& T,
                              const amrex::FArrayBox& Rho,
                              const amrex::Box&       box,
                              int              sCompR,
                              int              sCompY,
                              int              sCompT,
                              int              sCompC) const;

    void molarConcToMoleFrac (amrex::FArrayBox&       X,
                              const amrex::FArrayBox& C,
                              const amrex::Box&       box,
                              int              sCompC,
                              int              sCompX) const;
    //
    // Normalize mass fractions to prevent negatives.
    //
    void normalizeMassFrac (amrex::FArrayBox&       Ynorm,
                            const amrex::FArrayBox& Y,
                            const std::string&   excessSpecies,
                            const amrex::Box&       box,
                            int              sCompY,
                            int              sCompYnorm) const;
    //
    // Chemical Diffusivities.
    //
    void getMixAveragedRhoDiff (amrex::FArrayBox&       rhoD,
                                const amrex::FArrayBox& Y,
                                const amrex::FArrayBox& T,
                                amrex::Real             Patm,
                                const amrex::Box&       box,
                                int              sCompY,
                                int              sCompT,
                                int              sCompRD) const;
    //
    // Viscosity.
    //
    void getMixShearVisc (amrex::FArrayBox&       eta,
                          const amrex::FArrayBox& Y,
                          const amrex::FArrayBox& T,
                          const amrex::Box&       box,
                          int              sCompY,
                          int              sCompT,
                          int              sCompE) const;
    
    void getElementMoles (amrex::FArrayBox&       C_elt,
                          const std::string&   name,
                          const amrex::FArrayBox& C,
                          const amrex::Box&       box,
                          int              sCompC,
                          int              sCompC_elt) const;
    //
    // Optically thin radiation model.
    //
    void getOTradLoss_TDF (amrex::FArrayBox&       Qloss,
                           const amrex::FArrayBox& T,
                           const amrex::FArrayBox& X,
                           const amrex::Real       Patm,
                           const amrex::Real       T_bg,
                           const amrex::Box&       box,
                           int              sCompX,
                           int              sCompT,
                           int              sCompQ) const;
    //
    // H - to - T solve parameter access.
    //
    amrex::Real getHtoTerrMAX () const;
    void setHtoTerrMAX (amrex::Real err);
    int getHtoTiterMAX () const;
    void setHtoTiterMAX (int err);    
    //
    // Handy functions.
    //
    static amrex::Array<int> encodeStringForFortran(const std::string& astr);
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
	      const amrex::Array<std::pair<int,amrex::Real> > rwl,
	      const ChemDriver* CD);
        
        Edge (const std::string& n1,
	      const std::string& n2,
	      int reac,
	      amrex::Real weight,
	      const ChemDriver* CD);
        
        int equivSign (const Edge& rhs) const;

        void combine (const Edge& rhs, int sgn);

        bool touchesSp(const std::string& rhs) const;
        
        void reverse();

        const amrex::Array<std::pair<int,amrex::Real> >& rwl () const;

        const std::string left() const;

        const std::string right() const;

        bool operator< (const Edge& rhs) const;

        const amrex::Array<std::pair<int,amrex::Real> > RateWeightList() const {return RWL;}

        const ChemDriver* CD() const {return cd;}

    private:
        std::string sp1, sp2;
        amrex::Array<std::pair<int,amrex::Real> > RWL; // RateWeightList, each pair is (reac,coeff)
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
    
        amrex::Real awt(); // non-const because embedded lazy eval

        int size() const;
    
    private:
        void FillAtomicWeights ();
        std::map<std::string,int> mEltCnts;
        static std::map<std::string,amrex::Real> AtomicWeight;
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

    amrex::Array<std::string> mSpeciesNames;
    amrex::Array<std::string> mElementNames;
    amrex::Real               mHtoTerrMAX;
    int                mHtoTiterMAX;
    amrex::Array<amrex::Real>        mTmpData;
    int mMaxreac, mMaxspec, mMaxelts, mMaxord, mMaxthrdb, mMaxtp, mMaxsp, mMaxspnml;
    amrex::Array<int>         mNu;
    amrex::Array<int> reaction_map, reaction_rev_map;
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
const amrex::Array<std::string>&
ChemDriver::speciesNames () const
{
    return mSpeciesNames;
}

inline
const amrex::Array<std::string>&
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
amrex::Real
ChemDriver::getHtoTerrMAX () const
{
    return mHtoTerrMAX;
}

inline
void
ChemDriver::setHtoTerrMAX (amrex::Real err)
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
const amrex::Array<int>&
ChemDriver::reactionMap() const
{
  return reaction_map;
}

inline
const amrex::Array<int>&
ChemDriver::reactionReverseMap() const
{
  return reaction_rev_map;
}

#endif /*_ChemDriver_H_*/
