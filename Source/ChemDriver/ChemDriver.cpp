#include <winstd.H>
#include <iostream>

#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <sstream>

const Real HtoTerrMAX_DEF  = 1.e-8;
const int  HtoTiterMAX_DEF = 20;
const Real Tmin_trans_DEF  = 0.;
const ChemDriver::TRANSPORT transport_DEF = ChemDriver::CD_EG;

namespace
{
    bool initialized = false;
    ChemDriver::TRANSPORT transport = transport_DEF;
    static int number_of_species = -1;

    void ChemDriver_Finalize() {
      initialized = false;
      transport = transport_DEF;
      number_of_species = -1;
    }
}

extern "C" {
  int FORT_USINGEG() {
    return (int)(transport == ChemDriver::CD_EG);
  }

  int FORT_USINGMC() {
    return (int)(transport == ChemDriver::CD_TRANLIB);
  }
}

ChemDriver::ChemDriver ()
    :
    mHtoTerrMAX(HtoTerrMAX_DEF),
    mHtoTiterMAX(HtoTiterMAX_DEF)
{
    FORT_SETTMINTRANS(&Tmin_trans_DEF);

    //FORT_SETVERBOSEVODE();

    if (!initialized)
    {
        initOnce();
        BoxLib::ExecOnFinalize(ChemDriver_Finalize);
        initialized = true;
    }
    mTmpData.resize(mHtoTiterMAX);
}

ChemDriver::~ChemDriver ()
{
    FORT_FINALIZECHEM();
}

extern "C" {
  double* GetParamPtr(int                reaction_id,
                      REACTION_PARAMETER param_id,
                      int                species_id,
                      int                get_default);
  void ResetAllParametersToDefault();
};

ChemDriver::Parameter::Parameter(int                _reaction_id,
                                 REACTION_PARAMETER _param_id,
                                 int                _species_id)
  :reaction_id(_reaction_id),
   param_id(_param_id),
   species_id(_species_id)
{
  rp = GetParamPtr(reaction_id,param_id,species_id,0);
  rdefp = GetParamPtr(reaction_id,param_id,species_id,1);
}

const Real&
ChemDriver::Parameter::Parameter::Value() const
{
  return *rp;
}

Real&
ChemDriver::Parameter::Parameter::Value()
{
  return *rp;
}

Real
ChemDriver::Parameter::Parameter::DefaultValue() const
{
  return *rdefp;
}

void
ChemDriver::Parameter::Parameter::ResetToDefault()
{
  *rp = *rdefp;
}

std::string
ChemDriver::Parameter::GetParamString() const
{
  if (param_id == FWD_A)          {return "FWD_A";}
  else if (param_id == FWD_BETA)  {return "FWD_BETA";}
  else if (param_id == FWD_EA)    {return "FWD_EA";}
  else if (param_id == LOW_A)     {return "LOW_A";}
  else if (param_id == LOW_BETA)  {return "LOW_BETA";}
  else if (param_id == LOW_EA)    {return "LOW_EA";}
  else if (param_id == REV_A)     {return "REV_A";}
  else if (param_id == REV_BETA)  {return "REV_BETA";}
  else if (param_id == REV_EA)    {return "REV_EA";}
  else if (param_id == TROE_A)    {return "TROE_A";}
  else if (param_id == TROE_TS)   {return "TROE_TS";}
  else if (param_id == TROE_TSS)  {return "TROE_TSS";}
  else if (param_id == TROE_TSSS) {return "TROE_TSSS";}
  else if (param_id == SRI_A)     {return "SRI_A";}
  else if (param_id == SRI_B)     {return "SRI_B";}
  else if (param_id == SRI_C)     {return "SRI_C";}
  else if (param_id == SRI_D)     {return "SRI_D";}
  else if (param_id == SRI_E)     {return "SRI_E";}
  else if (param_id == THIRD_BODY) {return "THIRD_BODY";}
  else{
      BoxLib::Abort("Unknown reaction parameter");
  }

}
std::ostream&
ChemDriver::Parameter::operator<<(std::ostream& os) const
{
  // FIXME: Hard-coded list of parameters invites errors later...
  // RG: Moved it to its own routine so I could use it elsewhere,
  // inviting the same error
  std::string param_str;
  os << "(r: " << reaction_id << " p: ";
  param_str = GetParamString();
  os << param_str;
  if( param_id == THIRD_BODY ){
    os << " s: " << species_id;
  }
  os << " v: " << Value() << " DEF: " << DefaultValue() << ")";
  return os;
}

std::ostream& operator<< (std::ostream&  os, const ChemDriver::Parameter& param)
{
  return param.operator<<(os);
}

void
ChemDriver::Parameter::operator=(Real new_value)
{
  Real& val = Value();
  val = new_value;
}

void
ChemDriver::ResetAllParamsToDefault()
{
  ResetAllParametersToDefault();
}


void
ChemDriver::SetTransport (const ChemDriver::TRANSPORT& tran_in)
{
  if (initialized) {
    BoxLib::Abort("Must set transport model before constructing the ChemDriver");
  }
  transport = tran_in;
}

ChemDriver::TRANSPORT
ChemDriver::Transport ()
{
  return transport;
}

static void modify_parameters(ChemDriver& cd)
{
  /*
    To modify a chem parameter, p1 ... pn, add the following to the ParmParsed input, e.g.
    chem.parameters = p1 p2 ...
    chem.p1.type = FWD_A
    chem.p1.reaction_id = 0
    chem.p2.type = FWD_EA
    chem.p2.reaction_id = 0
    chem.values = 6.e17 6000
    ...
   */

  // FIXME: Redo this to avoid explicitly defining this map twice...
  std::map<std::string,REACTION_PARAMETER> PTypeMap;
  PTypeMap["FWD_A"]      = FWD_A;
  PTypeMap["FWD_BETAA"]  = FWD_BETA;
  PTypeMap["FWD_EA"]     = FWD_EA;
  PTypeMap["LOW_A"]      = LOW_A;
  PTypeMap["LOW_BETAA"]  = LOW_BETA;
  PTypeMap["LOW_EA"]     = LOW_EA;
  PTypeMap["REV_A"]      = REV_A;
  PTypeMap["REV_BETAA"]  = REV_BETA;
  PTypeMap["REV_EA"]     = REV_EA;
  PTypeMap["TROE_A"]     = TROE_A;
  PTypeMap["TROE_A"]     = TROE_A;
  PTypeMap["TROE_TS"]    = TROE_TS;
  PTypeMap["TROE_TSS"]   = TROE_TSS;
  PTypeMap["TROE_TSSS"]  = TROE_TSSS;
  PTypeMap["SRI_A"]      = SRI_A;
  PTypeMap["SRI_B"]      = SRI_B;
  PTypeMap["SRI_C"]      = SRI_C;
  PTypeMap["SRI_D"]      = SRI_D;
  PTypeMap["SRI_E"]      = SRI_E;
  PTypeMap["THIRD_BODY"] = THIRD_BODY;

  ParmParse pp("chem");
  Array<std::string> parameters;
  int np = pp.countval("parameters");
  if (np>0) {
    pp.getarr("parameters",parameters,0,np);
  }

  PArray<ChemDriver::Parameter> p(np,PArrayManage);
  std::map<int,PArray<ChemDriver::Parameter> > dependent_parameters;

  Array<Real> values(np);
  for (int i=0; i<np; ++i) {
    std::string prefix = parameters[i];
    ParmParse ppp(std::string("chem."+prefix).c_str());
    int reaction_id; ppp.get("reaction_id",reaction_id);
    if (reaction_id<0 || reaction_id > cd.numReactions()) {
      BoxLib::Abort("Reaction ID invalid");
    }

    std::string type; ppp.get("type",type);
    std::map<std::string,REACTION_PARAMETER>::const_iterator it = PTypeMap.find(type);
    if (it == PTypeMap.end()) {
      BoxLib::Abort("Unrecognized reaction parameter");
    }

    int id = -1;
    if (type == "THIRD_BODY") {
      std::string tb_name; ppp.get("tb_name",tb_name);
      id = cd.index(tb_name);
      BL_ASSERT(id >= 0);
    }

    p.set(i,new ChemDriver::Parameter(reaction_id,it->second,id));
    values[i] = p[i].DefaultValue();

    int ndp = ppp.countval("dependent_parameters");
    if (ndp > 0) {
      dependent_parameters[i].resize(ndp);
      Array<std::string> dpnames; ppp.getarr("dependent_parameters",dpnames,0,ndp);
      for (int j=0; j<ndp; ++j) {

	std::string dplist = "chem." + prefix + ".dependent_parameters." + dpnames[j];
	ParmParse pppd(dplist.c_str());

	int dpreaction_id; pppd.get("reaction_id",dpreaction_id);
	if (dpreaction_id<0 || dpreaction_id > cd.numReactions()) {
	  BoxLib::Abort("Dependent reaction ID invalid");
	}

	std::string dptype; pppd.get("type",dptype);
	std::map<std::string,REACTION_PARAMETER>::const_iterator it = PTypeMap.find(dptype);
	if (it == PTypeMap.end()) {
	  BoxLib::Abort("Unrecognized dependent reaction parameter");
	}

	int did = -1;
	if (dptype == "THIRD_BODY") {
	  std::string dtb_name; pppd.get("tb_name",dtb_name);
	  did = cd.index(dtb_name);
	  BL_ASSERT(id >= 0);
	}

	dependent_parameters[i].set(j, new ChemDriver::Parameter(dpreaction_id,it->second,did));
      }
    }
  }

  int nv = pp.countval("parameter_values");
  if (nv!=0) {
    BL_ASSERT(nv == np);
    pp.getarr("parameter_values",values,0,np);
    for (int i=0; i<np; ++i) {
      p[i].Value() = values[i];
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "************** Modified chem parameter \"" << parameters[i] << "\": " << p[i] << std::endl;
      }
      for (int j=0; j<dependent_parameters[i].size(); ++j) {
	dependent_parameters[i][j].Value() = values[i];
	if (ParallelDescriptor::IOProcessor()) {
	  std::cout << "                     dependent parameter: " << dependent_parameters[i][j] << std::endl;
	}
      }
    }
  }
}

void
ChemDriver::initOnce ()
{
    FORT_INITCHEM();
    getSpeciesNames();
    getElementNames();
    modify_parameters(*this);
    FORT_GETCKDIMPARAMS(&mMaxreac, &mMaxspec, &mMaxelts,  &mMaxord,
                        &mMaxthrdb, &mMaxtp,  &mMaxsp,    &mMaxspnml);
    getStoichCoeffs();

    ParmParse pp("ht");

    int  v_itol = 1;
    Real v_rtol = 1.e-10;
    Real v_atol = 1.e-10;

    pp.query("vode_itol",v_itol);
    pp.query("vode_rtol",v_rtol);
    pp.query("vode_atol",v_atol);

    BL_ASSERT(v_rtol > 0);
    BL_ASSERT(v_atol > 0);
    BL_ASSERT(v_itol == 1 || v_itol == 2);

    FORT_SETVODETOLS(&v_rtol,&v_atol,&v_itol);

    int  v_maxcyc = -1;

    pp.query("vode_max_subcycles",v_maxcyc);
    if (v_maxcyc > 0) {
      set_max_vode_subcycles(v_maxcyc);
    }

    reaction_map.resize(numReactions());
    FORT_GET_REACTION_MAP(reaction_map.dataPtr());
    reaction_rev_map.resize(numReactions());
    for (int i=0; i<reaction_map.size(); ++i) {
      reaction_rev_map[reaction_map[i]] = i;
    }
    number_of_species = numSpecies();
}

void
ChemDriver::set_verbose_vode()
{
    FORT_SETVERBOSEVODE();
}

void
ChemDriver::set_max_vode_subcycles(int maxcyc)
{
    FORT_SETVODESUBCYC(&maxcyc);
}

void
ChemDriver::set_species_Yscales(const std::string& scalesFile)
{
    Array<int> file = encodeStringForFortran(scalesFile);
    int len         = file.size();
    FORT_SETSPECSCALY(file.dataPtr(),&len);
}

void
ChemDriver::getSpeciesNames()
{
    int max_len = FORT_GETCKMAXNAMELEN();
    int* coded = new int[max_len];
    int Nspec = FORT_GETCKNUMSPEC();
    mSpeciesNames.clear();
    mSpeciesNames.resize(Nspec);
    for (int i=0; i<Nspec; ++i)
    {
        int ifort = i+1;
	int len = FORT_GETCKSPECNAME(&ifort, coded);
	mSpeciesNames[i] = decodeStringFromFortran(coded,len);
    }
    delete [] coded;
}

void
ChemDriver::getElementNames()
{
    const int max_len = 3;
    int* coded = new int[max_len];
    int Nelt = FORT_GETCKNUMELT();
    
    mElementNames.clear();
    mElementNames.resize(Nelt);
    for (int i=0; i<Nelt; ++i)
    {
        int ifort = i+1;
	int len = FORT_GETCKELTNAME(&ifort, coded);
	mElementNames[i] = decodeStringFromFortran(coded,len);
    }
    delete [] coded;
}

void
ChemDriver::getStoichCoeffs()
{
    mNu.resize(mMaxspec * mMaxreac);
    FORT_SETNU(mNu.dataPtr(),mNu.size());
}

Array<int>
ChemDriver::reactionsWithXonL(const std::string& specName) const
{
    const int idx = index(specName) + 1;
    int Nreacs = -1;
    Array<int> reactions(numReactions());
    FORT_FINDLHS(reactions.dataPtr(),&Nreacs,&idx);
    reactions.resize(Nreacs);
    for (int i=0; i<reactions.size(); ++i)
        reactions[i]--;
    return reactions;
}

Array<int>
ChemDriver::reactionsWithXonR(const std::string& specName) const
{
    const int idx = index(specName) + 1;
    int Nreacs = -1;
    Array<int> reactions(numReactions());
    FORT_FINDRHS(reactions.dataPtr(),&Nreacs,&idx);
    reactions.resize(Nreacs);
    for (int i=0; i<reactions.size(); ++i)
        reactions[i]--;
    return reactions;
}

int
ChemDriver::numberOfElementXinSpeciesY(const std::string& eltX,
                                          const std::string& spcY) const
{
    const int eltID = indexElt(eltX);
    const int spID  = index(spcY);
    return FORT_CKELTXINSPY(&eltID,&spID);
}

Array<std::pair<std::string,int> >
ChemDriver::specCoeffsInReactions(int reacIdx) const
{
    // TODO(rgrout) This may be less error prone if the reacIdx was the
    // reaction index in the chemkin input file instead of the array offset
    // I have applied the mapping in another function that uses this,
    // so if the mapping is put in here fix reactionStringBuild also
    Array<int> KI(mMaxsp);
    Array<int> NU(mMaxsp);
    int Nids = 0;
    const int fortReacIdx = reacIdx + 1;
    FORT_CKINU(&Nids,KI.dataPtr(),&mMaxsp,NU.dataPtr(),&mMaxsp,&fortReacIdx,mNu.dataPtr());
    Array<std::pair<std::string,int> > result(Nids);
    for (int i=0; i<Nids; ++i)
        result[i] = std::pair<std::string,int>(mSpeciesNames[KI[i]-1],NU[i]);
    return result;
}

std::string
ChemDriver::reactionString(int reacIdx) const
{
    const int max_len = 72;
    const int fortReacIdx = reacIdx + 1;
    int* coded = new int[max_len];
    int len = FORT_CKSYMR(&fortReacIdx, coded);
    std::string line = decodeStringFromFortran(coded,len);
    delete [] coded;
    return line;
}

/*
 * ChemDriver::reactionStringBuild
 * Build a reaction rate string from participating species
 * need this because FORT_CKSYMR isn't always implemented
 */
std::string
ChemDriver::reactionStringBuild(int reacIdx) const
{
    std::stringstream sline;
    Array<std::pair<std::string,int> > spcCoeffs; 

    const Array<int> rxnmap = reactionMap();
    int rxn_array_index = rxnmap[reacIdx];
    spcCoeffs = specCoeffsInReactions(rxn_array_index);
    Array<std::pair<std::string,int> > reactants, products; 

    for (int j = 0; j < spcCoeffs.size(); ++j) {
        if( spcCoeffs[j].second < 0 ) {
            reactants.push_back(spcCoeffs[j]);
        }
        else {
            products.push_back(spcCoeffs[j]);
        }
    }
    for (int j = 0; j < reactants.size(); j++) {
        if (reactants[j].second < -1) {
            sline << -1*reactants[j].second << reactants[j].first;
        }
        else {
            sline << reactants[j].first;
        }
        if (j < reactants.size()-1) {
            sline << " + ";
        }
    }
    sline << " = ";
    for (int j = 0; j < products.size(); j++) {
        if (products[j].second > 1) {
            sline << products[j].second << products[j].first;
        }
        else {
            sline << products[j].first;
        }
        if (j < products.size()-1) {
            sline << " + ";
        }
    }
    return sline.str();
}


/*
 * ChemDriver::printReactions
 * Dump out a table of reactions in this mechanism
 */
void
ChemDriver::printReactions() const
{
    int nr = numReactions();
    for (int j = 0; j < nr; ++j) {
        std::cout << "(R" << j << ") " << reactionStringBuild(j) << std::endl;
    }
}


Array<Real>
ChemDriver::speciesMolecWt() const
{
    Array<Real> mwt(numSpecies());
    FORT_GETCKMWT(mwt.dataPtr());
    return mwt;
}

Array<Real>
ChemDriver::elementAtomicWt() const
{
    Array<Real> awt(numElements());
    CD_MWT(awt.dataPtr());
    return awt;
}

extern "C" {
  void CD_MWT(Real* mwt)
  {
    if (!initialized) {
      BoxLib::Abort("Must construct the ChemDriver prior to calling CD_MWT");
    }
    FORT_GETCKMWT(mwt);
  }
  void CD_XTY_PT(const Real* X, Real* Y)
  {
    if (!initialized) {
      BoxLib::Abort("Must construct the ChemDriver prior to calling CD_XTY_PT");
    }
    static Array<int> idx(BL_SPACEDIM,0);
    static int* p = idx.dataPtr();
    FORT_MOLETOMASS(p, p, X, ARLIM(p), ARLIM(p), Y, ARLIM(p), ARLIM(p));
  }

  void CD_YTX_PT(const Real* Y, Real* X)
  {
    if (!initialized) {
      BoxLib::Abort("Must construct the ChemDriver prior to calling CD_YTX_PT");
    }
    static Array<int> idx(BL_SPACEDIM,0);
    static int* p = idx.dataPtr();
    FORT_MASSTOMOLE(p, p, Y, ARLIM(p), ARLIM(p), X, ARLIM(p), ARLIM(p));
  }
}

Array<Real>
ChemDriver::massFracToMoleFrac(const Array<Real>& Y) const
{
    int nc = Y.size();
    BL_ASSERT(nc = numSpecies());
    Array<Real> X(nc);
    CD_YTX_PT(Y.dataPtr(), X.dataPtr());
    return X;
}

Array<Real>
ChemDriver::moleFracToMassFrac(const Array<Real>& X) const
{
    int nc = X.size();
    BL_ASSERT(nc==numSpecies());
    Array<Real> Y(nc);
    CD_XTY_PT(X.dataPtr(),Y.dataPtr());
    return Y;
}

void
ChemDriver::normalizeMassFrac(FArrayBox&       Ynorm,
                                 const FArrayBox& Y,
                                 const std::string&   excessSpecies,
                                 const Box&       box,
                                 int              sCompY,
                                 int              sCompYn) const
{
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompYn+numSpecies() <= Ynorm.nComp());
    
    const Box& mabx = Y.box();
    const Box& mobx = Ynorm.box();
    
    const Box& ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;

    int xsID = index(excessSpecies) + 1;
    BL_ASSERT(xsID > 0);
    FORT_NORMMASS(ovlp.loVect(), ovlp.hiVect(), &xsID,
                  Y.dataPtr(sCompY),      ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                  Ynorm.dataPtr(sCompYn), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

int
ChemDriver::numReactions() const
{
    return FORT_GETCKNUMREAC();
}

void
ChemDriver::molarProduction(FArrayBox&       Q,
			       const std::string&   specName,
			       const FArrayBox& C,
			       const FArrayBox& T,
			       const Box&       box,
			       int              sCompC,
			       int              sCompT,
			       int              sCompQ) const
{
    const int Nreacs = numReactions();
    BL_ASSERT(Q.nComp() >= sCompQ + Nreacs);
    BL_ASSERT(C.nComp() >= sCompC + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = C.box();
    const Box& mbbx = T.box();
    const Box& mobx = Q.box();
    const Box& ovlp = box & mabx & mbbx & mobx;

    const int idx = index(specName) + 1; // to fortran indexing

    FORT_MOLPROD(ovlp.loVect(), ovlp.hiVect(), &idx,
		 Q.dataPtr(sCompQ), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
		 C.dataPtr(sCompC), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		 T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()));

}

void
ChemDriver::fwdRevReacRatesGivenXTP(FArrayBox&        FwdK,
                                       FArrayBox&        RevK,
                                       const Array<int>& rxnIDs,
                                       const FArrayBox&  X,
                                       const FArrayBox&  T,
                                       Real              Patm,
                                       const Box&        box,
                                       int               sCompX,
                                       int               sCompT,
                                       int               sCompFwdK,
                                       int               sCompRevK) const
{
    const int Nspec = numSpecies();
    const int Nreac = rxnIDs.size();
    BL_ASSERT(FwdK.nComp() >= sCompFwdK+Nreac);
    BL_ASSERT(RevK.nComp() >= sCompRevK+Nreac);
    BL_ASSERT(X.nComp() >= sCompX + Nspec);
    BL_ASSERT(T.nComp() > sCompT);
    
    const Box& mabx = X.box();
    const Box& mbbx = T.box();
    const Box& moabx = FwdK.box();
    const Box& mobbx = RevK.box();
    
    const Box& ovlp = box & mabx & mbbx & moabx & mobbx;
    if( ! ovlp.ok() ) return;
    
    FORT_FRrateXTP(ovlp.loVect(), ovlp.hiVect(),
                   X.dataPtr(sCompX),        ARLIM(mabx.loVect()),  ARLIM(mabx.hiVect()),
                   T.dataPtr(sCompT),        ARLIM(mbbx.loVect()),  ARLIM(mbbx.hiVect()),
                   FwdK.dataPtr(sCompFwdK), ARLIM(moabx.loVect()), ARLIM(moabx.hiVect()),
                   RevK.dataPtr(sCompRevK), ARLIM(mobbx.loVect()), ARLIM(mobbx.hiVect()),
                   &Patm,rxnIDs.dataPtr(),&Nreac);
}
    
void
ChemDriver::heatRelease(FArrayBox&       Q,
                           const FArrayBox& Y,
                           const FArrayBox& T,
                           Real             Patm,
                           const Box&       box,
                           int              sCompY,
                           int              sCompT,
                           int              sCompQ) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(Q.nComp() > sCompQ);
    BL_ASSERT(Y.nComp() >= sCompY + Nspec);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mobx = Q.box();
    
    const Box& ovlp = box & mabx & mbbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_HTRLS(ovlp.loVect(), ovlp.hiVect(),
               Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
               T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
               Q.dataPtr(sCompQ), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
               &Patm);
}

void
ChemDriver::reactionRateY(FArrayBox&       Ydot,
                             const FArrayBox& Y,
                             const FArrayBox& T,
                             Real             Patm,
                             const Box&       box,
                             int              sCompY,
                             int              sCompT,
                             int              sCompYdot) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(Ydot.nComp() >= sCompYdot + Nspec);
    BL_ASSERT(Y.nComp() >= sCompY + Nspec);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mobx = Ydot.box();
    
    const Box& ovlp = box & mabx & mbbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_RRATEY(ovlp.loVect(), ovlp.hiVect(),
                Y.dataPtr(sCompY),       ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                T.dataPtr(sCompT),       ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
                Ydot.dataPtr(sCompYdot), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
                &Patm);
}

#ifdef LMC_SDC
void
ChemDriver::reactionRateRhoY(FArrayBox&       RhoYdot,
                             const FArrayBox& RhoY,
                             const FArrayBox& RhoH,
                             const FArrayBox& T,
                             const Box&       box,
                             int              sCompRhoY,
                             int              sCompRhoH,
                             int              sCompT,
                             int              sCompRhoYdot) const
{
    const int Nspec = numSpecies();
    BL_ASSERT(RhoYdot.nComp() >= sCompRhoYdot + Nspec);
    BL_ASSERT(RhoY.nComp() >= sCompRhoY + Nspec);
    BL_ASSERT(RhoH.nComp() > sCompRhoH);
    BL_ASSERT(T.nComp() > sCompT);

    const Box& mabx = RhoY.box();
    const Box& mbbx = RhoH.box();
    const Box& mcbx = T.box();
    const Box& mobx = RhoYdot.box();
    
    const Box& ovlp = box & mabx & mbbx & mcbx & mobx;
    if( ! ovlp.ok() ) return;
    
    FORT_RRATERHOY(ovlp.loVect(), ovlp.hiVect(),
                   RhoY.dataPtr(sCompRhoY),       ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
                   RhoH.dataPtr(sCompRhoH),       ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
                   T.dataPtr(sCompT),             ARLIM(mcbx.loVect()), ARLIM(mcbx.hiVect()),
                   RhoYdot.dataPtr(sCompRhoYdot), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()) );
}
#endif

void
ChemDriver::massFracToMoleFrac(FArrayBox&       X,
				  const FArrayBox& Y,
				  const Box&       box,
				  int              sCompY,
				  int              sCompX) const
{
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    
    const Box& mabx = Y.box();
    const Box& mobx = X.box();
    
    const Box& ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSTOMOLE(ovlp.loVect(), ovlp.hiVect(),
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::moleFracToMassFrac(FArrayBox&       Y,
				  const FArrayBox& X,
				  const Box&       box,
				  int              sCompX,
				  int              sCompY) const
{
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    
    const Box& mobx = X.box();
    const Box& mabx = Y.box();
    
    const Box& ovlp = box & mobx & mabx;
    if (!ovlp.ok())
	return;
    
    FORT_MOLETOMASS(ovlp.loVect(), ovlp.hiVect(),
		    X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()),
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()));
}

void
ChemDriver::massFracToMolarConc(FArrayBox&       C,
				   const FArrayBox& Y,
				   const FArrayBox& T,
				   Real             Patm,
				   const Box&       box,
				   int              sCompY,
				   int              sCompT,
				   int              sCompC) const
{
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompT+1 <= T.nComp());
    
    const Box& mobx = C.box();
    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    
    const Box& ovlp = box & mabx & mbbx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSTP_TO_CONC(ovlp.loVect(), ovlp.hiVect(), &Patm,
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
		    C.dataPtr(sCompC), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::massFracToMolarConc(FArrayBox&       C,
				   const FArrayBox& Y,
				   const FArrayBox& T,
				   const FArrayBox& Rho,
				   const Box&       box,
				   int              sCompR,
				   int              sCompY,
				   int              sCompT,
				   int              sCompC) const
{
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Y.nComp());
    BL_ASSERT(sCompT+1 <= T.nComp());
    BL_ASSERT(sCompR <= Rho.nComp());

    const Box& mobx = C.box();
    const Box& mabx = Y.box();
    const Box& mbbx = T.box();
    const Box& mrbx = Rho.box();

    const Box& ovlp = box & mabx & mbbx & mobx & mrbx;
    if (!ovlp.ok())
	return;
    
    FORT_MASSR_TO_CONC(ovlp.loVect(), ovlp.hiVect(), 
		    Y.dataPtr(sCompY), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		    T.dataPtr(sCompT), ARLIM(mbbx.loVect()), ARLIM(mbbx.hiVect()),
		    Rho.dataPtr(sCompR),ARLIM(mrbx.loVect()),ARLIM(mrbx.hiVect()),
		    C.dataPtr(sCompC), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

void
ChemDriver::molarConcToMoleFrac(FArrayBox&       X,
                                   const FArrayBox& C,
                                   const Box&       box,
                                   int              sCompC,
                                   int              sCompX) const
{
    BL_ASSERT(sCompX+numSpecies() <= X.nComp());
    BL_ASSERT(sCompC+numSpecies() <= C.nComp());
    
    const Box& mobx = X.box();
    const Box& mabx = C.box();
    
    const Box& ovlp = box & mabx & mobx;
    if (!ovlp.ok())
	return;
    
    FORT_CONC_TO_MOLE(ovlp.loVect(), ovlp.hiVect(),
		      C.dataPtr(sCompC), ARLIM(mabx.loVect()), ARLIM(mabx.hiVect()),
		      X.dataPtr(sCompX), ARLIM(mobx.loVect()), ARLIM(mobx.hiVect()));
}

Array<int>
ChemDriver::encodeStringForFortran(const std::string& astr)
{
    long length = astr.size();
    Array<int> result(length);
    for (int i = 0; i < length; ++i)
        result[i] = astr[i];
    return result;
}

std::string
ChemDriver::decodeStringFromFortran(const int* coded,
				       int        length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}

#include "iostream"
using std::cout;
using std::endl;
bool
ChemDriver::solveTransient(FArrayBox&        Ynew,
                           FArrayBox&        Tnew,
                           const FArrayBox&  Yold,
                           const FArrayBox&  Told,
                           FArrayBox&        FuncCount,
                           const Box&        box,
                           int               sCompY,
                           int               sCompT,
                           Real              dt,
                           Real              Patm,
                           FArrayBox*        chemDiag,
                           bool              use_stiff_solver) const
{
    BL_ASSERT(sCompY+numSpecies() <= Ynew.nComp());
    BL_ASSERT(sCompY+numSpecies() <= Yold.nComp());
    BL_ASSERT(sCompT < Tnew.nComp());
    BL_ASSERT(sCompT < Told.nComp());
    
    BL_ASSERT(Ynew.box().contains(box) && Yold.box().contains(box));
    BL_ASSERT(Tnew.box().contains(box) && Told.box().contains(box));

    const int do_diag  = (chemDiag!=0);
    Real*     diagData = do_diag ? chemDiag->dataPtr() : 0;
    const int do_stiff = (use_stiff_solver);
    int success = FORT_CONPSOLV(box.loVect(), box.hiVect(),
				Ynew.dataPtr(sCompY), ARLIM(Ynew.loVect()), ARLIM(Ynew.hiVect()),
				Tnew.dataPtr(sCompT), ARLIM(Tnew.loVect()), ARLIM(Tnew.hiVect()),
				Yold.dataPtr(sCompY), ARLIM(Yold.loVect()), ARLIM(Yold.hiVect()),
				Told.dataPtr(sCompT), ARLIM(Told.loVect()), ARLIM(Told.hiVect()),
				FuncCount.dataPtr(),
				ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
				&Patm, &dt, diagData, &do_diag, &do_stiff);
    return success > 0;
}

#ifdef LMC_SDC
bool
ChemDriver::solveTransient_sdc(FArrayBox&        rhoYnew,
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
                               bool              use_stiff_solver) const
{
    BL_ASSERT(sComprhoY+numSpecies() <= rhoYnew.nComp());
    BL_ASSERT(sComprhoY+numSpecies() <= rhoYold.nComp());
    BL_ASSERT(sComprhoH < rhoHnew.nComp());
    BL_ASSERT(sComprhoH < rhoHold.nComp());
    BL_ASSERT(sCompT    < Tnew.nComp());
    BL_ASSERT(sCompT    < Told.nComp());
    
    BL_ASSERT(rhoYnew.box().contains(box) && rhoYold.box().contains(box));
    BL_ASSERT(rhoHnew.box().contains(box) && rhoHold.box().contains(box));
    BL_ASSERT(Tnew.box().contains(box)    && Told.box().contains(box));

    const int do_diag  = (chemDiag!=0);
    Real*     diagData = do_diag ? chemDiag->dataPtr() : 0;
    const int do_stiff = (use_stiff_solver);

    int success = FORT_CONPSOLV_SDC(box.loVect(), box.hiVect(),
				    rhoYnew.dataPtr(sComprhoY), ARLIM(rhoYnew.loVect()),   ARLIM(rhoYnew.hiVect()),
				    rhoHnew.dataPtr(sComprhoH), ARLIM(rhoHnew.loVect()),   ARLIM(rhoHnew.hiVect()),
				    Tnew.dataPtr(sCompT),       ARLIM(Tnew.loVect()),      ARLIM(Tnew.hiVect()),
				    rhoYold.dataPtr(sComprhoY), ARLIM(rhoYold.loVect()),   ARLIM(rhoYold.hiVect()),
				    rhoHold.dataPtr(sComprhoH), ARLIM(rhoHold.loVect()),   ARLIM(rhoHold.hiVect()),
				    Told.dataPtr(sCompT),       ARLIM(Told.loVect()),      ARLIM(Told.hiVect()),
				    const_src.dataPtr(0),       ARLIM(const_src.loVect()), ARLIM(const_src.hiVect()),
				    FuncCount.dataPtr(),        ARLIM(FuncCount.loVect()), ARLIM(FuncCount.hiVect()),
				    &dt, diagData, &do_diag, &do_stiff);
    return success > 0;
}
#endif

void
ChemDriver::getMixAveragedRhoDiff(FArrayBox&       rhoD,
                                     const FArrayBox& Y,
                                     const FArrayBox& T,
                                     Real             Patm,
                                     const Box&       box,
                                     int              sCompY,
                                     int              sCompT,
                                     int              sCompRD) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(rhoD.nComp() >= sCompRD+numSpecies());
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(rhoD.box().contains(box));

    int do_temp    = 0;
    int do_VelVisc = 0;

    FORT_MIXAVG_RHODIFF_TEMP(box.loVect(), box.hiVect(),
                             rhoD.dataPtr(sCompRD),ARLIM(rhoD.loVect()),ARLIM(rhoD.hiVect()),
                             T.dataPtr(sCompT),    ARLIM(T.loVect()),   ARLIM(T.hiVect()),
                             Y.dataPtr(sCompY),    ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                             &Patm,&do_temp,&do_VelVisc);


}

void
ChemDriver::getMixShearVisc(FArrayBox&       eta,
                               const FArrayBox& Y,
                               const FArrayBox& T,
                               const Box&       box,
                               int              sCompY,
                               int              sCompT,
                               int              sCompE) const
{
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    BL_ASSERT(eta.nComp() >= sCompE);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(eta.box().contains(box));
    
    FORT_MIX_SHEAR_VISC(box.loVect(), box.hiVect(),
                        eta.dataPtr(sCompE),ARLIM(eta.loVect()),ARLIM(eta.hiVect()),
                        T.dataPtr(sCompT),  ARLIM(T.loVect()),  ARLIM(T.hiVect()),
                        Y.dataPtr(sCompY),  ARLIM(Y.loVect()),  ARLIM(Y.hiVect()));
}

void
ChemDriver::getRhoGivenPTY(FArrayBox&       Rho,
			      Real             Patm,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompT,
			      int              sCompY,
			      int              sCompR) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_RHOfromPTY(box.loVect(), box.hiVect(),
		    Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
		    &Patm);
}

void
ChemDriver::getRhoGivenPvTY(FArrayBox&       Rho,
			      const FArrayBox& P,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompP,
			      int              sCompT,
			      int              sCompY,
			      int              sCompR) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(P.nComp() > sCompP);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(P.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_RHOfromPvTY(box.loVect(), box.hiVect(),
		    Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
		    P.dataPtr(sCompP),   ARLIM(P.loVect()),   ARLIM(P.hiVect()));
}

void
ChemDriver::getPGivenRTY(FArrayBox&       p,
			    const FArrayBox& Rho,
			    const FArrayBox& T,
			    const FArrayBox& Y,
			    const Box&       box,
			    int              sCompR,
			    int              sCompT,
			    int              sCompY,
			    int              sCompP) const
{
    BL_ASSERT(p.nComp() > sCompP);
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(p.box().contains(box));
    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_PfromRTY(box.loVect(), box.hiVect(),
		  p.dataPtr(sCompP),   ARLIM(p.loVect()),   ARLIM(p.hiVect()),
		  Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		  T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		  Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()));
}

void
ChemDriver::getTGivenPRY(FArrayBox&       T,
                         Real             Patm,
                         const FArrayBox& Rho,
                         const FArrayBox& Y,
                         const Box&       box,
                         int              sCompR,
                         int              sCompY,
                         int              sCompT) const
{
    BL_ASSERT(Rho.nComp() > sCompR);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(Rho.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_TfromPRY(box.loVect(), box.hiVect(),
		  T.dataPtr(sCompT),   ARLIM(T.loVect()),   ARLIM(T.hiVect()),
		  Rho.dataPtr(sCompR), ARLIM(Rho.loVect()), ARLIM(Rho.hiVect()),
		  Y.dataPtr(sCompY),   ARLIM(Y.loVect()),   ARLIM(Y.hiVect()),
                  &Patm);
}

void
ChemDriver::getCpmixGivenTY(FArrayBox&       cpmix,
			       const FArrayBox& T,
			       const FArrayBox& Y,
			       const Box&       box,
			       int              sCompT,
			       int              sCompY,
			       int              sCompCp) const
{
    BL_ASSERT(cpmix.nComp() > sCompCp);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(cpmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_CPMIXfromTY(box.loVect(), box.hiVect(),
		     cpmix.dataPtr(sCompCp),ARLIM(cpmix.loVect()), ARLIM(cpmix.hiVect()),
		     T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
		     Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getCvmixGivenTY(FArrayBox&       cvmix,
			       const FArrayBox& T,
			       const FArrayBox& Y,
			       const Box&       box,
			       int              sCompT,
			       int              sCompY,
			       int              sCompCv) const
{
    BL_ASSERT(cvmix.nComp() > sCompCv);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(cvmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_CVMIXfromTY(box.loVect(), box.hiVect(),
		     cvmix.dataPtr(sCompCv),ARLIM(cvmix.loVect()), ARLIM(cvmix.hiVect()),
		     T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
		     Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getHmixGivenTY(FArrayBox&       hmix,
			      const FArrayBox& T,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompT,
			      int              sCompY,
			      int              sCompH) const
{
    BL_ASSERT(hmix.nComp() > sCompH);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(hmix.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_HMIXfromTY(box.loVect(), box.hiVect(),
		    hmix.dataPtr(sCompH),ARLIM(hmix.loVect()), ARLIM(hmix.hiVect()),
		    T.dataPtr(sCompT),   ARLIM(T.loVect()),    ARLIM(T.hiVect()),
		    Y.dataPtr(sCompY),   ARLIM(Y.loVect()),    ARLIM(Y.hiVect()));
}

void
ChemDriver::getMwmixGivenY(FArrayBox&       mwmix,
			      const FArrayBox& Y,
			      const Box&       box,
			      int              sCompY,
			      int              sCompMw) const
{
    BL_ASSERT(mwmix.nComp() > sCompMw);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());

    BL_ASSERT(mwmix.box().contains(box));
    BL_ASSERT(Y.box().contains(box));
    
    FORT_MWMIXfromY(box.loVect(), box.hiVect(),
		    mwmix.dataPtr(sCompMw),ARLIM(mwmix.loVect()), ARLIM(mwmix.hiVect()),
		    Y.dataPtr(sCompY),     ARLIM(Y.loVect()),     ARLIM(Y.hiVect()));
}

void
ChemDriver::getCpGivenT(FArrayBox&       cp,
			   const FArrayBox& T,
			   const Box&       box,
			   int              sCompT,
			   int              sCompCp) const
{
    BL_ASSERT(cp.nComp() >= sCompCp + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    BL_ASSERT(cp.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    
    FORT_CPfromT(box.loVect(), box.hiVect(),
		 cp.dataPtr(sCompCp), ARLIM(cp.loVect()), ARLIM(cp.hiVect()),
		 T.dataPtr(sCompT),   ARLIM(T.loVect()),  ARLIM(T.hiVect()));
}

void
ChemDriver::getHGivenT(FArrayBox&       h,
			  const FArrayBox& T,
			  const Box&       box,
			  int              sCompT,
			  int              sCompH) const
{
    BL_ASSERT(h.nComp() >= sCompH + numSpecies());
    BL_ASSERT(T.nComp() > sCompT);

    BL_ASSERT(h.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    
    FORT_HfromT(box.loVect(), box.hiVect(),
		h.dataPtr(sCompH), ARLIM(h.loVect()), ARLIM(h.hiVect()),
		T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()));
}

int
ChemDriver::getTGivenHY(FArrayBox&       T,
                        const FArrayBox& H,
                        const FArrayBox& Y,
                        const Box&       box,
                        int              sCompH,
                        int              sCompY,
                        int              sCompT,
                        const Real&      errMAX) const
{
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(H.nComp() > sCompH);
    BL_ASSERT(Y.nComp() >= sCompY+numSpecies());
    
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(H.box().contains(box));
    BL_ASSERT(Y.box().contains(box));

    Real solveTOL = (errMAX<0 ? mHtoTerrMAX : errMAX);
    
    return FORT_TfromHY(box.loVect(), box.hiVect(),
			T.dataPtr(sCompT), ARLIM(T.loVect()), ARLIM(T.hiVect()),
			H.dataPtr(sCompH), ARLIM(H.loVect()), ARLIM(H.hiVect()),
			Y.dataPtr(sCompY), ARLIM(Y.loVect()), ARLIM(Y.hiVect()),
			&solveTOL,&mHtoTiterMAX,mTmpData.dataPtr());
}

void
ChemDriver::getElementMoles(FArrayBox&       C_elt,
                               const std::string&   name,
                               const FArrayBox& C,
                               const Box&       box,
                               int              sCompC,
                               int              sCompC_elt) const
{
    BL_ASSERT(C.nComp() >= sCompC+numSpecies());
    BL_ASSERT(C_elt.nComp() > sCompC_elt);

    Array<int> name_enc = encodeStringForFortran(name);
    const int name_len = name_enc.size();

    FORT_GETELTMOLES(name_enc.dataPtr(), &name_len,
                     box.loVect(), box.hiVect(),
                     C_elt.dataPtr(sCompC_elt),
                     ARLIM(C_elt.loVect()),ARLIM(C_elt.hiVect()),
                     C.dataPtr(sCompC),ARLIM(C.loVect()),ARLIM(C.hiVect()));
}

Real
ChemDriver::getRuniversal() const
{
    return FORT_RUNIV();
}

Real
ChemDriver::getP1atm_MKS() const
{
  return FORT_P1ATMMKS();
}

void
ChemDriver::getOTradLoss_TDF(FArrayBox&       Qloss,
                                const FArrayBox& T,
                                const FArrayBox& X,
                                const Real       Patm,
                                const Real       T_bg,
                                const Box&       box,
                                int              sCompX,
                                int              sCompT,
                                int              sCompQ) const
{
    BL_ASSERT(Qloss.nComp() > sCompQ);
    BL_ASSERT(T.nComp() > sCompT);
    BL_ASSERT(X.nComp() >= sCompX + numSpecies());

    BL_ASSERT(Qloss.box().contains(box));
    BL_ASSERT(T.box().contains(box));
    BL_ASSERT(X.box().contains(box));

    FORT_OTrad_TDF(box.loVect(), box.hiVect(),
                   Qloss.dataPtr(sCompQ), ARLIM(Qloss.loVect()), ARLIM(Qloss.hiVect()),
                   T.dataPtr(sCompT),     ARLIM(T.loVect()),     ARLIM(T.hiVect()),
                   X.dataPtr(sCompX),     ARLIM(X.loVect()),     ARLIM(X.hiVect()),
                   &Patm, &T_bg);
}

ChemDriver::Edge::Edge (const std::string& n1,
			const std::string& n2,
			const Array<std::pair<int,Real> > rwl,
			const ChemDriver* CD)
  : sp1(n1), sp2(n2), RWL(rwl), cd(CD) {}

ChemDriver::Edge::Edge (const std::string& n1,
			const std::string& n2,
			int reac,
			Real weight,
			const ChemDriver* CD)
  : sp1(n1), sp2(n2), cd(CD) { RWL.push_back(std::pair<int,Real>(reac,weight)); }

int
ChemDriver::Edge::equivSign (const Edge& rhs) const
{   
    BL_ASSERT(cd == rhs.CD());
    BL_ASSERT(cd != 0);
    if ( (sp1 == rhs.sp1) && (sp2 == rhs.sp2) )
        return 1;
    
    else
        if ( (sp1 == rhs.sp2) && (sp2 == rhs.sp1) )
            return -1;
    
    return 0;
}

void
ChemDriver::Edge::combine (const Edge& rhs, const int sgn)
{
    BL_ASSERT(cd == rhs.CD());
    BL_ASSERT(cd != 0);
    if (sgn!=0)
    {            
        int oldSize = RWL.size();
        int delSize = rhs.RWL.size();
        RWL.resize(oldSize + delSize);
        for (int i=0; i<delSize; ++i)
            RWL[oldSize+i] = std::pair<int,Real>(rhs.RWL[i].first,sgn*rhs.RWL[i].second);
    }
}

bool
ChemDriver::Edge::touchesSp(const std::string& rhs) const
{
    return ( (sp1 == rhs) || (sp2 == rhs) );
}

void
ChemDriver::Edge::reverse()
{
    for (int i=0; i<RWL.size(); ++i)
        RWL[i].second = - RWL[i].second;
    std::swap(sp1,sp2);
}

const Array<std::pair<int,Real> >&
ChemDriver::Edge::rwl () const
{
    return RWL;
}

const std::string
ChemDriver::Edge::left() const
{
    return sp1;
}

const std::string
ChemDriver::Edge::right() const
{
    return sp2;
}

bool
ChemDriver::Edge::operator<(const Edge& rhs) const
{
    return (this->left() == rhs.left() ?  this->right() < rhs.right() : this->left() < rhs.left());
}

std::ostream& operator<< (std::ostream& os, const ChemDriver::Edge& e)
{
  const ChemDriver* cd = e.CD();
  BL_ASSERT(cd != 0);
  const Array<int>& rvmap = cd->reactionReverseMap();
  os << e.sp1 << " " << e.sp2 << " ";
  for (int i=0; i<e.RWL.size(); ++i) {
    const std::pair<int,Real>& p=e.RWL[i];
    os << rvmap[p.first] << ":" << p.second;
    if (i<e.RWL.size()-1) os << " ";
  }
  return os;
}

ChemDriver::Group::Group (const std::map<std::string,int>& eltCnts)
{
    mEltCnts = eltCnts;
}

ChemDriver::Group::Group (const ChemDriver::Group& rhs)
{
    mEltCnts = rhs.mEltCnts;
}

ChemDriver::Group
ChemDriver::Group::operator- (const ChemDriver::Group& rhs) const
{
    ChemDriver::Group val(*this);
    for (std::map<std::string,int>::const_iterator rit = rhs.mEltCnts.begin();rit!=rhs.mEltCnts.end();++rit)
    {
        std::map<std::string,int>::iterator fit = val.mEltCnts.find( rit->first );
        if ( fit != val.mEltCnts.end() )
        {
            int res = fit->second - rit->second;
            if (res != 0)
                fit->second = res;
            else
                val.mEltCnts.erase(fit);
        }
        else
            val.mEltCnts[rit->first] = -rit->second;
    }
    return val;
}

ChemDriver::Group
ChemDriver::Group::operator* (int rhs) const
{
    ChemDriver::Group val(*this);
    for (std::map<std::string,int>::iterator vit = val.mEltCnts.begin();vit!=val.mEltCnts.end();++vit)
        vit->second = rhs*vit->second;
    return val;
}

int
ChemDriver::Group::operator[](const std::string& id) const
{
    std::map<std::string,int>::const_iterator fit = mEltCnts.find( id );
    if (fit!=mEltCnts.end())
        return fit->second;
    return 0;
}

bool
ChemDriver::Group::sameSign () const
{
    if (mEltCnts.size()>0)
    {
        std::map<std::string,int>::const_iterator it = mEltCnts.begin();
        if (it->second < 0) {
            for (++it; it!=mEltCnts.end(); ++it)
                if (it->second > 0)
                    return false;
        }
        else
        {
            for (++it; it!=mEltCnts.end(); ++it)
                if (it->second < 0)
                    return false;
        }
    }
    return true;
}

bool
ChemDriver::Group::contains (const std::string& id) const
{
    return mEltCnts.find(id)!=mEltCnts.end();
}
    
Real
ChemDriver::Group::awt ()
{
    if (AtomicWeight.size() == 0)
        FillAtomicWeights(); // read up static atomic weight data
    Real myWt = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it)
    {
        Real thisAwt = AtomicWeight[it->first];
        myWt += std::abs(it->second) * thisAwt;
    }
    return myWt;
}

int
ChemDriver::Group::size() const
{
    int mySize = 0;
    for (std::map<std::string,int>::const_iterator it = mEltCnts.begin(); it!=mEltCnts.end(); ++it)
        mySize += std::abs(it->second);
    return mySize;
}

ChemDriver::Group operator* (int n, const ChemDriver::Group& g)
{
    return g.operator*(n);
}

std::map<std::string,Real> ChemDriver::Group::AtomicWeight; // Static initializer

void
ChemDriver::Group::FillAtomicWeights ()
{
    AtomicWeight["H"] = 1.00797;
    AtomicWeight["HE"] = 4.00260;
    AtomicWeight["LI"] = 6.93900;
    AtomicWeight["BE"] = 9.01220;
    AtomicWeight["B"] = 10.81100;
    AtomicWeight["C"] = 12.01115;
    AtomicWeight["N"] = 14.00670;
    AtomicWeight["O"] = 15.99940;
    AtomicWeight["F"] = 18.99840;
    AtomicWeight["NE"] = 20.18300;
    AtomicWeight["NA"] = 22.98980;
    AtomicWeight["MG"] = 24.31200;
    AtomicWeight["AL"] = 26.98150;
    AtomicWeight["SI"] = 28.08600;
    AtomicWeight["P"] = 30.97380;
    AtomicWeight["S"] = 32.06400;
    AtomicWeight["CL"] = 35.45300;
    AtomicWeight["AR"] = 39.94800;
    AtomicWeight["K"] = 39.10200;
    AtomicWeight["CA"] = 40.08000;
    AtomicWeight["SC"] = 44.95600;
    AtomicWeight["TI"] = 47.90000;
    AtomicWeight["V"] = 50.94200;
    AtomicWeight["CR"] = 51.99600;
    AtomicWeight["MN"] = 54.93800;
    AtomicWeight["FE"] = 55.84700;
    AtomicWeight["CO"] = 58.93320;
    AtomicWeight["NI"] = 58.71000;
    AtomicWeight["CU"] = 63.54000;
    AtomicWeight["ZN"] = 65.37000;
    AtomicWeight["GA"] = 69.72000;
    AtomicWeight["GE"] = 72.59000;
    AtomicWeight["AS"] = 74.92160;
    AtomicWeight["SE"] = 78.96000;
    AtomicWeight["BR"] = 79.90090;
    AtomicWeight["KR"] = 83.80000;
    AtomicWeight["RB"] = 85.47000;
    AtomicWeight["SR"] = 87.62000;
    AtomicWeight["Y"] = 88.90500;
    AtomicWeight["ZR"] = 91.22000;
    AtomicWeight["NB"] = 92.90600;
    AtomicWeight["MO"] = 95.94000;
    AtomicWeight["TC"] = 99.00000;
    AtomicWeight["RU"] = 101.07000;
    AtomicWeight["RH"] = 102.90500;
    AtomicWeight["PD"] = 106.40000;
    AtomicWeight["AG"] = 107.87000;
    AtomicWeight["CD"] = 112.40000;
    AtomicWeight["IN"] = 114.82000;
    AtomicWeight["SN"] = 118.69000;
    AtomicWeight["SB"] = 121.75000;
    AtomicWeight["TE"] = 127.60000;
    AtomicWeight["I"] = 126.90440;
    AtomicWeight["XE"] = 131.30000;
    AtomicWeight["CS"] = 132.90500;
    AtomicWeight["BA"] = 137.34000;
    AtomicWeight["LA"] = 138.91000;
    AtomicWeight["CE"] = 140.12000;
    AtomicWeight["PR"] = 140.90700;
    AtomicWeight["ND"] = 144.24000;
    AtomicWeight["PM"] = 145.00000;
    AtomicWeight["SM"] = 150.35000;
    AtomicWeight["EU"] = 151.96000;
    AtomicWeight["GD"] = 157.25000;
    AtomicWeight["TB"] = 158.92400;
    AtomicWeight["DY"] = 162.50000;
    AtomicWeight["HO"] = 164.93000;
    AtomicWeight["ER"] = 167.26000;
    AtomicWeight["TM"] = 168.93400;
    AtomicWeight["YB"] = 173.04000;
    AtomicWeight["LU"] = 174.99700;
    AtomicWeight["HF"] = 178.49000;
    AtomicWeight["TA"] = 180.94800;
    AtomicWeight["W"] = 183.85000;
    AtomicWeight["RE"] = 186.20000;
    AtomicWeight["OS"] = 190.20000;
    AtomicWeight["IR"] = 192.20000;
    AtomicWeight["PT"] = 195.09000;
    AtomicWeight["AU"] = 196.96700;
    AtomicWeight["HG"] = 200.59000;
    AtomicWeight["TL"] = 204.37000;
    AtomicWeight["PB"] = 207.19000;
    AtomicWeight["BI"] = 208.98000;
    AtomicWeight["PO"] = 210.00000;
    AtomicWeight["AT"] = 210.00000;
    AtomicWeight["RN"] = 222.00000;
    AtomicWeight["FR"] = 223.00000;
    AtomicWeight["RA"] = 226.00000;
    AtomicWeight["AC"] = 227.00000;
    AtomicWeight["TH"] = 232.03800;
    AtomicWeight["PA"] = 231.00000;
    AtomicWeight["U"] = 238.03000;
    AtomicWeight["NP"] = 237.00000;
    AtomicWeight["PU"] = 242.00000;
    AtomicWeight["AM"] = 243.00000;
    AtomicWeight["CM"] = 247.00000;
    AtomicWeight["BK"] = 249.00000;
    AtomicWeight["CF"] = 251.00000;
    AtomicWeight["ES"] = 254.00000;
    AtomicWeight["FM"] = 253.00000;
    AtomicWeight["D"] = 002.01410;
    AtomicWeight["E"] = 5.48578E-4;
}

std::ostream& operator<< (std::ostream& os, const ChemDriver::Group& g)
{
    os << "Group < ";
    for (std::map<std::string,int>::const_iterator it=g.mEltCnts.begin(); it!=g.mEltCnts.end(); ++it)
        os << it->first << ":" << it->second << " ";
    os << ">";
    return os;
}

std::list<ChemDriver::Edge>
ChemDriver::getEdges (const std::string& trElt, int PrintVerbose, int HackSplitting) const
{
    std::list<Edge> edges;
    std::map<std::string,Group> groups;
    const Array<std::string>& spNames = speciesNames();
    const Array<std::string>& eltNames = elementNames();
    for (int i=0; i<spNames.size(); ++i)
    {
        const std::string sp = spNames[i];
        std::map<std::string,int> eltCnts;
        for (int j=0; j<eltNames.size(); ++j)
        {
            const std::string elt = eltNames[j];
            int m = numberOfElementXinSpeciesY(elt,sp);
            if (m>0)
                eltCnts[elt] = m;

        }
        groups[sp] = Group(eltCnts);
    }
    for (int r=0; r<numReactions(); ++r)
    {
        const Array<std::pair<std::string,int> >& coeffs = specCoeffsInReactions(r);

        std::map<std::string,int> net;
        for (int i=0; i<coeffs.size(); ++i)
        {
            const std::string& sp=coeffs[i].first;
            const int co=coeffs[i].second;
            std::map<std::string,int>::iterator it = net.find(sp);
            if (it!=net.end())
            {
                net[sp] += co;
                if (net[sp]==0)
                    net.erase(it);
            }
            else
            {
                net[sp] = co;
            }
        }

        std::map<std::string,int> reac, prod;

        for (std::map<std::string,int>::const_iterator it = net.begin(); it!=net.end(); ++it)
        {
            const std::string& sp = it->first;
            const int n = it->second;
            if (n < 0)
            {
                if (numberOfElementXinSpeciesY(trElt,sp)>0)
                    reac[sp] = -n;
            }
            else
            {
                if (numberOfElementXinSpeciesY(trElt,sp)>0)
                    prod[sp] = n;
            }
        }

        int LR = reac.size();
        int LP = prod.size();
        
        if (LR==0 || LP==0)  // no tr-containing species in reac
            continue;
        
	int rr = reaction_rev_map[r];
	if (PrintVerbose>0) {
	  cout << rr+1 << ": " << reactionStringBuild(rr) << endl;
	}
        if ((LR == 1) || (LP == 1)) // exactly one tr-containing species either side
        {
            for (std::map<std::string,int>::const_iterator rit=reac.begin(); rit!=reac.end(); ++rit)
            {
                const std::string& spcr = rit->first;
                const int cor = rit->second;

                for (std::map<std::string,int>::const_iterator pit=prod.begin(); pit!=prod.end(); ++pit)
                {
                    const std::string& spcp = pit->first;
                    const int cop = pit->second;
                    int w = std::min(cor*groups[spcr][trElt],cop*groups[spcp][trElt]);
		    if (PrintVerbose>0) {
		      cout << "    " << spcr << " => " << spcp << " w: " << w << endl;
		    }
                    edges.push_back(Edge(spcr,spcp,r,w,this));
                }
            }
            continue;
        }

        if ( (LR==2) && (LP==2) )  // two tr-containing species each side
        {
            std::map<std::string,int>::const_iterator r0 = reac.begin();
            std::map<std::string,int>::const_iterator r1 = r0; r1++;
            std::map<std::string,int>::const_iterator p0 = prod.begin();
            std::map<std::string,int>::const_iterator p1 = p0; p1++;
            
            const std::string& rs0 = r0->first;
            const std::string& rs1 = r1->first;
            const std::string& ps0 = p0->first;
            const std::string& ps1 = p1->first;
            
            int rc0 = r0->second;
            int rc1 = r1->second;
            int pc0 = p0->second;
            int pc1 = p1->second;

            Group b0(pc0 * groups[ps0] - rc0 * groups[rs0]);
            Group b1(pc1 * groups[ps1] - rc0 * groups[rs0]);
            int pick = 0;

            // HACK
            if ((HackSplitting==1) && (trElt=="H") && (rr==61)) // Intended for GRI30
            {
                pick = 0;
                if (PrintVerbose>0) {
                    cout << "decomposition of reaction " << rr+1 << "...resorting to HACK!!" << endl;
		}
            }
            else
            {
                if (b0.sameSign() && b1.sameSign())
                {
                    if (b1.size() < b0.size())
                    {
                        pick = 1;
                        if (PrintVerbose>0) {
			  cout << "decomposition of reaction " << rr+1 << " resorting to atom cnt" << endl;
			}
                    }
                    if (b1.size() == b0.size())
                    {
                        if (b0.awt() > b1.awt())
                        {
                            pick = 1;
                            if (PrintVerbose>0) {
			      cout << "decomposition of reaction " << rr+1 << " ...resorting to tot awt" << endl;
			    }
                        }
                        else if ( (PrintVerbose>0) && (b0.awt() < b1.awt()) )
                        {
                            cout << "decomposition of reaction " << rr+1 << " ...resorting to tot awt" << endl;
                        }
                    }
                }
                else
                {
		  if (b1.sameSign()) {
		    pick = 1;
		  }
                }
            }
            
            int nR0 = rc0*groups[rs0][trElt];
            int nR1 = rc1*groups[rs1][trElt];
            int nP0 = pc0*groups[ps0][trElt];
            int nP1 = pc1*groups[ps1][trElt];
            
            if (pick == 0)
            {
	      if (PrintVerbose>0) {
		cout << "  choosing to move " << b0 << " rather than " << b1 << endl;
	        cout << "    " << rs0 << " => " << ps0 << " w: " << std::min(nR0,nP0) << endl;
	      }
	      edges.push_back(Edge(rs0,ps0,r,std::min(nR0,nP0),this));
	      if (nP0<nR0) {
		if (PrintVerbose>0) {
		  cout << "    " << rs0 << " => " << ps1 << " w: " << nR0-nP0 << endl;
		}
		edges.push_back(Edge(rs0,ps1,r,nR0-nP0,this));
	      }
              
	      if (PrintVerbose>0) {
	        cout << "    " << rs1 << " => " << ps1 << " w: " << std::min(nR1,nP1) << endl;
	      }
	      edges.push_back(Edge(rs1,ps1,r,std::min(nR1,nP1),this));
	      if (nR0<nP0) {
		if (PrintVerbose>0) {
		  cout << "    " << rs1 << " => " << ps0 << " w: " << nP0-nR0 << endl;
		}
		edges.push_back(Edge(rs1,ps0,r,nP0-nR0,this));
	      }
            }
            else
            {
	      if (PrintVerbose>0) {
		cout << "  choosing to move " << b1 << " rather than " << b0 << endl;
	        cout << "    " << rs0 << " => " << ps1 << " w: " << std::min(nR0,nP1) << endl;
	      }
	      edges.push_back(Edge(rs0,ps1,r,std::min(nR0,nP1),this));
	      if (nP1<nR0) {
		if (PrintVerbose>0) {
		  cout << "    " << rs0 << " => " << ps0 << " w: " << nR0-nP1 << endl;
		}
		edges.push_back(Edge(rs0,ps0,r,nR0-nP1,this));
	      }

	      if (PrintVerbose>0) {
	        cout << "    " << rs1 << " => " << ps0 << " w: " << std::min(nR1,nP0) << endl;
	      }
	      edges.push_back(Edge(rs1,ps0,r,std::min(nR1,nP0),this));
	      if (nR0<nP1) {
		if (PrintVerbose>0) {
		  cout << "    " << rs1 << " => " << ps1 << " w: " << nP1-nR0 << endl;
		}
		edges.push_back(Edge(rs1,ps1,r,nP1-nR0,this));
	      }
            }
            continue;
        }

        if (LR==2 && LP==3)
        {
            std::map<std::string,int> reace; //to count the number of elements in given reactant species sp
            
            for (std::map<std::string,int>::const_iterator it = reac.begin(); it!=reac.end(); ++it)
            {
              const std::string& sp = it->first;
            
              for (int j=0; j<eltNames.size(); ++j)
              {
                const std::string elt = eltNames[j];
                int m = numberOfElementXinSpeciesY(elt,sp);
                reace[sp] += m;
              }
               // cout <<"CALCULATED VALUE=" << reace[sp] << endl;
                        
            }
            std::map<std::string,int>::const_iterator re0 = reace.begin();
            std::map<std::string,int>::const_iterator re1 = re0; re1++;

            int rce0 = re0->second;
            int rce1 = re1->second;

            std::map<std::string,int>::const_iterator r0 = reac.begin();
            std::map<std::string,int>::const_iterator r1 = r0; r1++;
            std::map<std::string,int>::const_iterator p0 = prod.begin();
            std::map<std::string,int>::const_iterator p1 = p0; p1++;
            std::map<std::string,int>::const_iterator p2 = p1; p2++;

            const std::string& rs0 = r0->first;
            const std::string& rs1 = r1->first;
            const std::string& ps0 = p0->first;
            const std::string& ps1 = p1->first;
            const std::string& ps2 = p2->first;
            int rc0 = r0->second;
            int rc1 = r1->second;
            int pc0 = p0->second;
            int pc1 = p1->second;
            int pc2 = p2->second;

            if(rce0<rce1)
            {
              Group b0(pc0 * groups[ps0] - rc0 * groups[rs0]);
              //cout << b0.awt() << endl;
              Group b1(pc1 * groups[ps1] - rc0 * groups[rs0]);
              //cout << b1.awt() << endl;
              Group b2(pc2 * groups[ps2] - rc0 * groups[rs0]);
              //cout << b2.awt() << endl;
              if(b0.awt()<b1.awt() && b0.awt()<b2.awt())
              {
                // link r0 to p0
                    int w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps0 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs1,ps2,r,w,this));
              }
              else if(b1.awt()<b0.awt() && b1.awt()<b2.awt())
              {
                // link r0 to p1
                    int w = std::min(rc0*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps1 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs1,ps2,r,w,this));
              }
              else
              {
                // link r0 to p2
                    int w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps2 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps2,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs1,ps0,r,w,this));
              }
            }
           else
            {
              Group b0(pc0 * groups[ps0] - rc1 * groups[rs1]);
              //cout << b0.awt() << endl;
              Group b1(pc1 * groups[ps1] - rc1 * groups[rs1]);
              //cout << b1.awt() << endl;
              Group b2(pc2 * groups[ps2] - rc1 * groups[rs1]);
              //cout << b2.awt() << endl;
                
              if(b0.awt()<b1.awt() && b0.awt()<b2.awt())
              {
                // link r1 to p0
                    int w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps0 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs0,ps2,r,w,this));
              }
              else if(b1.awt()<b0.awt() && b1.awt()<b2.awt())
              {
                // link r1 to p1
                    int w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps1 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs0,ps2,r,w,this));
              }
              else
              {
                // link r1 to p2
                    int w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps2 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps2,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs0,ps0,r,w,this));
              }
           }

           continue;
         }
         if (LR==2 && LP==4)
        {
            std::map<std::string,int> reace; //to count the number of elements in given reactant species sp
            
            for (std::map<std::string,int>::const_iterator it = reac.begin(); it!=reac.end(); ++it)
            {
              const std::string& sp = it->first;
            
              for (int j=0; j<eltNames.size(); ++j)
              {
                const std::string elt = eltNames[j];
                int m = numberOfElementXinSpeciesY(elt,sp);
                reace[sp] += m;
              }
               // cout <<"CALCULATED VALUE=" << reace[sp] << endl;
                        
            }
            std::map<std::string,int>::const_iterator re0 = reace.begin();
            std::map<std::string,int>::const_iterator re1 = re0; re1++;

            int rce0 = re0->second;
            int rce1 = re1->second;

            std::map<std::string,int>::const_iterator r0 = reac.begin();
            std::map<std::string,int>::const_iterator r1 = r0; r1++;
            std::map<std::string,int>::const_iterator p0 = prod.begin();
            std::map<std::string,int>::const_iterator p1 = p0; p1++;
            std::map<std::string,int>::const_iterator p2 = p1; p2++;
            std::map<std::string,int>::const_iterator p3 = p2; p3++;

            const std::string& rs0 = r0->first;
            const std::string& rs1 = r1->first;
            const std::string& ps0 = p0->first;
            const std::string& ps1 = p1->first;
            const std::string& ps2 = p2->first;
            const std::string& ps3 = p3->first;
            int rc0 = r0->second;
            int rc1 = r1->second;
            int pc0 = p0->second;
            int pc1 = p1->second;
            int pc2 = p2->second;
            int pc3 = p3->second;

            if(rce0<rce1)
            {
              Group b0(pc0 * groups[ps0] - rc0 * groups[rs0]);
              Group b1(pc1 * groups[ps1] - rc0 * groups[rs0]);
              Group b2(pc2 * groups[ps2] - rc0 * groups[rs0]);
              Group b3(pc3 * groups[ps3] - rc0 * groups[rs0]);             

              if(b0.awt()<b1.awt() && b0.awt()<b2.awt() && b0.awt()<b3.awt())
              {
                // link r0 to p0
                    int w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps0 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs1,ps2,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs1,ps3,r,w,this));
              }
              else if(b1.awt()<b0.awt() && b1.awt()<b2.awt() && b1.awt()<b3.awt())
              {
                // link r0 to p1
                    int w = std::min(rc0*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps1 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs1,ps2,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs1,ps3,r,w,this));
              }
              else if(b2.awt()<b0.awt() && b2.awt()<b1.awt() && b2.awt()<b3.awt())
              {
                // link r0 to p2
                    int w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps2 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps2,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs1,ps3,r,w,this));
              }
              else
              {
                // link r0 to p3
                    int w = std::min(rc0*groups[rs0][trElt],pc3*groups[ps3][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs0 << " => " << ps3 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs0,ps3,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs1,ps2,r,w,this));
              }
            }
           else
            {
              Group b0(pc0 * groups[ps0] - rc1 * groups[rs1]);
              Group b1(pc1 * groups[ps1] - rc1 * groups[rs1]);
              Group b2(pc2 * groups[ps2] - rc1 * groups[rs1]);
              Group b3(pc3 * groups[ps3] - rc1 * groups[rs1]);
                
              if(b0.awt()<b1.awt() && b0.awt()<b2.awt() && b0.awt()<b3.awt())
              {
                // link r1 to p0
                    int w = std::min(rc1*groups[rs1][trElt],pc0*groups[ps0][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps0 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps0,r,w,this));
                   w = std::min(rc1*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs0,ps2,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs0,ps3,r,w,this));
              }
              else if(b1.awt()<b0.awt() && b1.awt()<b2.awt() && b1.awt()<b3.awt())
              {
                // link r1 to p1
                    int w = std::min(rc1*groups[rs1][trElt],pc1*groups[ps1][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps1 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs0,ps2,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs0,ps3,r,w,this));
              }
              else if(b2.awt()<b0.awt() && b2.awt()<b1.awt() && b2.awt()<b3.awt())
              {
                // link r1 to p2
                    int w = std::min(rc1*groups[rs1][trElt],pc2*groups[ps2][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps2 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps2,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc3*groups[ps3][trElt]);
                   edges.push_back(Edge(rs0,ps3,r,w,this));
              }
             else
              {
                // link r1 to p3
                    int w = std::min(rc1*groups[rs1][trElt],pc3*groups[ps3][trElt]);
                    if (PrintVerbose>0) {
                      cout << "    " << rs1 << " => " << ps3 << " w: " << w << endl;
                    }
                    edges.push_back(Edge(rs1,ps3,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc1*groups[ps1][trElt]);
                   edges.push_back(Edge(rs0,ps1,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc0*groups[ps0][trElt]);
                   edges.push_back(Edge(rs0,ps0,r,w,this));
                   w = std::min(rc0*groups[rs0][trElt],pc2*groups[ps2][trElt]);
                   edges.push_back(Edge(rs0,ps2,r,w,this));
              }
           }

           continue;
         }
       

	cout << "Cannot decompose rxn: " << rr+1 << ": " << reactionStringBuild(rr) << endl;
    }

    // Uniquify/combine edges
    std::list<Edge> uEdges;
    while (edges.size()!=0)
    {
        const Edge oe(edges.front()); edges.pop_front();
        bool foundIt = false;
        // There's probably a faster way to search this list...
        for (std::list<Edge>::iterator it=uEdges.begin(); !foundIt && it!=uEdges.end(); ++it)
        {
            int sgn = it->equivSign(oe);
            if (sgn!=0)
            {
                it->combine(oe,sgn);
                foundIt = true;
            }
        }
        if (!foundIt)
            uEdges.push_back(oe);
    }
    
    return uEdges;
}
