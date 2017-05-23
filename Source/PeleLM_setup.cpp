//
// Note: define TEMPERATURE if you want variables T and rho*h, h = c_p*T,in the 
//       State_Type part of the state
//       whose component indices are Temp and RhoH
//       whose evoution equations are 
//       \pd (rho h)/\pd t + diver (\rho U h) = div k grad T
//       rho DT/dt = div k/c_p grad T
//       define RADIATION only if TEMPERATURE is also defined and the evolution equations
//       for T and rho*h are
//       \pd (rho h)/\pd t + diver (\rho U h) = div k grad T - div q_rad
//       rho DT/dt = div k/c_p grad T - 1/c_p div q_rad
//
//       The equation for temperature is used solely for computing div k grad T
//
//       Note: The reasons for using an auxiliary T equations are two fold:
//             1) solving the C-N difference equations for 
//                \pd (rho h)/\pd t + diver (\rho U h) = div k/c_p grad h
//                does not easily fit into our framework as of 10/31/96
//             2) the boundary condition for rho h at a wall is ill defined
//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
//
// see variableSetUp on how to use or not use these types in the state
//

#include <algorithm>
#include <cstdio>
#include <iomanip>

#include <PeleLM.H>
#include <RegType.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ErrorList.H>
#include <Prob_F.H>
#include <DERIVE_F.H>
#include <AMReX_FArrayBox.H>
#include <NAVIERSTOKES_F.H>
#include <PeleLM_F.H>
#include <AMReX_Utility.H>
#include <NS_error_F.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)      \
  const int* fablo = (fab).loVect();            \
  const int* fabhi = (fab).hiVect();            \
  Real* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)     \
  const int* fablo = (fab).loVect();            \
  const int* fabhi = (fab).hiVect();            \
  const Real* fabdat = (fab).dataPtr();

static Box the_same_box (const Box& b)    { return b;                 }
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }
static Box the_nodes (const Box& b) { return amrex::surroundingNodes(b); }

static bool do_group_bndry_fills = false;

//
// Components are  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall.
//

static
int
norm_vel_bc[] =
{
//  INT_DIR, EXT_DIR, HOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, EXT_DIR, EXT_DIR
};

static
int
tang_vel_bc[] =
{
//  INT_DIR, EXT_DIR, HOEXTRAP, REFLECT_EVEN, HOEXTRAP, EXT_DIR
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, HOEXTRAP, EXT_DIR
};

static
int
scalar_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

//
// Here we use slipwall to impose neumann condition, not to represent a real wall
//
static
int
temp_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
press_bc[] =
{
  INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
};

static
int
rhoh_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR
};

static
int
divu_bc[] =
{
  INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
dsdt_bc[] =
{
  INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
reflect_bc[] =
{
  INT_DIR,REFLECT_EVEN,REFLECT_EVEN,REFLECT_EVEN,REFLECT_EVEN,REFLECT_EVEN
};

static
int
species_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_x_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,norm_vel_bc[lo_bc[0]]);
  bc.setHi(0,norm_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,norm_vel_bc[lo_bc[1]]);
  bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_z_vel_bc (BCRec&       bc,
              const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
  bc.setLo(2,norm_vel_bc[lo_bc[2]]);
  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

static
void
set_scalar_bc (BCRec&       bc,
               const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,scalar_bc[lo_bc[i]]);
    bc.setHi(i,scalar_bc[hi_bc[i]]);
  }
}

static
void
set_reflect_bc (BCRec&       bc,
                const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,reflect_bc[lo_bc[i]]);
    bc.setHi(i,reflect_bc[hi_bc[i]]);
  }
}

static
void
set_temp_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,temp_bc[lo_bc[i]]);
    bc.setHi(i,temp_bc[hi_bc[i]]);
  }
}

static
void
set_pressure_bc (BCRec&       bc,
                 const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  int i;
  for (i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,press_bc[lo_bc[i]]);
    bc.setHi(i,press_bc[hi_bc[i]]);
  }
}

static
void
set_rhoh_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,rhoh_bc[lo_bc[i]]);
    bc.setHi(i,rhoh_bc[hi_bc[i]]);
  }
}

static
void
set_divu_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,divu_bc[lo_bc[i]]);
    bc.setHi(i,divu_bc[hi_bc[i]]);
  }
}

static
void
set_dsdt_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,dsdt_bc[lo_bc[i]]);
    bc.setHi(i,dsdt_bc[hi_bc[i]]);
  }
}

static
void
set_species_bc (BCRec&       bc,
                const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,species_bc[lo_bc[i]]);
    bc.setHi(i,species_bc[hi_bc[i]]);
  }
}

typedef StateDescriptor::BndryFunc BndryFunc;

extern "C"
{
//
// Function called by BCRec for user-supplied boundary data.
//
  typedef void (*ChemBndryFunc_FortBndryFunc)(Real* data, ARLIM_P(lo), ARLIM_P(hi),
                                              const int* dom_lo, const int* dom_hi,
                                              const Real* dx, const Real* grd_lo,
                                              const Real* time, const int* bc,
                                              const int* stateID);
}

class ChemBndryFunc
  :
  public BndryFunc
{
public:
  //
  // Bogus constructor.
  //
  ChemBndryFunc()
    :
    m_func(0),
    m_stateID(-1) {}
  //
  // Constructor.
  //
  ChemBndryFunc (ChemBndryFunc_FortBndryFunc  inFunc,
                 const std::string& stateName)
    :
    m_func(inFunc),
    m_stateName(stateName)
    {
      m_stateID = getStateID(m_stateName);
      BL_ASSERT(m_stateID >= 0);
    }
  //
  // Another Constructor which sets "regular" and "group" fill routines..
  //
  ChemBndryFunc (ChemBndryFunc_FortBndryFunc  inFunc,
                 const std::string&           stateName,
                 BndryFuncDefault             gFunc)
    :
    BndryFunc(gFunc,gFunc),
    m_func(inFunc),
    m_stateName(stateName)
    {
      m_stateID = getStateID(m_stateName);
      BL_ASSERT(m_stateID >= 0);
    }
  //
  // Destructor.
  //
  virtual ~ChemBndryFunc ()  override {}
    
  virtual StateDescriptor::BndryFunc* clone () const override
    {
      //
      // Bitwise copy ok here, no copy ctr required.
      //
      return new ChemBndryFunc(*this);
    }
  //
  // Fill boundary cells the "regular" way.
  // The other virtual function in BndryFunc will
  // give us the appropriate call for "group" fills.
  //
  virtual void operator () (Real* data, const int* lo, const int* hi,
                            const int* dom_lo, const int* dom_hi,
                            const Real* dx, const Real* grd_lo,
                            const Real* time, const int* bc) const override
    {
      BL_ASSERT(m_func != 0);
      m_func(data,ARLIM(lo),ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc,&m_stateID);
    }
  //
  // Fill boundary cells using "group" function.
  //
  virtual void operator () (Real* data, const int* lo, const int* hi,
                            const int* dom_lo, const int* dom_hi,
                            const Real* dx, const Real* grd_lo,
                            const Real* time, const int* bc, int ng) const override
    {
      BndryFunc::operator()(data, lo, hi, dom_lo, dom_hi, dx, grd_lo, time, bc, ng);
    }
  //
  // Access.
  //
  int getStateID () const              { return m_stateID;   }
  const std::string& getStateName () const { return m_stateName; }
  ChemBndryFunc_FortBndryFunc getBndryFunc () const  { return m_func;      }
    
protected:

  static int getStateID (const std::string& stateName)
    {
      const Array<std::string>& names = PeleLM::getChemSolve().speciesNames();
      for (int i=0; i<names.size(); i++)
        if (names[i] == stateName)
          return i;
      return -1;
    }
    
private:

  ChemBndryFunc_FortBndryFunc m_func;
  std::string                 m_stateName;
  int                         m_stateID;
};

//
// Indices of fuel and oxidizer -- ParmParsed in & used in a couple places.
//
std::string PeleLM::fuelName        = "CH4";
Array<std::string> PeleLM::consumptionName(1);
static std::string oxidizerName           = "O2";
static std::string productName            = "CO2";

void
PeleLM::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  int i;

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    phys_bc.setLo(dir,SlipWall);
    phys_bc.setHi(dir,SlipWall);
  }

  Initialize();
  BCRec bc;
  //
  // Set state variable Id's (Density and velocities set already).
  //
  int counter   = Density;
  int RhoH      = -1;
  int FirstSpec = -1;
  int Trac      = -1;
  int RhoRT     = -1;

  FirstSpec = ++counter;
  nspecies  = getChemSolve().numSpecies();
  counter  += nspecies - 1;
  RhoH = ++counter;
  Trac = ++counter;
  Temp = ++counter;
#ifndef BL_RHORT_IN_TRACER
  RhoRT = ++counter;
#endif
  NUM_STATE = ++counter;
  NUM_SCALARS = NUM_STATE - Density;

  const Array<std::string>& names = getChemSolve().speciesNames();

  amrex::Print() << nspecies << " Chemical species interpreted:\n { ";
  for (int i = 0; i < nspecies; i++)
    amrex::Print() << names[i] << ' ' << ' ';
  amrex::Print() << '}' << '\n' << '\n';

  //
  // Send indices of fuel and oxidizer to fortran for setting prob data in common block
  //
  ParmParse pp("ns");
  pp.query("fuelName",fuelName);
  consumptionName[0] = fuelName;
  if (int nc = pp.countval("consumptionName"))
  {
    consumptionName.resize(nc);
    pp.getarr("consumptionName",consumptionName,0,nc);
  }
  pp.query("oxidizerName",oxidizerName);
  pp.query("productName",productName);
  pp.query("do_group_bndry_fills",do_group_bndry_fills);

  //
  // Set scale of chemical components, used in ODE solves
  //
  std::string speciesScaleFile; pp.query("speciesScaleFile",speciesScaleFile);
  if (! speciesScaleFile.empty())
  {
    amrex::Print() << "  Setting scale values for chemical species\n\n";
    getChemSolve().set_species_Yscales(speciesScaleFile);
  }
  int verbose_vode=0; pp.query("verbose_vode",verbose_vode);
  if (verbose_vode!=0)
    getChemSolve().set_verbose_vode();

  int fuelID = getChemSolve().index(fuelName);
  int oxidID = getChemSolve().index(oxidizerName);
  int prodID = getChemSolve().index(productName);

  amrex::Print() << " fuel name " << fuelName << std::endl;
  amrex::Print() << " index for fueld and oxidizer " << fuelID << " " << oxidID << std::endl;

  FORT_SET_PROB_SPEC(&fuelID, &oxidID, &prodID, &nspecies);
  //
  // Get a species to use as a flame tracker.
  //
  std::string flameTracName = fuelName;
  pp.query("flameTracName",flameTracName);    
  //
  // **************  DEFINE VELOCITY VARIABLES  ********************
  //
  desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),StateDescriptor::Point,1,NUM_STATE,
                         &cell_cons_interp);

  Array<BCRec>       bcs(BL_SPACEDIM);
  Array<std::string> name(BL_SPACEDIM);

  set_x_vel_bc(bc,phys_bc);
  bcs[0]  = bc;
  name[0] = "x_velocity";
  desc_lst.setComponent(State_Type,Xvel,"x_velocity",bc,BndryFunc(FORT_XVELFILL));

  set_y_vel_bc(bc,phys_bc);
  bcs[1]  = bc;
  name[1] = "y_velocity";
  desc_lst.setComponent(State_Type,Yvel,"y_velocity",bc,BndryFunc(FORT_YVELFILL));

#if(BL_SPACEDIM==3)
  set_z_vel_bc(bc,phys_bc);
  bcs[2]  = bc;
  name[2] = "z_velocity";
  desc_lst.setComponent(State_Type,Zvel,"z_velocity",bc,BndryFunc(FORT_ZVELFILL));
#endif
  //
  // To enable "group" operations on filling velocities, we need to
  // overwrite the first component specifing how to do "regular"
  // and "group" fill operations.
  //
  if (do_group_bndry_fills)
  {
    desc_lst.setComponent(State_Type,
                          Xvel,
                          name,
                          bcs,
                          BndryFunc(FORT_XVELFILL,FORT_VELFILL));
  }
  //
  // **************  DEFINE SCALAR VARIABLES  ********************
  //
  // Set range of combination limit to include rho, rhoh and species, if they exist
  //
  //int combinLimit_lo = static_cast<int>(Density);
  //int combinLimit_hi = std::max(combinLimit_lo, RhoH);
  set_scalar_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Density,"density",bc,BndryFunc(FORT_DENFILL),
                        &cell_cons_interp);
  //
  // **************  DEFINE RHO*H  ********************
  //
  set_rhoh_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,RhoH,"rhoh",bc,BndryFunc(FORT_RHOHFILL),
                        &cell_cons_interp);
  //
  // **************  DEFINE TEMPERATURE  ********************
  //
  set_temp_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Temp,"temp",bc,BndryFunc(FORT_TEMPFILL));
  //
  // ***************  DEFINE SPECIES **************************
  //
  bcs.resize(nspecies);
  name.resize(nspecies);

  set_species_bc(bc,phys_bc);

  for (i = 0; i < nspecies; i++)
  {
    bcs[i]  = bc;
    name[i] = "rho.Y(" + names[i] + ")";

    desc_lst.setComponent(State_Type,
                          FirstSpec+i,
                          name[i].c_str(),
                          bc,
                          ChemBndryFunc(FORT_CHEMFILL,names[i]),
                          &cell_cons_interp);
  }
  //
  // To enable "group" operations on filling species, we need to
  // overwrite the first component specifing how to do "regular"
  // and "group" fill operations.
  //
//    if (do_group_bndry_fills)
  if (false)
  {        
    desc_lst.setComponent(State_Type,
                          FirstSpec,
                          name,
                          bcs,
                          ChemBndryFunc(FORT_CHEMFILL,names[0],FORT_ALLCHEMFILL),
                          &cell_cons_interp);
  }
  //
  // ***************  DEFINE TRACER and RhoRT **************************
  //
  // Force BCs to be REFLECT_EVEN for RhoRT ghost cells in UGRADP.
  // ADVFILL is ok for this, if all BC's are REFLECT_EVEN (ie, no EXT_DIR)
  //
  if (RhoRT >= 0)
  {
    set_scalar_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,BndryFunc(FORT_ADVFILL));

    set_reflect_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,RhoRT,"RhoRT",bc,BndryFunc(FORT_ADVFILL));
  }
  else
  {
    set_reflect_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,BndryFunc(FORT_ADVFILL));
  }

  advectionType.resize(NUM_STATE);
  diffusionType.resize(NUM_STATE);
  is_diffusive.resize(NUM_STATE);
  visc_coef.resize(NUM_STATE);
  //
  // Assume everything is diffusive and then change it if it is not.
  //
  for (i = 0; i < NUM_STATE; i++)
  {
    advectionType[i] = NonConservative;
    diffusionType[i] = RhoInverse_Laplacian_S;
    is_diffusive[i]  = true;
  }

  if (do_mom_diff == 1)
    for (int d = 0; d < BL_SPACEDIM; d++)
      advectionType[Xvel+d] = Conservative;

  is_diffusive[Density] = false;

  if (RhoRT > 0)
    is_diffusive[RhoRT] = false;

  if (Trac > 0)
  {
    if (RhoRT > 0)
    {
      advectionType[Trac] = NonConservative;
      diffusionType[Trac] = Laplacian_S;
      if (trac_diff_coef <= 0.0)
        is_diffusive[Trac] = false;
    }
    else
    {
      is_diffusive[Trac] = false;
    }
  }

  advectionType[Density] = Conservative;
  diffusionType[Density] = Laplacian_SoverRho;
  advectionType[Temp] = NonConservative;
  diffusionType[Temp] = RhoInverse_Laplacian_S;
  advectionType[RhoH] = Conservative;
  diffusionType[RhoH] = Laplacian_SoverRho;

  for (int i = 0; i < nspecies; ++i)
  {
    advectionType[FirstSpec + i] = Conservative;
    diffusionType[FirstSpec + i] = Laplacian_SoverRho;
  }

  if (is_diffusive[Density])
    amrex::Abort("PeleLM::variableSetUp(): density cannot diffuse");
  //
  // ---- pressure
  //
  // the #if 1 makes this a simpler problem CAR
#if 1
  desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                         StateDescriptor::Interval,1,1,
                         &node_bilinear_interp);

  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));
#else
  desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                         StateDescriptor::Point,1,1,
                         &node_bilinear_interp,true);

  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(FORT_PRESFILL));

  //
  // ---- time derivative of pressure
  //
  Dpdt_Type = desc_lst.size();
  desc_lst.addDescriptor(Dpdt_Type,IndexType::TheNodeType(),
                         StateDescriptor::Interval,1,1,
                         &node_bilinear_interp);
  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Dpdt_Type,Dpdt,"dpdt",bc,BndryFunc(FORT_PRESFILL));
#endif
  //
  // ---- right hand side of divergence constraint.
  //
  int ngrow;
  //
  // stick Divu_Type on the end of the descriptor list.
  //
  Divu_Type = desc_lst.size();
  ngrow = 1;
  desc_lst.addDescriptor(Divu_Type,IndexType::TheCellType(),StateDescriptor::Point,ngrow,1,
                         &cell_cons_interp);
  set_divu_bc(bc,phys_bc);
  desc_lst.setComponent(Divu_Type,Divu,"divu",bc,BndryFunc(FORT_DIVUFILL));
  //
  // Stick Dsdt_Type on the end of the descriptor list.
  //
  Dsdt_Type = desc_lst.size();
	    
  ngrow = 0;
  desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),StateDescriptor::Point,ngrow,1,
                         &cell_cons_interp);
  set_dsdt_bc(bc,phys_bc);
  desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,BndryFunc(FORT_DSDTFILL));
  //
  // Add in the fcncall tracer type quantity.
  //
  FuncCount_Type = desc_lst.size();
  desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),StateDescriptor::Point,0, 1, &cell_cons_interp);
  desc_lst.setComponent(FuncCount_Type, 0, "FuncCount", bc, BndryFunc(FORT_DQRADFILL));
  rhoydotSetUp();
  //
  // rho_temp
  //
  derive_lst.add("rho_temp",IndexType::TheCellType(),1,FORT_DERMPRHO,the_same_box);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Temp,1);
  //
  // enthalpy
  //
  derive_lst.add("enthalpy",IndexType::TheCellType(),1,FORT_DERDVRHO,the_same_box);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,RhoH,1);
  //
  // Species mass fractions.
  //
  for (i = 0; i < nspecies; i++)
  {
    const std::string name = "Y("+names[i]+")";
    derive_lst.add(name,IndexType::TheCellType(),1,FORT_DERDVRHO,the_same_box);
    derive_lst.addComponent(name,desc_lst,State_Type,Density,1);
    derive_lst.addComponent(name,desc_lst,State_Type,FirstSpec + i,1);
  }
  //
  // Species mole fractions
  //
  Array<std::string> var_names_molefrac(nspecies);
  for (i = 0; i < nspecies; i++)
    var_names_molefrac[i] = "X("+names[i]+")";
  derive_lst.add("molefrac",IndexType::TheCellType(),nspecies,
                 var_names_molefrac,FORT_DERMOLEFRAC,the_same_box);
  derive_lst.addComponent("molefrac",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("molefrac",desc_lst,State_Type,FirstSpec,nspecies);

  //
  // Species concentrations
  //
  Array<std::string> var_names_conc(nspecies);
  for (i = 0; i < nspecies; i++)
    var_names_conc[i] = "C("+names[i]+")";
  derive_lst.add("concentration",IndexType::TheCellType(),nspecies,
                 var_names_conc,FORT_DERCONCENTRATION,the_same_box);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,
                          FirstSpec,nspecies);

  if (nspecies > 0)
  {
    //
    // rho-sum rhoY.
    //
    derive_lst.add("rhominsumrhoY",IndexType::TheCellType(),1,FORT_DERRHOMINUSSUMRHOY,the_same_box);
    derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,Density,1);
    for (i = 0; i < nspecies; i++)
    {
      const int comp = FirstSpec + i;
      derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,comp,1);
    }
  }
  //
  // Sum rhoYdot
  //
  derive_lst.add("sumRhoYdot",IndexType::TheCellType(),1,FORT_DERSUMRHOYDOT,the_same_box);
  for (i = 0; i < nspecies; i++)
  {
    derive_lst.addComponent("sumRhoYdot",desc_lst,RhoYdot_Type,i,1);
  }
  //
  // **************  DEFINE DERIVED QUANTITIES ********************
  //
  // Divergence of velocity field.
  //
  derive_lst.add("diveru",IndexType::TheCellType(),1,FORT_DERMGDIVU,grow_box_by_one);
  derive_lst.addComponent("diveru",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // average pressure
  //
  derive_lst.add("avg_pressure",IndexType::TheCellType(),1,FORT_DERAVGPRES,
                 the_nodes);
  derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in X direction.
  //
  derive_lst.add("gradpx",IndexType::TheCellType(),1,FORT_DERGRDPX,the_nodes);
  derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in Y direction.
  //
  derive_lst.add("gradpy",IndexType::TheCellType(),1,FORT_DERGRDPY,the_nodes);
  derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);

#if (BL_SPACEDIM == 3)
  //
  // Pressure gradient in Z direction.
  //
  derive_lst.add("gradpz",IndexType::TheCellType(),1,FORT_DERGRDPZ,the_nodes);
  derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
#endif
  //
  // Magnitude of vorticity.
  //
  derive_lst.add("mag_vort",IndexType::TheCellType(),1,FORT_DERMGVORT,grow_box_by_one);
  derive_lst.addComponent("mag_vort",desc_lst,State_Type,Xvel,BL_SPACEDIM);
#ifdef DO_LMC_FORCE
  //
  // forcing - used to calculate the rate of injection of energy
  //
  derive_lst.add("forcing",IndexType::TheCellType(),1,FORT_DERFORCING,the_same_box);
  derive_lst.addComponent("forcing",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("forcing",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // forcex - used to put the forcing term in the plot file
  //
  derive_lst.add("forcex",IndexType::TheCellType(),1,FORT_DERFORCEX,the_same_box);
  derive_lst.addComponent("forcex",desc_lst,State_Type,Density,1);
  //    derive_lst.addComponent("forcex",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // forcey - used to put the forcing term in the plot file
  //
  derive_lst.add("forcey",IndexType::TheCellType(),1,FORT_DERFORCEY,the_same_box);
  derive_lst.addComponent("forcey",desc_lst,State_Type,Density,1);
  //    derive_lst.addComponent("forcey",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // forcez - used to put the forcing term in the plot file
  //
  derive_lst.add("forcez",IndexType::TheCellType(),1,FORT_DERFORCEZ,the_same_box);
  derive_lst.addComponent("forcez",desc_lst,State_Type,Density,1);
#endif
  //    derive_lst.addComponent("forcez",desc_lst,State_Type,Xvel,BL_SPACEDIM);

  std::string curv_str = "mean_progress_curvature";
  derive_lst.add(curv_str,IndexType::TheCellType(),1,&DeriveRec::GrowBoxByOne);
    
#ifdef PARTICLES
  //
  // The particle count at this level.
  //
  derive_lst.add("particle_count",IndexType::TheCellType(),1,
                 FORT_DERNULL,the_same_box);
  derive_lst.addComponent("particle_count",desc_lst,State_Type,Density,1);
  //
  // The total # of particles at our level or above.
  //
  derive_lst.add("total_particle_count",IndexType::TheCellType(),1,
                 FORT_DERNULL,the_same_box);
  derive_lst.addComponent("total_particle_count",desc_lst,State_Type,Density,1);
  //
  // Force all particles to be tagged.
  //
  //err_list.add("total_particle_count",1,ErrorRec::Special,FORT_PART_CNT_ERR);
#endif
  //
  // **************  DEFINE ERROR ESTIMATION QUANTITIES  *************
  //
  const int nGrowErr = 1;
  err_list.add("temp", nGrowErr, ErrorRec::Special, BL_FORT_PROC_CALL(FORT_TEMPERROR,fort_temperror));
  err_list.add("mag_vort", nGrowErr, ErrorRec::Special, BL_FORT_PROC_CALL(FORT_MVERROR,fort_mverror));
  err_list.add("tracer", nGrowErr, ErrorRec::Special, BL_FORT_PROC_CALL(FORT_ADVERROR,fort_adverror));
  //
  // Tag region of interesting chemistry.
  //
  const int idx = getChemSolve().index(flameTracName);
  if (idx >= 0)
  {
    amrex::Print() << "Flame tracer will be " << flameTracName << '\n';
    const std::string name = "Y("+flameTracName+")";
    err_list.add(name,nGrowErr,ErrorRec::Special,FORT_FLAMETRACERROR);
  }
}

class PLMBld
  :
  public LevelBld
{
  virtual void variableSetUp () override;
  virtual void variableCleanUp () override;
  virtual AmrLevel *operator() () override;
  virtual AmrLevel *operator() (Amr&            papa,
                                int             lev,
                                const Geometry& level_geom,
                                const BoxArray& ba,
                                const DistributionMapping& dm,
                                Real            time) override;
};

PLMBld HTbld;

LevelBld*
getLevelBld ()
{
  return &HTbld;
}

void
PLMBld::variableSetUp ()
{
  PeleLM::variableSetUp();
}

void
PLMBld::variableCleanUp ()
{
  PeleLM::variableCleanUp();
}

AmrLevel*
PLMBld::operator() ()
{
  return new PeleLM;
}

AmrLevel*
PLMBld::operator() (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    const DistributionMapping& dm,
                    Real            time)
{
    return new PeleLM(papa, lev, level_geom, ba, dm, time);
}

static
int
rhoydot_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_rhoydot_bc (BCRec&       bc,
                const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  for (int i = 0; i < BL_SPACEDIM; i++)
  {
    bc.setLo(i,rhoydot_bc[lo_bc[i]]);
    bc.setHi(i,rhoydot_bc[hi_bc[i]]);
  }
}

void
PeleLM::rhoydotSetUp()
{
  RhoYdot_Type       = desc_lst.size();
  const int ngrow = 1;
  const int nrhoydot = getChemSolve().numSpecies();

  amrex::Print() << "RhoYdot_Type, nrhoydot = " << RhoYdot_Type << ' ' << nrhoydot << '\n';

  desc_lst.addDescriptor(RhoYdot_Type,IndexType::TheCellType(),
                         StateDescriptor::Point,ngrow,nrhoydot,
                         &lincc_interp);
	
  //const StateDescriptor& d_cell = desc_lst[State_Type];
  const Array<std::string>& names   = getChemSolve().speciesNames();

  BCRec bc;	
  set_rhoydot_bc(bc,phys_bc);
  for (int i = 0; i < nrhoydot; i++)
  {
    const std::string name = "I_R[rhoY("+names[i]+")]";
    desc_lst.setComponent(RhoYdot_Type, i, name.c_str(), bc,
                          BndryFunc(FORT_RHOYDOTFILL), &lincc_interp, 0, nrhoydot-1);
  }
}


