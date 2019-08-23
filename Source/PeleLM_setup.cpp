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

#ifdef AMREX_USE_SUNDIALS_3x4x 
#include <actual_Creactor.h>
#endif

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
// Components are  Interior, Inflow, Outflow, Symmetry, &
// SlipWallAdiab, NoSlipWallAdiab, SlipWallIsoTherm, NoSlipWallIsoTherm.
//

static
int
norm_vel_bc[] =
{
  // 
INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD
};

static
int
tang_vel_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_ODD, REFLECT_EVEN, REFLECT_ODD
};

static
int
scalar_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
temp_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

static
int
press_bc[] =
{
  //INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP
  INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
rhoh_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

static
int
divu_bc[] =
{
  INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
dsdt_bc[] =
{
  INT_DIR, EXT_DIR, EXT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR 
};

static
int
reflect_bc[] =
{
  INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
int
species_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
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
      Vector<std::string> names;
      PeleLM::getSpeciesNames(names);
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
std::string PeleLM::productName     = "CO2";
Vector<std::string> PeleLM::consumptionName(1);
static std::string oxidizerName     = "O2";


void
PeleLM::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    phys_bc.setLo(dir,SlipWall);
    phys_bc.setHi(dir,SlipWall);
  }

  Initialize();

  // Initialize the runtime parameters for any of the external code
  init_extern();

  /* PelePhysics */
  amrex::Print() << " Initialization of network, reactor and transport \n";
  init_network();

#ifdef _OPENMP
#pragma omp parallel
#endif  
  reactor_init(&cvode_iE,&cvode_ncells);

  init_transport(use_tranlib);

  BCRec bc;
  //
  // Set state variable Id's (Density and velocities set already).
  //
  int counter   = Density;
//  int RhoH      = -1;
//  int FirstSpec = -1;
//  int Trac      = -1;
//  int RhoRT     = -1;

  first_spec = ++counter;
  pphys_get_num_spec(&nspecies);
  nreactions = pphys_numReactions();
  counter  += nspecies - 1;
  RhoH = ++counter;
  Trac = ++counter;
  Temp = ++counter;
#ifndef BL_RHORT_IN_TRACER
  RhoRT = ++counter;
#endif
  NUM_STATE = ++counter;
  NUM_SCALARS = NUM_STATE - Density;

  getSpeciesNames(spec_names); 

  amrex::Print() << nreactions << " Reactions in mechanism \n";
  amrex::Print() << nspecies << " Chemical species interpreted:\n { ";
  for (int i = 0; i < nspecies; i++)
    amrex::Print() << spec_names[i] << ' ' << ' ';
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

  // Fill spec_scalY that is not used anywhere anymore: FIXME 
  //if (! speciesScaleFile.empty())
  //{
  //  amrex::Print() << "  Setting scale values for chemical species\n\n";
  //  getChemSolve().set_species_Yscales(speciesScaleFile);
  //}
  int verbose_vode=0; pp.query("verbose_vode",verbose_vode);
  if (verbose_vode!=0)
    pphys_set_verbose_vode();

  int fuelID = getSpeciesIdx(fuelName);
  int oxidID = getSpeciesIdx(oxidizerName);
  int prodID = getSpeciesIdx(productName);
  int bathID = getSpeciesIdx("N2");

  amrex::Print() << " fuel name " << fuelName << std::endl;
  amrex::Print() << " index for fuel and oxidizer " << fuelID << " " << oxidID << std::endl;
  amrex::Print() << " index for bath " << bathID << std::endl;

  int dm = BL_SPACEDIM;
  int flag_active_control = 0;
  
  if (PeleLM::flag_active_control){
    flag_active_control = 1;}
  else {flag_active_control = 0;}
  
  set_prob_spec(dm,DefaultGeometry().ProbLo(),DefaultGeometry().ProbHi(),
                &bathID, &fuelID, &oxidID, &prodID, &nspecies,flag_active_control);
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

  Vector<BCRec>       bcs(BL_SPACEDIM);
  Vector<std::string> name(BL_SPACEDIM);

  set_x_vel_bc(bc,phys_bc);
  bcs[0]  = bc;
  name[0] = "x_velocity";
  desc_lst.setComponent(State_Type,Xvel,"x_velocity",bc,BndryFunc(xvel_fill));

  set_y_vel_bc(bc,phys_bc);
  bcs[1]  = bc;
  name[1] = "y_velocity";
  desc_lst.setComponent(State_Type,Yvel,"y_velocity",bc,BndryFunc(yvel_fill));

#if(BL_SPACEDIM==3)
  set_z_vel_bc(bc,phys_bc);
  bcs[2]  = bc;
  name[2] = "z_velocity";
  desc_lst.setComponent(State_Type,Zvel,"z_velocity",bc,BndryFunc(zvel_fill));
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
                          BndryFunc(xvel_fill,vel_fill));
  }
  //
  // **************  DEFINE SCALAR VARIABLES  ********************
  //
  // Set range of combination limit to include rho, rhoh and species, if they exist
  //
  //int combinLimit_lo = static_cast<int>(Density);
  //int combinLimit_hi = std::max(combinLimit_lo, RhoH);
  set_scalar_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Density,"density",bc,BndryFunc(den_fill),
                        &cell_cons_interp);
  //
  // **************  DEFINE RHO*H  ********************
  //
  set_rhoh_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,RhoH,"rhoh",bc,BndryFunc(rhoh_fill),
                        &cell_cons_interp);
  //
  // **************  DEFINE TEMPERATURE  ********************
  //
  set_temp_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Temp,"temp",bc,BndryFunc(temp_fill));
  //
  // ***************  DEFINE SPECIES **************************
  //
  bcs.resize(nspecies);
  name.resize(nspecies);

  set_species_bc(bc,phys_bc);

  for (int i = 0; i < nspecies; i++)
  {
    bcs[i]  = bc;
    name[i] = "rho.Y(" + spec_names[i] + ")";

    desc_lst.setComponent(State_Type,
                          first_spec+i,
                          name[i].c_str(),
                          bc,
                          ChemBndryFunc(chem_fill,spec_names[i]),
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
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,BndryFunc(adv_fill));

    set_reflect_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,RhoRT,"RhoRT",bc,BndryFunc(adv_fill));
  }
  else
  {
    set_reflect_bc(bc,phys_bc);
    desc_lst.setComponent(State_Type,Trac,"tracer",bc,BndryFunc(adv_fill));
  }

  advectionType.resize(NUM_STATE);
  diffusionType.resize(NUM_STATE);
  is_diffusive.resize(NUM_STATE);
  visc_coef.resize(NUM_STATE);
  //
  // Assume everything is diffusive and then change it if it is not.
  //
  for (int i = 0; i < NUM_STATE; i++)
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
    advectionType[first_spec + i] = Conservative;
    diffusionType[first_spec + i] = Laplacian_SoverRho;
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
  desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(press_fill));
#else
  desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                         StateDescriptor::Point,1,1,
                         &node_bilinear_interp,true);

  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,BndryFunc(press_fill));

  //
  // ---- time derivative of pressure
  //
  Dpdt_Type = desc_lst.size();
  desc_lst.addDescriptor(Dpdt_Type,IndexType::TheNodeType(),
                         StateDescriptor::Interval,1,1,
                         &node_bilinear_interp);
  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Dpdt_Type,Dpdt,"dpdt",bc,BndryFunc(press_fill));
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
  desc_lst.setComponent(Divu_Type,Divu,"divu",bc,BndryFunc(divu_fill));
  //
  // Stick Dsdt_Type on the end of the descriptor list.
  //
  Dsdt_Type = desc_lst.size();
	    
  ngrow = 0;
  desc_lst.addDescriptor(Dsdt_Type,IndexType::TheCellType(),StateDescriptor::Point,ngrow,1,
                         &cell_cons_interp);
  set_dsdt_bc(bc,phys_bc);
  desc_lst.setComponent(Dsdt_Type,Dsdt,"dsdt",bc,BndryFunc(dsdt_fill));
  //
  // Add in the fcncall tracer type quantity.
  //
  FuncCount_Type = desc_lst.size();
  desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),StateDescriptor::Point,0, 1, &cell_cons_interp);
  desc_lst.setComponent(FuncCount_Type, 0, "FuncCount", bc, BndryFunc(dqrad_fill));
  rhoydotSetUp();
  //
  // rho_temp
  //
  derive_lst.add("rho_temp",IndexType::TheCellType(),1,dermprho,the_same_box);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Temp,1);
  //
  // enthalpy
  //
  derive_lst.add("enthalpy",IndexType::TheCellType(),1,derdvrho,the_same_box);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,RhoH,1);

  //
  // Molecular Weight
  //
  derive_lst.add("molweight",IndexType::TheCellType(),1,dermolweight,the_same_box);
  derive_lst.addComponent("molweight",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("molweight",desc_lst,State_Type,first_spec,nspecies);
  
  //
  // Group Species Rho.Y (for ploting in plot file)
  //
  Vector<std::string> var_names_rhoY(nspecies);
  for (int i = 0; i < nspecies; i++)
    var_names_rhoY[i] = "rho.Y("+spec_names[i]+")";
  derive_lst.add("rhoY",IndexType::TheCellType(),nspecies,
                 var_names_rhoY,derRhoY,the_same_box);
  derive_lst.addComponent("rhoY",desc_lst,State_Type,first_spec,nspecies);
  
  //
  // Individual Species mass fractions (for error tag with tracer)
  //
  for (int i = 0; i < nspecies; i++)
  {
    const std::string chname = "Y("+spec_names[i]+")";
    derive_lst.add(chname,IndexType::TheCellType(),1,derdvrho,the_same_box);
    derive_lst.addComponent(chname,desc_lst,State_Type,Density,1);
    derive_lst.addComponent(chname,desc_lst,State_Type,first_spec + i,1);
  }
  //
  // Group Species mass fractions (for ploting in plot file)
  //
  Vector<std::string> var_names_massfrac(nspecies);
  for (int i = 0; i < nspecies; i++)
    var_names_massfrac[i] = "Y("+spec_names[i]+")";
  derive_lst.add("mass_fractions",IndexType::TheCellType(),nspecies,
                 var_names_massfrac,dermassfrac,the_same_box);
  derive_lst.addComponent("mass_fractions",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mass_fractions",desc_lst,State_Type,first_spec,nspecies);

  //
  // Species mole fractions
  //
  Vector<std::string> var_names_molefrac(nspecies);
  for (int i = 0; i < nspecies; i++)
    var_names_molefrac[i] = "X("+spec_names[i]+")";
  derive_lst.add("mole_fractions",IndexType::TheCellType(),nspecies,
                 var_names_molefrac,dermolefrac,the_same_box);
  derive_lst.addComponent("mole_fractions",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mole_fractions",desc_lst,State_Type,first_spec,nspecies);

  //
  // Species concentrations
  //
  Vector<std::string> var_names_conc(nspecies);
  for (int i = 0; i < nspecies; i++)
    var_names_conc[i] = "C("+spec_names[i]+")";
  derive_lst.add("concentration",IndexType::TheCellType(),nspecies,
                 var_names_conc,derconcentration,the_same_box);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,
                          first_spec,nspecies);

  //
  // Derive transport coefficients
  //
  Vector<std::string> var_names_transp_coeff(nspecies+2);
  for (int i = 0; i < nspecies; i++)
    var_names_transp_coeff[i] = "D_Y("+spec_names[i]+")";
  var_names_transp_coeff[nspecies] = "Lambda";
  var_names_transp_coeff[nspecies+1] = "Mu";
  derive_lst.add("cc_transport_coeffs",IndexType::TheCellType(),nspecies+2,
                 var_names_transp_coeff,dertransportcoeff,the_same_box);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,
                          first_spec,nspecies);




  if (nspecies > 0)
  {
    //
    // rho-sum rhoY.
    //
    derive_lst.add("rhominsumrhoY",IndexType::TheCellType(),1,drhomry,the_same_box);
    derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,Density,1);
    for (int i = 0; i < nspecies; i++)
    {
      const int comp = first_spec + i;
      derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,comp,1);
    }
  }
  //
  // Sum rhoYdot
  //
  derive_lst.add("sumRhoYdot",IndexType::TheCellType(),1,dsrhoydot,the_same_box);
  for (int i = 0; i < nspecies; i++)
  {
    derive_lst.addComponent("sumRhoYdot",desc_lst,RhoYdot_Type,i,1);
  }
  //
  // **************  DEFINE DERIVED QUANTITIES ********************
  //
  // Divergence of velocity field.
  //
  derive_lst.add("diveru",IndexType::TheCellType(),1,dermgdivu,grow_box_by_one);
  derive_lst.addComponent("diveru",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // average pressure
  //
  derive_lst.add("avg_pressure",IndexType::TheCellType(),1,deravgpres,
                 the_nodes);
  derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in X direction.
  //
  derive_lst.add("gradpx",IndexType::TheCellType(),1,dergrdpx,the_nodes);
  derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in Y direction.
  //
  derive_lst.add("gradpy",IndexType::TheCellType(),1,dergrdpy,the_nodes);
  derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);

#if (BL_SPACEDIM == 3)
  //
  // Pressure gradient in Z direction.
  //
  derive_lst.add("gradpz",IndexType::TheCellType(),1,dergrdpz,the_nodes);
  derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
#endif
  //
  // Magnitude of vorticity.
  //
  derive_lst.add("mag_vort",IndexType::TheCellType(),1,dermgvort,grow_box_by_one);
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
    
#ifdef AMREX_PARTICLES
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
  //err_list.add("total_particle_count",1,ErrorRec::Special,part_cnt_err);
#endif

  //
  // Dynamically generated error tagging functions
  //
  std::string amr_prefix = "amr";
  ParmParse ppamr(amr_prefix);
  Vector<std::string> refinement_indicators;
  ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
  for (int i=0; i<refinement_indicators.size(); ++i)
  {
    std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];

    ParmParse ppr(ref_prefix);
    Real min_time =  0; ppr.query("start_time",min_time);
    Real max_time = -1; ppr.query("end_time",max_time);
    int max_level = -1;  ppr.query("max_level",max_level);
    
    if (ppr.countval("value_greater")) {
      Real value; ppr.get("value_greater",value);
      std::string field; ppr.get("field_name",field);
      err_list.add(field.c_str(),0,ErrorRec::Special,
                   LM_Error_Value(valgt_error,value,min_time,max_time,max_level));
    }
    else if (ppr.countval("value_less")) {
      Real value; ppr.get("value_less",value);
      std::string field; ppr.get("field_name",field);
      err_list.add(field.c_str(),0,ErrorRec::Special,
                   LM_Error_Value(vallt_error,value,min_time,max_time,max_level));
    }
    else if (ppr.countval("vorticity_greater")) {
      Real value; ppr.get("vorticity_greater",value);
      std::string field; field="mag_vort";
      err_list.add(field.c_str(),0,ErrorRec::Special,
                   LM_Error_Value(magvort_error,value,min_time,max_time,max_level));
    }
    else if (ppr.countval("adjacent_difference_greater")) {
      Real value; ppr.get("adjacent_difference_greater",value);
      std::string field; ppr.get("field_name",field);
      err_list.add(field.c_str(),1,ErrorRec::Special,
                   LM_Error_Value(diffgt_error,value,min_time,max_time,max_level));
    }
    else if (ppr.countval("in_box_lo")) {
      std::vector<Real> box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
      ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
      ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
      RealBox realbox(&(box_lo[0]),&(box_hi[0]));
      err_list.add("dummy",1,ErrorRec::Special,
                   LM_Error_Value(box_error,realbox,min_time,max_time,max_level));
    }
    else {
      Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
    }
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
  int nrhoydot;
  pphys_get_num_spec(&nrhoydot);

  amrex::Print() << "RhoYdot_Type, nrhoydot = " << RhoYdot_Type << ' ' << nrhoydot << '\n';

  desc_lst.addDescriptor(RhoYdot_Type,IndexType::TheCellType(),
                         StateDescriptor::Point,ngrow,nrhoydot,
                         &lincc_interp);
	
  //const StateDescriptor& d_cell = desc_lst[State_Type];

  BCRec bc;	
  set_rhoydot_bc(bc,phys_bc);
  for (int i = 0; i < nrhoydot; i++)
  {
    const std::string name = "I_R[rhoY("+spec_names[i]+")]";
    desc_lst.setComponent(RhoYdot_Type, i, name.c_str(), bc,
                          BndryFunc(rhoYdot_fill), &lincc_interp, 0, nrhoydot-1);
  }
}


