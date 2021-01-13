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
#include <AMReX_FArrayBox.H>
#include <AMReX_Utility.H>
#include <EOS.H>
#include <Transport.H>

#include <PeleLM_derive.H>
#include <IndexDefines.H>
#include <reactor.h>

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
static Box grow_box_by_two (const Box& b) { return amrex::grow(b,2); }
static Box the_nodes (const Box& b) { return amrex::surroundingNodes(b); }


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

//FIXME -- Pulled from IAMR, needs update for PLM
static int norm_gradp_bc[] =
{
  INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_ODD, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int tang_gradp_bc[] =
{
  INT_DIR, FOEXTRAP, FOEXTRAP, REFLECT_EVEN, FOEXTRAP, FOEXTRAP, FOEXTRAP, FOEXTRAP
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
#if (AMREX_SPACEDIM == 3)
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
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

#if (AMREX_SPACEDIM == 3)
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
  for (int i = 0; i < AMREX_SPACEDIM; i++)
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
  for (int i = 0; i < AMREX_SPACEDIM; i++)
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
  for (int i = 0; i < AMREX_SPACEDIM; i++)
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
  for (i = 0; i < AMREX_SPACEDIM; i++)
  {
    bc.setLo(i,press_bc[lo_bc[i]]);
    bc.setHi(i,press_bc[hi_bc[i]]);
  }
}

static 
void 
set_gradpx_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_gradp_bc[lo_bc[0]]);
    bc.setHi(0,norm_gradp_bc[hi_bc[0]]);
    bc.setLo(1,tang_gradp_bc[lo_bc[1]]);
    bc.setHi(1,tang_gradp_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_gradp_bc[lo_bc[2]]);
    bc.setHi(2,tang_gradp_bc[hi_bc[2]]);
#endif
}

static
void
set_gradpy_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_gradp_bc[lo_bc[0]]);
    bc.setHi(0,tang_gradp_bc[hi_bc[0]]);
    bc.setLo(1,norm_gradp_bc[lo_bc[1]]);
    bc.setHi(1,norm_gradp_bc[hi_bc[1]]);
#if (BL_SPACEDIM == 3)
    bc.setLo(2,tang_gradp_bc[lo_bc[2]]);
    bc.setHi(2,tang_gradp_bc[hi_bc[2]]);
#endif
}

#if (BL_SPACEDIM == 3)
static
void
set_gradpz_bc (BCRec&       bc,
	       const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_gradp_bc[lo_bc[0]]);
    bc.setHi(0,tang_gradp_bc[hi_bc[0]]);
    bc.setLo(1,tang_gradp_bc[lo_bc[1]]);
    bc.setHi(1,tang_gradp_bc[hi_bc[1]]);
    bc.setLo(2,norm_gradp_bc[lo_bc[2]]);
    bc.setHi(2,norm_gradp_bc[hi_bc[2]]);
}
#endif

static
void
set_rhoh_bc (BCRec&       bc,
             const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < AMREX_SPACEDIM; i++)
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
  for (int i = 0; i < AMREX_SPACEDIM; i++)
  {
    bc.setLo(i,divu_bc[lo_bc[i]]);
    bc.setHi(i,divu_bc[hi_bc[i]]);
  }
}

static
void
set_species_bc (BCRec&       bc,
                const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < AMREX_SPACEDIM; i++)
  {
    bc.setLo(i,species_bc[lo_bc[i]]);
    bc.setHi(i,species_bc[hi_bc[i]]);
  }
}

//
// Indices of fuel and oxidizer -- ParmParsed in & used in a couple places.
//
std::string PeleLM::fuelName        = "CH4";
std::string PeleLM::productName     = "CO2";
Vector<std::string> PeleLM::consumptionName(1);
static std::string oxidizerName     = "O2";

//
// Get EB-aware interpolater when needed
//
#ifdef AMREX_USE_EB  
  static auto& cc_interp = eb_cell_cons_interp;
#else
  static auto& cc_interp = cell_cons_interp;
#endif


void
PeleLM::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  prob_parm.reset(new ProbParm{});
  ac_parm.reset(new ACParm{});

  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    phys_bc.setLo(dir,SlipWall);
    phys_bc.setHi(dir,SlipWall);
  }

  Initialize();

  amrex::Print() << " Initialization of reactor... \n";
  int reactor_type = 2;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif  
{
#ifdef USE_SUNDIALS_PP
  SetTolFactODE(relative_tol_chem,absolute_tol_chem);
#endif
#ifdef AMREX_USE_GPU
  reactor_info(reactor_type,ncells_chem);
#else
  reactor_init(reactor_type,ncells_chem);
#endif
}

  amrex::Print() << " Initialization of EOS (CPP)... \n";
  EOS::init();
  amrex::Print() << " Initialization of Transport (CPP)... \n";
  transport_init();

  BCRec bc;

  //
  // Set state variable Id's (Density and velocities set already).
  //

  first_spec = DEF_first_spec;
  RhoH = DEF_RhoH;
  Temp = DEF_Temp;
  RhoRT = DEF_RhoRT;
  NUM_STATE = DEF_NUM_STATE;
  NUM_SCALARS = DEF_NUM_SCALARS;

  EOS::speciesNames(spec_names);

  amrex::Print() << NUM_REACTIONS << " Reactions in mechanism \n";
  amrex::Print() << NUM_SPECIES << " Chemical species interpreted:\n { ";
  for (int i = 0; i < NUM_SPECIES; i++)
    amrex::Print() << spec_names[i] << ' ' << ' ';
  amrex::Print() << '}' << '\n' << '\n';

  //
  // Send indices of fuel and oxidizer to fortran for setting prob data in common block
  //
  ParmParse ppns("ns");
  ppns.query("fuelName",fuelName);
  consumptionName[0] = fuelName;
  if (int nc = ppns.countval("consumptionName"))
  {
    consumptionName.resize(nc);
    ppns.getarr("consumptionName",consumptionName,0,nc);
  }
  ppns.query("oxidizerName",oxidizerName);
  ppns.query("productName",productName);

  //
  // Set scale of chemical components, used in ODE solves
  //
  std::string speciesScaleFile; ppns.query("speciesScaleFile",speciesScaleFile);

  amrex::Print() << " fuel name " << fuelName << std::endl;
  
  //
  // Get a species to use as a flame tracker.
  //
  std::string flameTracName = fuelName;
  ppns.query("flameTracName",flameTracName);    
  //
  // **************  DEFINE VELOCITY VARIABLES  ********************
  //
  bool state_data_extrap = false;
  bool store_in_checkpoint = true;
  desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,1,NUM_STATE,
			 &cc_interp,state_data_extrap,store_in_checkpoint);

  amrex::StateDescriptor::BndryFunc pelelm_bndryfunc(pelelm_cc_ext_fill);
  pelelm_bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.


  Vector<BCRec>       bcs(AMREX_SPACEDIM);
  Vector<std::string> name(BL_SPACEDIM);

  set_x_vel_bc(bc,phys_bc);
  bcs[0]  = bc;
  name[0] = "x_velocity";

  set_y_vel_bc(bc,phys_bc);
  bcs[1]  = bc;
  name[1] = "y_velocity";

#if(AMREX_SPACEDIM==3)
  set_z_vel_bc(bc,phys_bc);
  bcs[2]  = bc;
  name[2] = "z_velocity";
#endif

  desc_lst.setComponent(State_Type,
                        Xvel,
                        name,
                        bcs,
                        pelelm_bndryfunc);

  //
  // **************  DEFINE SCALAR VARIABLES  ********************
  //
  // Set range of combination limit to include rho, rhoh and species, if they exist
  //
  set_scalar_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Density,"density",bc,pelelm_bndryfunc);

  //
  // **************  DEFINE RHO*H  ********************
  //
  set_rhoh_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,RhoH,"rhoh",bc,pelelm_bndryfunc);
  //
  // **************  DEFINE TEMPERATURE  ********************
  //
  set_temp_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,Temp,"temp",bc,pelelm_bndryfunc);
  //
  // ***************  DEFINE SPECIES **************************
  //
  bcs.resize(NUM_SPECIES);
  name.resize(NUM_SPECIES);

  set_species_bc(bc,phys_bc);

  for (int i = 0; i < NUM_SPECIES; i++)
  {
    bcs[i]  = bc;
    name[i] = "rho.Y(" + spec_names[i] + ")";
  }
  desc_lst.setComponent(State_Type,
                        first_spec,
                        name,
                        bcs,
                        pelelm_bndryfunc);
          
  //
  // ***************  DEFINE TRACER and RhoRT **************************
  //
  // Force BCs to be REFLECT_EVEN for RhoRT ghost cells in UGRADP.
  // ADVFILL is ok for this, if all BC's are REFLECT_EVEN (ie, no EXT_DIR)
  //
  set_reflect_bc(bc,phys_bc);
  desc_lst.setComponent(State_Type,RhoRT,"RhoRT",bc,pelelm_bndryfunc);

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

  advectionType[Density] = Conservative;
  diffusionType[Density] = Laplacian_SoverRho;
  advectionType[Temp] = NonConservative;
  diffusionType[Temp] = RhoInverse_Laplacian_S;
  advectionType[RhoH] = Conservative;
  diffusionType[RhoH] = Laplacian_SoverRho;

  for (int i = 0; i < NUM_SPECIES; ++i)
  {
    advectionType[first_spec + i] = Conservative;
    diffusionType[first_spec + i] = Laplacian_SoverRho;
  }

  if (is_diffusive[Density])
    amrex::Abort("PeleLM::variableSetUp(): density cannot diffuse");
  //
  // ---- pressure
  //
  desc_lst.addDescriptor(Press_Type,IndexType::TheNodeType(),
                         StateDescriptor::Interval,1,1,
                         &node_bilinear_interp);

  amrex::StateDescriptor::BndryFunc pelelm_nodal_bf(pelelm_dummy_fill);
  pelelm_nodal_bf.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.


  set_pressure_bc(bc,phys_bc);
  desc_lst.setComponent(Press_Type,Pressure,"pressure",bc,pelelm_nodal_bf);
  //
  // ---- grad P
  //
  desc_lst.addDescriptor(Gradp_Type,IndexType::TheCellType(),
			 StateDescriptor::Interval,gradp_grow,AMREX_SPACEDIM,
			 &cc_interp,state_data_extrap,store_in_checkpoint);
  amrex::StateDescriptor::BndryFunc gradp_bf(pelelm_dummy_fill);
  gradp_bf.setRunOnGPU(true);
  
  bcs.resize(BL_SPACEDIM);
  name.resize(BL_SPACEDIM);
  
  set_gradpx_bc(bc,phys_bc);
  bcs[0]  = bc;
  name[0] = "gradpx";
  
  set_gradpy_bc(bc,phys_bc);
  bcs[1]  = bc;
  name[1] = "gradpy";
  
#if(AMREX_SPACEDIM==3)
  set_gradpz_bc(bc,phys_bc);
  bcs[2]  = bc;
  name[2] = "gradpz";
#endif
  
  desc_lst.setComponent(Gradp_Type, Gradpx, name, bcs, gradp_bf);

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
                         &cc_interp);

  // EXT_DIR ghost cells should never actually get used, just need something computable 
  amrex::StateDescriptor::BndryFunc pelelm_divu_bf(pelelm_dummy_fill);
  pelelm_divu_bf.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

  set_divu_bc(bc,phys_bc);
  desc_lst.setComponent(Divu_Type,Divu,"divu",bc,pelelm_divu_bf);
  //
  // Add in the fcncall tracer type quantity.
  //
  FuncCount_Type = desc_lst.size();
  desc_lst.addDescriptor(FuncCount_Type, IndexType::TheCellType(),StateDescriptor::Point,0, 1, &cc_interp);
  
  amrex::StateDescriptor::BndryFunc pelelm_dummy_bf(pelelm_dummy_fill);
  pelelm_dummy_bf.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

  desc_lst.setComponent(FuncCount_Type, 0, "FuncCount", bc, pelelm_dummy_bf);

  rhoydotSetUp();
  //
  // rho_temp
  //
  derive_lst.add("rho_temp",IndexType::TheCellType(),1,pelelm_dermprho,the_same_box);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("rho_temp",desc_lst,State_Type,Temp,1);
  //
  // enthalpy
  //
  derive_lst.add("enthalpy",IndexType::TheCellType(),1,pelelm_derdvrho,the_same_box);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("enthalpy",desc_lst,State_Type,RhoH,1);

  //
  // Molecular Weight
  //
  derive_lst.add("molweight",IndexType::TheCellType(),1,pelelm_dermolweight,the_same_box);
  derive_lst.addComponent("molweight",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("molweight",desc_lst,State_Type,first_spec,NUM_SPECIES);
  
  //
  // Mixture heat capacity
  //
  derive_lst.add("cpmix",IndexType::TheCellType(),1,pelelm_dercpmix,the_same_box);
  derive_lst.addComponent("cpmix",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("cpmix",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("cpmix",desc_lst,State_Type,first_spec,NUM_SPECIES);
  
  //
  // Group Species Rho.Y (for ploting in plot file)
  //
  Vector<std::string> var_names_rhoY(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++)
    var_names_rhoY[i] = "rho.Y("+spec_names[i]+")";
  derive_lst.add("rhoY",IndexType::TheCellType(),NUM_SPECIES,
                 var_names_rhoY,pelelm_derRhoY,the_same_box);
  derive_lst.addComponent("rhoY",desc_lst,State_Type,first_spec,NUM_SPECIES);
  
  //
  // Individual Species mass fractions (for error tag with tracer)
  //
  for (int i = 0; i < NUM_SPECIES; i++)
  {
    const std::string chname = "Y("+spec_names[i]+")";
    derive_lst.add(chname,IndexType::TheCellType(),1,pelelm_derdvrho,the_same_box);
    derive_lst.addComponent(chname,desc_lst,State_Type,Density,1);
    derive_lst.addComponent(chname,desc_lst,State_Type,first_spec + i,1);
  }
  //
  // Group Species mass fractions (for ploting in plot file)
  //
  Vector<std::string> var_names_massfrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++)
    var_names_massfrac[i] = "Y("+spec_names[i]+")";
  derive_lst.add("mass_fractions",IndexType::TheCellType(),NUM_SPECIES,
                 var_names_massfrac,pelelm_dermassfrac,the_same_box);
  derive_lst.addComponent("mass_fractions",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mass_fractions",desc_lst,State_Type,first_spec,NUM_SPECIES);

  //
  // Species mole fractions
  //
  Vector<std::string> var_names_molefrac(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++)
    var_names_molefrac[i] = "X("+spec_names[i]+")";
  derive_lst.add("mole_fractions",IndexType::TheCellType(),NUM_SPECIES,
                 var_names_molefrac,pelelm_dermolefrac,the_same_box);
  derive_lst.addComponent("mole_fractions",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mole_fractions",desc_lst,State_Type,first_spec,NUM_SPECIES);

  //
  // Species concentrations
  //
  Vector<std::string> var_names_conc(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++)
    var_names_conc[i] = "C("+spec_names[i]+")";
  derive_lst.add("concentration",IndexType::TheCellType(),NUM_SPECIES,
                 var_names_conc,pelelm_derconcentration,the_same_box);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("concentration",desc_lst,State_Type,
                          first_spec,NUM_SPECIES);

  //
  // Derive transport coefficients
  //
  Vector<std::string> var_names_transp_coeff(NUM_SPECIES+2);
  for (int i = 0; i < NUM_SPECIES; i++)
    var_names_transp_coeff[i] = "D_Y("+spec_names[i]+")";
  var_names_transp_coeff[NUM_SPECIES] = "Lambda";
  var_names_transp_coeff[NUM_SPECIES+1] = "Mu";
  derive_lst.add("cc_transport_coeffs",IndexType::TheCellType(),NUM_SPECIES+2,
                 var_names_transp_coeff,pelelm_dertransportcoeff,the_same_box);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("cc_transport_coeffs",desc_lst,State_Type,
                          first_spec,NUM_SPECIES);




  if (NUM_SPECIES > 0)
  {
    //
    // rho-sum rhoY.
    //
    derive_lst.add("rhominsumrhoY",IndexType::TheCellType(),1,pelelm_drhomry,the_same_box);
    derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,Density,1);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
      const int comp = first_spec + i;
      derive_lst.addComponent("rhominsumrhoY",desc_lst,State_Type,comp,1);
    }
  }
  //
  // Sum rhoYdot
  //
  derive_lst.add("sumRhoYdot",IndexType::TheCellType(),1,pelelm_dsrhoydot,the_same_box);
  for (int i = 0; i < NUM_SPECIES; i++)
  {
    derive_lst.addComponent("sumRhoYdot",desc_lst,RhoYdot_Type,i,1);
  }
  //
  // **************  DEFINE DERIVED QUANTITIES ********************
  //
  //
  // average pressure
  //
  derive_lst.add("avg_pressure",IndexType::TheCellType(),1,pelelm_deravgpres,the_nodes);
  derive_lst.addComponent("avg_pressure",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in X direction.
  //
  derive_lst.add("gradpx",IndexType::TheCellType(),1,pelelm_dergrdpx,the_nodes);
  derive_lst.addComponent("gradpx",desc_lst,Press_Type,Pressure,1);
  //
  // Pressure gradient in Y direction.
  //
  derive_lst.add("gradpy",IndexType::TheCellType(),1,pelelm_dergrdpy,the_nodes);
  derive_lst.addComponent("gradpy",desc_lst,Press_Type,Pressure,1);

#if (AMREX_SPACEDIM == 3)
  //
  // Pressure gradient in Z direction.
  //
  derive_lst.add("gradpz",IndexType::TheCellType(),1,pelelm_dergrdpz,the_nodes);
  derive_lst.addComponent("gradpz",desc_lst,Press_Type,Pressure,1);
#endif
  //
  // Magnitude of vorticity.
  //
  derive_lst.add("mag_vort",IndexType::TheCellType(),1,pelelm_mgvort,grow_box_by_two);
  derive_lst.addComponent("mag_vort",desc_lst,State_Type,Xvel,AMREX_SPACEDIM);

#ifdef DO_LMC_FORCE
  //
  // forcing - used to calculate the rate of injection of energy
  //
  derive_lst.add("forcing",IndexType::TheCellType(),1,DeriveFunc3D(FORT_DERFORCING),the_same_box);
  derive_lst.addComponent("forcing",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("forcing",desc_lst,State_Type,Xvel,AMREX_SPACEDIM);
  //
  // forcex - used to put the forcing term in the plot file
  //
  derive_lst.add("forcex",IndexType::TheCellType(),1,DeriveFunc3D(FORT_DERFORCEX),the_same_box);
  derive_lst.addComponent("forcex",desc_lst,State_Type,Density,1);
  //    derive_lst.addComponent("forcex",desc_lst,State_Type,Xvel,AMREX_SPACEDIM);
  //
  // forcey - used to put the forcing term in the plot file
  //
  derive_lst.add("forcey",IndexType::TheCellType(),1,DeriveFunc3D(FORT_DERFORCEY),the_same_box);
  derive_lst.addComponent("forcey",desc_lst,State_Type,Density,1);
  //    derive_lst.addComponent("forcey",desc_lst,State_Type,Xvel,BL_SPACEDIM);
  //
  // forcez - used to put the forcing term in the plot file
  //
  derive_lst.add("forcez",IndexType::TheCellType(),1,DeriveFunc3D(FORT_DERFORCEZ),the_same_box);
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

  derive_lst.add("mixfrac_only",IndexType::TheCellType(),1,pelelm_dermixfrac,the_same_box);
  derive_lst.addComponent("mixfrac_only",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mixfrac_only",desc_lst,State_Type,first_spec,NUM_SPECIES);

  Vector<std::string> mix_and_diss(2);
  mix_and_diss[0] = "mixture_fraction";
  mix_and_diss[1] = "Scalar_diss";
  derive_lst.add("mixfrac",IndexType::TheCellType(),2,mix_and_diss,pelelm_dermixanddiss,grow_box_by_one,&cc_interp);
  derive_lst.addComponent("mixfrac",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("mixfrac",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("mixfrac",desc_lst,State_Type,first_spec,NUM_SPECIES);

  derive_lst.add("HeatRelease",IndexType::TheCellType(),1,pelelm_dhrr,the_same_box);
  derive_lst.addComponent("HeatRelease",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("HeatRelease",desc_lst,RhoYdot_Type,0,NUM_SPECIES);

  derive_lst.add("CMA",IndexType::TheCellType(),4,pelelm_dcma,the_same_box);
  derive_lst.addComponent("CMA",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("CMA",desc_lst,State_Type,first_spec,NUM_SPECIES);
  derive_lst.addComponent("CMA",desc_lst,State_Type,Temp,1);
  derive_lst.addComponent("CMA",desc_lst,RhoYdot_Type,0,NUM_SPECIES);


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
    RealBox realbox;
    if (ppr.countval("in_box_lo")) {
      std::vector<Real> box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
      ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
      ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
      realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
    }

    AMRErrorTagInfo info;

    if (realbox.ok()) {
      info.SetRealBox(realbox);
    }
    if (ppr.countval("start_time") > 0) {
      Real min_time; ppr.get("start_time",min_time);
      info.SetMinTime(min_time);
    }
    if (ppr.countval("end_time") > 0) {
      Real max_time; ppr.get("end_time",max_time);
      info.SetMaxTime(max_time);
    }
    if (ppr.countval("max_level") > 0) {
      int max_level; ppr.get("max_level",max_level);
      info.SetMaxLevel(max_level);
    }
    
    if (ppr.countval("value_greater")) {
      Real value; ppr.get("value_greater",value);
      std::string field; ppr.get("field_name",field);
      errtags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
    }
    else if (ppr.countval("value_less")) {
      Real value; ppr.get("value_less",value);
      std::string field; ppr.get("field_name",field);
      errtags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
    }
    else if (ppr.countval("vorticity_greater")) {
      Real value; ppr.get("vorticity_greater",value);
      const std::string field="mag_vort";
      errtags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
    }
    else if (ppr.countval("adjacent_difference_greater")) {
      Real value; ppr.get("adjacent_difference_greater",value);
      std::string field; ppr.get("field_name",field);
      errtags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
    }
    else if (realbox.ok())
    {
      errtags.push_back(AMRErrorTag(info));
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
  int nrhoydot = NUM_SPECIES;

  amrex::Print() << "RhoYdot_Type, nrhoydot = " << RhoYdot_Type << ' ' << nrhoydot << '\n';

  desc_lst.addDescriptor(RhoYdot_Type,IndexType::TheCellType(),
                         StateDescriptor::Point,ngrow,nrhoydot,
                         &cc_interp);


  amrex::StateDescriptor::BndryFunc pelelm_bndryfunc(pelelm_dummy_fill);
  pelelm_bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

	
  //const StateDescriptor& d_cell = desc_lst[State_Type];

  BCRec bc;	
  set_rhoydot_bc(bc,phys_bc);
  for (int i = 0; i < nrhoydot; i++)
  {
    const std::string name = "I_R[rhoY("+spec_names[i]+")]";
    desc_lst.setComponent(RhoYdot_Type, i, name.c_str(), bc,
                          pelelm_bndryfunc, &cc_interp, 0, nrhoydot-1);
  }
}
