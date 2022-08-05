#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <IndexDefines.H>
#include <pelelm_prob.H>
#include <PeleLM.H>


using namespace amrex;

struct PeleLMdummyFill
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int dcomp,
    const int /*numcomp*/,
    amrex::GeometryData const& geom,
    const amrex::Real /*time*/,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    
    const int* bc = bcr->data();

    // Shouldn't actually ever use this, just need something computable.
    // Set to some ridiculous value so we know if it does get used.
    amrex::Real s_ext[1] = {1.2345e40};

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

	 dest(iv, dcomp) = s_ext[0];
    
    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

	 dest(iv, dcomp) = s_ext[0];
    }


    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

      dest(iv, dcomp) = s_ext[0];

    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

	 dest(iv, dcomp) = s_ext[0];
    }

#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

	 dest(iv, dcomp) = s_ext[0];

    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

	  dest(iv, dcomp) = s_ext[0];
    }
#endif
  }
};

struct PeleLMNodalFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    // do something for external Dirichlet (BCType::ext_dir)
    amrex::Abort("PeleLMNodalFillExtDir: Need to write fill for external Dirichlet (BCType::ext_dir)");
  }
};

struct PeleLMFaceFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    // do something for external Dirichlet (BCType::ext_dir)
    amrex::Abort("PeleLMFaceFillExtDir: Need to write fill for external Dirichlet (BCType::ext_dir)");
  }
};

struct PeleLMCCFillExtDir
{
  ProbParm const* lprobparm; 
  ACParm const* lacparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  int m_do_turbinflow = 0;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDir(ProbParm const* a_prob_parm, ACParm const* a_ac_parm,
                               pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
                               int a_do_turbinflow)
      : lprobparm(a_prob_parm), lacparm(a_ac_parm), lpmfdata(a_pmf_data), m_do_turbinflow(a_do_turbinflow) {}

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int dcomp,
    const int numcomp,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int orig_comp) const
  {
    // do something for external Dirichlet (BCType::ext_dir)

    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
    prob_lo[0] + (iv[0] + 0.5) * dx[0], prob_lo[1] + (iv[1] + 0.5) * dx[1],
    prob_lo[2] + (iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    amrex::Real s_ext[DEF_NUM_STATE] = {0.0};

    bool iv_is_on_extdir = false;

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, 1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    }


    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, +1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    }

#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, +1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

         if (m_do_turbinflow && orig_comp == Xvel) {
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);
         iv_is_on_extdir = true;
    }
#endif

    if (iv_is_on_extdir) {
        if (orig_comp == Xvel){
          // IAMR sometimes call on a per velocity component basis
          // So don't loop over all the components if only one is requested
          for (int n = 0; n < std::min(AMREX_SPACEDIM,numcomp); n++) {
            dest(iv, dcomp + n) = s_ext[Xvel+n];
          }
        }
        else if (orig_comp == Yvel){
            dest(iv, dcomp) = s_ext[Yvel];
        }
#if AMREX_SPACEDIM == 3
        else if (orig_comp == Zvel){
            dest(iv, dcomp) = s_ext[Zvel];
        }
#endif
        else if (orig_comp == Density){
            dest(iv, dcomp) = s_ext[Density];
        }
        else if (orig_comp == DEF_first_spec){
          for (int n = 0; n < NUM_SPECIES; n++) {
            dest(iv, dcomp + n) = s_ext[DEF_first_spec+n];
          }
        }
        else if (orig_comp == DEF_RhoH){
           dest(iv, dcomp) = s_ext[DEF_RhoH];
        }
        else if (orig_comp == DEF_Temp){
           dest(iv, dcomp) = s_ext[DEF_Temp];
        }
        else if (orig_comp == DEF_RhoRT){
           dest(iv, dcomp) = 0.0;
        } else if (orig_comp == DEF_first_passive && DEF_NUM_PASSIVE > 0){
            for (int n = 0; n < DEF_NUM_PASSIVE; ++n) {
              dest(iv, dcomp + n) = s_ext[DEF_first_passive+n];
            }
        }
    }
 }
};


// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the desciptor set up in CNS::variableSetUp.

void pelelm_cc_ext_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
    // Get device pointers
    ProbParm const* lprobparm = PeleLM::prob_parm_d;
    ACParm const* lacparm = PeleLM::ac_parm_d;
    pele::physics::PMF::PmfData::DataContainer const* lpmfdata = PeleLM::pmf_data.getDeviceData();

    // Fill turbulence data if requested
    if (PeleLM::turb_inflow.is_initialized() && scomp < AMREX_SPACEDIM) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(),dir) & bx);
            if (bcr[1].lo()[dir]==EXT_DIR && bndryBoxLO.ok())
            {
                // Create box with ghost cells and set them to zero
                amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
                int Grow = 4;     // Being conservative
                for(int n=0;n<AMREX_SPACEDIM;n++) {
                    growVect[n] = Grow;
                }
                growVect[dir] = 0;
                amrex::Box modDom = geom.Domain();
                modDom.grow(growVect);
                auto bndryBoxLO_ghost = amrex::Box(amrex::adjCellLo(modDom,dir,Grow) & bx);
                data.setVal<amrex::RunOn::Host>(0.0,bndryBoxLO_ghost,Xvel,AMREX_SPACEDIM);
                PeleLM::turb_inflow.add_turb(bndryBoxLO, data, 0, geom, time, dir, amrex::Orientation::low);
            }

            auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(),dir) & bx);
            if (bcr[1].hi()[dir]==EXT_DIR && bndryBoxHI.ok())
            {
                //Create box with ghost cells and set them to zero
                amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
                int Grow = 4;
                for(int n=0;n<AMREX_SPACEDIM;n++) {
                    growVect[n] = Grow;
                }
                growVect[dir] = 0;
                amrex::Box modDom = geom.Domain();
                modDom.grow(growVect);
                auto bndryBoxHI_ghost = amrex::Box(amrex::adjCellHi(modDom,dir,Grow) & bx);
                data.setVal<amrex::RunOn::Host>(0.0,bndryBoxHI_ghost,Xvel,AMREX_SPACEDIM);
                PeleLM::turb_inflow.add_turb(bndryBoxHI, data, 0, geom, time, dir, amrex::Orientation::high);
            }
        }
    }

    GpuBndryFuncFab<PeleLMCCFillExtDir> gpu_bndry_func(PeleLMCCFillExtDir{lprobparm,lacparm,lpmfdata,PeleLM::turb_inflow.is_initialized()});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}


void pelelm_face_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

// The GpuBndryFuncFab routine is not yet supporting Edge Face data
// but it is called here so that it will automatically work when amrex support will be ok

        GpuBndryFuncFab<PeleLMFaceFillExtDir> gpu_bndry_func(PeleLMFaceFillExtDir{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}


void pelelm_press_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

        GpuBndryFuncFab<PeleLMNodalFillExtDir> gpu_bndry_func(PeleLMNodalFillExtDir{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}

void pelelm_dummy_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

        GpuBndryFuncFab<PeleLMdummyFill> gpu_bndry_func(PeleLMdummyFill{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

}
