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

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDir(ProbParm const* a_prob_parm, ACParm const* a_ac_parm,
                               pele::physics::PMF::PmfData::DataContainer const* a_pmf_data)
      : lprobparm(a_prob_parm), lacparm(a_ac_parm), lpmfdata(a_pmf_data) {}

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int dcomp,
    const int /*numcomp*/,
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

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif
         //
         // Fill s_ext with the EXT_DIR BC value
         // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
         //
         bcnormal(x, s_ext, idir, 1, time, geom, *lprobparm, *lacparm, lpmfdata);

         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }

    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif

         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);

         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }
    }


    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif

         bcnormal(x, s_ext, idir, +1, time, geom, *lprobparm, *lacparm, lpmfdata);

         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }

    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif

         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);

         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }

    }

#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif

         bcnormal(x, s_ext, idir, +1, time, geom, *lprobparm, *lacparm, lpmfdata);
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }

    } else if (
       (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
       (iv[idir] > domhi[idir])) {

#ifdef PELE_USE_TURBINFLOW
         // fill s_ext with turbulent velocity data so that the user can modify it.
         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             s_ext[Xvel+n] = dest(iv, dcomp + n);
           }
         }
#endif

         bcnormal(x, s_ext, idir, -1, time, geom, *lprobparm, *lacparm, lpmfdata);

         if (orig_comp == Xvel){
           for (int n = 0; n < AMREX_SPACEDIM; n++) {
             dest(iv, dcomp + n) = s_ext[Xvel+n];
           }
         }
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
         }

    }
#endif
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
        ProbParm const* lprobparm = PeleLM::prob_parm_d;
        ACParm const* lacparm = PeleLM::ac_parm_d;
        pele::physics::PMF::PmfData::DataContainer const* lpmfdata = PeleLM::pmf_data.getDeviceData();

#ifdef PELE_USE_TURBINFLOW
        if (PeleLM::prob_parm->do_turb && scomp < AMREX_SPACEDIM)  {
          ProbParm *probparmDD = PeleLM::prob_parm_d;
          ProbParm *probparmDH = PeleLM::prob_parm;

          // Copy problem parameter structs to host
          amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

          for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {

            auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(),dir) & bx);
            if (bcr[1].lo()[dir]==EXT_DIR && bndryBoxLO.ok())
            {
              // Create box with ghost cells and set them to zero
              amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
              int Grow = 4;     // Being conservative
              for(int n=0;n<AMREX_SPACEDIM;n++)
                growVect[n] = Grow;
              growVect[dir] = 0;
              amrex::Box modDom = geom.Domain();
              modDom.grow(growVect);
              auto bndryBoxLO_ghost = amrex::Box(amrex::adjCellLo(modDom,dir,Grow) & bx);
              data.setVal<amrex::RunOn::Host>(0.0,bndryBoxLO_ghost,Xvel,AMREX_SPACEDIM);

              add_turb(bndryBoxLO, data, 0, geom, time, dir, amrex::Orientation::low, probparmDH->tp);
              probparmDH->turb_ok[dir] = true;
            }

            auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(),dir) & bx);
            if (bcr[1].hi()[dir]==EXT_DIR && bndryBoxHI.ok())
            {
              //Create box with ghost cells and set them to zero
              amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
              int Grow = 4;
              for(int n=0;n<AMREX_SPACEDIM;n++)
                growVect[n] = Grow;
              growVect[dir] = 0;
              amrex::Box modDom = geom.Domain();
              modDom.grow(growVect);
              auto bndryBoxHI_ghost = amrex::Box(amrex::adjCellHi(modDom,dir,Grow) & bx);
              data.setVal<amrex::RunOn::Host>(0.0,bndryBoxHI_ghost,Xvel,AMREX_SPACEDIM);

              add_turb(bndryBoxHI, data, 0, geom, time, dir, amrex::Orientation::high, probparmDH->tp);
              probparmDH->turb_ok[dir+AMREX_SPACEDIM] = true;
            }
          }

          // Copy problem parameter structs back to device
          amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
        }
#endif

        GpuBndryFuncFab<PeleLMCCFillExtDir> gpu_bndry_func(PeleLMCCFillExtDir{lprobparm,lacparm,lpmfdata});

        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

#ifdef PELE_USE_TURBINFLOW
        if (PeleLM::prob_parm->do_turb && scomp < AMREX_SPACEDIM)  {

          ProbParm *probparmDD = PeleLM::prob_parm_d;
          ProbParm *probparmDH = PeleLM::prob_parm;

          // Copy problem parameter structs to host
          amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

          for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            if (probparmDH->turb_ok[dir]) {
              probparmDH->turb_ok[dir] = false;
            }
            if (probparmDH->turb_ok[dir+AMREX_SPACEDIM]) {
              probparmDH->turb_ok[dir+AMREX_SPACEDIM] = false;
            }
          }

          // Copy problem parameter structs back to device
          amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
        }
#endif
	  
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
