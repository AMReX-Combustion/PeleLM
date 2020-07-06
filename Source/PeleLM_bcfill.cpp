#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <IndexDefines.H>

using namespace amrex;

struct PeleLMFillExtDirGPU
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)

//amrex::Print() << "\n  HELLO from PeleLMFillExtDir  \n ";
//amrex::Print() << "\n  dcomp = " << dcomp << " numpmp = " << numcomp << "\n" ;
amrex::Print() << "\n  bcomp = " << bcomp << " orig_comp = " << orig_comp << "\n" ;

    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* prob_hi = geom.ProbHi();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + (iv[0] + 0.5) * dx[0], prob_lo[1] + (iv[1] + 0.5) * dx[1],
      prob_lo[2] + (iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    int NVAR = 1;
    if (orig_comp == Xvel){
#if AMREX_SPACEDIM < 3
      NVAR = 2;
#else
      NVAR = 3;
#endif
    }
    else if (orig_comp == DEF_first_spec){ 
     NVAR = NUM_SPECIES;
    }
    


amrex::Print() << "\n  NVAR = " << NVAR << "\n";

    amrex::Real s_int[NVAR] = {0.0};
    amrex::Real s_ext[NVAR] = {0.0};

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
//      bcnormal(x, s_int, s_ext, idir, +1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
//      bcnormal(x, s_int, s_ext, idir, -1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }



/*

    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
//      bcnormal(x, s_int, s_ext, idir, -1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domlo[idir], n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
      (iv[idir] > domhi[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domhi[idir], n);
      }
//      bcnormal(x, s_int, s_ext, idir, -1, time, geom);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#endif

*/

        }
};

struct PeleLMFillEdges
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
        }
};



void PeleLM_PressFill_CPU (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    // do something for external Dirichlet (BCType::ext_dir) if there are
}



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

amrex::Print() << "\n \n HELLO from pelelm_cc_ext_fill \n \n ";

//    if (Gpu::inLauncihRegion()) {
        GpuBndryFuncFab<PeleLMFillExtDirGPU> gpu_bndry_func(PeleLMFillExtDirGPU{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    } else {
        // Without EXT_DIR (e.g., inflow), we can pass a nullptr
//        CpuBndryFuncFab cpu_bndry_func(nullptr);
//        cpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    }
}


void pelelm_fillEdges (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

amrex::Print() << "\n \n HELLO from pelelm_fillEdges \n \n ";

//    if (Gpu::inLauncihRegion()) {
        GpuBndryFuncFab<PeleLMFillEdges> gpu_bndry_func(PeleLMFillEdges{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    } else {
        // Without EXT_DIR (e.g., inflow), we can pass a nullptr
//        CpuBndryFuncFab cpu_bndry_func(nullptr);
//        cpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    }
}


void pelelm_press_fill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{

amrex::Print() << "\n \n HELLO from pelelm_press_fill \n \n ";




CpuBndryFuncFab  cpu_bndry_func(PeleLM_PressFill_CPU);
//  cpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);

//    if (Gpu::inLauncihRegion()) {
//        GpuBndryFuncFab<PeleLMFillExtDir> gpu_bndry_func(PeleLMFillExtDir{});
//        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    } else {
        // Without EXT_DIR (e.g., inflow), we can pass a nullptr
//        CpuBndryFuncFab cpu_bndry_func(nullptr);
//        cpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
//    }
}
