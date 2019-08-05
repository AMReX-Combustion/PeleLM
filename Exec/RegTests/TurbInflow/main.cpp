#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>

using namespace amrex;

extern "C" {
  void fillinit(const int* turbin_enc, const int* turbin_len);

  void fillvel(const int* lo, const int* hi,
               Real* dat, const int* dlo, const int* dhi,
               const Real* dx, const Real* plo);
}

Vector<int>
encodeStringForFortran(const std::string& astr)
{
  long length = astr.size();
  Vector<int> result(length);
  for (int i = 0; i < length; ++i)
    result[i] = astr[i];
  return result;
}

static
void
Extend (FArrayBox& xfab,
        FArrayBox& vfab,
        const Box& domain)
{
  Box tbx = vfab.box();

  tbx.setBig(0, domain.bigEnd(0) + 3);

  const int ygrow = BL_SPACEDIM==3 ? 3 : 1;

  tbx.setBig(1, domain.bigEnd(1) + ygrow);

  xfab.resize(tbx,1);

  Box orig_box = vfab.box();
  vfab.shift(0, 1);
  vfab.shift(1, 1);
  xfab.copy(vfab); // (0,0)

  vfab.shift(0, domain.length(0)); // (1,0)
  xfab.copy(vfab);
  vfab.shift(1, domain.length(1)); // (1,1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0)); // (0,1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0)); // (-1,1)
  xfab.copy(vfab);
  vfab.shift(1, -domain.length(1)); // (-1,0)
  xfab.copy(vfab);
  vfab.shift(1, -domain.length(1)); // (-1,-1)
  xfab.copy(vfab);
  vfab.shift(0, domain.length(0)); // (0,-1)
  xfab.copy(vfab);
  vfab.shift(0, domain.length(0)); // (1,-1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0) - 1);
  vfab.shift(1,  domain.length(1) - 1);

  if (vfab.box() != orig_box) Abort("Oops, something bad happened");
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string TurbDir("Turb");

    if (ParallelDescriptor::IOProcessor())
      if (!UtilCreateDirectory(TurbDir, 0755))
        CreateDirectoryFailed(TurbDir);

    std::string Hdr = TurbDir; Hdr += "/HDR";
    std::string Dat = TurbDir; Dat += "/DAT";

    std::ofstream ifsd, ifsh;

    ifsh.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsh.good())
      FileOpenFailed(Hdr);

    ifsd.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsd.good())
      FileOpenFailed(Dat);

    Box box_turb(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(63,63,63)));
    RealBox rb_turb({D_DECL(0,0,0)},
                    {D_DECL(1,1,1)});
    int coord_turb(0);
    Array<int,BL_SPACEDIM> per_turb = {D_DECL(1,1,1)};
    Geometry geom_turb(box_turb,rb_turb,coord_turb,per_turb);
    const Real* dx_turb = geom_turb.CellSize();

    //
    // Write the first part of the Turb header.
    // Note that this is solely for periodic style inflow files.
    //
    Box box_turb_io(box_turb);
    box_turb_io.setBig(0, box_turb.bigEnd(0) + 3);
    box_turb_io.setBig(1, box_turb.bigEnd(1) + 3);
    box_turb_io.setBig(2, box_turb.bigEnd(2) + 1);

    ifsh << box_turb_io.length(0) << ' '
         << box_turb_io.length(1) << ' '
         << box_turb_io.length(2) << '\n';

    ifsh << rb_turb.length(0) + 2*dx_turb[0] << ' '
         << rb_turb.length(1) + 2*dx_turb[1] << ' '
         << rb_turb.length(2)                << '\n';

    ifsh << per_turb[0] << ' ' << per_turb[1] << ' ' << per_turb[2] << '\n';

    // Create a field to shove in
    FArrayBox vel_turb(box_turb,BL_SPACEDIM);
    Array4<Real> const& fab = vel_turb.array();
    AMREX_PARALLEL_FOR_3D ( box_turb, i, j, k,
    {
      Real x = (i+0.5)*dx_turb[0] + rb_turb.lo()[0];
      Real y = (j+0.5)*dx_turb[1] + rb_turb.lo()[1];
      Real z = (k+0.5)*dx_turb[2] + rb_turb.lo()[2];
      Real blobr = 0.25 * (rb_turb.hi()[0] - rb_turb.lo()[0]);
      Real blobx = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
      Real bloby = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
      Real blobz = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
      Real r = std::sqrt((x-blobx)*(x-blobx) +
                         (y-bloby)*(y-bloby) +
                         (z-blobz)*(z-blobz));
      fab(i,j,k,0) = r <= blobr ? 1 : 0;
      fab(i,j,k,1) = 2.;
      fab(i,j,k,2) = 3.;
    });

#if 0
    BoxArray ba(box_turb);
    DistributionMapping dm(ba);
    MultiFab vel(ba,dm,BL_SPACEDIM,0);
    vel[0].copy(vel_turb);
    VisMF::Write(vel,"TEST");
#endif

    // Dump field as a "turbulence file"
    IntVect sm = box_turb.smallEnd();
    IntVect bg = box_turb.bigEnd();
    int dir = BL_SPACEDIM - 1;
    FArrayBox xfab,TMP;
    //
    // We work on one cell wide Z-planes.
    // We first do the lo BL_SPACEDIM plane.
    // And then all the other planes in xhi -> xlo order.
    //
    for (int d = 0; d < BL_SPACEDIM; ++d)
    {
      bg[dir] = sm[dir];
      {
        Box bx(sm,bg);
        TMP.resize(bx,1);
        TMP.copy(vel_turb,bx,d,bx,0,1);
        Extend(xfab, TMP, box_turb);
        ifsh << ifsd.tellp() << std::endl;
        xfab.writeOn(ifsd);
      }
      for (int i = box_turb.bigEnd(dir); i >= box_turb.smallEnd(dir); i--)
      {
        sm[dir] = i;
        bg[dir] = i;
        Box bx(sm,bg);
        TMP.resize(bx,1);
        TMP.copy(vel_turb,bx,d,bx,0,1);
        Extend(xfab, TMP, box_turb);
        ifsh << ifsd.tellp() << std::endl;
        xfab.writeOn(ifsd);
      }
    }

    // Now that turbulence file written, read from it to fill the result
    auto TurbDirEnc = encodeStringForFortran(TurbDir);
    int len = TurbDirEnc.size();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      fillinit(&(TurbDirEnc[0]),&len);
    }
    
    Box box_res(IntVect(D_DECL(0,0,0)),
                IntVect(D_DECL(255,255,127)));
    RealBox rb_res({D_DECL(-1,-1,0)},
                   {D_DECL(1,1,1)});
    int coord_res(0);
    Array<int,BL_SPACEDIM> per_res = {D_DECL(1,1,1)};
    Geometry geom_res(box_res,rb_res,coord_res,per_res);
    Box domain = geom_res.Domain();
    BoxArray ba_res(domain);
    ba_res.maxSize(16);
    DistributionMapping dmap_res(ba_res);
    MultiFab mf_res(ba_res,dmap_res,BL_SPACEDIM,0);

    mf_res.setVal(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf_res,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& box = mfi.tilebox();
      fillvel(BL_TO_FORTRAN_BOX(box),
              BL_TO_FORTRAN_ANYD(mf_res[mfi]),
              geom_res.CellSize(), geom_res.ProbLo());
    }

    VisMF::Write(mf_res,"JUNK");
  }
  Finalize();
}
