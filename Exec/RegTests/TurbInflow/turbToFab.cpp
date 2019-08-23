#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_FArrayBox.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string pltfile("plt");
    
    std::string flctfile("Turb_s3d_1E4"); pp.query("flctfile",flctfile);

    std::string hdr = flctfile; hdr += "/HDR";

    std::ifstream ifs;

    ifs.open(hdr.c_str(), std::ios::in);

    if (!ifs.good())
      amrex::FileOpenFailed(hdr);

    int kmax[3];
    for (int i = 0; i < 3; i++)
      ifs >> kmax[i];

    Real d[3];
    ifs >> d[0] >> d[1] >> d[2];

    int per[3];
    ifs >> per[0] >> per[1] >> per[2];
    int turb_ncomp =3;

    Vector<size_t> offset(kmax[2]*turb_ncomp,0);

    for (int i = 0; i < offset.size(); i++)
      ifs >> offset[i];

    std::string dat = flctfile; dat += "/DAT";

    std::ifstream ifsd;

    ifsd.open(dat.c_str(), std::ios::in);

    if (!ifs.good())
      amrex::FileOpenFailed(dat);

    Print() << "size: " << kmax[0] << " " << kmax[1] << " " << kmax[2] << std::endl;

    Box bx(IntVect(D_DECL(0,0,0)),
           IntVect(D_DECL(kmax[0]-1,kmax[1]-1,kmax[2]-1)));
    FArrayBox data(bx,3); data.setVal(-10);
    
    for (int comp=0; comp<3; ++comp) {
      for (int k=0; k<kmax[2]; ++k) {
        int plane = k+1;
        int ncomp = comp+1;
        const long start = offset[(ncomp - 1)*kmax[2] + (plane - 1)];

        ifsd.seekg(start, std::ios::beg);

        if (!ifs.good())
          amrex::Abort("getplane(): seekg() failed");

        FArrayBox fab;
        fab.readFrom(ifsd);
        
        Box thisBox(fab.box());
        if (k==0) {
          thisBox.setSmall(2,kmax[2]-1);
          thisBox.setBig(2,kmax[2]-1);
        }
        else
        {
          thisBox.setSmall(2,kmax[2]-k-1);
          thisBox.setBig(2,kmax[2]-k-1);
        }
        
        data.copy(fab,fab.box(),0,thisBox,ncomp-1,1);
      }
    }
    std::ofstream ofs(flctfile+".fab");
    data.writeOn(ofs);
    ofs.close();
  }
  Finalize();
}
