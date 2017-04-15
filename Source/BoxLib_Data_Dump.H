#ifndef _BoxLib_Data_Dump_H_
#define _BoxLib_Data_Dump_H_

#include <iostream>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

extern "C" {
    
  // dump multifab to file, general
  bool writeMFgen (const MultiFab* mf,
                   const char*     file,
                   int             comp,
                   bool            abort)
  {
    int sc = (comp<0 ? 0 : comp);
    int nc = (comp<0 ? mf->nComp() : 1);
    MultiFab tmf(mf->boxArray(), nc, 0);
    MultiFab::Copy(tmf,*mf,sc,0,nc,0);
        
    std::string FullPath = file;
    if (FullPath!="" && FullPath[FullPath.length()-1] != '/')
      FullPath += '/';
        
    if (ParallelDescriptor::IOProcessor())
      if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
        BoxLib::CreateDirectoryFailed(FullPath);
        
    static const std::string MultiFabBaseName("MultiFab");
    FullPath += MultiFabBaseName;
	
    VisMF::Write(tmf,FullPath);
        
    if (abort)
      BoxLib::Abort("DEBUG ABORT::writeMF");
        
    return true;
  }
    
  // dump fab to file, general
  bool writeFABgen (const FArrayBox* fab,
                    const char*      file,
                    int              comp,
                    bool             abort)
  {
    int sc = (comp<0 ? 0 : comp);
    int nc = (comp<0 ? fab->nComp() : 1);
    FArrayBox tfab(fab->box(), nc);
    tfab.copy(*fab,sc,0,nc);
        
    std::ofstream ofs;
    ofs.open(file);
    if (!ofs.good())
      std::cerr << "Cannot open file: " << file << std::endl;
        
    tfab.writeOn(ofs);
        
    if (abort)
      BoxLib::Abort("DEBUG ABORT::writeFab");
        
    return true;
  }

  // dump multifab component to file, no abort
  bool writeMFc (const MultiFab* mf,
                 const char*     file,
                 int             comp)
  {
    return writeMFgen(mf,file,comp,false);
  }
    
  // dump multifab all components to file, no abort
  bool writeMF (const MultiFab* mf,
                const char*     file)
  {
    return writeMFgen(mf,file,-1,false);
  }
    
  // dump multifab one components to file, abort
  bool writeMFca (const MultiFab* mf,
                  const char*     file,
                  int             comp)
  {
    return writeMFgen(mf,file,comp,true);
  }
    
  // dump multifab all components to file, abort
  bool writeMFa (const MultiFab* mf,
                 const char*     file)
  {
    return writeMFgen(mf,file,-1,true);
  }
    
  // dump fab component to file, no abort
  bool writeFABc (const FArrayBox* fab,
                  const char*      file,
                  int              comp=-1)
  {
    return writeFABgen(fab,file,comp,false);
  }
    
  // dump fab all components to file, no abort
  bool writeFAB (const FArrayBox* fab,
                 const char*      file)
  {
    return writeFABgen(fab,file,-1,false);
  }
    
  // dump fab one components to file, abort
  bool writeFABca (const FArrayBox* fab,
                   const char*      file,
                   int              comp)
  {
    return writeFABgen(fab,file,comp,true);
  }
    
  // dump fab all components to file, abort
  bool writeFABa (const FArrayBox* fab,
                  const char*      file)
  {
    return writeFABgen(fab,file,-1,true);
  }
    
  // call amrvis to show a multifab on comp
  bool showMFgen(const MultiFab* mf,
                 int             comp)
  {
    std::string file = "DEBUGMF";
    writeMFgen(mf,file.c_str(),comp,false);
#if BL_SPACEDIM==2
    std::string cmd = "amrvis2d -mf " + file + "/MultiFab &";
#else
    std::string cmd = "amrvis3d -mf " + file + "/MultiFab &";
#endif
    int ierr = system(cmd.c_str());
    return bool(ierr);
  }
    
  // call amrvis to show a multifab, all comps
  bool showMF (const MultiFab* mf)
  {
    return showMFgen(mf,-1);
  }
    
  // call amrvis to show a fab
  bool showFABgen(const FArrayBox* fab,
                  int              comp)
  {
    std::string file = "DEBUGFAB";
    writeFABgen(fab,file.c_str(),comp,false);
#if BL_SPACEDIM==2
    std::string cmd = "amrvis2d -fab " + file + " &";
#else
    std::string cmd = "amrvis3d -fab " + file + " &";
#endif
    int ierr = system(cmd.c_str());
    return bool(ierr);
  }

  // call amrvis to show a fab
  bool showFAB(const FArrayBox* fab)
  {
    return showFABgen(fab,-1);
  }
    
}

void NeverCalledEXT()
{
  writeMFgen(0,"",-1,false);
  writeFABgen(0,"",-1,false);
  writeMFc(0,"",-1);
  writeMF(0,"");
  writeMFca(0,"",-1);
  writeMFa(0,"");
  writeFABc(0,"",-1);
  writeFAB(0,"");
  writeFABca(0,"",-1);
  writeFABa(0,"");
  showMFgen(0,-1);
  showFABgen(0,-1);
  showMF(0);
  showFAB(0);
}





#endif
