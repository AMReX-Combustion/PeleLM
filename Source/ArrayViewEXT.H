
#ifndef _ArrayViewEXT_H_
#define _ArrayViewEXT_H_

#include <iostream>
#include <cstdio>
#include <cfloat>

#include <AMReX_MultiFab.H>

#ifdef (BL_USE_FLOAT)
#  define FP_CLASS fp_classf
#else
#  define FP_CLASS fp_class
#endif

const int MAXSTRINGWIDTH = 14;
const std::string formatString("%8.4e");

bool BadFPClass(Real fpNumber)
{
  switch(FP_CLASS(fpNumber))
  {
  case FP_SNAN:           // Signaling NaN (Not-a-Number)
  case FP_QNAN:           // Quiet NaN (Not-a-Number)
  case FP_POS_INF:        // +infinity
  case FP_NEG_INF:        // -infinity
  case FP_POS_DENORM:     // positive denormalized
  case FP_NEG_DENORM:     // negative denormalized
    return true;
      
  case FP_POS_NORM:       // positive normalized
  case FP_NEG_NORM:       // negative normalized
  case FP_POS_ZERO:       // +0.0 (positive zero)
  case FP_NEG_ZERO:       // -0.0 (negative zero)
    return false;
      
  default:
    return true;  // ?
  } // end switch
}


std::string FPMax(const FArrayBox& fab, const Box& inBox, int nComp)
{
  Box box = inBox & fab.box();
  BL_ASSERT(box.ok());
  Real dataMax = fab(box.smallEnd(), nComp);
  std::string dataString;
  bool badVal = false;
  for (IntVect idx=box.smallEnd(); !badVal && idx<=box.bigEnd(); box.next(idx))
  {
    Real val = fab(idx,nComp);
    if (BadFPClass(val)) {
      dataString = "*** Bad Float Value ***";
      badVal = true;
    } else {
      dataMax = Max(val, dataMax);
    }
  }
  if (!badVal)
  {
    char tmp[MAXSTRINGWIDTH];
    if (snprintf(tmp, MAXSTRINGWIDTH, formatString.c_str(), dataMax) >= MAXSTRINGWIDTH)
      amrex::Abort("FPMax: buffer too small");
    dataString = tmp;
  }
  return "Max: " + dataString;
}

std::string FPMin(const FArrayBox& fab, const Box& inBox, int nComp)
{
  Box box = inBox & fab.box();
  BL_ASSERT(box.ok());
  Real dataMin = fab(box.smallEnd(), nComp);
  std::string dataString;
  bool badVal = false;
  for (IntVect idx=box.smallEnd(); !badVal && idx<=box.bigEnd(); box.next(idx))
  {
    Real val = fab(idx,nComp);
    if (BadFPClass(val)) {
      dataString = "*** Bad Float Value ***";
      badVal = true;
    } else {
      dataMin = Min(val, dataMin);
    }
  }
  if (!badVal)
  {
    char tmp[MAXSTRINGWIDTH];
    if (snprintf(tmp, MAXSTRINGWIDTH, formatString.c_str(), dataMin) >= MAXSTRINGWIDTH)
      amrex::Abort("FPMin: buffer too small");
    dataString = tmp;
  }
  return "Min: " + dataString;
}


extern "C" {
    
  // dump fab to tmp file
  bool WriteFab(const FArrayBox* fab, const char* file)
  {
    ofstream os(file,ios::out);
    if (!os)
    {
      cout << "writeFab: Cannot write to output file!" << endl;
      return false;
    }
    fab->writeOn(os);
    return true;
  };
    
  // dump multifab to tmp file
  bool WriteMultiFab(const MultiFab* mf, const char* file)
  {
    ofstream os(file,ios::out);
    if (!os)
    {
      cout << "writeMultiFab: Cannot write to output file!" << endl;
      return false;
    }
    mf->writeOn(os);
    return true;
  };

  // dump multifab stats to cout
  bool MultiFabStats(const MultiFab* mf, int n)
  {
    cout << "MultiFab stats:" << '\n';
    int fComp, eComp;
    fComp = eComp = n;
    if (n < 0)
    {
      fComp = 0;
      eComp = mf->nComp() - 1;
    }
    for (int k=fComp; k<=eComp; ++k)
    {
      cout << "  comp: " << k;
	    
      for (int i=0; i<mf->length(); i++)
      {
        const FArrayBox& fab = (*mf)[i];
        const Box& box = mf->boxArray()[i];
        cout << " grid " << i << ":  valid:: "
             << FPMin(fab, box, k) << ", " 
             << FPMax(fab, box, k) << ", total::"
             << FPMin(fab, fab.box(), k) << ", " 
             << FPMax(fab, fab.box(), k) << '\n';
      }
    }
    return true;
  }
}

void NeverCalledEXT()
{
  WriteFab( NULL, "DORK" );
  WriteMultiFab( NULL, "DORK" );
  MultiFabStats( NULL, -1 );
}


// More functions

extern "C"
{
  void PrintBoxArray (const BoxArray& ba);
  void PrintBoxDomain (const BoxDomain& bd);
  void PrintBndryRegister (const BndryRegister& br);
  void PrintTagBox (const TAGBOX& tb);
  void PrintTagBoxArray (const TagBoxArray& tba);
  void TagBoxCount (const TAGBOX& tb);
  void TagBoxArrayCount (const TagBoxArray& tba);
}

void PrintBoxArray (const BoxArray& ba) { cout << ba << endl; }

void PrintBoxDomain (const BoxDomain& bd) { cout << bd << endl; }

void
PrintBndryRegister (const BndryRegister& br)
{
  cout << "printing br" << endl;
  for (OrientationIter oitr; oitr; oitr++)
  {
    const FabSet& fs = br[oitr()];
    for (int i=0; i<fs.length(); i++)
    {
      cout << "(face, number) = ( "
           << oitr() << ", " << i << ")" << endl;
      cout << fs[i] << endl;
    }
  }
}
    
void
PrintTagBox (const TAGBOX& tb)
{
  const Box& bx = tb.box();

  long count = 0;

  cout << "TagBox: box = " << bx << ":\n";

  for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
  {
    if (!(tb(p) == TAGBOX::CLEAR))
    {
      count++;
      cout << p << ' ';
    }
  }

  cout << "Total tagged cells = " << count << endl;
}

void
TagBoxCount (const TAGBOX& tb)
{
  const Box& bx = tb.box();

  long count = 0;

  for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
  {
    if (!(tb(p) == TAGBOX::CLEAR))
    {
      count++;
    }
  }

  cout << "Total tagged cells = " << count << endl;
}

void
PrintTagBoxArray (const TagBoxArray& tba)
{
  long count = 0;

  cout << "TagBoxArray:\n";

  for (int i = 0; i < tba.length(); i++)
  {
    const Box& bx = tba[i].box();

    cout << "\ti = " << i << ", box = " << bx << ":\n";

    for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
    {
      if (!(tba[i](p) == TAGBOX::CLEAR))
      {
        count++;
        cout << p << ' ';
      }
    }

    cout << '\n';
  }

  cout << "Total tagged cells = " << count << endl;
}

void
TagBoxArrayCount (const TagBoxArray& tba)
{
  long count = 0;

  for (int i = 0; i < tba.length(); i++)
  {
    const Box& bx = tba[i].box();

    for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p))
    {
      if (!(tba[i](p) == TAGBOX::CLEAR))
      {
        count++;
      }
    }
  }

  cout << "Total tagged cells = " << count << endl;
}

#endif
