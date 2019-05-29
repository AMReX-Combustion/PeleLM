
#include <MyProb_F.H>

module probspec_module

  implicit none
  
  private
  
  public :: set_prob_spec, set_Y_from_Phi
            

contains

  subroutine set_prob_spec(fuel, oxid, prod, numspec) &
                                bind(C, name="set_prob_spec")
 
      implicit none
      integer fuel, oxid, prod, numspec
      
  end subroutine set_prob_spec

!---------------------------------

  subroutine set_Y_from_Phi(phi,Yt)bind(C, name="set_Y_from_Phi")
  
      use chem_driver, only: get_spec_name
  
      implicit none
#include <probdata.H>
#include <cdwrk.H>
#include <conp.H>
      REAL_T a, phi
      REAL_T Xt(maxspec), Yt(maxspec)
      
  end subroutine set_Y_from_Phi

end module probspec_module
