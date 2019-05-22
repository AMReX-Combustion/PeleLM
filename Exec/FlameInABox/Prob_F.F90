
#include <MyProb_F.H>

module probspec_module

  implicit none
  
  private
  
  public :: set_prob_spec, set_Y_from_Phi
            

contains

  subroutine set_prob_spec(fuel, oxid, prod, numspec) &
                                bind(C, name="set_prob_spec")
 
      implicit none
#include <cdwrk.H>
#include <probdata.H>
      integer fuel, oxid, prod, numspec
      fuelID = fuel + 1
      oxidID = oxid + 1
      prodID = prod + 1

      if (numspec .ne. Nspec) then
         call bl_pd_abort('number of species not consistent')
      endif
      
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
      integer n
      character*(maxspnml) name
      do n=1,Nspec
         Xt(n) = zero
      enddo
      
!     Set "a" for computing X from phi
!     hc + a.O2 -> b.CO2 + c.H2O
      
      call get_spec_name(name,fuelID)

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'H2') then
         a = .5d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'CH3OCH3') then
         a = 3.0d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      Xt(oxidID) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * Xt(oxidID) / a
      Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)
      
  end subroutine set_Y_from_Phi

end module probspec_module
