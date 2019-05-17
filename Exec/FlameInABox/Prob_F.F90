
#include <MyProb_F.H>

module probspec_module

  use fuego_chemistry

  implicit none
  
  private
  
  public :: set_prob_spec, set_Y_from_Phi
            

contains

  subroutine set_prob_spec(bath, fuel, oxid, prod, numspec) &
                                bind(C, name="set_prob_spec")
 
      use network,  only: nspec

      implicit none

#include <probdata.H>

      integer bath, fuel, oxid, prod, numspec

      fuelID = fuel + 1
      oxidID = oxid + 1
      prodID = prod + 1
      if (bath .le. 0) then
         call bl_pd_abort('no N2 species present in mechanism')
      endif
      bathID = bath + 1

      if (numspec .ne. Nspec) then
         call bl_pd_abort('number of species not consistent')
      endif
      
  end subroutine set_prob_spec

!---------------------------------

  subroutine set_Y_from_Phi(phi,Yt)bind(C, name="set_Y_from_Phi")
  
      use PeleLM_F, only: pphys_get_spec_name2
      use network,  only: nspec
  
      implicit none

#include <probdata.H>
#include <cdwrk.H>
#include <conp.H>

      REAL_T a, phi
      REAL_T Xt(nspec), Yt(nspec)
      integer n
      character*(maxspnml) name
      do n=1,nspec
         Xt(n) = zero
      enddo
      
!     Set "a" for computing X from phi
!     hc + a.O2 -> b.CO2 + c.H2O
      
      print * , "In 'set_Y_from_Phi' before pphys_get_spec_name2 ?", nspec, fuelID
      call pphys_get_spec_name2(name,fuelID)
      print *, "In 'set_Y_from_Phi' after pphys_get_spec_name2 ?", name

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'H2') then
         print *, " I've found my fuel ", name
         a = .5d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'CH3OCH3') then
         a = 3.0d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      print *, "In 'set_Y_from_Phi' bef call to CKXTY", oxidID, fuelID, bathID
      Xt(oxidID) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * Xt(oxidID) / a
      Xt(bathID)    = 1.d0 - Xt(fuelID) - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)
      
  end subroutine set_Y_from_Phi

end module probspec_module
