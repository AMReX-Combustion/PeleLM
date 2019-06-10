
#include <MyProb_F.H>

module probspec_module

  use fuego_chemistry

  implicit none
  
  private
  
  public :: set_prob_spec
            

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


end module probspec_module
