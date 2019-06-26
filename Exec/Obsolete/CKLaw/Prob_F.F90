
#include <MyProb_F.H>

module probspec_module

  implicit none
  
  private
  
  public :: set_prob_spec            

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

end module probspec_module
