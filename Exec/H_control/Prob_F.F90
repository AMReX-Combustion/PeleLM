
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

  subroutine set_Y_from_Phi(phi,Yt) bind(C,name="set_Y_from_Phi")

    use chem_driver, only: get_spec_name

    implicit none
#include <probdata.H>
#include <cdwrk.H>
#include <conp.H>
    REAL_T a, phi
    REAL_T alpha,beta,gamma,delt,factor
    REAL_T Xt(maxspec), Yt(maxspec)
    integer iO2,iH2,iCH4, iC12H26
    integer n, len
    character*(maxspnml) name

    len = len_trim(probtype)

    write(6,*)" should not be here"
    stop

    if (probtype(1:len).eq.BL_PROB_PREMIXED_FREE) then

       iO2 = -1; iH2 = -1; iCH4 = -1

       do n=1,Nspec
          call get_spec_name(name,n)
          if (name .eq. 'N2' ) iN2 = n
          if (name .eq. 'O2' ) iO2 = n
          if (name .eq. 'H2' ) iH2 = n
          if (name .eq. 'CH4' ) iCH4 = n
          if (name .eq. 'NC12H26' ) iC12H26 = n
       enddo

       do n = 1,Nspec
          Xt(n) = 0.d0
       end do

       alpha = H2_frac
       beta = 1.d0 - H2_frac
       gamma = (0.5d0*alpha + 2.d0*beta) / phi_in
       delt = gamma*.79d0/.21d0
       factor = alpha+beta+gamma+delt
       if (iH2  > 0) Xt(iH2) = alpha / factor
       if (iCH4 > 0) Xt(iCH4) = beta / factor
       if (iO2  > 0) Xt(iO2) = gamma / factor
       if (iN2  > 0) Xt(iN2) = delt / factor

    else

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
       else if (name .eq. 'NC12H26') then
          a = 18.50d0
       else
          call bl_abort('setupbc: Unknown fuel type')
       end if

       Xt(oxidID) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
       Xt(fuelID) = phi * Xt(oxidID) / a
       Xt(iN2)    = 1.d0 - Xt(fuelID) - Xt(oxidID)

    endif

    CALL CKXTY (Xt, Yt)
  end subroutine set_Y_from_Phi
end module probspec_module
