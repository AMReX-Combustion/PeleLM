#include <AMReX_REAL.H>

module main

  use amrex_fort_module, only : dp_t=>amrex_real, dim=>amrex_spacedim
  implicit none
  private
  public :: fillvel, fillinit

contains
  subroutine fillinit(turbin_enc,turbin_len) bind(C, name="fillinit")
    use turbinflow_module, only : init_turbinflow
    implicit none
    integer, intent(in) :: turbin_len, turbin_enc(turbin_len)

    integer maxlen,i
    parameter (maxlen=256)
    character*(:), allocatable ::  turbin
    logical :: turb_is_cgs = .FALSE.
    
    if (turbin_len .gt. maxlen) then
       call bl_abort('turbin file name too long')
    end if

    allocate(character(len=turbin_len) :: turbin)
    do i = 1, turbin_len
       turbin(i:i) = char(turbin_enc(i))
    end do

    call init_turbinflow(turbin, turb_is_cgs)
    deallocate(turbin)

  end subroutine fillinit

  subroutine fillvel(lo, hi, dat, dlo, dhi, dx, plo) bind(C, name="fillvel")

    use turbinflow_module

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(dp_t), intent(inout) :: dat(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),dim)
    real(dp_t), intent(in) :: dx(dim), plo(dim)

    real(dp_t) :: xx(dlo(1):dhi(1)), yy(dlo(2):dhi(2)), zz, vfluc(dlo(1):dhi(1),dlo(2):dhi(2),3)
    integer :: i,j,k

    do i = dlo(1),dhi(1)
       xx(i) = (float(i)+.5)*dx(1)+plo(1)
    enddo
    do j = dlo(2),dhi(2)
       yy(j) = (float(j)+.5)*dx(2)+plo(2)
    enddo

    do k=dlo(3), dhi(3)
       zz = (float(k)+0.5)*dx(3) + plo(3)
       dat(dlo(1):dhi(1),dlo(2):dhi(2),k,1:3) = 0.d0
       vfluc = 0.d0
       call get_turbstate(dlo(1),dlo(2),dhi(1),dhi(2),xx,yy,zz,vfluc)
       do j = dlo(2),dhi(2)
          do i = dlo(1),dhi(1)
             if (xx(i).lt.-0.5 .or. xx(i).gt.0.5 .or. yy(j).lt.-0.5 .or. yy(j).gt.0.5) then
                vfluc(i,j,:) = 0.d0
             endif
          enddo
       enddo
       dat(dlo(1):dhi(1),dlo(2):dhi(2),k,1:3) = vfluc(dlo(1):dhi(1),dlo(2):dhi(2),1:3)
    enddo
    
  end subroutine fillvel
end module main
