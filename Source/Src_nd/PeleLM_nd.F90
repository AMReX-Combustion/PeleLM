#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif 

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "PPHYS_CONSTANTS.H"
#include "mechanism.h"

module PeleLM_nd

  use amrex_fort_module,  only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_filcc_module, only : amrex_filccn

  implicit none

  private

  public :: part_cnt_err, &
            valgt_error, vallt_error, magvort_error, diffgt_error, &
            FORT_AVERAGE_EDGE_STATES

contains



! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on whether or not
! ::: they contain any particles.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

   subroutine part_cnt_err(tag, t_lo, t_hi,&
                           set, clear,&
                           var, v_lo, v_hi,&
                           lo, hi, nd, domlo, domhi,&
                           delta, xlo, problo, time, level) &
                           bind(C, name="part_cnt_err")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nd) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( var(i,j,k,1) .gt. zero) tag(i,j,k) = set
               !if (var(i,j,k,1) .gt. zero) print *,'TAGGING ',i,j,k
            end do
         end do
      end do

   end subroutine part_cnt_err

!=========================================================
!  Error tagging function : greater than 
!=========================================================

   subroutine valgt_error ( tag, t_lo, t_hi,&
                            set, clear,&
                            var, v_lo, v_hi,&
                            lo, hi, nd, domlo, domhi,&
                            delta, xlo, problo, time, level,value) &
                            bind(C, name="valgt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), var(i,j,k,1)>value)
            end do
         end do
      end do

   end subroutine valgt_error

!=========================================================
!  Error tagging function : lower than 
!=========================================================

   subroutine vallt_error ( tag, t_lo, t_hi,&
                            set, clear,&
                            var, v_lo, v_hi,&
                            lo, hi, nd, domlo, domhi,&
                            delta, xlo, problo, time, level,value) &
                            bind(C, name="vallt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), var(i,j,k,1)<value)
            end do
         end do
      end do

   end subroutine vallt_error

!=========================================================
!  Error tagging function : vorticity mag 
!  Value depends on level
!=========================================================

   subroutine magvort_error ( tag, t_lo, t_hi,&
                              set, clear,&
                              var, v_lo, v_hi,&
                              lo, hi, nd, domlo, domhi,&
                              delta, xlo, problo, time, level,value) &
                              bind(C, name="magvort_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      integer :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j,k) = merge(set, tag(i,j,k), ABS(var(i,j,k,1)) >= value*2.d0**level)
            end do
         end do
      end do

   end subroutine magvort_error

!=========================================================
!  Error tagging function : diff adj. cells greater than 
!=========================================================

   subroutine diffgt_error ( tag, t_lo, t_hi,&
                             set, clear,&
                             var, v_lo, v_hi,&
                             lo, hi, nd, domlo, domhi,&
                             delta, xlo, problo, time, level,value) &
                             bind(C, name="diffgt_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      integer :: v_lo(3), v_hi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, nd, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T,  dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),1) :: var
      REAL_T :: delta(3), xlo(3), problo(3), time, value

      REAL_T  :: axp, axm, ayp, aym, azp, azm, aerr
      integer :: i, j, k

      axp = 0.0d0
      axm = 0.0d0
      ayp = 0.0d0
      aym = 0.0d0
      azp = 0.0d0
      azm = 0.0d0

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               axp = ABS(var(i+1,j,k,1) - var(i,j,k,1))
               axm = ABS(var(i-1,j,k,1) - var(i,j,k,1))
               ayp = ABS(var(i,j+1,k,1) - var(i,j,k,1))
               aym = ABS(var(i,j-1,k,1) - var(i,j,k,1))
#if ( AMREX_SPACEDIM == 3 )
               azp = ABS(var(i,j,k+1,1) - var(i,j,k,1))
               azm = ABS(var(i,j,k-1,1) - var(i,j,k,1))
#endif
               aerr = MAX(azp,MAX(azm,MAX(axp,MAX(axm,MAX(ayp,aym)))))

               if ( aerr>=value ) then
                  tag(i,j,k) = set
               endif
            end do
         end do
      end do

   end subroutine diffgt_error

!=========================================================
!  Error tagging function : tag box 
!=========================================================

   subroutine box_error ( tag, t_lo, t_hi,&
                          set, clear,&
                          boxlo, boxhi, lo, hi,&
                          domlo, domhi,&
                          delta, xlo, problo, time, level) &
                          bind(C, name="box_error")

      implicit none

      integer :: t_lo(3), t_hi(3)
      REAL_T :: boxlo(3), boxhi(3)
      integer :: lo(3), hi(3), domlo(3), domhi(3)
      integer :: set, clear, level
      integer, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: tag
      REAL_T :: delta(3), xlo(3), problo(3), time

      integer :: i, j, k
      REAL_T  :: x, y, z

      do k = lo(3), hi(3)
#if ( AMREX_SPACEDIM == 3 )
         z = (float(k)+.5)*delta(3)+problo(3)
         if (z.ge.boxlo(3) .and. z.le.boxhi(3)) then
#endif
            do j = lo(2), hi(2)
               y = (float(j)+.5)*delta(2)+problo(2)
               if (y.ge.boxlo(2) .and. y.le.boxhi(2)) then
                  do i = lo(1), hi(1)
                     x = (float(i)+.5)*delta(1)+problo(1)

                     if (x.ge.boxlo(1) .and. x.le.boxhi(1) &
                          .and. y.ge.boxlo(2) .and. y.le.boxhi(2) &
#if ( AMREX_SPACEDIM == 3 )
                          .and. z.ge.boxlo(3) .and. z.le.boxhi(3) &
#endif
                          ) then
                        tag(i,j,k) = set
                     endif
                  end do
               endif
            end do
#if ( AMREX_SPACEDIM == 3 )
         endif
#endif
      end do

   end subroutine box_error

!=========================================================
!  This routine averages the mac face velocities for makeforce at half time
!=========================================================

   subroutine FORT_AVERAGE_EDGE_STATES( vel, v_lo, v_hi,&
                                        umacx, ux_lo, ux_hi,&
                                        umacy, uy_lo, uy_hi,&
#if ( AMREX_SPACEDIM == 3 )
                                        umacz, uz_lo, uz_hi,&
#endif
                                        getForceVerbose)&
                                        bind(C, name="FORT_AVERAGE_EDGE_STATES")

      implicit none

      integer :: v_lo(3), v_hi(3)
      integer :: ux_lo(3), ux_hi(3)
      integer :: uy_lo(3), uy_hi(3)
#if ( AMREX_SPACEDIM == 3 )
      integer :: uz_lo(3), uz_hi(3)
#endif
      integer :: getForceVerbose
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), dim) :: vel
      REAL_T, dimension(ux_lo(1):ux_hi(1),ux_lo(2):ux_hi(2),ux_lo(3):ux_hi(3)) :: umacx
      REAL_T, dimension(uy_lo(1):uy_hi(1),uy_lo(2):uy_hi(2),uy_lo(3):uy_hi(3)) :: umacy
#if ( AMREX_SPACEDIM == 3 )
      REAL_T, dimension(uz_lo(1):uz_hi(1),uz_lo(2):uz_hi(2),uz_lo(3):uz_hi(3)) :: umacz
#endif

      REAL_T  :: velmin(3)
      REAL_T  :: velmax(3)
      integer :: isioproc

      integer :: i, j, k, n

      do n = 1, dim
         velmin(n) = 1.d234
         velmax(n) = -1.d234
      enddo

      do k = v_lo(3), v_hi(3)
         do j = v_lo(2), v_hi(2)
            do i = v_lo(1), v_hi(1)
               vel(i,j,k,1) = half*(umacx(i,j,k)+umacx(i+1,j,k))
               vel(i,j,k,2) = half*(umacy(i,j,k)+umacy(i,j+1,k))
#if ( AMREX_SPACEDIM == 3 )
               vel(i,j,k,3) = half*(umacz(i,j,k)+umacz(i,j,k+1))
#endif
               do n = 1, dim
                  velmin(n) = min(velmin(n),vel(i,j,k,n))
                  velmax(n) = max(velmax(n),vel(i,j,k,n))
               enddo
            enddo
         enddo
      enddo

      if (getForceVerbose.gt.0) then
         call bl_pd_is_ioproc(isioproc)
         if (isioproc.eq.1) then
            do n = 1, dim
               write (6,*) "mac velmin (",n,") = ",velmin(n)
               write (6,*) "mac velmax (",n,") = ",velmax(n)
            enddo
         endif
      endif

   end subroutine FORT_AVERAGE_EDGE_STATES

end module PeleLM_nd
