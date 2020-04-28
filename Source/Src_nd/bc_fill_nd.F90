#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "mechanism.h"

module bc_fill_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_filcc_module, only : amrex_filccn
  use mod_Fvar_def, only : domnlo
  use user_defined_fcts_nd_module, only : bcfunction

  implicit none
  
  private
  
  public :: den_fill, temp_fill, rhoh_fill, chem_fill, &
            xvel_fill, yvel_fill, zvel_fill, vel_fill, &
            adv_fill, press_fill

contains

! ::: -----------------------------------------------------------
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data and that all non-interior cells have
! ::         have been filled with a large real number.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: den      <=> density array
! ::: lo,hi     => index extent of den array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	             corner of den array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: -----------------------------------------------------------

   subroutine den_fill (den, d_lo, d_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="den_fill")
      
      implicit none

! In/Out      
      integer :: d_lo(3), d_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: den

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( d_lo, d_hi, den, d_lo, d_hi, 1, domlo, domhi, delta, xlo, bc)
                       

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(d_lo(1)<domlo(1))) then
         do i = d_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = d_lo(3), d_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = d_lo(2), d_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(d_hi(1)>domhi(1))) then
         do i = domhi(1)+1, d_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = d_lo(3), d_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = d_lo(2), d_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(d_lo(2)<domlo(2))) then
         do j = d_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = d_lo(3), d_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = d_lo(1), d_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(d_hi(2)>domhi(2))) then
         do j = domhi(2)+1, d_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = d_lo(3),d_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = d_lo(1), d_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(d_lo(3)<domlo(3))) then
         do k = d_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = d_lo(2), d_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = d_lo(1), d_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(d_hi(3)>domhi(3))) then
         do k = domhi(3)+1, d_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = d_lo(2), d_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = d_lo(1), d_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  den(i,j,k) = rho
               enddo
            enddo
         enddo
      endif
#endif   
#endif

  end subroutine den_fill

! ::: -----------------------------------------------------------
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data and that all non-interior cells have
! ::         have been filled with a large real number.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: temp     <=> temperature array
! ::: lo,hi     => index extent of temp array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	             corner of den array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: -----------------------------------------------------------

   subroutine temp_fill (temp, t_lo, t_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc)&
                         bind(C, name="temp_fill")
      
      implicit none

! In/Out      
      integer :: t_lo(3), t_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: temp

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( t_lo, t_hi, temp, t_lo, t_hi, 1, domlo, domhi, delta, xlo, bc)
                       

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(t_lo(1)<domlo(1))) then
         do i = t_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = t_lo(3), t_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = t_lo(2), t_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(t_hi(1)>domhi(1))) then
         do i = domhi(1)+1, t_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = t_lo(3), t_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = t_lo(2), t_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(t_lo(2)<domlo(2))) then
         do j = t_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = t_lo(3), t_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = t_lo(1), t_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(t_hi(2)>domhi(2))) then
         do j = domhi(2)+1, t_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = t_lo(3),t_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = t_lo(1), t_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(t_lo(3)<domlo(3))) then
         do k = t_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = t_lo(2), t_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = t_lo(1), t_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(t_hi(3)>domhi(3))) then
         do k = domhi(3)+1, t_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = t_lo(2), t_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = t_lo(1), t_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  temp(i,j,k) = T
               enddo
            enddo
         enddo
      endif
#endif   
#endif

  end subroutine temp_fill

! ::: -----------------------------------------------------------
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data and that all non-interior cells have
! ::         have been filled with a large real number.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: rhoh     <=> rhoh array
! ::: lo,hi     => index extent of rhoh array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	             corner of den array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: -----------------------------------------------------------

   subroutine rhoh_fill (rhoh, r_lo, r_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc)&
                         bind(C, name="rhoh_fill")
      
      implicit none

! In/Out      
      integer :: r_lo(3), r_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: rhoh

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( r_lo, r_hi, rhoh, r_lo, r_hi, 1, domlo, domhi, delta, xlo, bc)
                       

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(r_lo(1)<domlo(1))) then
         do i = r_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = r_lo(2), r_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(r_hi(1)>domhi(1))) then
         do i = domhi(1)+1, r_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = r_lo(2), r_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(r_lo(2)<domlo(2))) then
         do j = r_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(r_hi(2)>domhi(2))) then
         do j = domhi(2)+1, r_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = r_lo(3),r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(r_lo(3)<domlo(3))) then
         do k = r_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = r_lo(2), r_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(r_hi(3)>domhi(3))) then
         do k = domhi(3)+1, r_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = r_lo(2), r_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoh(i,j,k) = rho*h
               enddo
            enddo
         enddo
      endif
#endif   
#endif

  end subroutine rhoh_fill

! ::: -----------------------------------------------------------
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data and that all non-interior cells have
! ::         have been filled with a large real number.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: rhoY     <=> rho*Y array
! ::: lo,hi     => index extent of rhoY array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	             corner of den array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: stateID   => id index of state being filled
! ::: -----------------------------------------------------------

   subroutine chem_fill (rhoY, r_lo, r_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc, id)&
                         bind(C, name="chem_fill")
      
      implicit none

! In/Out      
      integer :: r_lo(3), r_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      integer :: id
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: rhoY

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( r_lo, r_hi, rhoY, r_lo, r_hi, 1, domlo, domhi, delta, xlo, bc)
                       

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(r_lo(1)<domlo(1))) then
         do i = r_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = r_lo(2), r_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(r_hi(1)>domhi(1))) then
         do i = domhi(1)+1, r_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = r_lo(2), r_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(r_lo(2)<domlo(2))) then
         do j = r_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = r_lo(3), r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(r_hi(2)>domhi(2))) then
         do j = domhi(2)+1, r_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = r_lo(3),r_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(r_lo(3)<domlo(3))) then
         do k = r_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = r_lo(2), r_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(r_hi(3)>domhi(3))) then
         do k = domhi(3)+1, r_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = r_lo(2), r_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = r_lo(1), r_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  rhoY(i,j,k) = rho*Yl(id)
               enddo
            enddo
         enddo
      endif
#endif   
#endif

  end subroutine chem_fill

! ::: ===========================================================
! Fill all chem species at once
! ::: ===========================================================

   subroutine all_chem_fill (rhoY, r_lo, r_hi, &
                             domlo, domhi, delta, &
                             xlo, time, bc)&
                             bind(C, name="all_chem_fill")
      
      implicit none

! In/Out
      integer :: r_lo(3), r_hi(3)
      integer :: bc(dim,2,NUM_SPECIES)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NUM_SPECIES) :: rhoY

! Local
      integer :: n

      do n = 1, NUM_SPECIES
         call chem_fill (rhoY(:,:,:,n), r_lo, r_hi, &
                         domlo, domhi, delta, xlo, time, bc(1,1,n), n-1)
      enddo

   end subroutine all_chem_fill

! ::: ===========================================================
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: vel       <= x velocity array
! ::: lo,hi     => index extent of vel array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of rho array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: ===========================================================

   subroutine xvel_fill (vel, v_lo, v_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc)&
                         bind(C, name="xvel_fill")
      
      implicit none

! In/Out      
      integer :: v_lo(3), v_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( v_lo, v_hi, vel, v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc)

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(v_lo(1)<domlo(1))) then
         do i = v_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(v_hi(1)>domhi(1))) then
         do i = domhi(1)+1, v_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(v_lo(2)<domlo(2))) then
         do j = v_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(v_hi(2)>domhi(2))) then
         do j = domhi(2)+1, v_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3),v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(v_lo(3)<domlo(3))) then
         do k = v_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(v_hi(3)>domhi(3))) then
         do k = domhi(3)+1, v_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(1)
               enddo
            enddo
         enddo
      endif
#endif   
#endif
      
   end subroutine xvel_fill

! ::: ===========================================================
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: vel       <= x velocity array
! ::: lo,hi     => index extent of vel array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of rho array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: ===========================================================

   subroutine yvel_fill (vel, v_lo, v_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc)&
                         bind(C, name="yvel_fill")
      
      implicit none

! In/Out      
      integer :: v_lo(3), v_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( v_lo, v_hi, vel, v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc)

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(v_lo(1)<domlo(1))) then
         do i = v_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(v_hi(1)>domhi(1))) then
         do i = domhi(1)+1, v_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(v_lo(2)<domlo(2))) then
         do j = v_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(v_hi(2)>domhi(2))) then
         do j = domhi(2)+1, v_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3),v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(v_lo(3)<domlo(3))) then
         do k = v_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(v_hi(3)>domhi(3))) then
         do k = domhi(3)+1, v_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(2)
               enddo
            enddo
         enddo
      endif
#endif   
#endif
      
   end subroutine yvel_fill

! ::: ===========================================================
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: vel       <= x velocity array
! ::: lo,hi     => index extent of vel array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of rho array
! ::: time      => problem evolution time
! ::: bc	       => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: ===========================================================

   subroutine zvel_fill (vel, v_lo, v_hi, &
                         domlo, domhi, delta, &
                         xlo, time, bc)&
                         bind(C, name="zvel_fill")
      
      implicit none

! In/Out      
      integer :: v_lo(3), v_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

! Local      
      REAL_T  :: x(3)
      REAL_T  :: vel_bc(3), rho, Yl(0:NUM_SPECIES-1), T, h

      integer :: i, j, k

      call amrex_filccn ( v_lo, v_hi, vel, v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc)

      ! X-low boundary
      if ((bc(1,1)==EXT_DIR).and.(v_lo(1)<domlo(1))) then
         do i = v_lo(1), domlo(1)-1
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif
      
      ! X-hi boundary
      if ((bc(1,2)==EXT_DIR).and.(v_hi(1)>domhi(1))) then
         do i = domhi(1)+1, v_hi(1)
            x(1) = (float(i)+.5)*delta(1)+domnlo(1)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do j = v_lo(2), v_hi(2)
                  x(2) = (float(j)+.5)*delta(2)+domnlo(2)
                  call bcfunction(x, delta, 1, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif    

#if ( AMREX_SPACEDIM >=2 )      
      ! Y-low boundary
      if ((bc(2,1)==EXT_DIR).and.(v_lo(2)<domlo(2))) then
         do j = v_lo(2), domlo(2)-1
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3), v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif    
      
      ! Y-hi boundary
      if ((bc(2,2)==EXT_DIR).and.(v_hi(2)>domhi(2))) then
         do j = domhi(2)+1, v_hi(2)
            x(2) = (float(j)+.5)*delta(2)+domnlo(2)
            do k = v_lo(3),v_hi(3)
               x(3) = (float(k)+.5)*delta(3)+domnlo(3)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 2, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif

#if ( AMREX_SPACEDIM == 3 )      
      ! Z-low boundary
      if ((bc(3,1)==EXT_DIR).and.(v_lo(3)<domlo(3))) then
         do k = v_lo(3), domlo(3)-1
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, 1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif    
      
      ! Z-hi boundary
      if ((bc(3,2)==EXT_DIR).and.(v_hi(3)>domhi(3))) then
         do k = domhi(3)+1, v_hi(3)
            x(3) = (float(k)+.5)*delta(3)+domnlo(3)
            do j = v_lo(2), v_hi(2)
               x(2) = (float(j)+.5)*delta(2)+domnlo(2)
               do i = v_lo(1), v_hi(1)
                  x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                  call bcfunction(x, delta, 3, -1, time, .true., &
                                  vel_bc, rho, Yl, T, h)
                  vel(i,j,k) = vel_bc(3)
               enddo
            enddo
         enddo
      endif
#endif   
#endif
      
   end subroutine zvel_fill

!===========================================================
! Fill all velocity components at once.
!===========================================================

   subroutine vel_fill (vel, v_lo, v_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="vel_fill")
      
      implicit none

! In/Out      
      integer :: v_lo(3), v_hi(3)
      integer :: bc(dim,2,3)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3) :: vel
      
      call xvel_fill (vel(:,:,:,1), v_lo, v_hi, & 
                      domlo, domhi, delta, xlo, time, bc(1,1,1))

#if ( AMREX_SPACEDIM >=2 )      
      call yvel_fill (vel(:,:,:,2), v_lo, v_hi, &
                      domlo, domhi, delta, xlo, time, bc(1,1,2))

#if ( AMREX_SPACEDIM ==3 )      
      call zvel_fill (vel(:,:,:,3), v_lo, v_hi, &
                      domlo, domhi, delta, xlo, time, bc(1,1,3))
#endif   
#endif

   end subroutine vel_fill

! ::: ===========================================================
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.  
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data and that all non-interior cells have
! ::         have been filled with a large real number.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: adv      <=  advected quantity array
! ::: lo/hi     => index extent of adv array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of adv array
! ::: time      => problem evolution time
! ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: ===========================================================

   subroutine adv_fill (adv, a_lo, a_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="adv_fill")
      
      implicit none

! In/Out      
      integer :: a_lo(3), a_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3)) :: adv

! Local      
      integer :: i, j, k

      call amrex_filccn ( a_lo, a_hi, adv, a_lo, a_hi, 1, domlo, domhi, delta, xlo, bc)

      ! Dim-low boundary
      if ((bc(dim,1)==EXT_DIR).and.(a_lo(dim)<domlo(dim))) then
#if ( AMREX_SPACEDIM == 2 )      
         do j = a_lo(2), domlo(2)-1
            do k = a_lo(3), a_hi(3)
#elif ( AMREX_SPACEDIM == 3 )      
         do k = a_lo(3), domlo(3)-1
            do j = a_lo(2), a_hi(2)
#endif
               do i = a_lo(1), a_hi(1)
                  adv(i,j,k) = 0.0d0
               enddo
            enddo
         enddo
      endif    

   end subroutine adv_fill

! ::: ===========================================================
! ::: This routine is called during a filpatch operation when
! ::: the patch to be filled falls outside the interior
! ::: of the problem domain.  You are requested to supply the
! ::: data outside the problem interior in such a way that the
! ::: data is consistant with the types of the boundary conditions
! ::: you specified in the C++ code.
! ::: 
! ::: NOTE:  you can assume all interior cells have been filled
! :::        with valid data.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: p        <=  pressure array
! ::: DIMS(p)   => index extent of p array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of rho array
! ::: time      => problem evolution time
! ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
! ::: ===========================================================

   subroutine press_fill (p, p_lo, p_hi, &
                          domlo, domhi, delta, &
                          xlo, time, bc)&
                          bind(C, name="press_fill")

      implicit none

! In/Out
      integer :: p_lo(3), p_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p

! Local
      integer :: i, j, k
      integer :: ilo, ihi, jlo, jhi, klo, khi
      logical :: fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi
      logical :: per_xlo, per_xhi, per_ylo, per_yhi, per_zlo, per_zhi

      fix_xlo = (p_lo(1) < domlo(1)) .and. (bc(1,1) /= INT_DIR)
      per_xlo = (p_lo(1) < domlo(1)) .and. (bc(1,1) == INT_DIR)
      fix_xhi = (p_hi(1) > domhi(1)) .and. (bc(1,2) /= INT_DIR)
      per_xhi = (p_hi(1) > domhi(1)) .and. (bc(1,2) == INT_DIR)
#if ( AMREX_SPACEDIM >= 2 )      
      fix_ylo = (p_lo(2) < domlo(2)) .and. (bc(2,1) /= INT_DIR)
      per_ylo = (p_lo(2) < domlo(2)) .and. (bc(2,1) == INT_DIR)
      fix_yhi = (p_hi(2) > domhi(2)) .and. (bc(2,2) /= INT_DIR)
      per_yhi = (p_hi(2) > domhi(2)) .and. (bc(2,2) == INT_DIR)
#if ( AMREX_SPACEDIM == 3 )      
      fix_zlo = (p_lo(3) < domlo(3)) .and. (bc(3,1) /= INT_DIR)
      per_zlo = (p_lo(3) < domlo(3)) .and. (bc(3,1) == INT_DIR)
      fix_zhi = (p_hi(3) > domhi(3)) .and. (bc(3,2) /= INT_DIR)
      per_zhi = (p_hi(3) > domhi(3)) .and. (bc(3,2) == INT_DIR)
#endif
#endif

      ilo = max(p_lo(1),domlo(1))
      ihi = min(p_hi(1),domhi(1))
      jlo = p_lo(2)
      klo = p_lo(3)
      jhi = p_hi(2)
      khi = p_hi(3)
#if ( AMREX_SPACEDIM >= 2 )      
      jlo = max(p_lo(2),domlo(2))
      jhi = min(p_hi(2),domhi(2))
#if ( AMREX_SPACEDIM == 3 )      
      klo = max(p_lo(3),domlo(3))
      khi = min(p_hi(3),domhi(3))
#endif
#endif

!***************
!  SETTING BL_XLO
!***************

      if (fix_xlo) then
         do i = p_lo(1), domlo(1)-1
            do k = p_lo(3), p_hi(3)
               do j = p_lo(2), p_hi(2)
                  p(i,j,k) = p(ilo,j,k)
               end do 
            end do
         end do

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_ylo) then
            do i = p_lo(1), domlo(1)-1
                 do j = p_lo(2), domlo(2)-1
                    do k = klo, khi
                       p(i,j,k) = p(ilo,jlo,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jlo,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jlo,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jlo,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jlo,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_yhi) then
            do i = p_lo(1), domlo(1)-1
                 do j = domhi(2)+1, p_hi(2)
                    do k = klo, khi
                       p(i,j,k) = p(ilo,jhi,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jhi,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jhi,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jhi,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jhi,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (fix_zlo) then
            do i = p_lo(1), domlo(1)-1
                 do j = jlo, jhi
                    do k = p_lo(3), domlo(3)-1
                       p(i,j,k) = p(ilo,j,klo)
                    end do
                 end do
            end do
              if (per_ylo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,j,klo)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,j,klo)
                       end do
                    end do
                 end do
              end if

         end if

         if (fix_zhi) then
            do i = p_lo(1), domlo(1)-1
                 do j = jlo, jhi
                    do k = domhi(3)+1, p_hi(3)
                       p(i,j,k) = p(ilo,j,khi)
                    end do
                 end do
            end do
              if (per_ylo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,j,khi)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,j,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif
 
#if ( AMREX_SPACEDIM >= 2 )      
         if (per_ylo) then
               do i = p_lo(1), domlo(1)-1
                  do k = klo,khi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do i = p_lo(1), domlo(1)-1
                  do k = klo,khi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
#endif
 
#if ( AMREX_SPACEDIM == 3 )      
         if (per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = jlo,jhi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = jlo,jhi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
#endif

      end if     ! End if on fix_xlo

!*****************************************************************************
! SETTING BL_XHI
!*****************************************************************************

      if (fix_xhi) then
         do i = domhi(1)+1, p_hi(1)
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ihi,j,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_ylo) then
            do i = domhi(1)+1, p_hi(1)
                 do j = p_lo(2), domlo(2)-1
                    do k = klo, khi
                       p(i,j,k) = p(ihi,jlo,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jlo,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jlo,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jlo,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jlo,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_yhi) then
            do i = domhi(1)+1, p_hi(1)
                 do j = domhi(2)+1, p_hi(2)
                    do k = klo, khi
                       p(i,j,k) = p(ihi,jhi,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jhi,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jhi,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jhi,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jhi,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (fix_zlo) then
            do i = domhi(1)+1, p_hi(1)
                 do j = jlo, jhi
                    do k = p_lo(3), domlo(3)-1
                       p(i,j,k) = p(ihi,j,klo)
                    end do
                 end do
            end do
              if (per_ylo) then
               do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,j,klo)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
               do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,j,klo)
                       end do
                    end do
                 end do
              end if
         end if

         if (fix_zhi) then
            do i = domhi(1)+1, p_hi(1)
                 do j = jlo, jhi
                    do k = domhi(3)+1, p_hi(3)
                       p(i,j,k) = p(ihi,j,khi)
                    end do
                 end do
            end do
              if (per_ylo) then
               do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,j,khi)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
               do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,j,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (per_ylo) then
             do i = domhi(1)+1, p_hi(1)
                  do k = klo,khi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
             do i = domhi(1)+1, p_hi(1)
                  do k = klo,khi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (per_zlo) then
             do i = domhi(1)+1, p_hi(1)
                  do j = jlo,jhi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
              do i = domhi(1)+1, p_hi(1)
                  do j = jlo,jhi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (per_ylo .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

      end if   ! End if on fix_xhi

!*****************************************************************************
! SETTING BL_YLO
!*****************************************************************************

#if ( AMREX_SPACEDIM >= 2 )
      if (fix_ylo) then
         do j = p_lo(2), domlo(2)-1
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jlo,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM == 3 )
         if (fix_zlo) then
            do j = p_lo(2), domlo(2)-1
                 do k = p_lo(3), domlo(3)-1
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jlo,klo)
                    end do
                 end do
            end do
            if (per_xlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
         end if

         if (fix_zhi) then
            do j = p_lo(2), domlo(2)-1
                 do k = domhi(3)+1, p_hi(3)
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jlo,khi)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jlo,khi)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jlo,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

         if (per_xlo) then
               do j = p_lo(2), domlo(2)-1
                  do k = klo,khi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = p_lo(2), domlo(2)-1
                  do k = klo,khi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

#if ( AMREX_SPACEDIM == 3 )
         if (per_zlo) then
               do j = p_lo(2), domlo(2)-1
                  do i = ilo,ihi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = p_lo(2), domlo(2)-1
                  do i = ilo,ihi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
#endif

      end if      ! End if on ylo
 
!*****************************************************************************
! SETTING BL_YHI
!*****************************************************************************

      if (fix_yhi) then
         do j = domhi(2)+1, p_hi(2)
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jhi,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM == 3 )
         if (fix_zlo) then
            do j = domhi(2)+1, p_hi(2)
                 do k = p_lo(3), domlo(3)-1
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jhi,klo)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(i,jhi,klo)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(i,jhi,klo)
                       end do
                    end do
                 end do
              end if
         end if

         if (fix_zhi) then
            do j = domhi(2)+1, p_hi(2)
                 do k = domhi(3)+1, p_hi(3)
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jhi,khi)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jhi,khi)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jhi,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

         if (per_xlo) then
               do j = domhi(2)+1, p_hi(2)
                  do k = klo,khi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = domhi(2)+1, p_hi(2)
                  do k = klo,khi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

#if ( AMREX_SPACEDIM == 3 )
         if (per_zlo) then
               do j = domhi(2)+1, p_hi(2)
                  do i = ilo,ihi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = domhi(2)+1, p_hi(2)
                  do i = ilo,ihi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
#endif

      end if     ! End if on yhi

!*****************************************************************************
! SETTING BL_ZLO
!*****************************************************************************

#if ( AMREX_SPACEDIM == 3 )
      if (fix_zlo) then
         do k = p_lo(3), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,klo)
               end do
            end do
         end do

         if (per_xlo) then
               do k = p_lo(3), domlo(3)-1
                  do j = jlo,jhi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = p_lo(3), domlo(3)-1
                  do j = jlo,jhi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = ilo,ihi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = ilo,ihi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = p_lo(1), domlo(1)-1
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = p_lo(1), domlo(1)-1
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = domhi(1)+1, p_hi(1)
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = domhi(1)+1, p_hi(1)
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

      end if            

!*****************************************************************************
! SETTING BL_ZHI
!*****************************************************************************

      if (fix_zhi) then
         do k = domhi(3)+1, p_hi(3)
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,khi)
               end do
            end do
       end do

         if (per_xlo) then
               do k = domhi(3)+1, p_hi(3)
                  do j = jlo,jhi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = domhi(3)+1, p_hi(3)
                  do j = jlo,jhi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = ilo,ihi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = ilo,ihi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = p_lo(1), domlo(1)-1
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = p_lo(1), domlo(1)-1
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = domhi(1)+1, p_hi(1)
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = domhi(1)+1, p_hi(1)
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

      end if            
#endif
#endif

   end subroutine press_fill

end module bc_fill_nd_module
