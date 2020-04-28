
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

#include "mechanism.h"

module prob_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  use fuego_chemistry

  implicit none

  private
  
  public :: amrex_probinit, init_data

contains

! ::: -----------------------------------------------------------
! ::: This routine is called at problem initialization time
! ::: and when restarting from a checkpoint file.
! ::: The purpose is (1) to specify the initial time value
! ::: (not all problems start at time=0.0) and (2) to read
! ::: problem specific data from a namelist or other input
! ::: files and possibly store them or derived information
! ::: in FORTRAN common blocks for later use.
! ::: 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: init      => TRUE if called at start of problem run
! :::              FALSE if called from restart
! ::: strttime <=  start problem with this time variable
! ::: 
! ::: -----------------------------------------------------------

   subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)


      use PeleLM_F,  only: pphys_getP1atm_MKS
      use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
      use mod_Fvar_def, only : pamb
      use extern_probin_module, only: const_viscosity, const_bulk_viscosity, const_conductivity, const_diffusivity, mks_unit
      use probdata_module

      implicit none
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      ! Local
      double precision, dimension(:), allocatable :: data
      integer(kind=8) :: nx, ny, nz
      integer i

      namelist /fortin/ T_mean, P_mean, iname, binfmt, restart, lambda0, &
                        reynolds_lambda0, inres, uin_norm, urms0

                         
      namelist /heattransin/ pamb

!
!      Build `probin' filename -- the name of file containing fortin namelist.
!
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)
      REAL_T :: Yt(NUM_SPECIES), Density_mean, P_mean_in
      integer, parameter :: out_unit=20

      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
!         call bl_abort('probinit called with init ne 1')
      end if

      if (namlen .gt. maxlen) then
         call bl_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         namlen = 6
         probin(1:namlen) = 'probin'
      else
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif


    ! set namelist defaults here
    iname = ""
    binfmt = .false.
    restart = .false.
    lambda0 = 0.5d0
    reynolds_lambda0 = 100.0d0
    inres = 0
    uin_norm = 1.0d0
    urms0 = 1.0d0

! Initial pressure and temperature
    P_mean = 101325.0
    T_mean = 300.d0
    Yt(1) = 0.233d0
    Yt(2) = 0.767d0

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
!     Set defaults
      pamb = pphys_getP1atm_MKS()


      read(untin,fortin)
      
      read(untin,heattransin)
 
      close(unit=untin)

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
      end if


    P_mean_in = P_mean * 10.0d0   ! BEcause the chemkin routine is in CGS
    CALL CKRHOY(P_mean_in,T_mean,Yt,Density_mean)
    Density_mean = Density_mean * 1000.0d0  ! To convert back to MKS

    ! Wavelength associated to Taylor length scale
    k0 = 2.d0/lambda0

   ! Initial density, velocity, and material properties
    tau  = lambda0 / urms0
    const_bulk_viscosity = 0.d0
    const_diffusivity = 0.d0
    const_viscosity = Density_mean * urms0 * lambda0 / reynolds_lambda0 
    const_conductivity = 0.0d0
    mks_unit = .true.

    ! Write this out to file (might be useful for postprocessing)
    if ( amrex_pd_ioprocessor() ) then
       open(unit=out_unit,file="ic.txt",action="write",status="replace")
       write(out_unit,*)"lambda0, k0, rho0, urms0, tau, p0, T0, mu, Reynolds"
       write(out_unit,*) lambda0, "," , k0, "," , Density_mean, "," , urms0, "," , tau, "," , &
            P_mean, "," , T_mean, "," , const_viscosity, "," ,   reynolds_lambda0 
       close(out_unit)
    endif 

    ! Load velocity fields from file. Assume data set ordered in Fortran
    ! format and reshape the data accordingly. One thing to keep in mind
    ! is that this contains the entire input data. We will interpolate
    ! this data later to just match our box. Another assumption is that
    ! the input data is a periodic cube. If the input cube is smaller
    ! than our domain size, the cube will be repeated throughout the
    ! domain (hence the mod operations in the interpolation).

    if (restart) then
       if ( amrex_pd_ioprocessor() ) then
          write(*,*)"Skipping input file reading and assuming restart."
       endif
    else
       nx = int8(inres)
       ny = int8(1)
       nz = int8(1)
       if (dim .ge. 2) then
          ny = int8(inres)
          if (dim .ge. 3) then
             nz = int8(inres)
          endif
       endif

       allocate(data(0:nx*ny*nz*6-1))
       allocate(xinput(0:nx-1,0:ny-1,0:nz-1))
       ! allocate(yinput(0:nx-1,0:ny-1,0:nz-1))
       ! allocate(zinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(uinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(vinput(0:nx-1,0:ny-1,0:nz-1))
       allocate(winput(0:nx-1,0:ny-1,0:nz-1))
       if (binfmt) then
          call read_binary(iname, nx, ny, nz, data)
       else
          call read_csv(iname, nx, ny, nz, data)
       endif

       uinput = urms0 / uin_norm * reshape(data(3::6), (/nx, ny, nz/))
       vinput = urms0 / uin_norm * reshape(data(4::6), (/nx, ny, nz/))
       winput = urms0 / uin_norm * reshape(data(5::6), (/nx, ny, nz/))
       xinput = reshape(data(0::6), (/nx, ny, nz/))
       ! yinput = reshape(data(1::6), (/nx, ny, nz/))
       ! zinput = reshape(data(2::6), (/nx, ny, nz/))


       ! Get the xarray table and the differences.
       allocate(xarray(0:nx-1))
       allocate(xdiff(0:nx-1))
       xarray(0:nx-1) = xinput(:,0,0)
       xdiff(:nx-2) = xarray(1:) - xarray(:nx-2)
       xdiff(nx-1) = xarray(nx-1) - xarray(nx-2)

       ! Dimensions of the input box.
       Linput = maxval(xinput(:,0,0)) + 0.5d0*xdiff(nx-1)

       ! Deallocate some stuff
       deallocate(data)
    endif


  end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  The velocity field you
! ::: provide does not have to be divergence free and the pressure
! ::: field need not be set.  A subsequent projection iteration
! ::: will define aa divergence free velocity field along with a
! ::: consistant pressure.
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nscal     => number of scalar quantities.  You should know
! :::              this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::              ghost region).
! ::: -----------------------------------------------------------

   subroutine init_data(level, time, lo, hi, nscal, &
                        vel, scal, s_lo, s_hi, press, p_lo, p_hi, &
                        delta, xlo, xhi) &
                        bind(C, name="init_data")

      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, domnlo


      use probdata_module

      implicit none

! In/Out
      integer, intent(in) :: level, nscal
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: s_lo(3), s_hi(3)
      integer, intent(in) :: p_lo(3), p_hi(3)
      REAL_T, intent(in)  :: xlo(3), xhi(3)
      REAL_T, intent(in)  :: time, delta(3)
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),dim), intent(out) :: vel
      REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal), intent(out) :: scal
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)), intent(out) :: press

! Local
      REAL_T :: x, y, z, Yl(NUM_SPECIES), Patm
      REAL_T :: dx
      integer :: i, j, k, nspec
      double precision :: xmod, ymod, zmod
      integer :: m, mp1, n, np1, p, pp1
      double precision :: rr, s, t, uinterp, vinterp, winterp
       double precision :: f0, f1, f2, f3, f4, f5, f6, f7


      do k = lo(3), hi(3)
         z = (float(k)+.5d0)*delta(3)+domnlo(3)
         zmod = mod(z,Linput)
         call locate(xarray, inres, zmod, p)
         pp1 = mod(p+1,inres)
         t = (zmod - xarray(p)) / xdiff(p) 
      
         do j = lo(2), hi(2)
            y = (float(j)+.5d0)*delta(2)+domnlo(2)
            ymod = mod(y,Linput)
            call locate(xarray, inres, ymod, n)
            np1 = mod(n+1,inres)
            s = (ymod - xarray(n)) / xdiff(n)

            do i = lo(1), hi(1)
               x = (float(i)+.5d0)*delta(1)+domnlo(1)
               xmod = mod(x,Linput)
               call locate(xarray, inres, xmod, m)
               mp1 = mod(m+1,inres)
               rr = (xmod - xarray(m)) / xdiff(m)

             if (dim .eq. 1) then
                f0 = (1-rr)
                f1 = rr
                uinterp = uinput(m,0,0) * f0 + &
                     uinput(mp1,0,0) * f1
                vinterp = 0.0d0
                winterp = 0.0d0
             elseif (dim .eq. 2) then
                ! Factors for bilinear interpolation
                f0 = (1-rr) * (1-s)
                f1 = rr * (1-s)
                f2 = (1-rr) * s
                f3 = rr * s
                uinterp = uinput(m,n,0) * f0 + &
                     uinput(mp1,n,0) * f1 + &
                     uinput(m,np1,0) * f2 + &
                     uinput(mp1,np1,0) * f3
                vinterp = vinput(m,n,0) * f0 + &
                     vinput(mp1,n,0) * f1 + &
                     vinput(m,np1,0) * f2 + &
                     vinput(mp1,np1,0) * f3
                winterp = 0.0d0
             elseif (dim .eq. 3) then
 


               ! Factors for trilinear interpolation
                f0 = (1-rr) * (1-s) * (1-t)
                f1 = rr * (1-s) * (1-t)
                f2 = (1-rr) * s * (1-t)
                f3 = (1-rr) * (1-s) * t
                f4 = rr * (1-s) * t
                f5 = (1-rr) * s * t
                f6 = rr * s * (1-t)
                f7 = rr * s * t
 
               uinterp = uinput(m,n,p) * f0 + &
                     uinput(mp1,n,p) * f1 + &
                     uinput(m,np1,p) * f2 + &
                     uinput(m,n,pp1) * f3 + &
                     uinput(mp1,n,pp1) * f4 + &
                     uinput(m,np1,pp1) * f5+ &
                     uinput(mp1,np1,p) * f6 + &
                     uinput(mp1,np1,pp1) * f7
                vinterp = vinput(m,n,p) * f0 + &
                     vinput(mp1,n,p) * f1 + &
                     vinput(m,np1,p) * f2 + &
                     vinput(m,n,pp1) * f3 + &
                     vinput(mp1,n,pp1) * f4 + &
                     vinput(m,np1,pp1) * f5+ &
                     vinput(mp1,np1,p) * f6 + &
                     vinput(mp1,np1,pp1) * f7
                winterp = winput(m,n,p) * f0 + &
                     winput(mp1,n,p) * f1 + &
                     winput(m,np1,p) * f2 + &
                     winput(m,n,pp1) * f3 + &
                     winput(mp1,n,pp1) * f4 + &
                     winput(m,np1,pp1) * f5+ &
                     winput(mp1,np1,p) * f6 + &
                     winput(mp1,np1,pp1) * f7
             endif


               scal(i,j,k,Temp) = T_mean
               Yl(1) = 0.233d0
               Yl(2) = 0.767d0

               do nspec = 1,NUM_SPECIES
                  scal(i,j,k,FirstSpec+nspec-1) = Yl(nspec)
               end do

               vel(i,j,k,1) = uinterp
               vel(i,j,k,2) = vinterp
               vel(i,j,k,3) = winterp

            end do

         end do
      end do

      Patm = P_mean / pphys_getP1atm_MKS()

      call pphys_RHOfromPTY(lo,hi, &
                            scal(:,:,:,Density),   s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi, &
                            Patm)
      call pphys_HMIXfromTY(lo,hi, &
                            scal(:,:,:,RhoH),      s_lo, s_hi, &
                            scal(:,:,:,Temp),      s_lo, s_hi, &
                            scal(:,:,:,FirstSpec), s_lo, s_hi)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do nspec = 0,NUM_SPECIES-1
                  scal(i,j,k,FirstSpec+nspec) = scal(i,j,k,FirstSpec+nspec)*scal(i,j,k,Density)
               enddo
               scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
            enddo
         enddo
      enddo

   end subroutine init_data

  ! ::: -----------------------------------------------------------
  ! ::: Read a binary file
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: iname => filename
  ! ::: nx    => input resolution
  ! ::: ny    => input resolution
  ! ::: nz    => input resolution
  ! ::: data  <= output data
  ! ::: -----------------------------------------------------------
  subroutine read_binary(iname,nx,ny,nz,data)

    implicit none

    character(len=255), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny, nz
    double precision, intent(out) :: data(0:nx*ny*nz*6-1)

    integer, parameter :: in_unit=1
    integer :: ios = 0
    integer(kind=8) :: chunk
    integer(kind=8) :: i
    integer(kind=8), parameter :: i_one = 1
    integer(kind=8), parameter :: i_six = 6

    open(in_unit,file=trim(iname), access='stream', form='unformatted', status='old', action='read')

    ! We have to read the file in chunks in case it is bigger than (2^31-1) bytes
    chunk = nx*ny*i_six
    do i = 0, nz-i_one
       read(in_unit,iostat=ios)data(i*chunk:(i+i_one)*chunk-i_one)
    enddo
    close(in_unit)

    if (ios .ne. 0) then
       write(*,*)'Error in binary input file read. Exiting with read error', ios
       stop 99
    endif

  end subroutine read_binary

  ! ::: -----------------------------------------------------------
  ! ::: Read a csv file
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: iname => filename
  ! ::: nx    => input resolution
  ! ::: ny    => input resolution
  ! ::: nz    => input resolution
  ! ::: data  <= output data
  ! ::: -----------------------------------------------------------
  subroutine read_csv(iname,nx,ny,nz,data)

    implicit none

    character(len=255), intent(in) :: iname
    integer(kind=8), intent(in) :: nx, ny, nz
    double precision, intent(out) :: data(0:nx*ny*nz*6-1)

    integer :: i
    integer, parameter :: in_unit=1
    integer :: nlines = 0, ios = 0

    ! Get number of lines in file
    open(in_unit,file=trim(iname), access='sequential', form='formatted', status='old', action='read')
    read(in_unit,*) ! skip header
    nlines = 0
    do
       read(in_unit,*,iostat=ios)
       if (ios .ne. 0) exit
       nlines = nlines + 1
    enddo
    ios = 0

    ! Quick sanity check
    if (nlines .ne. nx*ny*nz) then
       write(*,'("Number of lines in the input file (=",I0,") does not ")')nlines
       write(*,'("  match the input resolution (n=",I0,") in the probin file")')nx
       stop 99
    endif

    ! Read the data from the file
    rewind(in_unit)
    read(in_unit,*) ! skip header
    do i = 0, nlines-1
       read(in_unit, *, iostat=ios)data(i*6:(i+1)*6-1)
       if (ios .ne. 0) then
          write(*,*)'Error in CSV input file read. Exiting with read error', ios
          stop 99
       endif
    enddo
    close(in_unit)

  end subroutine read_csv

  ! ::: -----------------------------------------------------------
  ! ::: Search for the closest index in an array to a given value
  ! ::: using the bisection technique.
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: xtable(0:n-1) => array to search in (ascending order)
  ! ::: n             => number of elements in array
  ! ::: x             => x location
  ! ::: idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
  ! ::: -----------------------------------------------------------
  subroutine locate(xtable, n, x, idxlo)

    implicit none

    double precision, intent(in) :: xtable(0:n-1)
    integer, intent(in) :: n
    double precision, intent(in) :: x
    integer, intent(out) :: idxlo

    ! Local variables
    integer :: idxhi, idxmid
    logical :: notdone

    ! If x is out of bounds, return boundary index
    if (x >= xtable(n-1)) then
       idxlo=n-1
       return
    elseif (x <= xtable(0)) then
       idxlo=0
       return
    endif

    ! Make sure the search array is increasing
    if (xtable(0) > xtable(n-1)) then
       write(*,'("Error in locate: non ascending input search array.")')
       stop 99
    endif

    ! Do the bisection
    idxlo = 0
    idxhi = n-1
    notdone = .true.
    do while (notdone)
       if (idxhi-idxlo <= 1) then
          notdone = .false.
       else
          idxmid = (idxhi+idxlo)/2
          if (x >= xtable(idxmid)) then
             idxlo = idxmid
          else
             idxhi = idxmid
          endif
       endif
    enddo
    return
  end subroutine locate

  ! ::: -----------------------------------------------------------
  ! ::: Calculate the Taylor length scale in x-direction.
  ! :::    lambda = <u^2>/<(du/dx)^2>
  ! :::
  ! ::: INPUTS/OUTPUTS:
  ! :::
  ! ::: u(n,n,n) => input velocity array
  ! ::: x(n,n,n) => input coordinate array
  ! ::: n        => size of array
  ! ::: lambda  <=> lambda
  ! ::: -----------------------------------------------------------
  subroutine calculate_taylor_microscale(u, x, n, lambda)

    implicit none

    double precision, intent(in) :: u(n,n,n)
    double precision, intent(in) :: x(n,n,n)
    integer, intent(in) :: n
    double precision, intent(out) :: lambda

    ! Local variables
    integer :: i
    double precision :: dudx2 = 0, u2 = 0

    ! Square of velocity
    u2 = sum(u(:,:,:)*u(:,:,:))

    ! Calculate the gradients. Assume periodicity for the edges.
    do i = 2, n-1
       dudx2 = dudx2 + sum((u(i+1,:,:) - u(i-1,:,:))**2 / (x(i+1,:,:) - x(i-1,:,:))**2)
    enddo
    dudx2 = dudx2 + sum((u(2,:,:) - u(n,:,:))**2 / (x(2,:,:) - x(n,:,:))**2)
    dudx2 = dudx2 + sum((u(1,:,:) - u(n-1,:,:))**2 / (x(1,:,:) - x(n-1,:,:))**2)

    lambda = sqrt(u2/dudx2)

    return
  end subroutine calculate_taylor_microscale




end module prob_nd_module
