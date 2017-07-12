
#include <MyPROB_F.H>

      subroutine setupbc()
      implicit none
#include <cdwrk.H>
#include <conp.H>
#include <bc.H>
#include <probdata.H>
#include <htdata.H>
      
      REAL_T Patm, FORT_P1ATMMKS
      REAL_T Xt(maxspec), Yt(maxspec)
      integer zone, n, getZone, fuelZone, airZone
      integer b(SDIM)
      integer num_zones_defined
      integer iO2,iCH3OCH3,len
      character*(maxspnml) name
      data  b / 1, 1 /
      
      Patm = pamb / FORT_P1ATMMKS()
      num_zones_defined = 0
      len = len_trim(probtype)

      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then

c     A diffusion flame
         fuelZone = getZone(0.d0, domnlo(2))
         airZone  = getZone(domnlo(1), domnlo(2))
         num_zones_defined = 2

c     Fuel
c     find O2 and DME
         iO2 = -1
         do n=1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'O2') iO2 = n
         end do
         if (iO2.eq.-1)
     &        write(6,*) '.....warning: no O2 in chemistry species list'

         iN2 = -1
         do n=1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'N2') iN2 = n
         end do
         if (iN2.eq.-1)
     &        write(6,*) '.....warning: no N2 in chemistry species list'

         iCH3OCH3 = -1
         do n=1,Nspec
            call get_spec_name(name,n)
            if (name .eq. 'CH3OCH3') iCH3OCH3 = n
         end do
         if (iCH3OCH3.eq.-1)
     &        write(6,*) '.....warning: no CH3OCH3 in chemistry species list'

c      set everything to zero
         do n = 1,Nspec
            Xt(n) = 0.d0
         end do 

c     fuel shall be half O2, half H2 by volume
         Xt(iCH3OCH3) = phi_in
         Xt(iN2) = 1.d0-Xt(iCH3OCH3)

c     convert mole frac to mass frac
         CALL CKXTY (Xt, IWRK, RWRK, Yt)

c     set bc's
         do n=1,Nspec
            Y_bc(n-1,fuelZone) = Yt(n)
         end do
         T_bc(fuelZone) = T_in
         u_bc(fuelZone) = 0.d0
         v_bc(fuelZone) = V_in

c     Air 
         do n=1,Nspec
            Xt(n) = zero
         enddo
         Xt(iN2) = 0.79d0
         Xt(iO2) = 0.21

         CALL CKXTY (Xt, IWRK, RWRK, Yt)         
         do n=1,Nspec
            Y_bc(n-1,airZone) = Yt(n)
         end do
         
         T_bc(airZone) = T_co
         u_bc(airZone) = 0.d0
         v_bc(airZone) = V_co

      else
         
         call bl_pd_abort()

      endif

      do zone=1,num_zones_defined
c     Set density and hmix consistent with data

         call FORT_RHOfromPTY(b, b, 
     &                        rho_bc(zone), DIMARG(b), DIMARG(b),
     &                        T_bc(zone),   DIMARG(b), DIMARG(b),
     &                        Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)
         call FORT_HMIXfromTY(b, b, 
     &                        h_bc(zone),   DIMARG(b), DIMARG(b),
     &                        T_bc(zone),   DIMARG(b), DIMARG(b),
     &                        Y_bc(0,zone), DIMARG(b), DIMARG(b))
      enddo
      bcinit = .true.
      end

c ::: -----------------------------------------------------------
      
      integer function getZone(x, y)
      implicit none
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      REAL_T x, y
      integer len

      getZone = BL_VOLUME
      len     = len_trim(probtype)
         
      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then
         
         if (x .ge. -splitx) then
            getZone = BL_AMBIENT
         else
            getZone = BL_FUELPIPE
         endif
         
      else
         call bl_pd_abort()
      endif
      end
      
c ::: -----------------------------------------------------------
      
      subroutine bcfunction(RegionID,x,y,time,u,v,rho,Yl,T,h,dx,getuv)

      implicit none

      integer RegionID
      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(SDIM)
      REAL_T  Patm, FORT_P1ATMMKS
      logical getuv
      integer dimloc(SDIM)

#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

      integer n, getZone, airZone,fuelZone, len
      REAL_T eta, sigma, pi
      data dimloc / 1,  1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if
      
      len = len_trim(probtype)

      if (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) then

         if (RegionID .eq. BL_YLO) then
            fuelZone = getZone(0.d0, domnlo(2))
            airZone  = getZone(domnlo(1), domnlo(2))

            sigma = 2.5d0*xfrontw*splitx
            eta = 0.5d0 * ( tanh((x + splitx)/sigma) 
     &                    - tanh((x - splitx)/sigma))
                  
            do n=1,Nspec
              Yl(n-1) = Y_bc(n-1,airZone)*(1.d0-eta)
     &               + eta*Y_bc(n-1,fuelZone)
            enddo
            T = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)
c  need to define rho and h from rest of state
c           rho = rho_bc(airZone)*(1.d0-eta) + eta*rho_bc(fuelZone)
c           h = h_bc(airZone)*(1.d0-eta) + eta*h_bc(fuelZone)
      Patm = pamb / FORT_P1ATMMKS()
 
         call FORT_RHOfromPTY(dimloc, dimloc,
     &                        rho, DIMARG(dimloc), DIMARG(dimloc),
     &                        T,   DIMARG(dimloc), DIMARG(dimloc),
     &                        Yl, DIMARG(dimloc), DIMARG(dimloc), Patm)
         call FORT_HMIXfromTY(dimloc, dimloc,
     &                        h,   DIMARG(dimloc), DIMARG(dimloc),
     &                        T,   DIMARG(dimloc), DIMARG(dimloc),
     &                        Yl, DIMARG(dimloc), DIMARG(dimloc))


            if (getuv .eqv. .TRUE.) then

              u = u_bc(airZone)*(1.d0-eta) + eta*u_bc(fuelZone)
              v = v_bc(airZone)*(1.d0-eta) + eta*v_bc(fuelZone)

              if (.not. forceInflow) then
                 ! sinusoidal variation of inflow
                 pi = 4.d0*atan(1.d0)
                 v = v + 100.d0*eta*sin(2.d0*pi*x/(domnhi(1)-domnlo(1)))*sin(2.d0*pi*time/1.d-5)
              end if

            endif

         else
            write(6,*) 'No bcfunction instruction for RegionID = ', RegionID
            call bl_pd_abort(' ')
         endif
      else
         write(6,*) 'No boundary condition for probtype = ', probtype(1:len)
         write(6,*) 'Available: '
         write(6,*) '            ',BL_PROB_JET_DIFFUSION
         call bl_pd_abort(' ')
      endif
      end

c ::: -----------------------------------------------------------
      
      subroutine FORT_INITDATANEWMECH(level,time,lo,hi,nscal,
     &     vel,scal,DIMS(state),press,DIMS(press),
     &     delta,xlo,xhi)
      implicit none
      integer  level, nscal
      integer  lo(SDIM), hi(SDIM)
      integer  DIMDEC(state)
      integer  DIMDEC(press)
      REAL_T   xlo(SDIM), xhi(SDIM)
      REAL_T   time, delta(SDIM)
      REAL_T   vel(DIMV(state),SDIM)
      REAL_T   scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))
 
#include <cdwrk.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>
 
      integer i, j, n
      REAL_T Patm, FORT_P1ATMMKS
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            scal(i,j,Trac) = zero
         end do
      end do
 
      Patm = pamb / FORT_P1ATMMKS()
      call FORT_RHOfromPTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state),
     &     Patm)
      call FORT_HMIXfromTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo
 
      end

c ::: -----------------------------------------------------------
c ::: This routine is called at problem setup time and is used
c ::: to initialize data on each grid.  The velocity field you
c ::: provide does not have to be divergence free and the pressure
c ::: field need not be set.  A subsequent projection iteration
c ::: will define aa divergence free velocity field along with a
c ::: consistant pressure.
c ::: 
c ::: NOTE:  all arrays have one cell of ghost zones surrounding
c :::        the grid interior.  Values in these cells need not
c :::        be set here.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: level     => amr level of grid
c ::: time      => time at which to init data             
c ::: lo,hi     => index limits of grid interior (cell centered)
c ::: nscal     => number of scalar quantities.  You should know
c :::		   this already!
c ::: vel      <=  Velocity array
c ::: scal     <=  Scalar array
c ::: press    <=  Pressure array
c ::: delta     => cell size
c ::: xlo,xhi   => physical locations of lower left and upper
c :::              right hand corner of grid.  (does not include
c :::		   ghost region).
c ::: -----------------------------------------------------------
      subroutine FORT_INITDATA(level,time,lo,hi,nscal,
     &	 	               vel,scal,DIMS(state),press,DIMS(press),
     &                         delta,xlo,xhi)
      implicit none
      integer    level, nscal
      integer    lo(SDIM), hi(SDIM)
      integer    DIMDEC(state)
      integer    DIMDEC(press)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     time, delta(SDIM)
      REAL_T     vel(DIMV(state),SDIM)
      REAL_T    scal(DIMV(state),nscal)
      REAL_T   press(DIMV(press))

#include <cdwrk.H>
#include <conp.H>
#include <htdata.H>
#include <bc.H>
#include <probdata.H>

      integer i, j, n, airZone, fuelZone, getZone, len
      REAL_T x, y, r, Patm, sigma, FORT_P1ATMMKS,eta

      len = len_trim(probtype)
           
      if ((probtype(1:len).eq.BL_PROB_JET_DIFFUSION)) then

         fuelZone = getZone(0.d0, domnlo(2))
         airZone  = getZone(domnlo(1), domnlo(2))

         sigma = 2.5d0*xfrontw*splitx

         do j = lo(2), hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)




c              if (blobr.lt.0) then

                  eta = 0.5d0 * (tanh((x + splitx)/sigma) 
     &                         - tanh((x - splitx)/sigma))
                  
                  do n=1,Nspec
                   scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*(1.d0-eta)
     &                  + eta*Y_bc(n-1,fuelZone)
                  enddo
                  scal(i,j,Temp) = T_bc(airZone)*(1.d0-eta) + eta*T_bc(fuelZone)
                  vel(i,j,1) = u_bc(airZone)*(1.d0-eta) + eta*u_bc(fuelZone)
                  vel(i,j,2) = v_bc(airZone)*(1.d0-eta) + eta*v_bc(fuelZone)
                  scal(i,j,Trac) = 0.d0
c                 
c                 eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
c                 do n=1,Nspec
c                    scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta
c    &                    + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
c                 enddo

c              else
               if (blobr.gt.0) then

                  eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
                  do n=1,Nspec
                     scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta
     &                    + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
                  enddo
               
                  
c     Superimpose blob of hot air
                  r = SQRT((x-blobx)**2 + (y-bloby)**2)
                  eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/Tfrontw))
                  do n=1,Nspec
                     scal(i,j,FirstSpec-1+n) = Y_bc(n-1,airZone)*eta
     &                    + (1.d0-eta)*scal(i,j,FirstSpec-1+n)
                  enddo
                  scal(i,j,Temp) = blobT*eta + (1.d0-eta)*scal(i,j,Temp)
c                 vel(i,j,1) = 0.d0
c                 vel(i,j,2) = v_bc(airZone)
                  
                  vel(i,j,1) = u_bc(airZone)*eta + (1.d0-eta)*u_bc(fuelZone)
                  vel(i,j,2) = v_bc(airZone)*eta + (1.d0-eta)*v_bc(fuelZone)
c                 scal(i,j,Trac) = 0.d0
                  
               endif

            enddo
         enddo
      endif

      Patm = pamb / FORT_P1ATMMKS()

      call FORT_RHOfromPTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state),
     &     Patm)

      call FORT_HMIXfromTY(lo,hi,
     &     scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),
     &     scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))

c     Update typical values
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               typVal_Y(n+1) = MAX(typVal_Y(n+1),scal(i,j,firstSpec+n))
            enddo
            typVal_Density = MAX(scal(i,j,Density),typVal_Density)
            typVal_Temp    = MAX(scal(i,j,Temp),   typVal_Temp)
            typVal_Trac    = MAX(scal(i,j,Trac),   typVal_Trac)
            typVal_RhoH = MAX(ABS(scal(i,j,RhoH)*scal(i,j,Density)),typVal_RhoH)
            do n = 1,BL_SPACEDIM
               typVal_Vel  = MAX(ABS(vel(i,j,n)),typVal_Vel)
            enddo
         enddo
      enddo
      do n = 1,Nspec
         typVal_Y(n) = MIN(MAX(typVal_Y(n),typVal_YMIN),typVal_YMAX)
      enddo

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,Nspec-1
               scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo
      end
      
c ::: -----------------------------------------------------------
c ::: This routine will zero out diffusivity on portions of the
c ::: boundary that are inflow, allowing that a "wall" block
c ::: the complement aperture
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: diff      <=> diffusivity on edges
c ::: DIMS(diff) => index extent of diff array
c ::: lo,hi      => region of interest, edge-based
c ::: domlo,hi   => index extent of problem domain, edge-based
c ::: dx         => cell spacing
c ::: problo     => phys loc of lower left corner of prob domain
c ::: bc         => boundary condition flag (on orient)
c :::                   in BC_TYPES::physicalBndryTypes
c ::: idir       => which face, 0=x, 1=y
c ::: isrz       => 1 if problem is r-z
c ::: id         => index of state, 0=u
c ::: ncomp      => components to modify
c ::: 
c ::: -----------------------------------------------------------
      subroutine FORT_ZEROVISC(diff,DIMS(diff),lo,hi,domlo,domhi,
     &                         dx,problo,bc,idir,isrz,id,ncomp)
      implicit none
      integer DIMDEC(diff)
      integer lo(SDIM), hi(SDIM)
      integer domlo(SDIM), domhi(SDIM)
      integer bc(2*SDIM)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(SDIM)
      REAL_T  problo(SDIM)
      
#include <probdata.H>
#include <cdwrk.H>
#include <htdata.H>
      integer i, j, n, Tid, RHid, YSid, YEid, ys, ye, len
      logical do_T, do_RH, do_Y
      REAL_T xl, xr, xh, y

      len = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) ) then
         Tid  = Temp      - id + SDIM
         RHid = RhoH      - id + SDIM
         YSid = FirstSpec - id + SDIM
         YEid = LastSpec  - id + SDIM
         
         do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
         do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
         ys = MAX(YSid,1)
         ye = MIN(YEid,ncomp)
         do_Y = (ye - ys + 1) .GE. 1
c     
c     Do species, Temp, rhoH
c     
         if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))
     &           .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
            y = DBLE(j)*dx(2)+domnlo(2)
            j = lo(2)
            do i = lo(1), hi(1)
               
               xl = DBLE(i)*dx(1)+domnlo(1) 
               xr = (DBLE(i)+1.d0)*dx(1)+domnlo(1) 
               xh = 0.5d0*(xl+xr)
                  
                  
c                 if (do_T)  diff(i,j,Tid ) = 0.d0
c                 if (do_RH) diff(i,j,RHid) = 0.d0
                  if (do_Y) then
                     do n=ys,ye
                        diff(i,j,n) = 0.d0
                     enddo
                  endif
                     
            end do
         endif
      end if
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: rho       => density array
c ::: DIMS(rho) => index extent of rho array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_DENERROR (tag,DIMS(tag),set,clear,
     &                          rho,DIMS(rho),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &			        problo,time,level)
      implicit none
      integer   DIMDEC(rho)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    rho(DIMV(rho), nvar)

#include <probdata.H>

      call bl_abort('DENERROR: should no be here')
      
      end

c ::: -----------------------------------------------------------

      subroutine FORT_FLAMETRACERROR (tag,DIMS(tag),set,clear,
     &                                ftrac,DIMS(ftrac),lo,hi,nvar,
     &                                domlo,domhi,dx,xlo,
     &			              problo,time,level)
      implicit none
      integer   DIMDEC(ftrac)
      integer   DIMDEC(tag)
      integer   lo(SDIM), hi(SDIM)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    ftrac(DIMV(ftrac), nvar)

      integer   i, j

#include <probdata.H>

      if (level.lt.max_trac_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),
     &              ftrac(i,j,1).gt.flametracval)
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of scalar array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_ADVERROR (tag,DIMS(tag),set,clear,
     &                          adv,DIMS(adv),lo,hi,nvar,
     &                          domlo,domhi,delta,xlo,
     &			        problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(adv)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    delta(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag)), len
      REAL_T    adv(DIMV(adv),nvar)

#include <probdata.H>

      len = len_trim(probtype)

      if ( (probtype(1:len).eq.BL_PROB_JET_DIFFUSION) )then
         call FORT_MVERROR(tag,DIMS(tag),set,clear,
     &                     adv,DIMS(adv),lo,hi,nvar,
     &                     domlo,domhi,delta,xlo,
     &                     problo,time,level)
      endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the
c ::: temperature gradient
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: temp      => density array
c ::: DIMS(temp)=> index extent of temp array
c ::: lo,hi     => index extent of grid
c ::: nvar      => number of components in rho array (should be 1)
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_TEMPERROR (tag,DIMS(tag),set,clear,
     &                          temperature,DIMS(temp),lo,hi,nvar,
     &                          domlo,domhi,dx,xlo,
     &                          problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(temp)
      integer   nvar, set, clear, level
      integer   domlo(SDIM), domhi(SDIM)
      integer   lo(SDIM), hi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    temperature(DIMV(temp),nvar)

      REAL_T    ax, ay, aerr
      integer   i, j, ng

#include <probdata.H>

      ng = min(ARG_H1(temp)-hi(1),ARG_H2(temp)-hi(2),
     &         lo(1)-ARG_L1(temp),lo(2)-ARG_L2(temp))

c      if (ng .lt. 1) then
c         write(6,*) "TEMPERR cannot compute gradient, ng = ",ng
c         call bl_abort(" ")
c      endif
c
c     ::::: refine where there is temperature gradient
c
      if (level .lt. max_temp_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = abs(temperature(i+1,j,1) - temperature(i,j,1))
               ay = abs(temperature(i,j+1,1) - temperature(i,j,1))
               ax = MAX(ax,abs(temperature(i,j,1) - temperature(i-1,j,1)))
               ay = MAX(ay,abs(temperature(i,j,1) - temperature(i,j-1,1)))
               aerr = max(ax,ay)
               tag(i,j) = cvmgt(set,tag(i,j),aerr.ge.tempgrad)
c              tag(i,j) = cvmgt(set,tag(i,j),temperature(i,j,1).lt.temperr)
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: magnitude of vorticity
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: vort      => array of vorticity values
c ::: DIMS(vor) => index extent of vort array
c ::: nvar      => number of components in vort array (should be 1)
c ::: lo,hi     => index extent of grid
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
      subroutine FORT_MVERROR (tag,DIMS(tag),set,clear,
     &                         vort,DIMS(vort),lo,hi,nvar,
     &                         domlo,domhi,dx,xlo,
     &			       problo,time,level)
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(vort)
      integer   nvar, set, clear, level
      integer   lo(SDIM), hi(SDIM)
      integer   domlo(SDIM), domhi(SDIM)
      REAL_T    dx(SDIM), xlo(SDIM), problo(SDIM), time
      integer   tag(DIMV(tag))
      REAL_T    vort(DIMV(vort),nvar)

      integer   i, j

#include <probdata.H>

      if (level .lt. max_vort_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               tag(i,j) = cvmgt(set,tag(i,j),
     &              ABS(vort(i,j,1)).ge.vorterr*2.d0**level)
            enddo
         enddo
      end if

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: den      <=  density array
c ::: DIMS(den) => index extent of den array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of den array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_DENFILL (den,DIMS(den),domlo,domhi,delta,
     &                         xlo,time,bc)
      implicit none

      integer DIMDEC(den), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  den(DIMV(den))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(den)
      lo(2) = ARG_L2(den)
      hi(1) = ARG_H1(den)
      hi(2) = ARG_H2(den)

      call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               den(i,j) = rho
            enddo
         enddo
      endif

      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data and that all non-interior cells have
c ::         have been filled with a large real number.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: adv      <=  advected quantity array
c ::: DIMS(adv) => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of adv array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_ADVFILL (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc)

      implicit none

      integer    DIMDEC(adv)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     delta(SDIM), xlo(SDIM), time
      REAL_T     adv(DIMV(adv))
      integer    bc(SDIM,2)

      integer    i,j
      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(adv)
      lo(2) = ARG_L2(adv)
      hi(1) = ARG_H1(adv)
      hi(2) = ARG_H2(adv)

      call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            do i = lo(1), hi(1)
               adv(i,j) = 0.0d0
            enddo
         enddo
      endif    

      end


c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: temp     <=  temperature array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_TEMPFILL (temp,DIMS(temp),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none

      integer DIMDEC(temp), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  temp(DIMV(temp))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(temp)
      lo(2) = ARG_L2(temp)
      hi(1) = ARG_H1(temp)
      hi(2) = ARG_H2(temp)

      call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               temp(i,j) = T
            enddo
         enddo
      endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: rhoh      <=  rho*h array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_RHOHFILL (rhoh,DIMS(rhoh),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none

      integer DIMDEC(rhoh), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoh(DIMV(rhoh))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      REAL_T  y, x
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoh)
      lo(2) = ARG_L2(rhoh)
      hi(1) = ARG_H1(rhoh)
      hi(2) = ARG_H2(rhoh)

      call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI,x,y,time,u,v,rho,Yl,T,h,delta,.false.)
               rhoh(i,j) = rho*h
            enddo
         enddo
      endif

      end
c
c Fill x & y velocity at once.
c
      subroutine FORT_VELFILL (vel,DIMS(vel),domlo,domhi,delta,
     &                         xlo,time,bc)

      implicit none
      integer DIMDEC(vel), bc(SDIM,2,SDIM)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  vel(DIMV(vel),SDIM)

      call FORT_XVELFILL (vel(ARG_L1(vel),ARG_L2(vel),1),
     & DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,1))

      call FORT_YVELFILL (vel(ARG_L1(vel),ARG_L2(vel),2),
     & DIMS(vel),domlo,domhi,delta,xlo,time,bc(1,1,2))

      end

c
c Fill all chem species at once
c
      subroutine FORT_ALLCHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                             xlo,time,bc)

      implicit none
#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

      integer DIMDEC(rhoY), bc(SDIM,2,Nspec)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY),Nspec)

      integer n
      
      do n=1,Nspec
         call FORT_CHEMFILL (rhoY(ARG_L1(rhoY),ARG_L2(rhoY),n),
     &        DIMS(rhoY),domlo,domhi,delta,xlo,time,bc(1,1,n),n-1)
      enddo
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: xvel     <=  x velocity array
c ::: lo,hi     => index extent of xvel array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_XVELFILL (xvel,DIMS(xvel),domlo,domhi,delta,
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(xvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  xvel(DIMV(xvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
      REAL_T  t_flct, dt_flct
      integer DIMDEC(uflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T uflct(:,:)
      allocatable uflct
#endif


      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(xvel)
      hi(1) = ARG_H1(xvel)
      lo(2) = ARG_L2(xvel)
      hi(2) = ARG_H2(xvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(uflct), loFlctArray, hiFlctArray)
         allocate(uflct(DIMV(uflct)))
         dt_flct = time - tbase_control
         t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
c         write(*,*) 'time,dt_flct,t_flct,tbase_control:',time,dt_flct,t_flct,tbase_control
         call INFL_FILL(FLCT_XVEL, DIMS(uflct), uflct, xlo, delta, t_flct,
     $                  bc, domnlo, domnhi)
      endif
#endif

      
      call filcc (xvel,DIMS(xvel),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 2) then
                  xvel(i,j) = xvel(i,j) + uflct(i,1)*turb_scale
               endif
#endif

            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               xvel(i,j) = u
            enddo
         enddo
      endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: yvel     <=  y velocity array
c ::: lo,hi     => index extent of yvel array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: -----------------------------------------------------------

      subroutine FORT_YVELFILL (yvel,DIMS(yvel),domlo,domhi,delta,
     &                          xlo,time,bc)
      implicit none
      integer DIMDEC(yvel), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM)
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  yvel(DIMV(yvel))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>

#if defined(BL_DO_FLCT)
#include <INFL_FORCE_F.H>
      REAL_T t_flct, dt_flct
      integer DIMDEC(vflct)
      integer loFlctArray(SDIM), hiFlctArray(SDIM)
      REAL_T vflct(:,:)
      allocatable vflct
#endif

      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(yvel)
      hi(1) = ARG_H1(yvel)
      lo(2) = ARG_L2(yvel)
      hi(2) = ARG_H2(yvel)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))

#if defined(BL_DO_FLCT)
      if (forceInflow) then
         do i = 1, SDIM
            loFlctArray(i) = lo(i)
            hiFlctArray(i) = hi(i)
         enddo
         loFlctArray(strmwse_dir) = 1
         hiFlctArray(strmwse_dir) = 1
         call SET_ARGS(DIMS(vflct), loFlctArray, hiFlctArray)
         allocate(vflct(DIMV(vflct)))
         dt_flct = time - tbase_control
         t_flct = zbase_control + V_in*dt_flct + dV_control*dt_flct**2
         call INFL_FILL(FLCT_YVEL, DIMS(vflct), vflct, xlo, delta, t_flct,
     $                  bc, domnlo, domnhi)
      endif
#endif

      
      call filcc (yvel,DIMS(yvel),domlo,domhi,delta,xlo,bc)

c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
#if defined(BL_DO_FLCT)
               if (forceLo .and. strmwse_dir .eq. 2) then
                  yvel(i,j) = yvel(i,j) + vflct(i,1)*turb_scale
               endif
#endif

            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, h, delta,.true.)
               yvel(i,j) = v
            enddo
         enddo
      endif
      end
      
c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.
c :::
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: rhoY      <= rho*Y (Y=mass fraction) array
c ::: lo,hi     => index extent of adv array
c ::: domlo,hi  => index extent of problem domain
c ::: delta     => cell spacing
c ::: xlo       => physical location of lower left hand
c :::              corner of temperature array
c ::: time      => problem evolution time
c ::: bc        => array of boundary flags bc(BL_SPACEDIM,lo:hi)
c ::: stateID   => id index of state being filled
c ::: -----------------------------------------------------------
      
      subroutine FORT_CHEMFILL (rhoY,DIMS(rhoY),domlo,domhi,delta,
     &                          xlo,time,bc,id )
      implicit none
      integer DIMDEC(rhoY), bc(SDIM,2)
      integer domlo(SDIM), domhi(SDIM), id
      REAL_T  delta(SDIM), xlo(SDIM), time
      REAL_T  rhoY(DIMV(rhoY))

#include <cdwrk.H>
#include <bc.H>
#include <probdata.H>
      
      integer i, j
      integer ilo, ihi, jlo, jhi
      REAL_T  y, x, hx
      REAL_T  u, v, rho, Yl(0:maxspec-1), T, h

      integer lo(SDIM), hi(SDIM)

      lo(1) = ARG_L1(rhoY)
      hi(1) = ARG_H1(rhoY)
      lo(2) = ARG_L2(rhoY)
      hi(2) = ARG_H2(rhoY)

      hx  = delta(1)
      ilo = max(lo(1),domlo(1))
      ihi = min(hi(1),domhi(1))
      jlo = max(lo(2),domlo(2))
      jhi = min(hi(2),domhi(2))
      
      call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)
      
c     NOTE:
c     In order to set Dirichlet boundary conditions in a mulitspecies
c     problem, we have to know all the state values, in a sense.  For
c     example, the total density rho = sum_l(rho.Yl).  So to compute any
c     rho.Yl, we need all Yl's...also need to evaluate EOS since we
c     really are specifying T and Yl's.  so, all this is centralized
c     here.  Finally, a layer of flexibilty is added to for the usual case
c     that the bc values may often be set up ahead of time.

      if (bc(1,1).eq.EXT_DIR.and.lo(1).lt.domlo(1)) then
         do i = lo(1), domlo(1)-1
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XLO, x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif
      
      if (bc(1,2).eq.EXT_DIR.and.hi(1).gt.domhi(1)) then
         do i = domhi(1)+1, hi(1)
            x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
            do j = lo(2), hi(2)
               y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
               call bcfunction(BL_XHI, x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    

      if (bc(2,1).eq.EXT_DIR.and.lo(2).lt.domlo(2)) then
         do j = lo(2), domlo(2)-1
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YLO, x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif    
      
      if (bc(2,2).eq.EXT_DIR.and.hi(2).gt.domhi(2)) then
         do j = domhi(2)+1, hi(2)
            y = (DBLE(j)+.5d0)*delta(2)+domnlo(2)
            do i = lo(1), hi(1)
               x = (DBLE(i)+.5d0)*delta(1)+domnlo(1)
               call bcfunction(BL_YHI, x, y, time, u, v, rho, Yl, T, h, delta,.false.)
               rhoY(i,j) = rho*Yl(id)
            enddo
         enddo
      endif
      
      end

c ::: -----------------------------------------------------------
c ::: This routine is called during a filpatch operation when
c ::: the patch to be filled falls outside the interior
c ::: of the problem domain.  You are requested to supply the
c ::: data outside the problem interior in such a way that the
c ::: data is consistant with the types of the boundary conditions
c ::: you specified in the C++ code.  
c ::: 
c ::: NOTE:  you can assume all interior cells have been filled
c :::        with valid data.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: p        <=  pressure array
c ::: DIMS(p)   => index extent of p array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of rho array
c ::: time      => problem evolution time
c ::: bc	=> array of boundary flags bc(BL_SPACEDIM,lo:hi) 
c ::: -----------------------------------------------------------

      subroutine FORT_PRESFILL (p,DIMS(p),domlo,domhi,dx,xlo,time,bc)
      implicit none
      integer    DIMDEC(p)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     dx(SDIM), xlo(SDIM), time
      REAL_T     p(DIMV(p))
      integer    bc(SDIM,2)

      integer    i, j
      integer    ilo, ihi, jlo, jhi
      logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi
      logical    per_xlo, per_xhi, per_ylo, per_yhi

      fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
      per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
      fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
      per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
      fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
      per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
      fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
      per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)

      ilo = max(ARG_L1(p),domlo(1))
      ihi = min(ARG_H1(p),domhi(1))
      jlo = max(ARG_L2(p),domlo(2))
      jhi = min(ARG_H2(p),domhi(2))
c
c     ::::: left side
c
      if (fix_xlo) then
         do i = ARG_L1(p), domlo(1)-1
            do j = jlo,jhi
               p(i,j) = p(ilo,j)
            end do
         end do
         if (fix_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,jlo)
               end do
            end do
         else if (per_ylo) then
            do i = ARG_L1(p), domlo(1)-1
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
         if (fix_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,jhi)
               end do
            end do
         else if (per_yhi) then
            do i = ARG_L1(p), domlo(1)-1
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ilo,j)
               end do
            end do
         end if
      end if
c
c     ::::: right side
c
      if (fix_xhi) then
         do i = domhi(1)+1, ARG_H1(p)
            do j = jlo,jhi
               p(i,j) = p(ihi,j)
            end do
	 end do
	 if (fix_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,jlo)
               end do
	    end do
	 else if (per_ylo) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = ARG_L2(p), domlo(2)-1
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
	 if (fix_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,jhi)
               end do
	    end do
	 else if (per_yhi) then
	    do i = domhi(1)+1, ARG_H1(p)
               do j = domhi(2)+1, ARG_H2(p)
                  p(i,j) = p(ihi,j)
               end do
	    end do
         end if
      end if
      
      if (fix_ylo) then
         do j = ARG_L2(p), domlo(2)-1
            do i = ilo, ihi
               p(i,j) = p(i,jlo)
            end do
	 end do
	 if (per_xlo) then
          do j = ARG_L2(p), domlo(2)-1
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
	 if (per_xhi) then
           do j = ARG_L2(p), domlo(2)-1
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jlo)
               end do
	    end do
         end if
      end if

      if (fix_yhi) then
         do j = domhi(2)+1, ARG_H2(p)
            do i = ilo, ihi
               p(i,j) = p(i,jhi)
            end do
	 end do
	 if (per_xlo) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = ARG_L1(p), domlo(1)-1
                  p(i,j) = p(i,jhi)
               end do
	    end do
         end if
	 if (per_xhi) then
	    do j = domhi(2)+1, ARG_H2(p)
               do i = domhi(1)+1, ARG_H1(p)
                  p(i,j) = p(i,jhi)
               end do
	    end do
         end if
      end if

      end

      subroutine FORT_RADLOSS(lo,hi,rad,DIMS(rad),
     &                        T,DIMS(T),Y,DIMS(Y),dx,Patm,time)
      implicit none

#include <cdwrk.H>
#include <probdata.H>

      integer DIMDEC(rad)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer lo(SDIM), hi(SDIM)
      REAL_T  rad(DIMV(rad))
      REAL_T  T(DIMV(T))
      REAL_T  Y(DIMV(Y),1)
      REAL_T  dx(SDIM), Patm, time

      integer i, j
      
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rad(i,j) = zero
         end do
      end do
      end

c
c
c ::: -----------------------------------------------------------
c
c     This routine add the forcing terms to the momentum equation
c
      subroutine FORT_MAKEFORCE(time,force,rho,
     &                          DIMS(istate),DIMS(state),
     &                          dx,xlo,xhi,gravity,scomp,ncomp)

      implicit none

      integer    DIMDEC(state)
      integer    DIMDEC(istate)
      integer    scomp, ncomp
      REAL_T     time, dx(SDIM)
      REAL_T     xlo(SDIM), xhi(SDIM)
      REAL_T     force  (DIMV(istate),scomp+1:scomp+ncomp)
      REAL_T     rho    (DIMV(state))
      REAL_T     gravity

#include <probdata.H>
#include <cdwrk.H>
#include <bc.H>

      integer i, j, n
      integer ilo, jlo
      integer ihi, jhi
      REAL_T  hx, hy
      integer isioproc
      integer nXvel, nYvel, nRho, nTrac

      call bl_pd_is_ioproc(isioproc)

      if (isioproc.eq.1 .and. pseudo_gravity.eq.1) then
         write(*,*) "pseudo_gravity::dV_control = ",dV_control
      endif

      hx = dx(1)
      hy = dx(2)

      ilo = istate_l1
      jlo = istate_l2
      ihi = istate_h1
      jhi = istate_h2

c     Assumes components are in the following order
      nXvel = 1
      nYvel = 2
      nRho  = 3
      nTrac = 4

      if (scomp.eq.0) then
         if (abs(gravity).gt.0.0001) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = gravity*rho(i,j)
               enddo
            enddo
c     else to zero
         else
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nXvel) = zero
                  force(i,j,nYvel) = zero
               enddo
            enddo
         endif
c     Add the pseudo gravity afterwards...
         if (pseudo_gravity.eq.1) then
            do j = jlo, jhi
               do i = ilo, ihi
                  force(i,j,nYvel) = force(i,j,nYvel) + dV_control*rho(i,j)
               enddo
            enddo
         endif
c     End of velocity forcing
      endif
      
      if ((scomp+ncomp).gt.BL_SPACEDIM) then
c     Scalar forcing
         do n = max(scomp+1,nRho), scomp+ncomp
            if (n.eq.nRho) then
c     Density
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else if (n.eq.nTrac) then
c     Tracer
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            else
c     Other scalar
               do j = jlo, jhi
                  do i = ilo, ihi
                     force(i,j,n) = zero
                  enddo
               enddo
            endif
         enddo
      endif

      end

