#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define CONPF_FILE conpFY
#define CONPJ_FILE conpJY

#if defined(BL_USE_FLOAT) || defined(BL_CRAY)
#define twothousand 2000
#define one100th    0.01
#define ten2minus19 1.e-19
#define million     1.e6
#define one2minus3  1.e-3
#else
#define twothousand 2000d0
#define one100th    0.01d0
#define ten2minus19 1.d-19
#define million     1.d6
#define one2minus3  1.d-3
#endif

#define SDIM 3

      subroutine FORT_NORMMASS(lo, hi, xsID,
     &                         Y, DIMS(Y), Ynorm, DIMS(YNORM))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(Ynorm)
      integer xsID
      REAL_T Y(DIMV(Y),*)
      REAL_T Ynorm(DIMV(Ynorm),*)

      integer i, j, k, n
      REAL_T sum

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sum = zero
               do n=1,Nspec
                  Ynorm(i,j,k,n) =  MAX( Y(i,j,k,n),zero)
                  sum = sum + Ynorm(i,j,k,n)
               end do
               Ynorm(i,j,k,xsID) = Y(i,j,k,xsID)+ one - sum
            end do
         end do
      end do
      end

      subroutine FORT_FRrateXTP(lo,hi,X,DIMS(X),T,DIMS(T),
     &                          FwdK,DIMS(FwdK),RevK,DIMS(RevK),
     &                          Patm,rxns,Nrxns)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(T)
      integer DIMDEC(FwdK)
      integer DIMDEC(RevK)
      integer Nrxns
      integer rxns(Nrxns)
      REAL_T X(DIMV(X),*)
      REAL_T T(DIMV(T))
      REAL_T FwdK(DIMV(FwdK),*)
      REAL_T RevK(DIMV(RevK),*)
      REAL_T Patm, scale

      REAL_T Xt(maxspec),FwdKt(maxreac),RevKt(maxreac)
      REAL_T sum, Yt(maxspec)
      integer i,j,k,n
      REAL_T P1atm,RU,RUC,Pdyne

      parameter(scale = million)

      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      Pdyne = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
#ifdef DO_JBB_HACK
               CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
               sum = zero
               do n=1,Nspec
                  Yt(n) =MAX( Yt(n),zero)
                  sum = sum+Yt(n)
               end do
               if (iN2 .gt. 0) then
                  Yt(iN2) = Yt(iN2)+one-sum
               endif
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
#endif
               CALL CKKFKR(Pdyne,T(i,j,k),Xt,IWRK(ckbi),RWRK(ckbr),FwdKt,RevKt)
               do n=1,Nrxns
                  FwdK(i,j,k,n) = FwdKt(rxns(n)+1)*scale
                  RevK(i,j,k,n) = RevKt(rxns(n)+1)*scale
               end do
            end do
         end do
      end do
      end

      subroutine FORT_HTRLS(lo,hi,Y,DIMS(Y),T,DIMS(T),
     &                      Q,DIMS(Q),Patm)
      implicit none

#include "cdwrk.H"
#include "conp.H"

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(Q)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T Q(DIMV(Q))
      REAL_T Patm

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      integer ndummy
      REAL_T tdummy,P1atm,RU,RUC
      REAL_T RHO_CGS, CPB_CGS, scal, HMIX_CGS, H_CGS(maxspec)

      ndummy = Nspec
      tdummy = 0.
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

#ifdef LMC_SDC
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               
               do n=1,Nspec
                  Zt(n) = Y(i,j,k,n)
               end do
               call CKHBMS(T(i,j,k),Zt,IWRK(ckbi),RWRK(ckbr),HMIX_CGS)
               
               Zt(Nspec+1) = HMIX_CGS * 1.d-4

               call CKRHOY(RWRK(NP),T(i,j,k), Zt(1), IWRK(ckbi), RWRK(ckbr), RHO_CGS)
               do n=1,Nspec
                  Zt(n) = Zt(n) * (RHO_CGS * 1.d3)
               end do
               c_0(1:Nspec+1) = zero
               rhoY_INIT(1:Nspec) = Zt(1:Nspec)
               rhoH_INIT = HMIX_CGS * 1.d-4
               T_cell = T(i,j,k)
               call conpFY_sdc(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               
               call CKHMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),H_CGS)
               Q(i,j,k) = 0.d0
               do n= 1, Nspec
                  Q(i,j,k) = Q(i,j,k) - Zdott(n) * H_CGS(n) * 1.d-4
               end do

            end do
         end do
      end do

#else

c     NOTE: scal converts result from assumed cgs to MKS (1 erg/s.cm^3 = .1 J/s.m^3)
      scal = tenth
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               CALL CKRHOY(RWRK(NP),Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),RHO_CGS)
               CALL CKCPBS(Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),CPB_CGS)
               Q(i,j,k) = Zdott(1) * scal * (RHO_CGS * 1.d3) * (CPB_CGS * 1.d-4)
            end do
         end do
      end do
#endif
      end

      subroutine FORT_RRATEY(lo,hi,Y,DIMS(Y),T,DIMS(T),
     &                       Ydot,DIMS(Ydot),Patm)
      implicit none

#include "cdwrk.H"
#include "conp.H"

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(Ydot)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T Ydot(DIMV(Ydot),*)
      REAL_T Patm

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      integer ndummy
      REAL_T tdummy,P1atm,RU,RUC, HMIX_CGS, RHO_CGS, RHOINV_MKS

      ndummy = Nspec
      tdummy = 0.
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

#ifdef LMC_SDC
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               
               do n=1,Nspec
                  Zt(n) = Y(i,j,k,n)
               end do
               call CKHBMS(T(i,j,k),Zt,IWRK(ckbi),RWRK(ckbr),HMIX_CGS)
               call CKRHOY(RWRK(NP),T(i,j,k), Zt(1), IWRK(ckbi), RWRK(ckbr), RHO_CGS)
               
               Zt(Nspec+1) = (HMIX_CGS * 1.d-4) * (RHO_CGS * 1.d3)
               
               do n=1,Nspec
                  Zt(n) = Zt(n) * (RHO_CGS * 1.d3)
               end do
               c_0(1:Nspec+1) = zero
               rhoY_INIT(1:Nspec) = Zt(1:Nspec)
               rhoH_INIT = (HMIX_CGS * 1.d-4) * (RHO_CGS * 1.d3)
               T_cell = T(i,j,k)
               call conpFY_sdc(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               
               RHOINV_MKS = 1.d0 / (RHO_CGS * 1.d3)
               do n=1,Nspec
                  Ydot(i,j,k,n) = Zdott(n) * RHOINV_MKS
               end do
            end do
         end do
      end do
#else
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               do n=1,Nspec
                  Ydot(i,j,k,n) = Zdott(n+1)
               end do
            end do
         end do
      end do
#endif
      end

#ifdef LMC_SDC
c     used by LMC
      subroutine FORT_RRATERHOY(lo,hi,RhoY,DIMS(RhoY),RhoH,DIMS(RhoH),T,DIMS(T),
     &                          RhoYdot,DIMS(RhoYdot))
      implicit none

#include "cdwrk.H"
#include "conp.H"

      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(RhoY)
      integer DIMDEC(RhoH)
      integer DIMDEC(T)
      integer DIMDEC(RhoYdot)
      REAL_T RhoY(DIMV(RhoY),Nspec)
      REAL_T RhoH(DIMV(RhoH))
      REAL_T T(DIMV(T))
      REAL_T RhoYdot(DIMV(RhoYdot),Nspec)

      REAL_T Zt(maxspec+1),Zdott(maxspec+1)
      integer i,j,k,n
      REAL_T TIME,P1atm,RU,RUC

      TIME = zero
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      !
      ! Note that c_0,c_1,rhoh_INIT & T_cell are thread-private.
      !
!$omp parallel do private(i,j,k,n,Zt,Zdott)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               c_0(1:Nspec+1) = zero
               c_1(1:Nspec+1) = zero
               rhoh_INIT      = RhoH(i,j,k)
               T_cell         = T(i,j,k)
               Zt(Nspec+1)    = RhoH(i,j,k)
               do n=1,Nspec
                  Zt(n) = RhoY(i,j,k,n)
               end do
               call conpFY_sdc(NEQ,TIME,Zt,Zdott,RWRK,IWRK)
               do n=1,Nspec
                  RhoYdot(i,j,k,n) = Zdott(n)
               end do
            end do
         end do
      end do
!$omp end parallel do
      end
#endif

      subroutine FORT_MASSTOMOLE(lo, hi, Y, DIMS(Y), X, DIMS(X))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(X)
      REAL_T Y(DIMV(Y),*)
      REAL_T X(DIMV(X),*)

      REAL_T Xt(maxspec), Yt(maxspec)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end

c     used by LMC
      subroutine FORT_MOLETOMASS(lo, hi, X, DIMS(X), Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(X)
      integer DIMDEC(Y)
      REAL_T X(DIMV(X),*)
      REAL_T Y(DIMV(Y),*)
      
      REAL_T Xt(maxspec), Yt(maxspec)
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
               CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
               do n = 1,Nspec
                  Y(i,j,k,n) = Yt(n)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MASSTP_TO_CONC(lo, hi, Patm,
     &                           Y, DIMS(Y), T, DIMS(T), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      REAL_T Patm
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      
      REAL_T Yt(maxspec), Ct(maxspec), RU, RUC, P1ATM, Ptmp
      integer i,j,k,n

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTCP(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*million
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MASSR_TO_CONC(lo, hi, Y, DIMS(Y), 
     &                              T, DIMS(T), RHO, DIMS(RHO), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(C)
      integer DIMDEC(RHO)
      REAL_T Y(DIMV(Y),*)
      REAL_T T(DIMV(T))
      REAL_T C(DIMV(C),*)
      REAL_T RHO(DIMV(RHO))

      REAL_T Yt(maxspec), Ct(maxspec), rhoScl
      integer i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               rhoScl = RHO(i,j,k)*one2minus3
               CALL CKYTCR(rhoScl,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*million
               end do
            end do
         end do
      end do
      end

      subroutine FORT_CONC_TO_MOLE(lo, hi,
     &                             C, DIMS(C), X, DIMS(X))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(C)
      integer DIMDEC(X)
      REAL_T C(DIMV(C),*)
      REAL_T X(DIMV(X),*)

      REAL_T Ct(maxspec), Xt(maxspec), scale
      integer i,j,k,n

      parameter(scale = one/million)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*scale
               end do
               CALL CKCTX(Ct,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end

      subroutine FORT_MOLPROD(lo, hi, id, 
     &                        Q, DIMS(Q), C, DIMS(C), T, DIMS(T) )
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM), id
      integer DIMDEC(Q)
      integer DIMDEC(C)
      integer DIMDEC(T)
      REAL_T Q(DIMV(Q),*)
      REAL_T C(DIMV(C),*)
      REAL_T T(DIMV(T))

      REAL_T Ct(maxspec), Qt(maxreac), Qkt(maxreac), millionth
      integer i,j,k,n

      millionth = one/million
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*millionth
               end do
               CALL CKQC(T(i,j,k),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
#ifdef MIKE
               CALL CKCONT(id,Qt,IWRK(ckbi),RWRK(ckbr),Qkt)
               do n = 1,Nreac
                  Q(i,j,k,n) = Qkt(n)*million
               end do
#else
               call bl_abort('FORT_MOLPROD: CKCONT not available')
#endif
            end do
         end do
      end do
      end

      subroutine FORT_GETELTMOLES(namenc, namlen, lo, hi,
     &                            Celt, DIMS(Celt), C, DIMS(C))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer namlen, maxlen
      integer namenc(namlen)
      integer DIMDEC(Celt)
      integer DIMDEC(C)
      REAL_T Celt(DIMV(Celt))
      REAL_T C(DIMV(C),*)
      integer thenames(maxelts*2)
      logical match
      integer i, j, k, theidx, n, lout
      integer NCF(Nelt,Nspec)
c     Find index of desired element
      CALL CKSYME(thenames,2)
      theidx = -1
      do i=1,Nelt
         match = .true.
         do j=1,namlen               
            if (namenc(j) .NE. thenames((i-1)*2+j)) match = .false.
         enddo
         if (match .eqv. .true.) theidx = i
      end do
      if (theidx.lt.0) then
         call bl_pd_abort()
      endif
c     Get the matrix of elements versus species
      call CKNCF(Nelt,IWRK,RWRK,NCF)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Celt(i,j,k) = zero
               do n = 1,Nspec
                  Celt(i,j,k) = Celt(i,j,k) +C(i,j,k,n)*NCF(theidx,n)
               end do
            end do
         end do
      end do
      end

c     use by LMC
      integer function FORT_CONPSOLV(lo, hi,
     &                         Ynew, DIMS(Ynew), Tnew, DIMS(Tnew),
     &                         Yold, DIMS(Yold), Told, DIMS(Told),
     &                         FuncCount, DIMS(FuncCount),
     &                         Patm,
     &                         dt,
     &                         diag,do_diag,do_stiff)
      implicit none

#include "cdwrk.H"
#include "conp.H"
#include "vode.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Yold)
      integer DIMDEC(Told)
      integer DIMDEC(Ynew)
      integer DIMDEC(Tnew)
      integer DIMDEC(FuncCount)
      integer do_diag, do_stiff
      REAL_T Yold(DIMV(Yold),*)
      REAL_T Told(DIMV(Told))
      REAL_T Ynew(DIMV(Ynew),*)
      REAL_T Tnew(DIMV(Tnew))
      REAL_T FuncCount(DIMV(FuncCount))
      REAL_T Patm, dt
      REAL_T diag(DIMV(FuncCount),*)
      
      integer ITOL, IOPT, ITASK, open_vode_failure_file
      parameter (IOPT=1, ITASK=1)
      REAL_T RTOL, ATOL(maxspec+1),ATOLEPS
      REAL_T spec_scalT, NEWJ_TOL
      parameter (spec_scalT=twothousand, NEWJ_TOL=one100th)
      external CONPF_FILE, CONPJ_FILE, open_vode_failure_file
      REAL_T TT1, TT2, RU, RUC, P1atm, sum, atoln
      integer i, j, k, m, MF, ISTATE, LOUTCK, lout, nfails

      integer nsub, node
      REAL_T dtloc, weight, TT1save, scale, tspecies(0:maxspec), nwdot(0:maxspec)
      REAL_T Ct(maxspec),Qt(maxreac),Ytemp(maxspec),Yres(maxspec)

      character*(maxspnml) name
      parameter (LOUTCK=6)

      logical newJ_triggered, bad_soln
      !
      ! Set to .true. if you want to see the 'nfails' output
      !
      logical, parameter :: verbose = .false.

      ITOL    = vode_itol
      RTOL    = vode_rtol
      ATOLEPS = vode_atol
      !
      ! Set molecular weights and pressure in area accessible by conpF
      !
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)

      if (do_stiff .eq. 1) then
         MF = 22  ! finite difference jacobian
      else
         MF = 10
      endif
      !
      ! Set up ATOL
      !
      if (ITOL.eq.2 .or. ITOL.eq.4) then
         ATOL(1) = ATOLEPS*typVal_Temp
         if (ATOL(1) <= 0.0d0) call bl_abort('typVal_Temp <= 0')
         do m=1,Nspec
            ATOL(m+1) = ATOLEPS*typVal_Y(m)
            if (ATOL(m+1) <= 0.0d0) call bl_abort('typVal_Y <= 0')
         end do
      else
         ATOL(1) = ATOLEPS
      endif

      if (do_diag.eq.1) then
         nsub  = nchemdiag
         dtloc = dt/nchemdiag
      else
         nsub  = 1
         dtloc = dt
      endif

      nfails = 0

      FORT_CONPSOLV = 1

!$omp parallel
      !
      ! Force recalculation of jacobian for each XYZ block.
      !
      FIRST = .TRUE.

      RWRK(NP) = Patm * P1atm

!$omp do private(i,j,k,m,TT1,ISTATE,sum,atoln,Ytemp,Yres,newJ_triggered)
!$omp&private(scale,bad_soln,Ct,Qt,node,weight,TT1save,TT2,tspecies)
!$omp&private(lout,name)
#ifdef REGRESSIONTEST
!$omp&schedule(static)
#else
!$omp&schedule(dynamic,1)
#endif
      do k=lo(3),hi(3)
         if (FORT_CONPSOLV .eq. 0) cycle
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               TT1    = zero
               ISTATE = 1
#ifdef SOLN_IS_1D
            if (i.ne.lo(1) .and. j.ne.lo(2)) then
               do m=1,Nspec
                  Ynew(i,j,k,m) = Ynew(lo(1),lo(2),k,m)
               end do
               Tnew(i,j,k) = Tnew(lo(1),lo(2),k)
            else
#endif
               tspecies(0) = Told(i,j,k)

#ifdef DO_JBB_HACK
#ifdef DO_JBB_HACK_TEMP
               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yold(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) Ytemp(iN2) = Ytemp(iN2)+one-sum
#else
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
            end do
#endif
            do m = 1,Nspec
               tspecies(m) = Ytemp(m)
            end do

#ifdef ALWAYS_NEW_J
            FIRST = .TRUE.
#endif
            newJ_triggered = FIRST

#ifdef TRIGGER_NEW_J
            if (.NOT. FIRST) then
               sum = zero
               do m=1,NEQ
                  scale = spec_scalT
                  if (m.ne.1) scale = spec_scalY(m-1)
                  sum = sum + ABS(tspecies(m-1)-YJ_SAVE(m))/scale
               end do
               if (sum .gt. NEWJ_TOL) then
                  FIRST = .TRUE.
                  newJ_triggered = .TRUE.
               endif
            endif
#endif
            if (do_diag.eq.1) then
               FuncCount(i,j,k) = 0
               CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
               CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
               do m=1,Nreac
                  diag(i,j,k,m) = diag(i,j,k,m)+half*dtloc*Qt(m)*million
               enddo
            endif

            do node = 1,nsub
               if (node.lt.nsub) then
                  weight = one
               else
                  weight = half
               endif

               TT1save = TT1
               TT2     = TT1 + dtloc

#if !defined(BL_USE_DOUBLE)
               CALL SVODE
#else 
               CALL DVODE
#endif
     &              (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &              ITASK, ISTATE, IOPT, RWRK(dvbr), dvr, IWRK(dvbi), dvi,
     &              CONPJ_FILE, MF, RWRK, IWRK)
               !
               ! If the step was bad, and we reused an old Jacobian,
               ! try again from scratch
               !
#if defined(TRIGGER_NEW_J) && defined(DO_JBB_HACK)
               if ((ISTATE .LT. 0) .and. (.NOT. newJ_triggered)) then
                  bad_soln = .FALSE.
                  do m=1,Nspec
                     if (tspecies(m) .lt. -1.D-6*spec_scalY(m))
     &                    bad_soln = .TRUE.
                  end do
                  if (bad_soln) then
                     FIRST  = .TRUE. 
                     TT1    = zero
                     ISTATE = 1
                     tspecies(0) = Told(i,j,k)
                     do m=1,Nspec
                        tspecies(m) = Ytemp(m)
                     end do
                     if (verbose) then
!$omp atomic
                        nfails = nfails + 1
!$omp end atomic
                     end if

#if !defined(BL_USE_DOUBLE)
                  CALL SVODE
#else 
                  CALL DVODE
#endif
     &           (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &              ITASK, ISTATE, IOPT, RWRK(dvbr), dvr, IWRK(dvbi), dvi,
     &              CONPJ_FILE, MF, RWRK, IWRK)
               end if
            end if
#endif
               TT1 = TT2

               if (do_diag.eq.1) then
                  CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
                  CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
                  do m=1,Nreac
                     diag(i,j,k,m) = diag(i,j,k,m)+weight*dtloc*Qt(m)*million
                  enddo
                  FuncCount(i,j,k) = FuncCount(i,j,k) + IWRK(dvbi+11)
               else
                  FuncCount(i,j,k) = IWRK(dvbi+11)
               endif

               if (verbose_vode .eq. 1) then
!$omp critical(output)
                  write(6,*) '......dvode done:'
                  write(6,*) ' last successful step size = ',RWRK(dvbr+10)
                  write(6,*) '          next step to try = ',RWRK(dvbr+11)
                  write(6,*) '   integrated time reached = ',RWRK(dvbr+12)
                  write(6,*) '      number of time steps = ',IWRK(dvbi+10)
                  write(6,*) '              number of fs = ',IWRK(dvbi+11)
                  write(6,*) '              number of Js = ',IWRK(dvbi+12)
                  write(6,*) '    method order last used = ',IWRK(dvbi+13)
                  write(6,*) '   method order to be used = ',IWRK(dvbi+14)
                  write(6,*) '            number of LUDs = ',IWRK(dvbi+18)
                  write(6,*) ' number of Newton iterations ',IWRK(dvbi+19)
                  write(6,*) ' number of Newton failures = ',IWRK(dvbi+20)
                  if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                     call get_spec_name(name,IWRK(dvbi+15))
                     write(6,*) '   spec with largest error = ',name
                  end if
!$omp end critical(output)
               end if

               if (ISTATE .LE. -1) then
!$omp critical(output)
                  call CONPF_FILE(NEQ, TT1, tspecies, nwdot, RWRK, IWRK)
                  lout = open_vode_failure_file()                  
                  write(lout,*)
                  write(lout,995) 'VODE Failed at (i,j,k) = (',i,',',j,
     &                 ',',k,'),   Return code = ',ISTATE
                  write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                write(lout,995) 'State ID, old, last, dY/dt, dY/dt*(dt)'
                  write(lout,996) 'T               ',
     &                 Told(i,j,k),tspecies(0),nwdot(0),nwdot(0)*(dt-TT1)
                  do m=1,Nspec
                     call get_spec_name(name,m)
                     write(lout,996) name,Yold(i,j,k,m),
     &                    tspecies(m),nwdot(m),nwdot(m)*(dt-TT1)
                  end do

 995              format(a,4(i4,a))
 996              format(a16,1x,4e30.22)
                  close(lout)
                  FORT_CONPSOLV = 0
!$omp end critical(output)
                  goto 800
               end if
             enddo

#ifdef DO_JBB_HACK_POST

               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
               Tnew(i,j,k) = tspecies(0)

               do m= 1,Nspec
                  Yres(m) = tspecies(m)
               end do

#if defined(DO_JBB_HACK)
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yold(i,j,k,m)+Yres(m)-Ytemp(m)
               end do
#elif  defined(DO_JBB_HACK_POST)
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yres(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yres(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) then
               Ytemp(iN2) = Ytemp(iN2)+one-sum
            endif
#else
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yres(m)
               end do
#endif

#ifdef SOLN_IS_1D
               endif
#endif
            end do
         end do
 800     continue
      end do
!$omp end do
!$omp end parallel

      if (verbose .and. nfails .gt. 0) then
         print*, '*** DVODE failures for last chem block: ', nfails; call flush(6)
      end if
      end
 
#ifdef LMC_SDC
c     used by LMC
      integer function FORT_CONPSOLV_SDC(lo, hi,
     &     rhoYnew,   DIMS(rhoYnew), 
     &     rhoHnew,   DIMS(rhoHnew),
     &     Tnew,      DIMS(Tnew),
     &     rhoYold,   DIMS(rhoYold), 
     &     rhoHold,   DIMS(rhoHold),
     &     Told,      DIMS(Told),
     &     const_src, DIMS(const_src),
     &     FuncCount, DIMS(FuncCount),
     &     dt,
     &     diag, do_diag, do_stiff)
      implicit none

#include "cdwrk.H"
#include "conp.H"
#include "vode.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(rhoYold)
      integer DIMDEC(rhoHold)
      integer DIMDEC(Told)
      integer DIMDEC(rhoYnew)
      integer DIMDEC(rhoHnew)
      integer DIMDEC(Tnew)
      integer DIMDEC(const_src)
      integer DIMDEC(FuncCount)
      integer do_diag, do_stiff
      REAL_T rhoYold(DIMV(rhoYold),*)
      REAL_T rhoHold(DIMV(rhoHold))
      REAL_T Told(DIMV(Told))
      REAL_T rhoYnew(DIMV(rhoYnew),*)
      REAL_T rhoHnew(DIMV(rhoHnew))
      REAL_T Tnew(DIMV(rhoHnew))
      REAL_T const_src(DIMV(const_src),1:Nspec+1)
      REAL_T FuncCount(DIMV(FuncCount))
      REAL_T dt
      REAL_T diag(DIMV(FuncCount),*)

      integer open_vode_failure_file
      external conpFY_sdc, CONPJ_FILE, open_vode_failure_file
      integer i, j, k, m, MF, ISTATE, lout, ITOL
      integer nsub, node, strang_fix, Niter, nfails
      character*(maxspnml) name
      REAL_T RTOL, ATOL(maxspec+2), ATOLEPS, TT1, TT2, RU, RUC, P1atm
      REAL_T Y(maxspec), Z(maxspec+2), ZP(maxspec+2), Yold(maxspec)
      REAL_T dtloc, weight, TT1save, rho, rhoInv
      REAL_T Ct(maxspec),Qt(maxreac), rsum
      REAL_T rhoYtemp(maxspec)
      REAL_T,  parameter :: HtoTerrMAX = BL_REAL_E(7.8,-12)
      integer, parameter :: HtoTiterMAX = 20, IOPT = 1, ITASK = 1
      REAL_T res(0:HtoTiterMAX-1), rhooldInv
      !
      ! Set to .true. if you want to see the 'nfails' output
      !
      logical, parameter :: verbose = .false.

      ITOL    = vode_itol
      RTOL    = vode_rtol
      ATOLEPS = vode_atol
      !
      ! Set molecular weights and pressure in area accessible by conpF
      !
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)

      if (do_stiff .eq. 1) then
         MF = 22  ! Backward difference solver.
      else
         MF = 10  ! Adams integrator
      endif

      if (ITOL.eq.2 .or. ITOL.eq.4) then
         ATOL(Nspec+1) = typVal_RhoH*ATOLEPS
         do m=1,Nspec
            ATOL(m) = ATOLEPS*typVal_Y(m)*typVal_Density
         end do
      else
         ATOL(1) = ATOLEPS
      end if

      if (do_diag.eq.1) then
         nsub  = nchemdiag
         dtloc = dt/nchemdiag
      else
         nsub  = 1
         dtloc = dt
      endif

      nfails = 0

      FORT_CONPSOLV_SDC = 1

!$omp parallel
      !
      ! Force recalculation of jacobian for each XYZ block.
      !
      FIRST = .TRUE.

      !
      ! The following are in /sdc_function/ :
      !
      !   c_0,c_1,rhoH_INIT,T_cell,rhoY_INIT,negative_Y_test
      !
!$omp do private(i,j,k,strang_fix,TT1,ISTATE,Y,node,Ct,Qt,res)
!$omp&private(rhoYtemp,Z,rho,rsum,rhoInv,lout,name,weight,NIter)
!$omp&private(TT1save,TT2,m,rhooldInv) 
#ifdef REGRESSIONTEST
!$omp&schedule(static)
#else
!$omp&schedule(dynamic,1)
#endif
      do k=lo(3),hi(3)
         if (FORT_CONPSOLV_SDC .eq. 0) cycle
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               TT1                = zero
               TT2                = dt
               strang_fix         = 0
               negative_Y_test    = 0
               ISTATE             = 1
               rhoYtemp(1:Nspec)  = rhoYold(i,j,k,1:Nspec)
               Z(1:Nspec)         = rhoYtemp(1:Nspec)
               rhoY_INIT(1:Nspec) = Z(1:Nspec)
               rhoH_INIT          = rhoHold(i,j,k)
               Z(Nspec+1)         = rhoH_INIT
               c_0(1:Nspec+1)     = const_src(i,j,k,1:Nspec+1)
               T_cell             = Told(i,j,k)
               FuncCount(i,j,k)   = 0

#ifdef ALWAYS_NEW_J
            FIRST = .TRUE.
#endif

               if (do_diag.eq.1) then
                  FuncCount(i,j,k) = 0
                  CALL CKYTCP(RWRK(NP),T_cell,Z(1),IWRK(ckbi),RWRK(ckbr),Ct)
                  CALL CKQC(T_cell,Ct,IWRK(ckbi),RWRK(ckbr),Qt)
                  do m=1,Nreac
                     diag(i,j,k,m) = diag(i,j,k,m)+half*dtloc*Qt(m)*million
                  enddo
               endif

               do node = 1,nsub
                  if (node.lt.nsub) then
                     weight = one
                  else
                     weight = half
                  endif

                  TT2     = TT1 + dtloc                  
                  TT1save = TT1

#if !defined(BL_USE_DOUBLE)
                  CALL SVODE
#else
                  CALL DVODE
#endif
     &                 (conpFY_sdc, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL,
     &                 ITASK, ISTATE, IOPT, RWRK(dvbr), dvr, IWRK(dvbi),
     &                 dvi, CONPJ_FILE, MF, RWRK, IWRK)

                  if (ISTATE .LE. -1 .or. negative_Y_test .eq. 1) then

                     strang_fix = 1

                     do m=1,Nspec
                        rhoYold(i,j,k,m) = rhoYold(i,j,k,m) + dt*const_src(i,j,k,m)
                     end do
                     rhoHold(i,j,k) = rhoHold(i,j,k) + dt*const_src(i,j,k,Nspec+1)

                     rho = zero
                     do m=1,Nspec
                        rho = rho+rhoYold(i,j,k,m)
                     end do
                     rhoInv = one / rho
                     rsum = zero
                     do m=1,Nspec
                        rhoYtemp(m) = MAX(rhoYold(i,j,k,m),zero)
                        rsum = rsum + rhoYtemp(m) * rhoInv
                     enddo
                     if (iN2 .gt. 0) then
                        rhoYtemp(iN2) = ( (rhoYtemp(iN2)*rhoInv + one - rsum) * rho )
                     endif

                     TT1                = zero
                     ISTATE             = 1
                     Z(1:Nspec)         = rhoYtemp(1:Nspec)
                     rhoY_INIT(1:Nspec) = Z(1:Nspec)
                     rhoH_INIT          = rhoHold(i,j,k)
                     Z(Nspec+1)         = rhoH_INIT
                     c_0(1:Nspec+1)     = zero
                     T_cell             = Told(i,j,k)
                     FIRST              = .TRUE.
                     FuncCount(i,j,k)   = 0
                     TT1save            = TT1
                     TT2                = TT1 + dtloc

                     if (verbose) then
!$omp atomic
                        nfails = nfails + 1
!$omp end atomic
                     end if

#if !defined(BL_USE_DOUBLE)
                     CALL SVODE
#else
                     CALL DVODE
#endif
     &                    (conpFY_sdc, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL,
     &                    ITASK, ISTATE, IOPT, RWRK(dvbr), dvr, IWRK(dvbi),
     &                    dvi, CONPJ_FILE, MF, RWRK, IWRK)

                     if (ISTATE .LE. -1) then
!$omp critical(output)
                        call conpFY_sdc(NEQ, TT1, Z, ZP, RWRK, IWRK)
                        lout = open_vode_failure_file()
                        write(lout,*)
                        write(lout,995) 'VODE Failed at (i,j,k) = (',i,
     &                       ',',j,',',k,'),   Return code = ',ISTATE
                        write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                        write(lout,995) 
     &                       'State ID,RY,RYp,dRH/dt,dRH/dt*(dt)'
                        name = 'RhoH'
                        write(lout,996) name,
     &                       RhoHold(i,j,k),Z(Nspec+1),
     &                       ZP(Nspec+1),ZP(Nspec+1)*(dt-TT1)
                        do m=1,Nspec
                           call get_spec_name(name,m)
                           write(lout,996) name,RhoYold(i,j,k,m),
     &                          Z(m),ZP(m),ZP(m)*(dt-TT1)
                        end do
                        write(lout,996) 'T_cell:',Z(Nspec+2)
                        rsum = zero
                        do m=1,Nspec
                           rsum = rsum + Z(m)
                        enddo
                        do m=1,Nspec
                           Y(m) = Z(m)/rsum
                        enddo
                        rsum = zero
                        do m=1,Nspec
                           rsum = rsum + RhoYold(i,j,k,m)
                        enddo
                        do m=1,Nspec
                           Yold(m) = RhoYold(i,j,k,m)/rsum
                        enddo
                        
                        write(lout,996)'Mass fracs: old, curr'
                        do m=1,Nspec
                           Y(m) = Z(m)/rsum
                           call get_spec_name(name,m)
                           write(lout,996) name,Yold(m),Y(m)
                        enddo
                        
 995                    format(a,3(i4,a))
 996                    format(a,1x,4e30.22)
                        close(lout)
                        FORT_CONPSOLV_SDC = 0
!$omp end critical(output)
                        goto 800
                     end if
                  end if
                  
                  TT1 = TT2

                  if (do_diag.eq.1) then
                     CALL CKYTCP(RWRK(NP),T_cell,Z(1),IWRK(ckbi),RWRK(ckbr),Ct)
                     CALL CKQC(T_cell,Ct,IWRK(ckbi),RWRK(ckbr),Qt)
                     do m=1,Nreac
                        diag(i,j,k,m) = diag(i,j,k,m)+weight*dtloc*Qt(m)*million
                     enddo
                  endif

                  FuncCount(i,j,k) = FuncCount(i,j,k) + IWRK(dvbi+11)
               enddo
               
               rhoHnew(i,j,k) = Z(Nspec+1)
               
               if (strang_fix .eq. 1) rhooldInv = one / rho

               rho = zero
               do m = 1,Nspec
                  rho = rho + Z(m)
               end do
               rhoInv = 1.d0 / rho

               if (strang_fix .eq. 1) then
                  !
                  ! We want Z to contain rhoYnew(i,j,k)
                  !
                  rhooldInv = one / rho

                  do m= 1,Nspec
                     Z(m) = Z(m) + (rhoYold(i,j,k,m) - rhoYtemp(m))*rhooldInv*rho
                  end do
               end if

               do m=1,Nspec
                  !
                  ! Finally update rhoYnew
                  !
                  rhoYnew(i,j,k,m) = Z(m)
               end do

               do m=1,Nspec
                  Y(m) = rhoYnew(i,j,k,m) * rhoInv
               enddo
               Tnew(i,j,k) = T_cell
               call FORT_TfromHYpt(Tnew(i,j,k),rhoHnew(i,j,k)*rhoInv,Y,HtoTerrMAX,HtoTiterMAX,res,Niter)

            end do
         end do
 800     continue
      end do
!$omp end do
!$omp end parallel

      if (verbose .and. nfails .gt. 0) then
         print*, '*** DVODE failures for last chem block: ', nfails; call flush(6)
      end if
      end
#endif


c     used by LMC
      subroutine FORT_BETA_WBAR(lo, hi, RD, DIMS(RD), RD_Wbar, DIMS(RD_Wbar), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RD)
      integer DIMDEC(RD_Wbar)
      integer DIMDEC(Y)
      REAL_T RD(DIMV(RD),*)
      REAL_T RD_Wbar(DIMV(RD_Wbar),*)
      REAL_T Y(DIMV(Y),*)

      integer i, j, k, n
      REAL_T Yt(maxspec), RHO, Wavg

      if (use_eg .eq. 1) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)

#ifdef LMC_SDC
                  RHO = 0.d0
                  do n=1,Nspec
                     RHO = RHO + Y(i,j,k,n)
                  end do

                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n) / RHO
                  end do
#else
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
#endif
                  CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)

                  do n=1,Nspec
                     RD_Wbar(i,j,k,n) = RD(i,j,k,n) * Yt(n) / Wavg
                  end do

               end do
            end do
         end do

      else if (use_mc .eq. 1) then

         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)

#ifdef LMC_SDC
                  RHO = 0.d0
                  do n=1,Nspec
                     RHO = RHO + Y(i,j,k,n)
                  end do

                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n) / RHO
                  end do
#else
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
#endif

                  CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
                  do n=1,Nspec
                     RD_Wbar(i,j,k,n) = RD(i,j,k,n) * Yt(n) / Wavg
                  end do

               end do
            end do
         end do

      endif

      end

c     used by LMC
      subroutine FORT_MIXAVG_RHODIFF_TEMP(lo, hi, RD, DIMS(RD), T,
     &     DIMS(T), Y, DIMS(Y), Patm, do_temp, do_VelVisc)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM), do_temp, do_VelVisc
      integer DIMDEC(RD)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RD(DIMV(RD),*)
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm

      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), Dt(maxspec)
      REAL_T SCAL, TSCAL, Wavg, RHO, Tt, invmwt(maxspec)
      REAL_T alpha, l1, l2, X(maxspec), CPMS(maxspec)

      parameter(SCAL = tenth, TSCAL = one / 100000.0D0)

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = one / invmwt(n)
      end do

      if (use_eg .eq. 1) then
!$omp parallel do private(i,j,k,n,Yt,Tt,alpha,Dt,Wavg)
!$omp&private(CPMS,X,RHO,l1,l2)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)

#ifdef LMC_SDC
                  RHO = 0.d0
                  do n=1,Nspec
                     RHO = RHO + Y(i,j,k,n)
                  end do

                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n) / RHO
                  end do
#else
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
#endif
                  Tt = MAX(T(i,j,k),TMIN_TRANS)
                  CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
                  CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL EGSPAR(Tt,X,Yt,CPMS,EGRWRK(1),EGIWRK(1))
                  CALL EGSV1(Ptmp,Tt,Yt,Wavg,EGRWRK(1),Dt)
#ifdef LMC_SDC
                  RHO = RHO*1.d-3
#else
                  CALL CKRHOY(Ptmp,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
#endif
                  do n=1,Nspec
                     RD(i,j,k,n) = RHO * Wavg * invmwt(n) * Dt(n) * SCAL
                  end do

                  if (do_temp .ne. 0) then
                     alpha = 1
                     CALL EGSL1(alpha,Tt,X,EGRWRK(1),l1)
                     alpha = -1
                     CALL EGSL1(alpha,Tt,X,EGRWRK(1),l2)
                     RD(i,j,k,Nspec+1) = half * (l1 + l2) * TSCAL
                  endif

                  if (do_VelVisc .ne. 0) then
                     CALL EGSE3(Tt,Yt,EGRWRK(1),RD(i,j,k,Nspec+2))
                     RD(i,j,k,Nspec+2) = RD(i,j,k,Nspec+2) * SCAL
                  endif

               end do
            end do
         end do
!$omp end parallel do
      else if (use_mc .eq. 1) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)

#ifdef LMC_SDC
                  RHO = 0.d0
                  do n=1,Nspec
                     RHO = RHO + Y(i,j,k,n)
                  end do

                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n) / RHO
                  end do
#else
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
#endif

                  Tt = MAX(T(i,j,k),TMIN_TRANS) 
                  CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL MCADIF(Ptmp,Tt,X,MCRWRK,Dt)
#ifdef LMC_SDC
                  RHO = RHO*1.d-3
#else
                  CALL CKRHOY(Ptmp,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
#endif
                  
                  do n=1,Nspec
                     RD(i,j,k,n) = RHO * Dt(n) * SCAL
                  end do
                  
                  if (do_temp .ne. 0) then
                     CALL MCACON(Tt,X,MCRWRK,l1)
                     RD(i,j,k,Nspec+1) = l1 * TSCAL
                  endif
                  
                  if (do_VelVisc .ne. 0) then
                     CALL MCAVIS(Tt,X,MCRWRK,RD(i,j,k,Nspec+2))
                     RD(i,j,k,Nspec+2) = RD(i,j,k,Nspec+2) * SCAL
                  endif
                  
               end do
            end do
         end do
      endif

      end

c     used LMC SDC
      subroutine FORT_MIX_SHEAR_VISC(lo, hi, eta, DIMS(eta),
     &                               T, DIMS(T), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(eta)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T eta(DIMV(eta))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T SCAL
      
      integer i, j, k, n
      REAL_T X(maxspec), Yt(maxspec), CPMS(maxspec), Tt

      parameter(SCAL = tenth)
      !
      !The following computes the mixture averaged shear viscosity 
      ! using EGLib. Note that SCAL converts assumed cgs units to
      ! MKS (1 g/cm.s = .1 kg/m.s)
      !
      if (use_eg .eq. 1) then
!$omp parallel do private(i,j,k,n,Yt,Tt,CPMS,X)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
                  Tt = MAX(T(i,j,k),TMIN_TRANS) 
                  CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL EGSPAR(Tt,X,Yt,CPMS,EGRWRK(1),EGIWRK(1))
                  CALL EGSE3(Tt,Yt,EGRWRK(1),eta(i,j,k))
                  eta(i,j,k) = eta(i,j,k) * SCAL
               end do
            end do
         end do
!$omp end parallel do
      else if (use_mc .eq. 1) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  do n=1,Nspec
                     Yt(n) = Y(i,j,k,n)
                  end do
                  Tt = MAX(T(i,j,k),TMIN_TRANS) 
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL MCAVIS(Tt,X,MCRWRK,eta(i,j,k))
                  eta(i,j,k) = eta(i,j,k) * SCAL
               end do
            end do
         end do
      endif

      end

c     used by LMC
      subroutine FORT_RHOfromPTY(lo, hi, RHO, DIMS(RHO), T, DIMS(T),
     &                           Y, DIMS(Y), Patm)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
c
      parameter(SCAL = one * 1000.0D0)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_RHOfromPvTY(lo, hi, RHO, DIMS(RHO), T, DIMS(T),
     &                           Y, DIMS(Y), P, DIMS(P))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      integer DIMDEC(P)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T P(DIMV(P))
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
c
      parameter(SCAL = one * 1000)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Ptmp = P(i,j,k) * P1ATM
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
      end
      
c     used LMC SDC
      subroutine FORT_PfromRTY(lo, hi, P, DIMS(P), RHO, DIMS(RHO),
     &                         T, DIMS(T), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(P)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T P(DIMV(P))
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), RHOt, SCAL, SCAL1
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
c           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
c
      parameter(SCAL = tenth, SCAL1 = tenth**3)

!$omp parallel do private(i,j,k,n,Yt,RHOt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               RHOt = RHO(i,j,k) * SCAL1
               CALL CKPY(RHOt,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),P(i,j,k))
               P(i,j,k) = P(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_TfromPRY(lo, hi, T, DIMS(T), RHO, DIMS(RHO),
     &                         Y, DIMS(Y), Patm)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(RHO)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T RHO(DIMV(RHO))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      REAL_T Patm
      
      integer i, j, k, n
      REAL_T RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL, Wavg, RHOt
c
c     NOTE: SCAL converts density (1 kg/m^3 = 1.e-3 g/cm^3)
c
      parameter(SCAL = tenth**3)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               RHOt = RHO(i,j,k) * SCAL
               T(i,j,k) = Ptmp / (RHOt * RU / Wavg)
            end do
         end do
      end do
      end
      
c     used by LMC
      subroutine FORT_CPMIXfromTY(lo, hi, CPMIX, DIMS(CPMIX), T, DIMS(T),
     &                            Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CPMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T CPMIX(DIMV(CPMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCPBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX(i,j,k))
               CPMIX(i,j,k) = CPMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do
      end
      
      subroutine FORT_CVMIXfromTY(lo, hi, CVMIX, DIMS(CVMIX), T, DIMS(T),
     &                            Y, DIMS(Y))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CVMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T CVMIX(DIMV(CVMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCVBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CVMIX(i,j,k))
               CVMIX(i,j,k) = CVMIX(i,j,k) * SCAL
            end do
         end do
      end do
      end
      
c     used by LMC
      subroutine FORT_HMIXfromTY(lo, hi, HMIX, DIMS(HMIX), T, DIMS(T),
     &                           Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(HMIX)
      integer DIMDEC(T)
      integer DIMDEC(Y)
      REAL_T HMIX(DIMV(HMIX))
      REAL_T T(DIMV(T))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k, n
      REAL_T Yt(maxspec), SCAL
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKHBMS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),HMIX(i,j,k))
               HMIX(i,j,k) = HMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end
      
c     used by LMC
      subroutine FORT_MWMIXfromY(lo, hi, MWMIX, DIMS(MWMIX), Y, DIMS(Y))
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(MWMIX)
      integer DIMDEC(Y)
      REAL_T MWMIX(DIMV(MWMIX))
      REAL_T Y(DIMV(Y),*)
      
      integer i, j, k,n
      REAL_T Yt(maxspec)
!
!     Returns mean molecular weight in kg/kmole
!
!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),MWMIX(i,j,k))
            end do
         end do
      end do
!$omp end parallel do

      end
      
      subroutine FORT_CPfromT(lo, hi, CP, DIMS(CP), T, DIMS(T))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(CP)
      integer DIMDEC(T)
      REAL_T CP(DIMV(CP),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, CPt(maxspec)
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
c
      parameter(SCAL = tenth**4)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKCPMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),CPt)
               do n=1,Nspec
                  CP(i,j,k,n) = CPt(n) * SCAL
               end do
            end do
         end do
      end do
      end
      
c     used by LMC
      subroutine FORT_HfromT(lo, hi, H, DIMS(H), T, DIMS(T))
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(H)
      integer DIMDEC(T)
      REAL_T H(DIMV(H),*)
      REAL_T T(DIMV(T))
      
      integer i, j, k, n
      REAL_T SCAL, Ht(maxspec)
c
c     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
c
      parameter(SCAL = tenth**4)

!$omp parallel do private(i,j,k,n,Ht)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),Ht)
               do n=1,Nspec
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do
!$omp end parallel do

      end

c     used by LMC
      integer function FORT_TfromHY(lo, hi, T, DIMS(T),
     &                              HMIX, DIMS(HMIX), Y, DIMS(Y),
     &                              errMax, NiterMAX, res)
      implicit none

#include "cdwrk.H"

      integer lo(SDIM), hi(SDIM)
      integer NiterMAX
      integer DIMDEC(T)
      integer DIMDEC(HMIX)
      integer DIMDEC(Y)
      REAL_T T(DIMV(T))
      REAL_T HMIX(DIMV(HMIX))
      REAL_T Y(DIMV(Y),*)
      REAL_T errMAX
      REAL_T res(0:NiterMAX-1)
      REAL_T Yt(maxspec), lres(0:NiterMAX-1)
      integer i,j,k,n,Niter,MAXiters

      MAXiters = 0

!$omp parallel do private(i,j,k,n,Yt,lres,Niter) reduction(max:MAXiters)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do

               call FORT_TfromHYpt(T(i,j,k),HMIX(i,j,k),Yt,errMax,NiterMAX,lres,Niter)

               if (Niter .lt. 0) then
                  write(6,*) 'T from h,y solve in FORT_TfromHY failed',Niter
                  call bl_abort(" ")
               end if
               
               if (Niter .gt. MAXiters) MAXiters = Niter

            end do
         end do
      end do
!$omp end parallel do

c     Set max iters taken during this solve, and exit
      FORT_TfromHY = MAXiters
      return
      end
c
c     Optically thin radiation model, specified at
c            http://www.ca.sandia.gov/tdf/Workshop/Submodels.html
c     
c     Q(T,species) = 4*sigma*SUM{pi*aP,i} *(T4-Tb4) 
c     
c     sigma=5.669e-08 W/m2K4 is the Steffan-Boltzmann constant, 
c     SUM{ } represents a summation over the species in the radiation calculation, 
c     pi is partial pressure of species i in atm (Xi times local pressure)
c     aP,i is the Planck mean absorption coefficient of species i, 1/[m.atm]
c     T is the local flame temperature (K)
c     Tb is the background temperature (300K or as spec. in expt)
c
c     For H2O and CO2,
c         aP = exp{c0 + c1*ln(T) + c2*{ln(T)}2 + c3*{ln(T)}3 + c4*{ln(T)}4} 
c     
c                            H2O                  CO2
c              c0       0.278713E+03         0.96986E+03
c              c1      -0.153240E+03        -0.58838E+03
c              c2       0.321971E+02         0.13289E+03
c              c3      -0.300870E+01        -0.13182E+02
c              c4       0.104055E+00         0.48396E+00
c     For CH4:
c     
c     aP,ch4 = 6.6334 - 0.0035686*T + 1.6682e-08*T2 + 2.5611e-10*T3 - 2.6558e-14*T4
c
c     For CO:   aP,co = c0+T*(c1 + T*(c2 + T*(c3 + T*c4)))
c
c           T <= 750                 else
c      
c         c0   4.7869              10.09       
c         c1  -0.06953             -0.01183    
c         c2   2.95775e-4          4.7753e-6   
c         c3  -4.25732e-7          -5.87209e-10
c         c4   2.02894e-10         -2.5334e-14 
c      
      
c     used by LMC
      subroutine FORT_OTrad_TDF(lo, hi, Qloss, DIMS(Qloss),
     &                          T, DIMS(T), X, DIMS(X), Patm, T_bg) 
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM)
      integer DIMDEC(Qloss)
      integer DIMDEC(T)
      integer DIMDEC(X)
      REAL_T Qloss(DIMV(Qloss))
      REAL_T T(DIMV(T))
      REAL_T X(DIMV(X),*)
      REAL_T Patm, T_bg
      
      character*(maxspnml) name      
      integer n, i, j, k, iH2O, iCO2, iCH4, iCO
      REAL_T aP, c0, c1, c2, c3, c4
      REAL_T T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4,Tb4,sigma

      parameter (sigma = 5.669D-08)
      
      iH2O = 0; iCO2 = 0; iCH4 = 0; iCO  = 0

      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .EQ. 'H20') iH2O = n
         if (name .EQ. 'CO2') iCO2 = n
         if (name .EQ. 'CH4') iCH4 = n
         if (name .EQ. 'CO')  iCO  = n
      end do
      
      Tb4 = T_bg**4

!$omp parallel do private(i,j,k,T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4)
!$omp&private(aP,c0,c1,c2,c3,c4)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               T1 = T(i,j,k)
               T2 = T1*T1
               T3 = T2*T1
               T4 = T3*T1

               lnT1 = zero; lnT2 = zero; lnT3 = zero; lnT4 = zero
               
               if ( (iH2O.gt.0) .or. (iCO2.gt.0) ) then
                  lnT1 = LOG(T1)
                  lnT2 = lnT1*lnT1
                  lnT3 = lnT2*lnT1
                  lnT4 = lnT3*lnT1
               end if
               
               aP = zero
               
               if ((iH2O.gt.0).and.(X(i,j,k,iH2O).gt.zero)) then
                  aP = aP + X(i,j,k,iH2O)*EXP(
     &                 + 0.278713D+03
     &                 - 0.153240D+03*lnT1
     &                 + 0.321971D+02*lnT2
     &                 - 0.300870D+01*lnT3
     &                 + 0.104055D+00*lnT4 )
               end if
               
               if ((iCO2.gt.0).and.(X(i,j,k,iCO2).gt.zero)) then            
                  aP = aP + X(i,j,k,iCO2)*EXP(
     &                 + 0.96986D+03
     &                 - 0.58838D+03*lnT1
     &                 + 0.13289D+03*lnT2
     &                 - 0.13182D+02*lnT3
     &                 + 0.48396D+00*lnT4 )
               end if

               if ((iCH4.gt.0).and.(X(i,j,k,iCH4).gt.zero)) then
                  aP = aP + X(i,j,k,iCH4)*
     &                 ( 6.6334D0
     &                 - 0.0035686D0 *T1
     &                 + 1.6682D-08*T2
     &                 + 2.5611D-10*T3
     &                 - 2.6558D-14*T4 )         
               end if
               
               if ((iCO.gt.0).and.(X(i,j,k,iCO).gt.zero)) then
                  if ( T1 .le. 750.0D0 ) then
                     c0 =  4.7869D0
                     c1 = -0.06953D0
                     c2 =  2.95775D-4
                     c3 = -4.25732D-7
                     c4 =  2.02894D-10
                  else
                     c0 =  10.09D0
                     c1 = -0.01183D0
                     c2 =  4.7753D-6
                     c3 = -5.87209D-10
                     c4 = -2.5334D-14
                  endif
                  aP = aP + X(i,j,k,iCO)*(c0 + c1*T1 + c2*T2 + c3*T3 + c4*T4)
               end if

               Qloss(i,j,k) = four*sigma*Patm*(T4-Tb4)*aP
               
            end do
         end do
      end do
!$omp end parallel do

      end
