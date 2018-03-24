#include "ChemDriver_F.H"
#include "AMReX_CONSTANTS.H"

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define three4th    0.75
#define onepoint27  1.27
#define point4      0.4
#define point67     0.67
#define point14     0.14
#define onetenthsnd 0.0001
#define ten2minus18 1.0e-18
#else
#define three4th    0.75d0
#define onepoint27  1.27d0
#define point4      0.4d0
#define point67     0.67d0
#define point14     0.14d0
#define onetenthsnd 0.0001d0
#define ten2minus18 1.0d-18
#endif

module chem_driver

  implicit none

  private
  
  public :: drv_get_reaction_map, set_Tmin_trans, SETVERBOSEVODE, &
            SETVODETOLS, set_vode_subcyc, set_spec_scal_Y, INITCHEM, &
            finalize_chem, GETCKMAXNAMELEN, GETCKDIMPARAMS, &
            FINDLHS, FINDRHS, SETNU, CKINU, CKELTXINSPY, GETCKNUMSPEC, &
            GETCKNUMELT, get_CK_num_reac, RUNIV, P1ATMMKS, GETCKELTNAME, &
            GETCKSPECNAME,CKSYMR,get_spec_name,get_spec_number,get_CKMWT, &
            FORT_GETCKAWT,conpFY,conpJY,conpFY_sdc,TfromeYpt,TfromHYpt, &
            open_vode_failure_file

contains

  subroutine drv_get_reaction_map(rmap)bind(C,name="drv_get_reaction_map")
  
    implicit none

#include "cdwrk.H"

    integer rmap(nReac)
    call GET_REACTION_MAP(rmap) 
  
  end subroutine drv_get_reaction_map

!------------------------------------
  
  subroutine set_Tmin_trans(TminTRANS)bind(C,name="set_Tmin_trans")
  
    implicit none
    REAL_T TminTRANS

#include "cdwrk.H"

    TMIN_TRANS = TminTRANS

  end subroutine set_Tmin_trans

!------------------------------------  
  
  subroutine SETVERBOSEVODE()bind(C, name="SETVERBOSEVODE")
    
    implicit none

#include "cdwrk.H"

    verbose_vode = 1

  end subroutine SETVERBOSEVODE

!------------------------------------
    
  subroutine SETVODETOLS(rtol,atol,itol)bind(C, name="SETVODETOLS")
  
    implicit none
    integer itol
    REAL_T rtol,atol

#include "cdwrk.H"

    vode_itol = itol
    vode_rtol = rtol
    vode_atol = atol

  end subroutine SETVODETOLS

!------------------------------------  
  
  subroutine set_vode_subcyc(maxcyc)bind(C, name="set_vode_subcyc")

    implicit none
    integer maxcyc

#include "cdwrk.H"

    max_vode_subcycles = maxcyc

  end subroutine set_vode_subcyc

!------------------------------------  
  
  subroutine set_spec_scal_Y(name, nlength)bind(C, name="set_spec_scal_Y")
    
    implicit none

#include "cdwrk.H"

    integer nlength, name(nlength), i, j, maxlen
    REAL_T val
    parameter (maxlen=256)
    character filet*(maxlen)
    character*(maxspnml) spname, spinname
!      
!     Convert encoded names to strings, and open file
!
    if (nlength.GT.maxlen) then
      call bl_abort('FORT_SETSPECSCAL: scale file name too long')
    end if
      
    do i = 1, nlength
      filet(i:i) = char(name(i))
    end do
    open(unit=51,status='OLD',form='FORMATTED', &
          file=filet(1:nlength),err=30)
      
 10 continue 
    read(51,*,end=20) spinname, val
    do j = 1,Nspec
      call get_spec_name(spname,j)
      if (spname .eq. spinname) then
        spec_scalY(j) = ABS(val)
      end if
    end do
    goto 10
 20 close(51)
    goto 40
 30 write(6,*) 'Trouble opening file = ',filet(1:nlength)
    call bl_abort(" ")
 40 continue
 
  end subroutine set_spec_scal_Y

!------------------------------------ 
  
  subroutine INITCHEM()bind(C, name="INITCHEM")
  
    implicit none

#include "cdwrk.H"
#include "conp.H"

    interface
      integer function FORT_USINGEG() bind(C,name="FORT_USINGEG")
        use iso_c_binding
      end function
         
      integer function FORT_USINGMC() bind(C,name="FORT_USINGMC")
        use iso_c_binding
       end function
    end interface
      
    integer n, RTOT, lout, egrlen, egilen, idummy(1)
    double precision rdummy(1)
    character*(maxspnml) name
    integer MAXFIT, NO, NFDIM, NT, NRANGE, NLITEMAX
    integer mcr, mci
    integer ierr, LOUTCK
    !
    ! Set a few default values
    !
    verbose_vode       = 0
    max_vode_subcycles = 15000
    spec_scalY         = one
    thickFacCH         = one
    !
    ! Get chemistry mechanism parameters.
    !
    CALL CKINIT()
    CALL CKINDX(idummy(1),rdummy(1),Nelt,Nspec,Nreac,Nfit)
    !
    ! Set up EGlib workspace.
    !
    ! When OPENMP give each thread its own space.
    !
    egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes &
            + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
    egilen = Nspec
    !
    ! Ditto for DVODE workspace.
    !
    NEQ   = Nspec + 1
    dvr   = 22 + 9*NEQ + 2*NEQ**2
    dvi   = 30 + NEQ
    dvdr  = Nspec + NEQ*2 + 2
    dvdbr = dvbr+dvr
    dvder = dvdbr+dvdr-1
    !      
    ! Set pointers into conp common blocks.
    !
    NP    = dvdbr
    NRHO  = NP  + 1
    NWT   = NRHO  + 1
    NWTI  = NWT + Nspec
    NZ    = NWTI + Nspec
    RTOT  = NZ  + NEQ - 1
    !
    ! Setup tranlib space
    !
    MAXFIT=7
    NO=4
    NFDIM=165
    NT=50
    NRANGE = MAXTP-1
    NLITEMAX=3
    mcr = Nspec*(19+2*NO+NO*NLITEMAX)+(NO+15)*Nspec**2
    mci = 4*Nspec + NLITEMAX
    LLINKMC = 44
    LOUTCK = 6

!$omp parallel
    allocate(RWRK(dvr+dvdr))
    allocate(IWRK(dvi))

!$omp single
    use_eg = 0
    use_mc = 0
!$omp end single
    if (FORT_USINGEG() .eq. 1) then
      allocate(EGRWRK(egrlen))
      allocate(EGIWRK(egilen))         
      call EGINICD(eg_nodes, lout, eg_IFLAG, eg_ITLS, &
                   EGRWRK, egrlen, EGIWRK, egilen)
!$omp single
      use_eg = 1
!$omp end single
    else if (FORT_USINGMC() .eq. 1) then
      allocate(MCRWRK(mcr))
      allocate(MCIWRK(mci))
      CALL MCINITCD (LOUTCK, mci, mcr, MCIWRK, MCRWRK, ierr)
      if (ierr .gt. 0) then
        WRITE(LOUTCK,*)' QUITTING BECAUSE MCINIT IFLAG = ', ierr
        call bl_abort(" ")
      end if
!$omp single
      use_mc = 1
!$omp end single
    else
      write(6,*) 'Unknown transport library'
      call bl_abort(" ")
    endif
    !
    !     Set IOPT=1 parameter settings for VODE
    !     They only really need to be set once.
    !     
    RWRK(dvbr+4) = 0
    RWRK(dvbr+5) = 0
    RWRK(dvbr+6) = 1.d-19
    IWRK(dvbi+4) = 0
    IWRK(dvbi+5) = max_vode_subcycles
    IWRK(dvbi+6) = 0
    !
    ! Set molecular weights where conpF can access'm.
    !
    CALL CKWT(IWRK(ckbi), RWRK(ckbr), RWRK(NWT))
    do n = 0,Nspec-1
      RWRK(NWTI+n) = 1.d0 / RWRK(NWT+n)
    enddo
!$omp end parallel

    if (RTOT .GT. dvder) then
      write(6,*) 'Memory layout bust, dvdr not big enough'
      write(6,*) RTOT, dvder
      call bl_abort(" ")
    end if
    !
    ! Find N2 in the list.
    !
    iN2 = -1
    do n = 1,Nspec
      call get_spec_name(name,n)
      if (name .eq. 'N2' ) iN2 = n
    end do
    if (iN2.eq.-1) &
        write(6,*) '.....warning: no N2 in chemistry species list'
  end subroutine INITCHEM

!----------------------------------
  
  subroutine finalize_chem()bind(C, name="finalize_chem")
 
    implicit none

#include "cdwrk.H"

    interface
      integer function FORT_USINGEG() bind(C,name="FORT_USINGEG")
        use iso_c_binding
      end function
         
      integer function FORT_USINGMC() bind(C,name="FORT_USINGMC")
        use iso_c_binding
      end function
    end interface

!$omp parallel
    
    if (FORT_USINGEG() .eq. 1) then
      deallocate(EGRWRK)
      deallocate(EGIWRK)
    else if (FORT_USINGMC() .eq. 1) then
      deallocate(MCRWRK)
      deallocate(MCIWRK)
    else
      write(6,*) 'Unknown transport library'
      call bl_abort(" ")
    endif
    deallocate(RWRK)
    deallocate(IWRK)
!$omp end parallel

    call CKFINALIZE();

  end subroutine finalize_chem

!------------------------------------  
  
  integer function GETCKMAXNAMELEN()bind(C, name="GETCKMAXNAMELEN")
  
    implicit none

#include "cdwrk.H"
    GETCKMAXNAMELEN = maxspnml

  end function GETCKMAXNAMELEN

!------------------------------------  
  
  subroutine GETCKDIMPARAMS(imaxreac, imaxspec, imaxelts, &
                            imaxord, imaxthrdb, imaxtp, imaxsp, &
                            imaxspnml)&
                            bind(C, name="GETCKDIMPARAMS")
                                    
                                    
    implicit none
    integer imaxreac, imaxspec, imaxelts, imaxord
    integer imaxthrdb, imaxtp, imaxsp, imaxspnml

#include "cdwrk.H"

    imaxreac = maxreac
    imaxspec = maxspec
    imaxelts = maxelts
    imaxord = 10
    imaxthrdb = maxthrdb
    imaxtp = maxtp
    imaxsp = maxsp
    imaxspnml = maxspnml
      
  end subroutine GETCKDIMPARAMS

!------------------------------------  
    
  subroutine FINDLHS(reactions, Nreacs, id)bind(C, name="FINDLHS")
  
    implicit none
    integer reactions(*), Nreacs, id

#ifdef MIKE
#include "cdwrk.H"

    integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
    Ndim = maxsp
    if ((id.le.0).or.(id.gt.Nspec)) then
      write(6,*) 'FINDLHS:  species id out of range: ',id
      call bl_abort(" ")
    end if
    Nreacs = 0
    do j=1,Nreac
      CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
      do n=1,Nids
        if ((KI(n).eq.id).and.(NU(n).lt.0)) then
               Nreacs = Nreacs + 1
               reactions(Nreacs) = j
        endif
      end do
    end do
#else
    call bl_abort("FINDLHS not implemented")
#endif
  end subroutine FINDLHS

!------------------------------------  
    
  subroutine FINDRHS(reactions, Nreacs, id)bind(C,name="FINDRHS")
  
    implicit none
    integer reactions(*), Nreacs, id

#ifdef MIKE
#include "cdwrk.H"
    
    integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
    Ndim = maxsp
    if ((id.le.0).or.(id.gt.Nspec)) then
      write(6,*) 'FINDRHS:  species id out of range: ',id
      call bl_abort(" ")
    end if
    Nreacs = 0
    do j=1,Nreac
      CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
      do n=1,Nids
        if ((KI(n).eq.id).and.(NU(n).gt.0)) then
               Nreacs = Nreacs + 1
               reactions(Nreacs) = j
        endif
      end do
    end do
#else
    call bl_abort("FINDRHS not implemented")
#endif
  end subroutine FINDRHS

!------------------------------------  
  
  subroutine SETNU(nu,lenNU)bind(C, name="SETNU")
  
    implicit none

#include "cdwrk.H"

    integer lenNU
    integer nu(maxreac,maxspec)

    if (lenNU .lt. maxreac*maxspec) then
      write(6,*) 'FORT_CKNU:  nu work array too small: '
      call bl_abort(" ")
    endif
    call CKNU(maxreac, IWRK(ckbi), RWRK(ckbr), nu)
 
  end subroutine SETNU

!------------------------------------  
  
  subroutine CKINU(Nids,KI,lenKI,NU,lenNU,rxnID,nuAll)bind(C, name="CKINU")
  
    implicit none

#include "cdwrk.H"

    integer lenKI,lenNU
    integer rxnID, Nids, KI(lenKI), NU(lenNU), nuAll(maxreac,maxspec)
    integer Ndim, k
    Ndim = MIN(lenKI,lenNU)
    if ((rxnID.le.0).or.(rxnID.gt.Nreac)) then
      write(6,*) 'CKINU:  reaction id out of range: ',rxnID
      call bl_abort(" ")
    end if
    if (Ndim.lt.maxsp) then
      call bl_abort('CKINU:  KI or NU not long enough')
    end if
    Nids = 0
    do k=1,Nspec
      if (nuAll(rxnID,k).ne.0) then
            Nids = Nids + 1
            KI(Nids) = k
            NU(Nids) = nuAll(rxnID,k)
      endif
    enddo

  end subroutine CKINU

!------------------------------------  
    
  integer function CKELTXINSPY(eltID, spID)bind(C, name="CKELTXINSPY")
  
    implicit none
      
#include "cdwrk.H"

    integer eltID, spID
    integer NCF(maxelts,maxspec)
    CALL CKNCF(maxelts, IWRK(ckbi), RWRK(ckbr), NCF)
    CKELTXINSPY = NCF(eltID+1,spID+1)
      
  end function CKELTXINSPY

!------------------------------------  
  
  integer function GETCKNUMSPEC()bind(C,name="GETCKNUMSPEC")
  
    implicit none
      
#include "cdwrk.H"

    GETCKNUMSPEC = Nspec
      
  end function GETCKNUMSPEC

!------------------------------------  
    
  integer function GETCKNUMELT()bind(C, name="GETCKNUMELT")
  
    implicit none
      
#include "cdwrk.H"

    GETCKNUMELT = Nelt
      
  end function GETCKNUMELT

!------------------------------------  
  
  integer function get_CK_num_reac()bind(C, name="get_CK_num_reac")
  
    implicit none
      
#include "cdwrk.H"

    get_CK_num_reac = Nreac
      
  end function get_CK_num_reac

!------------------------------------  
  
  double precision function RUNIV()bind(C, name="RUNIV")

    implicit none
    double precision Ruc, Pa

#include "cdwrk.H"

    call CKRP(IWRK(ckbi),RWRK(ckbr),RUNIV,Ruc,Pa)
!     1 erg/(mole.K) = 1.e-4 J/(kmole.K)
    RUNIV = RUNIV*1.d-4

  end function RUNIV

!------------------------------------  
  
  double precision function P1ATMMKS()bind(C, name="P1ATMMKS")

    implicit none
    double precision Ru, Ruc, Pa

#include "cdwrk.H"

    call CKRP(IWRK(ckbi),RWRK(ckbr),Ru,Ruc,Pa)
!     1 N/(m.m) = 0.1 dyne/(cm.cm)
    P1ATMMKS = Pa*1.d-1

  end function P1ATMMKS

!------------------------------------  
  
  integer function GETCKELTNAME(i, coded)bind(C, name="GETCKELTNAME")
  
    implicit none

#include "cdwrk.H"

    integer i
    integer coded(*)
    integer names(0:maxelts*2)
    integer nlen
    nlen = 2
    call CKSYME(names,nlen)
    coded(1) = names(2*(i-1)  )
    coded(2) = names(2*(i-1)+1)
    if (coded(2).eq.ICHAR(' ')) then
      GETCKELTNAME = 1
    else
      GETCKELTNAME = 2
    endif

  end function GETCKELTNAME

!------------------------------------  
  
  integer function GETCKSPECNAME(i, coded)bind(C, name="GETCKSPECNAME")
  
    implicit none
      
#include "cdwrk.H"

    integer i
    integer coded(*)
    integer names(maxspec*maxspnml)
    integer j, str_len
    str_len = 0
    call CKSYMS(names, maxspnml)
    do j = 1, maxspnml
      coded(j) = names(maxspnml*(i-1)+j)
    end do
    do j = 1, maxspnml
      if (coded(j).eq.ICHAR(' ')) then
        str_len = j
        exit
      endif 
    end do
    GETCKSPECNAME = str_len - 1
 
  end function GETCKSPECNAME

 !------------------------------------ 
  
  integer function CKSYMR(fortReacIdx, coded)bind(C, name="CKSYMR")
  
    implicit none
    integer fortReacIdx
    integer coded(*)

#ifdef MIKE
#include "cdwrk.H"

    character*(72) line 
    integer j, str_len, istr, iend, lout
    logical error

    error = .false.
    lout = 6
    call CKSYMR(fortReacIdx,lout,IWRK(ckbi),RWRK(ckbr), &
                CWRK(ckbc),str_len,line,error)
    if (error) then
      write(lout,*) 'Could not get reaction name for ',fortReacIdx
      call bl_abort(" ")
    end if
!      
!     Encode the name for transfer to C++
!
    istr = 1
    do while (line(istr:istr) .EQ. ' ')
      istr = istr + 1
    end do
    do j = 0, str_len-1
      coded(j+1) = ICHAR(line(istr+j:istr+j))
    end do
    CKSYMR = str_len
#else
    CKSYMR = 0
    call bl_abort("CKSYMR not implemented")
#endif
  end function CKSYMR

!------------------------------------  
    
  subroutine get_spec_name(name, j)
  
    implicit none
      
#include "cdwrk.H"

    integer i, j!, GETCKSPECNAME
    integer coded(maxspnml), len
    character*(maxspnml) name
    len = GETCKSPECNAME(j, coded)
    do i = 1, maxspnml
      name(i:i) = ' '
    end do
    do i = 1, len
      name(i:i) = char(coded(i))
    end do
    
  end subroutine get_spec_name

!------------------------------------  
  
  subroutine get_spec_number(name, j)
  
    implicit none

#include "cdwrk.H"

    integer j, n
    character*(*) name
    character*(maxspnml) locName
      
    j = -1
    do n = 1, Nspec
      call get_spec_name(locName, n)
      if (locName .EQ. name) j = n
    end do

  end subroutine get_spec_number

!------------------------------------  
  
  subroutine get_CKMWT(mwt)bind(C, name="get_CKMWT")
  
    implicit none
      
#include "cdwrk.H"

    REAL_T mwt(*)
!     Result in kg/kmole
    call CKWT(IWRK(ckbi),RWRK(ckbr),mwt)

  end subroutine get_CKMWT

!------------------------------------

  subroutine FORT_GETCKAWT(awt)

    implicit none
#include "cdwrk.H"

    REAL_T awt(*)
!     Result in kg/kmole
    call CKAWT(IWRK(ckbi),RWRK(ckbr),awt)

  end subroutine FORT_GETCKAWT

!------------------------------------  
  
  subroutine conpFY(N, TIME, Z, ZP, RPAR, IPAR)
  
    implicit none
      
#include "cdwrk.H"
#include "conp.H"

    REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
    integer N, IPAR(*)
      
    REAL_T RHO, CPB, SUM, H, WDOT, WT, THFAC
    integer K

    REAL_T CONC(maxspec), WDOTS(maxspec), ENTHALPY(maxspec)
!
!     Variables in Z are:  Z(1)   = T
!                          Z(K+1) = Y(K)
      
    CALL CKRHOY(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),RHO)
    CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)
    CALL CKYTCP(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CONC)

!
!     Get net production rates.  Compute such that production from -ve
!     reactants gives zero contrib.
!      
!      do k=1,Nspec
!         if (CONC(k) .lt. zero) write(6,*) '.....negative C',k
!      end do
      
#define MAKE_C_POS
#undef MAKE_C_POS
#ifdef MAKE_C_POS

    do k=1,Nspec
      CONC(k) = MAX(CONC(k),zero)
    end do 
#endif
      
    CALL CKWC(Z(1), CONC, IPAR(ckbi), RPAR(ckbr), WDOTS)
    CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), ENTHALPY)

!
!     Form governing equation
!

    THFAC = one / thickFacCH
    SUM = zero
    DO K = 1, Nspec
      H    = ENTHALPY(K)
      WDOT = WDOTS(K) * THFAC
      WT   = RPAR(NWT+K-1)
      ZP(K+1) = WDOT * WT / RHO
      SUM = SUM + H * WDOT * WT
    END DO
    ZP(1) = -SUM / (RHO*CPB)

#if 0
    print*, 'Z:'
    do k = 1, Nspec+1
      write(6,996) Z(K)
    end do
    print*, 'ZP:'
    do k = 1, Nspec+1
      write(6,996) ZP(K)
    end do
    print*, 'WDOT:'
    do k = 1, Nspec
      write(6,996) WDOTS(K)
    end do

996   format(e30.22)
#endif

  end subroutine conpFY

!------------------------------------  
  
  subroutine conpJY(N, TN, Y, SAVF, NFE, FTEM, ML, &
                       MU, PD, NRPD, RPAR, IPAR)
                       
    implicit none
      
    REAL_T SAVF
    REAL_T PD, RPAR(*), TN, Y
    dimension SAVF(*)
    integer N, NRPD, ML, MU, IPAR(*)
    dimension Y(N), PD(NRPD,N)
    REAL_T FTEM
    dimension FTEM(*)
    integer NFE
    
#include "cdwrk.H"
#include "conp.H"
     
    call bl_abort("conpJY: SHOULD NOT BE HERE!")
      
  END subroutine conpJY

!------------------------------------  
  
  subroutine conpFY_sdc(N, TIME, Z, ZP, RPAR, IPAR)

!
!     Variables in Z are:  Z(1:K) = rhoY(K) [MKS]
!                          Z(K+1) = RhoH    [MKS]

    implicit none

#include "cdwrk.H"
#include "conp.H"

    REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
    integer N, IPAR(*)
    
    REAL_T RHO_MKS, WDOT_CGS(maxspec), THFAC, HMIX_MKS, RINV_MKS
    REAL_T Y(maxspec), CONC_CGS(maxspec)
    integer K,Niter

    REAL_T,  parameter :: HtoTerrMAX = BL_REAL_E(7.8,-12)
    integer, parameter :: HtoTiterMAX = 20
    REAL_T res(0:HtoTiterMAX-1)

!
!     Variables in Z are:  Z(1:K) = rhoY(K)
!                          Z(K+1) = RhoH

    RHO_MKS = 0.d0
    do K=1,Nspec
      RHO_MKS = RHO_MKS + Z(K)
    enddo
    RINV_MKS = 1.d0 / RHO_MKS
      
    do K=1,Nspec
      CONC_CGS(K) = Z(K)*RPAR(NWTI+K-1)*1.d-3
      Y(K) = Z(K) * RINV_MKS
      if (Y(K) .lt. -1.d-3) then
        negative_Y_test = 1
      end if
    enddo

    HMIX_MKS = (rhoh_INIT + c_0(Nspec+1)*TIME) * RINV_MKS
    call TfromHYpt(T_cell,HMIX_MKS,Y,HtoTerrMAX,HtoTiterMAX,res,Niter)
    call CKWC(T_cell,CONC_CGS,IWRK,RWRK,WDOT_CGS)

    ZP(Nspec+1) = c_0(Nspec+1)
    THFAC = 1.d3/thickFacCH
    do k= 1, Nspec
      ZP(k) = WDOT_CGS(k) * RPAR(NWT+k-1) * THFAC + c_0(k)
    end do
      
  end subroutine conpFY_sdc

!------------------------------------  
  
  integer function TfromeYpt(T,ein,Y,errMax,NiterMAX,res)
  
    implicit none
      
#include "cdwrk.H"

    integer NiterMAX,Niter,n,NiterDAMP
    REAL_T T,ein,Y(*),errMAX,res(0:NiterMAX-1)
    REAL_T TMIN,TMAX,e

    parameter (TMIN=200, TMAX=5000)
    REAL_T T0,cv,de,temp
    REAL_T dT, etarg
    logical out_of_bounds, converged, soln_bad
    REAL_T e300,cv300,e6500,cv6500
    integer ihitlo,ihithi

    out_of_bounds(temp) = (temp.lt.TMIN) .or. (temp.gt.TMAX)

    NiterDAMP = NiterMAX
    if ((T.GE.TMIN).and.(T.LE.TMAX)) then
      T0 = T
    else
      T0 = half*(TMIN+TMAX)
      T = T0
    end if
    Niter = 0
    de = zero
    soln_bad = .FALSE.
    etarg = ein * BL_REAL_E(1.0,4)
    ihitlo = 0
    ihithi = 0

    CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)

    de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
    res(Niter) = de
    converged = de.le.errMAX

    do while ((.not.converged) .and. (.not.soln_bad))
      CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv)
      dT = (etarg - e)/cv
      if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
        T = TMAX
        ihithi = 1
      else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
        T = TMIN
        ihitlo = 1
      else
        T = T + dT
      end if
      soln_bad = out_of_bounds(T)
      if (soln_bad) then
        TfromeYpt = -1
        goto 100
      else
        CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)
        de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
        res(Niter) = de
        Niter = Niter + 1
      end if
      if (Niter .ge. NiterMAX) then
        TfromeYpt = -2
        goto 100
      endif
      converged = (de.le.errMAX) .or. (ABS(dT).le.errMAX)

      if((ihitlo.eq.1).and.(e.gt.etarg))then
        T = 300.d0
        CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e300)
        CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv300)
        T=300.d0+(etarg-e300)/cv300
        converged = .true.
      endif
      if((ihithi.eq.1).and.(e.lt.etarg))then
        T = 6500.d0
        CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e6500)
        CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv6500)
        T=6500.d0+(etarg-e6500)/cv6500
        converged = .true.
      endif

    end do

!     Set max iters taken during this solve, and exit
    TfromeYpt = Niter
    return

!     Error condition....dump state and bail out
 100  continue

      write(6,997) 'T from (e,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual = e - h + RT/Wbar [cgs]'
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)

  end  function TfromeYpt

 !------------------------------------ 
  
  subroutine TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)&
             bind(C, name="TfromHYpt")
             
      implicit none
      
#include "cdwrk.H"

      REAL_T T,Y(*),Hin,errMax

      integer NiterMAX,Niter,NiterDAMP,ihitlo,ihithi
      REAL_T  T0,cp,dH
      REAL_T  res(0:NiterMAX-1),dT, Htarg,HMIN,cpMIN,HMAX,cpMAX
      logical converged, soln_bad, stalled
      REAL_T  H, old_T, old_H, Tsec, Hsec

      integer, parameter :: Discont_NiterMAX = 100
      REAL_T,  parameter :: TMIN = 250.d0, TMAX = 5000.d0

      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T  = T0
      end if

      NiterDAMP = NiterMAX
      Niter     = 0
      soln_bad  = .FALSE.
      Htarg     = Hin * 1.d4
      ihitlo    = 0
      ihithi    = 0

      CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)

      old_T = T
      old_H = H

      dH         = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged  = dH.le.errMAX
      stalled    = .false.

      do while ((.not.converged) .and. (.not.stalled) .and. (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp)
         dT = (Htarg - H)/cp
         old_T = T
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = (T.lt.TMIN-one) .or. (T.gt.TMAX)
         if (soln_bad) then
            Niter = -1
            exit
         else
            old_H = H
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = min(dH,abs(dT))
            Niter = Niter + 1
         end if
         converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)
         if (Niter .ge. NiterMAX) then
            if(abs(T-1000.d0).le.1.d-3.and. dH.le.1.d-5)then
              converged = .true.
            else  
              Niter = -2
              exit
            endif
         endif

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),HMIN)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cpMIN)
            T=TMIN+(Htarg-HMIN)/cpMIN
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),HMAX)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cpMAX)
            T=TMAX+(Htarg-HMAX)/cpMAX
            converged = .true.
         endif

!c     If the iterations are failing, perhaps it is because the fits are discontinuous
!c     The following implements a modified secant method to hone in on a discontinity in h
!c     with T.  Discont_NiterMAX is fairly large because this process can be particularly
!c     slow to converge if the Htarg value happens to lay between the discontinuous function
!c     values.

         if (Niter .ge. NiterMAX) then
            do while (.not. stalled)
               dT = - (H - Htarg) * (old_T - T)/(old_H - H)
               Tsec = T + dT
               soln_bad = (Tsec.lt.TMIN-one) .or. (Tsec.gt.TMAX)
               if (soln_bad) then
                  Niter = -3
                  exit
               endif
               CALL CKHBMS(Tsec,Y,IWRK(ckbi),RWRK(ckbr),Hsec)
               if ( (Hsec-Htarg)*(Htarg-H) .gt. 0.d0 ) then
                  old_H = H
                  old_T = T
               endif
               H = Hsec
               T = Tsec
               stalled = (2*ABS(old_T-T)/(old_T+T).le.errMAX)
               Niter = Niter + 1
               if (Niter.gt.NiterMAX+Discont_NiterMAX) then
                  Niter = -2
                  exit
               endif
            enddo
            converged = .true.
         endif
      end do

      if (converged) return
!c
!c     Error condition....dump state and bail out
!c
#if 0
      print*, 'errMax: ', errMax
      write(6,997) 'T from (H,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'previous T = ',old_T
      write(6,998) 'target H = ',Htarg
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,NiterMAX-1
         write(6,998) '  ',res(n)
      end do
 997  format(a,3(i4,a))
 998  format(a,d21.12)
#endif
  end subroutine TfromHYpt
  
!------------------------------------  
  
  integer function open_vode_failure_file ()
  
      implicit none
      
      character*30 name, myproc
      integer lout,i,j,k,idx

!     Hardwire the unit number to 26 for the moment
      lout = 26 
      call bl_pd_myproc(i)
      write(myproc, *) i
      idx = 1 
      do j = 1, 30
         if (myproc(j:j) .ne. ' ') then
            idx = j
            goto 1
         end if 
      end do
 1    continue
      do k = 30, j+1, -1
         if (myproc(k:k) .ne. ' ') then
            goto 2
         end if
      end do
 2    continue
      write(name, '(2a)') 'vode.failed.', myproc(idx:k)
!      write(name, '(2a)') 'vode.failed.', myproc(idx:30)
      open(unit=lout, file=name, form='formatted', status='replace')
      open_vode_failure_file = lout
      
  end function open_vode_failure_file
      
end module chem_driver
