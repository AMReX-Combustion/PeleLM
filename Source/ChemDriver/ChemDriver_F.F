#include "ChemDriver_F.H"
#include "CONSTANTS.H"

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
      subroutine FORT_GET_REACTION_MAP(rmap)
      implicit none
#include "cdwrk.H"
      integer rmap(nReac)
      call GET_REACTION_MAP(rmap) 
      end

      subroutine FORT_SETTMINTRANS(TminTRANS)
      implicit none
      REAL_T TminTRANS
#include "cdwrk.H"
      TMIN_TRANS = TminTRANS
      end

      subroutine FORT_SETVERBOSEVODE()
      implicit none
#include "cdwrk.H"
      verbose_vode = 1
      end

      subroutine FORT_SETVODETOLS(rtol,atol,itol)
      implicit none
      integer itol
      REAL_T rtol,atol
#include "cdwrk.H"
      vode_itol = itol
      vode_rtol = rtol
      vode_atol = atol
      end

      subroutine FORT_SETVODESUBCYC(maxcyc)
      implicit none
      integer maxcyc
#include "cdwrk.H"
      max_vode_subcycles = maxcyc
      end

      subroutine FORT_SETSPECSCALY(name, nlength)
      implicit none
#include "cdwrk.H"
      integer nlength, name(nlength), i, j, maxlen
      REAL_T val
      parameter (maxlen=256)
      character filet*(maxlen)
      character*(maxspnml) spname, spinname
c      
c     Convert encoded names to strings, and open file
c
      if (nlength.GT.maxlen) then
         call bl_abort('FORT_SETSPECSCAL: scale file name too long')
      end if
      
      do i = 1, nlength
         filet(i:i) = char(name(i))
      end do
      open(unit=51,status='OLD',form='FORMATTED',
     &     file=filet(1:nlength),err=30)
      
 10   continue 
      read(51,*,end=20) spinname, val
      do j = 1,Nspec
         call get_spec_name(spname,j)
         if (spname .eq. spinname) then
            spec_scalY(j) = ABS(val)
         end if
      end do
      goto 10
 20   close(51)
      goto 40
 30   write(6,*) 'Trouble opening file = ',filet(1:nlength)
      call bl_abort(" ")
 40   continue 
      end

      subroutine FORT_INITCHEM()
      implicit none
#include "cdwrk.H"
#include "conp.H"
      integer n, RTOT, lout, egrlen, egilen, idummy(1)
      double precision rdummy(1)
      character*(maxspnml) name
      integer MAXFIT, NO, NFDIM, NT, NRANGE, NLITEMAX
      integer mcr, mci, FORT_USINGEG, FORT_USINGMC
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
      egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes
     &        + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
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
         call EGINICD(eg_nodes, lout, eg_IFLAG, eg_ITLS,
     &        EGRWRK, egrlen, EGIWRK, egilen)
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
      if (iN2.eq.-1)
     &     write(6,*) '.....warning: no N2 in chemistry species list'
      end

      subroutine FORT_FINALIZECHEM()
      implicit none
#include "cdwrk.H"
      integer FORT_USINGEG, FORT_USINGMC

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

      end

      SUBROUTINE EGINICD (NP, LOUT, IFLAG, ITLS, 
     &                    WEG, LWEG, IWEG, LIWEG)
C-----------------------------------------------------------------------
C
C     This subroutine initializes the pointers for the work arrays
C     WEG and IWEG and checks their length.
C     This subroutine should be called by the user once at the
C     beginning of the program.
C
C     Input
C     -----
C        NP        number of nodes
C        LOUT      output file number
C        IFLAG     flag for evaluating parameters and space allocation
C                  (see below)
C        ITLS      flag for space allocation (see below)
C        WEG       double precision work array for EGLIB
C        LWEG      length of WEG declared in main code
C        IWEG      integer work array for EGLIB
C        LIWEG     length of IWEG declared in main code
C
C        
C     The value of IFLAG and ITLS depends on the subroutines that
C     will be used as indicated by the following table
C
C
C     Subroutine     ITLS      IFLAG
C
C     EG*D(R)1         1         2
C     EG*D(R)2         1         2
C
C     EG*E1            0         1
C     EG*E2            1         2
C     EG*E3            1         3
C     EG*E4            1         3
C 
C     EG*K1            0         4
C     EG*K2            1         4
C     EG*K3            1         5
C     EG*K4            2         4
C     EG*K5            2         5
C     EG*K6            2         5
C  
C     EG*L1            0         1
C     EG*L2            1         6
C     EG*L3            1         7
C     EG*L4            2         6
C     EG*L5            2         7
C  
C     EG*LC1           1         7
C     EG*LC2           1         7
C     EG*LC3           2         7
C     EG*LC4           2         7
C  
C     EG*LTD(R)1       2         7
C     EG*LTD(R)2       2         7
C     EG*LTD(R)3       3         7
C     EG*LTD(R)4       3         7
C     EG*LTD(R)5       3         7
C     EG*LTD(R)6       3         7
C   
C     EG*TD(R)1        3         7
C  
C     EG*V(R)1         0         2
C
C
C     EGINI should be called with the highest possible values for
C     IFLAG and ITLS as read from the table.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
C     Check values for IFLAG and ITLS
C-----------------------------------------------------------------------
      IERROR = 0
      IF ( IFLAG .LT. 0 .OR. IFLAG .GT. 7 ) THEN
         WRITE(LOUT,'(1X,''IFLAG should be between 0 and 7'')')
         WRITE(LOUT,'(1X,''value read in EGini is'',2x,I5//)') IFLAG
         IERROR = 1
      ENDIF
      IF ( ITLS .LT. 0 .OR. ITLS .GT. 3 ) THEN
         WRITE(LOUT,'(1X,''ITLS should be between 0 and 3'')')
         WRITE(LOUT,'(1X,''value read in EGini is'',2x,I5//)') ITLS
         IERROR = 1
      ENDIF
      IF ( IERROR .EQ. 1 ) STOP
C-----------------------------------------------------------------------
C     Read the Linkeg file
C-----------------------------------------------------------------------
      LLEG   = 11
C-----------------------------------------------------------------------
c      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
c        READ (LLEG) NSLK, NO
c      CLOSE(UNIT=LLEG)

      call egtransetKK(NSLK)
      call egtransetNO(NO)

C-----------------------------------------------------------------------
C     Store IFLAG and the number of species in common 'eg.cmn'
C-----------------------------------------------------------------------
      JFLAG = IFLAG
C-----------------------------------------------------------------------
      NS = NSLK
      IF ( NS .LE. 1 ) THEN
         WRITE(LOUT,'(1X,''Error: the number of species must '',
     &                   '' be larger or equal to 2'')')
         STOP
      ENDIF
C-----------------------------------------------------------------------
C     Compute the size of the transport linear system.
C-----------------------------------------------------------------------
      NSS  = ITLS * NS
C-----------------------------------------------------------------------
C     NFIT is the degree for the polynomial fitting Aij -- Cij
C-----------------------------------------------------------------------
      NFIT = 7
      NAIJ = MAX0(NS*NS,NP)
C-----------------------------------------------------------------------
      IEGRU  = 1
      IEGPA  = IEGRU  + 1
      IFITA  = IEGPA  + 1
      IFITB  = IFITA  + NFIT * NS*NS
      IFITC  = IFITB  + NFIT * NS*NS
      IFITA0 = IFITC  + NFIT * NS*NS
      IFITB0 = IFITA0 + NFIT
      IFITC0 = IFITB0 + NFIT
      ICTAIJ = IFITC0 + NFIT
      ICTBIJ = ICTAIJ + NS*NS
      ICTCIJ = ICTBIJ + NS*NS
      IDLT1  = ICTCIJ + NS*NS
      IDLT2  = IDLT1  + NP
      IDLT3  = IDLT2  + NP
      IDLT4  = IDLT3  + NP
      IDLT5  = IDLT4  + NP
      IDLT6  = IDLT5  + NP
      IEGEPS = IDLT6  + NP
      IEGPOL = IEGEPS + NS
      IEGSIG = IEGPOL + NS
      IEGDIP = IEGSIG + NS
      IEPSIJ = IEGDIP + NS
      IEGCFD = IEPSIJ + NS*NS
      IEGCFE = IEGCFD + 4 * NS*NS
      IEGCFL = IEGCFE + 4 * NS
      IEGZRT = IEGCFL + 4 * NS
      IEGWT  = IEGZRT + NS 
      IAAA   = IEGWT  + NS
      IBBB   = IAAA   + NP
      IAUX   = IBBB   + NP
      IETA   = IAUX   + NS * NP
      IETALG = IETA   + NS * NP
      IXTR   = IETALG + NS * NP
      IYTR   = IXTR   + NS * NP
      IAIJ   = IYTR   + NS * NP
      IBIJ   = IAIJ   + NAIJ
      ICIJ   = IBIJ   + NAIJ
      IBIN   = ICIJ   + NAIJ
      ICINT  = IBIN   + (NS*(NS+1))/2 * NP
      ICXI   = ICINT  + NS * NP
      IEND   = ICXI   + NS * NP
C.....
      IF ( IFLAG .EQ. 1 ) THEN
         IEND = IBIN
      ELSEIF ( IFLAG .LE. 3 ) THEN
         IEND = ICINT
      ENDIF
C.....
      IDMI   = IEND
      IG     = IDMI   + NS*(ITLS*(ITLS+1))/2 * NP
      IAN    = IG     + (ITLS*NS*(ITLS*NS+1))/2 * NP
      IZN    = IAN    + NSS * NP
      IRN    = IZN    + NSS * NP
      ITEMP  = IRN    + NSS * NP
      IBETA  = ITEMP  + NSS * NP
      INEXT  = IBETA  + NSS * NP - 1
C.....
      IEGLIN = 1
      IINXT  = IEGLIN + NS - 1
C-----------------------------------------------------------------------
      ILOW = 0
      IF ( INEXT .GT. LWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of WEG should be '',
     &                   ''at least'',I12//)') INEXT
         ILOW = 1
      ENDIF
      IF ( IINXT .GT. LIWEG ) THEN
         WRITE(LOUT,'(//1X,''Error: the length of IWEG should be '',
     &                   ''at least'',I12//)') IINXT
         ILOW = 1
      ENDIF
      IF ( ILOW .EQ. 1 ) STOP
c     WRITE(LOUT,'(//1X,''The array WEG requires '',I12,
c    &                '' storage locations'')') INEXT
c     WRITE(LOUT,'(1X,''The array IWEG requires '',I12,
c    &                '' storage locations''//)') IINXT
C-----------------------------------------------------------------------
C     Store the universal gas constant and the atmospheric pressure
C     units: [erg/mol.K] for RU and [dyne/cm^2] for PA
C-----------------------------------------------------------------------
      WEG(IEGRU)  = 8.314D7
      WEG(IEGPA)  = 1.01325D6
C-----------------------------------------------------------------------
C     Read the Linkeg file
C-----------------------------------------------------------------------
      LLEG   = 11
      NONS   = NO*NS
      NONSNS = NO*NS*NS
C-----------------------------------------------------------------------
c
c     Set required data from funcs rather than Linkeg
c
      call egtransetWT(WEG(IEGWT))
      call egtransetEPS(WEG(IEGEPS))
      call egtransetSIG(WEG(IEGSIG))
      call egtransetDIP(WEG(IEGDIP))
      call egtransetPOL(WEG(IEGPOL))
      call egtransetZROT(WEG(IEGZRT))
      call egtransetNLIN(IWEG(IEGLIN))
      call egtransetCOFETA(WEG(IEGCFE))
      call egtransetCOFLAM(WEG(IEGCFL))
      call egtransetCOFD(WEG(IEGCFD))

c      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
c        READ (LLEG) NSLK, NO, (WEG(IEGWT+K-1), K=1, NS), 
c     &              (WEG(IEGEPS+K-1), K=1, NS), 
c     &              (WEG(IEGSIG+K-1), K=1, NS), 
c     &              (WEG(IEGDIP+K-1), K=1, NS), 
c     &              (WEG(IEGPOL+K-1), K=1, NS), 
c     &              (WEG(IEGZRT+K-1), K=1, NS), 
c     &              (IWEG(IEGLIN+K-1), K=1, NS), 
c     &              (WEG(IEGCFE+N-1), N=1, NONS), 
c     &              (WEG(IEGCFL+N-1), N=1, NONS), 
c     &              (WEG(IEGCFD+N-1), N=1, NONSNS)
c      CLOSE(UNIT=LLEG)
C-----------------------------------------------------------------------
      CALL LEVEPS (NS, WEG(IEGEPS), WEG(IEGSIG), WEG(IEGDIP), 
     &             WEG(IEGPOL), WEG(IEPSIJ) )
C-----------------------------------------------------------------------
C     Initialize the coefficients for fitting Aij, Bij and Cij
C-----------------------------------------------------------------------
      WEG(IFITA0    ) =  .1106910525D+01
      WEG(IFITA0 + 1) = -.7065517161D-02
      WEG(IFITA0 + 2) = -.1671975393D-01
      WEG(IFITA0 + 3) =  .1188708609D-01
      WEG(IFITA0 + 4) =  .7569367323D-03
      WEG(IFITA0 + 5) = -.1313998345D-02
      WEG(IFITA0 + 6) =  .1720853282D-03
C.....
      WEG(IFITB0    ) =  .1199673577D+01
      WEG(IFITB0 + 1) = -.1140928763D+00
      WEG(IFITB0 + 2) = -.2147636665D-02
      WEG(IFITB0 + 3) =  .2512965407D-01
      WEG(IFITB0 + 4) = -.3030372973D-02
      WEG(IFITB0 + 5) = -.1445009039D-02
      WEG(IFITB0 + 6) =  .2492954809D-03
C.....
      WEG(IFITC0    ) =  .8386993788D+00
      WEG(IFITC0 + 1) =  .4748325276D-01
      WEG(IFITC0 + 2) =  .3250097527D-01
      WEG(IFITC0 + 3) = -.1625859588D-01
      WEG(IFITC0 + 4) = -.2260153363D-02
      WEG(IFITC0 + 5) =  .1844922811D-02
      WEG(IFITC0 + 6) = -.2115417788D-03
C-----------------------------------------------------------------------
C     Evaluate Aij, Bij and Cij at the reference temperature of 1000K.
C-----------------------------------------------------------------------
      DDD = DLOG(1.0D3)
      DO J = 1, NS
        DO I = 1, NS
         IJ = (J-1) * NS + I-1
         TSLOG = DDD - WEG( IEPSIJ + IJ )
         T1 = TSLOG
         T2 = TSLOG * T1
         T3 = TSLOG * T2
         T4 = TSLOG * T3
         T5 = TSLOG * T4
         T6 = TSLOG * T5
         WEG(ICTAIJ+IJ) = WEG(IFITA0  )    + WEG(IFITA0+1)*T1
     1                    + WEG(IFITA0+2)*T2 + WEG(IFITA0+3)*T3
     2                    + WEG(IFITA0+4)*T4 + WEG(IFITA0+5)*T5
     3                    + WEG(IFITA0+6)*T6 
         WEG(ICTBIJ+IJ) = WEG(IFITB0  )    + WEG(IFITB0+1)*T1
     1                    + WEG(IFITB0+2)*T2 + WEG(IFITB0+3)*T3
     2                    + WEG(IFITB0+4)*T4 + WEG(IFITB0+5)*T5
     3                    + WEG(IFITB0+6)*T6 
         WEG(ICTCIJ+IJ) = WEG(IFITC0  )    + WEG(IFITC0+1)*T1
     1                    + WEG(IFITC0+2)*T2 + WEG(IFITC0+3)*T3
     2                    + WEG(IFITC0+4)*T4 + WEG(IFITC0+5)*T5
     3                    + WEG(IFITC0+6)*T6 
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     Evaluate FITA, FITB and FITC 
C-----------------------------------------------------------------------
      CALL EGABC ( NS, NFIT, WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &             WEG(IFITA0), WEG(IFITB0), WEG(IFITC0),
     &             WEG(IEPSIJ) )
C-----------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE MCINITCD (LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK,
     1           IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE MCINIT (LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK,
C                     IFLAG)
C  This subroutine reads the transport linkfile from the fitting code
C  and creates the internal storage and work arrays, IMCWRK(*) and
C  RMCWRK(*).  MCINIT must be called before any other transport
C  subroutine is called.  It must be called after the CHEMKIN package
C  is initialized.
C
C  INPUT
C  LOUT      - Integer scalar, formatted output file unit number.
C  LENIMC    - Integer scalar, minimum dimension of the integer
C              storage and workspace array IMCWRK(*);
C              LENIMC must be at least:
C              LENIMC = 4*KK + NLITE,
C              where KK is the total species count, and
C                    NLITE is the number of species with molecular
C                          weight less than 5.
C  LENRMC    - Integer scalar, minimum dimension of the real storage
C              and workspace array RMCWRK(*);
C              LENRMC must be at least:
C              LENRMC = KK*(19 + 2*NO + NO*NLITE) + (NO+15)*KK**2,
C              where KK is the total species count,
C                    NO is the order of the polynomial fits (NO=4),
C                    NLITE is the number of species with molecular
C                          weight less than 5.
C
C  OUTPUT
C  IMCWRK(*) - Integer workspace array; dimension at least LENIMC.
C  RMCWRK(*) - Real    workspace array; dimension at least LENRMC.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      include 'mccom.fh'
C
      DIMENSION IMCWRK(*), RMCWRK(*)
      CHARACTER*16 PRVERS, PRDATE, IFMT, RFMT, CFMT, LFMT
      PARAMETER
     1(CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
C
      LOGICAL IOK, ROK, LBIN
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C     The following number SMALL is used in the mixture diffusion
C     coefficient calculation; its use allows a smooth and well-
C     defined diffusion coefficient as the mixture approaches a pure
C     species, even though stictlyh speaking there does not exist a
C     diffusion coefficient in this case.  The value of SMALL should
C     be small relative to any species mole fraction of importance,
C     but large enough to be represented on the computer.
C
C*****SMALL 1) 64 bit floats
C      SMALL = 1.0E-50
C*****END SMALL 1) 64 bit floats
C*****SMALL 2) 32 bit floats
      SMALL = 1.0E-20
C*****END SMALL 2) 32 bit floats
C
C     Gas constant as reported in 1993 CRC, (J. Research of
C     National Bureau of Standards, 92, 85, 1987).
C     ( 8.314510(70)E+07 Joules mol-1 K-1)
C
      RU    = 8.314510E+07
C
C     Standard atmosphere (defined as an exact quantity)
C
      PATMOS= 1.01325E+06
C
C     Write version number
C
C*****precision > double
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C      PREC = 'SINGLE'
C*****END precision > single
C
      PRVERS ='4.3'
      PRDATE ='98/03/03'
C
c      WRITE (LOUT, '(/A, /1X,A, A, A, A, /A, /A, //)')
c     1 ' TRANLIB:  CHEMKIN-III MULTICOMPONENT TRANSPORT LIBRARY,',
c     2 PREC(1:CKLSCH(PREC)), ' PRECISION Vers. ',
c     3 PRVERS(1:CKLSCH(PRVERS)+1), PRDATE,
c     4 ' Copyright 1995, Sandia Corporation.',
c     5' The U.S. Government retains a limited license in this software.'
C
C     Read the problem size
c      CALL MCLEN (LINKMC, LOUT, LI, LR, IFLAG)

c     Set required data from funcs rather than tran.asc
      call egtransetLENIMC(LI)
      call egtransetLENRMC(LR)
      IFLAG = 0

      IOK = (LENIMC .GE. LI)
      ROK = (LENRMC .GE. LR)
C
      IF (.NOT.IOK .OR. .NOT.ROK) THEN
         IF (.NOT. IOK) WRITE (LOUT, 300) LI
         IF (.NOT. ROK) WRITE (LOUT, 350) LR
c         REWIND (LINKMC)
         IFLAG = 1
         RETURN
      ENDIF
C
c      REWIND LINKMC
C*****linkfile (transport) > binary
C      LBIN = .TRUE.
C*****END linkfile (transport) > binary
C*****linkfile (transport) > ascii
      LBIN = .FALSE.
C*****END linkfile (transport) > ascii
C
c      NREC = 1
c      IF (LBIN) THEN
c         READ (LINKMC, ERR=999) VERS
c         NREC = 2
c         READ (LINKMC, ERR=999) PRVERS
c         NREC = 3
c         READ (LINKMC, ERR=999) PREC
c         NREC = 4
c         READ (LINKMC, ERR=999) KERR
c         NREC = 5
c         READ (LINKMC, ERR=999) LI, LR, NO, NKK, NLITE
c         NREC = 6
c         READ (LINKMC, ERR=999) PATMOS
c      ELSE
c         READ (LINKMC, CFMT, ERR=999) VERS
c         NREC = 2
c         READ (LINKMC, CFMT, ERR=999) PRVERS
c         NREC = 3
c         READ (LINKMC, CFMT, ERR=999) PREC
c         NREC = 4
c         READ (LINKMC, LFMT, ERR=999) KERR
c         NREC = 5
c         READ (LINKMC, IFMT, ERR=999) LI, LR, NO, NKK, NLITE
c         NREC = 6
c         READ (LINKMC, RFMT, ERR=999) PATMOS
c      ENDIF
C
c     Do not bother to build functions to translate strings, assume this is correct
      VERS = '1.0'
      PRVERS = '3.11'
      PREC = 'DOUBLE'
      KERR = .false.
      call egtransetNO(NO)
      call egtransetKK(NKK)
      call egtransetNLITE(NLITE)
      call egtransetPATM(PATMOS)

      NK  = NO*NKK
      NK2 = NO*NKK*NKK
      K2  = NKK*NKK
      K3  = 3*NKK
      K32 = K3*K3
      NKT = NO*NKK*NLITE
C
C     APPORTION THE REAL WORKSPACE:
C
C     molecular weights for the species
      NWT  = 1
C     the epsilon/k well depth for the species
      NEPS = NWT + NKK
C     the collision diameter for the species
      NSIG = NEPS + NKK
C     the dipole moments for the species
      NDIP = NSIG + NKK
C     the polarizabilities for the species
      NPOL = NDIP + NKK
C     the rotational relaxation collision numbers
      NZROT= NPOL + NKK
C     the fit coefficients for conductivity
      NLAM = NZROT + NKK
C     the fit coefficients for viscosity
      NETA = NLAM + NK
      NDIF = NETA + NK
C     the fit coefficients for thermal diffusion ratio
      NTDIF= NDIF + NK2
C     mole fractions of the mixture
      NXX  = NTDIF + NO*NKK*NLITE
C     species viscosities
      NVIS = NXX + NKK
C     rotational relaxation collision numbers before Parker coffection
      NXI  = NVIS + NKK
C     species specific heats
      NCP  = NXI + NKK
C     rotational parts of the specific heats
      NCROT= NCP + NKK
C     internal parts of the specific heats
      NCINT= NCROT + NKK
C     the binary diffusion coefficients
      NBIND= NCINT + NKK
C     the matrix of reduced well depths
      NEOK = NBIND + K2
C     the matrix of reduced collision diameters
      NSGM = NEOK + K2
C     the matrix of A* collision integrals for each species pair
      NAST = NSGM + K2
C     the matrix of B* collision integrals for each species pair
      NBST = NAST + K2
C     the matrix of C* collision integrals for each species pair
      NCST = NBST + K2
C     the "L" matrix
      NXL  = NCST + K2
C     the right-hand sides of the linear system involving the
C     "L" matrix
      NR   = NXL + K32
C     the workspace needed by LINPACK to solve the "L" matrix linear
C     system
      NWRK = NR + K3
C      NTOT = NWRK + K3 - 1
C
C     APPORTION THE INTEGER WORKSPACE:
C
C     the indicators for the molecule linearity
      INLIN = 1
C     the species indices for the "light" species
      IKTDIF= INLIN + NKK
C     the pivot indices for LINPACK calls
      IPVT  = IKTDIF + NLITE
C      ITOT  = IPVT + K3 - 1
C
C     Read the data from the linkfile
C
c$$$      IF (LBIN) THEN
c$$$         NREC = 7
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
c$$$         NREC = 8
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
c$$$         NREC = 9
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
c$$$         NREC = 10
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
c$$$         NREC = 11
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
c$$$         NREC = 12
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
c$$$         NREC = 13
c$$$         READ (LINKMC, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
c$$$         NREC = 14
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
c$$$         NREC = 15
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
c$$$         NREC = 16
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
c$$$         NREC = 17
c$$$         READ (LINKMC, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
c$$$         NREC = 18
c$$$         READ (LINKMC, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
c$$$      ELSE
c$$$         NREC = 7
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
c$$$         NREC = 8
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
c$$$         NREC = 9
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
c$$$         NREC = 10
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
c$$$         NREC = 11
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
c$$$         NREC = 12
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
c$$$         NREC = 13
c$$$         READ (LINKMC, IFMT, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
c$$$         NREC = 14
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
c$$$         NREC = 15
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
c$$$         NREC = 16
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
c$$$         NREC = 17
c$$$         READ (LINKMC, IFMT, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
c$$$         NREC = 18
c$$$         READ (LINKMC, RFMT, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
c$$$         NREC = 19
c$$$      ENDIF

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NWT+N-1), N = 1, NKK)
      call egtransetWT(RMCWRK(NWT))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NEPS+N-1), N = 1, NKK)
      call egtransetEPS(RMCWRK(NEPS))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NSIG+N-1), N = 1, NKK)
      call egtransetSIG(RMCWRK(NSIG))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIP+N-1), N = 1, NKK)
      call egtransetDIP(RMCWRK(NDIP))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NPOL+N-1), N = 1, NKK)
      call egtransetPOL(RMCWRK(NPOL))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NZROT+N-1), N = 1, NKK)
      call egtransetZROT(RMCWRK(NZROT))

c      READ (LINKMC, IFMT, ERR=999) (IMCWRK(INLIN+N-1), N = 1, NKK)
      call egtransetNLIN(IMCWRK(INLIN))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NLAM+N-1), N = 1, NK)
      call egtransetCOFLAM(RMCWRK(NLAM))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NETA+N-1), N = 1, NK)
      call egtransetCOFETA(RMCWRK(NETA))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NDIF+N-1), N = 1, NK2)
      call egtransetCOFD(RMCWRK(NDIF))

c      READ (LINKMC, IFMT, ERR=999) (IMCWRK(IKTDIF+N-1), N = 1, NLITE)
      call egtransetKTDIF(IMCWRK(IKTDIF))

c      READ (LINKMC, RFMT, ERR=999) (RMCWRK(NTDIF+N-1), N = 1, NKT)
      call egtransetCOFTD(RMCWRK(NTDIF))

C
C     Set EPS/K and SIG for all I,J pairs
C
      CALL MCEPSG (NKK, RMCWRK(NEPS), RMCWRK(NSIG), RMCWRK(NDIP),
     1            RMCWRK(NPOL), RMCWRK(NEOK), RMCWRK(NSGM) )
C
  300 FORMAT (10X,'IMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
  350 FORMAT (10X,'RMCWRK MUST BE DIMENSIONED AT LEAST ', I5)
      RETURN
C  999 CONTINUE
c      WRITE (LOUT, *) ' Error reading Transport linkfile...'
c      REWIND (LINKMC)
      IFLAG = NREC
C
C     end of SUBROUTINE MCINITCD
      RETURN
      END

      integer function FORT_GETCKMAXNAMELEN()
      implicit none
#include "cdwrk.H"
      FORT_GETCKMAXNAMELEN = maxspnml
      end

      subroutine FORT_GETCKDIMPARAMS(imaxreac, imaxspec, imaxelts,
     &                               imaxord, imaxthrdb, imaxtp, imaxsp,
     &                               imaxspnml)
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
      end

      subroutine FORT_FINDLHS(reactions, Nreacs, id)
      implicit none
      integer reactions(*), Nreacs, id

#ifdef MIKE
#include "cdwrk.H"
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
      Ndim = maxsp
      if ((id.le.0).or.(id.gt.Nspec)) then
         write(6,*) 'FORT_FINDLHS:  species id out of range: ',id
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
      call bl_abort("FORT_FINDLHS not implemented")
#endif
      end

      subroutine FORT_FINDRHS(reactions, Nreacs, id)
      implicit none
      integer reactions(*), Nreacs, id

#ifdef MIKE
#include "cdwrk.H"
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
      Ndim = maxsp
      if ((id.le.0).or.(id.gt.Nspec)) then
         write(6,*) 'FORT_FINDRHS:  species id out of range: ',id
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
      call bl_abort("FORT_FINDRHS not implemented")
#endif
      end

      subroutine FORT_SETNU(nu,lenNU)
      implicit none
#include "cdwrk.H"
      integer lenNU
      integer nu(maxreac,maxspec)

      if (lenNU .lt. maxreac*maxspec) then
         write(6,*) 'FORT_CKNU:  nu work array too small: '
         call bl_abort(" ")
      endif
      call CKNU(maxreac, IWRK(ckbi), RWRK(ckbr), nu)
      end

      subroutine FORT_CKINU(Nids,KI,lenKI,NU,lenNU,rxnID,nuAll)
      implicit none
#include "cdwrk.H"
      integer lenKI,lenNU
      integer rxnID, Nids, KI(lenKI), NU(lenNU), nuAll(maxreac,maxspec)
      integer Ndim, k
      Ndim = MIN(lenKI,lenNU)
      if ((rxnID.le.0).or.(rxnID.gt.Nreac)) then
         write(6,*) 'FORT_CKINU:  reaction id out of range: ',rxnID
         call bl_abort(" ")
      end if
      if (Ndim.lt.maxsp) then
         call bl_abort('FORT_CKINU:  KI or NU not long enough')
      end if
      Nids = 0
      do k=1,Nspec
         if (nuAll(rxnID,k).ne.0) then
            Nids = Nids + 1
            KI(Nids) = k
            NU(Nids) = nuAll(rxnID,k)
         endif
      enddo
      end

      integer function FORT_CKELTXINSPY(eltID, spID)
      implicit none
#include "cdwrk.H"
      integer eltID, spID
      integer NCF(maxelts,maxspec)
      CALL CKNCF(maxelts, IWRK(ckbi), RWRK(ckbr), NCF)
      FORT_CKELTXINSPY = NCF(eltID+1,spID+1)
      end

      integer function FORT_GETCKNUMSPEC()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMSPEC = Nspec
      end

      integer function FORT_GETCKNUMELT()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMELT = Nelt
      end

      integer function FORT_GETCKNUMREAC()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMREAC = Nreac
      end

      double precision function FORT_RUNIV()
      implicit none
      double precision Ruc, Pa
#include "cdwrk.H"
      call CKRP(IWRK(ckbi),RWRK(ckbr),FORT_RUNIV,Ruc,Pa)
c     1 erg/(mole.K) = 1.e-4 J/(kmole.K)
      FORT_RUNIV = FORT_RUNIV*1.d-4
      end

      double precision function FORT_P1ATMMKS()
      implicit none
      double precision Ru, Ruc, Pa
#include "cdwrk.H"
      call CKRP(IWRK(ckbi),RWRK(ckbr),Ru,Ruc,Pa)
c     1 N/(m.m) = 0.1 dyne/(cm.cm)
      FORT_P1ATMMKS = Pa*1.d-1
      end

      integer function FORT_GETCKELTNAME(i, coded)
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
         FORT_GETCKELTNAME = 1
      else
         FORT_GETCKELTNAME = 2
      endif
      end

      integer function FORT_GETCKSPECNAME(i, coded)
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
      FORT_GETCKSPECNAME = str_len - 1
      end

      integer function FORT_CKSYMR(fortReacIdx, coded)
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
      call CKSYMR(fortReacIdx,lout,IWRK(ckbi),RWRK(ckbr),
     &     CWRK(ckbc),str_len,line,error)
      if (error) then
         write(lout,*) 'Could not get reaction name for ',fortReacIdx
         call bl_abort(" ")
      end if
c      
c     Encode the name for transfer to C++
c
      istr = 1
      do while (line(istr:istr) .EQ. ' ')
         istr = istr + 1
      end do
      do j = 0, str_len-1
         coded(j+1) = ICHAR(line(istr+j:istr+j))
      end do
      FORT_CKSYMR = str_len
#else
      FORT_CKSYMR = 0
      call bl_abort("FORT_CKSYMR not implemented")
#endif
      end

      subroutine get_spec_name(name, j)
      implicit none
#include "cdwrk.H"
      integer i, j, FORT_GETCKSPECNAME
      integer coded(maxspnml), len
      character*(maxspnml) name
      len = FORT_GETCKSPECNAME(j, coded)
      do i = 1, maxspnml
         name(i:i) = ' '
      end do
      do i = 1, len
         name(i:i) = char(coded(i))
      end do
      end

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
      end

      subroutine FORT_GETCKMWT(mwt)
      implicit none
#include "cdwrk.H"
      REAL_T mwt(*)
c     Result in kg/kmole
      call CKWT(IWRK(ckbi),RWRK(ckbr),mwt)
      end

      subroutine FORT_GETCKAWT(awt)
      implicit none
#include "cdwrk.H"
      REAL_T awt(*)
c     Result in kg/kmole
      call CKAWT(IWRK(ckbi),RWRK(ckbr),awt)
      end

      subroutine conpFY(N, TIME, Z, ZP, RPAR, IPAR)
      implicit none
#include "cdwrk.H"
#include "conp.H"

      REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer N, IPAR(*)
      
      REAL_T RHO, CPB, SUM, H, WDOT, WT, THFAC
      integer K

      REAL_T CONC(maxspec), WDOTS(maxspec), ENTHALPY(maxspec)
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
      
      CALL CKRHOY(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),RHO)
      CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)
      CALL CKYTCP(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CONC)
c
c     Get net production rates.  Compute such that production from -ve
c     reactants gives zero contrib.
c      
c      do k=1,Nspec
c         if (CONC(k) .lt. zero) write(6,*) '.....negative C',k
c      end do
      
#define MAKE_C_POS
#undef MAKE_C_POS
#ifdef MAKE_C_POS
      do k=1,Nspec
         CONC(k) = MAX(CONC(k),zero)
      end do 
#endif
      
      CALL CKWC(Z(1), CONC, IPAR(ckbi), RPAR(ckbr), WDOTS)
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), ENTHALPY)
C
C     Form governing equation
C
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

      END

      subroutine conpJY(N, TN, Y, SAVF, NFE, FTEM, ML,
     &                  MU, PD, NRPD, RPAR, IPAR)
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
      END

#ifdef LMC_SDC
      subroutine conpFY_sdc(N, TIME, Z, ZP, RPAR, IPAR)
C
C     Variables in Z are:  Z(1:K) = rhoY(K) [MKS]
C                          Z(K+1) = RhoH    [MKS]
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
C
C     Variables in Z are:  Z(1:K) = rhoY(K)
C                          Z(K+1) = RhoH
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
      call FORT_TfromHYpt(T_cell,HMIX_MKS,Y,HtoTerrMAX,HtoTiterMAX,res,Niter)
      call CKWC(T_cell,CONC_CGS,IWRK,RWRK,WDOT_CGS)

      ZP(Nspec+1) = c_0(Nspec+1)
      THFAC = 1.d3/thickFacCH
      do k= 1, Nspec
         ZP(k) = WDOT_CGS(k) * RPAR(NWT+k-1) * THFAC + c_0(k)
      end do
      END
#endif

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

c     Set max iters taken during this solve, and exit
      TfromeYpt = Niter
      return

c     Error condition....dump state and bail out
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

      end

      subroutine FORT_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)
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

c     If the iterations are failing, perhaps it is because the fits are discontinuous
c     The following implements a modified secant method to hone in on a discontinity in h
c     with T.  Discont_NiterMAX is fairly large because this process can be particularly
c     slow to converge if the Htarg value happens to lay between the discontinuous function
c     values.  
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
c
c     Error condition....dump state and bail out
c
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
      end
  
      integer function open_vode_failure_file ()
      implicit none
      character*30 name, myproc
      integer lout,i,j,k,idx

c     Hardwire the unit number to 26 for the moment
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
c      write(name, '(2a)') 'vode.failed.', myproc(idx:30)
      open(unit=lout, file=name, form='formatted', status='replace')
      open_vode_failure_file = lout
      end
