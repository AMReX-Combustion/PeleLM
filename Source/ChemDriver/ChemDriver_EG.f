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