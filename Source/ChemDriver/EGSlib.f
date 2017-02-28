C-----------------------------------------------------------------------
      SUBROUTINE EGSPAR ( T, X, Y, CP, WEG, IWEG )
C-----------------------------------------------------------------------
C
C     This subroutine initializes the thermomolecular
C     parameters that are needed in order to evaluate
C     the transport linear systems.
C     The parameters that have to be evaluated depend on the transport
C     coefficients that will be subsequently computed, as indicated
C     by JFLAG. This flag has been initialized when calling EGINI.
C
C     Input
C     -----
C        T         temperature
C        X(NS)     species mole fractions
C        Y(NS)     species mass fractions
C        CP(NS)    species heat capacities at constant pressure
C                  per unit mass
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION WEG(*), IWEG(*), X(*), Y(*)
      PARAMETER(SSS = 1.0D-16)
C-----------------------------------------------------------------------
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LSPAR ( JFLAG, WEG(IEGPA), 
     &       T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3), 
     &       WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &       WEG(IEGRU), NS, 
     &       WEG(IEGWT), WEG(IBIN), WEG(IETA), 
     &       WEG(IETALG), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &       WEG(ICTAIJ), WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &       WEG(ICINT), CP, WEG(ICXI), WEG(IEPSIJ), 
     &       WEG(IEGEPS), WEG(IEGCFD), WEG(IEGCFE), 
     &       WEG(IEGZRT), WEG(IEGDIP), IWEG(IEGLIN) )
C-----------------------------------------------------------------------
C     Add a small constant to the mole and mass fractions
C-----------------------------------------------------------------------
      AAA = 1.0D0 / DFLOAT(NS)
      DO I = 1, NS
         WEG ( IXTR + I - 1 ) = X(I) + SSS * ( AAA - X(I) )
         WEG ( IYTR + I - 1 ) = Y(I) + SSS * ( AAA - Y(I) )
      ENDDO
C-----------------------------------------------------------------------
C     AUX(i) = \sum_{j .ne. i} YTR(j)
C-----------------------------------------------------------------------
      CALL EGZERO ( NS, WEG(IAUX) )
      DO I = 1, NS
         DO J = I+1, NS
            WEG(IAUX + I-1) = WEG(IAUX + I-1) + WEG(IYTR + J-1)
            WEG(IAUX + J-1) = WEG(IAUX + J-1) + WEG(IYTR + I-1)
         ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LSPAR ( IFLAG, PATMOS, T, DLT1, DLT2, DLT3,
     &                   DLT4, DLT5, DLT6, RU, NS, WT, BIN, ETA,
     &                   ETALG, AIJ, BIJ, CIJ, CTAIJ, FITA, FITB, FITC,
     &                   CINT, CPMS, CXI, EPSIJ, 
     &                   EPS, COFD, COFE, ZROT, DIP, LIN )
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION BIN(*), ETA(NS), ETALG(NS), AIJ(*), CTAIJ(*),
     &          BIJ(*), CIJ(*), WT(*), 
     &          EPSIJ(*), LIN(NS), ZROT(NS), EPS(NS), DIP(NS)
      DIMENSION CPMS(NS), CINT(NS), CXI(NS)
      DIMENSION FITA(7,NS,NS), FITB(7,NS,NS), FITC(7,NS,NS)
C-----------------------------------------------------------------------
      PARAMETER(PI = 3.1415926535D0, PI32 = 5.5683D0)
      PARAMETER(PI32O2 = 2.7842D0, P2O4P2 = 4.4674D0)
C-----------------------------------------------------------------------
C     logarithm of the temperature and its powers
C-----------------------------------------------------------------------
      dlt1 = dlog ( t )
      dlt2 = dlt1 * dlt1
      dlt3 = dlt2 * dlt1
      dlt4 = dlt3 * dlt1
      dlt5 = dlt4 * dlt1
      dlt6 = dlt5 * dlt1
C-----------------------------------------------------------------------
      CALL EGSCOFE ( COFE, NS, DLT1, DLT2, DLT3, ETA, ETALG )
      IF ( IFLAG .LE. 1 ) RETURN
C-----------------------------------------------------------------------
      CALL EGSCOFD ( COFD, NS, DLT1, DLT2, DLT3, BIN )
      IF (IFLAG .LE. 2) RETURN
C-----------------------------------------------------------------------
C
C          DETERMINE A*, B*, AND C* FOR EACH SPECIES PAIR
C
C-----------------------------------------------------------------------
      IF ( IFLAG .EQ. 3 .OR. IFLAG .EQ. 5 ) THEN
         DO J = 1, NS
            DO I = J, NS
               IJ = (J-1) * NS + I
               FITA0 = FITA(1,I,J)
               FITA1 = FITA(2,I,J)
               FITA2 = FITA(3,I,J)
               FITA3 = FITA(4,I,J)
               FITA4 = FITA(5,I,J)
               FITA5 = FITA(6,I,J)
               FITA6 = FITA(7,I,J)
               AIJ(IJ) = FITA0 + FITA1*DLT1 + FITA2*DLT2 + FITA3*DLT3
     &                         + FITA4*DLT4 + FITA5*DLT5 + FITA6*DLT6 
               IF ( I .GT. J ) THEN
                  JI = (I-1) * NS + J
                  AIJ(JI) = AIJ(IJ)
               ENDIF
            ENDDO
         ENDDO
      ELSEIF ( IFLAG .EQ. 7 ) THEN
         DO J = 1, NS
            DO I = J, NS
               IJ = (J-1) * NS + I
C..............
               FITA0 = FITA(1,I,J)
               FITA1 = FITA(2,I,J)
               FITA2 = FITA(3,I,J)
               FITA3 = FITA(4,I,J)
               FITA4 = FITA(5,I,J)
               FITA5 = FITA(6,I,J)
               FITA6 = FITA(7,I,J)
C..............
               FITB0 = FITB(1,I,J)
               FITB1 = FITB(2,I,J)
               FITB2 = FITB(3,I,J)
               FITB3 = FITB(4,I,J)
               FITB4 = FITB(5,I,J)
               FITB5 = FITB(6,I,J)
               FITB6 = FITB(7,I,J)
C..............
               FITC0 = FITC(1,I,J)
               FITC1 = FITC(2,I,J)
               FITC2 = FITC(3,I,J)
               FITC3 = FITC(4,I,J)
               FITC4 = FITC(5,I,J)
               FITC5 = FITC(6,I,J)
               FITC6 = FITC(7,I,J)
C..............
               AIJ(IJ) = FITA0 + FITA1*DLT1 + FITA2*DLT2 + FITA3*DLT3
     &                         + FITA4*DLT4 + FITA5*DLT5 + FITA6*DLT6 
C
               BIJ(IJ) = FITB0 + FITB1*DLT1 + FITB2*DLT2 + FITB3*DLT3
     &                         + FITB4*DLT4 + FITB5*DLT5 + FITB6*DLT6 
C
               CIJ(IJ) = FITC0 + FITC1*DLT1 + FITC2*DLT2 + FITC3*DLT3
     &                         + FITC4*DLT4 + FITC5*DLT5 + FITC6*DLT6 
C..............
               IF ( I .GT. J ) THEN
                  JI = (I-1) * NS + J
                  AIJ(JI) = AIJ(IJ)
                  BIJ(JI) = BIJ(IJ)
                  CIJ(JI) = CIJ(IJ)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Test for IFLAG = 3
C-----------------------------------------------------------------------
      IF (IFLAG .EQ. 3) RETURN
C-----------------------------------------------------------------------
C         COMPUTE PARKER CORRECTION FOR ZROT
C         AND ALSO THE ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
C-----------------------------------------------------------------------
      DO 400 K = 1, NS
         IF (LIN(K) .EQ. 0) THEN
            CROT    = 0.0D0
            CINT(K) = 0.0D0
         ELSEIF (LIN(K) .EQ. 1) THEN
            WRU = WT(K) / RU
            CROT    = 1.0D0
            CINT(K) = CPMS(K) * WRU - 2.5D0
         ELSEIF (LIN(K) .EQ. 2) THEN
            WRU = WT(K) / RU
            CROT    = 1.5D0
            CINT(K) = CPMS(K) * WRU - 2.5D0
         ENDIF
C........
         DR   = EPS(K) / 298.0D0
         SQDR = DSQRT(DR)
         DR32 = SQDR*DR
         AAAA = (1.0D0 + PI32O2*SQDR + P2O4P2*DR + PI32*DR32) 
C........
         DD   = EPS(K) / T
         SQDD = DSQRT(DD)
         DD32 = SQDD*DD
         BBBB = (1.0D0 + PI32O2*SQDD + P2O4P2*DD + PI32*DD32) 
C........
         XI = ( AAAA / BBBB ) * MAX(1.0D0, ZROT(K))
         CXI(K) = CROT / ( XI * PI )
  400 CONTINUE
C-----------------------------------------------------------------------
C     Test for IFLAG = 4 or 5
C-----------------------------------------------------------------------
      IF (IFLAG .LE. 5) RETURN
C-----------------------------------------------------------------------
C        EVALUATE THE BINARY INTERNAL DIFFUSION COEFFICIENTS
C                     D_int,ii
C-----------------------------------------------------------------------
      DO I = 1, NS
         II = (I-1) * NS + I
         IF ( IFLAG .EQ. 7 ) THEN
            AAA = AIJ(II)
         ELSE
            AAA = CTAIJ(II)
         ENDIF
         IBIN = NS*(I-1) - (I*(I-1))/2 + I
         BIN(IBIN) = 5.0D0 * PATMOS * WT(I) / 
     &             ( 6.0D0 * RU * T * AAA * ETA(I) )
         IF ( DIP(I) .NE. 0.0D0 ) THEN
            BIN(IBIN) = BIN(IBIN) * ( 1.0D0 + 2.985D3 / (T*SQRT(T)) )
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSCOFD ( COF, NS, DLT, DLT2, DLT3, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COF(4,NS,*), BIN(*)
C-----------------------------------------------------------------------
C     This subroutine returns the array BIN formed by the
C     reciprocals of the binary diffusion coefficients at
C     atmospheric pressure.
C-----------------------------------------------------------------------
      DO 120 K=1, NS
         KP1 = K + 1
         DO 110 L=KP1, NS
            COF1 = COF(1,K,L)
            COF2 = COF(2,K,L)
            COF3 = COF(3,K,L)
            COF4 = COF(4,K,L)
            DLK = DEXP( -(COF1 + COF2*DLT + COF3*DLT2
     1                       + COF4*DLT3) )
            KL = NS*(K-1) - (K*(K-1))/2 + L
            BIN(KL) = DLK
110      CONTINUE
120   CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSCOFE ( COF, NS, DLT, DLT2, DLT3, ETA, ETALG )
C-----------------------------------------------------------------------
C     This subroutine returns the arrays ETA and ETALG formed by the
C     pure species shear viscosities and their logs'.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COF(4,NS), ETA(NS), ETALG(NS)
C-----------------------------------------------------------------------
      DO 120 K=1, NS
            COF1 = COF(1,K)
            COF2 = COF(2,K)
            COF3 = COF(3,K)
            COF4 = COF(4,K)
            ETALG(K) = COF1 + COF2*DLT + COF3*DLT2 + COF4*DLT3 
            ETA(K)   = DEXP( ETALG(K) )
120   CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSBIN (T, WEG, BIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        T           temperature
C        WEG         double precision work array for EGLIB
C
C     Output
C     ------
C        BIN(NS,NS)  binary diffusion coefficients at atmospheric 
C                    pressure for species pairs (i,j) with i<>j.
C
C-----------------------------------------------------------------------
      INCLUDE 'eg.cmn'
      DIMENSION WEG(*), BIN(NS,NS)
C-----------------------------------------------------------------------
      IND = 0
      DO I = 1, NS
         IND = IND + 1
         DO J = I+1, NS
            BIN(I,J) = 1.0D0 / WEG ( IBIN + IND )
            BIN(J,I) = BIN(I,J)
            IND = IND + 1
         ENDDO
      ENDDO
C-----------------------------------------------------------------------      
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSD1 (PRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS,NS)  flux diffusion matrix
C
C     
C     Two standard iterations are performed on matrix L_[00]
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSD1  ( NS, PRES, WEG(IXTR), WEG(IYTR), Y, WEG(IEGWT), 
     &               WW, WEG(IEGPA), D, WEG(IBIN), WEG(IG), 
     &               WEG(IDMI), WEG(ITEMP), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSD1 ( NS, PRES, XTR, YTR, Y, WT, WW, PATMOS, D,
     &                    BIN, G, DMI, TEMP, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      ITERMX = 2
      CALL EGSSI1 ( NS, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSD2 (PRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L_[00] 
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSD2  ( NS, PRES, WEG(IXTR), WEG(IYTR), Y, WEG(IEGWT), 
     &               WW, WEG(IEGPA), D, WEG(IBIN), WEG(IG), 
     &               WEG(IAN), WEG(ITEMP), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSD2 ( NS, PRES, XTR, YTR, Y, WT, WW, PATMOS, D,
     &                    BIN, G, AN, TEMP, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
C     Form the linear system
C-----------------------------------------------------------------------
      CALL EGSDDEC ( NS, XTR, YTR, G, BIN, TEMP )
C-----------------------------------------------------------------------
C     Form the matrix D
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NS, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C     Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSDR1 (T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     Two standard iterations are performed on matrix L_[00]
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C        D is also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSDR1 ( NS, T, WEG(IXTR), WEG(IYTR), Y, WEG(IEGWT), 
     &               WW, WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IG), WEG(IDMI), WEG(ITEMP), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSDR1 ( NS, TEMPER, XTR, YTR, Y, WT, WW, RU, PATMOS, 
     &                     D, BIN, G, DMI, TEMP, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      ITERMX = 2
      CALL EGSSI1 ( NS, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSDR2 (T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     We form the Cholski decomposition of matrix L_[00]
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C        D is also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSDR2 ( NS, T, WEG(IXTR), WEG(IYTR), Y, WEG(IEGWT), 
     &               WW, WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IG), WEG(IAN), WEG(ITEMP), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSDR2 ( NS, TEMPER, XTR, YTR, Y, WT, WW, RU, PATMOS, 
     &                     D, BIN, G, AN, TEMP, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
C     Form the linear system
C-----------------------------------------------------------------------
      CALL EGSDDEC ( NS, XTR, YTR, G, BIN, TEMP )
C-----------------------------------------------------------------------
C     Form the matrix D
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NS, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSE1 ( ALPHA, T, X, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 1
C        ITLS  = 0
C
C     Input
C     -----
C        ALPHA     parameter for the averaging formula
C        T         temperature
C        X(NS)     species mole fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        ETA       shear viscosity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate ETA.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSE1 ( NS, X, ETA, WEG(IETALG), ALPHA ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSE1 ( NS, X, ETAMA, ETALG, ALPHA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NS), ETALG(NS)
C-----------------------------------------------------------------------
      SUM = 0.0D0
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            SUM = SUM + X(I) * ETALG(I)
         ENDDO
         ETAMA = DEXP(SUM)
      ELSE
         DO I = 1, NS
            SUM = SUM + X(I) * DEXP ( ALPHA*ETALG(I) )
         ENDDO
         ETAMA = SUM ** (1.0D0/ALPHA)
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSE2 ( T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        ETA       shear viscosity 
C
C     
C     One CG iteration is performed on the matrix H
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSE2 ( NS, T, WEG(IXTR), WEG(IEGWT), ETA, WEG(IEGRU), 
     &              WEG(IEGPA),
     &              WEG(IBETA), WEG(IAN), WEG(IZN), WEG(IRN), 
     &              WEG(IG), WEG(IDMI), WEG(ITEMP),
     &              WEG(ICTAIJ), WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSE2 ( NS, T, X, WT, ETACG, RU, PATMOS, 
     &                    BETA, AN, ZN, RN, G, DMI, TEMP, 
     &                    AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMH ( NS, T, X, WT, RU, PATMOS, BETA, G, AIJ, ETA, BIN )
C-----------------------------------------------------------------------
      ITERMX = 1
      CALL DCOPY  ( NS, BETA, 1, RN, 1 )
      CALL EGSCG1 ( NS, G, DMI, AN, ZN, RN, TEMP, ITERMX )
      CALL EGSDOT ( NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSE3 ( T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 3
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        ETA       shear viscosity 
C
C     
C     One CG iteration is performed on the matrix H
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSE3 ( NS, T, WEG(IXTR), WEG(IEGWT), ETA, WEG(IEGRU), 
     &              WEG(IEGPA),
     &              WEG(IBETA), WEG(IAN), WEG(IZN), WEG(IRN), 
     &              WEG(IG), WEG(IDMI), WEG(ITEMP),
     &              WEG(IAIJ), WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSE3 ( NS, T, X, WT, ETACG, RU, PATMOS, 
     &                    BETA, AN, ZN, RN, G, DMI, TEMP, 
     &                    AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMH ( NS, T, X, WT, RU, PATMOS, BETA, G, AIJ, ETA, BIN )
C-----------------------------------------------------------------------
      ITERMX = 1
      CALL DCOPY  ( NS, BETA, 1, RN, 1 )
      CALL EGSCG1 ( NS, G, DMI, AN, ZN, RN, TEMP, ITERMX )
      CALL EGSDOT ( NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSE4 ( T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 3
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        ETA       shear viscosity
C
C     
C     We form the Choleski decomposition of matrix H 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSE4 ( NS, T, WEG(IXTR), WEG(IEGWT), ETA, WEG(IEGRU), 
     &              WEG(IEGPA),
     &              WEG(IBETA), WEG(IAN), WEG(IZN), WEG(IRN), 
     &              WEG(IG), WEG(IDMI), WEG(ITEMP),
     &              WEG(IAIJ), WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSE4 ( NS, T, X, WT, ETACG, RU, PATMOS, 
     &                    BETA, AN, ZN, RN, G, DMI, TEMP, 
     &                    AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMH ( NS, T, X, WT, RU, PATMOS, BETA, G, AIJ, ETA, BIN )
C-----------------------------------------------------------------------
      CALL EGSDEC ( NS, G, TEMP, IER )
      CALL DCOPY  ( NS, BETA, 1, AN, 1 )
      CALL EGSSOL ( NS, G, AN )
      CALL EGSDOT ( NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK1 ( ALPHA, T, X, WEG, VOLVIS )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 0
C
C     Input
C     -----
C        ALPHA     parameter for the averaging formula
C        T         temperature
C        X(NS)     species mole fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VOLVIS    volume viscosity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate VOLVIS.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK1 ( NS, X, VOLVIS, WEG(IETALG), WEG(ICXI),
     &              WEG(ICINT), ALPHA ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK1 ( NS, X, VOLVIS, ETALG, CXI, CINT, ALPHA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NS), ETALG(NS), CXI(NS), CINT(NS)
C-----------------------------------------------------------------------
      SXP = 0.0D0
      DO I = 1, NS
         IF ( CXI(I) .NE. 0.0D0 ) SXP = SXP + X(I)
      ENDDO
      SXP = 1.0D0 / SXP
      SUM = 0.0D0
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            IF ( CXI(I) .NE. 0.0D0 ) THEN
               CCC = CINT(I) / ( 1.5D0 + CINT(I) )
               VVV = ETALG(I) + DLOG ( 0.25D0*CCC*CCC/CXI(I) )
               SUM = SUM + SXP * X(I) * VVV
            ENDIF
         ENDDO
         VOLVIS = DEXP(SUM)
      ELSE
         DO I = 1, NS
            IF ( CXI(I) .NE. 0.0D0 ) THEN
               CCC = CINT(I) / ( 1.5D0 + CINT(I) )
               VVV = ETALG(I) + DLOG ( 0.25D0*CCC*CCC/CXI(I) )
               SUM = SUM + SXP * X(I) * DEXP ( ALPHA*VVV )
            ENDIF
         ENDDO
         VOLVIS = SUM ** (1.0D0/ALPHA)
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK2 (T, Y, WEG, VV01)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VV01      volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity \kappa_[01]
C
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK2 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV01, 
     &              WEG(IBIN), WEG(ICTAIJ), WEG(ICINT),
     &              WEG(IETA), WEG(ICXI), WEG(IG), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK2 ( NS, TEMPER, XTR, WT, RU, PATMOS, VV01, 
     &                    BIN, AIJ, CINT, ETA, CXI, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), CINT(NS), ETA(NS), XTR(NS), 
     &          WT(NS), CXI(*)
C-----------------------------------------------------------------------
      CALL EGSEMK01 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, CINT, ETA, CXI, G, BETA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      VV01 = 0.0D0
      DO I = 1, NS
         III  = NS*(I-1) - (I*(I-1))/2 + I
         VV01 = VV01 + BETA(I) * BETA(I) / G(III)
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK3 (T, Y, WEG, VV01)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VV01      volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity \kappa_[01]
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK3 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV01, 
     &              WEG(IBIN), WEG(IAIJ), WEG(ICINT), WEG(IETA), 
     &              WEG(ICXI), WEG(IG), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK3 ( NS, TEMPER, XTR, WT, RU, PATMOS, VV01, 
     &                    BIN, AIJ, CINT, ETA, CXI, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), CINT(NS), ETA(NS), XTR(NS), 
     &          WT(NS), CXI(*)
C-----------------------------------------------------------------------
      CALL EGSEMK01 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, CINT, ETA, CXI, G, BETA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      VV01 = 0.0D0
      DO I = 1, NS
         III  = NS*(I-1) - (I*(I-1))/2 + I
         VV01 = VV01 + BETA(I) * BETA(I) / G(III)
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK4 (T, Y, WEG, VV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VV        volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity 
C                  VV = \kappa_[s]^1 + \kappa_[01]
C     where \kappa_[s]^1 is obtained after one standard iteration
C     on the Schur complement K_[s].
C
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK4 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV, 
     &              WEG(IBIN), WEG(ICTAIJ), WEG(ICINT),
     &              WEG(IETA),  WEG(ICXI), WEG(IG), WEG(IBETA),
     &              WEG(IAN), WEG(ITEMP) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK4 ( NS, TEMPER, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, CINT, ETA, CXI, G, BETA, AN, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMK ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, CINT, ETA, CXI, G, BETA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGSEVK ( NS, VV, G, AN, BETA, TEMP )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK5 (T, Y, WEG, VV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VV        volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity 
C                  VV = \kappa_[s]^1 + \kappa_[01]
C     where \kappa_[s]^1 is obtained after one standard iteration
C     on the Schur complement K_[s].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK5 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV, 
     &              WEG(IBIN), WEG(IAIJ), WEG(ICINT), WEG(IETA), 
     &              WEG(ICXI), WEG(IG), WEG(IBETA),
     &              WEG(IAN), WEG(ITEMP) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK5 ( NS, TEMPER, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, CINT, ETA, CXI, G, BETA, AN, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMK ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, CINT, ETA, CXI, G, BETA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGSEVK ( NS, VV, G, AN, BETA, TEMP )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSK6 ( T, Y, WEG, VV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        VV        volume viscosity
C
C     
C     Direct inversion of the matrix K is performed by 
C     forming its LDL^t decomposition.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSK6 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV, 
     &              WEG(IBIN), WEG(IAIJ), WEG(ICINT), WEG(IETA), 
     &              WEG(ICXI), WEG(IG), WEG(IBETA),
     &              WEG(IAN), WEG(ITEMP) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSK6 ( NS, TEMPER, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, CINT, ETA, CXI, G, BETA, AN, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), CINT(NS), XTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEMK ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, CINT, ETA, CXI, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         II2 = NG*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         G(II1) = G(II1) + 1.5D0   * 1.5D0   * XTR(I) * XTR(I)
         G(II2) = G(II2) + 1.5D0   * CINT(I) * XTR(I) * XTR(I)
         G(II3) = G(II3) + CINT(I) * CINT(I) * XTR(I) * XTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            IJ2 = NG*(I-1) - (I*(I-1))/2 + J + NS
            JI2 = NG*(J-1) - (J*(J-1))/2 + I + NS
            IJ3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            G(IJ1) = G(IJ1) + 1.5D0   * 1.5D0   * XTR(I) * XTR(J)
            G(IJ2) = G(IJ2) + 1.5D0   * CINT(J) * XTR(I) * XTR(J)
            G(JI2) = G(JI2) + 1.5D0   * CINT(I) * XTR(I) * XTR(J)
            G(IJ3) = G(IJ3) + CINT(I) * CINT(J) * XTR(I) * XTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
      CALL EGSDOT ( NG, VV, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSL1 ( ALPHA, T, X, WEG, CON )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 1
C        ITLS  = 0
C
C     Input
C     -----
C        ALPHA     parameter for the averaging formula
C        T         temperature
C        X(NS)     species mole fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        CON       thermal conductivity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate CON.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSL1 ( NS, X, T, CON, WEG(IEGCFL), ALPHA ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSL1 ( NS, X, T, CONMA, COF, ALPHA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NS), COF(4,NS)
C-----------------------------------------------------------------------
      DLT  = DLOG(T)
      DLT2 = DLT * DLT
      DLT3 = DLT * DLT2
      SUM = 0.0D0
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            COF1 = COF(1,I)
            COF2 = COF(2,I)
            COF3 = COF(3,I)
            COF4 = COF(4,I)
            CONLG = COF1 + COF2*DLT + COF3*DLT2 + COF4*DLT3 
            SUM = SUM + X(I) * CONLG
         ENDDO
         CONMA = DEXP(SUM)
      ELSE
         DO I = 1, NS
            COF1 = COF(1,I)
            COF2 = COF(2,I)
            COF3 = COF(3,I)
            COF4 = COF(4,I)
            CONLG = COF1 + COF2*DLT + COF3*DLT2 + COF4*DLT3 
            SUM = SUM + X(I) * DEXP ( ALPHA*CONLG )
         ENDDO
         CONMA = SUM ** (1.0D0/ALPHA)
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSL2 (T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 6
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC        thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda_[e].
C
C     The collision integrals A_ij, B_ij and C_ij are assumed to be
C     temperature independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSL2 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), TC, WEG(IBIN), 
     &              WEG(ICTAIJ), WEG(ICTBIJ), WEG(ICTCIJ), 
     &              WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &              IWEG(IEGLIN),
     &              WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &              WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSL2 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, 
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMAE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 1
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG1 ( NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSL3 (T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC        thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda_[e].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSL3 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), TC, 
     &              WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &              WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &              IWEG(IEGLIN),
     &              WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &              WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSL3 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, 
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMAE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 1
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG1 ( NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSL4 (T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 6
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC        thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda.
C
C     The collision integrals A_ij, B_ij and C_ij are assumed to be
C     temperature independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSL4 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), TC, WEG(IBIN),
     &              WEG(ICTAIJ), WEG(ICTBIJ), WEG(ICTCIJ), 
     &              WEG(ICINT), WEG(IETA), WEG(ICXI),
     &              IWEG(IEGLIN), WEG(IG), WEG(IDMI), WEG(IAN),
     &              WEG(IZN),  WEG(IRN), WEG(ITEMP), WEG(IBETA))
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSL4 ( NS, TEMPER, XTR, WT, RU, PATMOS, TC, 
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMA ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 1
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSL5 (T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC        thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSL5 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &            WEG(IEGPA), TC, 
     &            WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &            WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &            WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &            WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSL5 ( NS, TEMPER, XTR, WT, RU, PATMOS, TC, 
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMA ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 1
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLCT1 (T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC         thermal conductivity 
C        CHIT(NS)   rescaled thermal diffusion ratios
C
C     
C     Three CG iterations are performed on the matrix \Lambda_[e].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLCT1 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &                WEG(IEGPA), TC, CHIT,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &                WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &                IWEG(IEGLIN),
     &                WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &                WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLCT1 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMAE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG1 ( NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL  ( NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVCT ( NS, CHIT, AN, XTR, BIN, CIJ, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLCT2 ( T, Y, WEG, IWEG, TC, CHIT )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC         thermal conducitvity
C        CHIT(NS)   rescaled thermal diffusion ratios
C
C     
C     We form the Choleski decomposition of matrix \Lambda_[e] 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLCT2 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &                WEG(IEGPA), TC, CHIT,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &                WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &                IWEG(IEGLIN),
     &                WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &                WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLCT2 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMAE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL  ( NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVCT ( NS, CHIT, AN, XTR, BIN, CIJ, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLCT3 (T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC         thermal conductivity 
C        CHIT(NS)   rescaled thermal diffusion ratios
C
C     
C     Three CG iterations are performed on the matrix \Lambda.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLCT3 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &                WEG(IEGPA), TC, CHIT,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &                WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &                IWEG(IEGLIN),WEG(IG), WEG(IDMI), WEG(IAN), 
     &                WEG(IZN), WEG(IRN), WEG(ITEMP), 
     &                WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLCT3 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMA ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL  ( NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVCT ( NS, CHIT, AN, XTR, BIN, CIJ, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLCT4 (T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        TC         thermal conducitvity
C        CHIT(NS)   rescaled thermal diffusion ratios
C
C     
C     We form the Choleski decomposition of matrix \Lambda 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLCT4 ( NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &             WEG(IEGPA), TC, CHIT,
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLCT4 ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMA ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
      NG = 2 * NS
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL  ( NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVCT ( NS, CHIT, AN, XTR, BIN, CIJ, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD1 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) thermal diffusion vector
C        D(NS,NS)  flux diffusion matrix
C
C
C     Two CG iterations are performed on the matrix L_[e] in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[e] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD1 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &        WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &        PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD1 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 2
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C.....
      ITERMX = 2
      CALL EGSSI2 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD2 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) thermal diffusion ratios
C        D(NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L_[e] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD2 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &        WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &        PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD2 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD3 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) thermal diffusion vector
C        D(NS,NS)  FIRST ORDER flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD3 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &             WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &             PTC, THETA, D,
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &             IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD3 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C.....
      ITERMX = 2
      CALL EGSSI1 ( NS, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD4 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) thermal diffusion ratios
C        D(NS,NS)  FIRST ORDER flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrices L and L_[00] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD4 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &                WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &                PTC, THETA, D,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &                WEG(ICIJ), WEG(ICINT), WEG(IETA), 
     &                WEG(ICXI), IWEG(IEGLIN), WEG(IG), 
     &                WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &                WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD4 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD5 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) thermal diffusion vector
C        D(NS,NS)  flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD5 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &             WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &             PTC, THETA, D,
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &             IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD5 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C.....
      ITERMX = 2
      CALL EGSSI3 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTD6 (PRES, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) thermal diffusion ratios
C        D(NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTD6 ( NS, PRES, T, WEG(IXTR), WEG(IYTR), Y,
     &                WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &                PTC, THETA, D,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &                WEG(ICIJ), WEG(ICINT), WEG(IETA), 
     &                WEG(ICXI), IWEG(IEGLIN), WEG(IG), 
     &                WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &                WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTD6 ( NS, PRES, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR1 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) rho * thermal diffusion vector
C        D(NS,NS)  rho * flux diffusion matrix
C
C
C     Two CG iterations are performed on the matrix L_[e] in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[e] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR1 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &        WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &        PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR1 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 2
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C.....
      ITERMX = 2
      CALL EGSSI2 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR2 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) rho * thermal diffusion ratios
C        D(NS,NS)  rho * flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L_[e] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR2 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &        WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &        PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &        IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR2 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEMLE ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR3 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) rho * thermal diffusion vector
C        D(NS,NS)  rho * FIRST ORDER flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR3 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &             WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &             PTC, THETA, D,
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &             IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR3 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C.....
      ITERMX = 2
      CALL EGSSI1 ( NS, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR4 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) rho * thermal diffusion ratios
C        D(NS,NS)  rho * FIRST ORDER flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrices L and L_[00] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR4 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &                WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &                PTC, THETA, D,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &                WEG(ICIJ), WEG(ICINT), WEG(IETA), 
     &                WEG(ICXI), IWEG(IEGLIN), WEG(IG), 
     &                WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &                WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR4 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR5 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       partial thermal conductivity
C        THETA(NS) rho * thermal diffusion vector
C        D(NS,NS)  rho * flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR5 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &             WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &             PTC, THETA, D,
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &             IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR5 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C.....
      ITERMX = 2
      CALL EGSSI3 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW, Y, 
     &             D, ITERMX, AUX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSLTDR6 (T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        PTC       thermal conducitvity
C        THETA(NS) rho * thermal diffusion ratios
C        D(NS,NS)  rho * flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSLTDR6 ( NS, T, WEG(IXTR), WEG(IYTR), Y,
     &                WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &                PTC, THETA, D,
     &                WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &                WEG(ICIJ), WEG(ICINT), WEG(IETA), 
     &                WEG(ICXI), IWEG(IEGLIN), WEG(IG), 
     &                WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &                WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSLTDR6 ( NS, TEMPER, XTR, YTR, Y, WT, WW,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &           LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      CALL DCOPY  ( NG, BETA, 1, AN, 1 )
      CALL EGSSOL ( NG, G, AN )
C-----------------------------------------------------------------------
      CALL EGSEVL ( NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
      CALL EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSTD1 (PRES, T, Y, WW, WEG, IWEG, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        THETA(NS) thermal diffusion vector
C        D(NS,NS)  flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSTD1 ( NS, PRES, T, Y,  WEG(IXTR), WEG(IYTR), 
     &               WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &               THETA, D,
     &               WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &               WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &               IWEG(IEGLIN), WEG(IG), 
     &               WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &               WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSTD1 ( NS, PRES, TEMPER, Y, XTR, YTR, WT, WW, 
     &           RU, PATMOS, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C       Evaluate theta 
C-----------------------------------------------------------------------
      ITERMX = 3
      NG = 3 * NS
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
C     Evaluate the diffusion matrix
C-----------------------------------------------------------------------
      CALL EGSEVDI1 ( NS, D, WT, WW, TEMP, XTR, YTR, BIN, AUX )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSTDR1 (T, Y, WW, WEG, IWEG, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        IWEG      integer work array for EGLIB
C
C     Output
C     ------
C        THETA(NS) rho * thermal diffusion vector
C        D(NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSTDR1 ( NS, T, Y,  WEG(IXTR), WEG(IYTR), 
     &       WEG(IEGWT), WW, WEG(IEGRU), WEG(IEGPA), 
     &       THETA, D,
     &       WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &       WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &       IWEG(IEGLIN),
     &       WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &       WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSTDR1 ( NS, TEMPER, Y, XTR, YTR, WT, WW, 
     &           RU, PATMOS, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGSEML ( NS, TEMPER, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &              LIN, G, BETA )
C-----------------------------------------------------------------------
C       Evaluate theta 
C-----------------------------------------------------------------------
      ITERMX = 3
      NG = 3 * NS
      CALL DCOPY  ( NG, BETA, 1, RN, 1 )
      CALL EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
C-----------------------------------------------------------------------
      CALL EGSEVT ( NS, THETA, AN, Y )
C-----------------------------------------------------------------------
C     Evaluate the diffusion matrix
C-----------------------------------------------------------------------
      CALL EGSEVDI1 ( NS, D, WT, WW, TEMP, XTR, YTR, BIN, AUX )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGSPRD ( NS, D, TEMP, Y )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      AA = PATMOS * WW / ( RU * TEMPER )
      CALL DSCAL ( NS, AA, THETA, 1 )
      CALL DSCAL ( NS*NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSV1 (PRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        PRES      pressure
C        T         temperature
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS)     flux diffusion coefficients
C
C     
C     The array D corresponds to the diagonal of the matrix
C     ~D_[00]^1 before projection.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSV1 ( NS, PRES, WEG(IXTR), WEG(IYTR), WEG(IEGWT), 
     &              WW, WEG(IEGPA), D, WEG(IBIN), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSV1 ( NS, PRES, XTR, YTR, WT, WW, PATMOS, 
     &                    D, BIN, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTR(NS), YTR(NS), WT(NS), D(NS), BIN(*), AUX(NS)
C-----------------------------------------------------------------------
      CALL EGZERO ( NS, D )
      DO I = 1, NS
         IP1 = I + 1
         DO J = IP1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            AAI = XTR(I) * BIN(IJ)
            AAJ = XTR(J) * BIN(IJ)
            D(I) = D(I) + AAJ
            D(J) = D(J) + AAI
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RWW = 1.0D0 / WW
      DO I = 1, NS
         D(I) = WT(I) * RWW * AUX(I) / D(I)
      ENDDO
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      AA = PATMOS / PRES
      CALL DSCAL ( NS, AA, D, 1 )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSVR1 (T, Y, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C
C     Output
C     ------
C        D(NS)     rho * flux diffusion coefficients
C
C        where rho is the density.
C
C     
C     The array D corresponds to the diagonal of the matrix
C     rho * ~D_[00]^1 before projection.
C     
C
C     Note
C     ----
C        The coefficients D are PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSVR1 ( NS, T, WEG(IXTR), WEG(IYTR), WEG(IEGWT), 
     &               WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSVR1 ( NS, TEMP, XTR, YTR, WT, RU, PATMOS, 
     &                     D, BIN, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTR(NS), YTR(NS), WT(NS), D(NS), BIN(*), AUX(NS)
C-----------------------------------------------------------------------
      CALL EGZERO ( NS, D )
      DO I = 1, NS
         IP1 = I + 1
         DO J = IP1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            AAI = XTR(I) * BIN(IJ)
            AAJ = XTR(J) * BIN(IJ)
            D(I) = D(I) + AAJ
            D(J) = D(J) + AAI
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      FAC = PATMOS / ( RU * TEMP )
      DO I = 1, NS
         D(I) = WT(I) * FAC * AUX(I) / D(I)
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSYV (PRES, Y, WW, WEG, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        PRES      pressure
C        Y(NS)     species mass fractions
C        WW        mean molecular weight
C        WEG       double precision work array for EGLIB
C        F(NS)     species diffusion driving forces
C
C     Output
C     ------
C        F(NS)     species flux diffusion velocities
C                  F_i <-- - sum_j \widetilde D_{ij} F_j  
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSYV ( NS, PRES, WEG(IXTR), WEG(IYTR), WEG(IEGWT), WW, 
     &              WEG(IEGPA), WEG(IG), WEG(ITEMP), WEG(IBIN), F )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSYV ( NS, PRES, XTR, YTR, WT, WW, PATMOS, G, 
     &                    TEMP, BIN, F )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION F(NS), WT(NS), XTR(NS)
C-----------------------------------------------------------------------
      CALL EGSDDEC ( NS, XTR, YTR, G, BIN, TEMP )
      CALL EGSSOL  ( NS, G, F )
      AA = - PATMOS / PRES / WW
      DO I = 1, NS
         F(I) = F(I) * WT(I) * XTR(I) * AA
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSRYV (T, Y, WEG, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        T         temperature
C        Y(NS)     species mass fractions
C        WEG       double precision work array for EGLIB
C        F(NS)     species diffusion driving forces
C
C     Output
C     ------
C        F(NS)     rescaled species flux diffusion velocities
C                  F_i <-- - rho * sum_j \widetilde D_{ij} F_j  
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGSRYV ( NS, WEG(IXTR), WEG(IYTR), WEG(IEGWT),
     &               WEG(IEGPA), WEG(IEGRU), T, WEG(IG), WEG(ITEMP), 
     &               WEG(IBIN), F )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGSRYV ( NS, XTR, YTR, WT, PATMOS, RU, TEMPER, 
     &                     G, TEMP, BIN, F )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION F(NS), WT(NS), XTR(NS)
C-----------------------------------------------------------------------
      CALL EGSDDEC ( NS, XTR, YTR, G, BIN, TEMP )
      CALL EGSSOL  ( NS, G, F )
      AA = - PATMOS / RU / TEMPER
      DO I = 1, NS
         F(I) = F(I) * WT(I) * XTR(I) * AA
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSDOT ( NG, DOT, AN, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION AN(NG), BETA(NG)
C-----------------------------------------------------------------------
      DOT = 0.0D0
      DO I = 1, NG
         DOT = DOT + AN(I) * BETA(I)
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVCT ( NS, CHIT, AN, XTR, BIN, CIJ, WT )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION CHIT(NS), AN(NS), XTR(NS), BIN(*), CIJ(NS,NS), 
     &          WT(NS) 
C-----------------------------------------------------------------------
      CALL EGZERO ( NS, CHIT )
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            TERM = 0.5D0 * BIN(IJ) * (6.0D0*CIJ(I,J) - 5.0D0 )
     &             * ( WT(I)*AN(J) - WT(J)*AN(I) ) / ( WT(I)+WT(J) )
            CHIT(I) = CHIT(I) + XTR(J) * TERM
            CHIT(J) = CHIT(J) - XTR(I) * TERM
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVD ( NS, NG, AN, YTR, G, XTR, WT, WW, D )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION AN(NG), YTR(NS), G(*), XTR(NS), WT(NS), D(NS,NS)
C-----------------------------------------------------------------------
      RWW = 1.0D0 / WW
      DO I = 1, NS
         CALL EGZERO ( NG, AN )
         CALL DCOPY ( NS, YTR, 1, AN, 1 )
         CALL DSCAL ( NS, -1.0D0, AN, 1 )
         AN(I) = AN(I) + 1.0D0
         CALL EGSSOL ( NG, G, AN )
         AAA = XTR(I) * WT(I) * RWW 
         DO J = 1, NS
            D(I,J) = AAA * AN(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVDI1 ( NS, D, WT, WW, TEMP, XTR, YTR, BIN, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION WT(NS), TEMP(NS), XTR(NS), YTR(NS), BIN(*), AUX(NS)
      DIMENSION D(NS,NS)
C-----------------------------------------------------------------------
      CALL EGZERO ( NS, TEMP )
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            TEMP(I) = TEMP(I) + XTR(J) * BIN(IJ)
            TEMP(J) = TEMP(J) + XTR(I) * BIN(IJ)
         ENDDO
      ENDDO
      DO I = 1, NS
         TEMP(I) = AUX(I) / TEMP(I)
      ENDDO
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            D(I,J) = WT(I) / WW * TEMP(I) * TEMP(J) * XTR(I) * BIN(IJ)
            D(J,I) = WT(J) / WW * TEMP(I) * TEMP(J) * XTR(J) * BIN(IJ)
         ENDDO
         D(I,I) = WT(I) / WW * TEMP(I) * ( 1.0D0 + YTR(I) )
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVK ( NS, VV, G, AN, BETA, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*), AN(*), TEMP(*)
C-----------------------------------------------------------------------
      NS2  = 2 * NS
      DO I = NS+1, NS2
         III  = NS2*(I-1) - (I*(I-1))/2 + I
         AN(I-NS) = BETA(I) / G(III)
      ENDDO
      DO I = 1, NS
         TEMP(I) = BETA(I)
         DO J = 1, NS
            IJ = NS2*(I-1) - (I*(I-1))/2 + J + NS
            TEMP(I) = TEMP(I) - G(IJ) * AN(J)
         ENDDO
      ENDDO
      VV01 = 0.0D0
      VVS  = 0.0D0
      DO I = 1, NS
         II = NS2*(I-1) - (I*(I-1))/2 + I 
         VV01 = VV01 + AN(I) * BETA(I+NS)
         VVS  = VVS  + TEMP(I) * TEMP(I) / G(II)
      ENDDO
      VV = VV01 + VVS
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVL ( NG, TC, AN, BETA, PATMOS, TEMPER )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION AN(NG), BETA(NG)
C-----------------------------------------------------------------------
      TC = 0.0D0
      DO I = 1, NG
         TC = TC + AN(I) * BETA(I)
      ENDDO
      TC = TC * PATMOS / TEMPER
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEVT ( NS, THETA, AN, Y )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION THETA(NS), AN(NS), Y(NS)
C-----------------------------------------------------------------------
      SUM = 0.0D0
      SOM = 0.0D0
      DO I = 1, NS
         SUM = SUM + Y(I) * AN(I)
         SOM = SOM + Y(I) 
      ENDDO
      SUM = SUM / SOM
      DO I = 1, NS
         THETA(I) =  SUM - AN(I)
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSPRD ( NS, D, TEMP, Y )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Input D and Y
C
C     Output D <-- P D P where P = I - Y*U / sum Y
C-----------------------------------------------------------------------
      DIMENSION D(NS,NS), TEMP(NS), Y(NS)
C-----------------------------------------------------------------------
      SOM = 0.0D0
      DO I = 1, NS
         SOM = SOM + Y(I)
      ENDDO
      RSOM = 1.0D0 / SOM
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         AAA = Y(J)
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(I,J) * AAA
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = TEMP(I) * RSOM
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA
         ENDDO
      ENDDO
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(J,I) 
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = Y(I) * RSOM
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA * TEMP(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSAXS ( N, A, X, B )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(*), X(N), B(N)
C-----------------------------------------------------------------------
C                             B == A.X
C
C     The matrix A is stored in symmetric form, i.e., instead of
C     A(i,j) we store A(n*(i-1) - i*(i-1)/2 + j)
C-----------------------------------------------------------------------
      DO 10 I = 1, N
         III = N*(I-1) - (I*(I-1))/2 + I
         B(I) = A(III) * X(I)
10    CONTINUE
      DO 30 I = 1, N
         DO 20 J = I+1, N
            III = N*(I-1) - (I*(I-1))/2 + J
            B(I) = B(I) + A(III)*X(J)
            B(J) = B(J) + A(III)*X(I)
20       CONTINUE
30    CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSCG1 ( NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), DMI(NG),  AN(NG), ZN(NG),  RN(NG),
     &          TEMP(NG)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C-----------------------------------------------------------------------
      NITER = 0
      DO I = 1, NG
         III = NG*(I-1) - (I*(I-1))/2 + I
         AN(I) = 0.0D0
         ZN(I) = 0.0D0
         DMI(I) = 1.0D0 / G(III)
      ENDDO
      BETAN = 0.0D0
      AAA = 0.0D0
      DO I = 1, NG
         AAA = AAA + DMI(I) * RN(I)*RN(I)
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NG
         ZN(I) = DMI(I)*RN(I) + BETAN*ZN(I)
      ENDDO
      CALL EGSAXS(NG, G, ZN, TEMP)
      BBB = VDDOT (NG, ZN, 1, TEMP, 1)
      DO I = 1, NG
         AN(I) = AN(I) + AAA/BBB*ZN(I)
         RN(I) = RN(I) - AAA/BBB*TEMP(I)
      ENDDO
      CCC = 0.0D0
      DO I = 1, NG
         CCC = CCC + DMI(I) * RN(I)*RN(I)
      ENDDO
      BETAN = CCC/AAA
      AAA   = CCC
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSCG2 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), DMI(3,NS),  AN(NG), ZN(NG),  RN(NG),
     &          TEMP(NG)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C
C     The preconditioner matrix DMI is stored in compact form
C     DMI(n,i) with n referring to
C                   1  2
C            DMI =  2  3
C-----------------------------------------------------------------------
      NITER = 0
      DO I = 1, NG
         AN(I) = 0.0D0
         ZN(I) = 0.0D0
      ENDDO
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         II2 = NG*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         DETM  = 1.0D0 / ( (G(II1)*G(II3) - G(II2)*G(II2)) )
C........
         DMI(1,I) =   G(II3) * DETM
         DMI(2,I) = - G(II2) * DETM
         DMI(3,I) =   G(II1) * DETM
C........
c         DMI(I,I) = 1.0D0/G(I,I)
c         DMI(2,I) = 0.0D0
c         DMI(3,I) = 1.0D0/G(I+NS,I+NS)
      ENDDO
      BETAN = 0.0D0
      AAA = 0.0D0
      DO I = 1, NS
         AAA = AAA + RN(I)*(DMI(1,I)*RN(I)+DMI(2,I)*RN(I+NS))
     &             + RN(I+NS)*(DMI(2,I)*RN(I)+DMI(3,I)*RN(I+NS))
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NS
         ZN(I) = DMI(1,I)*RN(I) + DMI(2,I)*RN(I+NS) + BETAN*ZN(I)
         ZN(I+NS) = DMI(2,I)*RN(I) + 
     &              DMI(3,I)*RN(I+NS) + BETAN*ZN(I+NS)
      ENDDO
      CALL EGSAXS(NG, G, ZN, TEMP)
      BBB = VDDOT (NG, ZN, 1, TEMP, 1)
      DO I = 1, NG
         AN(I) = AN(I) + AAA/BBB*ZN(I)
         RN(I) = RN(I) - AAA/BBB*TEMP(I)
      ENDDO
      CCC = 0.0D0
      DO I = 1, NS
         CCC = CCC + RN(I)*(DMI(1,I)*RN(I)+DMI(2,I)*RN(I+NS))
     &             + RN(I+NS)*(DMI(2,I)*RN(I)+DMI(3,I)*RN(I+NS))
      ENDDO
      BETAN = CCC/AAA
      AAA   = CCC
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSCG3 ( NS, NG, G, DMI, AN, ZN, RN, TEMP, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), DMI(6,NS),  AN(NG), ZN(NG),  RN(NG),
     &          TEMP(NG)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C
C     The preconditioner matrix DMI is stored in compact form
C     DMI(n,i) with n referring to
C                   1  2  3
C            DMI =  2  4  5
C                   3  5  6
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
      NITER = 0
      DO I = 1, NG
         AN(I) = 0.0D0
         ZN(I) = 0.0D0
      ENDDO
      DO I = 1, NS
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         DETM = 1.0D0/( G(II1) * G(II4) * G(II6)
     &               -  G(II1) * G(II5) * G(II5)
     &               -  G(II2) * G(II2) * G(II6) )
         DMI(1,I) = ( G(II4) * G(II6) - G(II5) * G(II5) ) * DETM
         DMI(2,I) = - G(II2) * G(II6) * DETM
         DMI(3,I) =   G(II2) * G(II5) * DETM
         DMI(4,I) =   G(II1) * G(II6) * DETM
         DMI(5,I) = - G(II5) * G(II1) * DETM
         DMI(6,I) = ( G(II1) * G(II4) - G(II2) * G(II2) ) * DETM
      ENDDO
C-----------------------------------------------------------------------
      BETAN = 0.0D0
      AAA = 0.0D0
      DO I = 1, NS
         FAC = RN(I) * (DMI(1,I)*RN(I) + DMI(2,I)*RN(I+NS) +
     &                  DMI(3,I)*RN(I+NS2))
     &       + RN(I+NS) * (DMI(2,I)*RN(I) + DMI(4,I)*RN(I+NS) + 
     &                     DMI(5,I)*RN(I+NS2))
     &       + RN(I+NS2) * (DMI(3,I)*RN(I) + DMI(5,I)*RN(I+NS) + 
     &                      DMI(6,I)*RN(I+NS2)) 
         AAA = AAA + FAC
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NS
         ZN(I) = DMI(1,I)*RN(I) + DMI(2,I)*RN(I+NS) +
     &           DMI(3,I)*RN(I+NS2) + BETAN*ZN(I)
         ZN(I+NS) = DMI(2,I)*RN(I) + DMI(4,I)*RN(I+NS) +
     &              DMI(5,I)*RN(I+NS2) + BETAN*ZN(I+NS)
         ZN(I+NS2) = DMI(3,I)*RN(I) + DMI(5,I)*RN(I+NS) +
     &               DMI(6,I)*RN(I+NS2) + BETAN*ZN(I+NS2)
      ENDDO
      CALL EGSAXS(NG, G, ZN, TEMP)
      BBB = VDDOT (NG, ZN, 1, TEMP, 1)
      DO I = 1, NG
         AN(I) = AN(I) + AAA/BBB*ZN(I)
         RN(I) = RN(I) - AAA/BBB*TEMP(I)
      ENDDO
      CCC = 0.0D0
      DO I = 1, NS
         FAC = RN(I) * (DMI(1,I)*RN(I) + DMI(2,I)*RN(I+NS) +
     &                  DMI(3,I)*RN(I+NS2))
     &       + RN(I+NS) * (DMI(2,I)*RN(I) + DMI(4,I)*RN(I+NS) + 
     &                     DMI(5,I)*RN(I+NS2))
     &       + RN(I+NS2) * (DMI(3,I)*RN(I) + DMI(5,I)*RN(I+NS) + 
     &                      DMI(6,I)*RN(I+NS2)) 
         CCC = CCC + FAC
      ENDDO
      BETAN = CCC/AAA
      AAA   = CCC
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSSI1 ( NS, G, DMI, TEMP, XTR, YTR, WT, WW,
     &                   Y, D, ITERMX, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*), DMI(NS), TEMP(NS), YTR(NS), Y(NS), D(NS,NS),
     &          XTR(NS), WT(NS), AUX(NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed.
C-----------------------------------------------------------------------
      CALL EGZERO ( NS * NS, D )
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         DMI(I) = AUX(I) / G(III)
         D(I,I) = DMI(I)
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            D(I,I) = DMI(I) * ( 1.0D0 + YTR(I) )
            DO J = I+1, NS
               IND = NS*(I-1) - (I*(I-1))/2 + J
               D(I,J) = - DMI(I) * DMI(J) * G(IND)
               D(J,I) = D(I,J)
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Project the resulting iterate
C-----------------------------------------------------------------------
      RWW = 1.0D0 / WW
      DO I = 1, NS
         FAC = XTR(I) * WT(I) * RWW
         DO J = 1, NS
            D(I,J) = FAC * D(I,J)
         ENDDO
      ENDDO
C.....
      SUMY = 0.0D0
      DO I = 1, NS
         SUMY = SUMY + Y(I)
      ENDDO
      RSUMY = 1.0D0 / SUMY
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         AAA = Y(J)
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(I,J) * AAA
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = TEMP(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA
         ENDDO
      ENDDO
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(J,I) 
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = Y(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA * TEMP(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN                                                           
      END                                                              
C=======================================================================
C=======================================================================
      SUBROUTINE EGSSI2 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW,
     &                   Y, D, ITERMX, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*), DMI(3,NS), TEMP(NG), YTR(NS), Y(NS), D(NS,NS),
     &          XTR(NS), WT(NS), AUX(NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      CALL EGZERO ( NS * NS, D )
      DO I = 1, NS
         FAC = 1.0D0 / AUX(I)
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         DETM = 1.0D0 / ( FAC*G(II1)*G(II4) - G(II2)*G(II2) )
         DMI(1,I) =   G(II4) * DETM
         DMI(2,I) = - G(II2) * DETM
         D(I,I) = DMI(1,I)
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            II1 = NS2*(I-1) - (I*(I-1))/2 + I
            II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
            II3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
            D(I,I) = 2.0D0 * DMI(1,I) - DMI(1,I)*DMI(1,I)*G(II1)
     &                                - DMI(1,I)*DMI(2,I)*G(II2)*2.0D0
     &                                - DMI(2,I)*DMI(2,I)*G(II3)
            DO J = I+1, NS
               IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
               IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
               JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
               IJ3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C..............
               D(I,J) = - ( DMI(1,I) * DMI(1,J) * G(IJ1)
     &                    + DMI(1,I) * DMI(2,J) * G(IJ2)
     &                    + DMI(2,I) * DMI(1,J) * G(JI2)
     &                    + DMI(2,I) * DMI(2,J) * G(IJ3) )
               D(J,I) = D(I,J)
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Project the resulting iterate
C-----------------------------------------------------------------------
      RWW = 1.0D0 / WW
      DO I = 1, NS
         FAC = XTR(I) * WT(I) * RWW
         DO J = 1, NS
            D(I,J) = FAC * D(I,J)
         ENDDO
      ENDDO
C.....
      SUMY = 0.0D0
      DO I = 1, NS
         SUMY = SUMY + Y(I)
      ENDDO
      RSUMY = 1.0D0 / SUMY
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         AAA = Y(J)
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(I,J) * AAA
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = TEMP(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA
         ENDDO
      ENDDO
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(J,I) 
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = Y(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA * TEMP(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN                                                           
      END                                                              
C=======================================================================
C=======================================================================
      SUBROUTINE EGSSI3 ( NS, NG, G, DMI, TEMP, XTR, YTR, WT, WW,
     &                   Y, D, ITERMX, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*), DMI(6,NS), TEMP(NG), YTR(NS), Y(NS), D(NS,NS),
     &          XTR(NS), WT(NS), AUX(NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
      CALL EGZERO ( NS * NS, D )
      DO I = 1, NS
         FAC = 1.0D0 / AUX(I)
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         DETM = 1.0D0/(FAC*G(II1)*G(II4)*G(II6)
     &               - FAC*G(II1)*G(II5)*G(II5)
     &               - G(II2)*G(II2)*G(II6))
         DMI(1,I) = ( G(II4) * G(II6) - G(II5) * G(II5) ) * DETM
         DMI(2,I) = - G(II2) * G(II6) * DETM
         DMI(3,I) =   G(II2) * G(II5) * DETM
         D(I,I) = DMI(1,I)
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            III = NS3*(I-1) - (I*(I-1))/2 + I
            D(I,I) = DMI(1,I) * ( 1.0D0 + G(III) * DMI(1,I)
     &                                  * YTR(I) / AUX(I) )
            DO J = I+1, NS
               IJ1 = NS3*(I-1) - (I*(I-1))/2 + J
               IJ2 = NS3*(I-1) - (I*(I-1))/2 + J + NS
               IJ3 = NS3*(I-1) - (I*(I-1))/2 + J + NS2
               IJ4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
               IJ5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS2
               IJ6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + J + NS2
               JI2 = NS3*(J-1) - (J*(J-1))/2 + I + NS
               JI3 = NS3*(J-1) - (J*(J-1))/2 + I + NS2
               JI5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + I + NS2
C..............
               D(I,J) = - ( DMI(1,I) * DMI(1,J) * G(IJ1)
     &                    + DMI(1,I) * DMI(2,J) * G(IJ2)
     &                    + DMI(1,I) * DMI(3,J) * G(IJ3)
     &                    + DMI(2,I) * DMI(1,J) * G(JI2)
     &                    + DMI(2,I) * DMI(2,J) * G(IJ4)
     &                    + DMI(2,I) * DMI(3,J) * G(IJ5)
     &                    + DMI(3,I) * DMI(1,J) * G(JI3)
     &                    + DMI(3,I) * DMI(2,J) * G(JI5)
     &                    + DMI(3,I) * DMI(3,J) * G(IJ6) )
               D(J,I) = D(I,J)
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Project the resulting iterate
C-----------------------------------------------------------------------
      RWW = 1.0D0 / WW
      DO I = 1, NS
         FAC = XTR(I) * WT(I) * RWW
         DO J = 1, NS
            D(I,J) = FAC * D(I,J)
         ENDDO
      ENDDO
C.....
      SUMY = 0.0D0
      DO I = 1, NS
         SUMY = SUMY + Y(I)
      ENDDO
      RSUMY = 1.0D0 / SUMY
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         AAA = Y(J)
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(I,J) * AAA
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = TEMP(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA
         ENDDO
      ENDDO
      CALL EGZERO ( NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            TEMP(I) = TEMP(I) + D(J,I) 
         ENDDO
      ENDDO
      DO I = 1, NS
         AAA = Y(I) * RSUMY
         DO J = 1, NS
            D(I,J) = D(I,J) - AAA * TEMP(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN                                                           
      END                                                              
C=======================================================================
C=======================================================================
      SUBROUTINE EGSDDEC ( NS, XTR, YTR, G, BIN, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Output
C     ------
C        G(*)      Choleski decomposition of the positive
C                  definite version of matrix L_[00]
C                  dimension G(*) at least NS*(NS+1)/2 
C
C-----------------------------------------------------------------------
      DIMENSION G(*), YTR(NS)
C-----------------------------------------------------------------------
      CALL EGSEML00 ( NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         G(II1) = G(II1) + YTR(I) * YTR(I)
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            G(IJ1) = G(IJ1) + YTR(I) * YTR(J)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGSDEC ( NG, G, TEMP, IER )
      IF ( IER .NE. 0 ) THEN
         WRITE(*,'(1X,''stopping in EGSDDEC with IER = '',i3)') IER
         STOP
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSDEC (N, A, W, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     A is stored in symmetric form
C     A(i,j) --> A(ij) with ij = N*(i-1) - i(i-1)/2 + j for j >= i
C-----------------------------------------------------------------------
      DIMENSION A(*), W(N)
      IER = 0
      DO  K = 1,N
         KM1 = K - 1
         KP1 = K + 1
         DO  J = 1, KM1
            jj = n*(j-1) - (j*(j-1))/2 + j
            kj = n*(j-1) - (j*(j-1))/2 + k
            W(J)= A(JJ)*A(KJ)
         ENDDO
         kk = n*(k-1) - (k*(k-1))/2 + k
         DO J = 1, KM1
            kj = n*(j-1) - (j*(j-1))/2 + k
            A(KK) = A(KK) - W(J)*A(KJ)
         ENDDO
         IF (A(KK) .EQ. 0.0D0) THEN
            WRITE(6, '(''SINGULAR MATRIX IN EGSDEC'')' )
            IER = K
            RETURN
         ENDIF
         FAC = 1.0D0/A(KK)
         DO J = 1, KM1
            DO I = KP1, N
               ik = n*(k-1) - (k*(k-1))/2 + i
               ij = n*(j-1) - (j*(j-1))/2 + i
               A(IK) = A(IK) - A(IJ)*W(J)
            ENDDO
         ENDDO
         DO I = KP1, N
            ik = n*(k-1) - (k*(k-1))/2 + i
            A(IK) = A(IK)/A(KK)
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSSOL (N, A, B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*), B(N)
      NM1 = N - 1
      DO J = 1, NM1
         JP1 = J + 1
         FAC = -B(J)
         DO  K = JP1, N
            kj = n*(j-1) - (j*(j-1))/2 + k
            B(K) = B(K) + A(KJ)*FAC
         ENDDO
      ENDDO
      DO J = 1, N
         jj = n*(j-1) - (j*(j-1))/2 + j
         B(J) = B(J)/A(JJ)
      ENDDO
      DO JB = 1, NM1
         J = N + 1 - JB
         JM1 = J - 1
         FAC = -B(J)
         DO K = 1, JM1
            jk = n*(k-1) - (k*(k-1))/2 + j
            B(K) = B(K) + A(JK)*FAC
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMH ( NS, T, X, WT, RU, PATMOS, 
     &                    BETA, G, AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NS), BIN(*), WT(NS), AIJ(NS,NS), ETA(NS)
      DIMENSION G(*), BETA(NS)
      PARAMETER (CCC = 5.0D0 / 3.0D0)
      DOUBLE PRECISION, ALLOCATABLE :: IWT(:)
      ALLOCATE(IWT(NS))
C-----------------------------------------------------------------------
C         EVALUATE THE RHS BETA
C-----------------------------------------------------------------------
      DO K = 1, NS
         BETA(K) = X(K)
         IWT(K) = 1.0D0 / WT(K)
      ENDDO
C-----------------------------------------------------------------------
C         EVALUATE THE MATRIX H
C
C     Note:    FAC * BIN = 2 W_{ij} / \eta_{ij} / A_{ij}
C-----------------------------------------------------------------------
      FAC = 6.0D0 * RU * T / ( 5.0D0 * PATMOS )
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         G(III) = X(I) * X(I) / ETA(I)
      ENDDO
      DO I = 1, NS
         IP1 = I + 1
         III = NS*(I-1) - (I*(I-1))/2 + I
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
            JJJ = NS*(J-1) - (J*(J-1))/2 + J
            AAA = X(I) * X(J) * FAC * BIN(IND) / ( WT(I) + WT(J) )
            G(IND) = AAA * ( AIJ(I,J) - CCC )
            G(III) = G(III) + AAA * 
     &                      ( AIJ(I,J)*WT(J)*IWT(I) + CCC )
            G(JJJ) = G(JJJ) + AAA * 
     &                      ( AIJ(I,J)*WT(I)*IWT(J) + CCC )
         ENDDO
      ENDDO
      DEALLOCATE(IWT)
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMK ( NS, TEMPER, X, WT, RU, PATMOS,
     &                    BIN, AIJ, CINT, ETA, CXI,
     &                    G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      CCC = 0.0D0
      DO I = 1, NS
         CCC = CCC + X(I) * CINT(I)
      ENDDO
      CV  = CCC + 1.5D0
      RCV = 1.0D0 / CV
      DO I = 1, NS
         BETA(I)      = X(I) * CCC * RCV
         BETA(I+NS)   = - X(I) * CINT(I) * RCV
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS2*(I-1)     - (I*(I-1))/2           + I
         II2 = NS2*(I-1)     - (I*(I-1))/2           + I + NS
         II3 = NS2*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         FAC = 4.0D0 * CXI(I) / ETA(I) * X(I) * X(I)
         G(II1) = FAC 
         G(II2) = - FAC 
         G(II3) = FAC 
      ENDDO
C-----------------------------------------------------------------------
      DEN = 2.4D0 * RU * TEMPER / PATMOS
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IBN = NS*(I-1) - (I*(I-1))/2 + J
C...........
            JJ1 = NS2*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS2*(J-1) - (J*(J-1))/2 + J + NS
            JJ3 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
            IJ3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            BB1  = X(I) * X(J) * BIN(IBN) * DEN * AIJ(I,J)
            BBB  = BB1 / ( WT(I) + WT(J) )
C ----- 10,10
            AAA = X(I) * X(J) * BIN(IBN) * 1.25D0 * DEN 
     &                 / ( WT(I) + WT(J) )
            BBJ = BBB * WT(J) / WT(I)
            BBI = BBB * WT(I) / WT(J)
            G(IJ1) = - AAA + BBB * ( CXI(I) + CXI(J) )
            G(II1) = G(II1) + AAA + BBJ * ( CXI(I) + CXI(J) )
            G(JJ1) = G(JJ1) + AAA + BBI * ( CXI(I) + CXI(J) )
C ----- 10,01 and 01,10
            G(IJ2) = - BB1 * CXI(J) / WT(I)
            G(JI2) = - BB1 * CXI(I) / WT(J)
            G(II2) = G(II2) - BB1 * CXI(I) / WT(I)
            G(JJ2) = G(JJ2) - BB1 * CXI(J) / WT(J)
C ----- 01,01
            BB2 = BB1 * ( WT(I) + WT(J) ) / WT(I) / WT(J)
            G(IJ3) = 0.0D0
            G(II3) = G(II3) + BB2 * CXI(I)
            G(JJ3) = G(JJ3) + BB2 * CXI(J)
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II3 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(I-NS) .EQ. 0.0D0) G(II3) = 1.0D0
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMK01 ( NS, TEMPER, X, WT, RU, PATMOS,
     &                      BIN, AIJ, CINT, ETA, CXI,
     &                      G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      CCC = 0.0D0
      DO I = 1, NS
         CCC = CCC + X(I) * CINT(I)
      ENDDO
      CV  = CCC + 1.5D0
      RCV = 1.0D0 / CV
      DO I = 1, NS
         BETA(I)   = - X(I) * CINT(I) * RCV
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K_[01]
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         FAC = 4.0D0 * CXI(I) / ETA(I) * X(I) * X(I)
         G(II1) = FAC 
      ENDDO
C-----------------------------------------------------------------------
      DEN = 2.4D0 * RU * TEMPER / PATMOS
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ1 = NS*(J-1) - (J*(J-1))/2 + J
            IJ1 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            BB1  = X(I) * X(J) * BIN(IJ1) * DEN * AIJ(I,J)
C ----- 01,01
            BB2 = BB1 * ( WT(I) + WT(J) ) / WT(I) / WT(J)
            G(IJ1) = 0.0D0
            G(II1) = G(II1) + BB2 * CXI(I)
            G(JJ1) = G(JJ1) + BB2 * CXI(J)
         ENDDO
      ENDDO
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(I) .EQ. 0.0D0) G(II1) = 1.0D0
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEML ( NS, TEMPER, X, WT, RU, PATMOS,
     &                    BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &                    LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), BIJ(NS,*), CIJ(NS,*), 
     &          CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*), LIN(*)
      PARAMETER(DDD = 2.0D1 / 3.0D0)
      PARAMETER(EEE = 4.0D0 / 1.5D1)
      PARAMETER(FFF = 2.0D1 / 3.0D0)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         BETA(I)      = 0.0D0
         BETA(I+NS)   = 2.5D0 * X(I)
         BETA(I+2*NS) = CINT(I) * X(I)
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS3*(I-1)     - (I*(I-1))/2           + I
         II2 = NS3*(I-1)     - (I*(I-1))/2           + I + NS
         II3 = NS3*(I-1)     - (I*(I-1))/2           + I + NS2
         II4 = NS3*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         II5 = NS3*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         G(II1) = 0.0D0
         G(II2) = 0.0D0
         G(II3) = 0.0D0
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         FAC = (5.0D0*PATMOS/(3.0D0*TEMPER)) * WT(I) / RU / ETA(I)
         IBN = NS*(I-1) - (I*(I-1))/2 + I
         G(II4) = FAC * X(I) * X(I) 
     &                   * ( 1.0D0 + 1.0D1 * CXI(I) / 3.0D0 )
         G(II5) = - FAC * X(I) * X(I) 
     &                   * 2.0D0 * CXI(I) 
         G(II6) = FAC * X(I) * X(I) 
     &                   * 1.2D0 * CXI(I)
     &                   + X(I) * X(I) * CINT(I) * BIN(IBN)
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
C........
         IP1 = I + 1
         DO J = IP1, NS
            IBN = NS*(I-1) - (I*(I-1))/2 + J
C...........
            JJ1 = NS3*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS3*(J-1) - (J*(J-1))/2 + J + NS
            JJ3 = NS3*(J-1) - (J*(J-1))/2 + J + NS2
            JJ4 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
            JJ5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS2
            JJ6 = NS3*(J+NS2-1) - ((J+NS2)*(J+NS2-1))/2 + J + NS2
C...........
            IJ1 = NS3*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS3*(I-1) - (I*(I-1))/2 + J + NS
            IJ3 = NS3*(I-1) - (I*(I-1))/2 + J + NS2
            IJ4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
            IJ5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS2
            IJ6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + J + NS2
C...........
            JI2 = NS3*(J-1) - (J*(J-1))/2 + I + NS
            JI3 = NS3*(J-1) - (J*(J-1))/2 + I + NS2
            JI5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + I + NS2
C...........
            WWI = WT(I) / ( WT(I) + WT(J) )
            WWJ = WT(J) / ( WT(I) + WT(J) )
C ----- 00,00
            AAA = X(I) * X(J) * BIN(IBN)
            G(IJ1) = - AAA
            G(II1) = G(II1) + AAA
            G(JJ1) = G(JJ1) + AAA
C ----- 00,10 and 10,00
            AAI = 0.5D0 * WWI * (6.0D0*CIJ(I,J) - 5.0D0)
            AAJ = 0.5D0 * WWJ * (6.0D0*CIJ(I,J) - 5.0D0)
            G(IJ2) = AAI * AAA
            G(JI2) = AAJ * AAA
            G(II2) = G(II2) - AAJ * AAA
            G(JJ2) = G(JJ2) - AAI * AAA
C ----- 00,01 and 01,00
            G(IJ3) = 0.0D0
            G(JI3) = 0.0D0
C ----- 10,10
            BB1 = DDD * AIJ(I,J) * (CXI(I)+CXI(J))
            BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J)
     &                - BB1 )
            BBJ = ( 7.5D0 * WWI * WWI + 
     &               (6.25D0 - 3.0D0*BIJ(I,J)) * WWJ * WWJ 
     &                + 4.0D0*AIJ(I,J) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
            BBI = ( 7.5D0 * WWJ * WWJ + 
     &               (6.25D0 - 3.0D0*BIJ(I,J)) * WWI * WWI 
     &                + 4.0D0*AIJ(I,J) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
            G(IJ4) = - BBB * AAA
            G(II4) = G(II4) + BBJ * AAA
            G(JJ4) = G(JJ4) + BBI * AAA
C ----- 10,01 and 01,10
            CCI = WWI * 4.0D0 * AIJ(I,J) * CXI(I)
            CCJ = WWJ * 4.0D0 * AIJ(I,J) * CXI(J)
            G(IJ5) = - CCJ * AAA
            G(JI5) = - CCI * AAA
            G(II5) = G(II5) - CCI * AAA
            G(JJ5) = G(JJ5) - CCJ * AAA
C ----- 01,01
            DDJ = AAA * CINT(I) 
     &          + AAA * 2.4D0 * AIJ(I,J) *(WT(I)/WT(J)) * CXI(I)
            DDI = AAA * CINT(J) 
     &          + AAA * 2.4D0 * AIJ(I,J) *(WT(J)/WT(I)) * CXI(J)
            G(IJ6) = 0.0D0
            G(II6) = G(II6) + DDJ
            G(JJ6) = G(JJ6) + DDI
         ENDDO
      ENDDO
      DO I = NS2+1, NS3
         II6 = NS3*(I-1) - (I*(I-1))/2 + I
         IF ( LIN(I-NS2) .EQ. 0) G(II6) = 1.0D0
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMA ( NS, TEMPER, X, WT, RU, PATMOS,
     &                    BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &                    LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), BIJ(NS,*), CIJ(NS,*), 
     &          CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*), LIN(*)
      PARAMETER(DDD = 2.0D1 / 3.0D0)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         BETA(I)    = 2.5D0 * X(I)
         BETA(I+NS) = CINT(I) * X(I)
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         FAC = (5.0D0*PATMOS/(3.0D0*TEMPER)) * WT(I) / RU / ETA(I)
         IBN = NS*(I-1) - (I*(I-1))/2 + I
         G(II4) = FAC * X(I) * X(I) 
     &                   * ( 1.0D0 + 1.0D1 * CXI(I) / 3.0D0 )
         G(II5) = - FAC * X(I) * X(I) 
     &                   * 2.0D0 * CXI(I) 
         G(II6) = FAC * X(I) * X(I) 
     &                   * 1.2D0 * CXI(I)
     &                   + X(I) * X(I) * CINT(I) * BIN(IBN)
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IBN = NS*(I-1) - (I*(I-1))/2 + J
C...........
            JJ4 = NS2*(J-1)    - (J*(J-1))/2         + J 
            JJ5 = NS2*(J-1)    - (J*(J-1))/2         + J + NS
            JJ6 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ4 = NS2*(I-1)    - (I*(I-1))/2         + J 
            IJ5 = NS2*(I-1)    - (I*(I-1))/2         + J + NS
            IJ6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI5 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            WWI = WT(I) / ( WT(I) + WT(J) )
            WWJ = WT(J) / ( WT(I) + WT(J) )
            AAA = X(I) * X(J) * BIN(IBN)
C ----- 10,10
            BB1 = DDD * AIJ(I,J) * (CXI(I)+CXI(J))
            BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J)
     &                - BB1 )
            BBJ = ( 7.5D0 * WWI * WWI + 
     &               (6.25D0 - 3.0D0*BIJ(I,J)) * WWJ * WWJ 
     &                + 4.0D0*AIJ(I,J) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
            BBI = ( 7.5D0 * WWJ * WWJ + 
     &               (6.25D0 - 3.0D0*BIJ(I,J)) * WWI * WWI 
     &                + 4.0D0*AIJ(I,J) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
            G(IJ4) = - BBB * AAA
            G(II4) = G(II4) + BBJ * AAA
            G(JJ4) = G(JJ4) + BBI * AAA
C ----- 10,01 and 01,10
            CCI = WWI * 4.0D0 * AIJ(I,J) * CXI(I)
            CCJ = WWJ * 4.0D0 * AIJ(I,J) * CXI(J)
            G(IJ5) = - CCJ * AAA
            G(JI5) = - CCI * AAA
            G(II5) = G(II5) - CCI * AAA
            G(JJ5) = G(JJ5) - CCJ * AAA
C ----- 01,01
            DDJ = AAA * CINT(I) 
     &          + AAA * 2.4D0 * AIJ(I,J) *(WT(I)/WT(J)) * CXI(I)
            DDI = AAA * CINT(J) 
     &          + AAA * 2.4D0 * AIJ(I,J) *(WT(J)/WT(I)) * CXI(J)
            G(IJ6) = 0.0D0
            G(II6) = G(II6) + DDJ
            G(JJ6) = G(JJ6) + DDI
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II6 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( LIN(I-NS) .EQ. 0) G(II6) = 1.0D0
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEML00 ( NS, X, BIN, G )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BIN(*), X(*)
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L_[00]
C
C     Note: the system matrix is evaluated at atmospheric pressure.
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         G(III) = 0.0D0
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         IP1 = I + 1
         III = NS*(I-1) - (I*(I-1))/2 + I
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
            JJJ = NS*(J-1) - (J*(J-1))/2 + J
            AAA = X(I) * X(J) * BIN(IND)
            G(IND) = - AAA
            G(III) = G(III) + AAA
            G(JJJ) = G(JJJ) + AAA
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMLE ( NS, TEMPER, X, WT, RU, PATMOS,
     &                     BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &                     LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), BIJ(NS,*), CIJ(NS,*), 
     &          CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*), LIN(*)
      PARAMETER(DDD = 4.0D0 / 3.0D0)
      PARAMETER(EEE = 4.0D0 / 1.5D1)
      PARAMETER(FFF = 2.0D1 / 3.0D0)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         BETA(I)      = 0.0D0
         BETA(I+NS)   = ( 2.5D0 + CINT(I) ) * X(I)
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L_[e]
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS2*(I-1)     - (I*(I-1))/2           + I
         II2 = NS2*(I-1)     - (I*(I-1))/2           + I + NS
         II4 = NS2*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         G(II1) = 0.0D0
         G(II2) = 0.0D0
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         FAC = (5.0D0*PATMOS/(3.0D0*TEMPER)) * WT(I) / RU / ETA(I)
         IBN = NS*(I-1) - (I*(I-1))/2 + I
         G(II4) = FAC * X(I) * X(I) 
     &                   * ( 1.0D0 + 8.0D0 * CXI(I) / 1.5D1 )
     &                + X(I) * X(I) * CINT(I) * BIN(IBN)
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IBN = NS*(I-1) - (I*(I-1))/2 + J
C...........
            JJ1 = NS2*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS2*(J-1) - (J*(J-1))/2 + J + NS
            JJ4 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
            IJ4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            WWI  = WT(I) / ( WT(I) + WT(J) )
            WWJ  = WT(J) / ( WT(I) + WT(J) )
            WIOJ = 3.0D0 * WT(I) / WT(J) - 2.0D0
            WJOI = 3.0D0 * WT(J) / WT(I) - 2.0D0
C ----- 00,00
            AAA = X(I) * X(J) * BIN(IBN)
            G(IJ1) = - AAA
            G(II1) = G(II1) + AAA
            G(JJ1) = G(JJ1) + AAA
C ----- 00,e and e,00
            AAI = 0.5D0 * WWI * (6.0D0*CIJ(I,J) - 5.0D0)
            AAJ = 0.5D0 * WWJ * (6.0D0*CIJ(I,J) - 5.0D0)
            G(IJ2) = AAI * AAA
            G(JI2) = AAJ * AAA
            G(II2) = G(II2) - AAJ * AAA
            G(JJ2) = G(JJ2) - AAI * AAA
C ----- e,e
            TTI = DDD * AIJ(I,J) * WIOJ * CXI(I)
            TTJ = DDD * AIJ(I,J) * WJOI * CXI(J)
            BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J)
     &                + TTI + TTJ )
            BBJ = 7.5D0 * WWI * WWI  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWJ * WWJ 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
     &            + EEE * AIJ(I,J) * CXI(I)
     &                          * WWI * WWJ * WIOJ * WIOJ
     &            + FFF * AIJ(I,J) * CXI(J) * WWI * WWJ
            BBI = 7.5D0 * WWJ * WWJ  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWI * WWI 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
     &            + EEE * AIJ(I,J) * CXI(J)
     &                          * WWI * WWJ * WJOI * WJOI
     &            + FFF * AIJ(I,J) * CXI(I) * WWI * WWJ
            DDJ = AAA * CINT(I) 
            DDI = AAA * CINT(J) 
            G(IJ4) = - BBB * AAA
            G(II4) = G(II4) + BBJ * AAA + DDJ
            G(JJ4) = G(JJ4) + BBI * AAA + DDI
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGSEMAE ( NS, TEMPER, X, WT, RU, PATMOS,
     &                     BIN, AIJ, BIJ, CIJ, CINT, ETA, CXI,
     &                     LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(*), BETA(*)
      DIMENSION BIN(*), AIJ(NS,*), BIJ(NS,*), CIJ(NS,*), 
     &          CINT(NS), ETA(NS), X(NS), 
     &          WT(NS), CXI(*), LIN(*)
      PARAMETER(DDD = 4.0D0 / 3.0D0)
      PARAMETER(EEE = 4.0D0 / 1.5D1)
      PARAMETER(FFF = 2.0D1 / 3.0D0)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         BETA(I) = ( 2.5D0 + CINT(I) ) * X(I)
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda_[e]
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II4 = NS*(I-1)  - (I*(I-1))/2   + I 
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         FAC = (5.0D0*PATMOS/(3.0D0*TEMPER)) * WT(I) / RU / ETA(I)
         IBN = NS*(I-1) - (I*(I-1))/2 + I
         G(II4) = FAC * X(I) * X(I) 
     &                   * ( 1.0D0 + 8.0D0 * CXI(I) / 1.5D1 )
     &                + X(I) * X(I) * CINT(I) * BIN(IBN)
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         II4 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ4 = NS*(J-1) - (J*(J-1))/2 + J
            IJ4 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            WWI  = WT(I) / ( WT(I) + WT(J) )
            WWJ  = WT(J) / ( WT(I) + WT(J) )
            WIOJ = 3.0D0 * WT(I) / WT(J) - 2.0D0
            WJOI = 3.0D0 * WT(J) / WT(I) - 2.0D0
            TTI = DDD * AIJ(I,J) * WIOJ * CXI(I)
            TTJ = DDD * AIJ(I,J) * WJOI * CXI(J)
            AAA = X(I) * X(J) * BIN(IJ4)
            BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J)
     &                + TTI + TTJ )
            BBJ = 7.5D0 * WWI * WWI  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWJ * WWJ 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
     &            + EEE * AIJ(I,J) * CXI(I)
     &                          * WWI * WWJ * WIOJ * WIOJ
     &            + FFF * AIJ(I,J) * CXI(J) * WWI * WWJ
            BBI = 7.5D0 * WWJ * WWJ  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWI * WWI 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
     &            + EEE * AIJ(I,J) * CXI(J)
     &                          * WWI * WWJ * WJOI * WJOI
     &            + FFF * AIJ(I,J) * CXI(I) * WWI * WWJ
            DDJ = AAA * CINT(I) 
            DDI = AAA * CINT(J) 
C...........
            G(IJ4) = - BBB * AAA
            G(II4) = G(II4) + BBJ * AAA + DDJ
            G(JJ4) = G(JJ4) + BBI * AAA + DDI
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C-----------------------------------------------------------------------

