      SUBROUTINE EGINI (NP, LOUT, IFLAG, ITLS, 
     &                  WEG, LWEG, IWEG, LIWEG )
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
      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
        READ (LLEG) NSLK, NO
      CLOSE(UNIT=LLEG)
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
      OPEN (UNIT=LLEG,STATUS='OLD',FORM='UNFORMATTED',FILE='Linkeg')
        READ (LLEG) NSLK, NO, (WEG(IEGWT+K-1), K=1, NS), 
     &              (WEG(IEGEPS+K-1), K=1, NS), 
     &              (WEG(IEGSIG+K-1), K=1, NS), 
     &              (WEG(IEGDIP+K-1), K=1, NS), 
     &              (WEG(IEGPOL+K-1), K=1, NS), 
     &              (WEG(IEGZRT+K-1), K=1, NS), 
     &              (IWEG(IEGLIN+K-1), K=1, NS), 
     &              (WEG(IEGCFE+N-1), N=1, NONS), 
     &              (WEG(IEGCFL+N-1), N=1, NONS), 
     &              (WEG(IEGCFD+N-1), N=1, NONSNS)
      CLOSE(UNIT=LLEG)
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
C=======================================================================
C=======================================================================
      SUBROUTINE EGABC ( NS, NFIT, FITA, FITB, FITC,
     &                   FITA0, FITB0, FITC0, EPSIJ )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION FITA(NFIT,NS,NS), FITB(NFIT,NS,NS), FITC(NFIT,NS,NS),
     &          FITA0(NFIT), FITB0(NFIT), FITC0(NFIT),
     &          EPSIJ(NS,NS)
C-----------------------------------------------------------------------
      DO J = 1, NS
         DO I = J, NS
            do m = 1, nfit
               SUMA = 0.0D0
               SUMB = 0.0D0
               SUMC = 0.0D0
               mm   = m - 1
               do k = mm, nfit-1
                  prod = 1.0d0
                  do l = 1, k-mm
                     prod = prod * (-epsij(i,j)) * dfloat(mm+l)
     &                           / dfloat(l)
                  enddo
                  SUMA = SUMA + FITA0(k+1) * PROD
                  SUMB = SUMB + FITB0(k+1) * PROD
                  SUMC = SUMC + FITC0(k+1) * PROD
               enddo
               FITA(m,I,J) = SUMA
               FITB(m,I,J) = SUMB
               FITC(m,I,J) = SUMC
            enddo
            IF ( I .GT. J ) THEN
               do m = 1, nfit
                  FITA(m,J,I) = FITA(m,I,J)
                  FITB(m,J,I) = FITB(m,I,J)
                  FITC(m,J,I) = FITC(m,I,J)
               enddo
            ENDIF
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE LEVEPS (NS, EPS, SIG, DIP, POL, EPSIJ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C-----------------------------------------------------------------------
      DIMENSION EPS(NS), SIG(NS), DIP(NS), POL(NS), EPSIJ(NS,NS)
      DATA PI/3.1415926535D0/, FAC/1.0D-12/, 
     &     DIPMIN/1.0D-20/, BOLTZ/1.38056D-16/
C-----------------------------------------------------------------------
C     This subroutine computes the Lennard-Jones potentials well
C     depths. It is called only once at the beginning.
C-----------------------------------------------------------------------
      DO 1000 J = 1, NS
         DO 1000 K = 1, J
           IF((DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN)) THEN
C-----------------------------------------------------------------------
C                K IS POLAR, J IS NONPOLAR
C-----------------------------------------------------------------------
              XI = 1.0D0 + 0.25D0*(POL(J)/SIG(J)**3) *
     1                     (FAC/BOLTZ) *
     2                     (DIP(K)**2/(EPS(K)*SIG(K)**3)) *
     3                      DSQRT(EPS(K)/EPS(J))
c             SGM(K,J) = 0.5D0 * (SIG(J)+SIG(K)) * XI**(-1.0D0/6.0D0)
c             SGM(J,K) = SGM(K,J)
              EPSIJ(K,J) = DSQRT(EPS(J)*EPS(K)) * XI**2
              EPSIJ(K,J) = DLOG ( EPSIJ(K,J) )
              EPSIJ(J,K) = EPSIJ(K,J)
          ELSE IF((DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN)) THEN
C-----------------------------------------------------------------------
C             J IS POLAR, K IS NONPOLAR
C-----------------------------------------------------------------------
              XI = 1.0D0 + 0.25D0*(POL(K)/SIG(K)**3) *
     1                     (FAC/BOLTZ) *
     2                     (DIP(J)**2/(EPS(J)*SIG(J)**3)) *
     3                      DSQRT(EPS(J)/EPS(K))
c             SGM(K,J) = 0.5D0 * (SIG(J)+SIG(K)) * XI**(-1.0D0/6.0D0)
c             SGM(J,K) = SGM(K,J)
              EPSIJ(K,J) = DSQRT(EPS(J)*EPS(K)) * XI**2
              EPSIJ(K,J) = DLOG ( EPSIJ(K,J) )
              EPSIJ(J,K) = EPSIJ(K,J)
          ELSE
C-----------------------------------------------------------------------
C              NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
C-----------------------------------------------------------------------
c             SGM(K,J) = 0.5D0 * (SIG(J) + SIG(K))
c             SGM(J,K) = SGM(K,J)
              EPSIJ(K,J) = DSQRT(EPS(J)*EPS(K))
              EPSIJ(K,J) = DLOG ( EPSIJ(K,J) )
              EPSIJ(J,K) = EPSIJ(K,J)
          ENDIF
1000  CONTINUE
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE  EGZERO ( NN, SX )                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)
      DO I = 1, NN                                                  
         SX(I) = 0.0D0
      ENDDO
      RETURN                                                           
      END                                                              
C=======================================================================
C=======================================================================

