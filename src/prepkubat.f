      SUBROUTINE PREPKUBAT(ND, NF, N, NDIM, XLAMBDA, ENTOTDENS, RNE,
     >                     POPNUM, ELEVEL, NCHARG, XJC, KONTNUP,
     >                     KONTLOW, SIGMAKI, FWEIGHT, WEIGHT, NOM,
     >                     ABXYZ, NFIRST, NLAST, INDNUP, INDLOW,
     >                     KODAT, ENLTE, CRATE, DCRATEDT,
     >                     NATOM, EION, EINST, COCO, KEYCBB,
     >                     KEYCBF, IONGRND, ALTESUM, MAXIND,
     >                     SIGMATHK, SEXPOK, EDGEK, MAXATOM, 
     >                     LASTKON, LASTIND, LASTINDAUTO, LASTINDALL, 
     >                     OPAROSS, TOLD, RADIUS,
     >                     RSTAR, OPALAMBDAMEAN, OPASMEANTC, 
     >                     MAXKONT, ADDCON1, ADDCON2, ADDCON3,
     >                     IGAUNT, ALPHA, SEXPO, MAXION, 
     >                     DTKUBAT, BCOLLIDONE, ZERO_RATES,
     >                     bKUBATDEBUG)
C***********************************************************************
C***  calculates the Q values for the THERMAL BALANCE method 
C***  cf. Kubat et al. (1999), Kubat (2001) 
C***  Using our notation, the formulae are given in the thesis by Andreas Sander. 
C***
C***  called by: STEAL
C***********************************************************************

      IMPLICIT NONE
                    
      INTEGER, INTENT(IN) :: ND, NF, N, NDIM, NATOM, MAXIND, MAXATOM,
     >                       MAXION, LASTKON, MAXKONT, LASTIND,
     >                       LASTINDAUTO, LASTINDALL
      INTEGER, DIMENSION(N), INTENT(IN) :: NOM
      INTEGER, DIMENSION(NDIM), INTENT(IN) :: IONGRND, NCHARG
      INTEGER, DIMENSION(MAXIND), INTENT(IN) :: INDNUP, INDLOW
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: KONTNUP, KONTLOW,
     >                                           KEYCBF
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
      INTEGER, DIMENSION(MAXATOM), INTENT(IN) :: KODAT
      CHARACTER(8), DIMENSION(MAXKONT), INTENT(IN) :: IGAUNT      

      REAL, INTENT(IN) :: RSTAR
      REAL, DIMENSION(NDIM,NDIM), INTENT(IN) :: EINST
      REAL, DIMENSION(4,MAXIND), INTENT(IN) :: COCO
      
      REAL, DIMENSION(4,NDIM), INTENT(IN) :: ALTESUM
      REAL, DIMENSION(NDIM), INTENT(IN) :: ELEVEL, WEIGHT, EION
      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ
      
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPNUM
      REAL, DIMENSION(ND), INTENT(IN) :: ENTOTDENS, RNE, OPAROSS, TOLD, 
     >                                   OPASMEANTC, OPALAMBDAMEAN, 
     >                                   RADIUS
      REAL, DIMENSION(NF), INTENT(IN) :: XLAMBDA, FWEIGHT
      REAL, DIMENSION(ND, NF), INTENT(IN) :: XJC
      REAL, DIMENSION(NF, LASTKON), INTENT(IN) :: SIGMAKI
      REAL, DIMENSION(MAXATOM, MAXION), INTENT(IN) :: 
     >                                        SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(MAXKONT), INTENT(IN) :: 
     >                ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3

      CHARACTER(4), DIMENSION(LASTIND), INTENT(IN) :: KEYCBB      
      
      REAL, DIMENSION(NDIM), INTENT(INOUT) :: ENLTE
      REAL, DIMENSION(NDIM, NDIM), INTENT(INOUT) :: CRATE, DCRATEDT
      
      REAL, DIMENSION(ND), INTENT(OUT) :: DTKUBAT      
      
      REAL, DIMENSION(ND) :: DELTAQ, DELTAQDT
      REAL, DIMENSION(ND, 2) :: QFF, QBF, QC, QFFDT, QBFDT, QCDT   !H is 1, C is 2
      
      LOGICAL, INTENT(IN) :: bKUBATDEBUG
      LOGICAL, DIMENSION(N, ND), INTENT(IN) :: ZERO_RATES
      LOGICAL, DIMENSION(N,N), INTENT(INOUT) :: bCOLLIDONE 
      LOGICAL :: bKSHELL, bPrintQTable
      
      INTEGER :: L, NUP, LOW, I, J, K, IND, KON, KONEDGE, NCHARGE, 
     >           NA, ISTATE, NOMJ
      REAL :: DELTALAM, AFFINT1, AFFINT2, AFFINT2d, PADD1,
     >        POPSUM1, POPSUM1d, POPSUM2, POPSUM2d, PADD2, GIII, 
     >        AFF, W3, W, ABF, ABFINT1, ABFINT2, ABFINT2d, RL, RCML,
     >        ENE, ETRANSIT, TL, EDGELAM, EDGE, ENTOTL, DTB, OPAMEAN,
     >        DGDT, GIIITP, AFFDGDT, AFFINT1d, AFFINT2dd, POPSUM2dd,
     >        WK, SIGMAK, TLmod, POPLOW, POPNUP, POPLJ, XJCLK
     
      INTEGER, EXTERNAL :: ISRCHFGT

      !Constants
      REAL, PARAMETER :: PI4  = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: HPLANCK = 6.625E-27     !PLANCK's CONSTANT (erg s)
      REAL, PARAMETER :: C1 = 1.4388             !C1 = H * C / K    ( CM * K )
      REAL, PARAMETER :: C2 = 3.9724E-16         !C2 = 2 * H * C    ( G * CM**3 / S**2 
      REAL, PARAMETER :: CLIGHT = 2.99792458E10      !Speed of Light in cm/s
      REAL, PARAMETER :: CFF = 1.370E-23         !COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100 )
      REAL, PARAMETER :: COM = 5.465E-11         !Com = a0**2 * sqrt( 8 * Pi * k / m )
      REAL, PARAMETER :: STEBOL = 5.6705E-5      !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5   !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI      
      REAL, PARAMETER :: BOLTZ = 1.38E-16        !BOLTZMANN CONSTANT (ERG/DEG)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      bPrintQTable = bKUBATDEBUG            !print out table if debug option is set
      
      bKSHELL = .FALSE.
      bksloop: DO NA=1, MAXATOM
        DO ISTATE=1, MAXATOM
          IF (SIGMATHK(NA,ISTATE) /= .0) THEN
            bKSHELL = .TRUE.
            EXIT bksloop
          ENDIF
        ENDDO
      ENDDO bksloop
      
      DO L=1, ND
        POPSUM1 = 0.
        POPSUM1d = 0.
        POPSUM2 = 0.
        POPSUM2d = 0.
        POPSUM2dd = 0.
        ENTOTL = ENTOTDENS(L)
        ENE = ENTOTDENS(L) * RNE(L)
        TL = TOLD(L)        
        DTKUBAT(L) = 0.
        QFF(L,1) = 0.
        QFF(L,2) = 0.
        QBF(L,1) = 0.
        QBF(L,2) = 0.
        QC(L,1) = 0.
        QC(L,2) = 0.
        
C****** FREE-FREE transitions ******************************************
        DO J=1, N
          AFFINT1 = 0.      
          AFFINT2 = 0.
          AFFINT1d = 0.
          AFFINT2d = 0.  
          AFFINT2dd = 0.
          NCHARGE = NCHARG(J)
          POPLJ = POPNUM(L, J)
          IF (ZERO_RATES(J, L))  POPLJ = 0.
          DO K=1, NF-1
            W = 1.E8/XLAMBDA(K)
            W3 = W*W*W         
            XJCLK = MAX(0., XJC(L,K))
            IF (NCHARGE > 0) THEN
              CALL GAUNTFF (GIII,NCHARGE,XLAMBDA(K),TL)
              CALL GAUNTFF (GIIITP,NCHARGE,XLAMBDA(K),1.01*TL)
              DGDT = (GIIITP - GIII) / (0.01 * TL)              
              AFF = CFF/W3/SQRT(TL) * NCHARGE**(2.) * GIII   !AlphaFF for Kubat method
              AFFDGDT = CFF/W3/SQRT(TL) * NCHARGE**(2.) * DGDT
            ELSE
              GIIITP = 0.
              GIII = 0.
              AFF = 0.
              AFFDGDT = 0.
            ENDIF
            AFFINT1 = AFFINT1 + FWEIGHT(K) * AFF * XJCLK
            AFFINT1d = AFFINT1d + FWEIGHT(K) * 
     >                                   AFFDGDT * XJCLK
            AFFINT2 = AFFINT2 + FWEIGHT(K) * AFF
     >        * ( XJCLK + C2 * W3) * EXP(-C1*W/TL)
            AFFINT2d = AFFINT2d + FWEIGHT(K) * AFF
     >        * ( XJCLK + C2 * W3) * EXP(-C1*W/TL) * W
            AFFINT2dd = AFFINT2dd + FWEIGHT(K) * AFFDGDT
     >        * ( XJCLK + C2 * W3) * EXP(-C1*W/TL)
          ENDDO
          POPSUM1 = POPSUM1 + ENTOTL*POPLJ * AFFINT1
          POPSUM1d = POPSUM1d + ENTOTL*POPLJ * AFFINT1d
          POPSUM2 = POPSUM2 + ENTOTL*POPLJ * AFFINT2
          POPSUM2d = POPSUM2d + ENTOTL*POPLJ * AFFINT2d
          POPSUM2dd = POPSUM2dd + ENTOTL*POPLJ * AFFINT2dd
        ENDDO
        QFF(L,1) = PI4 * ENE * POPSUM1        
        QFF(L,2) = PI4 * ENE * POPSUM2 
C        QFFDT(L,1) = - 1. / ( 2 * TL) * QFF(L,1) 
        QFFDT(L,1) = - 1. / ( 2 * TL) * QFF(L,1) + PI4 * ENE * POPSUM1d
        QFFDT(L,2) = - 1. / ( 2 * TL) * QFF(L,2) 
     >               + PI4 * ENE * (C1 / TL**2 * POPSUM2d + POPSUM2dd)

        POPSUM1 = 0.
        POPSUM2 = 0.
        POPSUM2d = 0.
C****** BOUND-FREE transitions *****************************************
        CALL LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,
     >               NOM,ABXYZ,NFIRST,NLAST,NATOM)
        CALL BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     >                XLAMBDA,ALPHA,SEXPO,ADDCON1,ADDCON2,ADDCON3,
     >                IGAUNT,KONTNUP,KONTLOW,LASTKON)
     
        DO KON=1, LASTKON
          NUP=KONTNUP(KON)
          LOW=KONTLOW(KON)
          POPLOW = POPNUM(L,LOW)
          POPNUP = POPNUM(L,NUP) 
          IF (ZERO_RATES(LOW,L)) POPLOW = 0.
          IF (ZERO_RATES(NUP,L)) THEN
            POPNUP = 0.
            POPLOW = 0.
          ENDIF
          IF (NCHARG(NUP) == NCHARG(LOW)) THEN
            STOP 'FATAL ERROR in PREPKUBAT (KON)'
          ENDIF
          EDGE = ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
          EDGELAM=1.E8/EDGE
          
          ABFINT1 = 0.
          ABFINT2 = 0.
          ABFINT2d = 0.
          KONEDGE = ISRCHFGT(NF,XLAMBDA,1,EDGELAM)-1
          DO K = 1, KONEDGE
            W = 1.E8/XLAMBDA(K)
            W3 = W*W*W
            ABF = SIGMAKI(K, KON)
C            XJCLK = MAX(0., XJC(L,K))
            XJCLK = XJC(L,K)
            IF (W < EDGE .AND. ABS(ABF) > 0.) THEN
              STOP 'FATAL ERROR in PREPKUBAT: SIGMAKI should be Zero!'
            ENDIF
            
            ABFINT1 = ABFINT1 
     >           + FWEIGHT(K) * ABF * XJCLK
     >               * (1. - EDGE / W)

            ABFINT2 = ABFINT2
     >           + FWEIGHT(K) * ABF * EXP(-C1*W/TL) *
     >             (XJCLK  + C2 * W3) *
     >             (1. -  EDGE / W)

            ABFINT2d = ABFINT2d
     >           + FWEIGHT(K) * ABF * EXP(-C1*W/TL) *
     >             (XJCLK  + C2 * W3) *
     >             (1. -  EDGE / W) * W
          ENDDO

          POPSUM1 = POPSUM1 + ENTOTL*POPLOW * ABFINT1
          PADD2 = ENTOTL*POPNUP
     >                     * ENLTE(LOW)/ENLTE(NUP) 
          POPSUM2 = POPSUM2 + PADD2 * ABFINT2
          POPSUM2d = POPSUM2d + C1 / (TL*TL) * PADD2 * ABFINT2d
     >                        - 3. / (2.*TL) * PADD2 * ABFINT2
     >                 - EDGE * C1 / (TL*TL) * PADD2 * ABFINT2
          
        ENDDO
        
C***  K-SHELL IONISATION (not carefully tested yet!)  ******************************
        IF (bKSHELL) THEN
          DO J=1, N 
            NOMJ = NOM(J)
            ISTATE = NCHARG(J) + 1
            IF (SIGMATHK(NOMJ,ISTATE) == 0.) CYCLE
            WK = 1.E8 / EDGEK(NOMJ,ISTATE)
            ABFINT1 = 0.
            ABFINT2 = 0.
            ABFINT2d = 0.
            DO K=1, NF
              IF (XLAMBDA(K) > WK) EXIT  !Stop if radiation is too soft for K-SHELL ionization
              W = 1.E8 / XLAMBDA(K)
              XJCLK = MAX(0., XJC(L,K))
              CALL KSIGMA (SIGMAK, SIGMATHK(NOMJ,ISTATE), 
     >                     EDGEK(NOMJ,ISTATE), W, SEXPOK(NOMJ,ISTATE))
              ABFINT1 = ABFINT1 
     >             + FWEIGHT(K) * SIGMAK * XJCLK
     >                 * (1. - WK / W)

              ABFINT2 = ABFINT2
     >             + FWEIGHT(K) * SIGMAK * EXP(-C1*W/TL) *
     >                 (XJCLK  + C2 * W3) *
     >                 (1. -  WK / W)
              ABFINT2d = ABFINT2d
     >           + FWEIGHT(K) * SIGMAK * EXP(-C1*W/TL) *
     >             (XJCLK  + C2 * W3) *
     >             (1. -  WK / W) * W
            ENDDO
            POPSUM1 = POPSUM1 + ENTOTL*POPNUM(L,LOW) * ABFINT1
            PADD2 = ENTOTL*POPNUM(L,NUP) 
     >                     * ENLTE(LOW)/ENLTE(NUP) 
            POPSUM2 = POPSUM2 + PADD2 * ABFINT2
            POPSUM2d = POPSUM2d + C1 / (TL*TL) *
     >         (ABFINT2d - WK * ABFINT2) * PADD2
          ENDDO
        ENDIF

        QBF(L,1) = PI4 * POPSUM1
        QBF(L,2) = PI4 * POPSUM2
        QBFDT(L,1) = 0.
        QBFDT(L,2) = PI4 * POPSUM2d 
        
        POPSUM1 = 0.
        POPSUM1d = 0.
        POPSUM2 = 0.
        POPSUM2d = 0.
C****** COLLISIONS *****************************************************
C***    In order to estimate the temperature derivative of the
C***    collisional rates, the subr. COLLI is called twice, the first time with
C***    a temperature being enhanced by 1%
        TLmod = 1.01*TL
        CALL COLLI (NDIM,N,ENLTE,TLmod,ENE,NCHARG,ELEVEL,EINST,DCRATEDT,
     >          EION,COCO,KEYCBB,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     >          INDNUP, INDLOW, LASTIND, LASTINDAUTO, LASTINDALL, 
     >          KONTNUP, KONTLOW, LASTKON, KEYCBF, IONGRND, MAXATOM)
        CALL COLLI (NDIM,N,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     >          EION,COCO,KEYCBB,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     >          INDNUP, INDLOW, LASTIND, LASTINDAUTO, LASTINDALL, 
     >          KONTNUP, KONTLOW, LASTKON, KEYCBF, IONGRND, MAXATOM)
        DO I=1, N
          DO J=1, N
            !calculate temperature derivative
            DCRATEDT(I,J) = (DCRATEDT(I,J) - CRATE(I,J)) / (TLmod - TL)
          ENDDO
        ENDDO
     
        bCOLLIDONE = .FALSE.
C***    Note: We need the pure OMEGA without boltzmann factor, therefore
C             always those CRATE is used which does not contain it. This
C             differs between line transitions and ionizations

C***    line transition collisions
        DO IND=1, LASTINDALL
C***       Special treatment required for DRTRANSITS - see RADNET
ccc        currently just neglected
           IF (IND .GT. LASTIND .AND. IND .LE. LASTINDAUTO) CYCLE
            NUP=INDNUP(IND)
            LOW=INDLOW(IND)
            W=ELEVEL(NUP)-ELEVEL(LOW)
            IF (NCHARG(NUP) /= NCHARG(LOW)) THEN
              STOP 'FATAL ERROR in PREPKUBAT (COLLI - IND)'
            ENDIF
            ETRANSIT = W * CLIGHT * HPLANCK
            
            PADD1 = ENTOTL * POPNUM(L,NUP) * ETRANSIT                  
            POPSUM1 = POPSUM1 + PADD1 * CRATE(NUP,LOW)
            POPSUM1d = POPSUM1d + PADD1 * DCRATEDT(NUP,LOW)
C            PADD1 = ENTOTL * POPNUM(L,NUP) * ENLTE(LOW)/ENLTE(NUP)
C     >                          * CRATE(LOW,NUP) * ETRANSIT
C            POPSUM1 = POPSUM1 + PADD1             
            
            
            PADD2 = ENTOTL * POPNUM(L,LOW) * ETRANSIT
            POPSUM2 = POPSUM2 + PADD2 * CRATE(LOW,NUP)
            POPSUM2d = POPSUM2d 
     >         + PADD2 * ENLTE(NUP)/ENLTE(LOW) * DCRATEDT(NUP,LOW)
     >         + PADD2 * (C1*W/(TL*TL)) * CRATE(LOW,NUP)
C            PADD2 = ENTOTL * POPNUM(L,LOW) * CRATE(LOW,NUP) * ETRANSIT
C            POPSUM2 = POPSUM2 + PADD2
C            POPSUM2d = POPSUM2d + (C1*W/(TL*TL)) * PADD2
            bCOLLIDONE(NUP,LOW) = .TRUE.
            bCOLLIDONE(LOW,NUP) = .TRUE.
        ENDDO        
C***    Additional iron line transition collisions which are not in the index numbering
C***    (due zero radiative cross sections, see PoWR-Memo 20100826.txt) 
        DO LOW=1, N  
          IF (NOM(LOW) == KODAT(26)) THEN 
            DO NUP=LOW+1, N
C***          Within the same ion?
              IF (NCHARG(LOW) == NCHARG(NUP)) THEN                  
C***            Only if transition has not been "used" in the collision part already
                IF (.NOT. bCOLLIDONE(LOW,NUP)) THEN
                  W=ELEVEL(NUP)-ELEVEL(LOW)
                  ETRANSIT = W * CLIGHT * HPLANCK
                  PADD1 = ENTOTL * POPNUM(L,NUP) * ETRANSIT
                  POPSUM1 = POPSUM1 + PADD1 * CRATE(NUP,LOW) 
                  POPSUM1d = POPSUM1d + PADD1 * DCRATEDT(NUP,LOW)
C                  PADD1 = ENTOTL * POPNUM(L,NUP) 
C     >                          * CRATE(NUP,LOW) * ETRANSIT
C                  PADD1 = ENTOTL * POPNUM(L,NUP) *ENLTE(LOW)/ENLTE(NUP)
C     >                          * CRATE(LOW,NUP) * ETRANSIT
C                  POPSUM1 = POPSUM1 + PADD1
C                  POPSUM1d = POPSUM1d + (1./2./TL) * PADD1 
C     >              - ENTOTL * POPNUM(L,NUP) * ETRANSIT * ENE * COM 
C     >                * C1 * W / TL**(3./2.) * WEIGHT(LOW)/WEIGHT(NUP)
C     >                * ENLTE(LOW)/ENLTE(NUP)

                  PADD2 = ENTOTL * POPNUM(L,LOW) * ETRANSIT
                  POPSUM2 = POPSUM2 + PADD2 * CRATE(LOW,NUP) 
                  POPSUM2d = POPSUM2d 
     >              + PADD2 * ENLTE(NUP)/ENLTE(LOW) * DCRATEDT(NUP,LOW)
     >              + PADD2 * (C1*W/(TL*TL)) * CRATE(LOW,NUP)
C                  PADD2 = ENTOTL * POPNUM(L,LOW)            
C     >                          * CRATE(LOW,NUP) * ETRANSIT
C                  PADD2 = ENTOTL * POPNUM(L,LOW)            
C     >                          * CRATE(LOW,NUP) * ETRANSIT
C                  POPSUM2 = POPSUM2 + PADD2
C                  POPSUM2d = POPSUM2d + (1./2./TL+C1*W/(TL*TL)) * PADD2 
C     >              - ENTOTL * POPNUM(L,LOW) * ETRANSIT * ENE * COM 
C     >                * C1 * W / TL**(3./2.) * WEIGHT(LOW)/WEIGHT(NUP) 

                  bCOLLIDONE(NUP,LOW) = .TRUE.
                  bCOLLIDONE(LOW,NUP) = .TRUE.
                ENDIF
              ENDIF              
            ENDDO
          ENDIF
        ENDDO        

C***    collisonal ionization        
        DO KON=1, LASTKON
            NUP=KONTNUP(KON)
            LOW=KONTLOW(KON)
            POPLOW = POPNUM(L,LOW)
            POPNUP = POPNUM(L,NUP) 
            IF (ZERO_RATES(LOW,L)) POPLOW = 0.
            IF (ZERO_RATES(NUP,L)) THEN
              POPNUP = 0.
              POPLOW = 0.
            ENDIF
            IF (NCHARG(NUP) == NCHARG(LOW)) THEN
              STOP 'FATAL ERROR in PREPKUBAT (COLLI - KON)'
            ENDIF
            EDGE = ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
C            EDGE = ELEVEL(NUP)-ELEVEL(LOW)
            
C            ETRANSIT = EDGE * HPLANCK
            ETRANSIT = EDGE * CLIGHT * HPLANCK
            
            PADD1 = ENTOTL * POPNUP * ETRANSIT
            POPSUM1 = POPSUM1 + PADD1 * CRATE(NUP,LOW)
            POPSUM1d = POPSUM1d 
     >        + PADD1 * ENLTE(LOW)/ENLTE(NUP) * DCRATEDT(LOW,NUP)
     >        - PADD1 * (3./(2.*TL) + C1*EDGE/(TL*TL)) * CRATE(NUP,LOW)
            
C            PADD1 = ENTOTL * POPNUP * CRATE(NUP,LOW) * ETRANSIT
C            PADD1 = ENTOTL * POPNUP * CRATE(LOW,NUP) 
C     >                        * ENLTE(LOW)/ENLTE(NUP) * ETRANSIT
C            POPSUM1 = POPSUM1 + PADD1
C            POPSUM1d = POPSUM1d - PADD1 * (1./TL)

            PADD2 = ENTOTL * POPLOW * ETRANSIT
            POPSUM2 = POPSUM2 + PADD2 * CRATE(LOW,NUP)
            POPSUM2d = POPSUM2d + PADD2 * DCRATEDT(LOW,NUP)

C            PADD2 = ENTOTL * POPLOW * CRATE(LOW,NUP) * ETRANSIT
C            POPSUM2 = POPSUM2 + PADD2
C            POPSUM2d = POPSUM2d + PADD2 * (1./(2.*TL) + C1*EDGE/(TL*TL))
C            POPSUM2d = POPSUM2d + PADD2 * (1./TL)
            bCOLLIDONE(NUP,LOW) = .TRUE.
            bCOLLIDONE(LOW,NUP) = .TRUE.
        ENDDO        
        
        
C****   All Qs are in erg / cm**3 / s  (energy density)
        QC(L,1) = POPSUM1
        QC(L,2) = POPSUM2
C        QCDT(L,1) = 0.      !ignores T-dependence of OMEGA
C        QCDT(L,2) = 0.
        QCDT(L,1) = POPSUM1d 
        QCDT(L,2) = POPSUM2d 
C        QCDT(L,1) = POPSUM1d * ENE / 2. * COM / SQRT(TL)
C        QCDT(L,2) = POPSUM2d * ENE / 2. * COM / SQRT(TL)

        DELTAQ(L) = ( QFF(L,1) + QBF(L,1) + QC(L,1) 
     >              - QFF(L,2) - QBF(L,2) - QC(L,2) ) 
        
        !simple approach
C        DTKUBAT(L) = ( QFF(L,1) + QBF(L,1) + QC(L,1) 
C     >               - QFF(L,2) - QBF(L,2) - QC(L,2) ) 
C     >               *4./3.* C1/C2 / 1.E3
        !Planck comparison approach
C        DTKUBAT(L) = PI4*DELTAQ(L) / (16.*STEBOL* OPAROSS(L)* TL**3)

        OPAMEAN = OPASMEANTC(L) - OPALAMBDAMEAN(L)
        DTB = PI4 / 
     >       (16. * TL**(3.) * STEBOL * OPAMEAN * 1000.)


        !Newton-Raphson approach
        DELTAQDT(L) = ( QFFDT(L,1) + QBFDT(L,1) + QCDT(L,1) 
     >                - QFFDT(L,2) - QBFDT(L,2) - QCDT(L,2) )
        DTKUBAT(L) = - DELTAQ(L) / DELTAQDT(L) / 1000.      !transform to kK for plots
        
      ENDDO
            
      IF (bPrintQTable) THEN
        WRITE (hCPR,*)
        WRITE (hCPR,*) 'TESTAUSGABE ZUR KUBAT-METHODE:'
        WRITE (hCPR,'(5X,4(3(A10,2X),3X),A7)') 
     >    '   QFF_H  ', '   QFF_C  ', ' DELTA QFF', 
     >    '   QBF_H  ', '   QBF_C  ', ' DELTA QBF', 
     >    '   QC_H   ', '   QC_C   ', ' DELTA QC ', 
     >    '  DELTA Q ', '  DDQDT   ', '   N-R    ', 'Planck '
        DO L=1, ND
          OPAMEAN = OPASMEANTC(L)
c          DTB = PI4 / 
c     >       (16. * TOLD(L)**(3.) * STEBOL * OPAMEAN * 1000.)
          IF (L > 1. .AND. L < ND) THEN
            DTB = RSTAR * 0.5 * (RADIUS(L-1) - RADIUS(L+1)) / 
     >         (4. * TOLD(L)**(3.) * STEBOLDPI * 1000.)
          ELSEIF (L==ND) THEN 
            DTB = RSTAR * 0.5 * (RADIUS(ND-1) - RADIUS(ND)) / 
     >         (4. * TOLD(L)**(3.) * STEBOLDPI * 1000.)            
          ELSE
            DTB = RSTAR * 0.5 * (RADIUS(1) - RADIUS(2)) / 
     >         (4. * TOLD(L)**(3.) * STEBOLDPI * 1000.)            
          ENDIF
          WRITE (hCPR,FMT='(I3,3(3(2X,G10.3),2X,"#"),4(2X,G10.3))') 
     >      L, QFF(L,1), QFF(L,2), QFF(L,1)- QFF(L,2), 
     >      QBF(L,1), QBF(L,2), QBF(L,1) - QBF(L,2),
     >      QC(L,1), QC(L,2), QC(L,1)-QC(L,2),      
     >      DELTAQ(L),  DELTAQDT(L), 
     >      - DELTAQ(L) / DELTAQDT(L) / 1000.,
C     >    (2./3.*(QFF(L,2)+QBF(L,2)+QC(L,2))/BOLTZ/ENE)/1000.
     >      DELTAQ(L) * DTB
        ENDDO
        WRITE (hCPR,*) 'maxmin: ', MAXVAL(DTKUBAT), MINVAL(DTKUBAT)
        WRITE (hCPR,*) 'TESTAUSGABE ZUR KUBAT-METHODE (ABLEITUNGEN):'
        DO L=1, ND
          WRITE (hCPR,FMT='(I3,3(3(2X,G10.3),2X,"#"),3(2X,G10.3))') 
     >      L, QFFDT(L,1), QFFDT(L,2), QFFDT(L,1)- QFFDT(L,2), 
     >      QBFDT(L,1), QBFDT(L,2), QBFDT(L,1) - QBFDT(L,2),
     >      QCDT(L,1), QCDT(L,2), QCDT(L,1)-QCDT(L,2),
     >      QFFDT(L,1) + QBFDT(L,1) + QCDT(L,1), 
     >      QFFDT(L,2) + QBFDT(L,2) + QCDT(L,2),
     >      QFFDT(L,1) + QBFDT(L,1) + QCDT(L,1) - 
     >       ( QFFDT(L,2) + QBFDT(L,2) + QCDT(L,2) )
        ENDDO
      ENDIF
      
      
      RETURN
      END
      
