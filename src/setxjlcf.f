      SUBROUTINE SETXJLCF (LASTIND,INDLOW,INDNUP,XRED,XBLUE,OPACIND,
     $   SCNEIND,SCOLIND,SLNEW,SLOLD,OPAL,XJLAPP,
     $   NF,XLAMBDA,SCNEW,OPAC,OPALOLD,ITNEL,LASTINDAUTO,
     >   NFL,PHI,PWEIGHT,NDIM,EINST,ELEVEL,EN,WEIGHT,ND,XJL,ENTOTL,
     >   RSTAR,VDOPUNIT,DELTAX,XMAX,L,TL,TLOLD,NOTEMP,
     $   IBLENDS,MAXLAP,XLAMZERO,BETA,PHIL,NBLENDS,KRUDAUT,MAXAUTO, 
     >   WFELOW,WFENUP,EN1,BDIAG)
C***********************************************************************
C***  CALLED FROM: SUBROUTINE COMA
C***  CALCULATE LINE RADIATION FIELD WITH APPROXIMATE LAMBDA OPERATOR TERMS
C***  USING THE "CORE FRACTION" METHOD FROM HAMANN (1985, 1986)
C***  THIS SUBROUTINE ALSO PROVIDES  (BY CALLING LIOP) THE LINE OPACITIES
C***  AND THE NEW LINE SOURCE FUNCTION
C***  ATTENTION: XJLAPP AND XJL ARE NOT DEFINED FOR RUDIMENTAL LINES!
C***          -  OPAL AND SLNEW ARE NOT DEFINED FOR RUDIMENTAL LINES
C***              AND FOR ZERO SCHARMER CORE (IF NO TEMP. CORRECTIONS)
C***          - FN USES THE SCRATCH ARRAY ATEST
C***********************************************************************
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: L, NF, NFL, ITNEL, MAXLAP, LASTIND, 
     >                       LASTINDAUTO, ND, NDIM, MAXAUTO
      REAL, INTENT(IN) :: VDOPUNIT, DELTAX, ENTOTL, RSTAR, TL, TLOLD, 
     >                    XMAX
 
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NDIM) :: ELEVEL, EN, EN1, WEIGHT
      REAL, DIMENSION(ND, LASTINDAUTO) :: XJL
      REAL, DIMENSION(NF) :: XLAMBDA, SCNEW, OPAC
      REAL, DIMENSION(NFL) :: PHI, PWEIGHT
      REAL, DIMENSION(LASTIND) :: XRED, XBLUE, XLAMZERO, 
     >                            SCOLIND, SCNEIND, SLNEW, SLOLD,
     >                            OPAL, OPALOLD, OPACIND, XJLAPP
      INTEGER, DIMENSION(LASTIND) :: INDNUP, INDLOW, NBLENDS
      INTEGER, DIMENSION(MAXLAP,LASTIND) :: IBLENDS
      INTEGER, DIMENSION(MAXAUTO) :: KRUDAUT
      REAL, DIMENSION(MAXLAP) :: BETA
      REAL, DIMENSION(NFL,MAXLAP) :: PHIL
      REAL, DIMENSION(ND,LASTIND) :: WFELOW, WFENUP

      LOGICAL :: NOTEMP, OLDOPAL, BNOCORE
      LOGICAL, DIMENSION(LASTIND) :: BDIAG

      INTEGER, PARAMETER :: NFLMAX = 50
      REAL, DIMENSION(NFLMAX) :: FN

      
      REAL, EXTERNAL :: BNUE
      REAL :: XLAM, DPLOW, DPNUP, XKL, XREL, XSHIFT, 
     >        DELTASC, DELTASL, XKLMIN, XKLMAX, Q, P, FCC, FCL,
     >        BETA1, WAVENUM, ETAL
      
      INTEGER :: IND, IND1, NUP, LOW, LB, KL, KLA, KLB, KLMIN, KLMAX

      
C***  WPI = SQUARE ROOT OF PI
      REAL, PARAMETER :: WPI    = 1.772453851
C***  CLIGHT = VELOCITY OF LIGHT IN KM/SEC
      REAL, PARAMETER :: CLIGHT = 2.99792458E5

      IF (NFL .GT. NFLMAX) THEN
        WRITE (0,*) 'ERROR: NFL > NFLMAX ', NFL, NFLMAX
        STOP 'FATAL ERROR IN SETXJLCF'
      ENDIF

      OLDOPAL = (L.GT.1) .AND. (ITNEL.EQ.1)

C***  LOOP OVER ALL LINE INDICES TO PREPARE OPACITIES, ETC. ---------
      DO 1 IND=1,LASTIND
      IF (XRED(IND) .GE. XBLUE(IND)) THEN
        BNOCORE=.TRUE.
      ELSE
        BNOCORE=.FALSE.
      ENDIF

      LOW=INDLOW(IND)
      NUP=INDNUP(IND)

C***  FOR RUDIMENTAL LINES: OPAL, SLNEW    NOT DEFINED
      IF (EINST(LOW,NUP) .EQ. -2.) THEN
         OPAL(IND)=TRANSFER('UNDEF', OPAL(IND))
         SLNEW(IND)=TRANSFER('UNDEF', SLNEW(IND))
         GOTO 1
      ENDIF

      XJLAPP(IND)=XJL(L,IND)

C***  OPAL, ETAL, SLNEW NOT CALCULATED IF ZERO LINE CORE AND "NOTEMP"
C***  NOTE THAT IRON LINES are always like lines with "core"
      IF (BNOCORE  .AND.  NOTEMP .AND. .NOT. BDIAG(IND))
     >    THEN
         OPAL(IND)=TRANSFER('UNDEF', OPAL(IND))
         SLNEW(IND)=TRANSFER('UNDEF', SLNEW(IND))
         GOTO 1
      ENDIF

      XLAM=XLAMZERO(IND)

      IF (OLDOPAL) THEN
         OPAL(IND) = OPALOLD(IND)
      ELSE
         CALL LIOP (EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     >         1,XLAM, ENTOTL, EN,RSTAR,OPAL(IND),ETAL,VDOPUNIT)
      ENDIF

C***  LASER SECURITY
      IF (OPAL(IND) .LE. .0) THEN
            SLNEW(IND)=.0
            XRED(IND)=.0
            XBLUE(IND)=.0
            OPAL(IND)=.0
            GOTO 1
            ENDIF

      IF (OLDOPAL) THEN
         SLNEW(IND) = SLOLD(IND)
         ELSE
         SLNEW(IND)=ETAL/OPAL(IND)
         ENDIF

      IF (BNOCORE) GOTO 1

      
C *** OUTPUT WARNING
C      IF (BDIAG(IND)) 
C     > write(0,*)'WARNING: Backgroundvalues calculated for Iron'
  
C***  BACKGROUND CONTINUUM OPACITY AND SOURCE FUNCTION ARE OBTAINED
C***  BY INTERPOLATION AT THE LINE FREQUENCIES
      CALL LIPO (OPACIND(IND),XLAM,OPAC,XLAMBDA,NF)
      IF (OLDOPAL) THEN
         SCNEIND(IND) = SCOLIND(IND)
      ELSE
         WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
         CALL XRUDI (SCNEIND(IND),WAVENUM,SCNEW,XLAMBDA,1,NF,1)
      ENDIF

   1  CONTINUE

C***  ENDLOOP  FOR PREPARATION -------------------------------------------

      IF (OLDOPAL) GOTO 102

C***  LOOP OVER ALL LINE INDICES TO CALCULATE XJLAPP ---------
      DO 101 IND=1,LASTIND
      IF (XRED(IND) .GE. XBLUE(IND)) THEN
        BNOCORE=.TRUE.
      ELSE
        BNOCORE=.FALSE.
      ENDIF

      LOW=INDLOW(IND)
      NUP=INDNUP(IND)

C**   NO LINE CORE -> XJLAPP = XJL REMAINS UNCHANGED
C***  IRON: EXCEPT FOR IRON-LINES
      IF (BNOCORE .AND. .NOT. BDIAG(IND)) GOTO 11

C***  IRON LINES
      IF (BDIAG(IND)) THEN
         DPLOW = EN(LOW) - EN1(LOW)
         DPNUP = EN(NUP) - EN1(NUP)
         XJLAPP(IND) = XJL(L, IND) +
     >                  WFELOW(L, IND)*DPLOW + WFENUP(L, IND)*DPNUP
         GOTO 101
      ENDIF
C *** OUTPUT WARNING
      IF (BDIAG(IND)) 
     > write(0,*)'WARNING: Iron zu weit'


C***  RUDIMENTAL LINES ARE NOT CALCULATED -> XJLAPP IS UNDEFINED
      IF (EINST(LOW,NUP) .EQ. -2.) GOTO 101

      XLAM=XLAMZERO(IND)

C***  LINE CORE INTEGRAL

C***  CONVERT XRED AND XBLUE TO CORRESPONDING LINE FREQUENCY INDICES
C***  NOTE THAT THE LINE FREQUENCIES ARE INDEXED IN FALLING SEQUENCE
      KLMAX=(XMAX-XRED(IND))/DELTAX + 1.
      KLMIN=(XMAX-XBLUE(IND))/DELTAX + 1.999999999
      KLA = MAX (KLMIN-1,  1 )
      KLB = MIN (KLMAX+1, NFL)

C***  DO OLD VERSION, IF CURRENT LINE HAS NO OVERLAPS
      IF (NBLENDS(IND) .EQ. 1) THEN 

      BETA1 = OPACIND(IND)/OPAL(IND)
      FCL=.0
      FCC=.0
      DO 16 KL=KLMIN, KLMAX
      P = 1. / (1. + BETA1 / PHI(KL))
      FCL = FCL+P*PWEIGHT(KL)
      FCC = FCC+(1.-P)*PWEIGHT(KL)
   16 CONTINUE
C***  CORRECT FOR THE INTEGRATION STEP XRED,X(KLMAX)
      IF (KLMAX .LT. NFL) THEN
         XKLMAX=XMAX-(KLMAX-1)*DELTAX
         Q=(XKLMAX-XRED(IND))/DELTAX
         P=1. / (1. + BETA1 / PHI(KLMAX))
         FCL=FCL+    P *PWEIGHT(KLMAX)*0.5*(Q-1.)
         FCC=FCC+(1.-P)*PWEIGHT(KLMAX)*0.5*(Q-1.)
         P=1. / (1. + BETA1 / PHI(KLMAX+1))
         FCL=FCL+    P *PWEIGHT(KLMAX+1)*0.5*Q
         FCC=FCC+(1.-P)*PWEIGHT(KLMAX+1)*0.5*Q
         ENDIF
C***  CORRECT FOR THE INTEGRATION STEP X(KLMIN),XBLUE
      IF (KLMIN .GT. 1) THEN
         XKLMIN=XMAX-(KLMIN-1)*DELTAX
         Q=(XBLUE(IND)-XKLMIN)/DELTAX
         P=1. / (1. + BETA1 / PHI(KLMIN))
         FCL=FCL+    P *PWEIGHT(KLMIN)*0.5*(Q-1.)
         FCC=FCC+(1.-P)*PWEIGHT(KLMIN)*0.5*(Q-1.)
         P=1. / (1. + BETA1 / PHI(KLMIN-1))
         FCL=FCL+    P *PWEIGHT(KLMIN-1)*0.5*Q
         FCC=FCC+(1.-P)*PWEIGHT(KLMIN-1)*0.5*Q
         ENDIF
      DELTASC=SCNEIND(IND)-SCOLIND(IND)
      DELTASL=SLNEW(IND)-SLOLD(IND)
      XJLAPP(IND)=XJLAPP(IND)+FCL*DELTASL+FCC*DELTASC

      ELSE

C***  THIS BRANCH HANDLES OVERLAP TREATMENT ---------------------------------

C***  Prepare FOR ALL OVERLAPPING LINES: 
C***  - Beta = LINE / CONTINUUM OPACITY (INVERSELY TO PREVIOUS DEFINITION!)
      DO 210 LB=1, NBLENDS(IND)
       IND1 = IBLENDS(LB,IND)
       IF (XRED(IND1) .GE. XBLUE(IND1)) GOTO 210
       BETA(LB) = OPAL(IND1) / OPACIND(IND)
C***  - PHIL = DISPLACED PROFILE FUNCTIONS
      IF (IND .EQ. IND1) THEN
         DO 250 KL=KLA, KLB
         PHIL(KL,LB) = PHI(KL)
  250    CONTINUE
         ELSE
C***     SHIFT OF CONSIDERED LINE RELATIVE TO BLENDING LINE, IN DOPPLER UNITS
         XSHIFT = (XLAMZERO(IND1)/XLAMZERO(IND) - 1. ) * CLIGHT / VDOPUNIT
         DO 220 KL=KLA, KLB
         XKL = XMAX - (KL-1) * DELTAX
         XREL = XKL + XSHIFT 
         PHIL(KL,LB) = EXP(-XREL*XREL) / WPI
  220    CONTINUE
         ENDIF
  210 CONTINUE


C***  STORE FACTOR 1/N  (INTEGRAND IN FCC) 
C***    N = 1 + SUM OVER (PHI(LB) * BETA(LB)  
C***    AS FUNCTION OF FREQUENCY
      DO 230 KL=KLA, KLB
        FN(KL)=1.
        DO 240 LB=1, NBLENDS(IND)
          IND1 = IBLENDS(LB,IND)
          IF (XRED(IND1) .GE. XBLUE(IND1)) GOTO 240
          FN(KL) = FN(KL) + PHIL(KL,LB) * BETA(LB)
  240     CONTINUE
        FN(KL)=1./FN(KL)
  230   CONTINUE


C***  INTEGRATE COEFFICIENT OF THE CONTINUUM SOURCE FUNCTION--> FCC
      FCC=0.0
      DO 26 KL=KLMIN, KLMAX
      FCC = FCC + FN(KL) * PWEIGHT(KL)
   26 CONTINUE
C***  CORRECT FOR THE INTEGRATION STEP XRED,X(KLMAX)
      IF (KLMAX .LT. NFL) THEN
         XKLMAX=XMAX-(KLMAX-1)*DELTAX
         Q=(XKLMAX-XRED(IND))/DELTAX
         FCC = FCC + FN(KLMAX)   * PWEIGHT(KLMAX)   * 0.5 * (Q-1.)
     >             + FN(KLMAX+1) * PWEIGHT(KLMAX+1) * 0.5 * Q
         ENDIF
C***  CORRECT FOR THE INTEGRATION STEP X(KLMIN),XBLUE
      IF (KLMIN .GT. 1) THEN
         XKLMIN=XMAX-(KLMIN-1)*DELTAX
         Q=(XBLUE(IND)-XKLMIN)/DELTAX
         FCC=FCC + FN(KLMIN)   * PWEIGHT(KLMIN)   * 0.5 * (Q-1.)
     >           + FN(KLMIN-1) * PWEIGHT(KLMIN-1) * 0.5 * Q
         ENDIF
      DELTASC = SCNEIND(IND)-SCOLIND(IND)
      XJLAPP(IND) = XJLAPP(IND) + FCC*DELTASC

C***  INTEGRATE OVER LINES
      DO 28 LB=1, NBLENDS(IND)
      IND1 = IBLENDS(LB,IND)
      IF (XRED(IND1) .GE. XBLUE(IND1)) GOTO 28
      FCL=0.0
        DO 27 KL=KLMIN,KLMAX
         FCL=FCL + PHIL(KL,LB) * FN(KL) * PWEIGHT(KL)
   27    CONTINUE
C***  CORRECT FOR THE INTEGRATION STEP XRED,X(KLMAX)
        IF (KLMAX .LT. NFL) THEN
           XKLMAX=XMAX-(KLMAX-1)*DELTAX
           Q=(XKLMAX-XRED(IND))/DELTAX
           FCL=FCL+ PHIL(KLMAX,LB)   * FN(KLMAX)   * PWEIGHT(KLMAX)
     >                                             * 0.5 * (Q-1.)
     >            + PHIL(KLMAX+1,LB) * FN(KLMAX+1) * PWEIGHT(KLMAX+1)
     >                                             * 0.5 * Q
           ENDIF
C***  CORRECT FOR THE INTEGRATION STEP X(KLMIN),XBLUE
        IF (KLMIN .GT. 1) THEN
           XKLMIN=XMAX-(KLMIN-1)*DELTAX
           Q=(XBLUE(IND)-XKLMIN)/DELTAX
           FCL= FCL + PHIL(KLMIN,LB)   * FN(KLMIN)   * PWEIGHT(KLMIN)
     >                                               * 0.5 * (Q-1.)
     >              + PHIL(KLMIN-1,LB) * FN(KLMIN-1) * PWEIGHT(KLMIN-1)
     >                                               * 0.5 * Q
           ENDIF

C***  Sum up XJLAPP: ADD THE CONTRIBUTION OF BLENDING LINE IND1
      DELTASL = SLNEW(IND1) - SLOLD(IND1)
      XJLAPP(IND) = XJLAPP(IND) + FCL * BETA(LB) * DELTASL

28    CONTINUE

      ENDIF


   11 CONTINUE

C***  ONLY IF TEMPERATURE CORRECTIONS ARE APPLIED:
C***  ADD SPECIAL TERM AT INNER BOUNDARY, WHICH ACCOUNTS FOR THE
C***    DIRECT TEMPERATURE-DEPENDENCE OF THE RADIATION FIELD
C***    VIA THE BOUNDARY CONDITION
      IF (.NOT. NOTEMP  .AND.  L .EQ. ND) THEN 
         XJLAPP(IND)=XJLAPP(IND)+0.5*(BNUE(XLAM,TL)-BNUE(XLAM,TLOLD))
         ENDIF

 101    CONTINUE
C***  ENDLOOP  ---------------------------------------------------------

  102 CONTINUE

C***  XJL --> XJLAPP FOR DR-TRANSITIONS
      DO 10 IND = LASTIND+1, LASTINDAUTO
         IF (KRUDAUT(IND-LASTIND) .EQ. 1) GOTO 10
         XJLAPP(IND) = XJL(L,IND)
   10 CONTINUE

      RETURN
      END
