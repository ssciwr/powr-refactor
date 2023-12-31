      SUBROUTINE DLIOP (I,ENTOTL,DOPAL,DETAL,XRED,XBLUE,VDOPUNIT,RSTAR,
     >       N,NDIM,EINST,WEIGHT,ELEVEL,LASTIND,INDLOW,INDNUP, EN, 
     >       LASTINDAUTO, bLAMAPPCOLI, bUSEALO, iBLOCKINVERSION)
C*******************************************************************************
C***  DERIVATIVES OF LINE OPACITIES AND EMISSIVITIES WITH RESPECT TO EN(I)
C***  AT CURRENT DEPTH POINT
C***  CALLED FROM: SUBROUTINE COMA
C***  ATTENTION: IF (NOTEMP), THESE DERIVATIVES ARE ONLY CALCULATED FOR
C***             NON-ZERO SCHARMER CORES!
C***           - DERIVATIVES NOT CALCULATED FOR RUDIMENTAL LINES!
C*******************************************************************************
 
      DIMENSION EINST(NDIM,NDIM),WEIGHT(NDIM),ELEVEL(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION DETAL(2),DOPAL(2),XRED(2),XBLUE(2), EN(2)
      LOGICAL, DIMENSION(LASTINDAUTO) :: bUSEALO
      LOGICAL :: bLAMAPPCOLI
      INTEGER :: iBLOCKINVERSION

C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  PI8 = 8*PI
      DATA PI8 /25.1327412288 /

C***  FOR EXACT COMPATIBILITY, CALCULATE C3 = H*C/(4*PI) :
      C3=C2/PI8

C***  DERIVATIVES WITH RESPECT TO T OR EL.DENSITY ARE ZERO:
      IF (I .EQ. N+1  .OR.  I .EQ. N+2) THEN
         DO 2 IND=1,LASTIND
         DETAL(IND)=0.
         DOPAL(IND)=0.
    2    CONTINUE
         RETURN
         ENDIF
  
C***  LOOP OVER ALL LINE TRANSITIONS  ----------------------------------
      DO 1 IND=1,LASTIND
      LOW=INDLOW(IND)
      NUP=INDNUP(IND)

C***  NO DERIVATIVES CALCULATED, IF NOT NECESSARY:
      IF (((.NOT. bLAMAPPCOLI) .AND. 
     >     ((XRED(IND) .GE. XBLUE(IND)) .OR.
     >      EINST(LOW,NUP) .EQ. -2. ))
     >   .OR. (bLAMAPPCOLI .AND. 
ccc   No derivates needed in case of rudimental lines (or sometimes in XJLCF method)
     >   (.NOT. bUSEALO(IND) .OR. EINST(LOW,NUP) .EQ. -2. ) ) ) THEN
         DETAL(IND)=TRANSFER('UNDEF', DETAL(IND))
         DOPAL(IND)=TRANSFER('UNDEF', DOPAL(IND))
         GOTO 1
      ENDIF

C***  ZERO DERIVATIVES, IF TRANSISTION LEVELS ARE NOT CONCERNED:
      IF (I .NE. LOW .AND. I .NE. NUP) THEN
         DETAL(IND)=.0
         DOPAL(IND)=.0
         GOTO 1
         ENDIF

      XLAMCM=1./(ELEVEL(NUP)-ELEVEL(LOW))
C***  DND= DELTA-NUE-DOPPLER (HERTZ)
      DND=VDOPUNIT*1.E5/XLAMCM
      EMINDU=EINST(NUP,LOW)*XLAMCM*XLAMCM/PI8*RSTAR

      IF (I .EQ. LOW) THEN
C***     DERIVATIVE WITH RESPECT TO LOWER LEVEL
         DETAL(IND)=.0
         ABSORP=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
         DOPAL(IND)=ENTOTL*ABSORP/DND
      ELSE

C***     DERIVATIVE WITH RESPECT TO UPPER LEVEL
         EMSPON=C3*RSTAR*EINST(NUP,LOW)/XLAMCM
C***     Set emissivities zero if both levels are equal (=POPMIN)
         IF (EN(LOW) .EQ. EN(NUP)) THEN
            DETAL(IND)=.0
            DOPAL(IND)=.0
         ELSE
            DETAL(IND)= ENTOTL*EMSPON/DND
            DOPAL(IND)=-ENTOTL*EMINDU/DND
         ENDIF
      ENDIF

    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
