      SUBROUTINE RMODCON (ND,NDDIM,RADIUS,NP,NPDIM,P,Z,ENTOT,T,RNE,NF,
     $             NFDIM,MODHIST,MAXHIST,LAST,ALTESUM,XLAMBDA,
     $             TEFF,NOTEMP,XJCOLD,HEDDI,EDDI,NCOLIP,FWEIGHT,
     $             KEY,POPNUM,RSTAR,MODHEAD,JOBNUM,NEXTK,N,NDIM,
     $             MAXXDAT,XDATA, DENSCON, FILLFAC, OPARND, 
     >             POPMIN, ZERO_RATES, bNoARAD, bNoXJC)
C***********************************************************************
C***  READING OF THE MODEL FILE, CALLED FROM WRCONT ****************************
C***********************************************************************
      INTEGER, INTENT(IN) :: NDDIM, NPDIM, NDIM
      INTEGER, INTENT(INOUT) :: ND, NP, N

      DIMENSION XDATA(MAXXDAT)
      DIMENSION HEDDI(NF)
      DIMENSION XJCOLD(2),EDDI(2)
      REAL, DIMENSION(NDDIM-1) :: ARAD

      REAL, DIMENSION(NDDIM*NDIM) :: POPNUM
      LOGICAL, DIMENSION(NDIM*NDDIM) :: ZERO_RATES
      
      LOGICAL :: NOTEMP, bNoARAD, bNoXJC

      REAL, INTENT(INOUT) :: POPMIN      
      
c***  tiefenabh. clumping
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC
 
C *** NOTE: ASSURE 64-BIT TYPE FOR USE IN WRITMS
      CHARACTER(8) :: NAME
      
      bNoXJC = .FALSE.

      CALL OPENMS (3, IDUMMY, IDUMMY, 1, IERR)
C***  IN CASE OF FIRST WRCONT-JOB:  SET OPTION "NOTEMP"
      CALL READMS (3,JOBNUM,1,'JOBNUM  ', IERR)
      JOBNUM=JOBNUM+1
      CALL READMS(3,ND,1,'ND      ', IERR)
      IF(ND.GT.NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,RADIUS,ND,'R       ', IERR)
      CALL READMS (3,NP,1,'NP      ', IERR)
      IF (NP.GT.NPDIM) THEN
            CALL REMARK ('TOO MANY IMPACT PARAMETER POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,P,NP,         'P       ', IERR)
      CALL READMS (3,Z,ND*NP,      'Z       ', IERR)
      CALL READMS (3,ENTOT,ND,     'ENTOT   ', IERR)
      DENSCON(1:ND) = 0.0
      CALL READMS(3,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               STOP 'Error in depth-dep. clumping during STEAL'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3,XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK ('ERROR WHEN READING XDATA FROM MODEL FILE')
C***     XFILL EQ 0.
         XDATA(1) = 0.  
      ENDIF
      
C***  Determine if COLI has already been run (at least once)
      CALL READMS (3,ARAD,ND-1, 'ARAD    ', IERR)
      IF (IERR == -10) THEN
        bNoARAD = .TRUE.
      ELSE 
        bNoARAD = .FALSE.
      ENDIF
      
      CALL READMS (3,TEFF,1,'TEFF    ',IERR)
      CALL READMS (3,T,ND,       'T       ', IERR)
      CALL READMS (3,RNE,ND,     'RNE     ', IERR)
      CALL READMS(3,NF,1,        'NF      ', IERR)
      IF(NF.GT.NFDIM) THEN
            CALL REMARK ('TOO MANY FREQUENCY POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS(3,XLAMBDA,NF,  'XLAMBDA ', IERR)
      CALL READMS (3,FWEIGHT,NF, 'FWEIGHT ', IERR)
      CALL READMS (3,KEY,NF,     'KEY     ', IERR)
      CALL READMS (3,POPNUM,ND*N,'POPNUM  ', IERR)
      CALL READMS (3,POPMIN,1,   'POPMIN  ', IERR)
      IF (IERR < 0) THEN
C***  POPMIN not on MODEL file => read from CARDS      
        POPMIN = 1.E-100
      ENDIF
      CALL READMS (3,RSTAR,1,    'RSTAR   ', IERR)
      OPARND = 0.       
      CALL READMS (3, OPARND, 1, 'OPARND  ', IERR)

C***  Flags for the POPMIN levels
      CALL READMS (3,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
      IF (IERR == -10) THEN
        DO I=1, N*ND
          ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF 
 
      ND3=3*ND
      DO 15 K=1,NF
      WRITE (NAME,'(A3,I4,A1)') 'XJC',K, ' '
      CALL READMS (3,XJCOLD(1+ND*(K-1)),ND,NAME, IERR)
      IF (IERR == -10) THEN
C***    no XJC on model file
        bNoXJC = .TRUE.
        GOTO 15
      ENDIF
      IF (.NOT. NOTEMP) THEN
        IF (K <= 999) THEN
          WRITE (NAME,'(A4,I3,A1)') 'EDDI',K, ' '
        ELSE
          WRITE (NAME,'(A4,I4)') 'EDDI',K
        ENDIF
        CALL READMS (3,EDDI,ND3,NAME, IERR)
        HEDDI(K)=EDDI(ND3)
      ENDIF
   15 CONTINUE
 
      CALL READMS (3,MODHEAD,13,   'MODHEAD ', IERR)
      CALL READMS (3,NCOLIP,1,     'NCOLIP  ', IERR)
      CALL READMS (3,NEXTK,1,      'NEXTK   ', IERR)
      CALL WRITMS (3,JOBNUM,1,     'JOBNUM  ',-1, IDUMMY, IERR)
      CALL READMS (3,LAST,1,       'MODHIST ', IERR)
      IF (LAST.GT.MAXHIST) THEN
            CALL REMARK ('MODEL HISTORY DIMENSION TOO SMALL')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,MODHIST,LAST, 'MODHIST ', IERR)
 
      RETURN
      END
