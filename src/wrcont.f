C***  MAIN PROGRAM WRCONT  *****************************************************
      SUBROUTINE WRCONT
C*******************************************************************************
C***  THIS PROGRAM IS TO SOLVE THE CONTINUOUS RADIATION TRANSFER
C***  WITH GIVEN POPULATION NUMBERS
C*******************************************************************************
 
      IMPLICIT NONE
 
C***  DEFINE ARRAY DIMENSIONS ******************************************
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM =          26
      INTEGER, PARAMETER :: NDIM    =        1560
      INTEGER, PARAMETER :: NFDIM   = 2*NDIM + 400
      INTEGER, PARAMETER :: MAXIND  =       45000
      INTEGER, PARAMETER :: MAXFEIND  =       1500
      INTEGER, PARAMETER :: MAXKONT =     NFDIM/2
      INTEGER, PARAMETER :: MAXKODR =        NDIM
      INTEGER, PARAMETER :: NDDIM   =          89
      INTEGER, PARAMETER :: NPDIM   =          94
      INTEGER, PARAMETER :: MAXHIST =        4000 
      INTEGER, PARAMETER :: MAXXDAT =          10 
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 3200 
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      COMMON /COMAUTO/ LOWAUTO, WAUTO, EAUTO, AAUTO, IONAUTO, KRUDAUT
      CHARACTER*10 LEVUPAUTO(MAXAUTO), LEVAUTO(MAXAUTO)

      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(NDDIM) :: IWARN
      INTEGER, DIMENSION(NFDIM) :: KEY
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, EN
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      REAL, DIMENSION(NDDIM) :: OPA, ETA, THOMSON, RADIUS, ENTOT, T,
     >                          RNE, XJC
      REAL, DIMENSION(3,NDDIM) :: EDDI
      REAL, DIMENSION(4,NDIM) :: ALTESUM
      REAL, DIMENSION(4,MAXIND) :: COCO
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NPDIM) :: P, A, W, WX
      REAL, DIMENSION(NPDIM, NPDIM) :: BX

C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY ]
      REAL, DIMENSION(NPDIM, NPDIM) :: B
      REAL, DIMENSION(NPDIM) :: C

      REAL, DIMENSION(NFDIM) :: XLAMBDA, EMFLUX, EMCOLI, FWEIGHT, HEDDI
      REAL, DIMENSION(NDDIM, NFDIM) :: XJCOLD
      REAL, DIMENSION(NDDIM, NPDIM) :: U, Z
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C      INTEGER, PARAMETER :: INDEXMAX = 4000000, NFEREADMAX = 100000    !standard vd100
C      INTEGER, PARAMETER :: INDEXMAX = 5000000, NFEREADMAX = 150000    !vd50
C      INTEGER, PARAMETER :: INDEXMAX = 15000000, NFEREADMAX = 400000   !vd20
      INTEGER, PARAMETER :: INDEXMAX = 100000000, NFEREADMAX = 600000     !custom hydro
C      INTEGER, PARAMETER :: INDEXMAX = 12000000, NFEREADMAX = 300000   !Goetz

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      COMMON /IRON/ FEDUMMY, INDRB, INDRF, SIGMAFE,
     >              IFRBSTA, IFRBEND, IFENUP, IFELOW, SIGMAINT
      
      INTEGER :: L, N, IWARNES, LASTIND, LASTCHAR, IFROM, ITO, NATOM,
     >           NAUTO, LASTKON, LASTFE, JOBNUM, ND, NP, NF,
     >           LAST, NEXTK, LSOPA, LSINT, IFLUX, JOBMAX, IVERS, K,
     >           NCOLIP, LASTK, KDONE, IDUMMY, ISTATS, KWORDS, KBLOCKS,
     >           IERR, IWARNES2, N_WITH_DRLEVELS 
     
      REAL :: VDOPFE, XLAM0FE, TEFF, RSTAR, DTDR, TOTOUT, TOTIN, DXFE,
     >        BCORE, OPARND, TNOW, TUSED, TLEFT, TPERK, FLUXIN, DBDR,
     >        DUMMY, TBEGIN, POPMIN

      LOGICAL BFEMODEL

      LOGICAL :: NOTEMP, BUNLU, BKUDRITZKI, bCLXJC, bNoARAD, bNoXJC

      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXION) :: SIGMATHK, SEXPOK, EDGEK

      CHARACTER(MAXHIST*8) :: MODHIST
      CHARACTER(255) :: HISTENTRY
      CHARACTER(100) :: MODHEAD
      CHARACTER(48) :: BUFFER48
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8) :: NAME, BUFFER8
 
      INTEGER, EXTERNAL :: IDX
 
C***  Array indicating POPMIN levels (flagged by steal)
      LOGICAL, DIMENSION(NDIM,NDDIM) :: ZERO_RATES

C***  To Store Feautrier Matrices ST(94*95,89)
      REAL, DIMENSION((NPDIM+1)*NPDIM,NDDIM) :: ST
      LOGICAL :: BELIFI

C***  tiefenabh. clumping nach goetz
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !write to MODEL file      
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
C***  Write Link Data (Program Version) to CPR file
      WRITE (hCPR,'(2A)') '>>> WRCONT started: Program Version from '
     >                 ,LINK_DATE
      WRITE (hCPR,'(4A)') '>>> created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL

      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

C***  Initialize BELIFI; To Store Feautrier Matrices in Memory is now default
      BELIFI = .FALSE.

C***  Initialize XJC and other vectors
      XJC = 0.
      EDDI = 0.
      EMFLUX = 0.
      
C***  INITIALIZE COUNTER IWARNES, COUNTING SMALL EDDI'S (F) IN ELIMIN
      IWARNES = 0
      IWARNES2 = 0

      CALL       DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >                  'WRCONT', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL, 
     >             LEVUPAUTO, LEVAUTO, N_WITH_DRLEVELS, MAXION)
 
C***  WARNING: THIS DEFINITION IS DIFFERENT TO THAT DONE IN STEAL
C***  UNSOELD-LUCY-PROCEDURE IS ACTIVE ==> NOTEMP = FALSE
      IF (BUNLU) THEN
        NOTEMP = .FALSE.
      ENDIF

      BKUDRITZKI = .FALSE.

C***  READING OF THE MODEL FILE ****************************************
      CALL       RMODCON (ND,NDDIM,RADIUS,NP,NPDIM,P,Z,ENTOT,T,RNE,NF,
     $             NFDIM,MODHIST,MAXHIST,LAST,ALTESUM,XLAMBDA,
     >             TEFF,NOTEMP,XJCOLD,HEDDI,EDDI,NCOLIP,FWEIGHT,
     >             KEY,POPNUM,RSTAR,MODHEAD,JOBNUM,NEXTK,N,NDIM,
     $             MAXXDAT,XDATA, DENSCON, FILLFAC, OPARND, 
     >             POPMIN, ZERO_RATES, bNoARAD, bNoXJC)
      WRITE(hCPR,'(A,I7)') '>>> This is job number ', JOBNUM

C***  IN CASE OF FIRST WRCONT-JOB: SET OPTION "NOTEMP"
c      IF (JOBNUM .EQ. 2 .AND. BKUDRITZKI) THEN
      IF (JOBNUM .EQ. 2) THEN
        NOTEMP = .TRUE.
      ELSE
        NOTEMP = .FALSE.
      ENDIF
      
      IF (bNoXJC) THEN
C***    IF no XJC is found in the model, we cannot call DIFDTDR      
        NOTEMP = .TRUE.
      ENDIF
      
      CALL POPMIN_NULLING (ZERO_RATES, POPNUM, POPMIN, ND, N)

C***  DECODING INPUT OPTIONS *******************************************
      CALL       DECON (LSOPA,LSINT,IFLUX,JOBMAX,MODHIST, 
     >                  BUNLU, bCLXJC, IVERS, POPMIN)

      TOTIN=.0
      TOTOUT=.0
      IF (NEXTK .GT. 1) THEN
          CALL READMS (3,EMFLUX,NF,'EMFLUX  ', IERR)
          CALL READMS (3,TOTIN,1,  'TOTIN   ', IERR)
          CALL READMS (3,TOTOUT,1, 'TOTOUT  ', IERR)
      ENDIF
 
C***  MASS STORAGE ON FILE 7 FOR THE FEAUTRIER MATRICES
C***  ARRAY SIZE: SUM OVER 'LL' (SEE SUBR. ELIMIN) FROM L=1 TO ND-1
C***  NUMBER OF WORDS
      KWORDS=(ND-1)*((NP+2)*(6*NP+6)-ND*(6*NP+9)+ND*(2*ND-1))/6
C***  NUMBER OF BLOCKS (1 BLOCK = 512 WORDS) IN MEMORY
      KBLOCKS=(KWORDS+511)/512+1
C***  SET MEMORY SIZE (ISTATS=1: COLLECT STATISTICS)
      ISTATS=1
c      CALL WOPEN (7,KBLOCKS,ISTATS)
      IF (BELIFI) THEN
        CALL OPENMS (7, IDUMMY, IDUMMY, 0, IERR)
      ENDIF
 
C***  CALCULATION OF "DTDR" FOR USE IN SUBR. DIFFUS
C***  (ONLY IF DEFAULT OPTION "TEMPERATURE CORRECTION" IS SET)
      IF (.NOT. NOTEMP) CALL DIFDTDR (DTDR,TEFF,XJCOLD,HEDDI,T(ND),
     >                  RADIUS(ND),ND,EN,POPNUM,POPMIN,RNE(ND),
     >                  ENTOT(ND),RSTAR,NDIM,N,LEVEL,NCHARG,
     $                  WEIGHT,ELEVEL,EION,EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,NOM,NF,XLAMBDA,FWEIGHT,
     $                  MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     $                  KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC, 
     >                  BKUDRITZKI, OPARND)

      IF (OPSYS .EQ. 'CRAY') THEN
        CALL SECOND(TBEGIN)
      ENDIF

C***  SOLUTION OF THE TRANSFER EQUATION AT EACH FREQUENCY-POINT ********
      DO L=1, ND
         ENTOT(L) = ENTOT(L) * DENSCON(L)
      ENDDO
      DO 1 K=NEXTK,NF
      CALL COOP (XLAMBDA(K),ND,T,RNE,POPNUM,POPMIN,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $           NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $           ALPHA,SEXPO,
     $           ADDCON1, ADDCON2, ADDCON3, 
     $           IGAUNT,
     $           SIGMATHK,SEXPOK,EDGEK,0,DUMMY,DUMMY,RADIUS,
     $           KONTNUP,KONTLOW,LASTKON,XDATA)
      DO L=1, ND
         OPA(L) = OPA(L) * FILLFAC(L)
         ETA(L) = ETA(L) * FILLFAC(L)
      ENDDO
      CALL DIFFUS (XLAMBDA(K),T,RADIUS,ND,BCORE,DBDR,DTDR,TEFF,NOTEMP)
      IF (LSOPA.GT.0)
     $ CALL PRIOPA (XLAMBDA(K),K,ND,LSOPA,RADIUS,
     $ OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)
      CALL ELIMIN (XLAMBDA(K),EMFLUX(K),FLUXIN,U,Z,
     $          A,B,C,W,BX,WX,XJC,RADIUS,P,BCORE,DBDR,
     $          OPA,ETA,THOMSON,EDDI,ND,NP,NPDIM,ENTOT,K,
     $           IWARNES, IWARNES2, ST, BELIFI, IVERS)
C***  INTEGRATION OF THE TOTAL INCIDENT AND EMERGENT FLUX
      TOTIN=TOTIN+FLUXIN*FWEIGHT(K)
      TOTOUT=TOTOUT+EMFLUX(K)*FWEIGHT(K)
 
C***  WRITING ON THE MODEL FILE
      IF (bNoARAD .OR. (.NOT. bCLXJC)) THEN
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL WRITMS(3,XJC,ND,NAME,-1, IDUMMY, IERR)
      ENDIF
      IF (K <= 999) THEN
        WRITE (NAME, '(A4, I3, A1)') 'EDDI', K, ' '
      ELSE
        WRITE (NAME, '(A4, I4)') 'EDDI', K
      ENDIF
      CALL WRITMS (3,EDDI,3*ND,NAME,-1, IDUMMY, IERR)
 
C***  ABORT, IF NOT SUFFICIENT TIME LEFT FOR NEXT FREQUENCY POINT
      LASTK=K
      IF (K.EQ.NF) GOTO 1
      KDONE=K+1-NEXTK
      IF (OPSYS .EQ. 'CRAY') THEN
        CALL SECOND(TNOW)
        TUSED=TNOW-TBEGIN
        CALL TREMAIN (TLEFT)
        TPERK=TUSED/FLOAT(KDONE)
      ELSE
        TLEFT = 20.
        TPERK = 1.
      ENDIF
      IF (TLEFT.LT.TPERK*10.0) THEN
            CALL REMARK ('CP TIME NEARLY ELAPSED')
C!!!            ASSIGN 2 TO LABEL
C!!!            CALL REMARKF (LABEL,LASTK,NF)
            WRITE (0,2) LASTK,NF
    2       FORMAT (I3,' OF ',I3,' FREQUENCY POINTS COMPLETED')
            GOTO 6
            ENDIF
    1 CONTINUE
 
    6 CONTINUE

C***  OUTPUT OF WARNING FOR SMALL EDDIS (F)
      IF (IWARNES .GT. 0) THEN
        WRITE (hCPR,*) 'WRCONT > ELIMIN: ', IWARNES, ' WARNINGS:'
        WRITE (hCPR,*) 'Either:  Eddies (F) smaller than 0.01 '
        WRITE (hCPR,*) 'or    :  Jcont .le. 0' 
      ENDIF
      IF (IWARNES2 .GT. 0) THEN
        WRITE (hCPR,*) 'WRCONT > ELIMIN: ', IWARNES2, ' WARNINGS:'
        WRITE (hCPR,*) 'Eddies (F) larger than 1.0 '
      ENDIF

C***  UPDATING THE MODEL HISTORY
      !LASTCHAR = LAST * 8
      !IFROM = LASTCHAR + 1
      !ITO = LASTCHAR + 40
      IF (NOTEMP) THEN
        WRITE(UNIT=BUFFER48, FMT=8) 
     >                   JOBNUM,IVERS,LASTK, 'NOTEMP  '
      ELSE
        WRITE(UNIT=BUFFER48, FMT=8) 
     >                   JOBNUM,IVERS,LASTK, 'TEMP    '
      ENDIF
    8 FORMAT ('/',I7,'. WRCONT (V.',I2,')  LASTK=',I4,1X,A8)
      CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,48,BUFFER48)
   
      IF (LASTK == NF) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') 'COMPLETE'
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF
 
C***  UPDATING THE MODEL FILE
      CALL WRITMS (3,MODHIST,MAXHIST,    'MODHIST ',-1, IDUMMY, IERR)
C!!!      IF (.NOT. NOTEMP) CALL WRITMS 
C!!!     >  (3,TEFF,1,'TEFF    ',-1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,EMFLUX,NF,'EMFLUX  ',-1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,TOTIN,1,  'TOTIN   ',-1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,TOTOUT,1, 'TOTOUT  ',-1, IDUMMY, IERR)
      IF (LASTK .EQ. NF) THEN
            NEXTK=1
         ELSE
            NEXTK=LASTK+1
         ENDIF
      CALL WRITMS (hMODEL,NEXTK,1,  'NEXTK   ',-1, IDUMMY, IERR)
C***  Write NCOLIP
      IF (NCOLIP /= -1) THEN
        NCOLIP = 0
        CALL WRITMS(hMODEL, NCOLIP,  1, 'NCOLIP  ', -1, IDUMMY, IERR)
      ENDIF

 
C***  PRINTOUTS
      IF (LSINT.GT.0 .AND. LASTK .EQ. NF)
     $ CALL PRIINT (XJC,EDDI,RADIUS,ND,XLAMBDA,NF,LSINT,JOBNUM,MODHEAD)
      IF (IFLUX.GT.0 .AND. LASTK .EQ. NF) THEN
        !Define dummy array with COLI flux
        DO K=1, NF
          EMCOLI(K) = 0.
        ENDDO
        CALL PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,JOBNUM,
     >          FWEIGHT,MODHEAD,KEY,EMCOLI)
      ENDIF
 
C***  DEFINITION OF SUBSEQUENT JOB
      IF (LASTK .EQ. NF) THEN
        CALL REMARK ('WRCONT: NEXTJOB=REPEAT')
        PRINT *,'WRCONT: NEXTJOB=REPEAT'
        CALL JSYMSET ('G1','REPEAT')
      ELSE
        PRINT *,'WRCONT: ',LASTK,' OF ',NF,' FREQUENCIES COMPLETED'
        PRINT *,'WRCONT: NEXTJOB=WRCONT'
        CALL REMARK ('WRCONT: NEXTJOB=WRCONT')
        CALL JSYMSET ('G1','WRCONT')
      ENDIF

C***  MAX. NUMBER OF JOBS EXCEEDED?
      IF (JOBNUM .GE. JOBMAX) THEN
         CALL REMARK ('WRCONT: MAX. NUMBER OF JOBS EXCEEDED')
         PRINT *,'WRCONT: MAX. NUMBER OF JOBS EXCEEDED'
         CALL JSYMSET ('G3','ENDJOB')
         ELSE
         CALL JSYMSET ('G3','MOREJOBS')
         ENDIF 

      CALL CLOSMS (hMODEL, IERR)
      IF (BELIFI) THEN
        CALL CLOSMS (7,IERR)
      ENDIF

      !write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)

      CALL JSYMSET ('G0','0')

      CALL STAMP (OPSYS, 'WRCONT', TIM1)

      STOP 'O.K.'
      END