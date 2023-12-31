      SUBROUTINE DECSTE (                                                   !Parameter count
     >     LSRAT,LSPOP,JOBMAX,EPSILON,REDUCE,IHIST,                         !6
     >     IFRRA, ITORA, IPRICC, IPRILC, LSEXPO, LSTAU, IFLUX,              !13
     >     IDAT, LEVELPL, MAXSETS, NPLOT, NDIM, NEWWRC,                     !19
     >     NGAMC,NGAMR,NGAML,NGAMD,AGAMC,AGAMR,AGAML,AGAMD,DELTAC,          !28
     >     TPLOT, TPLOTOPT, JOBNUM, NOEXTRAP, MODHIST,                      !33
     >     STHLP,NODATOM,NSCHAR,NOLAP,ITBR,ITMAX,OLDSTART,COMPO,            !41
     >     ELEMENT,SYMBOL,NATOM,NATOUT,BRRESET,                             !46
     >     DRLINES_CARD, NOPOP, BAUTO, BAUTO_ABORT, SMALLPOP, POPMIN,       !52
     >     BINBOX, POPRANG,                                                 !54
     >     COREX, BCOREX, VPLOT, BUNLU, UNLUTECLINE, BPRIUNLU,              !60
     >     ND, IPRINTZERORATES, PRILEVRA,                                   !63
     >     LEVELPLDEP, NPLOTDEP,                                            !65
     >     BITCONT, bTDIFFUS, BPGAHIST, AG, BAG, BPGAHISTE,                 !71
     >     PLOTOPT, MAXPLOTOPT, NPLOTOPT, OPC, BTALTER,                     !76
     >     BPLOCC, LPLOCC, KPLOCC, BRUDZERO,                                !80
     >     BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, IPLOT_XJLAPP, IPLOT_XJCAPP,    !85
     >     LPLOT_XJCAPP, NITER_PLOT_JAPP, BPLOTAPP, BNEWTONRESET,           !89
     >     XLAM_FINE_START, XLAM_FINE_END, BTRACE_POPNUM,                   !92
     >     BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX,                  !96
     >     BNUEFE, BXJLAPPCORE, WJCMIN, NKONV_THRESH,                       !100
     >     iBLOCKINVERSION, RSTAR, RMAX, TAUMAX, bENSURETAUMAX, TAUACC,     !106
     >     bTauStrict,                                                      !107
     >     RadiusGridParameters, VMIN, bThinWind, ThinCard, bSUMMARY,       !112   
     >     TaumaxCard, bHYDROSOLVE, bLateTau, HydroCard, AlphaCard,         !117
     >     DENSCON_FIX, DENSCON_LINE, MFORM, XMSTAR, GLOG, bUpdateMass,     !123
     >     WRTYPE, MODHEAD, bModHeadUpdate, bOLDVELO,                       !127
     >     bFULLHYDROSTAT, bGAMMARADMEAN, GEDDreduce, iZRType, iAMT,        !132
     >     CLUMP_SEP, MacroCard, OPALINE_SCALE, NBACKUP, GEFFKEY,           !137
     >     bHydroHelp, bHYSTloose, ICMMODE, bUCPP,                          !141
     >     bUpdateVT, VTURB, fTNDCOR, bTMTHOM, bKUBATDEBUG,                 !146
     >     bRELSCHARMERACC, bUseTWOPNT, FLUXEPS, HYSTACC, IHSSTATUS,        !151
     >     NDOUT, VTURB_LINE, ALOMIN, bLAMAPPCOLI, bDDVDOP,                 !156
     >     bFRACINV, bForceTDIFFUS, bNoFeTCORR, TRACEELEM, bINCADV,         !161
     >     ENTOT, RHO, iALOentry, bALOTri)                                  !165

C*******************************************************************************
C***  DECODING INPUT OPTIONS FOR MAIN PROGRAM "STEAL"
C*******************************************************************************

      IMPLICIT NONE

      INTEGER :: I, J, K, IN, MAXPLOTOPT, MAXSETS, NDIM, ND, L,
     >           NA, NATOM, NATOUT, ITBR, IPRINTZERORATES, NPAR,
     >           IND, ISTART, IERR, IHSSTATUS, NZ,
     >           LSTAU, LSRAT, LSPOP, LSEXPO, 
     >           IHIST, IFLUX, IDAT, IPRICC, IPRILC,
     >           JOBNUM, JOBMAX, NPLOT, NPLOTDEP, NPLOTOPT,
     >           NEWWRC, NSCHAR, ITMAX, MASSORIGIN, NKONV_THRESH,
     >           IDX, LPLOCC, KPLOCC, IFRRA, ITORA,
     >           MGC, MGR, MGL, MGD, MFORM, iAMT, ICMMODE, NDOUT,
     >           DELTAC, NITER_PLOT_JAPP, iBLOCKINVERSION, iZRType,
     >           IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NBACKUP

      INTEGER, INTENT(IN) :: iALOentry

      !INTEGER :: IFIX   !IFIX(x) ist ein Alias von INT(x) : konvertiert Real in Integer


      CHARACTER(256) :: KARTE, TPLOTOPT
      CHARACTER(1200) :: UNLUTECLINE, HydroCard
      CHARACTER(120) :: DENSCON_LINE, VTURB_LINE, DRLINES_CARD,
     >                 ThinCard, AlphaCard, MacroCard, TaumaxCard
      CHARACTER(256), DIMENSION(MAXPLOTOPT) :: PLOTOPT
      CHARACTER(10), DIMENSION(MAXSETS, NDIM) :: LEVELPL, LEVELPLDEP
      CHARACTER(10), DIMENSION(NATOM) :: ELEMENT
      LOGICAL, DIMENSION(NATOM) :: TRACEELEM
      CHARACTER(2), DIMENSION(NATOM) :: SYMBOL
      CHARACTER(10) :: PRILEVRA
      CHARACTER(8) :: OPC, GEFFKEY
      CHARACTER(2) :: WRTYPE
      CHARACTER(40) :: CURPAR
      CHARACTER(100) :: MODHEAD, HEADLINE

      REAL, DIMENSION(10) :: AGAMC,AGAMR,AGAML,AGAMD
      REAL, DIMENSION(7) :: AG
      REAL, DIMENSION(ND) :: VTURB, ENTOT, RHO

      INTEGER, DIMENSION(1) :: MODHIST
      INTEGER, DIMENSION(10) :: NGAMC,NGAMR,NGAML,NGAMD


      LOGICAL :: NOEXTRAP, TPLTAU, STHLP, NODATOM, NOLAP, OLDSTART,
     >           COMPO, DRNORUD, NOPOP, BAUTO, BAUTO_ABORT, NEWATOM,
     >           BCOREX, BUNLU, BPRIUNLU, bUCPP, bTMTHOM, BNEWTONRESET,
     >           bUseTWOPNT, BRUDZERO, RMAX_IN_RSUN, bRELSCHARMERACC,
     >           BITCONT, bTDIFFUS, BPGAHIST, BAG, BPGAHISTE, bDDVDOP,
     >           bForceTDIFFUS, bNoFeTCORR, bALOTri
      REAL :: RSTAR, RMAX, DENSCON_FIX, XL, XLC, XMSTAR, HYSTACC,
     >        EPSILON, REDUCE, Y0, BRRESET, COREX, ALOMIN,
     >        YMAX, YMIN, XMAX, XMIN, fTNDCOR,
     >        GAMMA, GAMMAC, GAMMAR, GAMMAL, GAMMAD, 
     >        ANG, ANGC, ANGR, ANGL, ANGD,
     >        SMALLPOP, SMALLPOP2, POPMIN, POPRANG, WJCMIN,
     >        XLAM_FINE_START, XLAM_FINE_END, CLUMP_SEP, OPALINE_SCALE,
     >        tempREAL, GEDDreduce, HYSTACCMIN, RONSET, RSTARorg
      REAL :: VFINAL, VMIN, VMINCAND, BETA, BETA2, BETA2FRACTION,
     >        VPAR1, VPAR2, VPAR1_2, VPAR2_2, RCON, HSCALE, VTURBND,
     >        TAUMAX, TAUACC, GLOG, POPMIN_OLD, VCON, FLUXEPS
      LOGICAL NOTEMP,TPLOT, BINBOX, VPLOT, bFRACINV, bINCADV
      LOGICAL BPLOCC, BNUEFE
      LOGICAL BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE,
     >        BPLOTAPP, BTRACE_POPNUM, BTALTER,
     >        BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX,
     >        bENSURETAUMAX, bHYDROSOLVE, bLateTau,
     >        bThinWind, bTauStrict, bSUMMARY, bUpdateMass, 
     >        bModHeadUpdate, bOLDVELO, bLAMAPPCOLI,
     >        bFULLHYDROSTAT, bGAMMARADMEAN, bKUBATDEBUG,
     >        bHydroHelp, bHYSTloose, bUpdateVT, bSMOCO
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters           !contains all RADIUS-GRID CARDS (for subr. RGRID)

      COMMON/VELPAR/ VFINAL,VMINCAND,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >            BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2


C***  Variablen zum Decodieren der Eingabezeilen
C***  Nach Einlesen einer Zeile und Parsen derselben
C***  enthaelt NPAR die Anzahl der Parameter in dieser Zeile
C***  ACTPAR(1)...ACTPAR(NPAR) enthalten die Parameter
      CHARACTER(40), DIMENSION(20) :: ACTPAR

C***  SOLAR RADIUS ( CM )
      REAL, PARAMETER :: RSUN = 6.96E10
      REAL, PARAMETER :: GCONST = 6.670E-8     !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33      !XMSUN = Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)



C***  DEFAULT SETTINGS ---------------------------------------------
C***  NO PRINT OPTIONS SET 
      LSTAU=-1
      LSRAT=-1
      LSPOP=-1
      LSEXPO=-1
      IHIST=-1
      IFLUX=-1
      IDAT=-1
      IPRICC=0
      IPRILC=0
      BPRIUNLU  = .FALSE.
      BPGAHIST  = .FALSE.
      BPGAHISTE = .FALSE.
      BAG = .FALSE.
      BGAMMARFIX = .FALSE. 
      BGAMMALFIX = .FALSE. 
      BGAMMACFIX = .FALSE. 
      BGAMMADFIX = .FALSE. 
C***  CORRECT LTE-SOURCEFUNCTION OF IRON-TRANSITIONS TO PLANCK VALUE 
      BNUEFE = .TRUE.
C***  JOBMAX EXCEEDED 
      JOBMAX=-1
C***  DEMANDED CONVERGENCE ACCURACY 
      EPSILON=.005
C***  NO REDUCTION OF CORRECTIONS
      REDUCE = 1.
C***  NO POPNUMBER PLOT 
      NPLOT=0
C***  NO DEPARTURE PLOT 
      NPLOTDEP = 0
C***  NO TEMPERATURE PLOT (DEFAULT)
      TPLOT=.FALSE.
      TPLTAU=.FALSE.
C***  NO VELOCITY PLOT
      VPLOT = .FALSE.
C***  NO PLOTOPTIONS STORED
      NPLOTOPT = 0
C***  NUMBER OF REPEAT ITERATIONS BETWEEN TWO WRCONT-JOBS 
C***   This default must agree with the definition in DECCOLI (COLI)
      NEWWRC=6
C***  NO SCHARMER AMPLIFICATION OF CORRECTIONS
      DO 21 I=1,10
      NGAMC(I)=9999999999
      NGAMR(I)=9999999999
      NGAML(I)=9999999999
      NGAMD(I)=9999999999
      AGAMC(I)=0.
      AGAMR(I)=0.
      AGAML(I)=0.
      AGAMD(I)=0.
   21 CONTINUE
      NGAMC(1)=0
      NGAMR(1)=0
      NGAML(1)=0
      NGAMD(1)=0
      MGC=1
      MGR=1
      MGL=1
      MGD=1
      DELTAC=1.0
      bLAMAPPCOLI = .FALSE.
      bALOTri = .FALSE.
C***  DR-RATES CALCULATED WITH CONTINUUM RADIATION FIELD
      DRLINES_CARD = ''
C***  RUDIMENTAL LINES ARE SET TO f=0 (NO RADIATIVE RATES)
      BRUDZERO = .FALSE.
C***  DISTANCE MODULUS 
      Y0=.0
C***  NO SUPPRESSION OF TEMPERATURE CORRECTIONS 
      NOTEMP=.FALSE.
C***  UNSOELD-LUCY PROCEDURE IS NOT STANDARD
      BUNLU = .FALSE.
C***  NO SUPPRESSION OF EXTRAPOLATION JOBS
C***    BCOREX IS USED TO SUPPRESS EXTRAPS UNLESS THE CORRECTION IS HIGHER
C***    THAN A TRESHOLD
      NOEXTRAP = .FALSE.
      BCOREX = .FALSE.
C***  NORMAL STEAL FUNCTION (NOT ONLY OUTPUT REQUEST)
      STHLP=.FALSE.
C***  Iteration of the Continuum at low JOBNUMs
      BITCONT = .TRUE.
C***  SUPPRESS OUTPUT OF ATOMIC DATA IN CASE OF CONVERGED MODEL
      NODATOM=.TRUE.
C***  SUPPRESS OUTPUT OF POPNUMBERS IN CASE OF CONVERGED MODEL
      NOPOP=.TRUE.
C***  OVERLAP OPTION NOT DISABLED
      NOLAP = .FALSE.
C***  BROYDEN BRANCH NOT ACTIVATED
      ITBR = 0
C***  CODE-NUMBER OF ELEMENT WHICH SHALL BE OMMITTED FROM SCHARMER TREATMENT
C***     DEFAULT: -1 MEANS NO VALID CODE NUMBER
      NSCHAR=-1
C***  MAXIMUM NUMBER OF NEWTON-RAPHSON ITERATIONS
      ITMAX = 50
C***  No blockwise inversion used for M^-1
      iBLOCKINVERSION = 0.
C***  STARTJOB (JOBNUM = 1): NO OLDSTART OPTION
      OLDSTART = .FALSE.
      NEWATOM = .FALSE.
C***  OUTPUT (PRINT DATOM, POP, RATE) ONLY FOR SPECIFIED ELEMENT
      NATOUT=0
C***  BRRESET IF NORM(Fnew) > NORM(Fold)
      BRRESET = 1.
      BNEWTONRESET = .TRUE.
C***  NO OUTPUT OF CHEMICAL COMPOSITION
      COMPO=.FALSE.
C***  AUTO MODIFY DEFAULT
      BAUTO = .TRUE.
C***  NO Job-Abort if Number of diverged Points is greater treshhold
      BAUTO_ABORT = .FALSE. 
C***  NO INBOX
      BINBOX = .FALSE.
C***  DIFFUSION APPROXIMATION is now standard (old TDIFFUS cards line)
      bTDIFFUS = .TRUE.
      bForceTDIFFUS = .FALSE.
      fTNDCOR = 0.
C***  Default: TAUMAX value does not refer to THOMSON scale
      bTMTHOM = .FALSE.
C***  Default: No debug output for thermal balance routine
      bKUBATDEBUG = .FALSE.
C***  Default: Absolute accuracy for Scharmer iteration (instead of relative)
      bRELSCHARMERACC = .FALSE.
C***  Default: Allow the use of the Two-point Broyden method
      bUseTWOPNT = .TRUE.
C***  Default: Do not invert DM matrix for fractional corrections
      bFRACINV = .FALSE.
C***  Default: Do not switch off temperature correction of COLI Fe rates
      bNoFeTCORR = .FALSE.
C***  Default: Do not consider advection terms in INITFCORR
      bINCADV = .FALSE.
      
C***  SMALLPOP: smaller POPNUMS are not accounted for covergence criterion
C***            of the inner iteration
      SMALLPOP = 1.E-12
C***  NDOUT: outermost depth point considered for convergence criterion 
C***         Default: Negative or zero value = automatic calculation in PRICORR
      NDOUT = -1

C***  POPMIN: smaller POPNUMs are set to this value 
      IF (POPMIN < 1.E-99) THEN
        POPMIN = 1.E-25
      ELSE 
        POPMIN_OLD = POPMIN
      ENDIF

C***  RANGE OF POPNUM PLOT
      POPRANG = -15.
C***  Initialize PLOT POP and DEPART LEVELS
      DO I=1, MAXSETS
         DO J=1, NDIM
            LEVELPL   (I,J) = '          '
            LEVELPLDEP(I,J) = '          '
         ENDDO
      ENDDO

C***  Radius and velocity defaults
      RMAX=-1.0
      RMAX_IN_RSUN = .FALSE.
      BETA2FRACTION = .0      

C***  Initialize radius grid lines
      DO i=1, 3, 1
        RadiusGridParameters(i) = ' '
      ENDDO

C***  Clumping stuff 
      DENSCON_FIX = 1.
      DENSCON_LINE = ' '
C***  Depth-dependent microturbulence
      VTURB_LINE = ' '
      
C***  Do not use new Integration of Approximate Radiation Fields
      BXJLAPPNEW  = .FALSE.
      BXJLAPPCORE = .FALSE.
      BXJCAPPNEW  = .FALSE.
C***  Do not use new Operator
      BNEWOPER   = .FALSE.
C***  Default values for extension of the fine grid
      XLAM_FINE_START =  100.
      XLAM_FINE_END   = 2000.
C***  No Plots of Approximate Radiation Fields
      IPLOT_XJLAPP = 1
      IPLOT_XJCAPP = 1
      LPLOT_XJCAPP = 1
      NITER_PLOT_JAPP = 2
      BPLOTAPP = .FALSE.

C***  Operator
      OPC = 'NONE'

C***  Alternate temperature corrections - perform them only every second time (if last COLI /= COLI+)
      BTALTER = .FALSE.

C***  Plot of Continuum Operator
      BPLOCC = .FALSE.

C***  Defaults for auxilary AUTO GAMMA Parameters
      AG(5) = 1.
      AG(6) = 1.
      AG(7) = 1.

C***  No Tracing of Popnumbers
      BTRACE_POPNUM = .FALSE.

C***  Minimum Scharmer Weight
      WJCMIN = 0.1    !continua
      ALOMIN = 0.1    !lines

C***  Threshold for maximum number of allowed diverged points in LINPOP
      NKONV_THRESH = 30      
      
c***  No printout of ZERO_RATES (levels which are omitted) in CPR file
      IPRINTZERORATES = 0

C***  new steal options from wrstart (added 08.08.2011 by ansander)
      bThinWind = .FALSE.
      bFULLHYDROSTAT = .FALSE.      !if true, the eddington gamma is build from a_rad instead of a_thom
      bGAMMARADMEAN = .FALSE.       !if true, uses mean value for GAMMA_rad instead of individual per depth point
      TAUMAX = .0
      TAUACC = 1.E-4
      bENSURETAUMAX = .FALSE.
      bTauStrict = .TRUE.
      GEDDreduce = 1.               !default: do not take a fraction of last GEDDRAD
      FLUXEPS = -1.
      bUCPP = .FALSE.
      HYSTACC = -0.05               !Accuracy for hydrostatic equation (negative = not forced for convergence)
      HYSTACCMIN = 0.001            !minimum input value for HYSTACC

C***  HYDRO PARAMATERS (added 17.08.2011 by ansander)           
      bHYDROSOLVE = .FALSE.     !Default: hydro branch is not used
      bLateTau = .FALSE.        !Default: Tau iteration is done directly after hydro
      bHYSTloose = .FALSE.
      
      UNLUTECLINE = 'UNLUTEC'
      HydroCard = 'HYDRO'
      
C***  Use Depth-dependent doppler profiles? (Default = no)
      bDDVDOP = .FALSE.
      
C***  Write model parameter summary file in this steal job
      bSUMMARY = .FALSE.

      NBACKUP = 0               !Default: No model backups during STEAL job
      iZRType = 1               !Default: Standard method for ZERO_RATES determination
                                !Values: 0 = no ZERO_RATES, 1 = standard/wrh, 2 = Goetz method
      iAMT = 0                  !AUTO MODIFY TEMPERATURE switch (0 = UNLU, 1 = TEMP, -1 = NOTEMP)                                
      ICMMODE = 1               !CORRMAX calculation modus (0 = consider everything, 1 = consider only converged)
      
      !how to calculate stellar mass (needed? perhaps not here but in wrstart)
      ! 0 = Mass-luminosity relation   
      ! 1 = direct CARDS input   
      ! 2 = log g
      MASSORIGIN = 0
      MFORM = 2         !New default: Use Goetz formula if MASSORIGIN = 0
      bUpdateMass = .FALSE.
      bModHeadUpdate = .FALSE.
      bUpdateVT = .FALSE.
      bOLDVELO = .FALSE.
      CLUMP_SEP = 0.            !default: no macroclumping
      OPALINE_SCALE = 1.        !default: no scaling of line opacities (changes a_lines)
      DO K=1, NATOM
        TRACEELEM(K) = .FALSE.
      ENDDO

C***  END OF DEFAULT SETTINGS  -----------------------------------------

      RSTARorg = RSTAR          !Backup input RSTAR (for update option)

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
    8 READ (1,4, END=99) KARTE
    4 FORMAT (A)

C***  COUNT Number of Parameters in Line KARTE 
C***  Skip Line, if line has zero parameters (i. e. Line is empty)
C***  Ignore more than 10 Parameter
      CALL SARGC(KARTE,NPAR)
      IF ( NPAR .LT. 1) GOTO 8
      IF ( NPAR .GT. 20) NPAR = 20

C***  Get actual parameters
      DO 20, I=1, NPAR
       CALL SARGV(KARTE,I,ACTPAR(I))
   20  CONTINUE
C***  Erase not used parameter fields
      DO 98, I=NPAR+1, 20
        ACTPAR(I) = '0.0'
98    CONTINUE


C***  NO BROYDEN IMMIDEATELY AFTER WRSTART (I.E. JOBNUM = 1)
      IF (ACTPAR(1) .EQ. 'BROYDEN' .AND. ACTPAR(2) .NE. 'RESET' .AND. 
C                         =======                        =====
     >   ACTPAR(2) .NE. 'TWOPOINT' .AND. JOBNUM .GT. 1) THEN
            READ (ACTPAR(2),18) XL
   18       FORMAT (G15.0)
            IF (XL .EQ. .0) THEN
              ITBR = 2
              ELSE
              ITBR = IFIX(XL)
              ENDIF
            GOTO 8
            ENDIF


      IF (ACTPAR(1) .EQ. 'BROYDEN' .AND. ACTPAR(2) .EQ. 'RESET' .AND. 
     >    JOBNUM .GT. 1) THEN
C                         =======                        =====      
            IF (ACTPAR(3) .EQ. 'INFINITY') THEN
              BNEWTONRESET = .FALSE.
            ELSE
              READ (ACTPAR(3),19) BRRESET
   19         FORMAT (F8.4)
              IF (BRRESET .LE. .0) THEN
                BRRESET = 1.
                ENDIF
              GOTO 8
            ENDIF
            ENDIF

      IF (ACTPAR(1) == 'BROYDEN' .AND. ACTPAR(2) == 'TWOPOINT') THEN
        IF (ACTPAR(3) == 'ON') THEN
          bUseTWOPNT = .TRUE.
        ELSEIF (ACTPAR(3) == 'OFF') THEN
          bUseTWOPNT = .FALSE.
        ENDIF
      ENDIF

      IF (ACTPAR(1) .EQ. 'OUTPUT') THEN
C                         ======
            STHLP=.TRUE.
            WRITE(0,*) 'DECSTE: OUTPUT ONLY MODE forced'

      ELSE IF (ACTPAR(1) .EQ. 'NEWWRC' ) THEN
C                              ======
            READ (ACTPAR(2),'(F10.0)',ERR=1) XL
            NEWWRC=IFIX(XL)

C***  PLOT OPTIONS *******************************
      ELSE IF (ACTPAR(1) .EQ. 'PLOT') THEN 
         IF (ACTPAR(2) .EQ. 'POP') THEN
C                            ===
            NPLOT=NPLOT+1
            IF (NPLOT .GT. NDIM ) GOTO 8
            DO 2 I=1,MAXSETS
            IF (NPAR .GE. I+2) LEVELPL(I,NPLOT) = ACTPAR(I+2)
    2       CONTINUE

         ELSE IF (ACTPAR(2) .EQ. 'DEPART') THEN
C                                 ======
            NPLOTDEP = NPLOTDEP + 1
            IF (NPLOTDEP .GT. NDIM ) GOTO 8
            DO I=1,MAXSETS
              IF (NPAR .GE. I+2) LEVELPLDEP(I,NPLOTDEP) = ACTPAR(I+2)
            ENDDO

         ELSE IF (ACTPAR(2)(1:1) .EQ. 'T') THEN
C                                      =
            !DECODE (80,15,KARTE) YMAX,YMIN,XMAX,XMIN
            TPLOTOPT = KARTE
            TPLOT=.TRUE.
C            IF (DKARTE(12:15) .EQ. 'TAUR') TPLTAU=.TRUE.

         ELSE IF (ACTPAR(2)(1:1) .EQ. 'V') THEN
C                                      =
            VPLOT=.TRUE.

         ELSE IF (ACTPAR(2) .EQ. 'JNUE'  .OR. 
     >            ACTPAR(2) .EQ. 'JLINE' .OR. 
     >            ACTPAR(2) .EQ. 'UNLU'  .OR. 
     >            ACTPAR(2) .EQ. 'HTOT'  .OR. 
     >            ACTPAR(2) .EQ. 'FLUX'  .OR. 
     >            ACTPAR(2)  ==  'ALPHA' .OR. 
     >            ACTPAR(2)  ==  'GAMMA' .OR. 
     >            ACTPAR(2)  ==  'SIGMAFE' .OR. 
     >            ACTPAR(2)  ==  'FGSTRAT' .OR. 
     >            ACTPAR(2)  ==  'FORCEMULT' .OR. 
     >            ACTPAR(2)  ==  'ACCELEM' .OR. 
     >            ACTPAR(2) .EQ. 'ACC'        ) THEN
            IF (NPLOTOPT .EQ. MAXPLOTOPT) THEN
               WRITE (0, '(A,/,A,/,A)') 'STEAL: WARNING ----  MORE ' //
     >           'PLOT OPTIONS REQUESTED THAN DIMENSIONED (MAXPLOTOPT)',
     >           'THE FOLLOWING OPTION IS IGNORED: ', KARTE(:IDX(KARTE))
            ELSE
               NPLOTOPT = NPLOTOPT + 1
               PLOTOPT(NPLOTOPT) = KARTE
            ENDIF
         ELSE IF (ACTPAR(2) .EQ. 'CCORE') THEN
            BPLOCC  = .TRUE.
            READ (ACTPAR(3),'(I8)') LPLOCC
            READ (ACTPAR(4),'(I8)') KPLOCC
            write (0,*) 'DECSTE: ', bplocc, lplocc, kplocc
         ELSE IF (ACTPAR(2) .EQ. 'XJLAPP') THEN
            BPLOTAPP = .TRUE.
            READ (ACTPAR(3),'(I8)') IPLOT_XJLAPP
            IF (NPAR .GT. 3) THEN
              IF (ACTPAR(4) .EQ. 'ITERATION') THEN
                READ (ACTPAR(5),'(I8)') NITER_PLOT_JAPP
              ENDIF
            ENDIF
         ELSE IF (ACTPAR(2) .EQ. 'XJCAPP') THEN
            BPLOTAPP = .TRUE.
            READ (ACTPAR(3),'(I8)') IPLOT_XJCAPP
            READ (ACTPAR(4),'(I8)') LPLOT_XJCAPP
            IF (NPAR .GT. 4) THEN
              IF (ACTPAR(5) .EQ. 'ITERATION') THEN
                READ (ACTPAR(6),'(I8)') NITER_PLOT_JAPP
              ENDIF
            ENDIF

         ENDIF

C********* PRINT OPTIONS *********************************
      ELSE IF (ACTPAR(1) .EQ. 'PRINT') THEN

         IF (ACTPAR(2) .EQ. 'RATES' ) THEN
C                            =====  
            LSRAT = 1
            IF (NPAR .GT. 2) THEN
               READ (ACTPAR(3),'(I10)', ERR=1) LSRAT
            ENDIF
            IFRRA = 0
            ITORA = 0
            PRILEVRA = ' '
            DO I = 4, NPAR-1
              IF (ACTPAR(I) .EQ. 'FROM')
     >            READ(ACTPAR(I+1),'(I10)',ERR=1) IFRRA
              IF (ACTPAR(I) .EQ. 'TO')
     >            READ(ACTPAR(I+1),'(I10)',ERR=1) ITORA
              IF (ACTPAR(I) .EQ. 'LEVEL')
     >            PRILEVRA = ACTPAR(I+1)
            ENDDO

         ELSE IF (ACTPAR(2) .EQ. 'HIST' ) THEN
C                                 ====
            IHIST=1
         ELSE IF (ACTPAR(2) .EQ. 'GAMMA' .AND. 
     >            ACTPAR(3) .EQ. 'HIST') THEN
C                                 =====
            BPGAHIST = .TRUE.
            IF (ACTPAR(4) .EQ. 'EXTEND') THEN
              BPGAHISTE = .TRUE.
            ENDIF
         ELSE IF (ACTPAR(2) .EQ. 'FLUX' ) THEN
C                                =====
            IFLUX=1
         ELSE IF (ACTPAR(2) .EQ. 'DATOM' ) THEN
C                                 ===== 
            IF (ACTPAR(3) == 'IFCONVERGED') THEN
              NODATOM = .FALSE.
            ELSE
              IDAT=1
            ENDIF
         ELSE IF (ACTPAR(2) .EQ. 'TAU' ) THEN
C                                 ===
            READ(ACTPAR(3),'(F10.0)',ERR=1) XL
            LSTAU=IFIX(XL)
            IF (LSTAU.EQ.0) LSTAU=1
         ELSE IF (ACTPAR(2) .EQ. 'POP' ) THEN
C                                 ===
            IF (ACTPAR(3) .EQ. 'IFCONVERGED') THEN
               NOPOP = .FALSE.
            ELSE
               READ (ACTPAR(3),'(F10.0)',ERR=1) XL
               LSPOP=IFIX(XL)
               IF (LSPOP.EQ.0) LSPOP=1
            ENDIF
         ELSE IF (ACTPAR(2) .EQ. 'EXPOP' ) THEN
C                                 =====
            READ(ACTPAR(3),'(F10.0)',ERR=1) XL
            LSEXPO=IFIX(XL)
            IF (LSEXPO .EQ. 0) LSEXPO=1
         ELSE IF (ACTPAR(2) .EQ. 'CCORE') THEN
C                                 =====
            IPRICC=1
         ELSE IF (ACTPAR(2) .EQ. 'LCORE') THEN
C                                 =====
            READ (ACTPAR(3),'(F10.0)',ERR=1) XLC
            IPRILC=IFIX(XLC)
         ELSE IF (ACTPAR(2) .EQ. 'ELEMENT') THEN
C                                 =======
            DO 11 NA=1,NATOM
            IF ((ACTPAR(3) .EQ. ELEMENT(NA)) .OR. 
     >          (ACTPAR(3) .EQ. SYMBOL(NA))) NATOUT=NA
   11       CONTINUE

         ELSE IF (ACTPAR(2)(1:5).EQ.'COMPO' ) THEN
C                                    =====   
            COMPO=.TRUE.

         ELSE IF (ACTPAR(2) .EQ. 'UNLU') THEN
            BPRIUNLU = .TRUE.

         ELSE IF ((ACTPAR(2) .EQ. 'SUMMARY') .OR. 
     >            (ACTPAR(2) .EQ. 'MODINFO')) THEN
            bSUMMARY = .TRUE.

         ELSE IF (ACTPAR(2) .EQ. 'ZERORATES') THEN
            IPRINTZERORATES = 1

         ENDIF

C***  END OF PRINT OPTIONS ****************************************

      ELSE IF (ACTPAR(1) .EQ. 'JOBMAX' ) THEN
C                              ======
            READ (ACTPAR(2),'(F10.0)',ERR=1) XL
            JOBMAX=IFIX(XL)
C            IF (JOBMAX >= 1000) THEN
              !PRINT *,'JOBMAX .GT. 1000: INVALID OPTION'
C              CALL REMARK ('WARNING: JOBMAX >= 1000 not recommended')
              !CALL REMARK ('JOBMAX .GT. 1000: INVALID OPTION')
             !STOP 'ERROR'
C            ENDIF

      ELSE IF (ACTPAR(1) .EQ. 'EPSILON' ) THEN
C                              =======
            READ(ACTPAR(2),'(F10.0)',ERR=1) EPSILON

      ELSE IF (ACTPAR(1) .EQ. 'FLUXEPS' ) THEN
C                              =======
            READ(ACTPAR(2),'(F10.0)',ERR=1) FLUXEPS

      ELSE IF (ACTPAR(1) .EQ. 'REDUCE' ) THEN
C                              ======
            READ (ACTPAR(2),'(F10.0)',ERR=1) REDUCE
            IF (REDUCE.EQ..0) REDUCE=.5

      ELSE IF (ACTPAR(1) .EQ. 'AUTO' .AND. ACTPAR(2) .EQ. 'GAMMA') THEN
C                              ====                        =====
            BAG = .TRUE.
            READ (ACTPAR(3), '(F10.0)') AG(1)
            READ (ACTPAR(4), '(F10.0)') AG(2)
            READ (ACTPAR(5), '(F10.0)') AG(3)
            READ (ACTPAR(6), '(F10.0)') AG(4)
            IF (NPAR .GE. 7) READ (ACTPAR(7), '(F10.0)') AG(5)
            IF (NPAR .GE. 8) READ (ACTPAR(8), '(F10.0)') AG(6)
            IF (NPAR .GE. 9) READ (ACTPAR(9), '(F10.0)') AG(7)
      ELSE IF (ACTPAR(1) .EQ. 'GAMMA') THEN
C                              =====
            MGC=MGC+1
            MGR=MGR+1
            MGL=MGL+1
            MGD=MGD+1
            IF (MGC.GT.10 .OR. MGR.GT.10 .OR. MGL.GT.10 .OR. MGD.GT.10)
     >        THEN
               CALL REMARK ('TOO MANY GAMMA OPTIONS')
               STOP 'ERROR'
               ENDIF
            READ (ACTPAR(2), '(F10.0)') GAMMA
            READ (ACTPAR(5), '(F10.0)') ANG
            NGAMC(MGC)=ANG
            NGAMR(MGR)=ANG
            NGAML(MGL)=ANG
            NGAMD(MGL)=ANG
            AGAMC(MGC)=GAMMA
            AGAMR(MGR)=GAMMA
            AGAML(MGL)=GAMMA
            AGAMD(MGL)=GAMMA
      ELSE IF (KARTE(:6) .EQ. 'GAMMAC' ) THEN
C                              ======
            BGAMMACFIX = .TRUE.
            MGC=MGC+1
            IF (MGC .GT. 10) GOTO 1
            READ (ACTPAR(2), '(F10.0)') GAMMAC
            READ (ACTPAR(5), '(F10.0)') ANGC
            NGAMC(MGC)=ANGC
            AGAMC(MGC)=GAMMAC
      ELSE IF (KARTE(:6) .EQ. 'GAMMAR' ) THEN
C                              ======
            BGAMMARFIX = .TRUE.
            MGR=MGR+1
            IF (MGR .GT. 10) GOTO 1
            READ (ACTPAR(2), '(F10.0)') GAMMAR
            READ (ACTPAR(5), '(F10.0)') ANGR
            NGAMR(MGR)=ANGR
            AGAMR(MGR)=GAMMAR
      ELSE IF (KARTE(:6) .EQ. 'GAMMAL' ) THEN
C                              ======
            BGAMMALFIX = .TRUE.
            MGL=MGL+1
            IF (MGL .GT. 10) GOTO 1
            READ (ACTPAR(2), '(F10.0)') GAMMAL
            READ (ACTPAR(5), '(F10.0)') ANGL
            NGAML(MGL)=ANGL
            AGAML(MGL)=GAMMAL
      ELSE IF (KARTE(:6) .EQ. 'GAMMAD' ) THEN
C                              ======
            BGAMMADFIX = .TRUE.
            MGD=MGD+1
            IF (MGD .GT. 10) GOTO 1
            READ (ACTPAR(2), '(F10.0)') GAMMAD
            READ (ACTPAR(5), '(F10.0)') ANGD
            NGAMD(MGL)=ANGD
            AGAMD(MGL)=GAMMAD

      ELSE IF (KARTE(:6) .EQ. 'DELTAC') THEN
C                              ======
            READ (ACTPAR(2), '(F10.0)', ERR=1) DELTAC

      ELSE IF (KARTE(:6) .EQ. 'DRLINE') THEN
C                              ======
            DRLINES_CARD = KARTE

      ELSE IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                              ========
            NOTEMP=.TRUE.
            IF (KARTE(42:46) .NE. ' ') THEN
              CALL DECNOT (MODHIST(1),MODHIST,KARTE,NOTEMP,'STEAL')
            ENDIF
      ELSE IF (KARTE(:7) .EQ. 'UNLUTEC') THEN
C                              =======
            UNLUTECLINE = TRIM(UNLUTECLINE) // KARTE(8:)
            BUNLU = .TRUE.

      ELSE IF (ACTPAR(1).EQ.'NO'.AND.ACTPAR(2)(1:6).EQ.'EXTRAP') THEN
C                            ==                         ======
            IF (NPAR .EQ. 2) THEN
              NOEXTRAP = .TRUE.
            ELSE IF (ACTPAR(3) .EQ. 'UNTIL' .AND. 
     >               ACTPAR(4) .EQ. 'JOB') THEN
              READ (ACTPAR(5), '(F20.0)') XL
              ISTART=IFIX(XL)
              IF (JOBNUM .LE. ISTART) NOEXTRAP=.TRUE.
            ELSE IF (ACTPAR(3) .EQ. 'UNTIL' .AND. 
     >               ACTPAR(4) .EQ. 'CORR') THEN
              READ (ACTPAR(5), '(F20.0)') COREX
              BCOREX = .TRUE.
            ENDIF

      ELSE IF (ACTPAR(1) .EQ. 'NSCHAR') THEN
C                              ======
            READ(ACTPAR(2),'(I10)',ERR=1) NSCHAR
 
      ELSE IF (ACTPAR(1) .EQ. 'ITMAX') THEN
C                              =====
            READ(ACTPAR(2),'(I10)',ERR=1) ITMAX
            
      ELSE IF (ACTPAR(1)(1:7) == 'FRACINV') THEN
C                                 =======
            bFRACINV = .TRUE.
      ELSE IF (ACTPAR(1) == 'SPLITINVERSION') THEN
C                            ==============
            iBLOCKINVERSION = 1
            IF (NPAR == 2) THEN
               IF (ACTPAR(2) == 'ION') THEN
                 iBLOCKINVERSION = 2
               ELSEIF (ACTPAR(2) == 'ATOM') THEN
                 iBLOCKINVERSION = 1
               ELSE
                  READ (ACTPAR(2), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                      iBLOCKINVERSION = IFIX(tempREAL)
                  ENDIF
               ENDIF
            ENDIF            

      ELSE IF (ACTPAR(1) == 'ZERORATES') THEN
C                            =========
C     @TODO:  Check for good keywords for the particular methods
            IF (NPAR == 2) THEN
               IF (ACTPAR(2) == 'NONE') THEN
                 iZRType = 0
               ELSEIF (ACTPAR(2) == 'PREDICT') THEN
                 iZRType = 1
               ELSEIF (ACTPAR(2) == 'LADDER') THEN
                 iZRType = 2
               ELSE
                  READ (ACTPAR(2), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                      iZRType = IFIX(tempREAL)
                  ENDIF
               ENDIF
            ENDIF                  

      ELSE IF (ACTPAR(1) == 'CMMODE') THEN
C                            ======
C     @TODO:  Check for good keywords for the particular methods
            IF (NPAR == 2) THEN
               IF (ACTPAR(2) == 'ALL') THEN
                 ICMMODE = 0
               ELSEIF (ACTPAR(2) == 'CONVERGED') THEN
                 ICMMODE = 1
C               ELSEIF (ACTPAR(2) == 'LADDER') THEN
C                 iZRType = 2
               ELSE
                  READ (ACTPAR(2), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                      ICMMODE = IFIX(tempREAL)
                  ENDIF
               ENDIF
            ENDIF                  
            
      ELSE IF (ACTPAR(1) .EQ. 'NO' .AND. ACTPAR(2) .EQ. 'OVERLAP') THEN
C                              ==                        =======
            NOLAP = .TRUE.

      ELSE IF (ACTPAR(1) .EQ. 'OLDSTART') THEN
C                              ========
            IF (JOBNUM == 1) OLDSTART = .TRUE.

      ELSE IF ((KARTE(:5) == 'OLD V') .OR.
C                             =====
     >         (KARTE(:9) == 'OLD STRAT')) THEN
C                             =========
            IF (JOBNUM == 1) bOLDVELO = .TRUE.

      ELSE IF (ACTPAR(1)(:5) .EQ. 'LEVEL') THEN
C                              ========
            NEWATOM = .TRUE.

      ELSE IF (ACTPAR(1) == 'LINPOP_THRESHOLD' ) THEN
C                            ================
            READ(ACTPAR(2),'(F10.0)',ERR=1) EPSILON
            IF (IERR == 0) THEN
                NKONV_THRESH = IFIX(tempREAL)
            ENDIF            

      ELSE IF ( (ACTPAR(1)(1:17) == 'STOP_AUTO_MODIFY') .OR. 
C                                    ================
     >          (ACTPAR(1)(1:10) == 'STOP_NOCON' ) ) THEN
C                                    ==========
            BAUTO_ABORT = .TRUE.
            
      ELSE IF (ACTPAR(1)(1:11) .EQ. 'NO_AUTO_MOD' ) THEN
C                                    ===========
            BAUTO = .FALSE.

      ELSE IF ( (ACTPAR(1)(1:11) == 'AUTO_MODIFY') .OR.
C                                    ===========
     >          (KARTE(1:11) == 'AUTO MODIFY') ) THEN
C                                ===========
            BAUTO = .TRUE.
            IF (ACTPAR(2) == 'MODIFY') THEN
              CURPAR = ACTPAR(3)
            ELSE
              CURPAR = ACTPAR(2)
            ENDIF
            SELECTCASE (CURPAR)
              CASE ('TEMP', 'TCOR', 'TCORR', 'TC')
                iAMT = 1
              CASE ('NOTEMP')
                iAMT = -1
            ENDSELECT
            
      ELSE IF (ACTPAR(1)(1:8) .EQ. 'SMALLPOP' ) THEN
C                                   ========
         READ (ACTPAR(2),'(F10.0)', ERR=1) SMALLPOP2
         IF (SMALLPOP2 .NE. SMALLPOP) THEN
            WRITE (*,'(A, G12.3, A, G12.3)') 
     >        'STEAL/DECSTE: WARNING: SMALLPOP=', SMALLPOP2, 
     >        ' chosen different from recommended default', SMALLPOP 
         ENDIF
         SMALLPOP = SMALLPOP2
C***     Optional parameter for outermost depth point considered in CORMAX
C***     (This is passed to STEAL->PRICORR)   added by ansander, 05-Oct-2016
         IF (ACTPAR(3) == 'NDOUT') THEN
           IF (ACTPAR(4) /= 'AUTO') THEN
             READ (ACTPAR(4),'(F10.0)', ERR=1, IOSTAT=IERR) tempREAL
             IF (IERR == 0) THEN
               NDOUT = tempREAL
             ENDIF
           ENDIF
         ENDIF

      ELSE IF (ACTPAR(1)(1:6) .EQ. 'POPMIN' ) THEN
C                                   =======
            READ (ACTPAR(2),'(F10.0)', ERR=1) POPMIN

      ELSE IF (ACTPAR(1)(1:10) .EQ. 'PLOT_INBOX' ) THEN
C                                    ==========
            BINBOX = .TRUE.

      ELSE IF (ACTPAR(1)(1:13) .EQ. 'PLOT_POPRANGE' ) THEN
C                                    =============
            READ (ACTPAR(2),'(F10.0)', ERR=1) POPRANG
            POPRANG = ALOG10(POPRANG)

      ELSE IF (KARTE(:22) == 'NO CONTINUUM ITERATION') THEN
C                             ======================
            BITCONT = .FALSE.

      ELSE IF (ACTPAR(1) == 'TDIFFCOR') THEN
C                              =======
            fTNDCOR = 1.
            IF (NPAR > 1) THEN
              READ (ACTPAR(2), '(F10.0)', IOSTAT=IERR) tempREAL
              IF (IERR == 0) THEN
                  fTNDCOR = tempREAL
              ENDIF
            ENDIF
            
      ELSE IF (ACTPAR(1) == 'TDIFFUS' .AND. ACTPAR(2) == 'ALWAYS') THEN
C                            =======
C***        Perform TDIFFUS corrections even if no other temperature
C***        corrections are applied
            bForceTDIFFUS = .TRUE.

      ELSE IF (ACTPAR(1) == 'NOTDIFFUS') THEN
C                              =======
            bTDIFFUS = .FALSE.

      ELSE IF (ACTPAR(1) .EQ. 'OPC') THEN
C                              ===
         OPC = ACTPAR(2)

      ELSE IF (ACTPAR(1) .EQ. 'TCORR' .AND. 
     >         ACTPAR(2) .EQ. 'ALTERNATE') THEN
C                              =========
         BTALTER = .TRUE.
         
      ELSE IF (ACTPAR(1) == 'TCORR' .AND. 
     >         ACTPAR(2) == 'ADVECTION') THEN
C                            =========
         bINCADV = .TRUE.

      ELSE IF (ACTPAR(1) == 'TCORR' .AND. 
C                            =====
     >         ACTPAR(2) == 'FERATES' .AND.
C                            =======
     >         ACTPAR(3) == 'OFF') THEN
C                            ===
         bNoFeTCORR = .TRUE.
         
      ELSE IF (ACTPAR(1) .EQ. 'RUDLINES_ZERORADRATES') THEN
C                              =====================
         BRUDZERO = .TRUE.

      ELSE IF (ACTPAR(1) .EQ. 'XJLAPP') THEN
C                              ======
         IF (ACTPAR(2) == 'COLI') THEN
C                          ====         
            IF (iALOentry > 0) THEN
              bLAMAPPCOLI = .TRUE.
            ELSE
              WRITE (hCPR,*) '* WARNING: ALO from COLI requested'
     >            // ' but not found -- fallback to default OPC'
            ENDIF
            DO I=3, NPAR
C***          Optional parameter: Specify ALOMIN = min "TAU" for using ALO
              IF (I < NPAR .AND. ACTPAR(3) == 'MIN') THEN
                READ (ACTPAR(4),'(F10.0)', ERR=1) ALOMIN
              ENDIF
C***          Optional parameter to request using the tri-diagonal ALO
              IF (ACTPAR(I)(:3) == 'TRI') THEN
                IF (iALOEntry < 3) THEN
                  WRITE (hCPR,*) '* WARNING: Tri-diagonal ALO requested'
     >            // ' but not found -- applying diagonal ALO instead'
                ELSE
                  bALOTri = .TRUE.
                ENDIF
              ENDIF
            ENDDO

         ELSE IF (ACTPAR(2) .EQ. 'NEW') THEN
C                                 ===
            BXJLAPPNEW = .TRUE.
         ELSE IF (ACTPAR(2) .EQ. 'CORE') THEN
C                                 ====
            BXJLAPPNEW  = .TRUE.
            BXJLAPPCORE = .TRUE.
         ENDIF

         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=1) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=1) XLAM_FINE_END
         ENDIF

      ELSE IF (ACTPAR(1) .EQ. 'XJCAPP' .AND. 
     >         ACTPAR(2) .EQ. 'NEW') THEN
C                              ===========
         IF (BXJLAPPCORE) THEN
          STOP
     >    'DECSTE: DO NOT USE "XJLAPP CORE" .AND. "XJCAPP NEW" TOGETHER'
         ENDIF
         BXJCAPPNEW = .TRUE.
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=1) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=1) XLAM_FINE_END
         ENDIF

      ELSE IF (ACTPAR(1) .EQ. 'XJAPP' .AND. 
     >         ACTPAR(2) .EQ. 'NEWOPERATOR') THEN
C                              ============
         BNEWOPER = .TRUE.

      ELSE IF (ACTPAR(1) .EQ. 'NOBNUEFE') THEN
C                              ========
         BNUEFE = .FALSE.

      ELSE IF (ACTPAR(1) == 'COLIPP') THEN
C                            ======
         bUCPP = .TRUE.
         
      ELSE IF (ACTPAR(1) == 'SCHARMERACC') THEN
C                            ===========

         IF (ACTPAR(2) == 'ABSOLUTE') THEN
C                          ========
           bRELSCHARMERACC = .FALSE.
           
         ELSEIF (ACTPAR(2) == 'RELATIVE') THEN
C                              ========
           bRELSCHARMERACC = .TRUE.
         ENDIF
         
      ELSE IF (ACTPAR(1) .EQ. 'TRACE' .AND. 
C                              =====
     >         ACTPAR(2) .EQ. 'POPNUM') THEN
C                              ======
         BTRACE_POPNUM = .TRUE.

      ELSE IF (ACTPAR(1) .EQ. 'WJCMIN') THEN
C                              ======
            READ(ACTPAR(2),'(E20.10)',ERR=1) WJCMIN

      ELSE IF (KARTE(:5) == 'RGRID') THEN
C                            =====
         RadiusGridParameters(1) = KARTE

      ELSE IF (ACTPAR(1) == 'SPECIAL_OUTER_POINTS') THEN
         RadiusGridParameters(2) = KARTE
         
      ELSE IF (ACTPAR(1) == 'SPECIAL_INNER_POINTS') THEN
         RadiusGridParameters(3) = KARTE         


      ELSE IF (KARTE(:6) == 'VELPAR') THEN
C                            ======
         VMINCAND = VMIN
         CALL DECVELPAR(KARTE, VFINAL, VMINCAND, BETA, RMAX)
         IF (VMIN < 0) THEN
C***        use VMIN from CARDS only if not stored in MODEL
            VMIN = VMINCAND
         ELSE
            VMINCAND = VMIN
         ENDIF

      ELSE IF (ACTPAR(1) == '2BETALAW') THEN
C                            ========
         READ (ACTPAR(3), '(F10.0)', ERR=1) BETA2 
         READ (ACTPAR(5), '(F10.0)', ERR=1) BETA2FRACTION  

      ELSE IF ((ACTPAR(1) == 'THIN') .OR.
C                             ====
     >         (KARTE(:9) == 'HYDROSTAT')) THEN
C                             =========
         bThinWind = .TRUE.
         ThinCard = KARTE
         IF (NPAR > 2) THEN
           DO i=3, NPAR
             SELECTCASE (ACTPAR(i))
               CASE ('LOOSE')
                 bHYSTloose = .TRUE.
               CASE ('FULL', 'FULLHD')
                 bFULLHYDROSTAT = .TRUE.
               CASE ('MEAN')
                 bGAMMARADMEAN = .TRUE.
               CASE ('REDUCEGAMMA', 'REDGAM', 'REDUCE', 'RG')
                 GEDDreduce = 0.1
                 IF (NPAR >= (i+1)) THEN
                   READ (ACTPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                   IF (IERR == 0) THEN
                     IF (tempREAL > 1. .OR. tempREAL < 0.)  THEN
                       WRITE (hCPR,'(A)') ' *** Error: invalid choice'
     >                    // ' of REDUCE parameter, must be '
     >                    // ' between 0 and 1!'
                       GOTO 1
                     ENDIF
                     GEDDreduce = tempREAL
                   ENDIF
                 ENDIF
               CASE ('ACC', 'EPS')
                 IF (NPAR >= (i+1)) THEN
                   READ (ACTPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                   IF (IERR == 0) THEN
                     HYSTACC = tempREAL
                   ENDIF          
                 ENDIF
               CASE ('NOCONVCRIT')
                 IHSSTATUS = -2
             ENDSELECT
           ENDDO
         ENDIF
         
C***     Parse accuracy input for hydrostatic equation
C***       negative value = not required for model convergence
C***       positive value = required for model convergence
C***     In addition, prevent unrealistically small accuracies  
         IF (HYSTACC > 0. .AND. IHSSTATUS > -2) THEN
           IHSSTATUS = 0
         ELSEIF (HYSTACC < 0. .OR. IHSSTATUS == -2) THEN
           IHSSTATUS = -1
         ENDIF
         HYSTACC = ABS(HYSTACC)
         IF (HYSTACC < HYSTACCMIN) THEN
           HYSTACC = HYSTACCMIN
           WRITE (hCPR,'(A,F8.5)') ' *** WARNING: Given accuracy for'
     >     // ' HYDROSTATIC EQUATION was too small and has been reset'
     >     // ' to ', HYSTACCMIN
         ENDIF
         
      ELSE IF (ACTPAR(1) == 'RMAX_IN_RSUN') THEN
         RMAX_IN_RSUN = .TRUE.

      ELSE IF (ACTPAR(1) == 'TAUMAX') THEN
C                            ======
         READ (ACTPAR(2), '(F10.0)', ERR=1) TAUMAX
         TaumaxCard = KARTE
C***     Only parameters which are also required for other routines
C***     than just ENSURETAUMAX are decoded here. The rest, such as 
C***     accuracy settings or numerical paramers, is decoded
C***     in ENSURETAUMAX itself.
         IF (NPAR > 2) THEN
          DO i=3, NPAR 
            SELECTCASE (ACTPAR(i))
              CASE ('FIX')
                bENSURETAUMAX = .TRUE.
                IF (NPAR >= (i+1)) THEN
                  READ (ACTPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF
                ENDIF
              CASE ('MIN')
                bTauStrict = .FALSE.
                bENSURETAUMAX = .TRUE.
              CASE ('EPS', 'ACC')
                IF (NPAR >= (i+1)) THEN
                  READ (ACTPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF                  
                ENDIF
              CASE ('REPS', 'RELEPS', 'RELACC')
C***            allows a relative specification of the accuracy
                IF (NPAR >= (i+1)) THEN
                  READ (ACTPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL * TAUMAX
                  ENDIF                  
                ENDIF
              CASE ('THOM', 'THOMSON')
                bTMTHOM = .TRUE.               
            ENDSELECT
          ENDDO

         ENDIF

      ELSE IF (ACTPAR(1) .EQ. 'TAUFIX') THEN
C                              ======
C     note: Deprecated CARD => use TAUMAX card
        WRITE (hCPR,*) '*** ERROR: DEPRECARED TAUFIX CARD FOUND ***'
        WRITE (hCPR,*) '--- Use TAUMAX card options instead! ---'
        STOP 'FATAL ERROR in DESCTE'
        
      ELSE IF (ACTPAR(1) .EQ. 'DENSCON') THEN
C                              =======         
         READ (ACTPAR(2), '(F10.0)',ERR=1) DENSCON_FIX
         DENSCON_LINE = KARTE

      ELSE IF (ACTPAR(1) == 'MACROCLUMP') THEN
C                            ==========         
         READ (ACTPAR(2), '(F10.0)',ERR=1) CLUMP_SEP
         MacroCard = KARTE

      ELSE IF (ACTPAR(1) == 'OPALINE') THEN
C                            =======         
         READ (ACTPAR(2), '(F10.0)',ERR=1) OPALINE_SCALE
         
      ELSE IF (ACTPAR(1) == 'HYDRO') THEN
C                            =====
        bHYDROSOLVE = .TRUE.
        HydroCard = TRIM(HydroCard) // KARTE(6:)
        DO I=2, NPAR
          SELECTCASE (ACTPAR(I))
            CASE ('LATETAU', 'LTAU', 'LT')
              bLateTau = .TRUE.
            CASE ('ANALYSE')
              bHydroHelp = .TRUE.
          ENDSELECT
        ENDDO
    
      ELSE IF (ACTPAR(1) .EQ. 'ALPHA') THEN
C                              =====
        AlphaCard = KARTE

      ELSE IF (ACTPAR(1) == 'WRTYPE' .OR. ACTPAR(1) == 'STARTYPE') THEN
        CALL SARGV (KARTE, 2, WRTYPE)
        IF (WRTYPE /= 'OB' .AND. WRTYPE /= 'WN' .AND.
     >      WRTYPE /= 'WC') THEN
            WRITE (0, *) '*** ERROR: Invalid choice of WRTYPE' 
            WRITE (0, *) '*** Valid types are: OB, WN, WC'
            STOP 'ERROR IN DECSTE'
        ENDIF

      ELSE IF (ACTPAR(1) == 'MLANGER') THEN
C                            =======
        MFORM = 1
      ELSE IF (ACTPAR(1) == 'MGOETZ') THEN
C                            =======
        MFORM = 2
      ELSE IF (ACTPAR(1) == 'MSTAR') THEN
C***    Fallback if MSTAR was not in model file
        IF ((XMSTAR < -90.) .AND. (GLOG < -90.)) THEN
          READ (ACTPAR(2), '(F10.0)') XMSTAR
          IF (XMSTAR <= 0.) THEN
            CALL REMARK('NO MASS INFORMATION FOUND IN CARDS')
            STOP 'ERROR IN STEAL'
          ELSE
            GLOG = ALOG10(GCONST * XMSTAR * XMSUN / RSTAR / RSTAR)
          ENDIF
        ENDIF
      ELSE IF (ACTPAR(1) == 'VTURB' .OR. ACTPAR(1) == 'VMIC') THEN
C                            =====                     ====      
        VTURB_LINE = KARTE
        
c      ELSE IF (ACTPAR(1) == 'VTURBUPDATE') THEN
c        bUpdateVT = .TRUE.
c        READ (ACTPAR(2), '(F10.0)', ERR=1) VTURBND
c        DO L=1, ND
c          VTURB(L) = VTURBND
c        ENDDO
      ELSE IF (ACTPAR(1) == 'MASSUPDATE') THEN
C                            ==========
C       !Changes the mass on the fly (either direct input or m-r relation)
C       !This is a debug option that should only be used for tests
        bUpdateMass = .TRUE.
        IF (NPAR > 1) THEN
          READ (ACTPAR(2), '(F10.0)', ERR=1) XMSTAR          
          GLOG = ALOG10(GCONST * XMSTAR * XMSUN / RSTAR / RSTAR)
        ELSE
          XMSTAR = 0.
        ENDIF
      ELSE IF (ACTPAR(1) == 'RSTARUPDATE') THEN
C                            ===========
C***    Manual change of RSTAR on the fly 
C***    This is only a debug option that should be used with care
C***    Input is expected in solar radii
        IF (NPAR > 1) THEN
          READ (ACTPAR(2), '(F10.0)', ERR=1) RSTAR
          RSTAR = RSTAR * RSUN
          GLOG = ALOG10(GCONST * XMSTAR * XMSUN / RSTAR / RSTAR)
C***      Scale total density with new radius
          DO L=1, ND
            ENTOT(L) = ENTOT(L) * (RSTARorg/RSTAR)**2
            RHO(L) = RHO(L) * (RSTARorg/RSTAR)**2
          ENDDO           
          bUpdateMass = .TRUE.
        ELSE
          WRITE (0,*) '*** ERROR: RSTARUPDATE NEEDS A PARAMETER'
          STOP '*** FATAL ERROR IN STEAL->DECSTE'
        ENDIF        

      ELSE IF (KARTE(:8) == 'HEADLINE') THEN
C                            ========
        HEADLINE = KARTE(10:)        
      ELSE IF (KARTE(:14) == 'RENEW HEADLINE') THEN
C                             ==============
        bModHeadUpdate = .TRUE.
      ELSE IF (ACTPAR(1) == 'DEBUGTB') THEN
C                            =======
        bKUBATDEBUG = .TRUE.
      ELSE IF (ACTPAR(1) == 'BACKUP') THEN
C                            ======
        NBACKUP = 50
        IF (NPAR > 1) THEN
          READ (ACTPAR(2), '(I10)', ERR=1) NBACKUP
        ENDIF
      ELSE IF (KARTE(:4) == 'VDOP') THEN
C                            ====
         IF (NPAR > 2) THEN
           IF (ACTPAR(3) == 'AUTO') bDDVDOP = .TRUE.
         ENDIF
      ELSE
C***  Check if the card refers to an element name
        CALL FINDCHARGE (ACTPAR(1), NZ)
C***    NZ > 0 means that an element of this name is known
        IF (NZ > 0) THEN
c          IFOUNDELEMENT = 0
          DO K=1, NATOM
            IF (ACTPAR(1) == ELEMENT(K)) THEN      
C***          ELEMENT IS NOT ONLY IN CARDS, BUT ALSO IN DATOM            
c              IFOUNDELEMENT = 1
C***          Look for keyword 'TRACE' -> this element will not be 
C***           considered for max/min correction calculation
              DO I=3, NPAR
                IF (ACTPAR(I) == 'TRACE') THEN
                  TRACEELEM(K) = .TRUE.
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      GOTO 8

C***  ERROR EXIT **************************************************
    1 CONTINUE
      WRITE (hCPR,*)'DECSTE: ERROR WHILE DECODING THE FOLLOWING LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'

C***  REGULAR EXIT: END-OF-FILE REACHED  **************************
   99 CONTINUE
      CLOSE (1)

      IF (RadiusGridParameters(1) == ' ' .AND. (.NOT. STHLP)) THEN
          WRITE (hCPR,'(A)') 'WARNING: RADIUS GRID NOT SPECIFIED' 
          !no grid is allowed if neither TAUMAX or HYDRO routine is called
          IF (bHYDROSOLVE .OR. bENSURETAUMAX) STOP 'ERROR'
      ENDIF

      IF (bHYDROSOLVE .OR. bENSURETAUMAX) THEN
        IF (RMAX < .0) THEN
          WRITE (hCPR,'(A)') 'RMAX NOT SPECIFIED'
          STOP 'ERROR'
        ENDIF
        IF (RMAX_IN_RSUN) THEN
            RMAX = RMAX / ( RSTAR / RSUN )
        ENDIF
      ENDIF

C***  Unsoeld-Lucy temperature correction UNLUTEC suppressed 
C**   ...   by the NOTEMP switch (note: NOTEMP only local in this routine!)
      IF (NOTEMP) BUNLU = .FALSE.
C**   ...   for WRSTART 
      IF (JOBNUM .LT. 15) BUNLU = .FALSE.
      IF (JOBNUM .LT. 15) NOTEMP = .TRUE.
C**   ...   for STEAL-Help (output only job)
      IF (STHLP) BUNLU  = .FALSE.
      IF (STHLP) NOTEMP = .TRUE.

      IF (.NOT. NOTEMP .AND. .NOT. BUNLU) THEN 
         WRITE (hCPR,*) '*** UNLUTEC Parameters missing!'
         STOP '*** FATAL ERROR detected in DECSTE'
      ENDIF
            
C***  Consistency check between GEFF keyword and HYDROSTATIC INTEGRATION
      IF (bThinWind .AND. 
     >       TRIM(GEFFKEY) /= '' .AND. TRIM(GEFFKEY) /= 'AUTO') THEN
        IF (bFULLHYDROSTAT .AND. TRIM(GEFFKEY) == 'THOM') THEN
          WRITE(hCPR,*) 'DECSTE: Inconsistency between HYDROSTATIC '
     >      // 'INTEGRATION card and intended meaning of LOG GEFF'
          WRITE(hCPR,*) 'You cannot use the FULL option as LOG GEFF '
     >      // 'refers only to the electron Gamma.'
          WRITE(hCPR,*) 'Remove the FULL option or restart the model '
          WRITE(hCPR,*) 'with RADFORCE=FULL on the LOG GEFF card.'
          STOP '*** FATAL ERROR in DECSTE ***'
        ELSEIF (.NOT. bFULLHYDROSTAT .AND. TRIM(GEFFKEY) == 'RAD') THEN
          WRITE(hCPR,*) 'DECSTE: Inconsistency between HYDROSTATIC '
     >      // 'INTEGRATION card and intended meaning of LOG GEFF'
          WRITE(hCPR,*) 'You have to use the FULL option as LOG GEFF '
     >      // 'refers to the full radiative Gamma.'
          WRITE(hCPR,*) 'Add the FULL option or restart the model with'
          WRITE(hCPR,*) 'RADFORCE=ELECTRON on the LOG GEFF card.'
          STOP '*** FATAL ERROR in DECSTE ***'
        ENDIF
      ENDIF
      

      IF (bModHeadUpdate) THEN
        MODHEAD(35:) = HEADLINE
      ENDIF
      
C*** If there are LEVEL cards (for ADAPTER), the first STEAL (JOBNUM=1)
C***  should not be suppressed, in order to provide reasonable POPNUMs
C***  for those levels which are not taken from OLD MODEL
      IF (NEWATOM) OLDSTART = .FALSE.

      RETURN
      END
