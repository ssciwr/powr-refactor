C*** INTERFACE LIBRARY ***********************
C
C
C
C*********************************************

C***  Interface for SHIFT calls
      INTERFACE SHIFT

        SUBROUTINE SHIFTREAL (RealArray,L,ND)
          INTEGER, INTENT(IN) :: L, ND
          REAL, DIMENSION(ND), INTENT(INOUT) :: RealArray
        END SUBROUTINE SHIFTREAL
        
        SUBROUTINE SHIFTINT (IntegerArray,L,ND)
          INTEGER, INTENT(IN) :: L, ND
          INTEGER, DIMENSION(ND), INTENT(INOUT) :: IntegerArray
        END SUBROUTINE SHIFTINT
        
        SUBROUTINE SHIFTSTRING (StringArray,L,ND)
          INTEGER, INTENT(IN) :: L, ND
          CHARACTER(*), DIMENSION(ND), INTENT(INOUT) :: StringArray
        END SUBROUTINE SHIFTSTRING

      END INTERFACE SHIFT

C***  Interface for unified initialization of beta law parameters
      INTERFACE INITVELBETAPAR
        SUBROUTINE INITVELBETAPAR (RMAX, RCONin, VCON, bUseMaxForPar1)
          REAL, INTENT(IN) :: RMAX, RCONin, VCON
          LOGICAL, INTENT(IN), OPTIONAL :: bUseMaxForPar1          
        END SUBROUTINE INITVELBETAPAR
      END INTERFACE INITVELBETAPAR

C***  Interface for unified SPLINPO calls
      INTERFACE SPLINPO
        SUBROUTINE SPLINPOX(F, X, FI, XI, N, SAFE, LPARAM,DFDX,D2FD2X)
          INTEGER, INTENT(IN) :: N          
          REAL, DIMENSION(N), INTENT(IN) :: XI, FI
          REAL, INTENT(OUT) :: F
          REAL, INTENT(IN) :: X
          LOGICAL, INTENT(IN), OPTIONAL :: SAFE
          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
          REAL, INTENT(OUT), OPTIONAL :: DFDX
          REAL, INTENT(OUT), OPTIONAL :: D2FD2X
        END SUBROUTINE SPLINPOX
      END INTERFACE SPLINPO

      INTERFACE SPLINPO_SAME_X
        SUBROUTINE SPLINPOX_SAME_X(F, X, FI, XI, N, SAFE, LPARAM,
     >                             DFDX,D2FD2X)
          INTEGER, INTENT(IN) :: N          
          REAL, DIMENSION(N), INTENT(IN) :: XI, FI
          REAL, INTENT(OUT) :: F
          REAL, INTENT(IN) :: X
          LOGICAL, INTENT(IN), OPTIONAL :: SAFE
          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
          REAL, INTENT(OUT), OPTIONAL :: DFDX
          REAL, INTENT(OUT), OPTIONAL :: D2FD2X
        END SUBROUTINE SPLINPOX_SAME_X
      END INTERFACE SPLINPO_SAME_X
      
C***  Interface for GRADQUAL calls
      INTERFACE GRADQUAL
        SUBROUTINE GRADQUAL(QUAL, GRAD, GRADCMP, 
     >                      Y, X, YI, XI, ND,
     >                      bSafeMode, LPARAM)
          INTEGER, INTENT(IN) :: ND
          REAL, DIMENSION(ND), INTENT(IN) :: YI, XI
          REAL, INTENT(OUT) :: QUAL, GRADCMP
          REAL, INTENT(IN) :: GRAD, Y, X     
          LOGICAL, INTENT(IN), OPTIONAL :: bSafeMode
          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
        END SUBROUTINE GRADQUAL
      END INTERFACE GRADQUAL
