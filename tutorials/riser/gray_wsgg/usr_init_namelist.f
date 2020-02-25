!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_INIT_NAMELIST                                      C
!  Purpose: initialize user_defined NAMELIST variables                 C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR_INIT_NAMELIST
!
      USE param 
      USE param1 
      IMPLICIT NONE
      INCLUDE 'usrnlst.inc' 
!
!
      PAFC = UNDEFINED
      PAA  = UNDEFINED
      ! default values for rad calculation
      RAD_ON = .false.
      RAD_EMIS_W = 0!UNDEFINED
      RAD_T_W = 0!UNDEFINED
      RAD_NQUAD = 1
      RAD_SKIP = 100
      RAD_NRR=0

      const_sfp = ZERO
      const_emisp= ZERO
       
      RETURN
      END SUBROUTINE USR_INIT_NAMELIST 
