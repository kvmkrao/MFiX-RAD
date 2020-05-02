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
      RAD_EMIS_W = UNDEFINED
      RAD_T_W = UNDEFINED
      RAD_NQUAD = 1
      RAD_SKIP = 100
      RAD_NRR=0
      RAD_RTE='P1'
      RAD_SPECTRAL='GRAY'
      RAD_WSGG='SMITH82'

      const_sfp  = ZERO
      const_emisp= ZERO
      const_absg = ZERO
       
      RETURN
      END SUBROUTINE USR_INIT_NAMELIST 
