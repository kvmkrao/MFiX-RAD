MODULE CALC_THERMO_DES_MOD

   USE DES_RADIATION_MOD, ONLY: CALC_AVGTS, DES_RADIATION
   USE RXNS_GS_DES1_MOD, ONLY: RXNS_GS_DES1

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_THERMO_DES                                        !
!                                                                      !
!  Purpose: This subroutine is called from DES routines. It calls      !
!  functions that calculate heat and mass transfer.                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_THERMO_DES

      USE compar
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE drag_gs_des1_mod, only: drag_gs_des1
      use conv_gs_des1_mod, only: conv_gs_des1
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE interpolation
      USE param1
      USE run

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: NP
! Functions
!---------------------------------------------------------------------//

! This is a quick work-around to keep the thermo routines from causes
! issues while the "check_data" routines are rewritten. Moving forward
! this routine should be split apart to avoid the particle loops for
! cold-flow, non-reacting cases.
      IF(.NOT.ENERGY_EQ .AND. .NOT.ANY_SPECIES_EQ) RETURN

! Calculate time dependent physical properties
      FORALL(NP=1:MAX_PIP, PARTICLE_STATE(NP)==NORMAL_PARTICLE) &
         DES_C_PS(NP) = CALC_CP_DES(NP)

      IF(DES_EXPLICITLY_COUPLED) THEN
!! Apply the convective heat transfer calculated by the gas phase.
         IF(CALC_CONV_DES) THEN
            WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
               Q_Source = Q_Source + CONV_Qs
         ENDIF
         IF(ANY_SPECIES_EQ) THEN
            WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
               Q_Source = Q_Source + RXNS_Qs
         ENDIF
         CALL DES_RADIATION
      ELSE
         IF(CALC_CONV_DES)CALL CONV_GS_DES1
         CALL DES_RADIATION
         IF(ANY_SPECIES_EQ) CALL RXNS_GS_DES1
      ENDIF

   END SUBROUTINE CALC_THERMO_DES

END MODULE CALC_THERMO_DES_MOD
