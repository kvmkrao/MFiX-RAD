MODULE DES_RADIATION_MOD

   USE compar, only: ijkstart3, ijkend3
   USE constant, only: pi
   USE derived_types, only: PIC
   USE des_thermo
   USE discretelement
   USE fldvar, only: ep_g, t_g
   USE functions, only: FLUID_AT
   USE functions, only: IS_NORMAL
   USE mfix_pic, only: MPPIC, DES_STAT_WT
   USE param1, only: one, zero
   USE particles_in_cell_mod, only: particles_in_cell
   USE run, only: ENERGY_EQ
   USE run, only: units
   use rad_fields, only: G,k_s
   USE geometry, only: vol
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_RADIATION                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_RADIATION
      use param, only : dimension_bc
      use rad_param, only : sigma 
      IMPLICIT NONE
      include 'usrnlst.inc'
      !character(30) rad_rte, rad_spectral 
      !common / USR_DATA_C / rad_rte, rad_spectral

! Passed variables
!---------------------------------------------------------------------//
! Global index of particle
      INTEGER :: NP
! Solid phase of particle I
      INTEGER :: lM
! Fluid cell index containing particle I
      INTEGER :: IJK

! Local variables
!---------------------------------------------------------------------//
! Environment temperature
      DOUBLE PRECISION :: Tenv
! SB constant TIMES particle surface area
      DOUBLE PRECISION :: SBx4Pi
      DOUBLE PRECISION :: Qrad
      INTEGER, PARAMETER :: lUpdateFreq=5
      INTEGER, SAVE :: lUpdate_avgTs=0
      DOUBLE PRECISION :: sbcon,gfield,abscoe
!......................................................................!


! Coupled simulations update the average solids temperature at the start
! of the DEM time march. Granular flow (no gas) simulations update the
! average solids temperature every "lUpdateFreq" time steps.
      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
         IF(MOD(lUpdate_avgTs,lUpdateFreq) == 0) THEN
            CALL PARTICLES_IN_CELL
            CALL CALC_avgTs
            lUpdate_avgTs = 0
         ELSE
            lUpdate_avgTs = lUpdate_avgTs + 1
         ENDIF
      ENDIF

      SBx4Pi = SB_CONST*4.0d0*Pi

      sbcon  = sigma 

      DO NP=1,MAX_PIP
         IF(IS_NORMAL(NP)) THEN
            lM = PIJK(NP,5)

               IJK = PIJK(NP,4)
               IF(FLUID_AT(IJK)) THEN
                  Tenv = EP_g(IJK)*T_g(IJK) + &
                     (ONE-EP_g(IJK))*avgDES_T_s(IJK)
               ELSE
                  Tenv = avgDES_T_s(IJK)
               ENDIF

               if (RAD_ON .eqv. .False.) then 
                   Qrad = DES_Em(lM)* SBx4Pi * &
                   (DES_RADIUS(NP)**2)*(Tenv**4 - DES_T_s(NP)**4)
               else if (RAD_SPECTRAL.eq.'GRAY-CONST') then
               ! Calculate the heat source/sink due to radiation.
                   Qrad = const_emisp*Pi*(DES_RADIUS(NP)**2) * (G(ijk,1) - 4.0d0*sigma*DES_T_s(NP)**4)

               else if (RAD_SPECTRAL.eq.'GRAY'.or. RAD_SPECTRAL.eq.'GRAY_WSGG') then 
               ! Calculate the heat source/sink due to radiation.
                   abscoe = k_s(ijk,lM,RAD_NQUAD)/PINC(IJK)*vol(ijk)
                   Qrad = abscoe * (G(ijk,1)  - 4.0d0*sbcon*DES_T_s(IJK)**4)
               end if

              !  Update the thermal source term.
              ! bug - needs to check 
               IF(MPPIC) THEN
                  Q_Source(NP) = Q_Source(NP) + Qrad !*DES_STAT_WT(NP)
              !bug - needs to check 
               ELSE
                  Q_Source(NP) = Q_Source(NP) + Qrad
               ENDIF
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE DES_RADIATION

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_avgTs                                             !
!  Author: J.Musser                                   Date: 06-NOV-12  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_avgTs

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER :: IJK
! Loop index for particles.
      INTEGER :: NP, lNP
! Sum of particle temperatures in fluid cell.
      DOUBLE PRECISION :: SUM_T_s
!---------------------------------------------------------------------//

      IF(.NOT.ENERGY_EQ) RETURN

! Loop over fluid cells.
      IJK_LP: DO IJK = IJKSTART3, IJKEND3

         avgDES_T_s(IJK) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP

         IF(PINC(IJK) > 0) THEN
! Initialize local solids temperature.
            SUM_T_s = ZERO
! Loop over all particles in cell IJK.
            lNP_LP: DO lNP = 1, PINC(IJK)
               NP = PIC(IJK)%p(lNP)
! Update the sum of particle temperatures in fluid cell IJK.
               IF(IS_NORMAL(NP)) SUM_T_s = SUM_T_s + DES_T_s(NP)
            ENDDO lNP_LP

! Average solids temperature in fluid cell IJK. The average method
! (over particles) will need changed for Hybrid model (mass? volume?).
            avgDES_T_s(IJK) = SUM_T_s/PINC(IJK)
         ENDIF
      ENDDO IJK_LP

      RETURN

   END SUBROUTINE CALC_avgTs

END MODULE DES_RADIATION_MOD
