!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SOLID INFORMATION                            !
!                                                                      !
!  Purpose:  							       !
!to find solid phase info(radius, void fraction, temperatue etc)       !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module rad_solid_species
use rad_param
use rad_config
use rad_fields
use rxns, only : species_s, dim_m, dim_n_all, dim_n_s
use fldvar, only : ep_s, D_p, T_s, X_s

implicit none
private

! ---- Data types -----
type solidParInfoType
  real(dp) :: radius
  real(dp) :: volFrac
  real(dp) :: T ! temperature
  integer :: totSpecies
  integer, dimension(MaxSolidSpecies) :: SpeciesId
  real(dp), dimension(MaxSolidSpecies) :: SpeciesMassFrac
end type solidParInfoType

! ---- Interfaces ------
public :: solidParInfoType
public :: rad_solid_species_init
public :: getSolidInfo
public :: getDesSolidInfo
public :: smax

interface rad_solid_species_init
    module procedure init1
end interface

! ---- Data members ----
! -- species data --
integer, dimension(DIM_M, DIM_N_S) :: radSpeciesId

logical, dimension(DIM_M, DIM_N_S) :: isRadSpecies


contains

subroutine init1(config)
    use rad_solid_prop, only : speciesNameToId
    implicit none
    type(configuration), intent(in) :: config
    integer :: i,m,sid
    isRadSpecies =.false.
    radSpeciesId = 0
    do m=1,SMAX
        do i =1, NMAX(m)
            sid=speciesNameToId(species_s(m,i), len(species_s(m,i)))
            radSpeciesId(m,i)=sid
            isRadSpecies(m,i) = sid.gt.0
        end do
    end do
end subroutine

subroutine getSolidInfo (cell, m, solidinfo)
    implicit none
    integer, intent(in) :: cell, m
    type(solidParInfoType), intent(out) :: solidinfo
    integer :: i
    solidinfo%volFrac = EP_s(cell, m)
    solidinfo%radius = D_p(cell, m)/2.d0
    solidinfo%T = T_s(cell, m)
    solidinfo%totSpecies = NMAX(m)
    do i=1, NMAX(m) ! mass mole fraction 
        solidinfo%speciesId(i) = radSpeciesId(m,i)
        solidinfo%speciesMassFrac(i) = X_s(cell,m,i)
    end do
end subroutine

subroutine getDesSolidInfo (cell, m, solidinfo)
      USE constant
      USE des_thermo
      USE discretelement
      USE fldvar
      USE param1
      USE toleranc
      use param, only: DIMENSION_N_s
      USE compar
      USE derived_types, only: PIC
      USE functions
      USE geometry
      USE indices
      USE physprop
      USE run, only: ENERGY_EQ
      use des_rxns, only: DES_X_s
      use rad_solid_prop, only : speciesNameToId
      USE particles_in_cell_mod, only: particles_in_cell
 
    implicit none
    integer, intent(in) :: cell, m
    type(solidParInfoType), intent(out) :: solidinfo
    integer :: i,lNP,NP

    INTEGER, PARAMETER :: lUpdateFreq=5
    INTEGER, SAVE :: lUpdate_avgTs=0 
    real(dp) :: sumrad,sum_T_s, Tenv
    real(dp) :: sumxs(NMAX(m))
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


    solidinfo%volFrac = EP_s(cell,m) 
    sumrad  = 0.0d0
    sum_T_s = 0.0d0
    sumxs   = 0.0d0 

    DO lNP = 1, PINC(cell)
          NP = PIC(cell)%p(lNP)
!         Update the sum of particle temperatures in fluid cell IJK.
          IF(IS_NORMAL(NP)) then 
            SUM_T_s = SUM_T_s + DES_T_s(NP)
            sumrad  = sumrad  + DES_RADIUS(NP) 
            do i=1,NMAX(m)
               sumxs(i)   = DES_X_s(NP,i) !sumxs(i)  + DES_X_s(NP,i)
            end do  
          END IF 
    ENDDO
 
    solidinfo%radius  = sumrad/max(1,PINC(cell)) 
    solidinfo%T       = SUM_T_s/max(1,PINC(cell)) 
    solidinfo%totSpecies = NMAX(m)
    do i=1, NMAX(m) 
        solidinfo%speciesId(i) = speciesNameToId(species_s(m,i), len(species_s(m,i))) 
        solidinfo%speciesMassFrac(i) =  sumxs(i) 
    end do
end subroutine

SUBROUTINE CALC_avgTs

      USE compar
      USE derived_types, only: PIC
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE functions
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE run, only: ENERGY_EQ

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

end module
