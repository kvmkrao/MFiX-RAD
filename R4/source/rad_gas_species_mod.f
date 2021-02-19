!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION GAS SPECIES                                  !
!                                                                      !
!  Purpose:  to find gas phase information                             !
!                                                                      !
!  Authors: V Kotteda, M Stoellinger                   Date: 1-july-19  !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module rad_gas_species
use rad_param
use rad_config
use rad_fields
use rxns, only : species_g, dim_m, dim_n_all, dim_n_g
use fldvar, only : ep_g, T_g, P_g, X_g
implicit none
private

! ---- Data types -----
type gasInfoType
    real(dp) :: p ! pressure unit bar = 10^5 N/m^2 (do the conversion based on MFIX unit choice)
    real(dp) :: T ! temperature unit Kelvin
    real(dp), dimension(nSp) :: C ! mole fractions
    real(dp) :: volFrac ! gas phase volume fraction
!    real(dp) :: sootfv ! soot volume fractions unit ppm
    real(dp) :: splen   
end type gasInfoType

! ---- Interfaces ------
public :: gasInfoType
public :: rad_gas_species_init
public :: getGasInfo

interface rad_gas_species_init
    module procedure init1
end interface

! ---- Data members ----
! -- species data --
integer, dimension(DIM_N_G) :: radSpeciesId

logical, dimension(DIM_N_G) :: isRadSpecies


contains

subroutine init1(config)
    implicit none
    type(configuration), intent(in) :: config
    integer :: i
    isRadSpecies =.false.
    radSpeciesId = 0
    do i =1, dim_n_g
        select case (trim(species_g(i)))
        case ("CO2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CO2
        case ("H2O")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2O
        case default
            isRadSpecies(i) = .false.
            radSpeciesId(i) = 0
        end select
    end do
end subroutine

subroutine getGasInfo (cell, gasinfo)
    use physprop, only: MW_AVG, MW_g
    use indices, only: i_of, j_of, k_of
    use geometry, only : vol, DX, DY, DZ, AXY,AXZ,AYZ
    implicit none
    integer, intent(in) :: cell
    type(gasInfoType), intent(out) :: gasinfo
    real(dp) :: xg
    integer :: i
    gasinfo%T = T_g(cell)
    gasinfo%P = P_g(cell)
    gasinfo%volFrac = ep_g(cell)
!    gasinfo%splen   = max(dx(i_of(cell)),dy(j_of(cell)),dz(k_of(cell)))
    gasinfo%splen   = 1.8d0*vol(cell)/(AXY(cell)+AXZ(cell)+AYZ(cell)) ! 3.6 vol/area
    do i=1, NMAX(0) ! mass to mole fraction
        if (isRadSpecies(i)) then
	    ! for gas species, gas mole fraction is used
            xg= X_g(cell,i)*MW_MIX_g(cell)/MW_AVG !MW_g(i)
            gasinfo%C(radSpeciesId(i)) = X_g(cell,i) !xg
        endif
    end do
end subroutine

end module rad_gas_species
