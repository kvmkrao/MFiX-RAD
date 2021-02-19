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
        case ("O3")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = O3
        case ("N2O")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = N2O
        case ("CO")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CO
        case ("CH4")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CH4
        case ("O2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = O2
        case ("NO")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = NO 
        case ("SO2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = SO2
        case ("NO2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = NO2
        case ("NH3")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = NH3
        case ("HNO3")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HNO3
        case ("OH")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = OH
        case ("HF")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HF
        case ("HCl")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HCl
        case ("HBr")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HBr
        case ("HI")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HI 
        case ("ClO")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = ClO
        case ("OCS")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = OCS
        case ("H2CO")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2CO
        case ("HOCl")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HOCl
        case ("N2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = N2
        case ("HCN")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HCN
        case ("CH3Cl")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CH3Cl
        case ("H2O2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2O2
        case ("C2H2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = C2H2
        case ("C2H6")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = C2H6
        case ("PH3")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = PH3
        case ("COF2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = COF2
        case ("SF6")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = SF6
        case ("H2S")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2S
        case ("HCOOH")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HCOOH 
        case ("HO2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HO2
        case ("O")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = O
        case ("ClONO20")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = ClONO20
        case ("NO+")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = NO
        case ("HOBr")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HOBr
        case ("C2H4")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = C2H4
        case ("CH3OH")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CH3OH
        case ("CH3Br")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CH3Br
        case ("CH3CN")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CH3CN
        case ("CF4")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CF4
        case ("C4H2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = C4H2
        case ("HC3N")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = HC3N
        case ("H2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2
        case ("CS")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CS
        case ("SO3")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = SO3
        case ("C2N2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = C2N2
        case ("COCl2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = COCl2
!HCl= 15, HBr= 16, HI= 17,   ClO= 18,   OCS= 19,  H2CO= 20, HOCl= 21, & 
!N2= 22,  HCN= 23, CH3Cl= 24,H2O2= 25,  C2H2= 26, C2H6= 27, PH3= 28, & 
!COF2= 29,SF6= 30, H2S= 31,  HCOOH= 32, HO2= 33,  O= 34,    ClONO20= 35, & 
!NO= 36,  HOBr= 37,C2H4= 38, CH3OH= 39, CH3Br= 40,CH3CN= 41,CF4= 42, & 
!C4H2= 43,HC3N= 44,H2= 45,   CS= 46,    SO3= 47,  C2N2= 48, COCl2= 49 
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
    write(*,*) "gasinfo", gasinfo%T, gasinfo%P 
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
