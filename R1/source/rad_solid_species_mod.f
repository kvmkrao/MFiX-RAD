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
    do i=1, NMAX(m) ! mass ole fraction
        solidinfo%speciesId(i) = radSpeciesId(m,i)
        solidinfo%speciesMassFrac(i) = X_s(cell,m,i)
    end do
end subroutine

end module
