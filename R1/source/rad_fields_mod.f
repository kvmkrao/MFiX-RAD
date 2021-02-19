module rad_fields
use rad_param
use rad_config
use rxns, only : ReactionRates, nRR ! to export radiative heat sources
implicit none

! ---- Interfaces ------
private :: init1, allocateFields1, allocateFields2, deallocateFields1

interface rad_fields_init
    module procedure init1
end interface

interface rad_fields_allocation
    module procedure allocateFields1, allocateFields2
end interface

interface rad_fields_deallocation
    module procedure deallocateFields1
end interface

! ---- Data members ----
logical :: radiationOn
! turn on radiation calculation

! -- spectral data --
integer :: nq
! number of spectral quadratures

real(dp), dimension(:), allocatable :: gq 
! quadrature points

real(dp), dimension(:), allocatable :: wq
! quadrature weights

! -- spectral fields data --
real(dp), dimension(:,:), allocatable :: k_g
! gas phase absorption coefficient, 
! dimensions 1=cells, 2=spectral points

real(dp), dimension(:,:,:), allocatable :: k_s
! solid phase absorption coefficient,
! dimensions 1=cells, 2=phases, 3=spectral points

real(dp), dimension(:,:,:), allocatable :: scat
! solid phase scattering coefficient,
! dimensions 1=cells, 2=phases, 3=spectral points

! -- radiation fields data --
real(dp), dimension(:,:), allocatable :: G
! spectral incident radiation, 
! dimensions 1=cells, 2=spectral points

real(dp), dimension(:,:,:), allocatable :: E
! spectral emission, 
! dimensions 1=cells, 2=phases, 2=spectral points

real(dp), dimension(:,:), allocatable :: Srad
! radiative heat sources in energy equations, 
! dimensions 1=cells, 2=phases

real(dp), dimension(dimension_bc) :: ew 
! wall emissivity

real(dp), dimension(dimension_bc) :: Tw
! wall temperature

! ---- Parameters ------
! model ids

contains

function nQuad() result (nq_)
    implicit none
    real(dp) ::  nq_
    nq_ = nq
end function

subroutine init1(config)
    implicit none
    type(configuration), intent(in) :: config
    call allocateFields1(config)
    ! input wall emissivity and temperature
end subroutine

subroutine allocateFields1(config)
    implicit none
    type(configuration), intent(in) :: config
    call allocateFields2(config%nQuad)
end subroutine

subroutine allocateFields2(nq_)
    implicit none
    integer, intent(in) :: nq_
    nq = nq_
    allocate (&
        gq(nq), &
        wq(nq), &
        k_g(DIMENSION_3, nq), &
        k_s(DIMENSION_3, SMAX, nq), &
        scat(DIMENSION_3, SMAX, nq), &
        G(DIMENSION_3, nq), &
        E(DIMENSION_3, 0:SMAX, nq), &
        Srad(DIMENSION_3, 0:SMAX) &
        )
    gq=0.d0
    wq=0.d0
    k_g=0.d0
    k_s=0.d0
    scat=0.d0
    G=0.d0
    E=0.d0
    Srad=0.d0
end subroutine

subroutine deallocateFields1()
    implicit none
    if (allocated(gq)) deallocate(gq)
    if (allocated(wq)) deallocate(wq)
    if (allocated(k_g)) deallocate(k_g)
    if (allocated(k_s)) deallocate(k_s)
    if (allocated(scat)) deallocate(scat)
    if (allocated(G)) deallocate(G)
    if (allocated(E)) deallocate(E)
    if (allocated(Srad)) deallocate(Srad)
end subroutine

end module rad_fields
