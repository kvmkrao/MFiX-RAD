!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION FIELDS                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 1-July-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

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

integer :: mphas 
real(dp), dimension(:), allocatable :: gq 
! quadrature points

real(dp), dimension(:), allocatable :: wq
! quadrature weights

! -- spectral fields data --
real(dp), dimension(:,:), allocatable :: k_g
real(dp), dimension(:,:), allocatable :: a_g
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
INTEGER RAD_NQUAD, RAD_SKIP, RAD_NRR
COMMON / USR_DATA_I / &
          RAD_NQUAD, RAD_SKIP, RAD_NRR

contains

function nQuad() result (nq_)
    implicit none
    real(dp) ::  nq_
    nq_ = rad_nquad !nq
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
    call allocateFields2(rad_nquad)
!    call allocateFields2(config%nQuad)
end subroutine

subroutine allocateFields2(nq_)
    use discretelement, only: DES_MMAX, discrete_element
    use mfix_pic, only: MPPIC
    USE physprop, only: smax, mmax
    implicit none
    integer, intent(in) :: nq_
    nq = rad_nquad !nq_


   mphas = smax
   if(discrete_element.or.MPPIC) mphas = DES_MMAX
   If(smax.lt.1) mphas = 1
    allocate (&
        gq(nq), &
        wq(nq), &
        k_g(DIMENSION_3, nq), &
        a_g(DIMENSION_3, nq), &
        k_s(DIMENSION_3, mphas, nq), &
        scat(DIMENSION_3, mphas, nq), &
        G(DIMENSION_3, nq), &
        E(DIMENSION_3, 0:mphas, nq), &
        Srad(DIMENSION_3, 0:mphas) &
        )
    gq=0.d0
    wq=0.d0
    k_g=0.d0
    a_g=0.d0
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
    if (allocated(a_g)) deallocate(a_g)
    if (allocated(k_s)) deallocate(k_s)
    if (allocated(scat)) deallocate(scat)
    if (allocated(G)) deallocate(G)
    if (allocated(E)) deallocate(E)
    if (allocated(Srad)) deallocate(Srad)
end subroutine

end module rad_fields
