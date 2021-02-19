!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION CONFIGURATION                                !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authors: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

module rad_config
use rad_param
implicit none
private

! ---- Data types -----
type configuration
    character(len=strLen) :: SpectralModelName
    character(len=strLen) :: RTEModelName
    integer :: nQuad
    integer :: skipSteps ! currently not functioning
    integer :: nrr
end type configuration

! ---- Interfaces -----
public :: configuration ! radiation configurations
public :: nullconfig
public :: to_upper

contains

subroutine nullconfig(config)
    type(configuration), intent(out) :: config
    config%SpectralModelName=''
    config%RTEModelName=''
    config%nQuad=0
    config%skipSteps=0
    config%nrr=0
end subroutine

subroutine to_upper(config)
    use rad_util, only : upper=>to_upper
    type(configuration), intent(inout) :: config
    config%SpectralModelName = upper(config%SpectralModelName)
    config%RTEModelName = upper(config%RTEModelName)
end subroutine
end module rad_config
