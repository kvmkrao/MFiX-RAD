!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION TRANSPORT                                    !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!  Commen:                                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


module rad_rte
use rad_param
use rad_config
use rad_spectral
implicit none
private

! ---- Interfaces ------
public :: rad_rte_init
public :: rad_rte_final
public :: rad_rte_calc
public :: rad_rte_modelName

interface rad_rte_init
    module procedure init1
end interface

interface rad_rte_final
    module procedure finalize1
end interface

interface rad_rte_calc
    module procedure calc
end interface

interface rad_rte_modelName
    module procedure modelName
end interface

! ---- Data members ----
character(len=strLen) :: modelName_
integer :: modelId_

! ---- Parameters ------
! model ids
integer, parameter :: unknown=0, p1=1

contains

subroutine init1(config)
    use rad_rte_p1, only : rad_rte_p1_init
    implicit none
    type(configuration), intent(in) :: config
    modelName_ = config%RTEModelName
    select case (modelName_)
        case ("P1")
            modelId_ = p1
            call rad_rte_p1_init(config)
        case default
    end select
end subroutine

subroutine finalize1()
    use rad_rte_p1, only : rad_rte_p1_final
    implicit none
    select case (modelId_)
        case (p1)
            call rad_rte_p1_final
        case default
    end select
end subroutine
subroutine calc()
    use rad_rte_p1, only : rad_rte_p1_calc
    implicit none
    select case (modelId_)
        case (p1)
            call rad_rte_p1_calc
        case default
    end select
end subroutine

function modelName()
    implicit none
    character(len=strLen) :: modelName
    modelName = modelName_
end function
end module rad_rte
