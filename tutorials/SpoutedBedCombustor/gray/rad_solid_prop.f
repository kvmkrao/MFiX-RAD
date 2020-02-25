!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SOLID PROPERTIRES                            !
!                                                                      !
!  Purpose:  							       !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

module rad_solid_prop
! solid species property database
use rad_param, only : dp
implicit none
integer, parameter, private :: NNU = 248

! ---- parameters -----
integer, parameter :: unknown=0, coal=201, char_=202, ash=203, lime=204
integer, parameter :: fixedCarbon=901, flyAsh=902
contains
function speciesNameToId(speciesName,strLen) result (sid)
    use rad_util, only : to_upper
    implicit none
    integer, intent(in) :: strLen
    character(strLen), intent(in) :: speciesName
    character(strLen) :: spName
    integer :: sid
    spName = to_upper(speciesName)
    select case (trim(spName))
        case ('COAL')
            sid = coal
        case ('CHAR')
            sid = char_
        case ('ASH')
            sid = ash
        case ('LIME')
            sid = 204
        case ('FIXED CARBON')
            sid = fixedCarbon
        case ('COAL ASH')
            sid = flyAsh
        case default
            sid = unknown
    end select
end 

function id2cmplxRefInd(sid, eta) result(m)
! interface function connects to higher levels
! given species id in dababase and wavenumber, return complex refractive index at this wavenumber
    implicit none
    integer, intent(in) :: sid ! species id in database
    real(dp), intent(in) :: eta ! wavenumber
    complex(dp) :: m  ! complex refractive index
    select case (sid)
        case (201) ! coal (use anthracite p395 table11.2 row 2 Modest 2003)
            m = cmplx(2.05d0, -0.54d0, dp)
        case (202) ! char (use carbon)
            m = cmplx(2.20d0, -1.12d0, dp)
        case (203) ! ash
            m = cmplx(1.5d0, -0.02d0, dp)
        !case (204) ! lime
        case (901) ! fixed carbon
            m = cmplx(2.20d0, -1.12d0, dp)
        case (902) ! fly ash
            m = cmplx(1.5d0, -0.02d0, dp)
        case default
            m=0
            print*, 'unknown sid', sid
    end select
end function id2cmplxRefInd

function genNBInterval() result(nbInt)
    real(dp), dimension(NNU,2) :: nbInt ! narrow band intervals
    ! row number equals to the total number of narrow-bands (248)
    ! first/second colume is the start/end of narrow band
    REAL(dp), dimension(5), parameter :: nbRes=(/10.0,25.0,50.0,100.0,250.0/) 
    ! narrowband resolution
    real(dp), dimension(6), parameter :: waveNumSeg = (/200, 300, 4000, 5000, 10000, 15000/) 
    ! narrowband resolution segments
    integer :: i, j, m, n
    real(dp), parameter :: eps = 1d-6  ! a small constant
    real(dp) :: nbl
    n=0
    nbl = waveNumSeg (1) ! narrowband lower bound 
    do i = 1, 5
        m = floor((waveNumSeg(i+1)- waveNumSeg(i))/nbRes(i)+eps) 
        ! total number of narrowbands over this spectral segment
        do j= 1, m
            n = n + 1 ! narrowband index
            nbInt(n, 1) = nbl
            nbl = nbl + nbRes(i)
            nbInt(n, 2) = nbl
        end do !j  
    end do ! i  
end function genNBInterval
end module
