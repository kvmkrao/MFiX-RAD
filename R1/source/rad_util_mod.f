module rad_util
! utility procedures
use rad_param
implicit none
private

public :: interp1 ! one-dimensional interpolation
public :: interp1es ! equal space one-dimensional interpolation
public :: interp1ess ! equal space one-dimensional interpolation scalar version
public :: ExtrapolateGhostCells ! extrapolate field variables to ghost cells
public :: bbem ! black body emissive power
public :: to_upper ! change a string to upper case

interface bbem
    module procedure bbem_s, bbem_v
end interface
contains

!C------------------------------------------
!C black body emission
!C------------------------------------------
subroutine bbem_s(E, T)
real(dp), intent(in) :: T	
real(dp), intent(out) :: E	

E = sigma * T**4d0 / Pi ! unit W.m^-2
if (UNITS=='CGS') E = E*1.d3 ! convert to erg.cm^-2.s^-1
end subroutine bbem_s

subroutine bbem_v(E, T)
real(dp), dimension(:), intent(in) :: T	
real(dp), dimension(:), intent(out) :: E	

E = sigma * T**4d0 / Pi ! unit W.m^-2
if (UNITS=='CGS') E = E*1.d3 ! convert to erg.cm^-2.s^-1
end subroutine bbem_v

SUBROUTINE interp1(nxy ,xx,yy, ni, xi,yi) 
!!$ given function xx(nxy)<=> yy(nxy), find yi(ni) correspondding to xi(ni)
!!$ suitable for interpolation among g's, and general interpolation
!!$   Input:
!!$     nxy    --  number of points in the original series y(x)
!!$     xx     --  original abscissas
!!$     yy     --  original function values
!!$     ni     --  number of points to be interpolated
!!$     xi     --  abscissas of interpolated points
!!$     yi     --  function values of interpolated points
!!$                                                     ---By Liangyu Wang
!! --- Added Notes on this program by Jian Cai ---------------------------
! (1) It is assumed that xx is in ascendant order, which is typical for
!     linear interpolation
! (2) xi is also assumed to be in ascendant order, which improves search
!     speed, but the original author failed to mension that
! (3) Original treatment for large outrange xi is constant extrapolation
!     in order to stablize calculation
! (4) Original treatment for small outrange xi is linear interpolation 
!     between (xx(1), yy(1)) and (0,0). This may cause error extrapolation
!     for negative xi
! Changes:
! (1) Comment out "close point approximation", because (a) it stops before
!     reaching the closest point; (b) it doesn't save much time at all
! (2) Replace "small interval treatment": instead of lower-bound the denominator
!     use "if" statement to detect small interval
! (3) Move constant extrapolation out of "i" loop and change shortcuts
!! ---End of Added Notes by Jian Cai ------------------------------------

IMPLICIT NONE
INTEGER,INTENT(in) :: nxy, ni
REAL(dp),DIMENSION(nxy),INTENT(in) :: xx, yy
REAL(dp),DIMENSION(ni),INTENT(in) :: xi
REAL(dp),DIMENSION(ni),INTENT(out) :: yi

INTEGER :: i, iq, ibgn
REAL(dp) :: dx

ibgn=1
DO iq= 1, ni
    ! constant extrapolate
    IF (xi(iq) < xx(1)) THEN
        yi(iq) = yy(1); cycle
!          EXIT
    ELSE IF (xi(iq) > xx(nxy)) THEN
        yi(iq:ni)=yy(nxy); exit
    END IF

    DO i=ibgn, nxy
        IF (xx(i) > xi(iq) ) THEN ! this means xx(i-1)<=xi(iq)<xx(i)
            dx = xx(i)-xx(i-1)
            if (abs(dx)<tiny(xx)) then
                yi(iq)=yy(i-1)
            else
                yi(iq)= yy(i-1) + ( yy(i)-yy(i-1) ) * ( xi(iq)-xx(i-1) ) /dx
                ibgn=i
                EXIT
            end if
        END IF
    END DO ! i
END DO ! iq
end subroutine

subroutine interp1es(nxy ,xx,yy, ni, xi,yi )
! given function xx(nxy)<=> yy(nxy), find yi(ni) corresponding to xi(ni)
! ONLY FOR EQUALLY SPACED GRID INTERPOLATION
implicit none
integer,intent(in) :: nxy, ni
real(dp),dimension(nxy),intent(in) :: xx, yy
real(dp),dimension(ni),intent(in) :: xi
real(dp),dimension(ni),intent(out) :: yi
 ! local variables
integer :: iq,i
real(dp) :: wx,dx
 
dx=xx(2)-xx(1)
 
do iq=1,ni
    if (xi(iq)<= xx(1)) then
        yi(iq)=yy(1) ! constant extrapolation
    else if (xi(iq)>= xx(nxy)) then
        yi(iq)=yy(nxy) ! constant extrapolation
    else
        i=int((xi(iq)-xx(1))/dx)+1
        wx=(xi(iq)-xx(i))/dx
        yi(iq)=wx*yy(i+1)+(1-wx)*yy(i)
    end if
end do 
end subroutine

function interp1ess(nxy ,xx,yy, xi) result (yi)
! given function xx(nxy)<=> yy(nxy), find yi(ni) corresponding to xi(ni)
! ONLY FOR EQUALLY SPACED GRID INTERPOLATION
implicit none
integer,intent(in) :: nxy
real(dp),dimension(nxy),intent(in) :: xx, yy
real(dp),intent(in) :: xi
real(dp) :: yi
 ! local variables
integer :: iq,i
real(dp) :: wx,dx
 
dx=xx(2)-xx(1)
 
    if (xi<= xx(1)) then
        yi=yy(1) ! constant extrapolation
    else if (xi>= xx(nxy)) then
        yi=yy(nxy) ! constant extrapolation
    else
        i=int((xi-xx(1))/dx)+1
        wx=(xi-xx(i))/dx
        yi=wx*yy(i+1)+(1-wx)*yy(i)
    end if
end function

subroutine ExtrapolateGhostCells(Phi)
use param
!use fldvar
use functions
use geometry
use indices
implicit none
real(dp), dimension(DIMENSION_3), intent(inout) :: Phi
    
integer :: ijk  	! loop variable over cells
integer :: M 	! loop variable over solid species
! ideally i should have dependency on species and temperature but will put them back later
! try uniform here first
integer :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP  ! neighbor indices
integer :: IJK_1, IJK_2  ! two neighbor indices for linear extrapolation
integer :: ind_1, ind_2, ind_0   ! directional indices of two real cells and ghost cell
! double precision :: coor_x(0:IMAX2), coor_y(0:JMAX2), coor_z(0:KMAX2) ! coordinate
! double precision :: YMIN, ZMIN
double precision :: dl_1, dl_2, dl_0 ! cell sizes at boundary for linear extrapolation
                                  ! dl_0 is ghost cell, dl_1 next to ghost, dl_2 next to dl_1
double precision, parameter :: zero=0.d0
! calculate gas phase absorption coefficients

!INCLUDE 'ep_s1.inc'
!INCLUDE 'function.inc'
!INCLUDE 'ep_s2.inc'

do IJK = 1, DIMENSION_3 ! loop to search ghost cells
  if (.NOT.FLUID_AT(IJK)) then  ! ghost cells
    if (DO_I) then
      IMJK = IM_OF(IJK)
      IPJK = IP_OF(IJK)
    else 
      IMJK = IJK
      IPJK = IJK
    end if
    if (DO_J) then
      IJMK = JM_OF(IJK)
      IJPK = JP_OF(IJK)
    else 
      IJMK = IJK
      IJPK = IJK
    end if

    if (DO_K) then
      IJKM = KM_OF(IJK)
      IJKP = KP_OF(IJK)
    else 
      IJKM = IJK
      IJKP = IJK
    end if

    ! search all neighbors, it is assumed that a ghost cell has only one real 
    ! cell neighbor
    if (DO_I.AND.FLUID_AT(IMJK)) then
       ! ghost cell at east boundary
       IJK_1 = IMJK
       IJK_2 = IM_OF(IJK_1) ! search toward west for second cell
       ind_0 = I_OF(IJK)
       ind_1 = I_OF(IJK_1)
       ind_2 = I_OF(IJK_2)
       dl_0 = DX(ind_0)
       dl_1 = DX(ind_1)
       dl_2 = DX(ind_2)
    else if (DO_I.AND.FLUID_AT(IPJK)) then
       ! ghost cell at west boundary
       IJK_1 = IPJK
       IJK_2 = IP_OF(IJK_1) ! search toward east for second cell
       ind_0 = I_OF(IJK)
       ind_1 = I_OF(IJK_1)
       ind_2 = I_OF(IJK_2)
       dl_0 = DX(ind_0)
       dl_1 = DX(ind_1)
       dl_2 = DX(ind_2)
    else if (DO_J.AND.FLUID_AT(IJMK)) then
       ! ghost cell at south boundary
       IJK_1 = IJMK
       IJK_2 = JM_OF(IJK_1) ! search toward north for second cell
       ind_0 = J_OF(IJK)
       ind_1 = J_OF(IJK_1)
       ind_2 = J_OF(IJK_2)
       dl_0 = DY(ind_0)
       dl_1 = DY(ind_1)
       dl_2 = DY(ind_2)
    else if (DO_J.AND.FLUID_AT(IJPK)) then
       ! ghost cell at north boundary
       IJK_1 = IJPK
       IJK_2 = JP_OF(IJK_1) ! search toward south for second cell
       ind_0 = J_OF(IJK)
       ind_1 = J_OF(IJK_1)
       ind_2 = J_OF(IJK_2)
       dl_0 = DY(ind_0)
       dl_1 = DY(ind_1)
       dl_2 = DY(ind_2)
    else if (DO_K.AND.FLUID_AT(IJKM)) then
       ! ghost cell at bottom boundary
       IJK_1 = IJKM
       IJK_2 = KM_OF(IJK_1) ! search toward top for second cell
       ind_0 = K_OF(IJK)
       ind_1 = K_OF(IJK_1)
       ind_2 = K_OF(IJK_2)
       dl_0 = DZ(ind_0)
       dl_1 = DZ(ind_1)
       dl_2 = DZ(ind_2)
    else if (DO_K.AND.FLUID_AT(IJKP)) then
       ! ghost cell at top boundary
       IJK_1 = IJKP
       IJK_2 = KP_OF(IJK_1) ! search toward bottom for second cell
       ind_0 = K_OF(IJK)
       ind_1 = K_OF(IJK_1)
       ind_2 = K_OF(IJK_2)
       dl_0 = DZ(ind_0)
       dl_1 = DZ(ind_1)
       dl_2 = DZ(ind_2)
    else ! corners
       IJK_1 = IJK
       IJK_2 = IJK
    end if
        
    if (FLUID_AT(IJK_1)) then
      if (FLUID_AT(IJK_2)) then
        ! regular ghost cells found two boundary cells
        Phi(IJK) = linearExtrapolation(Phi(IJK_1), Phi(IJK_2), dl_0, dl_1, dl_2)
      else
        ! search out of boundary
        print*, 'search out of boundary'
        stop
      end if
    else
      if (FLUID_AT(IJK_2)) then
        ! double layer of ghost cells?
        print*, 'Double layer of ghost cells?'
        stop
      else
        !corner cells
        Phi(IJK)= ZERO ! value does not matter
      end if
    end if


  end if ! if ghost cells
end do ! loop to search ghost cells

end subroutine ExtrapolateGhostCells

double precision function  linearExtrapolation(var_1, var_2, &
                dl_0, dl_1, dl_2)
double precision :: var_1, var_2, dl_0, dl_1, dl_2

linearExtrapolation = var_1 - (var_2-var_1)*(dl_0+dl_1)/(dl_1+dl_2)
end function linearExtrapolation


pure function to_upper (str) result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

end function to_upper

end module rad_util
