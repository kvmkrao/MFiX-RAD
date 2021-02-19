module rad
use rad_param
use rad_config
use rad_spectral
use rad_rte
use rad_fields, only: Srad, dimension_bc, smax
use rad_fields, only: k_g, k_s, G, E
implicit none
private

public :: rad_init ! inialization
public :: rad_final ! finalization, deallocation variables
public :: rad_calc ! calculate new radiative heat sources
public :: rad_write_src ! export rad sources
public :: S_Rpg, S_Rcg, S_Rps, S_Rcs ! rad heat sources for energy equation

type(configuration) :: config
contains

subroutine rad_init()
    use rad_fields, only: radiationOn, nq, ew, Tw
    implicit none
    include 'usrnlst.inc'
    radiationOn = RAD_ON

    !---M.S.: Here we set the wall bounday values
    Tw = RAD_T_W
    ew = RAD_EMIS_W
    
    ! create configuration data
    call nullconfig(config)
    config%SpectralModelName=RAD_SPECTRAL
    config%RTEModelName=RAD_RTE
    config%nQuad = nq
    config%skipSteps = RAD_SKIP
    config%nrr = RAD_NRR
    if (.not.radiationOn) return
    call to_upper(config)
    call rad_spectral_init(config)
    call rad_rte_init(config)
    
end subroutine rad_init

subroutine rad_final()
    use rad_fields, only: radiationOn
    if (.not.radiationOn) return
    call rad_spectral_final()
    call rad_rte_final()
end subroutine rad_final

subroutine rad_calc()
    use rad_fields, only: radiationOn
    implicit none
    if (.not.radiationOn) return
    ! spectral calculation
    call rad_spectral_calc()
    ! rte calculation
    call rad_rte_calc()
    ! source calculation
    call rad_spectral_srad()
end subroutine rad_calc

subroutine rad_write_src()
    use rad_fields, only: radiationOn, ReactionRates, nRR
    integer :: irr, m
    include 'usrnlst.inc'
    radiationOn = RAD_ON

    !    if (.not.radiationOn) return
    if (.not.radiationOn) then 
        ReactionRates = 0.0d0
!        write(*,*) "radiation is off"
    else 
!        write(*,*) "radiation is on"
    irr = config%nrr
!    write(*,*) "irr", irr, nrr 
    if (irr.eq.0) return
    if (irr.le.nRR) then
        ReactionRates(:,irr) = k_g(:,1) !G(:,1) !Srad(:,0) !Srad(:,0) G(:,1)
    endif

    !write(*,*)G(:,1)
!    write(*,*) "smax", smax
    do m=1, smax
        irr = m+config%nrr
        if (irr .le. nRR) then
            ReactionRates(:,irr) = k_s(:,m,1) !+k_g(:,1) !Srad(:,m)
        endif
    enddo
 
    do m=1,smax+1 ! +3 
       irr = m + config%nrr + smax 
        ReactionRates(:,irr) = Srad(:,m-1) 
    end do 

    do m=1,smax+1  ! +3 
       irr = m + config%nrr + 2*smax+1
        ReactionRates(:,irr) = E(:,m-1,1)
    end do

    do m=1,1       !+ 1 
       irr = m + config%nrr + 3*smax +2 
        ReactionRates(:,irr) = G(:,1) 
    end do     
    end if 
end subroutine

real(dp) function S_Rpg(ijk)
    implicit none
    integer ijk
    S_Rpg = 0.d0
end function

real(dp) function S_Rcg(ijk)
    use rad_fields, only: radiationOn
    implicit none
    integer ijk
    if (radiationOn) then
        S_Rcg = Srad(ijk,0)
    else
        S_Rcg = 0.d0
    endif
end function

real(dp) function S_Rps(ijk,m)
    implicit none
    integer ijk, m
    S_Rps = 0.d0
end function

real(dp) function S_Rcs(ijk, m)
    use rad_fields, only: radiationOn
    implicit none
    integer ijk, m
    if (radiationOn) then
        S_Rcs = Srad(ijk,m)
    else
        S_Rcs = 0.d0
    endif
end function

end module rad
