!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SPECTRAL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authors: V Kotteda, M Stoellinger                   Date: 1-july-19  !
!  Authors: V Kotteda                                  Date: 18-feb-20  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module rad_spectral
use rad_param
use rad_config
use rad_fields
use rad_gas_species, only : gasInfoType, rad_gas_species_init, getGasInfo
use rad_solid_species, only : solidParInfoType, rad_solid_species_init, getSolidInfo, smax
use rad_solid_species, only : getDesSolidInfo
implicit none
private

! ---- Interfaces ------
public :: rad_spectral_init
public :: rad_spectral_calc
public :: rad_spectral_srad
public :: rad_spectral_final

interface rad_spectral_init
    module procedure init1
end interface

interface rad_spectral_final
    module procedure finalize1
end interface

interface rad_spectral_calc
    module procedure calculateSpectralFields
end interface

interface rad_spectral_srad
    module procedure calcRadSources
end interface
! ---- Data members ----
character(len=strLen) :: modelName_
integer :: modelId_

! ---- Parameters ------
! model ids
integer, parameter :: unknown=0, gray=1
integer, parameter :: graycon=2,graywsgg=3, nongraywsgg=4 

contains

subroutine init1(config)
    use rad_spectral_gray, only : rad_spectral_gray_init
    implicit none
    type(configuration), intent(in) :: config
    modelName_ = config%SpectralModelName
    ! parse options
    select case (modelName_)
        case ("GRAY")
            modelId_ = gray
            call rad_spectral_gray_init(config)
        case("GRAY-CONST")
            modelId_ = graycon
            call rad_spectral_gray_init(config)
        case("GRAY-WSGG")
            modelId_ = graywsgg
            call rad_spectral_gray_init(config)
        CASE("NONGRAY-WSGG")
            modelId_ = nongraywsgg
            call rad_spectral_gray_init(config)
        case default
            modelId_ = unknown
    end select
    call rad_gas_species_init(config)
    call rad_solid_species_init(config)
end subroutine

subroutine finalize1()
    use rad_spectral_gray, only : rad_spectral_gray_final
    implicit none
    select case (modelId_)
        case(gray)
            call rad_spectral_gray_final()
        case(graycon) 
            call rad_spectral_gray_final()
        case(graywsgg)
            call rad_spectral_gray_final()
        CASE(nongraywsgg)
            call rad_spectral_gray_final()
        case default
    end select
end subroutine

subroutine calculateSpectralFields()
    use functions, only : fluid_at,wall_at, funijk
    use discretelement
    use mfix_pic, only: des_stat_wt, mppic, epg_p
    use functions, only: IS_NORMAL
    USE compar 
    use des_thermo, only : DES_T_s
    use fldvar, only: ep_g 
    use geometry, only : vol, DX, DY, DZ 
    use run, only: UNITS
    use indices, only: i_of, j_of, k_of
    use functions, only: funijk
    USE compar
    implicit none
    include 'usrnlst.inc'
    type(gasInfoType) :: gasinfo
    type(solidParInfoType), dimension(mphas) :: solidinfos
    real(dp),dimension(nq) :: kg          ! gas absorption coefficient cm^-1
    real(dp),dimension(mphas,nq) :: ks    ! solid absorption coefficient cm^-1
    real(dp),dimension(mphas,nq) :: scats ! solid scattering coefficient cm^-1
    real(dp),dimension(0:mphas,nq) :: emiss !spectral emission

    integer :: ijk, m, i, j, k, NP, lM, iq
    real(dp) :: p_conversion,cellvol,areabvol,sbcon
    real(dp), dimension(nq) :: ai,ki
    real(dp) :: tmp, rmrk,spl,ttr


    !---M.S-rad: initialize the gasinfo 
    gasinfo%T = 0.
    gasinfo%P = 0.
    gasinfo%C(:) = 0.
    gasinfo%volFrac = 0.
    !---M.S-rad: initialize the solidsinfos
    solidinfos%volFrac = 0.
    solidinfos%radius = 0.
    solidinfos%T = 0.
    solidinfos%totSpecies = 0
    do m=1, mphas
        solidinfos(m)%speciesId(:) = 0
        solidinfos(m)%speciesMassFrac(:) = 0.
    end do
    !---M.S-rad: Determine the pressure conversion from the current MFIX unit to bar
    if (UNITS=='CGS') then
        p_conversion = 1.01325e-6  ! 1atm = 1.01325e+6 barye 
    else if (UNITS=='SI') then
        p_conversion = 1.01325e-5  ! 1atm = 1.01325e+5 Pa 
    else 
        print*, "Radiative property calculation: unknown units: ", UNITS
        print*, "Known units are CGS and SI"
    endif

    sbcon = sigma  ! SI system
    if(UNITS=='CGS') sbcon = sigma * 1.d3

   if(modelId_.le.2) then  
       do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
              call getGasInfo(ijk, gasinfo)
              !---M.S-rad: convert pressure to bar
              gasinfo%P = gasinfo%P * p_conversion
              ! find absorption coeffient and emission term 
              call spectralgasCalc(gasinfo, kg, emiss) 
              ! put cell calculations into fields
              k_g(ijk,:) = kg
              E(ijk,0,:) = emiss(0,:)
          endif
       enddo
    end if 


    if(discrete_element.or.mppic) then
      if(modelId_.eq.2) then 
      ! gray constant DEM 
        do np=1,max_pip
          if(is_normal(np)) then
           lM = pijk(np,5)
             ijk = pijk(np,4)
               if(fluid_at(ijk)) then
                   cellvol        = vol(ijk)    
                   areabvol       = Pi*(des_radius(np)**2.0)/cellvol
                   k_s(ijk,lm,:)  = k_s(ijk,lm,:) + const_emisp*areabvol
                   E(ijk,lM,:)    = (sbcon/Pi)*(DES_T_s(NP))**4.d0
                   scat(ijk,lM,:) = scat(ijk,lM,:)+ (1.0d0-const_sfp)* (1.0d0-const_emisp)*areabvol
               end if! fluid_at loop ends here
          endif      ! normal   loop ends here
        end do       ! particle loop ends here
      else 
     ! gray  and non-gray DEM/PIC models 
         do ijk = ijkstart3, ijkend3
          if(fluid_at(ijk)) then
              if(pinc(ijk) > 0) then
                   do m=1, mphas
                      call getDesSolidInfo(ijk,m, solidinfos(m))
                   end do
                   call spectralsldCalc(solidinfos, ks, scats, emiss)
                   do m=1,mphas
                      k_s(ijk,m,:)  = ks(m,:)
                      E(ijk,m,:)    = emiss(m,:)
                      if(modelId_.eq.4) E(ijk,m,:) = E(ijk,m,:)*1.0/float(nq) ! Non-gray model 
                      scat(ijk,m,:) = scats(m,:)
                   end do
              end if
          end if
         end do
      end if 
   ! solid phase in TFM (gray and non-gray models) 
    else 
      do ijk = 1, dimension_3
        if (fluid_at(ijk)) then
            call getGasInfo(ijk, gasinfo)
            !---M.S-rad: convert pressure to bar
            gasinfo%P = gasinfo%P * p_conversion
              do m=1, smax
                 call getSolidInfo(ijk,m, solidinfos(m))
              end do
              call spectralsldCalc(solidinfos, ks, scats, emiss)
            ! put cell calculations into fields
             do m=1, smax
                   E(ijk,m,:)    = emiss(m,:) 
                   if(modelId_.eq.4) E(ijk,m,:) = E(ijk,m,:)*1.0/float(nq) ! Non-gray model 
                   k_s(ijk,m,:)  = ks(m,:)
                   scat(ijk,m,:) = scats(m,:)
             end do
        endif
      enddo
    end if 

    ! unit conversions
    if (UNITS=='CGS') then
        scat = scat/1.d2
    endif

    call extrapolation()

end subroutine

subroutine calcRadSources()
    use discretelement, only: discrete_element
    use mfix_pic,       only: mppic
    use fldvar,         only: T_g
    use functions,      only: funijk
    use geometry,       only: dx,dy,dz
    USE compar

    implicit none
    type(gasInfoType) :: gasinfo
    real(dp),dimension(1) :: T
    integer iq, m,ijk,i,j,k
    double precision :: sbcon
    include 'usrnlst.inc'
    radiationOn = RAD_ON 
    Srad = 0.d0

    sbcon = sigma  ! SI system
    if(UNITS=='CGS') sbcon = sigma * 1.d3

   if(modelId_.eq.3) then
     do ijk = 1, dimension_3
       do iq = 1, nq
          Srad(ijk,0) = Srad(ijk,0)+k_g(ijk,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,0,iq))
          do m=1, mphas
             Srad(ijk,m) = Srad(ijk,m)+k_s(ijk,m,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,m,iq))
          end do
       end do
      end do 
   else if(modelId_.eq.4) then 
     do ijk = 1, dimension_3
       do iq = 1, nq
          Srad(ijk,0) = Srad(ijk,0)+ k_g(ijk,iq)*(G(ijk,iq)-a_g(ijk,iq)*4.d0*sbcon*T_g(ijk)**4.0)  
          do m=1, mphas
             Srad(ijk,m) = Srad(ijk,m)+k_s(ijk,m,iq)*(G(ijk,iq)-4.d0*pi*E(ijk,m,iq))
          end do
       end do
      end do
   else 
       do iq = 1, nq
        Srad(:,0) = Srad(:,0)+k_g(:,iq)*(G(:,iq)-4.d0*pi*E(:,0,iq)) 
    end do
   end if 

   if(discrete_element.or.mppic) then
       do iq = 1, nq
           do m=1, mphas
               Srad(:,m) = Srad(:,m)+k_s(:,m,iq)*(G(:,iq)-4.d0*pi*E(:,m,iq))
           end do
       end do 
   else
       do iq = 1,nq   
           do m=1, smax
               Srad(:,m) = Srad(:,m)+k_s(:,m,iq)*(G(:,iq)-4.d0*pi*E(:,m,iq))
           end do
       end do
   end if 

      ! by far Srad units are W/(m^3) in SI or erg/(cm^3.s) in CGS
      ! C_pg uses J/kg.K in SI or cal/g.K in CGS
      ! 1 erg = 2.39005736d-8 cal
      if (UNITS=='CGS') Srad = Srad * 2.39005736d-8
    
end subroutine


subroutine spectralsldCalc(solidinfos, ks, scats, emiss)
    use rad_spectral_gray, only :rad_spectral_graysld_calc,rad_spectral_graysldconst_calc
    implicit none
    type(solidParInfoType),dimension(mphas), intent(in) :: solidinfos
    real(dp),dimension(mphas,nq), intent(out) :: ks
        ! solid absorption coefficient cm^-1
    real(dp),dimension(mphas,nq), intent(out) :: scats
        ! solid scattering coefficient cm^-1
    real(dp),dimension(0:mphas,nq), intent(out) :: emiss
        !spectral emission 

    select case(modelId_)
        case (gray)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case(graycon)
            call rad_spectral_graysldconst_calc(solidinfos, ks, scats, emiss)
        case(graywsgg)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case(nongraywsgg)
            call rad_spectral_graysld_calc(solidinfos, ks, scats, emiss)
        case default
    end select
end subroutine

subroutine spectralgasCalc(gasinfo, kg, emiss)
    use rad_spectral_gray, only :rad_spectral_graygas_calc,rad_spectral_graygasconst_calc
    implicit none
    type(gasInfoType), intent(in) :: gasinfo
    real(dp),dimension(nq), intent(out) :: kg
        ! gas absorption coefficient cm^-1
        ! gas emission coefficient cm^-1
    real(dp),dimension(0:mphas,nq), intent(out) :: emiss
        !spectral emission 

    select case(modelId_)
        case (gray)
            call rad_spectral_graygas_calc(gasinfo, kg, emiss)
        case(graycon)
            call rad_spectral_graygasconst_calc(gasinfo, kg, emiss)
        case default
    end select
end subroutine


subroutine extrapolation()
    use rad_util, only : extrapolateGhostCells
    implicit none
    integer :: iq, is
    do iq=1, nq
        call extrapolateGhostCells(k_g(:,iq))
        do is = 1, size(k_s, 3)
            call extrapolateGhostCells(k_s(:,is,iq))
            call extrapolateGhostCells(scat(:,is,iq))
        end do
    end do

end subroutine

end module rad_spectral
