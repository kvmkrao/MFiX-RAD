!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SPECTRAL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 1-july-19  !
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

! smith 1982    http://doi.org/10.1115/1.3245174
! P_H20/p_CO2 = 1 
real(dp), dimension(3) :: ki1
real(dp), dimension(3,4) :: bei1
! P_H20/p_CO2 = 2
real(dp), dimension(3) :: ki2
real(dp), dimension(3,4) :: bei2

! BORDBER2014    https://doi.org/10.1016/j.combustflame.2014.03.013  
real(dp), dimension(4,5,5) :: cijk
real(dp), dimension(4,5) :: dik
! JOHANSSON2011  https://doi.org/10.1016/j.combustflame.2011.02.001
real(dp), dimension(4,3) :: c1ik,c2ik,c3ik
real(dp), dimension(4)   :: k1j, k2j
! Dorigon2013    https://doi.org/10.1016/j.ijheatmasstransfer.2013.05.010 
real(dp), dimension(4) :: kpi1, kpi2
real(dp), dimension(4) :: b1i0, b1i1, b1i2,b1i3,b1i4 
real(dp), dimension(4) :: b2i0, b2i1, b2i2,b2i3,b2i4 
! TAYLOR74       https://doi.org/10.1016/0017-9310(74)90067-2 
real(dp), dimension(4)  :: kf2i, bf21i, bf22i 
real(dp), dimension(4)  :: kf1i, bf11i, bf12i

! smith1982 
data ki1  / 0.4303, 7.055, 178.1 /
data bei1 /5.150,  0.7749, 1.907, &
          -2.303,  3.399, -1.824, &
           0.9779,-2.297,  0.5608, &
          -1.494,  3.770, -0.5122 /
data ki2 / 0.4201, 6.516, 131.9 /
data bei2 / 6.508, -0.2504, 2.718, &
            -5.551, 6.112, -3.118, &
             3.029,-3.882, 1.221, &
            -5.353, 6.528,-1.612 /

!  TAYLOR74  
data kf1i  / 0.0, 0.91, 9.4, 130.0/
data kf2i  / 0.0, 0.69, 7.4, 80.0 /
data bf21i / 0.364, 0.266, 0.252, 0.118 /
data bf22i / 4.74, 7.19, -7.41, -4.52 /

data bf11i/ 0.4092, 0.284, 0.211,  0.0958/ 
data bf12i / 7.53,    2.58, -6.54, - 3.57/

! JOHANSSON2011 
data c1ik / 0.358 , 0.0731,-0.0466, &
            0.392, -0.212,  0.0191, &
            0.142, -0.0831, 0.0148, &
            0.0798,-0.0370, 0.0023/
data c2ik / -0.165,-0.0554, 0.0930, &
            -0.291, 0.644, -0.209, &
            0.348,  -0.294, 0.0662, &
            0.0866, -0.106,  0.0305 /
data c3ik / 0.0598,  0.0028,-0.0256, &
            0.0784, -0.197,  0.0662, &
            -0.122,  0.118, -0.0295, &
            -0.0127, 0.0169,-0.0051 /

data k1j / 0.055, 0.88, 10.0, 135 /
data k2j / 0.012,-0.021,-1.6, -35 /

! Dorigon2013
data kpi2 / 0.1921, 1.719, 11.37, 111.016/
data kpi1 / 0.01873, 1.723, 12.48,  144.9/

data b1i0 / 0.07197 , 0.11070, 0.20910, 0.07092/
data b1i1 / 87.24, 33.97, -6.423, 6.586/
data b1i2 / -96.90, -24.67, -3.200, -12.78/ 
data b1i3 / 46.51, 4.647, 1.718, 5.577/ 
data b1i4 / -79.17, -1.039, -2.105, -7.709/ 

data b2i0 / 0.05617, 0.14260, 0.13620, 0.12220/
data b2i1 / 78.44  , 17.95  , 25.74  , -2.327/
data b2i2 / -85.63 , -1.077 , -37.11 , -7.492/
data b2i3 / 42.46  , -6.971 , 15.75  , 4.275 /
data b2i4 / -74.4  , 17.74  , -22.67 , -6.608/

! Borddar et al. (2014), A line by line based weighted sum of gray gases model for inhomogeneous
! co2-h2o mixture in oxy-fired combustion 
data dik  / 0.0340429, 0.3509457, 4.5707400, 109.81690, &
            0.0652305, 0.7465138, 2.1680670,-50.923590, &
           -0.0463685,-0.5293090,-1.4989010, 23.432360, &
           0.0138684, 0.1594423,  0.4917165, -5.1638920, &
          -0.0014450,-0.0166326, -0.0542999, 0.4393889/

data cijk /0.7412956, 0.1552073, 0.2550242,-0.0345199, &
          -0.9412652, 0.6755648,-0.6065428, 0.4112046, &
           0.8531866,-1.1253940, 0.8123855,-0.5055995, &
          -0.3342806, 0.6040543,-0.4532290, 0.2317509, &
           0.0431436,-0.1105453, 0.0869309,-0.0375491, &
          -0.5244441,-0.4862117, 0.3805403, 0.2656726, &
           0.2799577, 1.4092710, 0.3494024,-0.5728350, &
           0.0823075,-0.5913199,-1.1020091, 0.4579559, &
           0.1474987,-0.0553385, 0.6784475,-0.1656759, &
          -0.0688622, 0.0464663,-0.1306996, 0.0229520, &
           0.5822860, 0.3668088,-0.4249709,-0.1225365, &
          -0.7672319,-1.3834490, 0.1853509, 0.2924490, &
           0.5289430, 0.9085441, 0.4046178,-0.2616436, &
          -0.4160689,-0.1733014,-0.3432603, 0.1052608, &
           0.1109773,-0.0016129, 0.0741446,-0.0160047, &
          -0.2096994,-0.1055508, 0.1429446, 0.0300151, &
           0.3204027, 0.4575210,-0.1013694,-0.0798076, &
          -0.2468463,-0.3334201,-0.0811822, 0.0764841, &
           0.1697627, 0.0791608, 0.0883088,-0.0321935, &
          -0.0420861,-0.0035398,-0.0202929, 0.0050463, &
           0.0242031, 0.0105857,-0.0157408,-0.0028205, &
          -0.0391017,-0.0501976, 0.0130244, 0.0079966, &
           0.0310940, 0.0384236, 0.0062981,-0.0079084, &
          -0.0204066,-0.0098934,-0.0084152, 0.0033870, &
           0.0049188, 0.0006121, 0.0020110,-0.0005364 /

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

    ! gas phase   (gray medium, WSGG model)  
     if(modelId_.eq.3) then
 
      do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
            call getGasInfo(ijk, gasinfo)
            !---M.S-rad: convert pressure to bar
            gasinfo%P = gasinfo%P * p_conversion
            rmrk       = gasinfo%C(H2O)/gasinfo%C(CO2)
            if(gasinfo%C(H2O).eq.0) rmrk = 0.0d0 
            spl        = rad_spl
            if(UNITS=='CGS') spl = spl*1.e-2 ! cm to m  
            tmp        = 0.0d0 
            ai         = 0.0d0
            ki         = 0.0d0 
            ttr        = gasinfo%T
            if(rmrk.le.0.0d0.or.gasinfo%C(CO2).eq.0.0d0) then 
               tmp = 0.0d0 
            else if(rad_wsgg.eq.'SMITH82') then
                if (rmrk.eq.2) then
                  do i=1,3
                    ki(i) = ki2(i)
                    ai(i) = bei2(i,1)*1.0D-1 + bei2(i,2)*1.0D-4*ttr+ &
                    bei2(i,3)*1.0D-7*ttr**2.0 +bei2(i,4)*1.0D-11*ttr**3.0
                    tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  end do
                else 
                  do i=1,3
                    ki(i) = ki1(i)
                    ai(i) = bei1(i,1)*1.0D-1 + bei1(i,2)*1.0D-4*ttr+ &
                    bei1(i,3)*1.0D-7*ttr**2.0 +bei1(i,4)*1.0D-11*ttr**3.0
                    tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  end do
                end if
            else if(rad_wsgg.eq.'BORDBAR14') then
                do i=1,4
                   ai(i)  = 0.0d0 
                   ki(i) = 0.0d0
                   do k=1,5   ! absorption co-efficient 
                      ki(i) = ki(i) + dik(i,k) * (rmrk)**(k-1)
                      do j=1,5 ! weighting factors 
                        ai(i)= ai(i)+cijk(i,j,k)*(rmrk**(k-1))*(ttr/1200.0)**(j-1)
                      end do
                   end do
                tmp= tmp + ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl)) 
                end do
            else if(rad_wsgg.eq.'JOHANSSON11') then
                do i=1,4
                  ai(i) = 0.0d0
                  ki(i) = k1j(i) +k2j(i) *rmrk
                  do k=1,3
                    ai(i)=ai(i)+(c1ik(i,k)+c2ik(i,k)*rmrk+c3ik(i,k)*rmrk*rmrk)*(ttr/1200.0)**(k-1)
                  end do
                  tmp= tmp+ai(i)*(1.0d0-exp(-ki(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
               end do
            else if (rad_wsgg.eq.'TAYLOR74') then  
              if(rmrk.eq.2) then 
                 do i=1,4
                    ai(i) = bf21i(i) + bf22i(i) * ttr* 1.0d-5
                    tmp   = tmp + ai(i)*(1.0-exp(-kf2i(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                    ki(i) = kf2i(i)
                 end do
              else 
                 do i=1,4
                    ai(i) = bf11i(i) + bf12i(i) * ttr* 1.0D-5
                    tmp   = tmp + ai(i)*(1.0d0-exp(-kf1i(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                    ki(i) = kf1i(i) 
                 end do
              end if

            else if(rad_wsgg.eq.'DORIGON13') then 
              if (rmrk.eq.2) then  ! 400 K to 2500 K , 0.001 atmm to 10 atm m 
                do i=1,4 
                  ai(i) = b2i0(i) + b2i1(i)*1.0e-5*ttr  + b2i2(i)*1.0e-8*ttr*ttr + b2i3(i)*1.0e-11*ttr**3.0 + b2i4(i) * 1.0E-15*ttr**4.0
                  tmp = tmp + ai(i)*(1.0-exp(-kpi2(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  ki(i) = kpi2(i)
               end do 
              else   ! rmrk.eq.1 
               do i=1,4
                  ai(i) = b1i0(i) + b1i1(i)*1.0e-5*ttr  + b1i2(i)*1.0e-8*ttr*ttr + b1i3(i)*1.0e-11*ttr**3.0 + b1i4(i) * 1.0E-15*ttr**4.0
                  tmp = tmp + ai(i)*(1.0-exp(-kpi1(i)*(gasinfo%P * 0.986923)*(gasinfo%C(H2O)+gasinfo%C(CO2))*spl))
                  ki(i) = kpi1(i)
               end do 
              end if
            end if  ! rad_wsgg loop end
            kg(1) =  -log(1.0d0 - tmp) / spl
            k_g(ijk,:) = kg(1)   
            E(ijk,0,:) = sbcon*(ttr**4.0)/Pi  
            if(UNITS=='CGS') k_g(ijk,:) = k_g(ijk,:) *1.0D-2 ! meters to centi meters
      end if     ! fluid_ijk loop  
     end do      ! do loop 
     else if (modelId_.eq.4) then
      do ijk = 1, dimension_3
          if (fluid_at(ijk)) then
          call getGasInfo(ijk, gasinfo)
          !---M.S-rad: convert pressure to bar
          gasinfo%P = gasinfo%P * p_conversion
          rmrk       = gasinfo%C(H2O)/gasinfo%C(CO2)
          if(gasinfo%C(H2O).eq.0) rmrk = 0.0d0 
          spl        = rad_spl
          if(UNITS=='CGS') spl = spl*1.e-2 ! cm to m 
          ttr        = gasinfo%T
            if(gasinfo%C(H2O).eq.0.0.or.gasinfo%C(CO2).eq.0) goto 999
            if(rad_wsgg.eq.'SMITH82') then
                if (rmrk.eq.2) then
                  do i=1,nq  
                    ki(i) = ki2(i)
                    ai(i) = bei2(i,1)*1.0D-1 + bei2(i,2)*1.0D-4*ttr+ &
                    bei2(i,3)*1.0D-7*ttr**2.0 +bei2(i,4)*1.0D-11*ttr**3.0
                  end do
                else 
                  do i=1,nq  
                    ki(i) = ki1(i)
                    ai(i) = bei1(i,1)*1.0D-1 + bei1(i,2)*1.0D-4*ttr+ &
                    bei1(i,3)*1.0D-7*ttr**2.0 +bei1(i,4)*1.0D-11*ttr**3.0
                  end do
                end if
            else if(rad_wsgg.eq.'BORDBAR14') then
                do i=1,nq 
                   ki(i) = 0.0d0 
                   ai(i) = 0.0d0
                   do k=1,5   ! absorption co-efficient 
                      ki(i) = ki(i) + dik(i,k) * (rmrk)**(k-1)
                      do j=1,5 ! weighting factors 
                        ai(i)= ai(i)+cijk(i,j,k)*(rmrk**(k-1))*(ttr/1200.0)**(j-1)
                      end do
                   end do
                end do
            else if(rad_wsgg.eq.'JOHANSSON11') then
                do i=1,nq 
                  ai(i) = 0.0d0 
                  ki(i) = k1j(i) +k2j(i) *rmrk
                  do k=1,3
                    ai(i)=ai(i)+(c1ik(i,k)+c2ik(i,k)*rmrk+c3ik(i,k)*rmrk*rmrk)*(ttr/1200.0)**(k-1)
                  end do
              end do
            else if (rad_wsgg.eq.'TAYLOR74') then  
              if(rmrk.eq.2) then 
                 do i=1,nq
                    ai(i) = bf21i(i) + bf22i(i) * ttr* 1.0d-5
                    ki(i) = kf2i(i)
                 end do
              else
                 do i=1,nq
                    ai(i) = bf11i(i) + bf12i(i) * ttr* 1.0d-5
                    ki(i) = kf1i(i)
                 end do
              end if  
            else if(rad_wsgg.eq.'DORIGON13') then 
              if (rmrk.eq.2) then  ! 400 K to 2500 K , 0.001 atmm to 10 atm m 
               do i=1,nq  
                 ai(i) = b2i0(i) + b2i1(i)*1.0e-5*ttr  + b2i2(i)*1.0e-8*ttr*ttr + b2i3(i)*1.0e-11*ttr**3.0 + b2i4(i) * 1.0E-15*ttr**4.0
                 ki(i) = kpi2(i)
               end do 
              else                 ! rmrk.eq.1 
               do i=1,nq 
                  ai(i) = b1i0(i) + b1i1(i)*1.0e-5*ttr  + b1i2(i)*1.0e-8*ttr*ttr + b1i3(i)*1.0e-11*ttr**3.0 + b1i4(i) * 1.0E-15*ttr**4.0
                  ki(i) = kpi1(i)
               end do 
              end if
            end if  ! rad_wsgg loop end
999	   continue 
           do i=1,nq 
              k_g(ijk,i) = ki(i) 
              a_g(ijk,i) = ai(i)
              E(ijk,0,i) = ai(i)*sbcon*(ttr**4.0)/Pi  
              if(UNITS=='CGS') k_g(ijk,i) = k_g(ijk,i) *1.0d-2 ! 1/meters to 1/cm 
          end do ! nq loop end
      end if     ! fluid_ijk loop  end  
     end do      ! do loop ijk end 

     else 
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
    use indices,        only: i_of, j_of, k_of
    use functions,      only: fluid_at, wall_at
    use functions,      only: funijk
    use geometry,       only: dx,dy,dz
    USE compar

    implicit none
    type(gasInfoType) :: gasinfo
    real(dp),dimension(1) :: T
    integer iq, m,ijk,i,j,k
    double precision :: sbcon
    real(dp), dimension(dimension_3)  :: tdis
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
