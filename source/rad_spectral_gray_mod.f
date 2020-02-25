!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SPECTRAL GRAY                                !
!                                                                      !
!  Purpose:  to find absorption, emission and scattering co-efficient  !
!                                                                      !
!  Authors: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module rad_spectral_gray
use rad_param
use rad_config
use rad_fields
use rad_gas_species, only : gasInfoType, rad_gas_species_init, getGasInfo
use rad_solid_spectral, only : solidParInfoType 
implicit none
private

! ---- Interfaces ------
public :: rad_spectral_gray_init
public :: rad_spectral_gray_final
public :: rad_spectral_graygas_calc
public :: rad_spectral_graysld_calc
public :: rad_spectral_graygasconst_calc
public :: rad_spectral_graysldconst_calc

interface rad_spectral_gray_init
    module procedure init1
end interface

interface rad_spectral_gray_final
    module procedure finalize1
end interface

interface rad_spectral_graygas_calc
    module procedure calcgas
end interface

interface rad_spectral_graysld_calc
    module procedure calcsld
end interface

interface rad_spectral_graygasconst_calc
    module procedure constgascalc
end interface

interface rad_spectral_graysldconst_calc
    module procedure constsldcalc
end interface

! ---- Data members ----
character(len=strLen) :: modelName_
integer :: modelId_
! ---- Parameters ------
! database
! Planck mean absorption coefficient of CO2
integer, parameter :: kpLen=126
real(dp), dimension(kpLen) :: kpT, kpCO2, kpH2O
data kpT   / 300, 320, 340, 360, 380, 400,&
             420, 440, 460, 480, 500, 520,&
             540, 560, 580, 600, 620, 640,&
             660, 680, 700, 720, 740, 760,&
             780, 800, 820, 840, 860, 880,&
             900, 920, 940, 960, 980,1000,&
            1020,1040,1060,1080,1100,1120,&
            1140,1160,1180,1200,1220,1240,&
            1260,1280,1300,1320,1340,1360,&
            1380,1400,1420,1440,1460,1480,&
            1500,1520,1540,1560,1580,1600,&
            1620,1640,1660,1680,1700,1720,&
            1740,1760,1780,1800,1820,1840,&
            1860,1880,1900,1920,1940,1960,&
            1980,2000,2020,2040,2060,2080,&
            2100,2120,2140,2160,2180,2200,&
            2220,2240,2260,2280,2300,2320,&
            2340,2360,2380,2400,2420,2440,&
            2460,2480,2500,2520,2540,2560,&
            2580,2600,2620,2640,2660,2680,&
            2700,2720,2740,2760,2780,2800/

data kpCO2 /0.259375,0.250671,0.246764,0.247224,0.251313,0.258158,&
            0.266878,0.276654,0.286778,0.296672,0.305889,0.314100,&
            0.321089,0.326724,0.330951,0.333771,0.335228,0.335395,&
            0.334368,0.332253,0.329164,0.325215,0.320517,0.315180,&
            0.309303,0.302980,0.296298,0.289334,0.282157,0.274829,&
            0.267406,0.259934,0.252456,0.245008,0.237619,0.230316,&
            0.223120,0.216049,0.209119,0.202339,0.195721,0.189270,&
            0.182992,0.176890,0.170966,0.165221,0.159654,0.154264,&
            0.149050,0.144008,0.139135,0.134429,0.129886,0.125500,&
            0.121269,0.117188,0.113252,0.109457,0.105799,0.102272,&
            0.0988730,0.0956280,0.0924685,0.0894237,0.0864895,0.0836616,&
            0.0809361,0.0783092,0.0757772,0.0733364,0.0709833,0.0687146,&
            0.0665270,0.0644172,0.0623824,0.0604196,0.0585259,0.0566987,&
            0.0549354,0.0532336,0.0515907,0.0500045,0.0484729,0.0469937,&
            0.0455648,0.0441845,0.0428507,0.0415618,0.0403159,0.0391116,&
            0.0379472,0.0368212,0.0357322,0.0346787,0.0336596,0.0326735,&
            0.0317192,0.0307956,0.0299015,0.0290359,0.0281978,0.0273861,&
            0.0265999,0.0258384,0.0251006,0.0243858,0.0236931,0.0230217,&
            0.0223710,0.0217402,0.0211287,0.0205357,0.0199608,0.0194031,&
            0.0188623,0.0183377,0.0178288,0.0173350,0.0168560,0.0163911,&
            0.0159400,0.0155022,0.0150772,0.0146647,0.0142642,0.0138754/

data kpH2O /0.370897,0.335415,0.305788,0.280695,0.259118,0.240292,&
            0.223648,0.208767,0.195337,0.183127,0.171965,0.161717,&
            0.152277,0.143561,0.135499,0.128031,0.121106,0.114678,&
            0.108706,0.103152,0.0979827,0.0931668,0.0886759,0.0844838,&
            0.0805666,0.0769022,0.0734704,0.0702529,0.0672326,0.0643942,&
            0.0617236,0.0592079,0.0568353,0.0545951,0.0524776,0.0504738,&
            0.0485755,0.0467754,0.0450666,0.0434429,0.0418987,0.0404303,&
            0.0390298,0.0376943,0.0364198,0.0352026,0.0340393,0.0329268,&
            0.0318620,0.0308424,0.0298654,0.0289287,0.0280302,0.0271677,&
            0.0263396,0.0255439,0.0247792,0.0240437,0.0233363,0.0226554,&
            0.0219998,0.0213753,0.0207665,0.0201797,0.0196139,0.0190682,&
            0.0185418,0.0180337,0.0175432,0.0170695,0.0166120,0.0161699,&
            0.0157426,0.0153295,0.0149301,0.0145438,0.0141700,0.0138083,&
            0.0134582,0.0131192,0.0127910,0.0124731,0.0121651,0.0118666,&
            0.0115773,0.0112969,0.0110249,0.0107612,0.0105054,0.0102573,&
            0.0100164,0.00978269,0.00955579,0.00933549,0.00912156,0.00891377,&
            0.00871191,0.00851578,0.00832519,0.00813994,0.00795986,0.00778477,&
            0.00761451,0.00744891,0.00728783,0.00713111,0.00697861,0.00683020,&
            0.00668574,0.00654510,0.00640817,0.00627483,0.00614495,0.00601845,&
            0.00589519,0.00577510,0.00565807,0.00554400,0.00543280,0.00532439,&
            0.00521868,0.00511559,0.00501505,0.00491696,0.00482127,0.00472790/
contains

subroutine init1(config)
    use rad_fields, only : rad_fields_allocation
    implicit none
    type(configuration), intent(in) :: config
    nq=1
    call rad_fields_allocation(nq)
    wq(:) = 1.0
    gq(:) = 1.0
end subroutine

subroutine finalize1()
    use rad_fields, only : rad_fields_deallocation
    call rad_fields_deallocation
end subroutine

subroutine calcgas(gasinfo,kg, emiss)
    use rad_util, only : bbem
    implicit none
    type(gasInfoType), intent(in) :: gasinfo
    real(dp), dimension(1), intent(out) :: kg 
    real(dp), dimension(0:mphas,1), intent(out) :: emiss
    integer m
    kg(1) = grayGasKp(gasinfo)
    call bbem(emiss(0,1), gasInfo%T)
end subroutine

subroutine constgascalc(gasinfo,kg, emiss)
    use rad_util, only : bbem
    implicit none
    include 'usrnlst.inc'
    type(gasInfoType), intent(in) :: gasinfo
    real(dp), dimension(1), intent(out) :: kg
    real(dp), dimension(0:SMAX,1), intent(out) :: emiss
    integer m
    kg(1) = const_absg 
    call bbem(emiss(0,1), gasInfo%T)
end subroutine

subroutine calcsld(solidinfos, ks, scats, emiss)
    use rad_util, only : bbem
    USE discretelement, only : discrete_element
    implicit none
    type(solidParInfoType), dimension(:), intent(in) :: solidinfos
    real(dp), dimension(mphas,1), intent(out) :: ks, scats
    real(dp), dimension(0:mphas,1), intent(out) :: emiss
    integer m
    do m=1, size(solidinfos)
         call graySolidProperties(solidinfos(m), ks(m,1), scats(m,1), emiss(m,1))
    end do
end subroutine

subroutine constsldcalc(solidinfos, ks, scats, emiss)
    use rad_util, only : bbem
    USE discretelement, only : discrete_element
    implicit none 
    include 'usrnlst.inc'
    type(solidParInfoType), dimension(:), intent(in) :: solidinfos
    real(dp), dimension(SMAX,1), intent(out) :: ks, scats
    real(dp), dimension(0:SMAX,1), intent(out) :: emiss
    integer m
    if(discrete_element) then 
      do m=1, size(solidinfos)
        call graySolidProperties(solidinfos(m), ks(m,1), scats(m,1), emiss(m,1))
      end do
    else 
      do m=1, size(solidinfos)
         call grayConstSolidProperties(solidinfos(m), ks(m,1), scats(m,1), emiss(m,1))
      end do
    end if 
end subroutine

function grayGasKp(info) result (kp)
    use rad_util, only : interp1ess
    implicit none
    type(gasinfotype), intent(in) :: info
    real(dp) :: kp
    real(dp) :: kco2,kh2o
    kco2 = interp1ess(kpLen,kpT, kpCO2, info%T)
    kh2o = interp1ess(kpLen,kpT, kpH2O, info%T)
    !---M.S-rad: the tabulated absorption coefficients are in cm^-1 atm^-1 and P is in bar
    kp = (kco2*info%C(CO2)+kh2o*info%C(H2O)) * (info%P * 0.986923) * info%volFrac
end function

subroutine graySolidProperties(solidInfo, ks, scats, emiss)
    use rad_solid_spectral, only : genPlanckMeanParticleAbsc,&
                                   genGrayParticleScat
    use rad_util, only : bbem
    USE discretelement, only : discrete_element
    implicit none
    type(solidParInfoType), intent(in) :: solidinfo
    real(dp), intent(out) :: ks, scats, emiss
    ks = genPlanckMeanParticleAbsc(solidInfo)
    scats = genGrayParticleScat(solidInfo)
    call bbem(emiss, solidInfo%T)
end subroutine


subroutine grayConstSolidProperties(solidInfo, ks, scats, emiss)
    use rad_solid_spectral, only : genPlanckMeanParticleAbsc,&
                                   genGrayParticleScat
    use rad_util, only : bbem
    use fldvar, only: ep_s, d_p, rop_s, ro_s, T_s, X_s, theta_m
    implicit none
    include 'usrnlst.inc'
    type(solidParInfoType), intent(in) :: solidinfo
    real(dp), intent(out) :: ks, scats, emiss

    ks =  (3.0d0/4.0d0)*solidInfo%volFrac * const_emisp/solidInfo%radius
    emiss = (sigma/Pi)*(solidInfo%T)**4.d0
    if (UNITS=='CGS') emiss = emiss*1.d3  
    scats = (3.0d0/4.0d0)*(1.0-const_sfp)*(1.0-const_emisp)/solidInfo%radius
end subroutine

end module rad_spectral_gray
