!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION PARAMETERS                                   !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authors: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!!** module rad_param
!!*****************************************************************
module rad_param
use param, only : dimension_bc, dimension_3
use physprop, only : MW_G, MW_MIX_g
use physprop, only : smax, nmax
use run, only : units
implicit none
! Parameters used to control radiation calculations.
! Occationally they may need to be adjusted.
! They could be converted to keywords if necessary

!-------------- Data type parameters --------------
integer, parameter :: dp = kind(1.d0)
! real number precision

integer, parameter :: strLen = 128
! default length of a string

!-------------- Thermochemical parameters ---------
integer, parameter :: nSp = 49
! recognized gas compositions 
! Currently two are recognized CO2 and H2O

integer, parameter :: CO2= 1, H2O= 2,  O3= 3,  N2O= 4,     CO= 5,    CH4= 6, & 
O2= 7,   NO= 8,   SO2= 9,   NO2= 10,   NH3= 11,  HNO3= 12, OH= 13,   HF= 14, & 
HCl= 15, HBr= 16, HI= 17,   ClO= 18,   OCS= 19,  H2CO= 20, HOCl= 21, & 
N2= 22,  HCN= 23, CH3Cl= 24,H2O2= 25,  C2H2= 26, C2H6= 27, PH3= 28, & 
COF2= 29,SF6= 30, H2S= 31,  HCOOH= 32, HO2= 33,  O= 34,    ClONO20= 35, & 
NOp= 36,  HOBr= 37,C2H4= 38, CH3OH= 39, CH3Br= 40,CH3CN= 41,CF4= 42, & 
C4H2= 43,HC3N= 44,H2= 45,   CS= 46,    SO3= 47,  C2N2= 48, COCl2= 49 


!!DATA (MOLID(I),I=0,NMOL)/ '   All', &
!!!           1        2        3        4        5        6        7
!!     '   H2O','   CO2','    O3','   N2O','    CO','   CH4','    O2', &
!!!           8        9       10       11       12       13       14
!!     '    NO','   SO2','   NO2','   NH3','  HNO3','    OH','    HF', &
!!!          15       16       17       18       19       20       21
!!     '   HCl','   HBr','    HI','   ClO','   OCS','  H2CO','  HOCl', &
!!!          22       23       24       25       26       27       28
!!     '    N2','   HCN',' CH3Cl','  H2O2','  C2H2','  C2H6','   PH3', &
!!!          29       30       31       32       33       34       35
!!     '  COF2','   SF6','   H2S',' HCOOH','   HO2','     O','ClONO20',&
!!!          36       37       38       39       40       41       42
!!     '   NO+','  HOBr','  C2H4',' CH3OH',' CH3Br',' CH3CN','   CF4',&
!!!          43       44       45       46       47       48       49
!!     '  C4H2','  HC3N','    H2','    CS','   SO3','  C2N2',' COCl2',&
!!!          50       51       52       53       54       55       56

! gas species id in gas mole fraction vector

integer, parameter :: GAS =0
! gas phase id in arrays for multiphase mixture

integer, parameter :: MaxSolidSpecies=10
! solid phase species array length

!-------------- Radiation spectral parameters -----
real(dp), parameter :: BETAMIN = 1.d-2
! minimum beta (extinction coefficient) for radiation calculation
! If the maximum extinction coefficient is below this number,
! the problem is considered as optically thin. Then absorption is
! negligible, and radiative heat source comes from emission only.

! real(dp), parameter :: BETALB = 1.d-5
! ! the lower bound value of extinction coefficient. The extinction
! ! coefficient in the flow field will be lower bounded by this
! ! number. This is because 1/Beta is used in P1 calculation and
! ! a small number of beta causes numerical instability

logical, parameter :: isSootGray=.true.  
! logical constant for gray soot in spectral calculation

!--------------- Physical constants --------------
real(dp), parameter :: Pi = 3.14159265359D0
real(dp), parameter :: sigma = 5.6704d-8 
! Stefan Boltzmann constant in SI unit W.m^(-2).K^(-4) SI unit
end module
