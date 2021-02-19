!!*****************************************************************
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
integer, parameter :: nSp = 2
! recognized gas compositions 
! Currently two are recognized CO2 and H2O

integer, parameter :: CO2 = 1, H2O = 2 
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
