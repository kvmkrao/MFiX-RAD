!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: RADIATION SOLID SPECTRAL                               !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Authos: V Kotteda, M Stoellinger                   Date: 19-Apr-19  !
!   Comments                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

module rad_solid_spectral
! (1) performs all calculation between complex refractive indices to narrowband constant absorption coefficient
! (2) will also have gray scattering coefficient
use rad_param, only : dp, strLen, MaxSolidSpecies, Pi, sigma
use rad_solid_species, only : solidParInfoType
implicit none
private ! by default private for all
integer :: totalRadPartSpecies ! total involved particle species
integer, parameter :: NNU = 248 ! 248 narrowbands from gas module
integer, dimension(:), allocatable :: radPartIdList
real(dp), dimension(:,:), allocatable :: abscParticle ! particle absorption coefficients over narrowbands
!---------------------------------------------------------------------------
! database narrowband kappaStar
type kappaStarDBInfo
  integer :: SpeciesId
  real(dp) :: radius
  real(dp) :: kappaStar(NNU)
end type kappaStarDBInfo
integer, parameter :: maxkappaStarDBSize = 10
integer :: kappaStarDBsize=0
type(kappaStarDBInfo), dimension(maxkappaStarDBSize) :: kappaStarDB

public :: solidParInfoType, & ! type for Solid particle information
          initSolidSpectralCal, & ! initialize calculation
          particleAbscStar,     & ! m, a, eta -> kappaStar
          particleAbscStarLBL2NB, & ! project LBL kappaStar to narrowband
          genNBparticleAbscStar, & ! repeat project LBL kappaStar to narrowband for all narrowbands
          genNBparticleAbsc, & ! get narrow band linear based absorption coefficient
          genLBLparticleAbsc, & ! LBL linear based absorption coefficient
	  genGrayParticleScat, & ! generatre gray particle scattering coefficient
	  genPlanckMeanParticleAbsc ! generate Planck mean absorption coefficient from LBL calculation
private :: BBFN
contains

function initSolidSpectralCal(totNumSpecies, idList) result (statInfo)
! initialize this module
! not finished
integer, intent(in) :: totNumSpecies
integer, intent(in), dimension(:) :: idList
integer :: statInfo
integer :: totalRadPartSpecies
totalRadPartSpecies = totNumSpecies
statInfo=0
allocate(abscParticle(NNU, totalRadPartSpecies), radPartIdList(totalRadPartSpecies), stat=statInfo)
end function initSolidSpectralCal

function particleAbscStar(m, a, eta) result (kappaStar)
! calculate volume fraction particle absorption coefficient based on (1) solid complex refractive index
! (2) particle size and (3) wavenumber
! algorithm described by Modest 2003
complex(dp), intent(in) :: m ! complex refractive index
real(dp), intent(in) :: a, eta ! particle radius (unit m), wavenumber (unit cm^-1)
real(dp) :: kappaStar ! volume fraction based particle absorption coefficient
! the final absorption coefficient will be kappaStar times volume fraction (of species)
real(dp) :: x, betaStar0 ! size parameter, absorption/extinction coefficient at small particle limit
real(dp), parameter :: pi2 = 6.283185307d0 ! 2*pi
real(dp), parameter :: C1 = -1d0/1.6d0
complex(dp) :: f
real(dp) :: t ! temp variable
x= pi2*a*eta*100.d0
f = (m*m-1.d0)/(m*m+2.d0)
betaStar0 = -4.d0* x* imag(f)
t = (betaStar0+2.3d0*betaStar0**3)**(-1.6d0) 
t = t + (1.66d0/(betaStar0**0.16d0))**(-1.6d0)
kappaStar = t**C1
end function particleAbscStar
function particleExtcStar(m, a, eta) result ( betaStar)
! calculate volume fraction particle extinction coefficient based on (1) solid complex refractive index
! (2) particle size and (3) wavenumber
! algorithm described by Modest 2003
complex(dp), intent(in) :: m ! complex refractive index
real(dp), intent(in) :: a, eta ! particle radius (unit m), wavenumber (unit cm^-1)
real(dp) :: betaStar ! volume fraction based particle absorption coefficient
! the final absorption coefficient will be kappaStar times volume fraction (of species)
real(dp) :: x, betaStar0 ! size parameter, absorption/extinction coefficient at small particle limit
real(dp), parameter :: pi2 = 6.283185307d0 ! 2*pi
real(dp), parameter :: C1 = -1d0/1.2d0
complex(dp) :: f
real(dp) :: t ! temp variable
x= pi2*a*eta*100.d0
f = (m*m-1.d0)/(m*m+2.d0)
betaStar0 = -4.d0* x* imag(f)
t = (betaStar0+6.78d0*betaStar0**3)**(-1.2d0) 
t = t + (3.09d0/(betaStar0**0.10d0))**(-1.2d0)
betaStar = t**C1
end function particleExtcStar

function particleAbscStarLBL2NB(sid, a, eta_min, eta_max) result(kappaStarNB)
! project LBL particle absorption coefficient into constant narrowband value
! sophiscated numerical integration shall be used but currently use trap
use rad_solid_prop, only : id2cmplxRefInd
integer, intent(in) :: sid ! solid species id in database
real(dp), intent(in) :: a, eta_min, eta_max ! particle size, minimum and maximum wavenumber of a narrow band
! a (m), eta (cm^-1)
real(dp) :: kappaStarNB ! volume fraction based particle absorption coefficient over narrowband
integer, parameter :: N =2 ! sample finite number of points for a mean value
  ! if N=2, sample two bounds
real(dp), dimension(N) :: eta, kappaStar
real(dp) :: deltaEta
integer :: i
complex(dp) :: m ! complex refractive index
deltaEta = (eta_max - eta_min) /real(N-1)
eta = (/ (eta_min+ i*deltaEta, i=0, N-1)/)
do i = 1, N
  m = id2cmplxRefInd(sid, eta(i))
  kappaStar(i) = particleAbscStar(m, a, eta(i))
end do
kappaStarNB = (sum(kappaStar, 1)-0.5d0*(kappaStar(1)+kappaStar(N)))/(N-1)
  ! trapizoidal rule for integration
end function particleAbscStarLBL2NB

function genNBparticleAbscStar (sid, a) result(kappaStarNBList)
! given a species id and particle size generate narrow band constant kappaStar
use rad_solid_prop, only : genNBInterval
integer, intent(in) :: sid ! solid species id in database
real(dp), intent(in) :: a ! particle size (m)
real(dp), dimension(NNU) :: kappaStarNBList
real(dp), dimension(NNU,2) :: nbInt ! narrow band intervals
integer :: n, iDB
! check if databased
do iDB = 1, kappaStarDBSize
  if ((kappaStarDB(iDB)%SpeciesId.eq.sid) .and. abs(a-kappaStarDB(iDB)%radius)<1.d-10) then
    ! load databased value
    kappaStarNBList = kappaStarDB(iDB)%kappaStar
    !print*, 'Read from kappaStarNBdb'
    return
  endif
enddo
! get narrowband intervals from gas module
kappaStarNBList =0.d0
nbInt = genNBInterval()
do n = 1, NNU
  kappaStarNBList(n) = particleAbscStarLBL2NB(sid, a, nbInt(n,1), nbInt(n,2))
end do !
! add new kappaStar to database
if (kappaStarDBsize<maxKappaStarDBSize) then
  kappaStarDBsize= kappaStarDBsize+1
  kappaStarDB(kappaStarDBSize)%radius = a
  kappaStarDB(kappaStarDBsize)%SpeciesId=sid
  kappaStarDB(kappaStarDBsize)%kappaStar=kappaSTarNBList
  print*, "Add to KappaStarNBDB", kappaStarDBSize, sid, a
endif
end function genNBparticleAbscStar

function genNBparticleAbsc(particle_info) result(kappaNBList)
! unit m^-1
type(solidParInfoType), intent(in) :: particle_info
real(dp), dimension(NNU) :: kappaNBList, kappaStarNBList
real(dp) :: fv, a, fA
integer :: sid, s
kappaNBList=0.d0
do s = 1, particle_info%totSpecies
  sid = particle_info%speciesId(s)
  if (sid.le.0) cycle
  a = particle_info%radius
  fv = particle_info%volFrac
  fA = 0.75d0*fv/a
  kappaStarNBList = genNBparticleAbscStar(sid, a)
  kappaNBList = kappaNBList + fA*kappaStarNBList*particle_info%SpeciesMassFrac(s) 
  ! weighted by mass fraction
end do
end function genNBparticleAbsc

function genLBLparticleAbsc(particle_info, eta) result(kappa)
! unit m^-1
use rad_solid_prop, only : id2cmplxRefInd
type(solidParInfoType), intent(in) :: particle_info
real(dp), dimension(:) :: eta  ! list of wavenumbers
real(dp) :: kappa(size(eta,1)), kappaStar(size(eta,1))
real(dp) :: a, fv, fA
complex(dp) :: m
integer :: i, s, sid

a = particle_info%radius
fv = particle_info%volFrac
fA = 0.75d0*fv/a
kappa=0.d0
do s = 1, particle_info%totSpecies
  sid = particle_info%speciesId(s)
  if (sid.le.0) cycle
  do i = 1, size(eta,1)
    m = id2cmplxRefInd(sid, eta(i))
    kappaStar(i) = particleAbscStar(m, a, eta(i))
  end do
  kappa =kappa + fA*kappaStar*particle_info%SpeciesMassFrac(s) 
  ! weighted by mass fraction
end do
! missing species vol frac
end function genLBLparticleAbsc

function genGrayParticleScat(particle_info) result(scat)
use rad_solid_prop, only : id2cmplxRefInd
type(solidParInfoType), intent(in) :: particle_info
real(dp) :: scat, eta, a, fv, fA
! assume gray value taken at the peak spectrum which depends on
! temperature, as of 2011/DEC/09 this is the only one function
! uses temperature.
integer :: sid ! solid species id in database
! a (m), eta (cm^-1), output scat in (m^-1)
real(dp) :: kappaStar, betaStar, sigmaStar
integer :: s ! species loop variable
complex(dp) :: m ! complex refractive index
real(dp), parameter :: C3=1.d4/2898.d0 ! eqn 1.16 Modest 2003 book
eta = C3 * particle_info%T
a = particle_info%radius
fv = particle_info%volFrac
fA = 0.75d0*fv/a
scat=0.d0
do s = 1, particle_info%totSpecies
  sid = particle_info%speciesId(s)
  if (sid.le.0) cycle
  m = id2cmplxRefInd(sid, eta)
  kappaStar = particleAbscStar(m, a, eta)
  betaStar = particleExtcStar(m, a, eta)
  sigmaStar = betaStar-kappaStar
  scat  =scat  + fA*sigmaStar*particle_info%SpeciesMassFrac(s) 
  ! weighted by mass fraction
end do
end function genGrayParticleScat

function genPlanckMeanParticleAbsc(particle_info) result(abscp)
! calculate planck mean particle absorption coefficient from LBL calculation
type(solidParInfoType), intent(in) :: particle_info
real(dp) :: abscp
real(dp),parameter :: etaMin=100.d0, etaMax=15000.d0, etaStep=100.d0
integer, parameter :: etaLen = ceiling((etaMax-etaMin)/etaStep)
integer :: i
real(dp), dimension(etaLen) :: eta, kappa, BB
real(dp) :: Q

eta = (/(etaMin+(i-1)*etaStep, i=1, etaLen)/)
kappa  = genLBLparticleAbsc(particle_info, eta)
do i = 1, etaLen
  Q = particle_info%T * 1.d4/eta(i)
  BB(i) = BBFN(Q)
enddo ! bb decreasing
if (BB(etaLen)>2.d-2) print*, "Maximum wavenumber too small for Planck function"
abscp = BB(etaLen)*kappa(etaLen)
do i= (etaLen-1),1, -1
  abscp = abscp + (BB(i)-BB(i+1)) * kappa(i) 
enddo
end function genPlanckMeanParticleAbsc

!-----------------------------------------------------------------------
! black body emission function 
REAL(dp) FUNCTION BBFN(Q) 
!!$     ********************************************************************
!!$     *  This subroutine calculates the fractional blackbody             *
!!$     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
!!$     ********************************************************************
    REAL(dp) :: PI,CC,C2,EPS,V,EX,M,EM,VM,BM,Q
    REAL(dp),PARAMETER :: x1= 21.d0, x2= 75671.d0

    IF (Q <x1) THEN        ! if Q is too small, return 0.0
       BBFN=0.d0; RETURN
    ELSEIF (Q>x2) THEN     ! if Q is too large, return 1.0
       BBFN= 1.d0; RETURN
    ENDIF   

    PI=4.D0*DATAN(1.D0)
    CC=1.5D1/PI**4
    C2=1.4388D4
    EPS=1.D-16

    V=C2/Q  
    EX=DEXP(V)

    M=0
    BBFN=0.D0
    EM=1.D0 
5   M=M+1   
    VM=M*V  
    BM=(6.D0+VM*(6.D0+VM*(3.D0+VM)))/M**4
    EM=EM/EX
    BBFN=BBFN+BM*EM
    IF(VM**3*EM.GT.EPS) GOTO 5
    BBFN=CC*BBFN
    RETURN  
  END FUNCTION BBFN
!------------------------------------------------------------------------

end module rad_solid_spectral

