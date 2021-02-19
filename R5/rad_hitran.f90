!***********************************************************************
!   SUBROUTINE Getdata
module rad_hitran

   IMPLICIT NONE
   INTEGER   :: i,isotp,j
   DOUBLE PRECISION :: dummy
   integer :: dummy1, dmy
   CHARACTER (len=256) :: PATHCO2,CH4File, H2OFile, CO2File
   double precision :: h1, h2, h3, h4, h5, h6 
   double precision :: hp1, hp2, hp3, hp4, hp5, hp6 
   double precision :: sint, bair, bself, elow, expb
   DOUBLE PRECISION,PARAMETER   :: cav=6.022136d23  ! Avogadro constant
   DOUBLE PRECISION,PARAMETER   :: hh =6.62607015d-27 !erg s    ! Planck  constant
   DOUBLE PRECISION,PARAMETER   :: cc =2.99792458d10  !cm s−1   ! speed of light 
   DOUBLE PRECISION,PARAMETER   :: kk =1.380649d-16   !erg K−1  ! Boltzmann constant
   DOUBLE PRECISION,PARAMETER   :: cc2=1.4387769      !cm K     ! second radiation constant
   DOUBLE PRECISION, parameter  :: cbz = 1.38064852d-23 !J/K     ! botlzmann constant 

! data(i,1) = h1  = wavenumber  v
! data(i,2) = h2 = intensity   S
! data(i,3) = h3 = b_air      gamma_air 
! data(i,4) = h4 = b_self     gamma_self
! data(i,5) = h5 = E''        lower state energy 
! data(i,6) = h6 = exponent for b   n_air 

!  read the data file 
!    CO2File='/home/vkotteda/rad_model/5ece4919_CO2.par'
!    CH4File='/home/vkotteda/rad_model/5ebd5178_ch4.par'
!    H2OFile='/home/vkotteda/rad_model/5ebd5203_H2O.par'

  contains 
  subroutine abs_co2_old(wvn, Pres, Temp, QT, QT0, absfm) 
  double precision, intent(in)  ::  wvn, Pres, Temp, QT, QT0  
  double precision, intent(out) :: absfm
  double precision :: patm,T0 
  double precision :: TT0,dens, line_width, sij_t,sij_t0, sijT, fiv, kv, denum  
  
  CO2File='/home/vkotteda/rad_model/5ece4919_CO2.par'
  patm  = Pres/1.01325d5
  T0    = 296.0d0 
  open(11,file=CO2FILE) 
   i=1
   DO 
   READ(11,151) dummy1,isotp,h1
   If(isotp/=1) CYCLE
   if(h1.gt.wvn-20.and.h1.le.wvn) then 
     READ(11,150) dummy1,isotp,h1,h2,dummy,h3, h4,h5,h6
     READ(11,150) dummy1,isotp,hp1,hp2,dummy,hp3, hp4,hp5,hp6 
      ! write(*,*) dummy1,isotp,hp1,hp2,dummy,hp3, hp4,hp5,hp6 
      if(hp1.ge.wvn) then 
      ! write(*,*) h1, hp1, h2, hp2 
      ! linear interpolation 
        sint  = h2 + (h2-hp2)*(wvn-h1)/(h1-hp1)
        bair  = h3 + (h3-hp3)*(wvn-h1)/(h1-hp1)
        bself = h4 + (h4-hp4)*(wvn-h1)/(h1-hp1)
        elow  = h5 + (h5-hp5)*(wvn-h1)/(h1-hp1)
        expb  = h6 + (h6-hp6)*(wvn-h1)/(h1-hp1)
!        write(*,*) dummy1, isotp, wvn, sint, dummy, bair,bself, elow,expb
      exit
      end if 
   exit 
   end if 
!   If(hdata(i,5)<0.d0) THEN
!    hdata(i,5)=0.d0
!    write(*,*) 'E"<0 at ', hdata(i,1)
!    Endif
!   If(hdata(i,4)<1.e-5) hdata(i,4)=5.d0*hdata(i,3)
   i=i+1
   ENDDO
150 FORMAT(i2,i1,f12.6,2e10.3,f5.4,f5.3,f10.4,f4.2) 
151 FORMAT(i2,i1,f12.6) 
  4 i=i-1

    TT0=T0/Temp
    dens  = Pres/((8314.4621/30.0d0)* Temp) ! density = p/ R_spefic T)  kg/m^3

    kv = 0.0d0
    line_width=(bself)*patm*TT0**expb
    sij_t  = QT*exp(-cc2*elow/Temp) * (1.0d0 - exp(-cc2* wvn/Temp))
    sij_t0 = QT0*exp(-cc2*elow/T0)  * (1.0d0 - exp(-cc2* wvn/T0))
    ! line intensity
    sijT   = sint* sij_t/ sij_t0
    fiv    = 1.0d0/line_width  !line_width /(PI*((wvnumber - wvn)**2.0d0+ line_width**2.0d0))
    kv     =  kv + sijT*fiv!*dens
    denum     =  2.479d19*(296/temp) !88 /(cbz * 1d6*Temp)
    absfm = absfm + kv
    write(*,*) kv, absfm, (101325*1/(cbz*273*1d6)) !!*293))

  END SUBROUTINE abs_co2_old



end module rad_hitran  
!************************************************************************
