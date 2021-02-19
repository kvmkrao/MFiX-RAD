MODULE rad_mc_calabs
use rad_gas_species, only : gasInfoType, rad_gas_species_init, getGasInfo
use param, only: dimension_3
use functions, only : funijk
use compar, only: istart1, jstart1, kstart1,iend1, jend1, kend1
IMPLICIT NONE
INTEGER  :: lu,ifg
INTEGER,PARAMETER  :: rows=1420000

DOUBLE PRECISION   :: hdata(rows,7)
integer :: MOL, NSO82,ISO85  
contains 

subroutine calabslbl(neta,Temp, eta, T_g,abscoef)
!subroutine calabslbl(gasinfo,wvnumber,absfm)
   use tipartitionsum
   use rad_hitran 
!   type(gasinfotype), intent(in) :: gasinfo
   integer          :: neta 
   double precision :: Pres,Temp,Patm,wvnumber,tc
   double precision :: xmfr(3)
   double precision :: eta(neta), absfm(neta)
   double precision :: halflw,gammai, fiv,line_width,QT,QT0,gi, TmI
   double precision :: abund, gj, massmol, plen 
   integer          :: it,icl,l,js,jj,lines,imax,ifg,im1,ip1,iq,ijk,kn
   integer          :: i1, j1, k1 
   character(len=126) :: hitempfile
   double precision :: abscoef(dimension_3,neta), T_g(dimension_3) 

! PRESSURE IN atm
   Pres = 1.0d0 !gasinfo%P 
   patm  = Pres/1.01325d5  ! pref = 1 atm 
   absfm = 0.0d0
!   Temp  = 300.0d0 !gasinfo%T 
   absfm = 0.0d0
   plen  = 1.0d0 
   abscoef = 0.0d0 
!C     GO TO PARTICULAR MOLECULE:  
!      IF(MOL.EQ.1) THEN
         MOL     = 1 
         abund   = 9.97317d-1
         QT0     = 1.7458E+02
         massmol = 18.010565
         NSO82   = 161 
         ISO     = 1
         TmI     = Tmax(MOL,ISO)
         CALL QT_H2O(Temp,ISO,gi,QT)
         hitempfile = '5ebd5203_H2O.par'
         call Getdata(hitempfile,lines) !eta range 0.072059 - 25710.825500
         do i1=istart1,iend1
           do j1=jstart1, jend1
             do k1=kstart1, kend1
               ijk = funijk(i1,j1,k1)
               Tc  = T_g(ijk) 
               do kn=1,neta
                if(eta(kn).gt.0.08.and.eta(kn).lt.25710) then
                  call abs_xg(lines, eta(kn), Pres, Temp, plen, 0.1d0, QT, QT0, absfm(kn))
                  abscoef(ijk,kn) = absfm(kn) 
                end if 
               end do 
             end do 
           end do 
         end do 
!      ELSE IF(MOL.EQ.2) THEN
         MOL     = 2 
         abund   = 9.84204E-01                                                                                  
         QT0     = 2.8609E+02
         massmol = 43.989830
         NSO82   =  626
         ISO     = 1 
         TmI     = Tmax(MOL,ISO)
         hitempfile = '5ece4919_CO2.par'
         call QT_CO2(Temp,ISO,gi,QT)
         call Getdata(hitempfile,lines) ! 158.301811 - 14073.811564
         do i1=istart1,iend1 
           do j1=jstart1, jend1
             do k1=kstart1, kend1
               ijk = funijk(i1,j1,k1)
               Tc  = T_g(ijk) 
               do kn=1,neta
                if(eta(kn).gt.159.0.and.eta(kn).lt.14073.0) then 
                  call abs_xg(lines, eta(kn), Pres, Temp, plen, 0.1d0, QT, QT0, absfm(kn))
                  abscoef(ijk,kn) = absfm(kn) 
                end if 
               end do 
             end do 
           end do 
         end do 

!      ELSE IF(MOL.EQ.3) THEN
        IF (MOL.EQ.3) THEN 
         abund   = 9.92901E-01                                                                                                            
         ISO     = 46
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_NO(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.9) THEN
         abund   = 9.45678E-01                                                                                                               
         QT0     = 6.3403E+03
         massmol = 63.961901
         ISO     = 626
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_SO2(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.10) THEN
         abund   = 9.91616E-01                                                                                                             
         QT0     = 1.3577E+04
         massmol = 45.992904
         ISO     = 646
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_NO2(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.11) THEN
         abund   = 9.95872E-01
         QT0     = 1.7252E+03
         massmol = 17.026549
         ISO     = 4111
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_NH3(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.12) THEN
         abund   = 9.89110E-01
         QT0     = 2.1393E+05
         massmol = 62.995644
         ISO     = 146
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_NH3(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.13) THEN
         abund   = 9.97473E-01
         QT0     = 8.0348E+01
         massmol = 17.002740
         ISO     = 61
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_OH(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.14) THEN
         abund   = 9.99844E-01
         QT0     = 4.1469E+01
         massmol = 20.006229
         ISO     = 19
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HF(Temp,ISO,gi,QT) 
 
      ELSE IF(MOL.EQ.15) THEN
         abund   = 7.57587E-01
         QT0     = 1.6065E+02
         massmol = 35.976678
         ISO     = 15
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HCL(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.16) THEN
         abund   = 5.06781E-01
         QT0     = 2.0017E+02
         massmol = 79.926160
         ISO     = 19
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HBR(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.17) THEN
         abund   = 9.99844E-01
         QT0     = 3.8899E+02
         massmol = 127.912297
         ISO     = 17
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HI(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.18) THEN
         abund   = 7.55908E-01
         QT0     = 3.2746E+03
         massmol = 50.963768
         ISO     = 56
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CLO(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.19) THEN
         abund   = 9.37395E-01
         QT0     = 1.2210E+03
         massmol = 59.966986
         ISO     = 622
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_OCS(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.20) THEN
         abund   = 9.86237E-01
         QT0     = 2.8445E+03
         massmol = 30.010565
         ISO     = 126
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_H2CO(Temp,ISO,gi,QT) 
 
      ELSE IF(MOL.EQ.21) THEN
         abund   = 7.55790E-01
         QT0     = 1.9275E+04
         massmol = 51.971593
         ISO     = 165
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HOCL(Temp,ISO,gi,QT) 
 
      ELSE IF(MOL.EQ.22) THEN
         abund   = 9.92687E-01
         QT0     = 4.6710E+02
         massmol = 28.006148
         ISO     = 44
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_N2(Temp,ISO,gi,QT) 
 
      ELSE IF(MOL.EQ.23) THEN
        abund    = 9.85114E-01
         QT0     = 8.9220E+02
         massmol = 27.010899
         ISO     = 124
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HCN(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.24) THEN
         abund   = 7.48937E-01
         QT0     = 5.7916E+04
         massmol = 49.992328
         ISO     = 215
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CH3CL(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.25) THEN
         abund   = 9.94952E-01
         QT0     = 9.8480E+03
         massmol = 34.005480
         ISO     = 1661
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_H2O2(Temp,ISO,gi,QT) 
 
      ELSE IF(MOL.EQ.26) THEN
         abund   = 9.77599E-01
         QT0     = 4.1245E+02
         massmol = 26.015650
         ISO     = 1221
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C2H2(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.27) THEN
         abund   = 9.76990E-01
         QT0     = 7.0883E+04
         massmol = 30.046950
         ISO     = 1221
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C2H6(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.28) THEN
         abund   = 9.99533E-01
         QT0     = 3.2494E+03
         massmol = 33.997238
         ISO     = 1111
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_PH3(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.29) THEN
         abund   = 9.86544E-01
         QT0     = 7.0028E+04
         massmol = 65.991722
         ISO     = 269
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_COF2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.30) THEN
         abund   = 0.950180E+00
         QT0     = 1.6233E+06
         massmol = 145.962492
         ISO     = 29
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_SF6(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.31) THEN
         abund   = 9.49884E-01
         QT0     = 5.0579E+02
         massmol = 33.987721
         ISO     = 121
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_H2S(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.32) THEN
         abund   = 9.83898E-01
         QT0     = 3.9133E+04                                                                                                                
         massmol = 46.005480
         ISO     = 126
         TmI = Tmax(MOL,ISO)
         CALL QT_HCOOH(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.33) THEN
         abund   = 9.95107E-01
         QT0     = 4.3004E+03
         massmol = 32.997655
         ISO     = 166
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HO2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.34) THEN 
!...not applicable to O
      gi=0
      QT = 0.

      ELSE IF(MOL.EQ.35) THEN
         abund   = .749570E+00
         QT0     = 4.7884E+06
         massmol = 96.956672
         ISO     = 5646
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_ClONO2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.36) THEN 
         abund   = 9.93974E-01
         QT0     = 3.1169E+02
         massmol = 29.997989
         ISO     = 46
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_NOp(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.37) THEN 
         abund   = 5.05579E-01
         QT0     = 2.8339E+04
         massmol = 95.921076
         ISO     = 169
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HOBr(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.38) THEN
         abund   = 9.77294E-01
         QT0     = 1.1042E+04
         massmol = 28.031300
         ISO     = 221
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C2H4(Temp,ISO,gi,QT)
 
      ELSE IF(MOL.EQ.39) THEN
         abund   = 9.85930E-01
         QT0     = 7.0570E+04
         massmol = 32.026215
         ISO     = 2161
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CH3OH(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.40) THEN
         abund   = 5.00995E-01
         QT0     = 8.3052E+04
         massmol = 93.941811
         ISO     = 219
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CH3Br(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.41) THEN
         abund   = 9.73866E-01
         QT0     = 8.8672E+04
         massmol = 41.026549
         ISO     = 2124
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CH3CN(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.42) THEN 
!c...no data for CF4
      gi=0
      QT = 0.
!      go to 100

      ELSE IF(MOL.EQ.43) THEN
         abund   = 9.55998E-01
         QT0     = 9.8190E+03
         massmol = 50.015650
         ISO     =  2211
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C4H2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.44) THEN
         abund   = 9.63346E-01
         QT0     = 2.4787E+04
         massmol = 51.010899
         ISO     = 1224
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_HC3N(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.45) THEN
         abund   = 9.99688E-01
         QT0     = 7.6712E+00
         massmol = 2.015650
         ISO     = 11
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_H2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.46) THEN
         abund   = 9.39624E-01
         QT0     = 2.5362E+02
         massmol = 43.971036
         ISO     = 22
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CS(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.47) THEN
         abund   = 9.43400E-01
         QT0     = 7.7833E+03
         massmol = 79.956820
         ISO     = 26
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_SO3(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.48) THEN
         abund   = 9.70752E-01
         QT0     = 1.5582E+04
         massmol = 52.006148
         ISO     = 4224
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C2N2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.49) THEN
         abund   = 5.66392E-01
         QT0     = 1.4800E+06
         massmol = 97.932620 
         ISO     = 2655
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_COCL2(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.50) THEN
         abund   = 9.97317d-1
         QT0     = 0.7458d2
         massmol = 18.010565
         ISO     = 161
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_SO(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.51) THEN
         abund   = 9.97317d-1
         QT0     = 0.7458d2
         massmol = 18.010565
         ISO     = 161
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_C3H4(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.52) THEN
         abund   = 9.97317d-1
         QT0     = 0.7458d2
         massmol = 18.010565
         ISO     = 161
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CH3(Temp,ISO,gi,QT)

      ELSE IF(MOL.EQ.53) THEN
         abund   = 8.92811E-01
         QT0     = 1.3526E+03
         massmol = 75.944140
         ISO     = 222
         TmI = Tmax(MOL,ISO)
         ! TOTAL INTERNAL PARTITION FUNCTION
         CALL QT_CS2(Temp,ISO,gi,QT)
      ENDIF     
   END 

   subroutine calabssnb(gasinfo, Temp, eta, abs_coef)
   use fldvar, only : ep_g, T_g, P_g, X_g
   use rad_param
   ! This function uses the Elsasser Narrow Band model (Tong 1992 with
   ! parameters from Edwards 1967) to calculate the total spectral
   ! absorption coefficient (m^-1 ) of a CO2/H20/N2 mixture given the
   ! current wavenumber (cm^-1), gas temperature (K) , total pressure (atm) ,
   ! mole fraction of CO2, and mole fraction of H20 - it is assumed
   ! that the remaining mole fraction consitutes the broadening N2 gas k.
   !The datails about the model can be found in 
   !Li, Weiming, Timothy W. Tong, Dean Dobranich, and Louis A. Gritzo. "A combined narrow-and wide-band model for computing the spectral absorption coefficient of CO2, CO, H2O, CH4, C2H2, and NO." Journal of Quantitative Spectroscopy and Radiative Transfer 54, no. 6 (1995): 961-970. (DOI: https://doi.org/10.1016/0022-4073(95)00126-6)

   implicit none
   include 'usrnlst.inc'
   ! Initialize Variables
   type(gasInfoType), intent(in) :: gasinfo
!   integer, intent(in) :: ijk 
   double precision, intent(in) ::   eta
   double precision, intent(out) :: abs_coef
   double precision :: Temp , Pres , F_CO2 , F_H2O, TT0,hcbkt, P0  
   !ch - planck constant    = 6.62607015×10−27 erg⋅s (CGS) 
   !cc - speed of sound     = 2.998 10^10 cm/s (CGS) 
   !ck - Boltzmann constant = 1.380649×10−16 erg⋅K−1
   double precision, parameter ::  ch = 6.6237d-27, ck = 1.3802d-16, &
   cc = 2.9979d10 , T0 = 100d0 
   integer :: i
   double precision, dimension(5) :: &
   eta_c_c = (/667d0,960d0,1060d0,2350d0,3715d0/), &
   n_c = (/0.7d0,0.8d0,0.8d0,0.8d0,0.65d0/) , &
   a_c = (/2d0,2d0,2d0,1d0,2d0/), &
   delta_c = (/100d0,20d0,30d0,30d0,1000d0/),&
   b_c, C1_c , C2_c , C3_c , C3_T0_c

   double precision, dimension(4) :: &
   eta_c_h = (/1600d0,3750d0,5350d0,7250d0/), &
   n_h = (/1.0d0,1.0d0,1.0d0,1.0d0/) , &
   delta_h = (/250d0,100d0,5d0,200d0/), & 
   a_h,b_h, C1_h , C2_h , C3_h , C3_T0_h

   double precision :: rho_CO2 , rho_H2O , Pb , Pe , v1_c , v3_c , &
   phi1_c , phi2_c , phi3_c , v1_h , v2_h , v3_h , phi_011_h , phi_101_h , &
   phi4_h , Bv , delta , Sc_d

   P0  = 1.0d0 !101325.0d0  ! reference pressure 
   b_c = 1.3d0       
   b_h = 5.0d0 
!   a_c = 2.0d0 ! a = 1 for asymmetric bands and 2 for symmetric bands 
   a_h = 2.0d0 
   ! Calculate Densities
   Pres = 1.0d0 ! 101325.0 !P_g(ijk) !gasinfo%P
!   Temp =          !T_g(ijk) !gasinfo%T
   F_CO2 = 0.1d0 !gasinfo%C(CO2)
   F_H2O = 0.2d0 !gasinfo%C(H2O)
   
!   write(*,*) pres, Temp, F_CO2, F_H2O
   !          MM*            P/     (0.0821 (L atm/mol K) * Temp)   ! 1lt = 1000 cm^3 
   rho_CO2 = 44.0095d0*F_CO2*(Pres/101325d0)/0.082057d0/Temp/1000d0 ! g/L or kg/m^3 !*1000d0  
   rho_H2O = 18.0153d0*F_H2O*(Pres/101325d0)/0.082057d0/Temp/1000d0 ! g/L or kg/m^3 !*1000d0 
   !write(*,*) rho_CO2, rho_H2O, Pres, F_CO2, F_H2O, Temp
   ! Calculate N2 Partial Pressure
   
   Pb = Pres *(1d0-F_CO2-F_H2O)
   
   !CO2 Constant Phi Functions
   v1_c = 1351d0
   v3_c = 2396d0
   hcbkt = ch*cc/(ck*Temp) 
   phi1_c = (1.0d0-exp(-hcbkt*(v3_c-v1_c))) * (exp(-hcbkt*v1_c) - &
             0.5d0*exp(-2d0*hcbkt*v1_c)) / ( (1d0-exp(-hcbkt*v1_c)) * & 
                                             (1d0-exp(-hcbkt*v3_c)) )
   phi2_c = (1d0 - exp(-hcbkt*(v1_c+v3_c)))/((1d0-exp(-hcbkt*v1_c)) *  &
                                             (1d0-exp(-hcbkt*v3_c)))
   phi3_c = 1.0d0 + 0.053d0*(Temp/T0)**1.5d0
   
   !write(*,*) phi1_c, phi2_c, phi3_c, Temp, T0 
   !H2O Constant Phi Functions
   v1_h = 3652d0
   v2_h = 1595d0
   v3_h = 3756d0
   phi_011_h = (1d0- exp(-hcbkt*(v2_h+v3_h)))/ ( (1d0- exp(-hcbkt*v1_h)) &
            *  (1d0- exp(-hcbkt*v2_h))            *   (1d0- exp(-hcbkt*v3_h)) )
   phi_101_h = (1d0- exp(-hcbkt*(v1_h+v3_h)))/ ( (1d0- exp(-hcbkt*v1_h)) &
            *  (1d0- exp(-hcbkt*v2_h))            *   (1d0- exp(-hcbkt*v3_h)) )
!   phi4_h = exp(-17.6d0*(Temp/T0)**(-0.5))
   
   !write(*,*) phi_011_h, phi_101_h, phi4_h
   !hc/k = 1.4388 cm K 
   !CO2 Constants
   TT0 = Temp/T0
   C1_c  = (/19d0, 0.76d0*phi1_c,0.76d0*phi1_c,110d0, 4d0*phi2_c/)
   C2_c  = (/6.9d0*(TT0)**0.5d0,1.6d0*(TT0)**0.5d0*C1_c(2)**0.5d0, &
             1.6d0*(TT0)**0.5d0*C1_c(3)**0.5d0,31d0*(TT0)**0.5d0,8.6d0*phi3_c/)
   C3_c  = (/12.9d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,& 
             11.5d0*(TT0)**0.5d0,24d0*(TT0)**0.5d0/)
!   C3_T0_c = (/12.9d0*(T0/T0)**0.5d0,12.4d0*(T0/T0)**0.5d0,12.4d0* &
!              (T0/T0)**0.5d0,11.5d0*(T0/T0)**0.5d0,24d0*(T0/T0)**0.5d0/)
   !H2O Constants
   C1_h = (/41.2d0,23.3d0,3d0*phi_011_h,2.5d0*phi_101_h/)
   C2_h = (/44d0,39d0,6d0*C1_h(3)**0.5d0,8d0*C1_h(4)**0.5d0/)
   C3_h = (/52d0*(TT0)**0.5d0,65d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0/)
!   C3_T0_h = (/52d0*(T0/T0)**0.5d0,65d0*(T0/T0)**0.5d0,46d0*(T0/T0)**0.5d0,46d0*(T0/T0)**0.5d0 /)
   ! Calculate CO2 Contribution to Absorption Coefficient
   
   !write(*,*) phi_011_h, phi_101_h, phi4_h
   abs_coef = 0d0
   do i = 1 , size(eta_c_c)
      Pe = ((Pb + (b_c(i)-1.0d0)*Pres*F_CO2)/ P0 )**n_c(i)  ! Pe 
      Bv = (C2_c(i)**2d0*Pe)/(4d0*C1_c(i)*C3_c(i)) ! beta 
      delta = delta_c(i) !30d0*C3_T0_c(i)                       ! delta 
      Sc_d  =C1_c(i)/C3_c(i)* exp(-(a_c(i)/C3_c(i))* abs(eta - eta_c_c(i)))
      abs_coef = abs_coef + rho_CO2*(Sc_d)* ( (sinh(pi*Bv/2d0))  / &
      (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta-eta_c_c(i))/delta) ) )
      !write(*,*) Pe, Bv, delta, C1_c(i), C3_c(i),Sc_d, pi, eta, eta_c_c(i), (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta - eta_c_c(i)) / delta)), abs_coef 
   end do
!   stop  
   ! Calculate H20 Contribution to Absorption Coefficient
   do i = 1 , size(eta_c_h)
      Pe = ((Pb + (b_h(i)-1.0d0)*Pres*F_H2O)/ P0)**n_h(i)
      Bv =  (C2_h(i)**2d0*Pe)/(4d0*C1_h(i)* C3_h(i))
      delta = delta_h(i) !30d0*C3_T0_h(i)
      Sc_d  =C1_h(i)/C3_h(i)* exp(-(a_h(i)/C3_h(i))* abs(eta - eta_c_h(i)))
      abs_coef = abs_coef + rho_H2O*(Sc_d)* ( (sinh(pi*Bv/2d0))  / &
      (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta-eta_c_h(i))/delta) ) )
   end do
   !write(*,*) "abs", abs_coef
   stop 
   end subroutine calabssnb

   SUBROUTINE Getdata(hirdfl,i)
   IMPLICIT NONE
   INTEGER   :: i,isotp,j
   INTEGER   :: dummy1
   DOUBLE PRECISION :: dummy
   CHARACTER (len=10) :: dummy2
   character (len=126) :: hirdfl 

! hdata(i,1) = wavenumber  v_ij (cm^-1)
! hdata(i,2) = intensity   S_ij (cm^-1/mole.cm^-2)
! hdata(i,3) = b_air      gamma_air  (cm^-1 atm^-1)
! hdata(i,4) = b_self     gamma_self (cm^-1 atm^-1)
! hdata(i,5) = E''        lower state energy (cm^-1) 
! hdata(i,6) = exponent for b   n_air 
! hdata(i,7) = delta  pressure shift   (cm^-1 atm^-1)

   open(101,file=hirdfl,status="unknown")  
   i=1
   DO
   READ(101,FMT=160,END=4) dummy1,isotp,hdata(i,1),hdata(i,2),dummy,hdata(i,3),    &
    hdata(i,4),hdata(i,5),hdata(i,6),hdata(i,7)
!   write(*,*) i,dummy1,isotp,hdata(i,1),hdata(i,2),dummy,(hdata(i,3),j=3,6)
!   If(hdata(i,1)<wvmin) CYCLE
!   If(hdata(i,2)>wvmax) GOTO 4
! Ignore if not the most abundent isotope
   If(isotp/=1) CYCLE
!   If(hdata(i,5)<0.d0) THEN
!      hdata(i,5)=0.d0
!      write(*,*) 'E"<0 at ', hdata(i,1)
!   Endif
   If(hdata(i,4)<1.e-5) hdata(i,4)=5.d0*hdata(i,3)
   i=i+1
   ENDDO
160 FORMAT(i2,i1,f12.6,2e10.3,f5.4,f5.3,f10.4,f4.2,f8.6)
!160 FORMAT(i2,i1,f12.6,es10.3,a10,f5.4,f5.3,f10.4,f4.2) commented Murali 
! The above format specification is for HITRAN 2008
! For HITRAN 96 use the following
!160 FORMAT(i3,f12.6,e10.3,a10,f5.4,f5.4,f10.4,f4.2)
  4 i=i-1
  close(101) 
  !5 RETURN
   END SUBROUTINE Getdata

 subroutine abs_xg(nlines,wvn, Pres, Temp, plen, xco2, QT, QT0, absfm)

  double precision ::  wvn, Temp, QT, QT0, xco2   
  double precision, intent(out):: absfm
  integer          :: nlines 
  double precision :: patm,T0, ypt, vsh,plen,Pres 
  double precision :: TT0,dens, line_width, sij_t,sij_t0, sijT, fiv, kv, denum  
  integer          :: i
! Input parameters
  double precision, parameter :: PI =3.1415926d0
  double precision, parameter :: cav=6.022136d23              ! Avogadro constant
  double precision, parameter :: hh =6.62607015d-27 !erg s    ! Planck  constant
  double precision, parameter :: cc =2.99792458d10  !cm s−1   ! speed of light 
  double precision, parameter :: kk =1.380649d-16   !erg K−1  ! Boltzmann constant
  double precision, parameter :: cc2=1.4387769      !cm K     ! second radiation constant
  double precision, parameter :: cbz = 1.38064852d-23 !J/K    ! botlzmann constant 

! The wavenumber of the spectral line transition (cm-1) in vacuum
! hdata(i,1) = wavenumber  v_ij (1/cm)

!The spectral line intensity (cm−1/(molecule·cm−2)) at Tref=296K
! hdata(i,2) = intensity   S_ij          ( cm^-1 /(molecule.cm^-2) at Trf = 296K 

!The air-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm
! hdata(i,3) = b_air      gamma_air      HWHM (cm^-1 /atm) p_ref = 1atm 

! The self-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm
! hdata(i,4) = b_self     gamma_self     

!The lower-state energy of the transition (cm-1)
! hdata(i,5) = E''

! The coefficient of the temperature dependence of the air-broadened half width
! hdata(i,6) = exponent for b   n_air    

!The pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to the vacuum transition wavenumber νij
!    hdata(i,7)  pressure shift 

  patm  = 1.0d0   !Pres/1.01325d5 in atm 
  T0    = 296.0d0 ! Kelvin 
  Pres  = 1.0d0   ! in atm 

  TT0=T0/Temp    
  ! column density U = 1.01325d-4 * plen * (xco2*patm) / (8.314*Temp)  
  ! conversion factor 1.01325d-4 
  dens  = 1.01325d-4*plen*(xco2*Pres)/(8.3144621* Temp) ! KiloMoles/cm^2
  !dens  = Pres*plen*(xco2*Pres)/(8.3144621* Temp) ! KiloMoles/cm^2

  kv = 0.0d0
  do i=1,nlines !,50 
!The Lorentzian (pressure-broadened) HWHM, γ(p,T) for a gas at pressure p (atm), temperature T (K) and partial pressure pself (atm) is calculated as:  
!see https://hitran.org/docs/definitions-and-units/
  ypt       = (TT0**hdata(i,6))* (hdata(i,3)*(1.0d0 - xco2)* Pres + hdata(i,4)*xco2*Pres)

! The pressure shift, δ, of the transition wavenumber leads to a shifted position ν∗ij given by
  vsh       = hdata(i,1) + hdata(i,7)*Pres 
!  line_width= (hdata(i,4))*patm* TT0**expb
  sij_t  = QT *exp(-cc2*hdata(i,5)/Temp)* (1.0d0 - exp(-cc2* hdata(i,1)/Temp))
  sij_t0 = QT0*exp(-cc2*hdata(i,5)/T0)  * (1.0d0 - exp(-cc2* hdata(i,1)/T0))
  ! line intensity
  sijT   = hdata(i,2)* sij_t/ sij_t0     !cm^-1/(moleculescm^-2) 
  fiv    = ypt /(PI*(ypt**2.0d0 + (wvn- vsh)**2.0d0)) ! cm or 1/cm^-1  
  kv     = kv + sijT*fiv*dens*6.023d26                ! NA = 6.02214076×10^23 mol−1
  end do 
  absfm = absfm + kv

  END SUBROUTINE abs_xg


   subroutine calabssnb_old(gasinfo, Temp, eta, abs_coef)
   use fldvar, only : ep_g, T_g, P_g, X_g
   use rad_param
   ! This function uses the Elsasser Narrow Band model (Tong 1992 with
   ! parameters from Edwards 1967) to calculate the total spectral
   ! absorption coefficient (m^-1 ) of a CO2/H20/N2 mixture given the
   ! current wavenumber (cm^-1), gas temperature (K) , total pressure (atm) ,
   ! mole fraction of CO2, and mole fraction of H20 - it is assumed
   ! that the remaining mole fraction consitutes the broadening N2 gas k.
   !The datails about the model can be found in 
   !Li, Weiming, Timothy W. Tong, Dean Dobranich, and Louis A. Gritzo. "A combined narrow-and wide-band model for computing the spectral absorption coefficient of CO2, CO, H2O, CH4, C2H2, and NO." Journal of Quantitative Spectroscopy and Radiative Transfer 54, no. 6 (1995): 961-970. (DOI: https://doi.org/10.1016/0022-4073(95)00126-6)

   implicit none
   include 'usrnlst.inc'
   ! Initialize Variables
   type(gasInfoType), intent(in) :: gasinfo
!   integer, intent(in) :: ijk 
   double precision, intent(in) ::   eta
   double precision, intent(out) :: abs_coef
   double precision :: Temp , Pres , F_CO2 , F_H2O, TT0,hcbkt, P0  
   !ch - planck constant    = 6.62607015×10−27 erg⋅s (CGS) 
   !cc - speed of sound     = 2.998 10^10 cm/s (CGS) 
   !ck - Boltzmann constant = 1.380649×10−16 erg⋅K−1
   double precision, parameter ::  ch = 6.6237d-27, ck = 1.3802d-16, &
   cc = 2.9979d10 , T0 = 100d0 
   integer :: i
   double precision, dimension(5) :: &
   eta_c_c = (/667d0,962d0,1064d0,2326d0,3703d0/), &
   n_c = (/0.7d0,0.8d0,0.8d0,0.8d0,0.65d0/) , &
   a_c = (/2d0,2d0,2d0,1d0,2d0/), &
   b_c, C1_c , C2_c , C3_c , C3_T0_c

   double precision, dimension(5) :: &
   eta_c_h = (/500d0,1587d0,3704d0,5348d0,7247d0/), &
   n_h = (/1.0d0,1.0d0,1.0d0,1.0d0,1.0d0/) , &
   a_h,b_h, C1_h , C2_h , C3_h , C3_T0_h

   double precision :: rho_CO2 , rho_H2O , Pb , Pe , v1_c , v3_c , &
   phi1_c , phi2_c , phi3_c , v1_h , v2_h , v3_h , phi_011_h , phi_101_h , &
   phi4_h , Bv , delta , Sc_d

   P0  = 1.0d0  ! reference pressure 
   b_c = 1.3d0       
   b_h = 5.0d0 
!   a_c = 2.0d0 ! a = 1 for asymmetric bands and 2 for symmetric bands 
   a_h = 2.0d0 
   ! Calculate Densities
   Pres = 1.0d0 !101325.0 !P_g(ijk) !gasinfo%P
!   Temp =          !T_g(ijk) !gasinfo%T
   F_CO2 = 0.1d0 !gasinfo%C(CO2)
   F_H2O = 0.2d0 !gasinfo%C(H2O)
   
!   write(*,*) pres, Temp, F_CO2, F_H2O
   !          MM*            P/     (0.0821 (L atm/mol K) * Temp)   ! 1lt = 1000 cm^3 
   rho_CO2 = 44.0095d0*F_CO2*(Pres)/0.082057d0/Temp*1000d0 ! g/L or kg/m^3 !*1000d0  
   rho_H2O = 18.0153d0*F_H2O*(Pres)/0.082057d0/Temp*1000d0 ! g/L or kg/m^3 !*1000d0 
   !write(*,*) rho_CO2, rho_H2O, Pres, F_CO2, F_H2O, Temp
   ! Calculate N2 Partial Pressure
   
   Pb = Pres *(1d0-F_CO2-F_H2O)
   
   !CO2 Constant Phi Functions
   v1_c = 1351d0
   v3_c = 2396d0
   hcbkt = ch*cc/(ck*Temp) 
   phi1_c = (1.0d0-exp(-hcbkt*(v3_c-v1_c))) * (exp(-hcbkt*v1_c) - &
             0.5d0*exp(-2d0*hcbkt*v1_c)) / ( (1d0-exp(-hcbkt*v1_c)) * & 
                                             (1d0-exp(-hcbkt*v3_c)) )
   phi2_c = (1d0 - exp(-hcbkt*(v1_c+v3_c)))/((1d0-exp(-hcbkt*v1_c)) *  &
                                             (1d0-exp(-hcbkt*v3_c)))
   phi3_c = 1.0d0 + 0.053d0*(Temp/T0)**1.5d0
   
   !write(*,*) phi1_c, phi2_c, phi3_c, Temp, T0 
   !H2O Constant Phi Functions
   v1_h = 3652d0
   v2_h = 1595d0
   v3_h = 3756d0
   phi_011_h = (1d0- exp(-hcbkt*(v2_h+v3_h)))/ ( (1d0- exp(-hcbkt*v1_h))/ &
               (1d0- exp(-hcbkt*v2_h))            *   (1d0- exp(-hcbkt*v3_h)) )
   phi_101_h = (1d0- exp(-hcbkt*(v1_h+v3_h)))/ ( (1d0- exp(-hcbkt*v1_h)) /&
               (1d0- exp(-hcbkt*v2_h))            *   (1d0- exp(-hcbkt*v3_h)) )
   phi4_h = exp(-17.6d0*(Temp/T0)**(-0.5))
   
   !write(*,*) phi_011_h, phi_101_h, phi4_h
   !hc/k = 1.4388 cm K 
   !CO2 Constants
   TT0 = Temp/T0
   C1_c  = (/19d0, 0.76d0*phi1_c,0.76d0*phi1_c,110d0, 4d0*phi2_c/)
   C2_c  = (/6.9d0*(TT0)**0.5d0,1.6d0*(TT0)**0.5d0*C1_c(2)**0.5d0, &
             1.6d0*(TT0)**0.5d0*C1_c(3)**0.5d0,31d0*(TT0)**0.5d0,8.6d0*phi3_c/)
   C3_c  = (/12.9d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,& 
             11.5d0*(TT0)**0.5d0,24d0*(TT0)**0.5d0/)
   C3_T0_c = (/12.9d0,12.4d0,12.4d0,11.5d0,24d0/)
   !H2O Constants
   C1_h = (/5200d0*phi4_h,41.2d0,23.3d0,3d0*phi_011_h,2.5d0*phi_101_h/)
   C2_h = (/3.8d0*C1_h(1)**0.5d0,44d0,39d0,6d0*C1_h(4)**0.5d0,8d0*C1_h(5)**0.5d0/)
   C3_h = (/28.5d0*(TT0)**0.5d0,52d0*(TT0)**0.5d0,65d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0/)
   C3_T0_h = (/28.5d0,52d0,65d0,46d0,46d0 /)
   ! Calculate CO2 Contribution to Absorption Coefficient
   
   !write(*,*) phi_011_h, phi_101_h, phi4_h
   abs_coef = 0d0
   do i = 1 , size(eta_c_c)
      Pe = ((Pb + b_c(i)*Pres*F_CO2)/ P0)**n_c(i)  ! Pe 
      Bv = (C2_c(i)**2d0*Pe)/(4d0*C1_c(i)*C3_c(i)) ! beta 
      delta = 30d0*C3_T0_c(i)                       ! delta 
      Sc_d  =C1_c(i)/C3_c(i)* exp(-(a_c(i)/C3_c(i))* abs(eta - eta_c_c(i)))
      abs_coef = abs_coef + rho_CO2*Sc_d*( (sinh(pi*Bv/2d0))  / &
      (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta-eta_c_c(i))/delta) ) )
      write(*,*) Pe, Bv, delta, C1_c(i), C3_c(i),Sc_d, pi, eta, eta_c_c(i), (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta - eta_c_c(i)) / delta)), abs_coef 
   end do
!   stop  
   ! Calculate H20 Contribution to Absorption Coefficient
   do i = 1 , size(eta_c_h)
      Pe = ((Pb + b_h(i)*Pres*F_H2O)/ P0)**n_h(i)
      Bv =  (C2_h(i)**2d0*Pe)/(4d0*C1_h(i)* C3_h(i))
      delta = 30d0*C3_T0_h(i)
      Sc_d  =C1_h(i)/C3_h(i)* exp(-(a_h(i)/C3_h(i))* abs(eta - eta_c_h(i)))
      abs_coef = abs_coef + rho_H2O*Sc_d* ( (sinh(pi*Bv/2d0))  / &
      (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta-eta_c_h(i))/delta) ) )
   end do
   write(*,*) "abs", abs_coef
!   stop 
   end subroutine calabssnb_old
!!   subroutine calabsfsk()
!!! CODE TO CALCULATE FSK DISTRIBUTION FOR THE FULL SPECTRUM [k(g) and a(T_i,g)]
!!! FOR MIXTURES OF CO2(ifg=1), H2O(ifg=2) AND CH4(ifg=3), WITH OR WITHOUT SOOT.
!!!
!!! Input parameters (set by changing values in front matter):
!!!   Tref:   reference temperature (K)
!!!   Tmin:   minimum temperature for a-function (K)
!!!   Tmax:   maximum temperature for a-function (K)
!!!   numT:   number of temperatures considered for a (equally spaced Tmin-->Tmax)
!!            1.6d0*(TT0)**0.5d0*C1_c(3)**0.5d0,31d0*(TT0)**0.5d0,8.6d0*phi3_c/)
!!   C3_c  = (/12.9d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,12.4d0* &
!!          (TT0)**0.5d0,11.5d0*(TT0)**0.5d0,24d0*(TT0)**0.5d0/)
!!   C3_T0_c = (/12.9d0*(TT0)**0.5d0,12.4d0*(TT0)**0.5d0,12.4d0* &
!!              (TT0)**0.5d0,11.5d0*(TT0)**0.5d0,24d0*(TT0)**0.5d0/)
!!   !H2O Constants
!!   C1_h = (/5200d0*phi4_h,41.2d0,23.3d0,3d0*phi_011_h,2.5d0*phi_101_h/)
!!   C2_h = (/3.8d0*C1_h(1)**0.5d0,44d0,39d0,6d0*C1_h(4)**0.5d0, &
!!              8d0*C1_h(5) **0.5d0/)
!!   C3_h = (/28.5d0*(TT0)**0.5d0,52d0*(TT0)**0.5d0,65d0*(TT0)**0.5d0, &
!!              46d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0/)
!!   C3_T0_h = (/28.5d0*(TT0)**0.5d0,52d0*(TT0)**0.5d0, &
!!            65d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0,46d0*(TT0)**0.5d0 /)
!!   ! Calculate CO2 Contribution to Absorption Coefficient
!!   
!!   !write(*,*) phi_011_h, phi_101_h, phi4_h
!!   abs_coef = 0d0
!!   do i = 1 , size(eta_c_c)
!!      Pe = ((Pb + b_c(i)*Pres*F_CO2)/ P0 )**n_c(i)  ! Pe 
!!      Bv =  (C2_c(i)**2d0*Pe)/(4d0*C1_c(i)*C3_c(i)) ! beta 
!!      delta = 30d0*C3_T0_c(i)                       ! delta 
!!      Sc_d = C1_c(i)/C3_c(i)* exp(-(a_c(i)/C3_c(i))* abs(eta - eta_c_c(i)))
!!      abs_coef = abs_coef + (rho_CO2*(Sc_d/delta)* sinh(pi*Bv/2d0))  / &
!!      (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta - eta_c_c(i))))
!!      write(*,*) Pe, Bv, delta, Sc_d, pi, eta, eta_c_c(i), (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta - eta_c_c(i)) / delta))
!!   end do
!!   
!!   
!!   ! Calculate H20 Contribution to Absorption Coefficient
!!   do i = 1 , size(eta_c_h)
!!      Pe = ((Pb + b_h(i)*Pres*F_CO2)/ P0)**n_h(i)
!!      Bv =  (C2_h(i)**2d0*Pe)/(4d0*C1_h(i)* C3_h(i))
!!      delta = 30d0*C3_T0_h(i)
!!      Sc_d = C1_h(i)/ C3_h(i)* exp(-(a_h(i)/ C3_h(i))*abs(eta - eta_c_h(i)))
!!     abs_coef = abs_coef + (rho_H2O*Sc_d* sinh(pi*Bv/2d0))! / &
!!   !   (cosh(pi*Bv/2d0) - cos(2d0*pi*(eta - eta_c_h(i)) / delta))
!!   end do
!!   end subroutine calabssnb

!   subroutine calabsfsk()
! CODE TO CALCULATE FSK DISTRIBUTION FOR THE FULL SPECTRUM [k(g) and a(T_i,g)]
! FOR MIXTURES OF CO2(ifg=1), H2O(ifg=2) AND CH4(ifg=3), WITH OR WITHOUT SOOT.
!
! Input parameters (set by changing values in front matter):
!   Tref:   reference temperature (K)
!   Tmin:   minimum temperature for a-function (K)
!   Tmax:   maximum temperature for a-function (K)
!   numT:   number of temperatures considered for a (equally spaced Tmin-->Tmax)
!   P:      total pressure of mixture (bar)
!   xmfr(ifg): mole fraction of component ifg
!   fvsoot: soot volume fraction (-)
!   nsoot,ksoot: soot complex index of refraction (assumed constant across spectrum)
!   wvnm_b: minimum wavenumber to be considered (cm-1)
!   wvnm_e: maximum wavenumber to be considered (cm-1)
!   wvnmst: wavenumber step, equally spaced (cm-1)
!   klmin:  minimum absorption coefficient evaluated for individual lines (because of
!           line overlap klmin should be 2-3 orders of magnitude lower than desired accuracy)
!   n_pwrk: number of k-boxes for k-distribution
!   pwr:    exponent for setting of k-box values; in equal steps of k^pwr
!   kdmin:  minimum k-value to be considered for k-distribution (values<1E-9 force kdmin=1E-9)
!   kdmax:  maximum k-value to be considered for k-distribution (kdmax<=0 sets kdmax=kmax,
!           ie, the maximum absorption coefficient found across the spectrum)
!     Warning: setting kdmin too low and/or kdmax too large will cause many empty k-boxes!
!   nq:     number of quadrature points for radiative calculations
! Switches:
!   iwr=0: single run for absco and k and a(g), no absco stored
!   iwr=1: calculate absco, then k and a(g), and absco stored
!   iwr=2: absco from file, calculate  k and a(g)
!   ipl=0: calculate linear absorption coefficient
!   ipl=1: calculate pressure-based absorption coefficient (only for single absorbing gas!)
!   ipr=1: single output file containing k,g(T,k),a(T,k) for all values of T and k
!   ipr=2: an additional data file containing wq,gq(Tref),kq(gq),aq(T,kq) for all values of T for the
!          nq quadrature points gq (and Gaussian weights wq)
! Notes:
!   If the absorption coefficient absco is calculated (iwr=0 or 1), usually the
!   LINEAR absorption coefficient is calculated (ipl=0); sometimes--FOR A PURE GAS ONLY--
!   the PRESSURE-BASED absorption coefficient is desired (ipl=1). If the pressure-based
!   absorption coefficient for a dilute gas is desired, set xmfr=1.d-3 (=0.1%)
!
! Files (names and paths may be changed):
!   xxx/co2.dat, xxx/h2o.dat, xxx/ch4.dat: HITRAN/HITEMP data files if required
!   absco.dat: output file for absorption coefficient absco if iwr=0 or iwr=1
!              input file for absorption coefficient absco if iwr=2
!   k+avsg.dat: output file (in Tecplot format), providing:
!               k, f(k), g0(k), numT*g(T,k), numT*a(T,g0)
!   DOUBLE PRECISION   :: data(rows,6)
!  END MODULE Key

!!!   PROGRAM Main
!!!   USE Key
!!!! Input parameters
!!!   INTEGER,PARAMETER :: numT=4,n_pwrk=5000,iwr=1,ipl=1,ipr=2,nq=12 ! J.Cai set nq=12 old nq=10
!!!   DOUBLE PRECISION,PARAMETER :: P=1.d0,Tref=2000d0,Tmin=000d0,Tmax=1500d0
!!!   DOUBLE PRECISION,PARAMETER :: xmfr(3)=(/0.0d0,1.0d-3,0.d0/)
!!!   DOUBLE PRECISION,PARAMETER :: klmin=1.d-9,pwr=0.1d0
!!!   DOUBLE PRECISION,PARAMETER :: fvsoot=0.d-6,nsoot=1.89d0,ksoot=0.92d0
!!!!
!!!   INTEGER,PARAMETER :: nabs=(wvnm_e-wvnm_b)/wvnmst+100,ncdim=100
!!!   DOUBLE PRECISION,PARAMETER :: molemass(3)=(/44.d0,18.d0,16.d0/)
!!!   DOUBLE PRECISION,PARAMETER :: hhh=6.626076d-34,ccc=2.997925d10
!!!   DOUBLE PRECISION,PARAMETER :: kkk=1.380658d-23,nnn=6.022136d23
!!!   DOUBLE PRECISION,PARAMETER   :: PI=3.1415926d0,T0=296.d0,P0=1.d0
!!!   DOUBLE PRECISION,PARAMETER   :: c1=3.7419d-12,c2=1.4388d0,sigma=5.67d-12
!!!   DOUBLE PRECISION  :: hck,hckt,hckt0,hcktt0,bk,ck,sootf,sootfc,patm,cnchk
!!!   DOUBLE PRECISION  :: wvnm(nabs),absc(0:nabs),abscof,ri,mass,RR0,VV0
!!!   DOUBLE PRECISION  :: kmin=1.d10,kmax=1.d-10,pwrk_min,pwrk_max,pwrk_step,temp
!!!   INTEGER           :: i,it,icl,j,l,js,jj,lines,imax,ifg,im1,ip1,iq,ic(0:nq)
!!!   INTEGER           :: number,iadd,iaddo=0,bin1st,binlst,numbins,number2
!!!   DOUBLE PRECISION  :: k(0:n_pwrk),pwrk(0:n_pwrk),af(numT,n_pwrk),as(numT,n_pwrk)
!!!   DOUBLE PRECISION  :: ki,kim1,kpwri,kpwrim1,dtot,dtoto,dlcl,T(0:numT)
!!!   DOUBLE PRECISION  :: ff(0:numT,0:n_pwrk+1),gg(0:numT,0:n_pwrk+1)
!!!   DOUBLE PRECISION  :: V,V0,deta,dk
!!!   DOUBLE PRECISION  :: wvnm_b2,wvnm_e2,wvnmst2,kplancklbl,kplanckfsk,kplancksum,kperr
!!!   DOUBLE PRECISION  :: intensity,line_width,klmax,TT0,cr,crcn,c1sigt4(0:numT),eb(0:numT),sum(0:numT)
!!!   CHARACTER         :: dummy*25
!!!! J.Cai added filenames start
!!!   character(256),parameter :: kvsgFile='kvsgh2o-0p-2000K.dat'
!!!   character(256),parameter :: kvsgqFile='kvsgqh2o-0p-2000K.dat'
!!!   character(256),parameter :: abscFile='absch2o-0p-2000K.dat'
!!!   character(256),parameter :: CO2File='/opt/HITRAN2008/02_hit08.par'
!!!   character(256),parameter :: CH4File='/opt/HITRAN2008/06_hit08.par'
!!!   character(256),parameter :: H2OFile='/opt/HITRAN2008/01_hit08.par'
!!!   !character(256),parameter :: H2OFile='/opt/HITRAN/hitran96/by_molec/01_hit96.par'
!!!! J.Cai added filenames end
!!!! J.Cai comment out starts
!!!!! Selection of g-values for numerical quadrature, using a Numerical Recipes routine
!!!!! If Numerical Recipes is not available, set nq=12, comment out the following 8 lines of code,
!!!!! and uncomment the 5-line REAL declaration following it
!!!!   REAL              :: gqs(nq),wqs(nq),kq(nq),aq(numt,nq),gq(nq),wq(nq),gaujac,alf=3.,bet=-.9,wsum
!!!!! Get quadrature coefficients from Numerical recipies
!!!!    wsum=0.
!!!!    CALL GAUJAC(gqs,wqs,nq,alf,bet)
!!!!          do iq=1,nq
!!!!            gq(iq)=0.5*(1.-gqs(iq))
!!!!            wq(iq)=wqs(iq)/(2.**(alf+bet+1)*gq(iq)**alf*(1.-gq(iq))**bet)
!!!!            wsum=wsum+wq(iq)
!!!!          enddo
!!!!! Correction to make sum(wq)=1
!!!!         wq=wq/wsum
!!!!! End quadrature coefficients from Numerical recipies
!!!! J.Cai comment out ends
!!!! Selection of precalculated g-values for numerical quadrature, for nq=12,alf=3.,bet=0.
!!!! J.Cai uncomment out starts
!!!    REAL              :: kq(nq),aq(numt,nq), &
!!!           gq(nq)=(/ 5.120075E-02,1.170678E-01,2.015873E-01,3.007074E-01,4.095012E-01,5.225285E-01,  &
!!!                     6.341280E-01,7.387071E-01,8.310236E-01,9.064499E-01,9.612060E-01,9.925594E-01/),&
!!!           wq(nq)=(/ 5.556622E-02,7.576839E-02,9.258290E-02,1.048306E-01,1.118451E-01,1.132605E-01,  &
!!!                     1.090012E-01,9.927844E-02,8.457905E-02,6.563999E-02,4.341329E-02,1.904792E-02/)
!!!! J.Cai uncomment out starts
!!!! PRESSURE IN atm
!!!   patm=P/1.01325
!!!! IF P-BASED absc is requested make sure it is not a mixture
!!!    IF(ipl==1) THEN
!!!        cnchk=xmfr(1)*xmfr(2)+xmfr(1)*xmfr(3)+xmfr(3)*xmfr(2)
!!!        IF(cnchk > 1.d-6) PAUSE 'ERROR: you cannot use pressure-based absorption coefficient for mixture!'
!!!    ENDIF
!!!! Open output files
!!!! Output file for k,g(T,k),a(T,k) for all numT temperatures and all n_pwrk k-values
!!!   OPEN(7,FILE=kvsgFile) ! J.Cai centralize filename
!!!! Header formatted for TECPLOT, for a numT of 4
!!!   write(7,6)
!!! 6 format('VARIABLES = k,g0,g1,g2,g3,g4,a1,a2,a3,a4')
!!!! Output file for wq,gq(Tref),kq(gq),aq(T,kq) for all numT temperatures and all nq g-values
!!!   IF(ipr==2) THEN
!!!       OPEN(8,FILE=kvsgqFile,STATUS='unknown') ! J.Cai centralize filename
!!!! Header formatted for readability, for a numT of 4
!!!        write(8,8)
!!!   ENDIF
!!! 8 format('    wq',9x,'gq',9x,'kq',8x,'aq1',8x,'aq2',8x,'aq3',8x,'aq4')
!!!! File containing absorption coefficient
!!!   !IF(iwr>0) OPEN(9,FILE='c:\absco\absch2o-0p-2000K.dat',STATUS='unknown')
!!!   IF(iwr>0) OPEN(9,FILE=abscFile,STATUS='unknown') ! J.Cai centralize filename
!!!
!!!   hck=hhh*ccc/kkk
!!!   hckt0=hck/T0
!!!   IF(iwr==2) THEN
!!!! READ (GAS-ONLY) ABSORPTION COEFFICIENT FROM FILE
!!!      read(9,97) dummy
!!!      read(9,98) number
!!!      read(9,92) dummy,wvnm_b2,wvnm_e2,wvnmst2
!!!      number2=(wvnm_e2-wvnm_b2)/wvnmst2+1
!!!      IF(number2 /= number) PAUSE 'bad data file'
!!!      read(9,*) (absc(i),i=1,number)
!!!      CLOSE(9)
!!!      DO i=1,number
!!!         wvnm(i)=wvnm_b2+(i-1)*wvnmst2
!!!      ENDDO
!!!
!!!   ELSE
!!!
!!!! CALCULATE NECESSARY WAVENUMBERS AND ABSORPTION COEFFICIENTS
!!!! WVNM(I): WAVENUMBER
!!!! ABSC(I): ABSORPTION COEFFICIENT; INITIALIZE FOR SOOT
!!!   number=(wvnm_e-wvnm_b)/wvnmst+1
!!!   IF(number>nabs) PAUSE 'increase nabs'
!!!   DO i=1,number
!!!    wvnm(i)=wvnm_b+(i-1)*wvnmst
!!!    absc(i)=0.d0
!!!   ENDDO
!!!
!!!  kplancksum=0.d0
!!!! Scan over gases
!!!  DO ifg=1,3
!!!  IF(xmfr(ifg)<.9d-3) CYCLE
!!!! Open HITEMP database
!!!  IF(ifg==1) THEN
!!!     lu=11
!!!     write(*,*) 'Reading CO2 data'
!!!     OPEN(lu,FILE=CO2File, action='READ')  ! J.Cai change path
!!!  ELSEIF(ifg==2) THEN
!!!     lu=12
!!!     write(*,*) 'Reading H2O data'
!!!     OPEN(lu,FILE=H2OFile, action='READ')  ! J.Cai change path
!!!  ELSE
!!!     lu=13
!!!     write(*,*) 'Reading CH4 data'
!!!     OPEN(lu,FILE=CH4File, action='READ') ! J.Cai change path
!!!  ENDIF
!!!! GET LINE INFORMATION FROM HITEMP (all lines wvnm_b.le.lmbda.le.wvnm_e)
!!!   CALL Getdata(lines)
!!!   CLOSE(lu)
!!!   write(*,*) 'Gas ',ifg,',  lines read: ', lines
!!!   mass=1.d3*molemass(ifg)/nnn
!!!   hckt=hck/Tref
!!!   hcktt0=hckt0-hckt
!!!! MULTIPLIER FOR PRESSURE-BASED ABSORPTION COEFFICIENT
!!!   cr=nnn/8.314d1/Tref
!!!! MULTIPLIER FOR LINEAR ABSORPTION COEFFICIENT
!!!   IF(ipl==0) cr=cr*xmfr(ifg)*P
!!!   IF(ifg==1) THEN
!!!    write(*,*) 'start of absc loop for co2'
!!!! VIBRATIONAL PARTITION FUNCTION FOR CO2
!!!    V0=(1d0-DEXP(-666d0*hckt0))**2*(1d0-DEXP(-2396d0*hckt0))*   &
!!!            (1d0-DEXP(-1351d0*hckt0))
!!!    V=(1d0-DEXP(-666d0*hckt))**2*(1d0-DEXP(-2396d0*hckt))*   &
!!!            (1d0-DEXP(-1351d0*hckt))
!!!    VV0=V/V0
!!!! ROTATIONAL PARTITION FUNCTION
!!!    TT0=T0/Tref
!!!    RR0=TT0
!!!  ELSEIF(ifg==2) THEN
!!!   write(*,*) 'start of absc loop for h2o'
!!!! VIBRATIONAL PARTITION FUNCTION FOR H2O
!!!    V0=(1d0-DEXP(-3652d0*hckt0))*(1d0-DEXP(-1595d0*hckt0))*   &
!!!            (1d0-DEXP(-3756d0*hckt0))
!!!    V=(1d0-DEXP(-3652d0*hckt))*(1d0-DEXP(-1595d0*hckt))*   &
!!!            (1d0-DEXP(-3756d0*hckt))
!!!    VV0=V/V0
!!!! ROTATIONAL PARTITION FUNCTION
!!!    TT0=T0/Tref
!!!    RR0=TT0**1.5
!!!  ELSE
!!!   write(*,*) 'start of absc loop for ch4'
!!!! VIBRATIONAL PARTITION FUNCTION FOR CH4
!!!    V0=(1d0-DEXP(-1306d0*hckt0))**3*(1d0-DEXP(-1526d0*hckt0))**2*   &
!!!            (1d0-DEXP(-2914d0*hckt0))*(1d0-DEXP(-3020d0*hckt0))**3
!!!    V=(1d0-DEXP(-1306d0*hckt))**3*(1d0-DEXP(-1526d0*hckt))**2*   &
!!!            (1d0-DEXP(-2914d0*hckt))*(1d0-DEXP(-3020d0*hckt))**3
!!!    VV0=V/V0
!!!! ROTATIONAL PARTITION FUNCTION
!!!    TT0=T0/Tref
!!!    RR0=TT0**1.5
!!!  ENDIF
!!!   VR0=VV0*RR0
!!!
!!!! Scan over lines
!!!    c1sigt4(0)=c1/(sigma*Tref**4)
!!!    DO l=1,lines
!!!! print out every 10000 lines so that user sees machine is working
!!!    IF(l-(l/10000)*10000 == 0) write(*,*) 'reading line #',l
!!!      line_width=(Data(l,4)*xmfr(ifg)+Data(l,3)*(1.d0-xmfr(ifg)))*patm*TT0**Data(l,6)
!!!! molecule-based intensity
!!!      intensity=Data(l,2)*VR0*DEXP(hcktt0*Data(l,5))
!!!! pressure-based or linear intensity
!!!      intensity=intensity*cr*(1d0-DEXP(-data(l,1)*hckt))/(1d0-DEXP(-data(l,1)*hckt0))
!!!! calculate Planck-mean absorption coefficient to ensure adequacy of spectral resolution
!!!      kplancksum=kplancksum+intensity*c1sigt4(0)*data(l,1)**3/(EXP(hck/Tref*data(l,1))-1.d0)
!!!! absorption coefficient at line center
!!!      klmax=intensity/(PI*line_width)
!!!      IF(klmax<klmin) CYCLE
!!!! Find wavenumber (subscript) closest to line center
!!!      icl=(data(l,1)-wvnm_b)/wvnmst+1
!!!! Scan over adjacent wavenumbers to see whether line makes contribution to absco
!!!        DO i=icl,number
!!!          deta=(data(l,1)-wvnm(i))/line_width
!!!          dk=klmax/(1.d0+deta*deta)
!!!          absc(i)=absc(i)+dk
!!!          IF(dk<klmin) EXIT
!!!        ENDDO
!!!        DO i=icl-1,1,-1
!!!          deta=(data(l,1)-wvnm(i))/line_width
!!!          dk=klmax/(1.d0+deta*deta)
!!!          absc(i)=absc(i)+dk
!!!          IF(dk<klmin) EXIT
!!!        ENDDO
!!!    ENDDO ! l
!!!  ENDDO ! ifg
!!!   write(*,*) 'end of absc loop'
!!!
!!!! WRITE ABSC TO DATA FILE IF DESIRED
!!!  IF(iwr==1) THEN
!!!      write(9,96)
!!!      write(9,98) number
!!!      write(9,95) wvnm_b,wvnm_e,wvnmst
!!!      write(9,99) (absc(i),i=1,number)
!!!      CLOSE(9)
!!!      ENDIF
!!!92 FORMAT(a1,3f12.5)
!!!95 FORMAT('#',3f12.5)
!!!96 FORMAT('variables = "absco"')
!!!97 FORMAT(a25)
!!!98 FORMAT('zone i=',i8)
!!!99 FORMAT(e14.5)
!!!   ENDIF ! iwr
!!!
!!!! ADD SOOT CONTRIBUTION TO ABSCO (AFTER STORING/RETRIEVING GAS-ONLY VALUES)
!!!  sootfc=sootf(fvsoot,nsoot,ksoot)
!!!  DO i=1,number
!!!    absc(i)=absc(i)+sootfc*wvnm(i)
!!!  ENDDO
!!!
!!!! CALCULATE K-DISTRIBUTION
!!!! pwrK_MAX AND pwrK_MIN DEFINES THE MAXIMUM AND MINIMUM pwr(K) VALUES.
!!!! pwrK_STEP IS THE pwr(K) INTERVAL WE USE TO SCAN THE SPECTRUM.
!!!! pwrK(I): pwr(K) VALUES
!!!! FF(I): F*DELTAK
!!!! GG(I): G FUNCTION (CUMULATIVE K-DISTRIBUTION FUNCTION)
!!!! Find maximum kappa value
!!!    DO i=1,number
!!!      IF(absc(i)>kmax) then
!!!        kmax=absc(i)
!!!        imax=i
!!!      endif
!!!    ENDDO
!!!    IF(ipl==0) THEN
!!!        write(*,93) kmax,wvnm(imax)
!!!    ELSE
!!!        write(*,94) kmax,wvnm(imax)
!!!    ENDIF
!!!93 FORMAT(' kmax =',f12.5,'cm-1 at ',f10.4,'cm-1')
!!!94 FORMAT(' kmax =',f12.5,'bar-1cm-1 at ',f10.4,'cm-1')
!!!   IF(kmax>kdmax .and. kdmax>0.d0) &
!!!        write(*,*) 'WARNING!!! kdmax is less than maximum absorption coefficient!!'
!!!    kmin=max(1.d-9,kdmin)
!!!    kmax=max(kdmax,kmax)
!!!    pwrk_min=kmin**pwr
!!!    pwrk_max=kmax**pwr
!!!    pwrk_step=(pwrk_max-pwrk_min)/REAL(n_pwrk-1)
!!!
!!!!*****************************************************************************
!!!! K-DISTRIBUTION
!!!   write(*,*) 'start of k-dist loop'
!!!
!!!   DO i=0,n_pwrk
!!!    pwrk(i)=REAL(i-1)*pwrk_step+pwrk_min
!!!    k(i)=pwrk(i)**(1./pwr)
!!!   ENDDO
!!!
!!!   c1sigt4(0)=c1/(sigma*Tref**4)*wvnmst
!!!   T(0)=Tref
!!!   DO it=1,numT
!!!    T(it)=max(300d0,Tmin+(Tmax-Tmin)*(it-1.)/(numT-1.))
!!!    c1sigt4(it)=c1/(sigma*T(it)**4)*wvnmst
!!!   ENDDO
!!!    absc(0)=2.*absc(1)-absc(2)
!!!    kpwri=absc(0)**pwr
!!!    ki=absc(0)
!!!    iaddo=max(1,Int((kpwri-pwrk_min)/pwrk_step+1))
!!!    absc(number+1)=2.*absc(number)-absc(number-1)
!!!
!!!   ff=0.d0
!!!! SCAN OVER SPECTRUM
!!!    sum=0.
!!!   DO i=1,number
!!!     DO it=0,numT
!!!       eb(it)=c1sigt4(it)*wvnm(i)**3/(EXP(hck/T(it)*wvnm(i))-1.d0)
!!!     ENDDO
!!!     sum=sum+eb
!!!       kpwrim1=kpwri
!!!       kim1=ki
!!!       kpwri=absc(i)**pwr
!!!       ki=absc(i)
!!!       iadd=max(1,Int((kpwri-pwrk_min)/pwrk_step+1))
!!!       dtoto=0.d0
!!!       dk=ki-kim1
!!!! New k-cell iadd > iaddo: calculate local detas for each k-cell
!!!     IF(iadd>iaddo) THEN
!!!        DO j=iaddo,iadd-1
!!!           dtot=(k(j+1)-kim1)/dk
!!!           dlcl=dtot-dtoto
!!!           dtoto=dtot
!!!           DO it=0,numT
!!!              ff(it,j)=ff(it,j)+eb(it)*dlcl
!!!           ENDDO
!!!        ENDDO
!!!        DO it=0,numT
!!!           ff(it,iadd)=ff(it,iadd)+eb(it)*(1.d0-dtot)
!!!        ENDDO
!!!! New k-cell iadd < iaddo: calculate local detas for each k-cell
!!!     ELSEIF(iadd<iaddo) THEN
!!!        DO j=iaddo,iadd+1,-1
!!!           dtot=(k(j)-kim1)/dk
!!!           dlcl=dtot-dtoto
!!!           dtoto=dtot
!!!           DO it=0,numT
!!!              ff(it,j)=ff(it,j)+eb(it)*dlcl
!!!           ENDDO
!!!        ENDDO
!!!        DO it=0,numT
!!!           ff(it,iadd)=ff(it,iadd)+eb(it)*(1.d0-dtot)
!!!        ENDDO
!!!! New k-cell iadd = iaddo: entire deta remains in cell
!!!     ELSE
!!!     DO it=0,numT
!!!        ff(it,iadd)=ff(it,iadd)+eb(it)
!!!     ENDDO
!!!     ENDIF
!!!     iaddo=iadd
!!!   ENDDO ! number
!!!! FLAG EMPTY k-BINS
!!!    DO i=1,n_pwrk
!!!      bin1st=i
!!!      IF(ff(0,i)>0.d0) EXIT
!!!    ENDDO
!!!    DO i=n_pwrk,1,-1
!!!      binlst=i
!!!      IF(ff(0,i)>0.d0) EXIT
!!!    ENDDO
!!!
!!!! CORRECT K-VALUES BY HALF A BIN UPWARDS
!!!   DO i=1,n_pwrk-1
!!!    k(i)=.5D0*(k(i)+k(i+1))
!!!   ENDDO
!!!
!!!! CALCULATE G FUNCTIONS
!!!    DO it=0,numT
!!!      gg(it,binlst)=1.d0
!!!      DO i=binlst,bin1st,-1
!!!        gg(it,i-1)=gg(it,i)-ff(it,i)
!!!! CONVERT FROM F*DELTAK TO F IF DESIRED (comment out if not)
!!!!       ff(it,i)=ff(it,i)/(k(i)-k(i-1))
!!!      ENDDO
!!!! EXTRAPOLATE TO binlst+1 TO ALLOW AF-CALCULATION AT binlst
!!!    ff(it,binlst+1)=ff(it,binlst-1)-2*ff(it,binlst)
!!!      gg(it,binlst+1)=2.d0-gg(it,binlst-1)
!!!! IF (ALMOST) ENTIRE SPECTRUM IS CONSIDERED, SUM(eb)>1 DUE TO QUADRATURE ERROR --> g NEED ADJUSTMENT
!!!      IF(gg(it,bin1st-1)<0.d0) THEN
!!!        DO i=binlst,bin1st-1,-1
!!!          gg(it,i)=(gg(it,i)-gg(it,bin1st-1))/(1.d0-gg(it,bin1st-1))
!!!        ENDDO
!!!      ENDIF
!!!! SET GG(i<bin1st-1)=0 TO ACCOUNT FOR NEGLECTED PARTS OF SPECTRUM
!!!      DO i=0,bin1st-1
!!!        gg(it,i)=0.d0
!!!      ENDDO
!!!! SET GG(i>binlst)=1 FOR COMPATIBILIITTYY
!!!      DO i=binlst+1,n_pwrk
!!!        gg(it,i)=1.0000000000001d0
!!!      ENDDO
!!!!-- set first bin always to zero !! J.Cai
!!!	  gg(it,0) = 0.d0  !! J.Cai
!!!    ENDDO ! it
!!!
!!!! CALCULATE A-FUNCTIONS
!!!    numbins=binlst-bin1st+1
!!!    DO it=1,numT
!!!      DO i=bin1st,binlst
!!!        af(it,i)=(ff(it,i+1)+4*ff(it,i)+ff(it,i-1))/(ff(0,i+1)+4*ff(0,i)+ff(0,i-1))
!!!      ENDDO
!!!      DO i=1,bin1st
!!!        af(it,i)=af(it,bin1st)
!!!      ENDDO
!!!      DO i=binlst+1,n_pwrk
!!!        af(it,i)=af(it,binlst)
!!!      ENDDO
!!!    ENDDO   ! it
!!!! do a little smoothing on raw af
!!!      as=af
!!!      DO j=1,25
!!!        DO it=1,numT
!!!          DO i=bin1st,binlst
!!!            im1=i-1+1/i
!!!            ip1=i+1-i/n_pwrk
!!!            as(it,i)=(af(it,ip1)+4*af(it,i)+af(it,im1))/6.d0
!!!          ENDDO
!!!        ENDDO
!!!        af=as
!!!      ENDDO
!!!   write(*,*) 'end of k-dist loop'
!!!
!!!! OUTPUT
!!!   DO i=1,n_pwrk
!!!    WRITE(7,50) k(i),(gg(it,i),it=0,numT),(af(it,i),it=1,numT)
!!!   ENDDO
!!!50 FORMAT(1p10d16.8)
!!!   CLOSE(7)
!!!
!!!! Find k- and a-values for selected quadrature point g-values gq
!!!    IF(ipr == 2) THEN
!!!! k values for g(Tref)
!!!      it=0
!!!        i=bin1st-1
!!!        ic(0)=bin1st
!!!        DO iq=1,nq
!!! 11       IF(gg(it,i+1)>gq(iq)) THEN
!!!            kq(iq)=k(i)+(k(i+1)-k(i))*(gq(iq)-gg(it,i))/(gg(it,i+1)-gg(it,i))
!!!            ic(iq)=max(i,bin1st) !! J.Cai
!!!          ELSE
!!!            i=i+1
!!!            GOTO 11
!!!          ENDIF
!!!        ENDDO ! iq
!!!        DO iq=1,nq-1
!!!          ic(iq)=(ic(iq)+ic(iq+1))/2
!!!        ENDDO
!!!        ic(nq)=binlst
!!!! Find smoothened a for quadrature points
!!!      DO it=1,numT
!!!        DO iq=1,nq
!!!		  IF (ic(iq).NE.ic(iq-1)) THEN !! J.Cai, DO loop later assumes this condition
!!!            aq(it,iq)=0.d0
!!!            DO i=ic(iq-1),ic(iq)-1
!!!              aq(it,iq)=aq(it,iq)+.5d0*(af(it,i)+af(it,i+1))*(gg(it,i+1)-gg(it,i))
!!!            ENDDO
!!!            aq(it,iq)=aq(it,iq)/(gg(it,ic(iq))-gg(it,ic(iq-1)))
!!!		  ELSE  !! J.Cai
!!!		    aq(it,iq) = 1.d0  !! J.Cai
!!!		  ENDIF !! J.Cai
!!!        ENDDO
!!!      ENDDO   ! it
!!!      DO i=1,nq
!!!         WRITE(8,60) wq(i),gq(i),kq(i),(aq(it,i),it=1,numT)
!!!      ENDDO
!!!      CLOSE(8)
!!! 60   FORMAT(1p10e11.4)
!!!    ENDIF
!!!
!!!!
!!!! CALCULATE KPLANCK BY LBL AND FSK TO ASCERTAIN QUALITY OF K-DISTRIBUTION
!!!! LBL
!!!    kplancklbl=0.d0
!!!    DO i=1,number
!!!      kplancklbl=kplancklbl+absc(i)*c1sigt4(0)*wvnm(i)**3/(EXP(hck/T(0)*wvnm(i))-1.d0)
!!!    ENDDO
!!!! FSK
!!!    kplanckfsk=0.d0
!!!    DO i=bin1st,n_pwrk
!!!      kplanckfsk=kplanckfsk+0.5d0*(k(i)+k(i-1))*(gg(0,i)-gg(0,i-1))
!!!    ENDDO
!!!    IF(iwr<2) THEN
!!!        kperr=1d2*(1.d0-kplancklbl/kplancksum)
!!!        write(*,70) kplancksum,kplancklbl,kperr
!!!70      FORMAT(' k_Pl-sum =',e12.4,',  k_Pl-lbl =',e12.4,'; error =',f7.2,'%')
!!!        IF(kperr > .5d0) write(*,*) 'spectral resolution questionable: decrease wvnst!'
!!!    ENDIF
!!!    kperr=1d2*(1.d0-kplanckfsk/kplancklbl)
!!!    write(*,80) kplancklbl,kplanckfsk,kperr
!!!80  FORMAT(' k_Pl-lbl =',e12.4,',  k_Pl-fsk =',e12.4,'; error =',f7.2,'%')
!!!    IF(kperr > .5d0) write(*,*) 'k-distribution questionable: increase n_pwrk or change pwr!'
!!!
!!!   END PROGRAM Main

!!!************************************************************************
!!   SUBROUTINE Getdata(i)
!!!   USE Key
!!   IMPLICIT NONE
!!   INTEGER   :: i,isotp
!!   INTEGER   :: dummy1
!!   CHARACTER :: dummy2*10
!!! data(i,1) = wavenumber
!!! data(i,2) = intensity
!!! data(i,3) = b_air
!!! data(i,4) = b_self
!!! data(i,5) = E''
!!! data(i,6) = exponent for b
!!
!!   i=1
!!   DO
!!   READ(lu,FMT=150,END=4) dummy1,isotp,data(i,1),data(i,2),dummy2,data(i,3),    &
!!   READ(lu,FMT=150,END=4) dummy1,isotp,data(i,1),data(i,2),dummy2,data(i,3),    &
!!    data(i,4),data(i,5),data(i,6)
!!   If(data(i,1)<wvnm_b) CYCLE
!!   If(data(i,1)>wvnm_e) GOTO 4
!!! Ignore if not the most abundent isotope
!!   If(isotp/=1) CYCLE
!!   If(data(i,5)<0.d0) THEN
!!    data(i,5)=0.d0
!!    write(*,*) 'E"<0 at ', data(i,1)
!!    Endif
!!   If(data(i,4)<1.e-5) data(i,4)=5.d0*data(i,3)
!!   i=i+1
!!   ENDDO
!!150 FORMAT(i2,i1,f12.6,es10.3,a10,f5.4,f5.3,f10.4,f4.2)
!!! The above format specification is for HITRAN 2008
!!! For HITRAN 96 use the following
!!!150 FORMAT(i3,f12.6,e10.3,a10,f5.4,f5.4,f10.4,f4.2)
!!  4 i=i-1
!!  5 RETURN
!!   END SUBROUTINE Getdata
!************************************************************************
    FUNCTION sootf(fv,n,k)
    DOUBLE PRECISION    :: fv,n,k,sootf,pi36=1.130973355d2,nsq,ksq
    nsq=n*n
    ksq=k*k
    sootf=pi36*n*k*fv/((nsq-ksq+2.d0)**2+4.d0*nsq*ksq)
    RETURN
    END
!************************************************************************

END MODULE rad_mc_calabs
