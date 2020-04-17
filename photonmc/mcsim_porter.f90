program main

integer, parameter:: nsize = 50
double precision, dimension(nsize,nsize,nsize):: ak,as,aqdot
integer i,j,k,npasses
npasses = 1

ak = 0.0d0    ! initialize absorption co-efficient  
as = 0.0d0    ! intialize  scattering co-effcient 
aqdot = 0.0d0 ! initialize heat flux 
call PMCSimulation(ak, as, aqdot, nsize)

do i=1,nsize
       do j=1,nsize 
         do k=1,nsize
           write(900,*) i, j, k , ak(i,j,k), as(i,j,k),aqdot(i,j,k)
         end do 
       end do 
end do 

end  program main


subroutine PMCSimulation(abc, aem,aq,nn) 
!Total number of photon packets 
!integer, parameter :: nphotons = 10
integer :: nphotons
double precision, parameter :: sigma = 5.67d-8
integer, parameter :: pi = 22.0d0/7.0d0
double precision, dimension(nn,nn,nn):: abc, aem, aq
double precision, dimension(nn,nn,nn):: tt, xco2, xh2o
double precision:: sigma_a, sigma_s, sigma_t 
integer, dimension (nn,nn,nn) :: bcid 

double precision:: xm, ym, zm, gsa, theta, phi,costheta, absijk,abstp,waven,eta,eta2
double precision:: rho,pres
double precision:: xi,yi,zi,xp, yp, zp 
double precision:: xl, yl, zl, pl 
double precision:: dcx, dcy, dcz,wt,dw
double precision:: r1, r2, r3, r4  
integer         :: nx, ny,nz, icount  
double precision:: dx, dy,dz, tin, dvol  
double precision, allocatable :: xloc(:),yloc(:),zloc(:)
integer         :: npart,npht
double precision:: dpart,cq, gasvolFrac
double precision:: getCosTheta
double precision:: Absorp_Coeff, em_surf
double precision :: rad,tc, ppr
integer :: Nrays, Nrays_tot
double precision:: ray_energy,PPR_crit_perc,E_abs
integer :: ix,iy,iz 

PPR_crit_perc = 0.1d0
em_surf = 1.0d0 
gsa = 0.0d0 !0.75d0            ! scattering anisotropy (g)  [-1, 1]  ! do not know the correct value 
tin = 700.0d0           ! To. do not know the To value.
xm = 0.0d0              ! x_min 
ym = 0.0d0              ! y_min 
zm = 0.0d0              ! z_min 
xl = 2.0d0              ! domain length along x-direction  
yl = 4.0d0              ! domain length along y-direction 
zl = 2.0d0              ! domain length along z-direction 

!   mesh parameters 
nx = nn                 ! no of nodes in x  
ny = nn                 ! no of nodes in y 
nz = nn                 ! no of nodes in z 
dx = xl/float(nx-1)     ! element size in x-direction      
dy = yl/float(ny-1)     ! element size in y-direction 
dz = zl/float(nz-1)     ! element size in z-direction 

allocate(xloc(nx))
allocate(yloc(nx))
allocate(zloc(nx))

do i=1,nx 
  do j=1,ny 
    do k=1,nz 
       xloc(i) =  xm + (i-1)*dx 
       yloc(j) =  ym + (j-1)*dy  !
       zloc(k) =  zm + (k-1)*dz
       rad = sqrt((xloc(i)-1d0)**2 + (yloc(j)-1d0)**2)
        if (rad < 1d0)then
            if (z_c <= 0.375d0)then
                tc = 400d0+1400d0*(zloc(k)/0.375d0)
            else
                tc = 1800d0-1000d0*(zloc(k)-0.375d0)/3.625d0
            end if
            tt(i,j,k) = (tc-800d0)*(1d0-3d0*rad**2+2d0*rad**3)+800d0
        else
            tt(i,j,k) = 800d0
        end if
       xco2(i,j,k) = 0.21 ! mole fraction of CO2
       xh2o(i,j,k) = 0.00 ! mole fraction of H2O
    end do 
  end do 
end do 

bcid = 0 
! front and back boundaries (xy plane) 
do i=1,nx 
  do j=1,ny
     tt(i,j,1)  = 400.0d0 
     tt(i,j,nz) = 800.0d0 
     bcid(i,j,1)  = 1 
     bcid(i,j,nz) = 2
  end do 
end do 

! left and right boundaries  (yz plane)
do j=1,ny 
 do k=1,nz 
  tt(1,j,k)  = 300.0d0 
  tt(nx,j,k) = 300.0d0 
  bcid(1,j,k)  = 3
  bcid(nx,j,k) = 4
 end do 
end do 

! top and bottom boundaries  (xz plane ) 
do i=1,nx 
 do k=1,nz 
  tt(i,1,k)   = 300.0d0 
  tt(i,ny,k)  = 300.0d0 
  bcid(1,j,k)  = 5 
  bcid(i,ny,k) = 6
 end do 
end do 

Nrays_tot = 0
do i=1,nx 
  do j=1,ny 
    do k=1,nz
     pres  = 1.0d0  !is this the correct value ?? 
     gasvolFrac = 1.0d0 
     npart = 10     !donot know the correct value. please chage this variable 
     dpart = 1.0d0  !do not know the correct value. please change this variable 
     rho   = 1.0d0  !do not know the correct value. please change this variable 
 !   call plank_mean(xco2(i,j,k), xh2o(i,j,k), pres,gasvolFrac, tt(i,j,k), absijk)
!     call plank_mean(0.21d0,       0.00d0,       1.0,1.0d0, 700.0, absijk)
!************************************************************************
! not sure about this loop. Please correct this
!      cq = (1.0d0/float(nphotons))*absijk*5.67d-8*tt(i,j,k)**4.0d0        !

      abstp     = 0.1d0 
      ppr       = abstp * 2d0 
      cq        = 4.0*abstp*sigma*tt(i,j,k)**4.0
      dvol      = dx*dy*dz  
      nphotons  = nint(cq*dvol/ppr)
      Nrays_tot = Nrays_tot + Nrays

!       xi  = xloc(i) + dx*rand() 
!       yi  = yloc(j) + dy*rand()
!       zi  = zloc(k) + dz*rand()
!       dcx = 0.0d0                  ! direction cosines
!       dcy = 0.0d0                  ! direction cosines
!       dcz = 1.0d0                  ! direction cosines
!       icount = 0
!       npht = 1 
      do icount 1, nphotons           ! photons loop 
       r1 = rand()
       xi  = xloc(i) + dx*rand() 
       yi  = yloc(j) + dy*rand()
       zi  = zloc(k) + dz*rand()
       theta = acos(1d0-2d0*rand())              !2*pi*r1    
!       if(cos(theta)>0.0d0) then 
          waven =  400.0d0 + (1.0d0/(10000.0-400.0))*r1       ! calculate wave number for a given rand()
          if (abs(xi) > xl)then
              write(*,*)"Invalid x Coordinate:", xi
              cycle
          elseif (abs(yi) > yl)then
              write(*,*) "Invalid y Coordinate:", yi
              cycle
          elseif (zi > zl .OR. zi < 0d0)then
              write(*,*)"Invalid z Coordinate:", zi
              cycle 
          end if 
          phi   =  2.0*pi*rand()
          call edge_length(xi,yi,zi,dx,dy,dx,nx,ny,nz,theta,phi,pl)
          dcx = cos(phi)*sin(theta)
          dcy = sin(phi)*sin(theta)
          dcz = cos(theta)
          ray_energy = PPR 
          do while(ray_energy > PPR*PPR_crit_perc/100d0) 
          xe = xi + sl *  dcx
          ye = yi + sl *  dcy
          ze = zi + sl *  dcz
          ix = floor(xe/dx) 
          iy = floor(ye/dy) 
          iz = floor(ze/dz) 
          !Absorb Fraction of Ray Energy to Cell
          abstp = 0.1d0 
          E_abs = (1d0-exp(-abstp*sl))*ray_energy
          ray_energy = ray_energy - E_abs
          abc(ix,iy,iz) = abc(ix,iy,iz) + E_abs

          !Check if Ray Hit a Surface
          if (abs(abs(xe)-xl) <= eps .OR. & ! xe <=eps
              abs(abs(ye)-yl) <= eps .OR. &
              abs(ze-zl) <= eps .OR. &
              abs(ze) <= eps)then

              !Absorb Fraction of Ray Energy to Surf
!              surfID = GetsurfID(x_ray,y_ray,z_ray)
              E_abs = em_surf*ray_energy
              aem(ix,iy,iz) = aem(ix,iy,iz) + E_abs
              ray_energy = ray_energy - E_abs
                !Reflect Diffusely (if Wall Emissivity /= 1)
                 if (em_surf == 1d0)then
                     EXIT 
                 end if
         !        call random_number(rand)
         !        call random_number(rand2)

         !        !faceID = surf_data%face(surfID,1)
         !        reflected_direction = Reflect(faceID,rand,rand2)
         !        theta = reflected_direction(1)
         !        phi = reflected_direction(2)
             end if

             !Find New Distance to Next Cell Edge in Current Direction
!             D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
             call edge_length(xe,ye,ze,dx,dy,dx,nx,ny,nz,theta,phi,pl) 
             !Find New Cell ID
             x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
             y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
             z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
             ix = x_mid/dx
             iy = y_mid/dy
             iz = z_mid/dz
!             cellID = GetcellID(x_mid,y_mid,z_mid)
         end do
 !Absorb Final Small Fraction of Ray Energy if Below Threshold
            aem(ix,iy,iz) = aem(ix,iy,iz) + ray_energy
            ray_energy = 0d0
        end do
    end do

do i=1,nx,nx-1
  do j=1,ny,ny-1
   do k=1,nz,nz-1 
      abstp     = 0.1d0 
      ppr       = abstp * 2d0  
      cq        = 4.0*abstp*sigma*tt(i,j,k)**4.0
      dvol      = dx*dy*dz  
      nphotons  = nint(cq*dvol/ppr)
      Nrays_tot = Nrays_tot + Nrays


               r2 = rand()
               write(*,*) abstp, r2 
                if(abstp>r2) then 
                  npht= npht +1
                  aa(ix,iy,iz) = npht
                  write(*,*) 
                end if 
             end if
          endif
       end if  !cos theta loop  
       call spin(dcx, dcy, dcz, gsa) ! compute new photon direction
!      dcz =  sin(eta2)*cos(theta)
! not sure about this loop. Please correct this 
!************************************************************************
      end do   !photon loop 
10 continue
    end do 
  end do 
end do 

deallocate(xloc)
deallocate(yloc) 
deallocate(zloc) 
return
end

!Combined Absorption Coefficient Function---------------------------------
    function Absorp_Coeff(eta,rho,T,P_abs,N_part,d_part) result(abs_coef)
!This function returns the combined absorption coefficient of a CO2 gas and carbon particle mixture (m^-1), using the Elsasser narrow-band model (Edwards 1967), given the current wave number (cm^-1), density (gm/m^3), temperature (K), partial pressure of CO2 (atm), number of carbon particles (m^-3), and diameter of particles (m). 
        implicit none
        !Initialize Variables
        integer, intent(in) :: N_part
        integer :: i
        double precision, parameter :: pi = 22.0/7.0
        double precision, intent (in) :: eta, rho, T, P_abs, d_part
        double precision, parameter :: h = 6.6237D-27, k = 1.3802D-16, &
            c = 2.9979D10
        double precision :: abs_coef, nu1 = 1351D0, nu3 = 2396D0, &
            phi1, phi2, phi3, T_ref = 100D0, P_ref = 1D0, Pe, delta, &
            beta, Sc_d
        double precision, dimension(5) :: &
            eta_c = (/667D0,960D0,1060D0, 2410D0,3660D0/), & 
            b = (/1.3D0,1.3D0,1.3D0,1.3D0,1.3D0/), &
            n = (/0.7D0,0.8D0,0.8D0,0.8D0,0.65D0/), &
            a = (/2D0,2D0,2D0,1D0,2D0/), C1, C2, C3, C3_ref

        !Calculate Phi(T) Functions
        phi1 = ((1D0-exp(-(h*c/k/T)*(nu3-nu1)))*(exp(-h*c*nu1/k/T)-0.5D0* &
            exp(-2D0*h*c*nu1/k/T)))/((1D0-exp(-h*c*nu1/k/T))*&
            (1D0-exp(-h*c*nu3/k/T)))
        phi2 = ((1D0-exp(-h*c*(nu1+nu3)/k/T)))/((1D0-exp(-h*c*nu1/k/T))* &
            (1D0-exp(-h*c*nu3/k/T)))
        phi3 = (1D0+0.053D0*(T/T_ref)**(1.5D0))

        !Calculate Constants
        C1 = (/19D0,0.76D0*phi1,0.76D0*phi1,110D0,4D0*phi2/)
        C2 = (/6.9D0 *(T/T_ref)**0.5D0,1.6D0*(T/T_ref)**0.5D0*C1(2)**0.5D0, & 
            1.6D0*(T/T_ref)**0.5D0*C1(3)**0.5D0,31D0*(T/T_ref)**0.5D0, &
            8.6D0*phi3/)
        C3 = (/12.9D0*(T/T_ref)**0.5D0,12.4D0*(T/T_ref)**0.5D0, &
            12.4D0*(T/T_ref)**0.5D0, 11.5D0*(T/T_ref)**0.5D0, &
            24D0*(T/T_ref)**0.5D0/)
        C3_ref = (/12.9D0*(T_ref/T_ref)**0.5D0,12.4D0*(T_ref/T_ref)**0.5D0, &
            12.4D0*(T_ref/T_ref)**0.5D0, 11.5D0*(T_ref/T_ref)**0.5D0,24D0* &
            (T_ref/T_ref)**0.5D0/)

        !Loop Over All Major Spectral Absorption Bands
        abs_coef = 0D0
        do i = 1,size(eta_c)
            Pe = (((1D0-P_abs)+b(i)*P_abs)/P_ref)**n(i)
            delta = 30D0*C3_ref(i)
            beta = C2(i)**2D0*Pe/4D0/C1(i)/C3(i)
            Sc_d = C1(i)/C3(i)*exp(-(a(i)/C3(i))*abs(eta-eta_c(i)))
            abs_coef = abs_coef+rho*Sc_d*sinh(PI*beta/2D0)/ &
                (cosh(PI*beta/2D0)-cos(2*PI*(eta-eta_c(i))/delta))
        end do

        !Add Contribution from Carbon Particle Absorption (Howell Eq. 15.3)
        !CURRENTLY ASSUMING ABSORPTION EFFICIENCY FACTOR = 1 FOR ALL ETA
!        abs_coef = abs_coef + N_part*PI*(d_part**2D0)/4D0
    end function Absorp_Coeff


!Compute the new photon direction (due to scattering event) 
subroutine spin(mu_x, mu_y, mu_z, g) 
double precision, parameter :: PI = 22.0d0/7.0d0 
double precision :: getCosTheta
double precision :: costheta 
double precision :: phi
double precision :: sintheta 
double precision :: sinphi 
double precision :: cosphi
double precision :: denom, muzcosphi, ux, uy, uz
double precision :: g  
double precision :: mu_x, mu_y, mu_z

costheta = getCosTheta(g)  !scattering angle, θ | Henyey-Greenstein phase function 
phi = 2.0d0 * PI * rand()      ! polar angle 
sintheta = sqrt(1.0d0 - costheta * costheta) ! sin(theta)
sinphi = sin(phi)
cosphi = cos(phi)

if (mu_z .eq. 1.0d0) then
    mu_x = sintheta * cosphi  
    mu_y = sintheta * sinphi
    mu_z = costheta
else if (mu_z .eq. -1.0d0) then  
    mu_x = sintheta * cosphi 
    mu_y = -sintheta * sinphi
    mu_z = -costheta
else 
    denom = sqrt(1.0d0 - mu_z * mu_z)
    muzcosphi = mu_z * cosphi 
    ux = sintheta * (mu_x * muzcosphi - mu_y * sinphi) / denom + mu_x * costheta
    uy = sintheta * (mu_y * muzcosphi + mu_x * sinphi) / denom + mu_y * costheta 
    uz = -denom * sintheta * cosphi + mu_z * costheta
    mu_x = ux
    mu_y = uy
    mu_z = uz
end if 

return 
end  

function getCosTheta(g) ! sampling the H-G scattering phase function 
double precision :: getCosTheta
double precision :: g 
double precision mu
!    integer,parameter :: seed = 86456
!    call srand(seed)
mu = (1.0 - g * g) / (1.0 - g + 2.0 * g * rand())     ! part3 in cosTheta 
if (g.eq.0.0d0) then 
      getCosTheta = 2.0 * rand()  - 1.0d0
else 
     getCosTheta = (1.0 + g * g - mu * mu) / (2.0 * g)
end if 
end  

subroutine plank_mean(gasCO2, gasH2O, gasp,  gasvolFrac, temp, kp)
implicit none
integer, parameter :: kpLen=126
double precision, dimension(kpLen) :: kpT, kpCO2, kpH2O
double precision, intent(in) :: gasCO2, gasH2O, gasp,  gasvolFrac, temp
double precision :: kco2, kh2o
double precision, intent(out) :: kp
double precision:: interp1ess

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

kco2 = interp1ess(kpLen,kpT, kpCO2, temp)
kh2o = interp1ess(kpLen,kpT, kpH2O, temp)
    !---M.S-rad: the tabulated absorption coefficients are in cm^-1 atm^-1 and P is in bar
kp = (kco2*gasCO2+kh2o*gasH2O) * (gasp* 0.986923) * gasvolFrac
return 
end 


function interp1ess(nxy ,xx,yy, xi) result (yi)
! given function xx(nxy)<=> yy(nxy), find yi(ni) corresponding to xi(ni)
! ONLY FOR EQUALLY SPACED GRID INTERPOLATION
implicit none
integer,intent(in) :: nxy
double precision,dimension(nxy),intent(in) :: xx, yy
double precision,intent(in) :: xi
double precision :: yi
 ! local variables
integer :: iq,i
double precision :: wx,dx

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

subroutine edge_length(xc,yc,zc,dx,dy,dx,nx,ny,nz,thet,phia,sl)
double precision, parameter ::  eps=1.0e-12
double precision :: xc, yc, zc, dx, dy, dz, sl 
double precision :: xl, yl, zl,xd,yd,zd
integer :: nx, ny, nz  
double precision :: thet, phia 
double precision :: dcx, dcy, dcz 
double precision, dimension(6) :: t = -1
double precision :: x_err,y_err,z_err
integer :: i, cellID

dcx = sin(thet)*cos(phia) 
dcy = sin(thet)*sin(phia) 
dcz = cos(thet) 

xl = 2.0d0
yl = 2.0d0
zl = 4.0d0

xd = xc 
yd = yc
zd = zc 
!Adjust Coordinates if Point is on Positive Surface
if (abs(xd- xl) < eps)then
    xd = xd - dx/2d0
end if
if (abs(yd - yl) < eps)then
    yd = yd - dy/2d0
end if    
if (abs(zd - zl) < eps)then
    zd = zd - dz/2d0
end if

x0 = dx*floor(xd/dx)
y0 = dy*floor(yd/dy)
z0 = dz*floor(zd/dz) 

!Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
x_err = abs(mod(xc+mod(nx,2)*dx/2d0,dx))
y_err = abs(mod(yc+mod(ny,2)*dy/2d0,dy))
z_err = abs(mod(zc,dz))
if (x_err <= eps .OR. dx-x_err <= eps)then
    if (dcx > 0)then
        x0 = xc
    else
        x0 = xc - dx
    end if
end if
if (y_err <= eps .OR. dy-y_err <= eps)then
    if (dcy > 0)then
        y0 = yc
    else
        y0 = yc - dy
    endif
end if
if (z_err <= eps .OR. dz-z_err <= eps)then
    if (dcz > 0)then
        z0 = zc
    else
        z0 = zc - dz
    end if
end if

!Find All t's
if (abs(dcx) >= 1d-15)then 
    t(1) = (x0-xc)/dcx
    t(2) = (x0+dx-xc)/dcx
end if
if (abs(dcy) >= 1d-15)then
    t(3) = (y0-yc)/dcy
    t(4) = (y0+dy-yc)/dcy
end if
if (abs(dcz) >= 1d-15)then          
    t(5) = (z0-zc)/dcz
    t(6) = (z0+dz-zc)/dcz
end if
!Find Distance to Edge (Smallest, Positive t)
sl = 1d15
do i = 1,6
    if (t(i) > 0d0 .and. t(i) < sl)then
        sl = t(i)
    end if
end do
end 

!    function Reflect(faceID,rand,rand2,dr1,dr2) result(reflected_direction)
!        !This function provides the new ray direction after a reflecting event occurs (reflected_direction = [theta',phi']), relative to the global cartesian coordinates, given the faceID at which the reflection occurs(-x=1, +x=2, -y=3, +y=4, -z=5, +z=6) and two random numbers seeded between 0 and 1 (rand and rand2). The reflecting surface is assumed to be diffuse. (This function can also be utilized when emmitting rays from surfaces not normal to the global coordinate system.)
!        implicit none
!        !Initialize Variables
!        integer, intent(in) :: faceID
!        double precision, intent(in) :: rand, rand2
!        double precision :: theta, phi, theta_p, phi_p, &
!        x_p, y_p, z_p, x_p_t, y_p_t, z_p_t
!        double precision, dimension(2) :: reflected_direction
!        double precision, dimension(6,2) :: face_angles
!        double precision, dimension(3,3) :: Q
!
!        !Define Angles (Theta, Phi) of Face Normal Vectors
!        face_angles(1,1) = PI/2D0
!        face_angles(1,2) = 0D0
!        face_angles(2,1) = PI/2D0
!        face_angles(2,2) = PI
!        face_angles(3,1) = PI/2D0
!        face_angles(3,2) = PI/2D0
!        face_angles(4,1) = PI/2D0
!        face_angles(4,2) = 3D0*PI/2D0
!        face_angles(5,1) = 0D0
!        face_angles(5,2) = 0D0
!        face_angles(6,1) = PI
!        face_angles(6,2) = 0D0
!        !Get Angles of Current Face Normal Vector
!        theta = face_angles(faceID,1)
!        phi = face_angles(faceID,2)
!
!        !Calculate Reflected Direction (Angles) Relative to Face Normal
!        theta_p = asin(sqrt(rand))
!        phi_p = 2D0*PI*rand2
!
!        !Calculate Reflected Direction (Vector) Relative to Face Normal
!        x_p = sin(theta_p)*cos(phi_p)
!        y_p = sin(theta_p)*sin(phi_p)
!        z_p = cos(theta_p)
!
!        !Define the Transformation Matrix
!        Q(1,1) = cos(phi)*cos(theta)
!        Q(1,2) = -sin(phi)
!        Q(1,3) = sin(theta)*cos(phi)
!        Q(2,1) = sin(phi)*cos(theta)
!        Q(2,2) = cos(phi)
!        Q(2,3) = sin(theta)*sin(phi)
!        Q(3,1) = -sin(theta)
!        Q(3,2) = 0
!        Q(3,3) = cos(theta)
!        !Transform Reflected Direction (Vector) to Global Coordinate System
!        x_p_t = x_p*Q(1,1)+y_p*Q(1,2)+z_p*Q(1,3)
!        y_p_t = x_p*Q(2,1)+y_p*Q(2,2)+z_p*Q(2,3)
!        z_p_t = x_p*Q(3,1)+y_p*Q(3,2)+z_p*Q(3,3)
!
!        !Back Calculate Reflected Direction (Angles) Relative to GCS
!        reflected_direction(1) = atan2(sqrt(x_p_t**2+y_p_t**2),z_p_t)
!        reflected_direction(2) = atan2(y_p_t,x_p_t)
!    end function Reflect

