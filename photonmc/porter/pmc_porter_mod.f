module rad_pmc_mod 
!use rad_gas_abs
!      use SET_GEOMETRY_PMC_mod, only :set_geometry_pmc
! Fluid grid cell dimensions and mesh size
USE geometry, only: DX, DY, DZ 
! Number of particles in the I/J/K direction
use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
USE param, only: dimension_3, dimension_m, north, south, east, west, top, bottom
USE geometry, only: IMIN2, IMAX2
!      USE param1, only: zero, one, half, small_number, undefined

use indices, only: i_of, j_of, k_of
! Global Parameters:
!---------------------------------------------------------------------//
use param1, only: ZERO

! Module procedures.
!---------------------------------------------------------------------//
use mpi_utility
use error_manager
USE geometry, only: imax, jmax, kmax
USE geometry, only: imax2, jmax2, kmax2
USE geometry, only: imax3, jmax3, kmax3, ijkmax3
USE geometry, only: x_max, x_min, y_max, y_min, z_max, z_min
! Fluid grid cell dimensions and mesh size
      USE functions, only: funijk
      use geometry, only: VOL, AYZ, AXZ, AXY
!     implicit none
!     integer, parameter :: imax = 17, jmax = 17, kmax = 34 
!     double precision, parameter :: x_max= 2d0, y_max= 2d0, z_max= 4d0
     double precision, dimension (:), allocatable  :: XEF  !(0:DIMENSION_I)
     double precision, dimension (:), allocatable  :: YNF  !(0:DIMENSION_J)
     double precision, dimension (:), allocatable  :: ZTF  !(0:DIMENSION_K)

!      use set_geometry_pmc_mod, only:XEC,YNC,ZTC  
!      double precision, dimension (:,:), allocatable :: coords 
      double precision, dimension (:), allocatable :: T_g 
      double precision, dimension (:), allocatable :: emise,abse,srad 
      double precision, parameter :: eps = 1d-12
      double precision, parameter :: PI = 4.d0*DATAN(1.d0)
      double precision, parameter :: stef_boltz = 5.670374d-8
      double precision :: dx1, dy1, dz1 
contains 

    subroutine rad_pmc_calc
!    USE geometry, only: imax, jmax, kmax
!    USE geometry, only: imax2, jmax2, kmax2
!    USE geometry, only: imax3, jmax3, kmax3, ijkmax3
!    USE geometry, only: x_max, x_min, y_max, y_min, z_max, z_min 
!!    USE geometry_pmc, only: XE, YN, ZT 
!! The East/North/Top face location of a given I/J/K index.
!!      use discretelement, only: XE, YN, ZT
!! Flag for 2D simulations.
      use geometry, only: NO_K
! The start and end indices of IJK loops
      use compar, only: IJKStart3, IJKEnd3
! The Upper and Loper indices covered by the current process.
      use compar, only: ISTART3, IEND3
      use compar, only: JSTART3, JEND3
      use compar, only: KSTART3, KEND3
! Fixed array sizes in the I/J/K direction
      use bc, only: bc_jj_ps, bc_defined, bc_type_enum, par_slip_wall
      use bc, only: bc_k_b, bc_k_t, bc_j_s, bc_j_n, bc_i_w, bc_i_e
      use bc, only: free_slip_wall, no_slip_wall, dimension_bc
      implicit none    
      double precision, parameter :: em_surf = 1d0
      double precision, parameter :: abs_coef = 0.1d0

    !Initialize Working Variables-----------------------------------------------
!     integer :: N_cells = imax*jmax*kmax
!     integer :: N_surfs = 2*(imax*jmax+imax*kmax+jmax*kmax)
      integer :: i,j,k,l,m,n, ijk, cellID, wlbcid !, faceID
      integer :: ip,jp,kp 
      integer :: nrays, nrays_tot,ncountr 
      double precision, parameter :: PPR_crit_perc = 0.1d0
      double precision :: PPR, pl 
      double precision :: dv, x_c, y_c, z_c, x_ray, y_ray, z_ray
      double precision :: Tc, rad, rand, rand2, theta, phi, ray_energy, x_mid, y_mid, z_mid
      double precision :: D_edge, E_abs, E_in_sum = 0d0, E_out_sum = 0d0, time1, time2
      integer :: inew, jnew,knew, ijkn
!    double precision, dimension (N_cells,3) :: coords = 0
!        type bc_type
!        integer,          dimension (N_surfs,1) :: node = 0
!        integer,          dimension (N_surfs,1) :: face = 0
!        double precision, dimension (N_surfs,1) :: glbnode = 0
!        double precision, dimension (N_surfs,3) :: gy = 0
!    end type
!    type(bc_type) :: bc

    double precision :: pres, xH2O, xCO2, volfrac,d_part, rho,eta 
    integer :: npart_cell
!    double precision, dimension (6) :: da
   ! double precision, dimension (2) :: reflected_direction
    character(len=1) :: c, bsp = char(8)
!    character(len=50) :: filename, grid_t, abs_coef_t, PPR_t
    integer :: IER 
    integer :: I1, I2, J1, J2, K1, K2
!    double precision, dimension (N_cells) :: T_g
!    double precision, dimension (N_cells,3) :: coords = 0
!    type cell_data_type
!        double precision, dimension (N_cells,3) :: gy = 0
!    end type
!    type(cell_data_type) :: cell_data
!    call CPU_TIME(time1)

!    nrays_tot = 10000000 
!    allocate(coords(dimension_3,3))
    allocate(T_g(dimension_3))
    allocate(emise(dimension_3))
    allocate(abse(dimension_3))
    allocate(srad(dimension_3))
    allocate( XEF (0:DIMENSION_I), STAT=IER )
    allocate( YNF (0:DIMENSION_J), STAT=IER )
    allocate( ZTF (0:DIMENSION_K), STAT=IER )

    !Calculate Size of Homogeneous Subregions---------------------------------- 
    dx1 = x_max/float(imax)
    dy1 = y_max/float(jmax)
    dz1 = z_max/float(kmax)
!    dv = dx1*dy1*dz1
!    da = (/dy1*dz1,dy1*dz1,dx1*dz1,dx1*dz1,dx1*dy1,dx1*dy1/)
    !write(*,*) imax,  jmax,kmax
!    write(*,*) imin,  jmin,kmin 
    !stop 
    !write(*,*) imin2, imax2,jmin2, jmax2, kmin2, kmax2
!    write(*,*) imin3, imax3,jmin3, jmax3, kmin3, kmax3
    !write(*,*) x_min, x_max, y_min, y_max,z_min, z_max 
!    stop 
    !call setgeometry_pmc(XEF, YNF, ZTF)  
!    do k=kmin2,kmax2
!     do j=jmin2,jmax2 
!       do i=imin2, imax2 
!          write(*,*)funijk(i,j,k),XEC(i), YNC(j), ZTC(k),dimension_3
!       end do 
!     end do 
!    end do 
!!    stop 
  
!    !Calculate Cell Center Coordinates------------------------------------------
!    do k = 1,kmax
!        do j = 1, jmax
!            do i = 1,imax
!                cellID = i+(j-1)*imax+(k-1)*imax*jmax
!                coords(cellID,1) = i*dx1-x_max/2d0-(dx1/2d0)
!                coords(cellID,2) = j*dy1-y_max/2d0-(dy1/2d0)
!                coords(cellID,3) = k*dz1-(dz1/2d0)
!            end do
!        end do
!    end do
    pres = 1.0d0 
    xH2O = 0.2d0 
    xCo2 = 0.1d0 
    volfrac = 1.0d0 
    d_part = 1.0d-3
    npart_cell = 100
    rho   = 1.0d0
!    T_g = 300.0d0  

!   Read total number of photons to consider 
    open(21,file='input.dat') 
    read(21,*) nrays_tot
    close(21) 
!    nrays_tot = 10000000 
   !Location of cell centers 	
   do k=1,kmax 
     do i=1,imax 
      do j=1,jmax 
        ijk = i+(j-1)*imax + (k-1)*imax*jmax 
        XEF(i) = x_min + (i-1)*dx(i)  + dx(i)/2.0
        YNF(j) = y_min + (j-1)*dy(j)  + dy(j)/2.0
        ZTF(k) = z_min + (k-1)*dz(k)  + dz(k)/2.0
        !write(*,*) ijk, i, XEF(i)
      end do 
     end do 
    end do 

    !Calculate Temperatures via cell centers 
    do k=kmin2, kmax
     do i=imin2, imax
       do j=jmin2, jmax
        ijk = i+(j-1)*imax + (k-1)*imax*jmax !funijk(i,j,k) 
        x_c = XEF(i) !+dx(i)/2.0d0
        y_c = YNF(j) !+dy(j)/2.0d0
        z_c = ZTF(k) !+dz(k)/2.0d0 
        rad = sqrt((x_c-1.0)**2.0 + (y_c-1.0d0)**2.0)
        if (rad < 1d0)then
            if (z_c <= 0.375d0)then
                Tc = 400d0+1400d0*(z_c/0.375d0)
            else
                Tc = 1800d0-1000d0*(z_c-0.375d0)/3.625d0
            end if
            T_g(ijk) = (Tc-800d0)*(1d0-3d0*rad**2+2d0*rad**3)+800d0
        else
            T_g(ijk) = 800d0
        end if
        if(k.eq.kmin2) T_g(ijk) = 400.0d0 
        if(k.ge.kmax)  T_g(ijk) = 800.0d0 
        if(j.eq.jmin2) T_g(ijk) = 300.0d0 
        if(j.ge.jmax)  T_g(ijk) = 300.0d0 
        if(i.eq.imin2) T_g(ijk) = 300.0d0 
        if(i.ge.imax)  T_g(ijk) = 300.0d0 
!        write(*,*) i, j, k, ijk, T_g(ijk), x_c, xec(i)
       end do
     end do
    end do 

    emise = 0.0d0 
    E_out_sum = 0.0d0
!    !Calculate Emissive Power at cell centers 
!    !Calculate total Emissive Power
    do k = kmin2, kmax
      do i = imin2, imax 
          do j=jmin2, jmax  !1,N_cells
          ijk =i+(j-1)*imax + (k-1)*imax*jmax ! funijk(i,j,k) 
!!        eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!        abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)  
          emise(ijk) = 4.0d0*abs_coef*stef_boltz*(T_g(ijk))**4.0
     !     write(*,*) emise(ijk),abs_coef, stef_boltz, T_g(ijk)
          ! dv = dx(i)*dy(j)*dz(k) 
           E_out_sum = E_out_sum + emise(ijk)*dx(i)*dy(j)*dz(k)
          end do
      end do 
    end do 

!    stop
     PPR = E_out_sum/float(nrays_tot)
!    !MONTE CARLO SIMULATION-----------------------------------------------------
!    !Loop Over All Cells
    ncountr = 0 
     do k=kmin2, kmax
       do i=imin2, imax 
          do j=jmin2, jmax
          ijk = i+(j-1)*imax + (k-1)*imax*jmax !funijk(i,j,k) 
!!        !Get Center Coordinates of Current Cell
          !x_c = XEF(i) 
          !y_c = YNF(j) 
          !z_c = ZTF(k) 
!!        !Determine Number of Rays to emit in the current cell 
          dv = dx(i)*dy(j)*dz(k) 
          nrays = nint(emise(ijk)*dv/PPR) 
          ncountr = ncountr + nrays

!        Absorb Energy due to rounding off the number of rays  at current cell
         abse(ijk) = abse(ijk) + emise(ijk)*dv-nrays*PPR  
         ip  = i
         jp  = j 
         kp  = k 
        !loop over all Rays emit from current cell
         if(nrays.lt.1) cycle 
         do l = 1, nrays
!            !Display Current Progress
            write(*,10,advance='no')(bsp,n=1,70),ijk,imax*jmax*kmax, l, nrays
            10 format(70A1,'Emitting Rays from Cell: ',I5,'/',I5, &
                '     Tracing Ray: ',I7,"/",I7)
!
!            !Pick Uniform Point of Emission
            call random_number(rand)
            x_ray = XEF(i)+rand*dx(ip) - dx(ip)/2.0  !dx(ip)
            call random_number(rand)
            y_ray = YNF(j)+rand*dy(jp) - dy(jp)/2.0 !dy(jp)
            call random_number(rand)
            z_ray = ZTF(k)+rand*dz(kp) -dz(kp)/2.0 !dz(kp)

!            !Pick Diffuse Direction of Emission
            call random_number(rand)
            theta = ACOS(1d0-2d0*rand)
            call random_number(rand)
            phi = 2d0*PI*rand  
            !write(*,*) i,j,k, x_ray,y_ray,z_ray
!test1
!!	    find new cell 
!            !Find Distance to Cell Edge in Current Direction
          !  call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
            call find_length(x_ray,y_ray,z_ray, theta, phi,pl) 
            !write(*,*) "pl", x_ray, y_ray, z_ray, theta, phi, pl
!
!            !Ray Tracing Until Ray Energy is Below Threshold Energy
             ray_energy  = PPR

            do while(ray_energy > PPR*PPR_crit_perc/100d0)
!!                !Move to Cell Edge
           !     call find_length(ip,jp,kp, x_ray,y_ray,z_ray, theta, phi,pl) 
                x_ray = x_ray + pl*sin(theta)*cos(phi)
                y_ray = y_ray + pl*sin(theta)*sin(phi)
                z_ray = z_ray + pl*cos(theta)

            !    call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                !call find_length(ip,jp,kp, x_ray,y_ray,z_ray, theta, phi,pl) 

!!                !Absorb Fraction of Ray Energy to Cell
!!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
!!!                constant absorption 
!!!                call graygasconst(T_g(cellID),kg) 
!!!                abs_coef = kg(1)
                  
      !          call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
                ray_energy = ray_energy - E_abs
               ! write(*,*) "cellid",  ijk
                abse(ijk) =abse(ijk) + E_abs
!!
!!                !Check if Ray Hit a Surface
                if (x_ray.le.x_min+eps .OR. x_ray.ge.x_max-eps  .OR. &
                    y_ray.le.y_min+eps .OR. y_ray.ge.y_max-eps  .OR. & 
                    z_ray.le.z_min+eps .OR. z_ray.ge.z_max-eps) then
!!                   
!!                   !Absorb Fraction of Ray Energy to Surf
                    call getwallbcnode(x_ray,y_ray,z_ray,wlbcid)
                    !call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                    E_abs       = em_surf*ray_energy 
                    abse(wlbcid)= abse(wlbcid) + E_abs
                    ray_energy  = ray_energy   - E_abs
!!
!!                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf .eq. 1.0d0)then
                        EXIT
                    end if
            !        exit
                    call random_number(rand)
                    call random_number(rand2)
                    call reflect(x_ray,y_ray,z_ray,rand,rand2,theta,phi)
                end if

             !   E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
             !   ray_energy = ray_energy - E_abs
            !    write(*,*) "before abse",  ip, jp, kp, ijk
                !abse(ijk) =abse(ijk) + E_abs

!!              Find new cell 
        !        call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk) 
!!              !Find New Distance to Next Cell Edge in Current Direction
                call find_length(x_ray,y_ray,z_ray, theta, phi,pl) 
!!                
!!                !Find New Cell ID
                x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
             !   x_mid = x_ray + pl*sin(theta)*cos(phi)
                y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
             !   y_mid = y_ray + pl*sin(theta)*sin(phi)
                z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
             !   z_mid = z_ray + pl*cos(theta)
                 call find_celln(x_mid, y_mid, z_mid, ip, jp, kp,ijk) 
              !  call getcellid(x_mid, y_mid, z_mid, dx1, dy1, dz1, ip, jp, kp, ijk) 
!!                cellID =  GetcellID(x_mid,y_mid,z_mid)
            end do

            !Absorb Final Small Fraction of Ray Energy if Below Threshold
          !  ijk=ip+(jp-1)*imax + (kp-1)*imax*jmax  !funijk(ip,jp,kp)
           ! write(*,*) ijk, ip, jp, kp, ray_energy
           ! abse(ijk) = abse(ijk) + ray_energy
           ! ray_energy = 0d0
         end do 
        end do
    end do
  end do 
!

! stop 
  do m = 1, DIMENSION_BC
    if (BC_DEFINED(m)) then
! The range of boundary cells
            I1 = BC_I_W(m)
            I2 = BC_I_E(m)
            J1 = BC_J_S(m)
            J2 = BC_J_N(m)
            K1 = BC_K_B(m)
            K2 = BC_K_T(m)
!         write(*,*) L, I1, I2, J1, J2, K1, K2
      do k=k1,k2
        do i=i1,i2
          do j=j1,j2
            ijk = i+(j-1)*imax + (k-1)*imax*jmax !funijk(i,j,k)

            if(i1.eq.I2) dv = dy(j)*dz(k)
            if(j1.eq.j2) dv = dx(i)*dz(k)
            if(k1.eq.k2) dv = dx(i)*dy(j)
        
!!          !Get Center Coordinates of Current Cell
            !x_c = XEF(i)
            !y_c = YNF(j)
            !z_c = ZTF(k)
!!          !Determine Number of Rays to emit in the current cell 
            !nrays = nint(emise(ijk)*vol(ijk)/PPR) 
            nrays = nint(emise(ijk)*dv/PPR)
            ncountr = ncountr + nrays

!           Absorb Energy left due to round off at current cell
            abse(ijk) = abse(ijk) + emise(ijk)*dv-nrays*PPR
            ip  = i
            jp  = j
            kp  = k
            if(nrays.lt.1) cycle
        !loop over all Rays emit from current cell
     !    if(nrays.lt.1) cycle
         do l = 1, nrays
!            !Display Current Progress
            write(*,10,advance='no')(bsp,n=1,70),ijk, imax*jmax, l, nrays
           ! 10 format(65A1,'Emitting Rays from Cell: ',I5,'/',I5, &
           !     '     Tracing Ray: ',I5,"/",I5)
!
!            !Pick Uniform Point of Emission

            call random_number(rand)
            x_ray = XEF(i)+rand*dx(ip) - dx(ip)/2.0  !dx(ip)
            call random_number(rand)
            y_ray = YNF(j)+rand*dy(jp) - dy(jp)/2.0 !dy(jp)
            call random_number(rand)
            z_ray = ZTF(k)+rand*dz(kp) -dz(kp)/2.0 !dz(kp)

!            !Pick Diffuse Direction of Emission
            call random_number(rand)
            theta = ACOS(1d0-2d0*rand)
            call random_number(rand)
            phi = 2d0*PI*rand
            !write(*,*) i,j,k, x_ray,y_ray,z_ray
!test1
!!          find new cell 
          !  call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
            call find_length(x_ray,y_ray,z_ray, theta, phi,pl)
            !write(*,*) "pl", x_ray, y_ray, z_ray, theta, phi, pl
!            !Find Distance to Cell Edge in Current Direction
!           D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
!
!            !Ray Tracing Until Ray Energy is Below Threshold Energy
             ray_energy  = PPR
             !PPR = 0.1d0
             !ray_energy= 1.0d0

            do while(ray_energy > PPR*PPR_crit_perc/100d0)
!!                !Move to Cell Edge

           !     call find_length(ip,jp,kp, x_ray,y_ray,z_ray, theta, phi,pl) 
           !     call find_length(ip,jp,kp, x_ray,y_ray,z_ray, theta, phi,pl) 
                x_ray = x_ray + pl*sin(theta)*cos(phi)
                y_ray = y_ray + pl*sin(theta)*sin(phi)
                z_ray = z_ray + pl*cos(theta)

            !    call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                !call find_length(ip,jp,kp, x_ray,y_ray,z_ray, theta, phi,pl) 

!!                !Absorb Fraction of Ray Energy to Cell
!!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
!!!                constant absorption 
!!!                call graygasconst(T_g(cellID),kg) 
!!!                abs_coef = kg(1)

      !          call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
                ray_energy = ray_energy - E_abs
               ! write(*,*) "cellid",  ijk
                abse(ijk) =abse(ijk) + E_abs
!!
!!                !Check if Ray Hit a Surface
                if (x_ray.le.x_min+eps .OR. x_ray.ge.x_max-eps  .OR. &
                    y_ray.le.y_min+eps .OR. y_ray.ge.y_max-eps  .OR. &
                    z_ray.le.z_min+eps .OR. z_ray.ge.z_max-eps) then
!!                   
!!                   !Absorb Fraction of Ray Energy to Surf
                    call getwallbcnode(x_ray,y_ray,z_ray,wlbcid)
                    !call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk)
                    E_abs      = em_surf*ray_energy
                    abse(wlbcid)  = abse(wlbcid) +  E_abs
                    ray_energy = ray_energy - E_abs
!!
!!                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf .eq. 1.0d0)then
                        EXIT
                    end if
            !        exit
                    call random_number(rand)
                    call random_number(rand2)
                    call reflect(x_ray,y_ray,z_ray,rand,rand2,theta,phi)
                end if

             !   E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
             !   ray_energy = ray_energy - E_abs
            !    write(*,*) "before abse",  ip, jp, kp, ijk
                !abse(ijk) =abse(ijk) + E_abs

!!              Find new cell 
        !        call find_celln(x_ray, y_ray, z_ray, ip, jp, kp,ijk) 
!!              !Find New Distance to Next Cell Edge in Current Direction
                call find_length(x_ray,y_ray,z_ray, theta, phi,pl)
!!                
!!                !Find New Cell ID
                x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
             !   x_mid = x_ray + pl*sin(theta)*cos(phi)
                y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
             !   y_mid = y_ray + pl*sin(theta)*sin(phi)
                z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
             !   z_mid = z_ray + pl*cos(theta)
                 call find_celln(x_mid, y_mid, z_mid, ip, jp, kp,ijk)
              !  call getcellid(x_mid, y_mid, z_mid, dx1, dy1, dz1, ip, jp, kp, ijk) 
!!                cellID =  GetcellID(x_mid,y_mid,z_mid)
            end do

            !Absorb Final Small Fraction of Ray Energy if Below Threshold
          !  ijk=ip+(jp-1)*imax + (kp-1)*imax*jmax  !funijk(ip,jp,kp)
           ! write(*,*) ijk, ip, jp, kp, ray_energy
           ! abse(ijk) = abse(ijk) + ray_energy
           ! ray_energy = 0d0

            end do
          end do
        end do
      end do

    end if
  end do


    srad = 0.0d0
!    !Calculate Source Terms in Cells--------------------------------------------
    do k=kmin2, kmax
      do j=jmin2, jmax
        do i=imin2, imax
          ijk = i+(j-1)*imax + (k-1)*imax*jmax !funijk(i,j,k)
          dv =  dx(i)*dy(j)*dz(k)
          srad(ijk) =  abse(ijk)/dv - emise(ijk)
          if(i.eq.imax/2.and.j.eq.jmax/2) write(11,*) XEF(i),YNF(j),ZTF(k),T_g(ijk),srad(ijk)
        end do 
      end do 
    end do  
!    end do
!
!    !Calculate Flux Terms in Surfs----------------------------------------------
!    do i = 1,N_surfs
!       faceID = bc%face(i,1)
!       bc%gy(i,3) = bc%energy(i,2)/da(faceID) - &
!                               bc%gy(i,1)
!    end do
!
!    !Save Centerline Data for Plots---------------------------------------------
!    write(grid_t,'(i2,a1,i2,a1,i2)') imax,'x',jmax,'x',kmax
!    write(abs_coef_t,'(1f5.3)') abs_coef
!    write(PPR_t,'(1f5.3)') PPR
!    filename = 'Porter'//'Grid_'//trim(grid_t)//'a_' &
!                //trim(abs_coef_t)//'PPR_'//trim(PPR_t)//'.dat'
!    open (unit = 1, file = filename)
!    do i = int(imax*jmax/2d0+0.5d0),int(imax*jmax*kmax-imax*jmax/2d0+0.5d0),imax*jmax
!       write (1,"(3f9.3,1f12.3,1f12.3)")  coords(i,1:3), &
!                                          T_g(i), &
!                                          cell_data%gy(i,3)
!    end do
!
!    !Display Final Results------------------------------------------------------
!
!    !Energy Balance Check
!    do i = 1,N_cells
!        E_out_sum = E_out_sum + cell_data%gy(i,1)*dv 
!        E_in_sum = E_in_sum + cell_data%gy(i,2)
!    end do
!    do i = 1,N_surfs
!        faceID = bc%face(i,1)
!        E_out_sum = E_out_sum + bc%gy(i,1)*da(faceID) 
!        E_in_sum = E_in_sum + bc%gy(i,2)
!    end do
!    write(*,*) New_Line(c)
!    write(*,"(a22,1f10.2,a4)") "Total Energy Emitted :", E_out_sum, " (W)"
!    write(*,"(a22,1f10.2,a4)") "Total Energy Absorbed:", E_in_sum, " (W)"
!
!    !Total Number of Rays Traced
!    write(*,*)
!    write(*,"(a28,1i9)") "Total Number of Rays Traced:", N_rays_tot
!
!    !CPU Time
!    call CPU_TIME(time2)
!    write(*,*)
!    write(*,"(a15,1f7.2,a4)") "Total CPU Time:", time2-time1, " (s)"
!
!    !Saved File
!    write(*,*)
!    write(*,"(a30,a46)") "Saving Centerline Results to: ", trim(filename)

!write(*,*)
!contains
end subroutine rad_pmc_calc

subroutine find_celln(ray_x,ray_y,ray_z, ix, jy, kz,ijknew)
          
          integer :: ix, jy , kz
          integer, intent(out) :: ijknew
          integer :: inew, jnew, knew,i,j,k
          double precision :: ray_x, ray_y, ray_z 
          double precision :: dx, dy, dz 

!        !Check For Valid Coordinates
!        if (ray_x > x_max .or. ray_x < x_min)then
!            print *, "Invalid x Coordinate:", ray_x
!            ijknew= 0
!            return
!        elseif (ray_y > y_max.or. ray_y < y_min)then
!            print *, "Invalid y Coordinate:", ray_y
!            ijknew = 0
!            return
!        elseif (ray_z > z_max.OR. ray_z < z_min )then
!            print *, "Invalid z Coordinate:", ray_z
!            ijknew  = 0
!            return
!        end if 
                
          dx = x_max/float(imax) 
          dy = y_max/float(jmax)
          dz = z_max/float(kmax) 

        !Adjust Coordinates if Point is on Positive Surface
        if (abs(ray_x - x_max) < eps)then
            ray_x = ray_x - dx/2d0
        end if
        if (abs(ray_y - y_max) < eps)then
            ray_y = ray_y - dy/2d0
        end if
        if (abs(ray_z - z_max) < eps)then
            ray_z = ray_z - dz/2d0
        end if


          inew = ceiling(ray_x/dx)
          jnew = ceiling(ray_y/dy)
          knew = ceiling(ray_z/dz)

!        !Check For Valid Coordinates
          if(ray_x .gt. x_max) inew = imax 
          if(ray_x .lt. x_min) inew = 1 
          if(ray_y .gt. y_max) jnew = jmax 
          if(ray_y .lt. y_min) jnew = 1 
          if(ray_z .gt. z_max) knew = kmax 
          if(ray_z .lt. z_min) knew = 1 

!          if(ray_x.lt.XEF(ix)) then
!               inew =ix-1
!               if(inew.le.1) inew = 1 
!          else if(ray_x.gt.XEF(ix)) then 
!               inew = ix+1
!               if(inew.ge.imax) inew = imax
!          else 
!               inew = ix 
!          end if
!
!          if(ray_y.lt.YNF(jy)) then
!               jnew = jy-1
!               if(jnew.le.1) jnew = 1
!          else if(ray_y.gt.YNF(jy)) then 
!               jnew = jy+1
!               if(jnew.ge.jmax) jnew = jmax
!          else 
!               jnew = jy 
!          end if
!
!            if(ray_z.lt.ZTF(kz)) then
!               knew = kz-1
!               if(knew.le.1) knew = 1
!             else if(ray_z.gt.ZTF(kz)) then 
!               knew = kz+1 
!               if(knew.ge.kmax) knew = kmax
!             else 
!               knew = kz 
!            end if
          
            ijknew = inew+(jnew-1)*imax + (knew-1)*imax*jmax !funijk(inew,jnew,knew) 
            ix = inew 
            jy = jnew 
            kz = knew 
!            write(*,*) "cell number", ix,jy,kz,XEF(ix), YNF(jy),ZTF(kz)
end subroutine find_celln 

!    !Current Cell ID From Coordinates Function----------------------------------
    subroutine getcellid(x,y,z,ip, jp, kp, cellID)
!        !This function returns the cell ID number that a point is in, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. If point is on cell face/edge/corner, cellID of the positive-side cell is returned (except for when on a positive surface). Cell ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z.
        implicit none
!        !Initialize Variables
        double precision :: x, y, z
        double precision :: dx1, dy1, dz1
        integer, intent(out) :: cellID
        integer :: ip, jp, kp 

       dx1 = x_max/float(imax) 
       dy1 = y_max/float(jmax) 
       dz1 = z_max/float(kmax) 

!        !Check For Valid Coordinates
        if (x > x_max .or. x < x_min)then
            print *, "Invalid x Coordinate:", x
            cellID = 0
            return
        elseif (y > y_max.or. y < y_min)then
            print *, "Invalid y Coordinate:", y
            cellID = 0
            return
        elseif (z > z_max.OR. z < z_min )then
            print *, "Invalid z Coordinate:", z
            cellID = 0
            return           
        end if
!
!        !Adjust Coordinates if Point is on Positive Surface
        if (abs(x - x_max) < eps)then
            x = x - dx1/2d0
        end if
        if (abs(y - y_max) < eps)then
            y = y - dy1/2d0
        end if    
        if (abs(z - z_max) < eps)then
            z = z - dz1/2d0
        end if
!
!        !Calculate Cell ID
        cellID = 1+floor(x/dx1)+floor(y/dy1)*imax+ &
                 floor(z/dz1)*imax*jmax
        ip = floor(x/dx1) 
        jp = floor(y/dy1)
        kp = floor(z/dz1) 
    end subroutine getcellid
!
!    !Current Surface ID From Coordinates Function-------------------------------
    subroutine getwallbcnode(x,y,z,noden)
!        !This function returns the surface ID number that a point is on, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Note that if point is not exactly on a surface, it is snapped to the nearest surface. Surface ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z (on the current face), with faces ordered as -x,x,-y,y,-z,z.
        implicit none
!        !Initialize Variables
        double precision :: x, y, z
        double precision :: dx, dy, dz
    !    double precision, dimension(6) :: err
        integer :: noden 
        integer :: inew,jnew,knew 

          dx = x_max/float(imax)
          dy = y_max/float(jmax)
          dz = z_max/float(kmax)

          inew = ceiling(x/dx)
          jnew = ceiling(y/dy)
          knew = ceiling(z/dz)

!        !Check For Valid Coordinates
          if(x-eps .ge. x_max) inew = imax
          if(x .le. x_min+eps) inew = 1
          if(y-eps .ge. y_max) jnew = jmax
          if(y .le. y_min+eps) jnew = 1
          if(z-eps .ge. z_max) knew = kmax
          if(z .le. z_min+eps) knew = 1

          noden = inew + (jnew-1)*inew + (knew-1)*imax*jmax
  !      dx1 = x_max/float(imax) 
  !      dy1 = y_max/float(jmax) 
  !      dz1 = z_max/float(kmax) 
  !       
  !      !Snap Coordinate to Nearest Surface
  !      if (x /= x_max.AND. x /=0d0 .and. y /= y_max .AND. y /=0d0 .and. &
  !          z /= z_max.AND. z /=0d0)then
  !          err(1) = abs(x)
  !          err(2) = abs(x-x_max)
  !          err(3) = abs(y)
  !          err(4) = abs(y-y_max)
  !          err(5) = abs(z)
  !          err(6) = abs(z-z_max)
  !          if (err(1) == minval(err))then
  !              x = 0d0
  !          elseif (err(2) == minval(err))then
  !              x = x_max
  !          elseif (err(3) == minval(err))then
  !              y = 0d0
  !          elseif (err(4) == minval(err))then
  !              y = y_max
  !          elseif (err(5) == minval(err))then
  !              z = 0d0    
  !          elseif (err(6) == minval(err))then
  !              z = z_max
  !          else
  !              print *, "Surface Snap Error"
  !          end if
  !      end if

  !      !Display Warning if Snap Distance was not Small
  !      if (minval(err) >= eps)then
  !          print *, "Warning: Surface Snap Distance >= eps"
  !      end if

  !      !Calculate Surface ID
  !      if (x == 0.0d0) then
  !          surfID = 1+floor(y/dy1)+floor(z/dz1)*jmax
  !      else if (x == x_max) then
  !          surfID = 1+floor(y/dy1)+floor(z/dz1)*jmax + jmax*kmax
  !      else if (y == 0.0d0) then
  !          surfID = 1+floor(x/dx1)+floor(z/dz1)*imax + 2*jmax*kmax
  !      else if (y == y_max) then
  !          surfID = 1+floor(x/dx1)+floor(z/dz1)*imax + 2*jmax*kmax &
  !                  + imax*kmax
  !      else if (z == 0d0) then
  !          surfID = 1+floor(x/dx1)+floor(y/dy1)*imax &
  !                  + 2*(jmax*kmax+imax*kmax)
  !      else if (z == z_max) then
  !          surfID = 1+floor(x/dx1)+floor(y/dy1)*imax &
  !                  + 2*(jmax*kmax+imax*kmax) + imax*jmax
  !      else 
  !          print *, "Surface ID Error"
  !          surfID = 0
  !      end if
    end subroutine getwallbcnode 

    !Distance to Cell Edge in Current Direction Function------------------------
     subroutine find_length(x,y,z,theta, phi, pl) 
!    function Get_D_edge(x,y,z,theta,phi) result(D_edge)
        !This function computes the distance from the current point (defined by x,y,z) to the bounding cell surface in the current direction of emission (defined by input cone angle theta [0,PI] and azimuthal angle phi [0,2*PI]). Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. Point may be on a cell face/edge/corner. More information on page 158 of Mahan's book.
        implicit none
        !Initialize Variables
        double precision, intent(in)  :: x,y,z,theta,phi
        double precision, intent(out) :: pl
        integer :: i,j,k
        integer :: ix,jy,kz,ijknew
        double precision :: x0,y0,z0,dl,dm,dn,x_err,y_err,z_err
        double precision :: dx,dy,dz
        double precision, dimension(6) :: t = -1 
        integer :: l, cellID,ijk
        !Find Direction Cosines
        dx = x_max/float(imax) 
        dy = y_max/float(jmax) 
        dz = z_max/float(kmax) 

        dl = sin(theta)*cos(phi)
        dm = sin(theta)*sin(phi)
        dn = cos(theta)
        !Find Cell Corner
       call find_celln(x,y,z, ix, jy, kz,ijknew) 
!        cellID = ijk  !GetcellID(x,y,z)
        x0 = XEF(ix) -dx/2.0d0 !i_of(ijk)) 
        y0 = YNF(jy) -dy/2.0d0 !_of(ijk)) 
        z0 = ZTF(kz) -dz/2.0d0 !_of(ijk))
        !write(*,*) "d_edge", x, y, z, ijknew, x0, y0, z0
       ! stop 
        !write(*,*) "find_legth", ijk, i_of(ijk), j_of(ijk), k_of(ijk), x0, y0, z0
        !write(*,*) "find_legth", i, j, k, x0, y0, z0
!        !Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
        x_err = abs(mod(x+mod(imax,2)*dx/2d0,dx))
        y_err = abs(mod(y+mod(jmax,2)*dy/2d0,dy))
        z_err = abs(mod(z,dz)) !(kz))) !_of(ijk))))
        !if (x_err <= eps .OR.dx(i_of(ijk))-x_err <= eps)then
        if (x_err <= eps .OR.dx-x_err <= eps)then
            if (dl > 0)then
                x0 = x
            else
                x0 = x - dx !dx(ix)  !dx(i_of(ijk))
            end if
        end if
!        if (y_err <= eps .OR.dy(j_of(ijk))-y_err <= eps)then
        if (y_err <= eps .OR.dy-y_err <= eps) then
            if (dm > 0)then
                y0 = y
            else
                y0 = y - dy !dy(jy) !dy(j_of(ijk))
            endif
        end if
!        if (z_err <= eps .OR.dz(k_of(ijk))-z_err <= eps)then
        if (z_err <= eps .OR.dz-z_err <= eps)then
            if (dn > 0)then
                z0 = z
            else
                z0 = z - dz !dz(kz) !dz(k_of(ijk))
            end if
        end if
        !Find All t's
        if (abs(dl) >= eps)then 
            t(1) = (x0-x)/dl
            !t(2) = (x0+dx(i_of(ijk))-x)/dl
            !t(2) = (x0+dx(ix)-x)/dl
            t(2) = (x0+dx-x)/dl
        end if
        if (abs(dl) >= eps)then
            t(3) = (y0-y)/dm
            !t(4) = (y0+dy(j_of(ijk)) -y)/dm
            t(4) = (y0+dy -y)/dm
        end if
        if (abs(dn) >= eps)then          
            t(5) = (z0-z)/dn
            !t(6) = (z0+dz(k_of(ijk))-z)/dn
            t(6) = (z0+dz-z)/dn
        end if

        !Find Distance to Edge (Smallest, Positive t)
        pl = 1d15
        do l = 1,6
            if (t(l) > 0d0 .AND. t(l) < pl)then
                pl = t(l)
            end if
        end do
     return 
     end subroutine find_length 
!    end function Get_D_edge
!
!    !Reflecting Direction Function----------------------------------------------
    subroutine  reflect(x_r,y_r,z_r,rand,rand2, theta, phi)  
!        !This function provides the new ray direction after a reflecting event occurs, relative to the global cartesian coordinates, given the faceID at which the reflection occurs(-x=1, +x=2, -y=3, +y=4, -z=5, +z=6) and two random numbers seeded between 0 and 1 (rand and rand2). The reflecting surface is assumed to be diffuse. (This function can also be utilized when emmitting rays from surfaces not normal to the global coordinate system.)
        implicit none
!        !Initialize Variables
        integer :: IDface
        integer :: i, j, k 
        double precision, intent(in) :: rand, rand2
        double precision :: theta, phi, theta_p, phi_p
        double precision ::  x_p, y_p, z_p, x_p_t, y_p_t, z_p_t
        double precision ::  x_r, y_r, z_r
!        double precision, dimension(2) :: reflected_direction
        double precision, dimension(6,2) :: face_angles
        double precision, dimension(3,3) :: Q
        double precision :: dx, dy, dz 
        
       
!        !Define Angles (Theta, Phi) of Face Normal Vectors
        face_angles(1,1) = PI/2D0
        face_angles(1,2) = 0D0
        face_angles(2,1) = PI/2D0
        face_angles(2,2) = PI
        face_angles(3,1) = PI/2D0
        face_angles(3,2) = PI/2D0
        face_angles(4,1) = PI/2D0
        face_angles(4,2) = 3D0*PI/2D0
        face_angles(5,1) = 0D0
        face_angles(5,2) = 0D0
        face_angles(6,1) = PI
        face_angles(6,2) = 0D0
!
        if(x_r .le.x_min+eps ) idface = 1 
        if(x_r .ge.x_max-eps ) idface = 2  
        if(y_r .le.y_min+eps ) idface = 3 
        if(y_r .ge.y_max-eps ) idface = 4 
        if(z_r .le.z_min+eps ) idface = 5 
        if(z_r .ge.z_max-eps ) idface = 6 
      !  if(i.eq.imin1)   idface = 1 
      !  if(i.ge.imax2-1) idface = 2
      !  if(j.eq.jmin1)   idface = 3
      !  if(j.ge.jmax2-1) idface = 4
      !  if(k.eq.kmin1)   idface = 5
      !  if(k.ge.kmax2-1) idface = 6
 

!        !Get Angles of Current Face Normal Vector
        theta = face_angles(idface,1)
        phi = face_angles(idface,2)
!
!        !Calculate Reflected Direction (Angles) Relative to Face Normal
        theta_p = asin(sqrt(rand))
        phi_p = 2D0*PI*rand2
!
!        !Calculate Reflected Direction (Vector) Relative to Face Normal
        x_p = sin(theta_p)*cos(phi_p)
        y_p = sin(theta_p)*sin(phi_p)
        z_p = cos(theta_p)
!
!        !Define the Transformation Matrix
        Q(1,1) = cos(phi)*cos(theta)
        Q(1,2) = -sin(phi)
        Q(1,3) = sin(theta)*cos(phi)
        Q(2,1) = sin(phi)*cos(theta)
        Q(2,2) = cos(phi)
        Q(2,3) = sin(theta)*sin(phi)
        Q(3,1) = -sin(theta)
        Q(3,2) = 0
        Q(3,3) = cos(theta)
!
!        !Transform Reflected Direction (Vector) to Global Coordinate System
        x_p_t = x_p*Q(1,1)+y_p*Q(1,2)+z_p*Q(1,3)
        y_p_t = x_p*Q(2,1)+y_p*Q(2,2)+z_p*Q(2,3)
        z_p_t = x_p*Q(3,1)+y_p*Q(3,2)+z_p*Q(3,3)

!        !Back Calculate Reflected Direction (Angles) Relative to GCS
        theta  = atan2(sqrt(x_p_t**2+y_p_t**2),z_p_t) 
        phi    = atan2(y_p_t,x_p_t)
    end subroutine reflect
!!Combined Absorption Coefficient Function---------------------------------                                                                                                                                   
!    function Absorp_Coeff(eta,rho,T,P_abs,N_part,d_part) result(abs_coef)
!!This function returns the combined absorption coefficient of a CO2 gas and carbon particle mixture (m^-1), using the Elsasser narrow-band model (Edwards 1967), given the current wave number (cm^-1), density (gm/m^3), temperature (K), partial pressure of CO2 (atm), number of carbon particles (m^-3), and diameter of particles (m). 
!        implicit none 
!        !Initialize Variables
!        integer, intent(in) :: N_part
!        integer :: i 
!        double precision, parameter :: pi = 22.0/7.0
!        double precision, intent (in) :: eta, rho, T, P_abs, d_part
!        double precision, parameter :: h = 6.6237D-27, k = 1.3802D-16, &
!            c = 2.9979D10
!        double precision :: abs_coef, nu1 = 1351D0, nu3 = 2396D0, &
!            phi1, phi2, phi3, T_ref = 100D0, P_ref = 1D0, Pe, delta, &
!            beta, Sc_d 
!        double precision, dimension(5) :: & 
!            eta_c = (/667D0,960D0,1060D0, 2410D0,3660D0/), & 
!            b = (/1.3D0,1.3D0,1.3D0,1.3D0,1.3D0/), &
!            n = (/0.7D0,0.8D0,0.8D0,0.8D0,0.65D0/), &
!            a = (/2D0,2D0,2D0,1D0,2D0/), C1, C2, C3, C3_ref
!
!        !Calculate Phi(T) Functions
!        phi1 = ((1D0-exp(-(h*c/k/T)*(nu3-nu1)))*(exp(-h*c*nu1/k/T)-0.5D0* &
!            exp(-2D0*h*c*nu1/k/T)))/((1D0-exp(-h*c*nu1/k/T))*&
!            (1D0-exp(-h*c*nu3/k/T)))
!        phi2 = ((1D0-exp(-h*c*(nu1+nu3)/k/T)))/((1D0-exp(-h*c*nu1/k/T))* &
!            (1D0-exp(-h*c*nu3/k/T)))
!        phi3 = (1D0+0.053D0*(T/T_ref)**(1.5D0))
!        !Calculate Constants
!        C1 = (/19D0,0.76D0*phi1,0.76D0*phi1,110D0,4D0*phi2/)
!        C2 = (/6.9D0 *(T/T_ref)**0.5D0,1.6D0*(T/T_ref)**0.5D0*C1(2)**0.5D0, & 
!            1.6D0*(T/T_ref)**0.5D0*C1(3)**0.5D0,31D0*(T/T_ref)**0.5D0, &
!            8.6D0*phi3/)
!        C3 = (/12.9D0*(T/T_ref)**0.5D0,12.4D0*(T/T_ref)**0.5D0, &
!            12.4D0*(T/T_ref)**0.5D0, 11.5D0*(T/T_ref)**0.5D0, &
!            24D0*(T/T_ref)**0.5D0/)
!        C3_ref = (/12.9D0*(T_ref/T_ref)**0.5D0,12.4D0*(T_ref/T_ref)**0.5D0, &
!            12.4D0*(T_ref/T_ref)**0.5D0, 11.5D0*(T_ref/T_ref)**0.5D0,24D0* &
!            (T_ref/T_ref)**0.5D0/)
!
!        !Loop Over All Major Spectral Absorption Bands
!        abs_coef = 0D0
!        do i = 1,size(eta_c)
!            Pe = (((1D0-P_abs)+b(i)*P_abs)/P_ref)**n(i)
!            delta = 30D0*C3_ref(i)
!            beta = C2(i)**2D0*Pe/4D0/C1(i)/C3(i)
!            Sc_d = C1(i)/C3(i)*exp(-(a(i)/C3(i))*abs(eta-eta_c(i)))
!            abs_coef = abs_coef+rho*Sc_d*sinh(PI*beta/2D0)/ &
!                (cosh(PI*beta/2D0)-cos(2*PI*(eta-eta_c(i))/delta))
!        end do
!
!        !Add Contribution from Carbon Particle Absorption (Howell Eq. 15.3)
!        !CURRENTLY ASSUMING ABSORPTION EFFICIENCY FACTOR = 1 FOR ALL ETA
!!        abs_coef = abs_coef + N_part*PI*(d_part**2D0)/4D0
!    end function Absorp_Coeff

      subroutine setgeometry_pmc(XEF, YNF, ZTF) 
      IMPLICIT NONE
! User specified dimension of the system (by default 2D, but if 3D system is
! desired then it must be explicitly specified)
      INTEGER, PARAMETER :: DIMN = 3

! Local Variables:
!---------------------------------------------------------------------//
! Gic loop indices
      INTEGER :: I, J, K
! Error Flag
      INTEGER :: IER
! Position of domain boundaries gally given as
!   (x_min, x_max, y_min, y_max, z_min, z_max)
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2
!......................................................................!
! Collect the error flags from all ranks. If all allocaitons were
! successful, do nothing. Otherwise, flag the error and abort.
    double precision :: XEF (0:DIMENSION_I)
    double precision :: YNF (0:DIMENSION_J)
    double precision :: ZTF (0:DIMENSION_K)



      CALL GLOBAL_ALL_SUM(IER)

! Set boundary edges.
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout
      EX2 = X_MAX;    WX1 = X_MIN ! East/West
      TY2 = Y_MAX;    BY1 = Y_MIN ! North/South
      NZ2 = Z_MAX;    SZ1 = Z_MIN ! Top/Bottom

! Initialize arrays.
      XEF(:) = ZERO
      YNF(:) = ZERO
      ZTF(:) = ZERO

! Set boundary edges.
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout
      EX2 = X_MAX;    WX1 = X_MIN ! East/West
      TY2 = Y_MAX;    BY1 = Y_MIN ! North/South
      NZ2 = Z_MAX;    SZ1 = Z_MIN ! Top/Bottom

! Initialize arrays.
      XEF(:) = ZERO
      YNF(:) = ZERO
      ZTF(:) = ZERO

! Each loop starts at 2 and goes to max+2 (i.e., imin1=2, imax2=imax+2)
! However, the indices range to include ghost cells (0-imax2) to avoid
! multiple if statements in particles_in_cell
      XEF(IMIN2-1) = X_MIN-DX(IMIN2)
      DO I = IMIN2, IMAX2
         XEF(I) = XEF(I-1) + DX(I) !+DX(I)/2.0d0 
        ! write(*,*) I, XEF(I)
      ENDDO

      YNF(JMIN2-1) = Y_MIN-DY(JMIN2)
      DO J  = JMIN2, JMAX2
         YNF(J) = YNF(J-1) + DY(J) !+DY(J)/2.0d0 
        ! write(*,*) j, YNF(j)
      ENDDO

      IF(DIMN.EQ.3) THEN
         ZTF(KMIN2-1) = Z_MIN-DZ(KMIN2)
         DO K = KMIN2, KMAX2
            ZTF(K) = ZTF(K-1) + DZ(K) !+ DZ(k)/2.0
       !     write(*,*) k, ZTF(k)
         ENDDO
      ENDIF

!      stop
      RETURN
   END SUBROUTINE setgeometry_pmc

end module rad_pmc_mod
