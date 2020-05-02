! photon monte carlo solver to find source terms due to thermal radiation 
! modified Apr 28, 2020 
module pmc_mod 
!use rad_gas_abs
!      use SET_GEOMETRY_PMC_mod, only :set_geometry_pmc
! Fluid grid cell dimensions and mesh size
USE geometry, only: DX, IMIN2, IMAX2
USE geometry, only: DY, JMIN2, JMAX2
USE geometry, only: DZ, KMIN2, KMAX2
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
!     implicit none
!     integer, parameter :: imax = 17, jmax = 17, kmax = 34 
!     double precision, parameter :: x_max= 2d0, y_max= 2d0, z_max= 4d0
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XEC  !(0:DIMENSION_I)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YNC  !(0:DIMENSION_J)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZTC  !(0:DIMENSION_K)

!      use set_geometry_pmc_mod, only:XEC,YNC,ZTC  
!      double precision, dimension (:,:), allocatable :: coords 
      double precision, dimension (:), allocatable :: T_g 
      double precision, dimension (:), allocatable :: engre,engra 
      double precision, parameter :: eps = 1d-12
      double precision, parameter :: PI = 4.d0*DATAN(1.d0)
      double precision, parameter :: stef_boltz = 5.670374d-8
!      Allocate( XEC (0:DIMENSION_I)) !, STAT=IER )
!      Allocate( YNC (0:DIMENSION_J)) !, STAT=IER )
!      Allocate( ZTC (0:DIMENSION_K)) !, STAT=IER )
contains 

    subroutine mcfurnace
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
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2
      USE functions, only: funijk 
      use geometry, only: VOL, AYZ, AXZ, AXY 
! Fixed array sizes in the I/J/K direction

     implicit none    
     double precision, parameter :: em_surf = 1d0
     double precision, parameter :: abs_coef = 0.1d0

    !Initialize Working Variables-----------------------------------------------
!    integer :: N_cells = imax*jmax*kmax
!    integer :: N_surfs = 2*(imax*jmax+imax*kmax+jmax*kmax)
    integer :: i,j,k,l, ijk, cellID, surfID, faceID
    integer :: nrays, nrays_tot,ncountr 
    double precision, parameter :: PPR_crit_perc = 0.1d0
    double precision :: PPR, pl 
    double precision :: dx1, dy1, dz1, dv, x_c, y_c, z_c, x_ray, y_ray, z_ray
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
    double precision, dimension (6) :: da
    double precision, dimension (2) :: reflected_direction
    character(len=1) :: c, bsp = char(8)
    character(len=50) :: filename, grid_t, abs_coef_t, PPR_t
    integer :: IER 
!    double precision, dimension (N_cells) :: T_g
!    double precision, dimension (N_cells,3) :: coords = 0
!    type cell_data_type
!        double precision, dimension (N_cells,3) :: gy = 0
!    end type
!    type(cell_data_type) :: cell_data
!    call CPU_TIME(time1)

    nrays_tot = 1e8 
!    allocate(coords(dimension_3,3))
    allocate(T_g(dimension_3))
    allocate(engre(dimension_3))
    allocate(engra(dimension_3))
    allocate( XEC (0:DIMENSION_I), STAT=IER )
    allocate( YNC (0:DIMENSION_J), STAT=IER )
    allocate( ZTC (0:DIMENSION_K), STAT=IER )

    !Calculate Size of Homogeneous Subregions---------------------------------- 
    dx1 = x_max/imax
    dy1 = y_max/jmax
    dz1 = z_max/kmax
!    dv = dx1*dy1*dz1
!    da = (/dy1*dz1,dy1*dz1,dx1*dz1,dx1*dz1,dx1*dy1,dx1*dy1/)
    write(*,*) x_max, y_max, z_max, imax, jmax, kmax 
    write(*,*) imin2, imax2,jmin2, jmax2, kmin2, kmax2
    write(*,*) imin3, imax3,jmin3, jmax3, kmin3, kmax3
   ! stop 
    call setgeometry_pmc(XEC, YNC, ZTC)  
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
    T_g = 300.0d0  
    !Calculate Cell Temperatures------------------------------------------------

    do k=kmin2, kmax2-1
     do i=imin2, imax2-1
       do j=jmin2, jmax2-1
        ijk = funijk(i,j,k) 
        x_c = XEC(i) !+dx(i)/2.0d0
        y_c = YNC(j) !+dy(j)/2.0d0
        z_c = ZTC(k) !+dz(k)/2.0d0 
        rad = sqrt(x_c**2 + y_c**2)
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
!        write(*,*) i, j, k, ijk, T_g(ijk), x_c, xec(i)
       end do
     end do
    end do 

!    stop 

!    
!    !Calculate Surf Center Coordinates------------------------------------------
!    !x-Faces
!    do k = 1,kmax
!        do j = 1,jmax
!            do i = 1,2
!                surfID = j+(k-1)*jmax+(i-1)*(jmax*kmax)
!                cellID = i+(j-1)*imax+(k-1)*imax*jmax
!                if(i.eq.2)cellID = imax+(j-1)*imax+(k-1)*imax*jmax
!                bc%node(surfID,1) = surfID
!                bc%face(surfID,1) = i
!                bc%glbnode(surfID,1) = cellID
!                T_g(cellID) = 300.0d0 
!            end do
!        end do
!    end do
!    !y-Faces
!    do k = 1,kmax
!        do i = 1,imax
!            do j = 1,2
!                surfID = i+(k-1)*imax+(imax*kmax)*(j-1)+2*(kmax*jmax)
!                cellID = i+(j-1)*imax+(k-1)*imax*jmax
!                if(j.eq.2)cellID = i+(jmax-1)*imax+(k-1)*imax*jmax
!                bc%node(surfID,1) = surfID
!                bc%face(surfID,1) = j+2
!                bc%glbnode(surfID,1) = cellID 
!                T_g(cellID) = 300.0d0 
!            end do
!        end do
!    end do
!    !z-Faces
!    do j = 1,jmax
!        do i = 1,imax
!            do k = 1,2
!                surfID = i+(j-1)*imax+(imax*jmax)*(k-1)+2*(kmax*jmax+imax*kmax)
!                cellID = i+(j-1)*imax+(k-1)*imax*jmax
!                if(k.eq.2)cellID = i+(j-1)*imax+(kmax-1)*imax*jmax
!                bc%node(surfID,1) = surfID
!                bc%face(surfID,1) = k+4
!                bc%glbnode(surfID,1) = cellID ! j*dx-x_enc/2d0-(dx/2d0)
!                T_g(cellID) = 300.0d0
!                if(k.eq.2) T_g(cellID) = 300.0d0
!            end do
!        end do
!    end do

   
     engre = 0.0d0 
!    !Calculate Cell Emissive Powers---------------------------------------------
    do k = kmin2, kmax2-1 
      do i = imin2, imax2-1 
          do j=1,jmin2, jmax2-1  !1,N_cells
          ijk = funijk(i,j,k) 
!!        eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!        abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)  
          engre(ijk) = 4d0*abs_coef*stef_boltz*(T_g(ijk))**4
          end do
      end do 
    end do 

!   Total gy 
    E_out_sum = 0.0d0
    do k=kmin2, kmax2-1 
      do i=imin2, imax2-1 
        do j=jmin2, jmax2-1    
           ijk = funijk(i,j,k)
           E_out_sum = E_out_sum + engre(ijk)*vol(ijk)
        end do
      end do 
    end do 

   PPR = E_out_sum/nrays_tot  
 !Calculate Surf Emissive Powers---------------------------------------------
!    do i = 1,N_surfs
!        bc%gy(i,1) = em_surf*stef_boltz*T_g(bc%glbnode(surfID,1))**4
!    end do
!
!    !MONTE CARLO SIMULATION-----------------------------------------------------
!    !Loop Over All Cells
    ncountr = 0 
     do k = kmin1, kmax-2 
       do i=imin1, imax2-1 
          do j=jmin1, jmax2-1
          ijk = funijk(i,j,k) 
!!        !Get Center Coordinates of Current Cell
          x_c = XEC(i) 
          y_c = YNC(j) 
          z_c = ZTC(k) 
!!        !Determine Number of Rays to emit in the current cell 
          nrays = nint(engre(ijk)*vol(ijk)/PPR) 
          ncountr = ncountr + nrays

!        Absorb rEnergy at current cell
         engra(ijk) = engra(ijk) + engre(ijk)*vol(ijk)-nrays*PPR  
!
        !loop over all Rays emit from current cell
         do l = 1, 10 !nrays
!            !Display Current Progress
!            write(*,10,advance='no')(bsp,k=1,65),i, N_cells, j, N_rays
!            10 format(65A1,'Emitting Rays from Cell: ',I5,'/',I5, &
!                '     Tracing Ray: ',I5,"/",I5)
!
!            !Pick Uniform Point of Emission
            call random_number(rand)
            x_ray = XEC(i)+rand*dx(i)
            call random_number(rand)
            y_ray = YNC(j)+rand*dy(j)
            call random_number(rand)
            z_ray = ZTC(k)+rand*dz(k)
!
!            !Reset cellID
            cellID = ijk
!
!            !Pick Diffuse Direction of Emission
            call random_number(rand) 
            theta = ACOS(1d0-2d0*rand)
            call random_number(rand)
            phi = 2d0*PI*rand  

!	    find new cell 
              if(x_ray.le.XEC(i)) then 
                 inew =i
               else 
                 inew = i+1 
                 if(inew.eq.imax2-1) inew = i
              end if 

              if(y_ray.le.YNC(j)) then
                 jnew =j
               else
                 jnew = j+1
                 if(jnew.eq.jmax2-1) jnew = j 
              end if
 
              if(z_ray.le.ZTC(k)) then
                 knew = k
               else
                 knew = k+1
                 if(knew.eq.kmax2-1) knew = k 
              end if 

             ijkn = funijk(inew,jnew,knew) 
!             call find_length(ijkn, x_ray,y_ray,z_ray, EC, YNC, ZTC, theta, phi,plen) 
             call find_length(ijkn, x_ray,y_ray,z_ray, theta, phi,pl) 

!            !Find Distance to Cell Edge in Current Direction
!            D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
!
!            !Ray Tracing Until Ray Energy is Below Threshold Energy
            ray_energy  = PPR
            do while(ray_energy > PPR*PPR_crit_perc/100d0)
!                !Move to Cell Edge
                x_ray = x_ray + pl*sin(theta)*cos(phi)
                y_ray = y_ray + pl*sin(theta)*sin(phi)
                z_ray = z_ray + pl*cos(theta)
!
!                !Absorb Fraction of Ray Energy to Cell
!!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
!!                constant absorption 
!!                call graygasconst(T_g(cellID),kg) 
!!                abs_coef = kg(1) 
                E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
                ray_energy = ray_energy - E_abs
                engra(ijkn) =engra(ijkn) + E_abs
!
!                !Check if Ray Hit a Surface
                if (abs(abs(x_ray)-x_max) <= eps .OR. &
                    abs(abs(y_ray)-y_max) <= eps .OR. &
                    abs(abs(z_ray)-z_max) <= eps .OR. &
                        abs(z_ray) <= eps)then
!
!                   !Absorb Fraction of Ray Energy to Surf
!                    surfID = GetsurfID(x_ray,y_ray,z_ray)
                    E_abs = em_surf*ray_energy 
!                   bc%gy(surfID,2) = bc%energy(surfID,2) &
!                                               + E_abs
                    ray_energy = ray_energy - E_abs
!
!                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf == 1d0)then
                        EXIT
                    end if
                    call random_number(rand)
                    call random_number(rand2)
!                    faceID = bc%face(surfID,1)
!                    reflected_direction = Reflect(faceID,rand,rand2)
!                    theta = reflected_direction(1)
!                    phi = reflected_direction(2)
                end if

!              find new cell 
                 if(x_ray.le.XEC(i)) then
                    inew =i
                  else
                    inew = i+1
                    if(inew.eq.imax2-1) inew = i
                 end if

                 if(y_ray.le.YNC(j)) then
                    jnew =j
                  else
                    jnew = j+1
                    if(jnew.eq.jmax2-1) jnew = j
                 end if

                 if(z_ray.le.ZTC(k)) then
                    knew = k
                  else
                    knew = k+1
                    if(knew.eq.kmax2-1) knew = k
                 end if

                ijkn = funijk(inew,jnew,knew) 

!                !Find New Distance to Next Cell Edge in Current Direction
                !D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
                call find_length(ijkn, x_ray,y_ray,z_ray, theta, phi,pl) 
!                
!                !Find New Cell ID
                x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
                y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
                z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
!                cellID =  GetcellID(x_mid,y_mid,z_mid)
            end do
!
!            !Absorb Final Small Fraction of Ray Energy if Below Threshold
            engra(ijkn) = engra(ijkn) + ray_energy
            ray_energy = 0d0
         end do 
        end do
    end do
  end do 
!
!    !Loop Over All Surfs
!    write(*,*)
!    do i = 1,N_surfs
!
!        !Get Center Coordinates and Face ID of Current Surf
!        x_c = coords((bc%glbnode(i,1)),1)
!        y_c = coords((bc%glbnode(i,1)),2)
!        z_c = coords((bc%glbnode(i,1)),3)
!        faceID = bc%face(i,1)
!
!        !Determine Number of Rays Current Cell Must Emit
!        N_rays = nint(bc%gy(i,1)*da(faceID)/PPR)
!        N_rays_tot = N_rays_tot + N_rays
!
!        !Absorb Roundoff Energy at Current Surf
!        bc%gy(i,2) = bc%energy(i,2) + bc%energy(i,1)*da(faceID)-N_rays*PPR  
!
!        !Emit All Rays from Current Surf
!        do j = 1,N_rays
!            !Display Current Progress
!            write(*,20,advance='no')(bsp,k=1,65),i, N_surfs, j, N_rays
!            20 format(65A1,'Emitting Rays from Surf: ',I5,'/',I5, &
!                '     Tracing Ray: ',I5,"/",I5)
!
!            !Pick Uniform Point of Emission
!            call random_number(rand)
!            x_ray = x_c-dx1/2d0+rand*dx1
!            call random_number(rand)
!            y_ray = y_c-dy1/2d0+rand*dy1
!            call random_number(rand)
!            z_ray = z_c-dz1/2d0+rand*dz1
!            if (abs(abs(x_c)-x_max/2D0) <= eps)then
!                x_ray = x_c
!            elseif(abs(abs(y_c)-y_max/2D0) <= eps)then
!                y_ray = y_c
!            elseif(abs(z_c-z_max) <= eps)then
!                z_ray = z_c
!            elseif(abs(z_c) <= eps)then
!                z_ray = z_c
!            end if
!
!            !Reset cellID
!            cellID = GetcellID(x_c,y_c,z_c)
!
!            !Pick Diffuse Direction of Emission
!            call random_number(rand)
!            call random_number(rand2)
!            faceID = bc%face(i,1)
!            reflected_direction = Reflect(faceID,rand,rand2)
!            theta = reflected_direction(1)
!            phi = reflected_direction(2)
!
!            !Find Distance to Cell Edge in Current Direction
!            D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
!
!            !Ray Tracing Until Ray Energy is Below Threshold Energy
!            ray_energy = PPR
!            do while(ray_energy > PPR*PPR_crit_perc/100d0)
!                !Move to Cell Edge
!                x_ray = x_ray + D_edge*sin(theta)*cos(phi)
!                y_ray = y_ray + D_edge*sin(theta)*sin(phi)
!                z_ray = z_ray + D_edge*cos(theta)
!
!                !Absorb Fraction of Ray Energy to Cell
!!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
!                E_abs = (1d0-exp(-abs_coef*D_edge))*ray_energy
!                ray_energy = ray_energy - E_abs
!                cell_data%gy(cellID,2) = cell_data%energy(cellID,2) + E_abs
!
!                !Check if Ray Hit a Surface
!                if (abs(abs(x_ray)-x_max/2d0) <= eps .OR. &
!                    abs(abs(y_ray)-y_max/2d0) <= eps .OR. &
!                    abs(z_ray-z_max) <= eps .OR. &
!                    abs(z_ray) <= eps)then
!
!                    !Absorb Fraction of Ray Energy to Surf
!                    surfID = GetsurfID(x_ray,y_ray,z_ray)
!                    E_abs = em_surf*ray_energy
!                    bc%gy(surfID,2) = bc%energy(surfID,2) &
!                                                + E_abs
!                    ray_energy  = ray_energy - E_abs
!
!                    !Reflect Diffusely (if Wall Emissivity /= 1)
!                    if (em_surf == 1d0)then
!                        EXIT
!                    end if
!                    call random_number(rand)
!                    call random_number(rand2)
!                    faceID = bc%face(surfID,1)
!                    reflected_direction = Reflect(faceID,rand,rand2)
!                    theta = reflected_direction(1)
!                    phi = reflected_direction(2)
!                end if
!
!                !Find New Distance to Next Cell Edge in Current Direction
!                D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
!                
!                !Find New Cell ID
!                x_mid = (x_ray + (x_ray + D_edge*sin(theta)*cos(phi)))/2d0
!                y_mid = (y_ray + (y_ray + D_edge*sin(theta)*sin(phi)))/2d0
!                z_mid = (z_ray + (z_ray + D_edge*cos(theta)))/2d0
!                cellID = GetcellID(x_mid,y_mid,z_mid)
!            end do
!
!            !Absorb Final Small Fraction of Ray Energy
!            cell_data%gy(cellID,2) = cell_data%energy(cellID,2) &
!                                        + ray_energy 
!            ray_energy  = 0d0
!        end do
!    end do
!
!    !Calculate Source Terms in Cells--------------------------------------------
!    do i = 1,N_cells
!       cell_data%gy(i,3) =  cell_data%energy(i,2)/dv - &
!                                cell_data%gy(i,1) 
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

write(*,*)
!contains
end subroutine mcfurnace 

!    !Current Cell ID From Coordinates Function----------------------------------
!    function GetcellID(x,y,z) result(cellID)
!        !This function returns the cell ID number that a point is in, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. If point is on cell face/edge/corner, cellID of the positive-side cell is returned (except for when on a positive surface). Cell ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z.
!        implicit none
!        !Initialize Variables
!        double precision :: x, y, z
!        integer :: cellID
!
!        !Check For Valid Coordinates
!        if (abs(x) > x_max/2d0)then
!            print *, "Invalid x Coordinate:", x
!            cellID = 0
!            return
!        elseif (abs(y) > y_max/2d0)then
!            print *, "Invalid y Coordinate:", y
!            cellID = 0
!            return
!        elseif (z > z_max.OR. z < 0d0)then
!            print *, "Invalid z Coordinate:", z
!            cellID = 0
!            return           
!        end if
!
!        !Adjust Coordinates if Point is on Positive Surface
!        if (abs(x - x_max/2d0) < eps)then
!            x = x - dx1/2d0
!        end if
!        if (abs(y - y_max/2d0) < eps)then
!            y = y - dy1/2d0
!        end if    
!        if (abs(z - z_max) < eps)then
!            z = z - dz1/2d0
!        end if
!
!        !Calculate Cell ID
!        cellID = 1+floor((x+x_max/2d0)/dx1)+floor((y+y_max/2d0)/dy1)*imax+ &
!                 floor(z/dz1)*imax*jmax
!    end function GetcellID
!
!    !Current Surface ID From Coordinates Function-------------------------------
!    function GetsurfID(x,y,z) result(surfID)
!        !This function returns the surface ID number that a point is on, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Note that if point is not exactly on a surface, it is snapped to the nearest surface. Surface ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z (on the current face), with faces ordered as -x,x,-y,y,-z,z.
!        implicit none
!        !Initialize Variables
!        double precision :: x, y, z
!        double precision, dimension(6) :: err
!        integer :: surfID
!
!        !Snap Coordinate to Nearest Surface
!        if (abs(x) /= x_max/2d0 .AND. abs(y) /= y_max/2d0 .AND. &
!            z /= z_max.AND. z /= 0d0)then
!            err(1) = abs(x+x_max/2d0)
!            err(2) = abs(x-x_max/2d0)
!            err(3) = abs(y+y_max/2d0)
!            err(4) = abs(y-y_max/2d0)
!            err(5) = abs(z)
!            err(6) = abs(z-z_max)
!            if (err(1) == minval(err))then
!                x = -x_max/2d0
!            elseif (err(2) == minval(err))then
!                x = x_max/2d0
!            elseif (err(3) == minval(err))then
!                y = -y_max/2d0
!            elseif (err(4) == minval(err))then
!                y = y_max/2d0
!            elseif (err(5) == minval(err))then
!                z = 0d0    
!            elseif (err(6) == minval(err))then
!                z = z_max
!            else
!                print *, "Surface Snap Error"
!            end if
!        end if
!
!        !Display Warning if Snap Distance was not Small
!        if (minval(err) >= eps)then
!            print *, "Warning: Surface Snap Distance >= eps"
!        end if
!
!        !Calculate Surface ID
!        if (x == -x_max/2d0) then
!            surfID = 1+floor((y+y_max/2d0)/dy1)+floor(z/dz1)*jmax
!        else if (x == x_max/2d0) then
!            surfID = 1+floor((y+y_max/2d0)/dy1)+floor(z/dz1)*jmax + jmax*kmax
!        else if (y == -y_max/2d0) then
!            surfID = 1+floor((x+x_max/2d0)/dx1)+floor(z/dz1)*imax + 2*jmax*kmax
!        else if (y == y_max/2d0) then
!            surfID = 1+floor((x+x_max/2d0)/dx1)+floor(z/dz1)*imax + 2*jmax*kmax &
!                    + imax*kmax
!        else if (z == 0d0) then
!            surfID = 1+floor((x+x_max/2)/dx1)+floor((y+y_max/2)/dy1)*imax &
!                    + 2*(jmax*kmax+imax*kmax)
!        else if (z == z_max) then
!            surfID = 1+floor((x+x_max/2)/dx1)+floor((y+y_max/2)/dy1)*imax &
!                    + 2*(jmax*kmax+imax*kmax) + imax*jmax
!        else 
!            print *, "Surface ID Error"
!            surfID = 0
!        end if
!    end function GetsurfID
!
    !Distance to Cell Edge in Current Direction Function------------------------
     subroutine find_length(ijk,x, y, z, theta, phi, pl) 
!    function Get_D_edge(x,y,z,theta,phi) result(D_edge)
        !This function computes the distance from the current point (defined by x,y,z) to the bounding cell surface in the current direction of emission (defined by input cone angle theta [0,PI] and azimuthal angle phi [0,2*PI]). Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. Point may be on a cell face/edge/corner. More information on page 158 of Mahan's book.
        implicit none
        !Initialize Variables
        double precision :: x,y,z,theta,phi
        double precision :: pl,x0,y0,z0,L,M,N,x_err,y_err,z_err
        double precision, dimension(6) :: t = -1
        integer :: i, cellID,ijk
        !Find Direction Cosines
        L = sin(theta)*cos(phi)
        M = sin(theta)*sin(phi)
        N = cos(theta)
        !Find Cell Corner
        cellID = ijk  !GetcellID(x,y,z)
        x0 = XEC(i_of(ijk)) 
        y0 = YNC(j_of(ijk)) 
        z0 = ZTC(k_of(ijk))
        write(*,*) "find_legth", ijk, x0, y0, z0
!        !Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
        x_err = abs(mod(x+mod(imax,2)*dx(i_of(ijk))/2d0,dx(i_of(ijk))) )
        y_err = abs(mod(y+mod(jmax,2)*dy(j_of(ijk))/2d0,dy(j_of(ijk))) )
        z_err = abs(mod(z,dz(k_of(ijk))))
        if (x_err <= eps .OR.dx(i_of(ijk))-x_err <= eps)then
            if (L > 0)then
                x0 = x
            else
                x0 = x - dx(i_of(ijk))
            end if
        end if
        if (y_err <= eps .OR.dy(j_of(ijk))-y_err <= eps)then
            if (M > 0)then
                y0 = y
            else
                y0 = y - dy(j_of(ijk))
            endif
        end if
        if (z_err <= eps .OR.dz(k_of(ijk))-z_err <= eps)then
            if (N > 0)then
                z0 = z
            else
                z0 = z - dz(k_of(ijk))
            end if
        end if
        !Find All t's
        if (abs(L) >= 1d-15)then 
            t(1) = (x0-x)/L
            t(2) = (x0+dx(i_of(ijk))-x)/L
        end if
        if (abs(M) >= 1d-15)then
            t(3) = (y0-y)/M
            t(4) = (y0+dy(j_of(ijk)) -y)/M
        end if
        if (abs(N) >= 1d-15)then          
            t(5) = (z0-z)/N
            t(6) = (z0+dz(k_of(ijk))-z)/N
        end if

        !Find Distance to Edge (Smallest, Positive t)
        pl = 1d15
        do i = 1,6
            if (t(i) > 0d0 .AND. t(i) < pl)then
                pl = t(i)
            end if
        end do
     return 
     end subroutine find_length 
!    end function Get_D_edge
!
!    !Reflecting Direction Function----------------------------------------------
!    function Reflect(faceID,rand,rand2) result(reflected_direction)
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
!
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
!
!        !Transform Reflected Direction (Vector) to Global Coordinate System
!        x_p_t = x_p*Q(1,1)+y_p*Q(1,2)+z_p*Q(1,3)
!        y_p_t = x_p*Q(2,1)+y_p*Q(2,2)+z_p*Q(2,3)
!        z_p_t = x_p*Q(3,1)+y_p*Q(3,2)+z_p*Q(3,3)
!
!        !Back Calculate Reflected Direction (Angles) Relative to GCS
!        reflected_direction(1) = atan2(sqrt(x_p_t**2+y_p_t**2),z_p_t) 
!        reflected_direction(2) = atan2(y_p_t,x_p_t)
!    end function Reflect
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

      subroutine setgeometry_pmc(XEC, YNC, ZTC) 
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
    double precision :: XEC (0:DIMENSION_I)
    double precision :: YNC (0:DIMENSION_J)
    double precision :: ZTC (0:DIMENSION_K)



      CALL GLOBAL_ALL_SUM(IER)

! Set boundary edges.
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout
      EX2 = X_MAX;    WX1 = X_MIN ! East/West
      TY2 = Y_MAX;    BY1 = Y_MIN ! North/South
      NZ2 = Z_MAX;    SZ1 = Z_MIN ! Top/Bottom

! Initialize arrays.
      XEC(:) = ZERO
      YNC(:) = ZERO
      ZTC(:) = ZERO

! Set boundary edges.
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout
      EX2 = X_MAX;    WX1 = X_MIN ! East/West
      TY2 = Y_MAX;    BY1 = Y_MIN ! North/South
      NZ2 = Z_MAX;    SZ1 = Z_MIN ! Top/Bottom

! Initialize arrays.
      XEC(:) = ZERO
      YNC(:) = ZERO
      ZTC(:) = ZERO

! Each loop starts at 2 and goes to max+2 (i.e., imin1=2, imax2=imax+2)
! However, the indices range to include ghost cells (0-imax2) to avoid
! multiple if statements in particles_in_cell
      XEC(IMIN2-1) = X_MIN-DX(IMIN2)
      DO I = IMIN2, IMAX2
         XEC(I) = XEC(I-1) + DX(I)
!         write(*,*) i, xec(i) 
      ENDDO

      YNC(JMIN2-1) = Y_MIN-DY(JMIN2)
      DO J  = JMIN2, JMAX2
         YNC(J) = YNC(J-1) + DY(J)
      ENDDO

      IF(DIMN.EQ.3) THEN
         ZTC(KMIN2-1) = Z_MIN-DZ(KMIN2)
         DO K = KMIN2, KMAX2
            ZTC(K) = ZTC(K-1) + DZ(K)
         ENDDO
      ENDIF

      RETURN
   END SUBROUTINE setgeometry_pmc

end module pmc_mod

