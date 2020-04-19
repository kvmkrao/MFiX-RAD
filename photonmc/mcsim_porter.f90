! photon monte carlo solver to find source terms due to thermal radiation 
! medium - gray gas 
!        - no particles 
! modified Apr 18, 2020 
program MC_FURNACE
    implicit none
    !Initialize Input Parameters------------------------------------------------
    integer, parameter :: imax = 17, jmax = 17, kmax = 34 
     double precision, parameter :: x_max= 2d0, y_max= 2d0, z_max= 4d0
     double precision, parameter :: em_surf = 1d0, abs_coef = 0.1d0

    !Initialize Working Variables-----------------------------------------------
    integer, parameter :: N_cells = imax*jmax*kmax, N_surfs = 2*(imax*jmax+imax*kmax+jmax*kmax)
    integer :: i,j,k, ijk, cellID, surfID, faceID, N_rays, N_rays_tot = 0
    double precision, parameter :: PPR = abs_coef*10d0, eps = 1d-12, & 
        PPR_crit_perc = 0.1d0, PI = 4.d0*DATAN(1.d0), stef_boltz = 5.670374d-8 
    double precision :: dx, dy, dz, dv, x_c, y_c, z_c, x_ray, y_ray, z_ray, &
        Tc, rad, rand, rand2, theta, phi, ray_energy, x_mid, y_mid, z_mid, &
        D_edge, E_abs, E_in_sum = 0d0, E_out_sum = 0d0, time1, time2
    double precision :: pres, xH2O, xCO2, volfrac,d_part, rho,eta 
    integer :: npart_cell
    double precision, dimension (6) :: da
    double precision, dimension (2) :: reflected_direction
    character(len=1) :: c, bsp = char(8)
    character(len=50) :: filename, grid_t, abs_coef_t, PPR_t
    double precision, dimension (N_cells) :: T_g
    double precision, dimension (N_cells,3) :: coords = 0
    type cell_data_type
        double precision, dimension (N_cells,3) :: energy = 0
    end type
    type(cell_data_type) :: cell_data
    type bc_type
        integer,          dimension (N_surfs,1) :: node = 0
        integer,          dimension (N_surfs,1) :: face = 0
        double precision, dimension (N_surfs,1) :: glbnode = 0
        double precision, dimension (N_surfs,3) :: energy = 0
    end type
    type(bc_type) :: bc
    call CPU_TIME(time1)

    !Calculate Size of Homogeneous Subregions---------------------------------- 
    dx = x_max/imax
    dy = y_max/jmax
    dz = z_max/kmax
    dv = dx*dy*dz
    da = (/dy*dz,dy*dz,dx*dz,dx*dz,dx*dy,dx*dy/)

    !Calculate Cell Center Coordinates------------------------------------------
    do k = 1,kmax
        do j = 1, jmax
            do i = 1,imax
                cellID = i+(j-1)*imax+(k-1)*imax*jmax
                coords(cellID,1) = i*dx-x_max/2d0-(dx/2d0)
                coords(cellID,2) = j*dy-y_max/2d0-(dy/2d0)
                coords(cellID,3) = k*dz-(dz/2d0)
            end do
        end do
    end do
    pres = 1.0d0 
    xH2O = 0.2d0 
    xCo2 = 0.1d0 
    volfrac = 1.0d0 
    d_part = 1.0d-3
    npart_cell = 100
    rho   = 1.0d0 
    !Calculate Cell Temperatures------------------------------------------------
    do i = 1,N_cells
        x_c = coords(i,1)
        y_c = coords(i,2)
        z_c = coords(i,3)
        rad = sqrt(x_c**2 + y_c**2)
        if (rad < 1d0)then
            if (z_c <= 0.375d0)then
                Tc = 400d0+1400d0*(z_c/0.375d0)
            else
                Tc = 1800d0-1000d0*(z_c-0.375d0)/3.625d0
            end if
            T_g(i) = (Tc-800d0)*(1d0-3d0*rad**2+2d0*rad**3)+800d0
        else
            T_g(i) = 800d0
        end if
    end do

    !Calculate Surf Center Coordinates------------------------------------------
    !x-Faces
    do k = 1,kmax
        do j = 1,jmax
            do i = 1,2
                surfID = j+(k-1)*jmax+(i-1)*(jmax*kmax)
                cellID = i+(j-1)*imax+(k-1)*imax*jmax
                if(i.eq.2)cellID = imax+(j-1)*imax+(k-1)*imax*jmax
                bc%node(surfID,1) = surfID
                bc%face(surfID,1) = i
                bc%glbnode(surfID,1) = cellID
                T_g(cellID) = 300.0d0 
            end do
        end do
    end do
    !y-Faces
    do k = 1,kmax
        do i = 1,imax
            do j = 1,2
                surfID = i+(k-1)*imax+(imax*kmax)*(j-1)+2*(kmax*jmax)
                cellID = i+(j-1)*imax+(k-1)*imax*jmax
                if(j.eq.2)cellID = i+(jmax-1)*imax+(k-1)*imax*jmax
                bc%node(surfID,1) = surfID
                bc%face(surfID,1) = j+2
                bc%glbnode(surfID,1) = cellID 
                T_g(cellID) = 300.0d0 
            end do
        end do
    end do
    !z-Faces
    do j = 1,jmax
        do i = 1,imax
            do k = 1,2
                surfID = i+(j-1)*imax+(imax*jmax)*(k-1)+2*(kmax*jmax+imax*kmax)
                cellID = i+(j-1)*imax+(k-1)*imax*jmax
                if(k.eq.2)cellID = i+(j-1)*imax+(kmax-1)*imax*jmax
                bc%node(surfID,1) = surfID
                bc%face(surfID,1) = k+4
                bc%glbnode(surfID,1) = cellID ! j*dx-x_enc/2d0-(dx/2d0)
                T_g(cellID) = 300.0d0
                if(k.eq.2) T_g(cellID) = 300.0d0
            end do
        end do
    end do

    !Calculate Cell Emissive Powers---------------------------------------------
    do i = 1,N_cells
!       eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!       abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)  
        cell_data%energy(i,1) = 4d0*abs_coef*stef_boltz*(T_g(i))**4
    end do

    !Calculate Surf Emissive Powers---------------------------------------------
    do i = 1,N_surfs
        bc%energy(i,1) = em_surf*stef_boltz*T_g(bc%glbnode(surfID,1))**4
    end do

    !MONTE CARLO SIMULATION-----------------------------------------------------

    !Loop Over All Cells
    do i = 1,N_cells

        !Get Center Coordinates of Current Cell
        x_c = coords(i,1)
        y_c = coords(i,2)
        z_c = coords(i,3)

        !Determine Number of Rays Current Cell Must Emit
        N_rays = nint(cell_data%energy(i,1)*dv/PPR)
        N_rays_tot = N_rays_tot + N_rays

        !Absorb Roundoff Energy at Current Cell
        cell_data%energy(i,2) = cell_data%energy(i,2) + &
                                cell_data%energy(i,1)*dv-N_rays*PPR  

        !Emit All Rays from Current Cell
        do j = 1,N_rays
            !Display Current Progress
            write(*,10,advance='no')(bsp,k=1,65),i, N_cells, j, N_rays
            10 format(65A1,'Emitting Rays from Cell: ',I5,'/',I5, &
                '     Tracing Ray: ',I5,"/",I5)

            !Pick Uniform Point of Emission
            call random_number(rand)
            x_ray = x_c-dx/2d0+rand*dx
            call random_number(rand)
            y_ray = y_c-dy/2d0+rand*dy
            call random_number(rand)
            z_ray = z_c-dz/2d0+rand*dz

            !Reset cellID
            cellID = i

            !Pick Diffuse Direction of Emission
            call random_number(rand) 
            theta = ACOS(1d0-2d0*rand)
            call random_number(rand)
            phi = 2d0*PI*rand  

            !Find Distance to Cell Edge in Current Direction
            D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)

            !Ray Tracing Until Ray Energy is Below Threshold Energy
            ray_energy = PPR
            do while(ray_energy > PPR*PPR_crit_perc/100d0)
                !Move to Cell Edge
                x_ray = x_ray + D_edge*sin(theta)*cos(phi)
                y_ray = y_ray + D_edge*sin(theta)*sin(phi)
                z_ray = z_ray + D_edge*cos(theta)

                !Absorb Fraction of Ray Energy to Cell
!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
                E_abs = (1d0-exp(-abs_coef*D_edge))*ray_energy
                ray_energy = ray_energy - E_abs
                cell_data%energy(cellID,2) = cell_data%energy(cellID,2) + E_abs

                !Check if Ray Hit a Surface
                if (abs(abs(x_ray)-x_max/2d0) <= eps .OR. &
                    abs(abs(y_ray)-y_max/2d0) <= eps .OR. &
                    abs(z_ray-z_max) <= eps .OR. &
                    abs(z_ray) <= eps)then

                    !Absorb Fraction of Ray Energy to Surf
                    surfID = GetsurfID(x_ray,y_ray,z_ray)
                    E_abs = em_surf*ray_energy
                    bc%energy(surfID,2) = bc%energy(surfID,2) &
                                                + E_abs
                    ray_energy = ray_energy - E_abs

                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf == 1d0)then
                        EXIT
                    end if
                    call random_number(rand)
                    call random_number(rand2)
                    faceID = bc%face(surfID,1)
                    reflected_direction = Reflect(faceID,rand,rand2)
                    theta = reflected_direction(1)
                    phi = reflected_direction(2)
                end if

                !Find New Distance to Next Cell Edge in Current Direction
                D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
                
                !Find New Cell ID
                x_mid = (x_ray + (x_ray + D_edge*sin(theta)*cos(phi)))/2d0
                y_mid = (y_ray + (y_ray + D_edge*sin(theta)*sin(phi)))/2d0
                z_mid = (z_ray + (z_ray + D_edge*cos(theta)))/2d0
                cellID = GetcellID(x_mid,y_mid,z_mid)
            end do

            !Absorb Final Small Fraction of Ray Energy if Below Threshold
            cell_data%energy(cellID,2) = cell_data%energy(cellID,2) &
                                        + ray_energy
            ray_energy = 0d0
        end do
    end do

    !Loop Over All Surfs
    write(*,*)
    do i = 1,N_surfs

        !Get Center Coordinates and Face ID of Current Surf
        x_c = coords((bc%glbnode(i,1)),1)
        y_c = coords((bc%glbnode(i,1)),2)
        z_c = coords((bc%glbnode(i,1)),3)
        faceID = bc%face(i,1)

        !Determine Number of Rays Current Cell Must Emit
        N_rays = nint(bc%energy(i,1)*da(faceID)/PPR)
        N_rays_tot = N_rays_tot + N_rays

        !Absorb Roundoff Energy at Current Surf
        bc%energy(i,2) = bc%energy(i,2) + bc%energy(i,1)*da(faceID)-N_rays*PPR  

        !Emit All Rays from Current Surf
        do j = 1,N_rays
            !Display Current Progress
            write(*,20,advance='no')(bsp,k=1,65),i, N_surfs, j, N_rays
            20 format(65A1,'Emitting Rays from Surf: ',I5,'/',I5, &
                '     Tracing Ray: ',I5,"/",I5)

            !Pick Uniform Point of Emission
            call random_number(rand)
            x_ray = x_c-dx/2d0+rand*dx
            call random_number(rand)
            y_ray = y_c-dy/2d0+rand*dy
            call random_number(rand)
            z_ray = z_c-dz/2d0+rand*dz
            if (abs(abs(x_c)-x_max/2D0) <= eps)then
                x_ray = x_c
            elseif(abs(abs(y_c)-y_max/2D0) <= eps)then
                y_ray = y_c
            elseif(abs(z_c-z_max) <= eps)then
                z_ray = z_c
            elseif(abs(z_c) <= eps)then
                z_ray = z_c
            end if

            !Reset cellID
            cellID = GetcellID(x_c,y_c,z_c)

            !Pick Diffuse Direction of Emission
            call random_number(rand)
            call random_number(rand2)
            faceID = bc%face(i,1)
            reflected_direction = Reflect(faceID,rand,rand2)
            theta = reflected_direction(1)
            phi = reflected_direction(2)

            !Find Distance to Cell Edge in Current Direction
            D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)

            !Ray Tracing Until Ray Energy is Below Threshold Energy
            ray_energy = PPR
            do while(ray_energy > PPR*PPR_crit_perc/100d0)
                !Move to Cell Edge
                x_ray = x_ray + D_edge*sin(theta)*cos(phi)
                y_ray = y_ray + D_edge*sin(theta)*sin(phi)
                z_ray = z_ray + D_edge*cos(theta)

                !Absorb Fraction of Ray Energy to Cell
!                eta = 400.0d0 + (1.0d0/(10000.0-400.0))*rand
!                abs_coef = Absorp_Coeff(eta,rho,T_g(cellID),pres,npart_cell,d_part)
                E_abs = (1d0-exp(-abs_coef*D_edge))*ray_energy
                ray_energy = ray_energy - E_abs
                cell_data%energy(cellID,2) = cell_data%energy(cellID,2) + E_abs

                !Check if Ray Hit a Surface
                if (abs(abs(x_ray)-x_max/2d0) <= eps .OR. &
                    abs(abs(y_ray)-y_max/2d0) <= eps .OR. &
                    abs(z_ray-z_max) <= eps .OR. &
                    abs(z_ray) <= eps)then

                    !Absorb Fraction of Ray Energy to Surf
                    surfID = GetsurfID(x_ray,y_ray,z_ray)
                    E_abs = em_surf*ray_energy
                    bc%energy(surfID,2) = bc%energy(surfID,2) &
                                                + E_abs
                    ray_energy = ray_energy - E_abs

                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf == 1d0)then
                        EXIT
                    end if
                    call random_number(rand)
                    call random_number(rand2)
                    faceID = bc%face(surfID,1)
                    reflected_direction = Reflect(faceID,rand,rand2)
                    theta = reflected_direction(1)
                    phi = reflected_direction(2)
                end if

                !Find New Distance to Next Cell Edge in Current Direction
                D_edge = Get_D_edge(x_ray,y_ray,z_ray,theta,phi)
                
                !Find New Cell ID
                x_mid = (x_ray + (x_ray + D_edge*sin(theta)*cos(phi)))/2d0
                y_mid = (y_ray + (y_ray + D_edge*sin(theta)*sin(phi)))/2d0
                z_mid = (z_ray + (z_ray + D_edge*cos(theta)))/2d0
                cellID = GetcellID(x_mid,y_mid,z_mid)
            end do

            !Absorb Final Small Fraction of Ray Energy
            cell_data%energy(cellID,2) = cell_data%energy(cellID,2) &
                                        + ray_energy
            ray_energy = 0d0
        end do
    end do

    !Calculate Source Terms in Cells--------------------------------------------
    do i = 1,N_cells
       cell_data%energy(i,3) =  cell_data%energy(i,2)/dv - &
                                cell_data%energy(i,1) 
    end do

    !Calculate Flux Terms in Surfs----------------------------------------------
    do i = 1,N_surfs
       faceID = bc%face(i,1)
       bc%energy(i,3) = bc%energy(i,2)/da(faceID) - &
                               bc%energy(i,1)
    end do

    !Save Centerline Data for Plots---------------------------------------------
    write(grid_t,'(i2,a1,i2,a1,i2)') imax,'x',jmax,'x',kmax
    write(abs_coef_t,'(1f5.3)') abs_coef
    write(PPR_t,'(1f5.3)') PPR
    filename = 'Porter'//'Grid_'//trim(grid_t)//'a_' &
                //trim(abs_coef_t)//'PPR_'//trim(PPR_t)//'.dat'
    open (unit = 1, file = filename)
    do i = int(imax*jmax/2d0+0.5d0),int(imax*jmax*kmax-imax*jmax/2d0+0.5d0),imax*jmax
       write (1,"(3f9.3,1f12.3,1f12.3)")  coords(i,1:3), &
                                          T_g(i), &
                                          cell_data%energy(i,3)
    end do

    !Display Final Results------------------------------------------------------

    !Energy Balance Check
    do i = 1,N_cells
        E_out_sum = E_out_sum + cell_data%energy(i,1)*dv 
        E_in_sum = E_in_sum + cell_data%energy(i,2)
    end do
    do i = 1,N_surfs
        faceID = bc%face(i,1)
        E_out_sum = E_out_sum + bc%energy(i,1)*da(faceID) 
        E_in_sum = E_in_sum + bc%energy(i,2)
    end do
    write(*,*) New_Line(c)
    write(*,"(a22,1f10.2,a4)") "Total Energy Emitted :", E_out_sum, " (W)"
    write(*,"(a22,1f10.2,a4)") "Total Energy Absorbed:", E_in_sum, " (W)"

    !Total Number of Rays Traced
    write(*,*)
    write(*,"(a28,1i9)") "Total Number of Rays Traced:", N_rays_tot

    !CPU Time
    call CPU_TIME(time2)
    write(*,*)
    write(*,"(a15,1f7.2,a4)") "Total CPU Time:", time2-time1, " (s)"

    !Saved File
    write(*,*)
    write(*,"(a30,a46)") "Saving Centerline Results to: ", trim(filename)

write(*,*)
contains

    !Display Cell Data Matrix Subroutine----------------------------------------
    subroutine Display_Cell_Data(min,max)
        !This subroutine displays a formatted view of the current data contained in the custom 'cell_data' matrix.
        implicit none
        !Initalize Variables
        integer :: min, max

        !Display Formatted Data
        do i = min,max
            if (mod(i-1,33) == 0)then
                print "(*(a1))", ("-", j=1,76)
                print "(a6,5x,a3,5x,a3,5x,a3,8x,a1,8x,a4,9x,a3,8x,a2)", &
                "CellID","x_c","y_c","z_c","T","Eout","Ein","dq"
                print "(*(a1))", ("-", j=1,76)
            end if
            print "(1x,1i4,2x,3f8.2,4x,1f6.1,3f11.1)", &
                coords(i,1:3), T_g(i),cell_data%energy(i,1:3)
!                cell_data%ID(i,1), coords(i,1:3), &
!                T_g(i),cell_data%energy(i,1:3)
        end do
    end subroutine Display_Cell_Data

    !Display Surf Data Matrix Subroutine----------------------------------------
    subroutine Display_bc(min,max)
        !This subroutine displays a formatted view of the current data contained in the custom 'bc' matrix.
        implicit none
        !Initalize Variables
        integer :: min, max

        !Display Formatted Data
        do i = min,max
            if (mod(i-1,33) == 0)then
                print "(*(a1))", ("-", j=1,76)
                print "(a6,3x,a4,4x,a3,5x,a3,5x,a3,6x,a1,7x,a4,7x,a3,8x,a2)", &
                "SurfID","Face","x_c","y_c","z_c","T","Eout","Ein","dq"
                print "(*(a1))", ("-", j=1,76)
            end if
            print "(1x,1i4,6x,1i1,1x,3f8.2,3x,1f5.1,3f10.1)", &
                bc%node(i,1), bc%face(i,1), &
                coords((bc%glbnode(i,1)),1:3), T_g(bc%glbnode(i,1)),&
                bc%energy(i,1:3)
        end do
    end subroutine Display_bc

    !Current Cell ID From Coordinates Function----------------------------------
    function GetcellID(x,y,z) result(cellID)
        !This function returns the cell ID number that a point is in, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. If point is on cell face/edge/corner, cellID of the positive-side cell is returned (except for when on a positive surface). Cell ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z.
        implicit none
        !Initialize Variables
        double precision :: x, y, z
        integer :: cellID

        !Check For Valid Coordinates
        if (abs(x) > x_max/2d0)then
            print *, "Invalid x Coordinate:", x
            cellID = 0
            return
        elseif (abs(y) > y_max/2d0)then
            print *, "Invalid y Coordinate:", y
            cellID = 0
            return
        elseif (z > z_max.OR. z < 0d0)then
            print *, "Invalid z Coordinate:", z
            cellID = 0
            return           
        end if

        !Adjust Coordinates if Point is on Positive Surface
        if (abs(x - x_max/2d0) < eps)then
            x = x - dx/2d0
        end if
        if (abs(y - y_max/2d0) < eps)then
            y = y - dy/2d0
        end if    
        if (abs(z - z_max) < eps)then
            z = z - dz/2d0
        end if

        !Calculate Cell ID
        cellID = 1+floor((x+x_max/2d0)/dx)+floor((y+y_max/2d0)/dy)*imax+ &
                 floor(z/dz)*imax*jmax
    end function GetcellID

    !Current Surface ID From Coordinates Function-------------------------------
    function GetsurfID(x,y,z) result(surfID)
        !This function returns the surface ID number that a point is on, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Note that if point is not exactly on a surface, it is snapped to the nearest surface. Surface ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z (on the current face), with faces ordered as -x,x,-y,y,-z,z.
        implicit none
        !Initialize Variables
        double precision :: x, y, z
        double precision, dimension(6) :: err
        integer :: surfID

        !Snap Coordinate to Nearest Surface
        if (abs(x) /= x_max/2d0 .AND. abs(y) /= y_max/2d0 .AND. &
            z /= z_max.AND. z /= 0d0)then
            err(1) = abs(x+x_max/2d0)
            err(2) = abs(x-x_max/2d0)
            err(3) = abs(y+y_max/2d0)
            err(4) = abs(y-y_max/2d0)
            err(5) = abs(z)
            err(6) = abs(z-z_max)
            if (err(1) == minval(err))then
                x = -x_max/2d0
            elseif (err(2) == minval(err))then
                x = x_max/2d0
            elseif (err(3) == minval(err))then
                y = -y_max/2d0
            elseif (err(4) == minval(err))then
                y = y_max/2d0
            elseif (err(5) == minval(err))then
                z = 0d0    
            elseif (err(6) == minval(err))then
                z = z_max
            else
                print *, "Surface Snap Error"
            end if
        end if

        !Display Warning if Snap Distance was not Small
        if (minval(err) >= eps)then
            print *, "Warning: Surface Snap Distance >= eps"
        end if

        !Calculate Surface ID
        if (x == -x_max/2d0) then
            surfID = 1+floor((y+y_max/2d0)/dy)+floor(z/dz)*jmax
        else if (x == x_max/2d0) then
            surfID = 1+floor((y+y_max/2d0)/dy)+floor(z/dz)*jmax + jmax*kmax
        else if (y == -y_max/2d0) then
            surfID = 1+floor((x+x_max/2d0)/dx)+floor(z/dz)*imax + 2*jmax*kmax
        else if (y == y_max/2d0) then
            surfID = 1+floor((x+x_max/2d0)/dx)+floor(z/dz)*imax + 2*jmax*kmax &
                    + imax*kmax
        else if (z == 0d0) then
            surfID = 1+floor((x+x_max/2)/dx)+floor((y+y_max/2)/dy)*imax &
                    + 2*(jmax*kmax+imax*kmax)
        else if (z == z_max) then
            surfID = 1+floor((x+x_max/2)/dx)+floor((y+y_max/2)/dy)*imax &
                    + 2*(jmax*kmax+imax*kmax) + imax*jmax
        else 
            print *, "Surface ID Error"
            surfID = 0
        end if
    end function GetsurfID

    !Distance to Cell Edge in Current Direction Function------------------------
    function Get_D_edge(x,y,z,theta,phi) result(D_edge)
        !This function computes the distance from the current point (defined by x,y,z) to the bounding cell surface in the current direction of emission (defined by input cone angle theta [0,PI] and azimuthal angle phi [0,2*PI]). Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. Point may be on a cell face/edge/corner. More information on page 158 of Mahan's book.
        implicit none
        !Initialize Variables
        double precision :: x,y,z,theta,phi
        double precision :: D_edge,x0,y0,z0,L,M,N,x_err,y_err,z_err
        double precision, dimension(6) :: t = -1
        integer :: i, cellID
        !Find Direction Cosines
        L = sin(theta)*cos(phi)
        M = sin(theta)*sin(phi)
        N = cos(theta)
        !Find Cell Corner
        cellID = GetcellID(x,y,z)
        x0 = coords(cellID,1) - dx/2d0
        y0 = coords(cellID,2) - dy/2d0
        z0 = coords(cellID,3) - dz/2d0
        !Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
        x_err = abs(mod(x+mod(imax,2)*dx/2d0,dx))
        y_err = abs(mod(y+mod(jmax,2)*dy/2d0,dy))
        z_err = abs(mod(z,dz))
        if (x_err <= eps .OR. dx-x_err <= eps)then
            if (L > 0)then
                x0 = x
            else
                x0 = x - dx
            end if
        end if
        if (y_err <= eps .OR. dy-y_err <= eps)then
            if (M > 0)then
                y0 = y
            else
                y0 = y - dy
            endif
        end if
        if (z_err <= eps .OR. dz-z_err <= eps)then
            if (N > 0)then
                z0 = z
            else
                z0 = z - dz
            end if
        end if
        !Find All t's
        if (abs(L) >= 1d-15)then 
            t(1) = (x0-x)/L
            t(2) = (x0+dx-x)/L
        end if
        if (abs(M) >= 1d-15)then
            t(3) = (y0-y)/M
            t(4) = (y0+dy-y)/M
        end if
        if (abs(N) >= 1d-15)then          
            t(5) = (z0-z)/N
            t(6) = (z0+dz-z)/N
        end if

        !Find Distance to Edge (Smallest, Positive t)
        D_edge = 1d15
        do i = 1,6
            if (t(i) > 0d0 .AND. t(i) < D_edge)then
                D_edge = t(i)
            end if
        end do
    end function Get_D_edge

    !Reflecting Direction Function----------------------------------------------
    function Reflect(faceID,rand,rand2) result(reflected_direction)
        !This function provides the new ray direction after a reflecting event occurs (reflected_direction = [theta',phi']), relative to the global cartesian coordinates, given the faceID at which the reflection occurs(-x=1, +x=2, -y=3, +y=4, -z=5, +z=6) and two random numbers seeded between 0 and 1 (rand and rand2). The reflecting surface is assumed to be diffuse. (This function can also be utilized when emmitting rays from surfaces not normal to the global coordinate system.)
        implicit none
        !Initialize Variables
        integer, intent(in) :: faceID
        double precision, intent(in) :: rand, rand2
        double precision :: theta, phi, theta_p, phi_p, &
        x_p, y_p, z_p, x_p_t, y_p_t, z_p_t
        double precision, dimension(2) :: reflected_direction
        double precision, dimension(6,2) :: face_angles
        double precision, dimension(3,3) :: Q

        !Define Angles (Theta, Phi) of Face Normal Vectors
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

        !Get Angles of Current Face Normal Vector
        theta = face_angles(faceID,1)
        phi = face_angles(faceID,2)

        !Calculate Reflected Direction (Angles) Relative to Face Normal
        theta_p = asin(sqrt(rand))
        phi_p = 2D0*PI*rand2

        !Calculate Reflected Direction (Vector) Relative to Face Normal
        x_p = sin(theta_p)*cos(phi_p)
        y_p = sin(theta_p)*sin(phi_p)
        z_p = cos(theta_p)

        !Define the Transformation Matrix
        Q(1,1) = cos(phi)*cos(theta)
        Q(1,2) = -sin(phi)
        Q(1,3) = sin(theta)*cos(phi)
        Q(2,1) = sin(phi)*cos(theta)
        Q(2,2) = cos(phi)
        Q(2,3) = sin(theta)*sin(phi)
        Q(3,1) = -sin(theta)
        Q(3,2) = 0
        Q(3,3) = cos(theta)

        !Transform Reflected Direction (Vector) to Global Coordinate System
        x_p_t = x_p*Q(1,1)+y_p*Q(1,2)+z_p*Q(1,3)
        y_p_t = x_p*Q(2,1)+y_p*Q(2,2)+z_p*Q(2,3)
        z_p_t = x_p*Q(3,1)+y_p*Q(3,2)+z_p*Q(3,3)

        !Back Calculate Reflected Direction (Angles) Relative to GCS
        reflected_direction(1) = atan2(sqrt(x_p_t**2+y_p_t**2),z_p_t) 
        reflected_direction(2) = atan2(y_p_t,x_p_t)
    end function Reflect
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
integer :: i!,iq
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


end program MC_FURNACE
