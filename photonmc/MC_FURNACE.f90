!Gray Furnace Problem
!David Tobin - April 2020
!Details in Porter Paper

program MC_FURNACE
    implicit none
    !Initialize Input Parameters------------------------------------------------
    integer, parameter :: Nx = 17, Ny = 17, Nz = 34 
    double precision, parameter :: x_enc = 2d0, y_enc = 2d0, z_enc = 4d0, &
        em_surf = 1d0, abs_coef = 0.1d0

    !Initialize Working Variables-----------------------------------------------
    integer, parameter :: N_cells = Nx*Ny*Nz, N_surfs = 2*(Nx*Ny+Nx*Nz+Ny*Nz)
    integer :: i,j,k, cellID, surfID, faceID, N_rays, N_rays_tot = 0
    double precision, parameter :: PPR = abs_coef*10d0, epsilon = 1d-12, & 
        PPR_crit_perc = 0.1d0, PI = 4.d0*DATAN(1.d0), stef_boltz = 5.670374d-8 
    double precision :: dx, dy, dz, dv, x_c, y_c, z_c, x_ray, y_ray, z_ray, &
        Tc, rad, rand, rand2, theta, phi, ray_energy, x_mid, y_mid, z_mid, &
        D_edge, E_abs, E_in_sum = 0d0, E_out_sum = 0d0, time1, time2
    double precision, dimension (6) :: da
    double precision, dimension (2) :: reflected_direction
    character(len=1) :: c, bsp = char(8)
    character(len=50) :: filename, grid_t, abs_coef_t, PPR_t
    type cell_data_type
        integer,          dimension (N_cells,1) :: ID = 0
        double precision, dimension (N_cells,3) :: coords = 0
        double precision, dimension (N_cells,1) :: temp = 0
        double precision, dimension (N_cells,3) :: energy = 0
    end type
    type(cell_data_type) :: cell_data
    type surf_data_type
        integer,          dimension (N_surfs,1) :: ID = 0
        integer,          dimension (N_surfs,1) :: face = 0
        double precision, dimension (N_surfs,3) :: coords = 0
        double precision, dimension (N_surfs,1) :: temp = 0
        double precision, dimension (N_surfs,3) :: energy = 0
    end type
    type(surf_data_type) :: surf_data
    call CPU_TIME(time1)

    !Calculate Size of Homogeneous Subregions---------------------------------- 
    dx = x_enc/Nx
    dy = y_enc/Ny
    dz = z_enc/Nz
    dv = dx*dy*dz
    da = (/dy*dz,dy*dz,dx*dz,dx*dz,dx*dy,dx*dy/)

    !Calculate Cell Center Coordinates------------------------------------------
    do i = 1,Nz
        do j = 1, Ny
            do k = 1,Nx
                cellID = k+(j-1)*Nx+(i-1)*Nx*Ny
                cell_data%ID(cellID,1) = cellID
                cell_data%coords(cellID,1) = k*dx-x_enc/2d0-(dx/2d0)
                cell_data%coords(cellID,2) = j*dy-y_enc/2d0-(dy/2d0)
                cell_data%coords(cellID,3) = i*dz-(dz/2d0)
            end do
        end do
    end do

    !Calculate Surf Center Coordinates------------------------------------------
    !x-Faces
    do i = 1,Nz
        do j = 1,Ny
            do k = 1,2
                surfID = j+(i-1)*Ny+(k-1)*(Ny*Nz)
                surf_data%ID(surfID,1) = surfID
                surf_data%face(surfID,1) = k
                surf_data%coords(surfID,1) = x_enc/2d0*(-1)**k
                surf_data%coords(surfID,2) = j*dy-y_enc/2d0-(dy/2d0)
                surf_data%coords(surfID,3) = i*dz-(dz/2d0)
            end do
        end do
    end do
    !y-Faces
    do i = 1,Nz
        do j = 1,Nx
            do k = 1,2
                surfID = j+(i-1)*Nx+(Nx*Nz)*(k-1)+2*(Nz*Ny)
                surf_data%ID(surfID,1) = surfID
                surf_data%face(surfID,1) = k+2
                surf_data%coords(surfID,1) = j*dx-x_enc/2d0-(dx/2d0)
                surf_data%coords(surfID,2) = y_enc/2d0*(-1)**k
                surf_data%coords(surfID,3) = i*dz-(dz/2d0)
            end do
        end do
    end do
    !z-Faces
    do i = 1,Ny
        do j = 1,Nx
            do k = 1,2
                surfID = j+(i-1)*Nx+(Nx*Ny)*(k-1)+2*(Nz*Ny+Nx*Nz)
                surf_data%ID(surfID,1) = surfID
                surf_data%face(surfID,1) = k+4
                surf_data%coords(surfID,1) = j*dx-x_enc/2d0-(dx/2d0)
                surf_data%coords(surfID,2) = i*dy-y_enc/2d0-(dy/2d0)
                surf_data%coords(surfID,3) = z_enc*(k-1)
            end do
        end do
    end do

    !Calculate Cell Temperatures------------------------------------------------
    do i = 1,N_cells
        x_c = cell_data%coords(i,1)
        y_c = cell_data%coords(i,2)
        z_c = cell_data%coords(i,3)
        rad = sqrt(x_c**2 + y_c**2)
        if (rad < 1d0)then
            if (z_c <= 0.375d0)then
                Tc = 400d0+1400d0*(z_c/0.375d0)
            else
                Tc = 1800d0-1000d0*(z_c-0.375d0)/3.625d0
            end if
            cell_data%temp(i,1) = (Tc-800d0)*(1d0-3d0*rad**2+2d0*rad**3)+800d0
        else
            cell_data%temp(i,1) = 800d0
        end if
    end do

    !Calculate Surf Temperatures------------------------------------------------
    do i = 1,N_surfs
        surf_data%temp(i,1) = 300d0
    end do

    !Calculate Cell Emissive Powers---------------------------------------------
    do i = 1,N_cells
        cell_data%energy(i,1) = 4d0*abs_coef*stef_boltz*(cell_data%temp(i,1))**4
    end do

    !Calculate Surf Emissive Powers---------------------------------------------
    do i = 1,N_surfs
        surf_data%energy(i,1) = em_surf*stef_boltz*(surf_data%temp(i,1))**4
    end do

    !MONTE CARLO SIMULATION-----------------------------------------------------

    !Loop Over All Cells
    do i = 1,N_cells

        !Get Center Coordinates of Current Cell
        x_c = cell_data%coords(i,1)
        y_c = cell_data%coords(i,2)
        z_c = cell_data%coords(i,3)

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
                E_abs = (1d0-exp(-abs_coef*D_edge))*ray_energy
                ray_energy = ray_energy - E_abs
                cell_data%energy(cellID,2) = cell_data%energy(cellID,2) + E_abs

                !Check if Ray Hit a Surface
                if (abs(abs(x_ray)-x_enc/2d0) <= epsilon .OR. &
                    abs(abs(y_ray)-y_enc/2d0) <= epsilon .OR. &
                    abs(z_ray-z_enc) <= epsilon .OR. &
                    abs(z_ray) <= epsilon)then

                    !Absorb Fraction of Ray Energy to Surf
                    surfID = GetsurfID(x_ray,y_ray,z_ray)
                    E_abs = em_surf*ray_energy
                    surf_data%energy(surfID,2) = surf_data%energy(surfID,2) &
                                                + E_abs
                    ray_energy = ray_energy - E_abs

                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf == 1d0)then
                        EXIT
                    end if
                    call random_number(rand)
                    call random_number(rand2)
                    faceID = surf_data%face(surfID,1)
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
        x_c = surf_data%coords(i,1)
        y_c = surf_data%coords(i,2)
        z_c = surf_data%coords(i,3)
        faceID = surf_data%face(i,1)

        !Determine Number of Rays Current Cell Must Emit
        N_rays = nint(surf_data%energy(i,1)*da(faceID)/PPR)
        N_rays_tot = N_rays_tot + N_rays

        !Absorb Roundoff Energy at Current Surf
        surf_data%energy(i,2) = surf_data%energy(i,2) + &
                                surf_data%energy(i,1)*da(faceID)-N_rays*PPR  

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
            if (abs(abs(x_c)-x_enc/2D0) <= epsilon)then
                x_ray = x_c
            elseif(abs(abs(y_c)-y_enc/2D0) <= epsilon)then
                y_ray = y_c
            elseif(abs(z_c-z_enc) <= epsilon)then
                z_ray = z_c
            elseif(abs(z_c) <= epsilon)then
                z_ray = z_c
            end if

            !Reset cellID
            cellID = GetcellID(x_c,y_c,z_c)

            !Pick Diffuse Direction of Emission
            call random_number(rand)
            call random_number(rand2)
            faceID = surf_data%face(i,1)
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
                E_abs = (1d0-exp(-abs_coef*D_edge))*ray_energy
                ray_energy = ray_energy - E_abs
                cell_data%energy(cellID,2) = cell_data%energy(cellID,2) + E_abs

                !Check if Ray Hit a Surface
                if (abs(abs(x_ray)-x_enc/2d0) <= epsilon .OR. &
                    abs(abs(y_ray)-y_enc/2d0) <= epsilon .OR. &
                    abs(z_ray-z_enc) <= epsilon .OR. &
                    abs(z_ray) <= epsilon)then

                    !Absorb Fraction of Ray Energy to Surf
                    surfID = GetsurfID(x_ray,y_ray,z_ray)
                    E_abs = em_surf*ray_energy
                    surf_data%energy(surfID,2) = surf_data%energy(surfID,2) &
                                                + E_abs
                    ray_energy = ray_energy - E_abs

                    !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf == 1d0)then
                        EXIT
                    end if
                    call random_number(rand)
                    call random_number(rand2)
                    faceID = surf_data%face(surfID,1)
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
       faceID = surf_data%face(i,1)
       surf_data%energy(i,3) = surf_data%energy(i,2)/da(faceID) - &
                               surf_data%energy(i,1)
    end do

    !Save Centerline Data for Plots---------------------------------------------
    write(grid_t,'(i2,a1,i2,a1,i2)') nx,'x',ny,'x',nz
    write(abs_coef_t,'(1f5.3)') abs_coef
    write(PPR_t,'(1f5.3)') PPR
    filename = 'MC_FURNACE'//',Grid='//trim(grid_t)//',a=' &
                //trim(abs_coef_t)//',PPR='//trim(PPR_t)//'.dat'
    open (unit = 1, file = filename)
    do i = int(Nx*Ny/2d0+0.5d0),int(Nx*Ny*Nz-Nx*Ny/2d0+0.5d0),Nx*Ny
       write (1,"(3f9.3,1f12.3,1f12.3)")  cell_data%coords(i,1:3), &
                                          cell_data%temp(i,1), &
                                          cell_data%energy(i,3)
    end do

    !Display Final Results------------------------------------------------------

    !Energy Balance Check
    do i = 1,N_cells
        E_out_sum = E_out_sum + cell_data%energy(i,1)*dv 
        E_in_sum = E_in_sum + cell_data%energy(i,2)
    end do
    do i = 1,N_surfs
        faceID = surf_data%face(i,1)
        E_out_sum = E_out_sum + surf_data%energy(i,1)*da(faceID) 
        E_in_sum = E_in_sum + surf_data%energy(i,2)
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
                cell_data%ID(i,1), cell_data%coords(i,1:3), &
                cell_data%temp(i,1),cell_data%energy(i,1:3)
        end do
    end subroutine Display_Cell_Data

    !Display Surf Data Matrix Subroutine----------------------------------------
    subroutine Display_Surf_Data(min,max)
        !This subroutine displays a formatted view of the current data contained in the custom 'surf_data' matrix.
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
                surf_data%ID(i,1), surf_data%face(i,1), &
                surf_data%coords(i,1:3), surf_data%temp(i,1),&
                surf_data%energy(i,1:3)
        end do
    end subroutine Display_Surf_Data

    !Current Cell ID From Coordinates Function----------------------------------
    function GetcellID(x,y,z) result(cellID)
        !This function returns the cell ID number that a point is in, given the x,y, and z coordinates of the point. Function must also have access to grid geometry/coordinate varaibles. Point must be inside the enclosure or on a surface. If point is on cell face/edge/corner, cellID of the positive-side cell is returned (except for when on a positive surface). Cell ID numbers increase incrementally starting at the most negative x,y, and z cell, and increasing first in x, then in y, then in z.
        implicit none
        !Initialize Variables
        double precision :: x, y, z
        integer :: cellID

        !Check For Valid Coordinates
        if (abs(x) > x_enc/2d0)then
            print *, "Invalid x Coordinate:", x
            cellID = 0
            return
        elseif (abs(y) > y_enc/2d0)then
            print *, "Invalid y Coordinate:", y
            cellID = 0
            return
        elseif (z > z_enc .OR. z < 0d0)then
            print *, "Invalid z Coordinate:", z
            cellID = 0
            return           
        end if

        !Adjust Coordinates if Point is on Positive Surface
        if (abs(x - x_enc/2d0) < epsilon)then
            x = x - dx/2d0
        end if
        if (abs(y - y_enc/2d0) < epsilon)then
            y = y - dy/2d0
        end if    
        if (abs(z - z_enc) < epsilon)then
            z = z - dz/2d0
        end if

        !Calculate Cell ID
        cellID = 1+floor((x+x_enc/2d0)/dx)+floor((y+y_enc/2d0)/dy)*Nx+ &
                 floor(z/dz)*Nx*Ny
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
        if (abs(x) /= x_enc/2d0 .AND. abs(y) /= y_enc/2d0 .AND. &
            z /= z_enc .AND. z /= 0d0)then
            err(1) = abs(x+x_enc/2d0)
            err(2) = abs(x-x_enc/2d0)
            err(3) = abs(y+y_enc/2d0)
            err(4) = abs(y-y_enc/2d0)
            err(5) = abs(z)
            err(6) = abs(z-z_enc)
            if (err(1) == minval(err))then
                x = -x_enc/2d0
            elseif (err(2) == minval(err))then
                x = x_enc/2d0
            elseif (err(3) == minval(err))then
                y = -y_enc/2d0
            elseif (err(4) == minval(err))then
                y = y_enc/2d0
            elseif (err(5) == minval(err))then
                z = 0d0    
            elseif (err(6) == minval(err))then
                z = z_enc
            else
                print *, "Surface Snap Error"
            end if
        end if

        !Display Warning if Snap Distance was not Small
        if (minval(err) >= epsilon)then
            print *, "Warning: Surface Snap Distance >= epsilon"
        end if

        !Calculate Surface ID
        if (x == -x_enc/2d0) then
            surfID = 1+floor((y+y_enc/2d0)/dy)+floor(z/dz)*Ny
        else if (x == x_enc/2d0) then
            surfID = 1+floor((y+y_enc/2d0)/dy)+floor(z/dz)*Ny + Ny*Nz
        else if (y == -y_enc/2d0) then
            surfID = 1+floor((x+x_enc/2d0)/dx)+floor(z/dz)*Nx + 2*Ny*Nz
        else if (y == y_enc/2d0) then
            surfID = 1+floor((x+x_enc/2d0)/dx)+floor(z/dz)*Nx + 2*Ny*Nz &
                    + Nx*Nz
        else if (z == 0d0) then
            surfID = 1+floor((x+x_enc/2)/dx)+floor((y+y_enc/2)/dy)*Nx &
                    + 2*(Ny*Nz+Nx*Nz)
        else if (z == z_enc) then
            surfID = 1+floor((x+x_enc/2)/dx)+floor((y+y_enc/2)/dy)*Nx &
                    + 2*(Ny*Nz+Nx*Nz) + Nx*Ny
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
        x0 = cell_data%coords(cellID,1) - dx/2d0
        y0 = cell_data%coords(cellID,2) - dy/2d0
        z0 = cell_data%coords(cellID,3) - dz/2d0

        !Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
        x_err = abs(mod(x+mod(Nx,2)*dx/2d0,dx))
        y_err = abs(mod(y+mod(Ny,2)*dy/2d0,dy))
        z_err = abs(mod(z,dz))
        if (x_err <= epsilon .OR. dx-x_err <= epsilon)then
            if (L > 0)then
                x0 = x
            else
                x0 = x - dx
            end if
        end if
        if (y_err <= epsilon .OR. dy-y_err <= epsilon)then
            if (M > 0)then
                y0 = y
            else
                y0 = y - dy
            endif
        end if
        if (z_err <= epsilon .OR. dz-z_err <= epsilon)then
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

end program MC_FURNACE
