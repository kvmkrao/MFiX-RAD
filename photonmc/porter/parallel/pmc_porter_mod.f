module rad_pmc_mod 
! Module procedures.
     use mpi_utility
     use error_manager
     use param, only: dimension_3
     use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
     use geometry, only: imax, jmax, kmax
     use geometry, only: imax1, imin1, imin2, imax2, imax3, imin3
     use geometry, only: jmax1, jmin1, jmin2, jmax2, jmax3, jmin3
     use geometry, only: kmax1, kmin1, kmin2, kmax2, kmax3, kmin3
     use geometry, only: DX, DY, DZ 
     use geometry, only: x_max, x_min, y_max, y_min, z_max, z_min
     use geometry, only: VOL, AYZ, AXZ, AXY
     use compar 
     use functions, only: funijk, funijk_gl,is_on_mype_owns 
     use parallel_mpi
     use sendrecvnode
     double precision, dimension (:), allocatable  :: XEF  
     double precision, dimension (:), allocatable  :: YNF  
     double precision, dimension (:), allocatable  :: ZTF  

      double precision, dimension (:), allocatable :: T_g 
      double precision, dimension (:), allocatable :: emise,abse,srad 
      double precision, dimension (:), allocatable :: absebc 
! Global Parameters:
      double precision, parameter :: eps = 1d-12
      double precision, parameter :: PI = 4.d0*DATAN(1.d0)
      double precision, parameter :: stef_boltz = 5.670374d-8
      double precision :: dx1, dy1, dz1 
      integer, dimension(:), allocatable :: iprs, jprs, kprs
      integer, dimension(:), allocatable :: ipre, jpre, kpre

      integer, dimension(:), allocatable :: iprsf, jprsf, kprsf
      integer, dimension(:), allocatable :: ipref, jpref, kpref

contains 

subroutine rad_pmc_calc
!! Flag for 2D simulations.
    use geometry, only: NO_K
! Fixed array sizes in the I/J/K direction
    use bc, only: bc_jj_ps, bc_defined, bc_type_enum, par_slip_wall
    use bc, only: bc_k_b, bc_k_t, bc_j_s, bc_j_n, bc_i_w, bc_i_e
    use bc, only: free_slip_wall, no_slip_wall, dimension_bc
    use functions, only: east_of, west_of, top_of, bottom_of, north_of, south_of,fluid_at,is_on_mype_plus2layers 
    use desgrid, only: IofPROC, JofPROC, KofPROC
    use desgrid, only: procIJK

     implicit none    
    double precision, parameter :: em_surf = 1.0d0
    double precision, parameter :: abs_coef = 0.1d0

    !Define Working Variables-----------------------------------------------
    integer :: i,j,k,l,m,n, ijk, cellID, wlbcid
    integer :: ip,jp,kp,idiff, jdiff, kdiff 
    integer :: nrays, nrays_tot,ncountr 
    double precision, parameter :: PPR_crit_perc = 0.1d0
    double precision :: PPR, pl 
    double precision :: dv, x_c, y_c, z_c, x_ray, y_ray, z_ray
    double precision :: Tc, rad, rand, rand2, theta, phi, ray_energy, x_mid, y_mid, z_mid
    double precision :: D_edge, E_abs, E_in_sum, E_out_sum, time1, time2
    double precision :: E_out_sum1
    integer :: inew, jnew,knew, ijkn,ijkd

    double precision :: pres, xH2O, xCO2, volfrac,d_part, rho,eta,tmp
    integer :: npart_cell

    character(len=70) :: msg
    character(len=1) :: c, bsp = char(8)
    integer :: IER,OWNER,iproc,irequest  
    integer :: I1, I2, J1, J2, K1, K2,li,lj,lk, li2,lj2,lk2,lijkproc
    integer :: ic0, ic1, ic2
    integer :: liproc_start, liproc_end
    integer :: ljproc_start, ljproc_end
    integer :: lkproc_start, lkproc_end
    integer :: liproc, ljproc, lkproc,itotalneigh
    integer :: indx(4),indxr(4) 
    ! status variable - tells the status of send/ received calls
    ! Needed for receive subroutine
    integer, dimension(MPI_STATUS_SIZE) :: status1

    call CPU_TIME(time1)

    allocate(T_g(dimension_3))
    allocate(emise(dimension_3))
    allocate(abse(dimension_3))
    allocate(absebc(dimension_3))
    allocate(srad(dimension_3))
    allocate(XEF (0:DIMENSION_I), STAT=IER )
    allocate(YNF (0:DIMENSION_J), STAT=IER )
    allocate(ZTF (0:DIMENSION_K), STAT=IER )
    allocate(iprs(numPEs)) 
    allocate(jprs(numPEs)) 
    allocate(kprs(numPEs)) 
    allocate(ipre(numPEs)) 
    allocate(jpre(numPEs)) 
    allocate(kpre(numPEs)) 

    allocate(iprsf(numPEs)) 
    allocate(jprsf(numPEs)) 
    allocate(kprsf(numPEs)) 
    allocate(ipref(numPEs)) 
    allocate(jpref(numPEs)) 
    allocate(kpref(numPEs)) 

    call init_random_seed
    !Calculate Size of Homogeneous Subregions---------------------------------- 
    dx1 = x_max/float(imax)
    dy1 = y_max/float(jmax)
    dz1 = z_max/float(kmax)

    iprs = 0 
    ipre = 0 
    jprs = 0 
    jpre = 0 
    kprs = 0 
    kpre = 0 

    do iproc = 1,numPEs
      if(myPE.eq.iproc-1) then  
        iprs(iproc) = istart1 
        ipre(iproc) = iend1
        jprs(iproc) = jstart1
        jpre(iproc) = jend1
        kprs(iproc) = kstart1
        kpre(iproc) = kend1
      end if
    end do

    call MPI_ALLREDUCE(iprs,iprsf, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
    call MPI_ALLREDUCE(ipre,ipref, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
    call MPI_ALLREDUCE(jprs,jprsf, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
    call MPI_ALLREDUCE(jpre,jpref, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
    call MPI_ALLREDUCE(kprs,kprsf, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
    call MPI_ALLREDUCE(kpre,kpref, numPEs, MPI_INTEGER,MPI_SUM , MPI_COMM_WORLD,IER)
!   Read total number of photons to consider 
    open(21,file='input.dat') 
    read(21,*) nrays_tot
    close(21)

   !Location of cell centers 	
    do k=kstart1,kend1
      do j=jstart1,jend1
        do i=istart1,iend1
          XEF(i) = x_min + (i-2)*dx(i)  + dx(i)/2.0
          YNF(j) = y_min + (j-2)*dy(j)  + dy(j)/2.0
          ZTF(k) = z_min + (k-2)*dz(k)  + dz(k)/2.0
        end do 
      end do 
    end do 

    !Calculate Temperatures via cell centers 
    do k=kstart1, kend1
      do i=istart1, iend1
        do j=jstart1, jend1
           ijk =  funijk(i,j,k)  
           x_c = XEF(i)
           y_c = YNF(j) 
           z_c = ZTF(k)
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
        end do
      end do
    end do 

    E_out_sum1 = 0.0d0
    !Calculate Emissive Power at cell centers 
    !Calculate total Emissive Power
    do k = kstart1, kend1
      do i = istart1, iend1
        do j=jstart1, jend1 
          if(k.eq.2)      E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dx(i)*dy(j) 
          if(k.eq.kmax+1) E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dx(i)*dy(j) 
          if(j.eq.2)      E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dx(i)*dz(k) 
          if(j.eq.jmax+1) E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dx(i)*dz(k) 
          if(i.eq.2)      E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dy(j)*dz(k) 
          if(i.eq.imax+1) E_out_sum1 = E_out_sum1 + em_surf*stef_boltz*300.0d0**4*dy(j)*dz(k) 
          ijk =funijk(i,j,k)
          emise(ijk)  = 4.0d0*abs_coef*stef_boltz*(T_g(ijk))**4.0
          abse(ijk)   = 0.0d0
          absebc(ijk) = 0.0d0
          E_out_sum1   = E_out_sum1 + emise(ijk)*dx(i)*dy(j)*dz(k)
        end do 
      end do 
    end do 

    CALL MPI_ALLREDUCE(E_out_sum1, E_out_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD,IER)
    PPR = E_out_sum/float(nrays_tot)
    
!   !MONTE CARLO SIMULATION-----------------------------------------------------
!   !Loop Over All Cells
    ncountr = 0 
    do k =kstart1, kend1
      do i=istart1, iend1 
        do j=jstart1, jend1
          ijk = funijk(i,j,k) 
!!        !Determine Number of Rays to emit in the current cell 
          dv = dx(i)*dy(j)*dz(k) 
          nrays = nint(emise(ijk)*dv/PPR) 
          ncountr = ncountr + nrays

!         Absorb Energy due to rounding off the number of rays  at current cell
          abse(ijk) = abse(ijk) + emise(ijk)*dv-nrays*PPR  
          !loop over all Rays emit from current cell
          if(nrays.lt.1) cycle 
          do l = 1, nrays
!            !Display Current Progress
            !write(*,10,advance='no')(bsp,n=1,70),ijk,imax*jmax*kmax, l, nrays
            10 format(70A1,'Emitting Rays from Cell: ',I5,'/',I5, &
                '     Tracing Ray: ',I7,"/",I7)
            !Pick Uniform Point of Emission
            call random_number(rand)
!            rand = 0.5 
            x_ray = XEF(i)+rand*dx(i) - dx(i)/2.0
            call random_number(rand)
!            rand = 0.6 
            y_ray = YNF(j)+rand*dy(j) - dy(j)/2.0
            call random_number(rand)
!            rand = 0.7 
            z_ray = ZTF(k)+rand*dz(k) -dz(k)/2.0 

           !Pick Diffuse Direction of Emission
            call random_number(rand)
!            rand = 0.5
            theta = ACOS(1d0-2d0*rand)
            call random_number(rand)
!            rand = 0.6
            phi = 2d0*PI*rand
!!	    find new cell 
            call find_length(x_ray,y_ray,z_ray, theta, phi,pl) 
!            !Ray Tracing Until Ray Energy is Below Threshold Energy
            ray_energy  = PPR
            indx(1) = i 
            indx(2) = j 
            indx(3) = k
            owner   = myPE
            do while(ray_energy > PPR*PPR_crit_perc/100d0)
                x_ray = x_ray + pl*sin(theta)*cos(phi)
                y_ray = y_ray + pl*sin(theta)*sin(phi)
                z_ray = z_ray + pl*cos(theta)

                E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
                ray_energy = ray_energy - E_abs
                ! Establish the OWNER of the BC
                if(IS_ON_myPE_owns(indx(1),indx(2),indx(3))) then 
                   ijk = funijk(indx(1),indx(2),indx(3))  
                   abse(ijk) =abse(ijk) + E_abs
                else if(myPE.ne.owner) then 
               !call MPI_SENDRECV(indx,3,MPI_INTEGER, owner,1,indxr,3,MPI_INTEGER,myPE,1, MPI_COMM_WORLD,status1,IER) 
                   call MPI_ISEND(indx, 3, MPI_INTEGER, owner, 1, MPI_COMM_WORLD,irequest,IER)
                   call MPI_ISEND(E_abs, 1,MPI_DOUBLE_PRECISION, owner, 112, MPI_COMM_WORLD,irequest,IER)
                   if(myPE.eq.owner) then
                       call MPI_RECV(indxr, 3,MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,status1,IER)
                       call MPI_RECV(E_abs, 1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 112, MPI_COMM_WORLD,status1,IER)
                       ijkn = funijk(indxr(1),indxr(2),indxr(3))
                       abse(ijkn) =abse(ijkn) + E_abs
                   end if 
                end if 

                !Check if Ray Hit a Surface
                if (x_ray.le.x_min+eps .OR. x_ray.ge.x_max-eps  .OR. &
                    y_ray.le.y_min+eps .OR. y_ray.ge.y_max-eps  .OR. & 
                    z_ray.le.z_min+eps .OR. z_ray.ge.z_max-eps) then
!!                    
!!                  !Absorb Fraction of Ray Energy to Surf
                    call getwallbcnode(x_ray,y_ray,z_ray,indx(1),indx(2),indx(3),owner)
                    E_abs       = em_surf*ray_energy
                    if(IS_ON_myPE_owns(indx(1),indx(2),indx(3))) then 
                        wlbcid = funijk(indx(1),indx(2),indx(3))
                        absebc(wlbcid) =abse(wlbcid) + E_abs
                    else if (mype.ne.owner) then 
                        call MPI_ISEND(indx, 3, MPI_INTEGER, owner, 2, MPI_COMM_WORLD,irequest, IER)                                         
                        call MPI_ISEND(E_abs, 1,MPI_DOUBLE_PRECISION, owner, 113, MPI_COMM_WORLD,irequest,IER)
                        if(myPE.eq.owner) then
                           call MPI_RECV(indxr, 3,MPI_INTEGER, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,status1,IER)
                           call MPI_RECV(E_abs, 1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 113, MPI_COMM_WORLD,status1,IER)
                           wlbcid = funijk(indxr(1),indxr(2),indxr(3))
                           absebc(wlbcid) =absebc(wlbcid) + E_abs
                        end if
                    end if

                    ray_energy     = ray_energy - E_abs
!!                  
!!                  !Reflect Diffusely (if Wall Emissivity /= 1)
                    if (em_surf .eq. 1.0d0) exit
                    call random_number(rand)
                 !   rand =  0.5
                    call random_number(rand2)
                 !   rand =  0.6
                    call reflect(x_ray,y_ray,z_ray,rand,rand2,theta,phi)
                end if

                call find_length(x_ray,y_ray,z_ray, theta, phi,pl) 
!!              !Find New Cell ID
                x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
                y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
                z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
                call find_celln(x_mid, y_mid, z_mid,indx(1),indx(2),indx(3),owner) 
            end do 

!             Absorb Final Small Fraction of Ray Energy if Below Threshold
            ijk=funijk(i,j,k) 
            abse(ijk) = abse(ijk) + ray_energy
          end do   ! rays 
        end do   !j
      end do  ! i 
    end do  ! k 

    DO M = 1, DIMENSION_BC
      IF (BC_DEFINED(M)) THEN
        I1 = max(2,BC_I_W(M))
        I2 = min(BC_I_E(M),imax+1)
        J1 = max(2,BC_J_S(M))
        J2 = min(BC_J_N(M),jmax+1)
        K1 = max(2,BC_K_B(M))
        K2 = min(BC_K_T(M),kmax+1)
  
        DO K = K1, K2
           DO J = J1, J2
              DO I = I1, I2
                 ijk = funijk(i,j,k) 
                 if(m.le.2) then 
                      dv = dx1*dy1               
                 else if(m.ge.5) then 
                      dv = dy1*dz1    
                 else
                      dv = dx1*dz1 
                 end if            
  !!             !Determine Number of Rays to emit in the current cell 
                 nrays = nint(em_surf*stef_boltz*(300.0d0**4)*dv/PPR)
                 ncountr = ncountr + nrays
  
  !              Absorb Energy left due to round off at current cell
                 indx(1) = i
                 indx(2) = j
                 indx(3) = k
                 if(nrays.lt.1) cycle
                !loop over all Rays emit from current cell
                 do l = 1, nrays
  !               !Display Current Progress
                   !write(*,10,advance='no')bsp,ijk,imax*jmax*kmax, l, nrays
  !
  !                !Pick Uniform Point of Emission
                   call random_number(rand)
                 !  rand = 0.5
                   if(m.eq.5) then 
                   x_ray = 0.0+rand*dx(indx(1))    - dx(indx(1))/2.0 
                   else if (m.eq.6) then 
                   x_ray = x_max+rand*dx(indx(1))  - dx(indx(1))/2.0
                   else 
                   x_ray = XEF(i)+rand*dx(indx(1)) - dx(indx(1))/2.0
                   end if 
  
                   call random_number(rand)
                  ! rand = 0.6
                   if(m.eq.3) then
                   y_ray = 0.0+rand*dy(indx(2))   - dy(indx(2))/2.0 
                   else if (m.eq.4) then
                   y_ray = y_max+rand*dy(indx(2)) - dy(indx(2))/2.0
                   else
                   y_ray = YNF(j)+rand*dy(indx(2)) - dy(indx(2))/2.0
                   end if
  
                   call random_number(rand)
                  ! rand = 0.7
                   if(m.eq.1) then
                   z_ray = 0.0+rand*dz(indx(3))   - dz(indx(3))/2.0  
                   else if (m.eq.2) then
                   z_ray = z_max+rand*dz(indx(3)) - dz(indx(3))/2.0
                   else
                   z_ray = ZTF(k)+rand*dz(indx(3)) - dz(indx(3))/2.0
                   end if
  
  !                !Pick Diffuse Direction of Emission
                   call random_number(rand)
                   ! rand = 0.5 
                   theta = ACOS(1d0-2d0*rand)
                   call random_number(rand)
                   ! rand = 0.6 
                   phi = 2d0*PI*rand
  !!               find new cell 
                   call find_length(x_ray,y_ray,z_ray, theta, phi,pl)
  !                !Ray Tracing Until Ray Energy is Below Threshold Energy
                   ray_energy  = PPR
  
                   do while(ray_energy > PPR*PPR_crit_perc/100d0)
  !!                 !Move to Cell Edge
                     x_ray = x_ray + pl*sin(theta)*cos(phi)
                     y_ray = y_ray + pl*sin(theta)*sin(phi)
                     z_ray = z_ray + pl*cos(theta)
  
                     E_abs = (1d0-exp(-abs_coef*pl))*ray_energy
                     ray_energy = ray_energy - E_abs
                     ! Establish the OWNER of the BC
                     if(IS_ON_myPE_owns(indx(1),indx(2),indx(3))) then
                          ijk = funijk(indx(1),indx(2),indx(3))  
                          abse(ijk) =abse(ijk) + E_abs
                     else if(myPE.ne.owner) then 
                          call MPI_ISEND(indx, 3, MPI_INTEGER, owner, 1, MPI_COMM_WORLD,irequest,IER)
                          call MPI_ISEND(E_abs, 1,MPI_DOUBLE_PRECISION, owner, 112, MPI_COMM_WORLD,irequest,IER)
                         if(myPE.eq.owner) then 
                           call MPI_RECV(indxr, 3,MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,status1,IER)
                           call MPI_RECV(E_abs, 1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 112, MPI_COMM_WORLD,status1,IER)
                           ijkn = funijk(indxr(1),indxr(2),indxr(3))
                           abse(ijkn) =abse(ijkn) + E_abs
                         end if 
                     end if 
  
  !!                 !Check if Ray Hit a Surface
                     if (x_ray.le.x_min+eps .OR. x_ray.ge.x_max-eps  .OR. &
                         y_ray.le.y_min+eps .OR. y_ray.ge.y_max-eps  .OR. &
                         z_ray.le.z_min+eps .OR. z_ray.ge.z_max-eps) then
  !!                      
  !!                     !Absorb Fraction of Ray Energy to Surf
                         call getwallbcnode(x_ray,y_ray,z_ray,indx(1),indx(2),indx(3),owner)
                         E_abs          = em_surf*ray_energy
                         if(IS_ON_myPE_owns(indx(1),indx(2),indx(3))) then
                             wlbcid = funijk(indx(1),indx(2),indx(3))
                             absebc(wlbcid) =abse(wlbcid) + E_abs
                         else if (mype.ne.owner) then
                             call MPI_ISEND(indx, 3, MPI_INTEGER, owner, 1, MPI_COMM_WORLD,irequest,IER)
                             call MPI_ISEND(E_abs, 1,MPI_DOUBLE_PRECISION, owner, 112, MPI_COMM_WORLD,irequest,IER)
                             if(myPE.eq.owner) then
                                call MPI_RECV(indxr, 3,MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,status1,IER)
                                call MPI_RECV(E_abs, 1,MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 112, MPI_COMM_WORLD,status1,IER)
                                wlbcid = funijk(indxr(1),indxr(2),indxr(3))
                                absebc(wlbcid) =abse(wlbcid) + E_abs
                             end if
                         end if
                         ray_energy     = ray_energy - E_abs
  !!
  !!                     !Reflect Diffusely (if Wall Emissivity /= 1)
                         if (em_surf .eq. 1.0d0) exit
                         call random_number(rand)
                         !rand = 0.5 
                         call random_number(rand2)
                         !rand2 = 0.6 
                         call reflect(x_ray,y_ray,z_ray,rand,rand2,theta,phi)
                     end if
  
  !!                 !Find New Distance to Next Cell Edge in Current Direction
                     call find_length(x_ray,y_ray,z_ray, theta, phi,pl)
  !!                 !Find New Cell ID
                     x_mid = (x_ray + (x_ray + pl*sin(theta)*cos(phi)))/2d0
                     y_mid = (y_ray + (y_ray + pl*sin(theta)*sin(phi)))/2d0
                     z_mid = (z_ray + (z_ray + pl*cos(theta)))/2d0
                     call find_celln(x_mid, y_mid, z_mid,indx(1),indx(2),indx(3),owner) 
                   end do  ! ray energy
                   ijk=funijk(i,j,k) 
                   abse(ijk) = abse(ijk) + ray_energy 
                end do !nray loop 
             end do ! j 
          end do  ! i
        end do ! k 
        
      end if
    end do

    srad = 0.0d0
!    !Calculate Source Terms in Cells--------------------------------------------
    do k =kstart1, kend1
      do i=istart1, iend1 
        do j=jstart1, jend1
          ijk = funijk(i,j,k) 
          dv  = dx(i)*dy(j)*dz(k)
          srad(ijk) =  abse(ijk)/dv - emise(ijk)
          if(i.eq.11.and.j.eq.11) write(11,*) XEF(i),YNF(j),ZTF(k),T_g(ijk),srad(ijk)
          if(i.eq.10.and.j.eq.10) write(10,*) XEF(i),YNF(j),ZTF(k),T_g(ijk),srad(ijk)
          if(i.eq.9.and.j.eq.9) write(9,*) XEF(i),YNF(j),ZTF(k),T_g(ijk),srad(ijk)
        end do 
      end do 
    end do  

    !CPU Time
    call CPU_TIME(time2)
    write(*,*)
    write(*,"(a15,1f7.2,a4)") "Total CPU Time:", time2-time1, " (s)"

    deallocate(T_g)
    deallocate(emise)
    deallocate(abse)
    deallocate(absebc)
    deallocate(srad)
    deallocate(XEF )
    deallocate(YNF )
    deallocate(ZTF )
    deallocate(iprs)
    deallocate(jprs)
    deallocate(kprs)
    deallocate(ipre)
    deallocate(jpre)
    deallocate(kpre)

    deallocate(iprsf)
    deallocate(jprsf)
    deallocate(kprsf)
    deallocate(ipref)
    deallocate(jpref)
    deallocate(kpref)

    call MPI_Barrier(MPI_COMM_WORLD,IER)
    return 
end subroutine rad_pmc_calc

!subroutine find_celln(ray_x,ray_y,ray_z, ix, jy, kz,ijknew)
subroutine find_celln(ray_x,ray_y,ray_z, ix,jy,kz,iproc)
     integer :: ix, jy , kz,iproc
     integer :: ijkn,nproc
     integer :: inew, jnew, knew,i,j,k
     double precision :: ray_x, ray_y, ray_z 

    !Adjust Coordinates if Point is on Positive Surface
     if (abs(ray_x - x_max) < eps)then
        ray_x = ray_x - dx1/2d0
     end if
     if (abs(ray_y - y_max) < eps)then
        ray_y = ray_y - dy1/2d0
     end if
     if (abs(ray_z - z_max) < eps)then
        ray_z = ray_z - dz1/2d0
     end if

     ix = max(2,ceiling(ray_x/dx1)+1)
     jy = max(2,ceiling(ray_y/dy1)+1)
     kz = max(2,ceiling(ray_z/dz1)+1)
     if(ix.gt.imax+1) ix = imax+1
     if(jy.gt.jmax+1) jy = jmax+1
     if(kz.gt.kmax+1) kz = kmax+1
     if(IS_ON_myPE_owns(ix, jy, kz)) then 
              iproc  = myPE 
     else
       do nproc=1,numPEs
           if(ix.ge.iprsf(nproc).and.ix.le.ipref(nproc).and.& 
              jy.ge.jprsf(nproc).and.jy.le.jpref(nproc).and.& 
              kz.ge.kprsf(nproc).and.kz.le.kpref(nproc)) then 
              iproc   = nproc-1
           end if 
       end do 
     end if 
end subroutine find_celln 

!Subroutine to the id of the cell on a wall boundary-------------------------------
subroutine getwallbcnode(x,y,z,inew,jnew,knew,owner)
!This subroutine returns the id the cell on a wall boundary for the given the x,y, and z coordinates of the point.
        implicit none
!       !Initialize Variables
        double precision :: x, y, z
        integer :: noden,OWNER,nproc 
        integer :: inew,jnew,knew 

        inew = max(2,ceiling(x/dx1)+1)
        jnew = max(2,ceiling(y/dy1)+1)
        knew = max(2,ceiling(z/dz1)+1)

        if(x-eps .ge. x_max) inew = imax+1
        if(x .le. x_min+eps) inew = 2
        if(y-eps .ge. y_max) jnew = jmax+1
        if(y .le. y_min+eps) jnew = 2
        if(z-eps .ge. z_max) knew = kmax+1
        if(z .le. z_min+eps) knew = 2

        if(inew.gt.imax+1) inew = imax+1
        if(jnew.gt.jmax+1) jnew = jmax+1
        if(knew.gt.kmax+1) knew = kmax+1

        if(IS_ON_myPE_owns(inew, jnew, knew)) then 
          owner = myPE
        else 
          do nproc=1,numPEs 
            if(knew.ge.kprsf(nproc).and.knew.le.kpref(nproc)) then 
              if(inew.ge.iprsf(nproc).and.inew.le.ipref(nproc)) then 
                if(jnew.ge.jprsf(nproc).and.jnew.le.jpref(nproc)) then 
                   owner   = nproc-1
                end if 
              end if 
            end if 
          end do 
        end if 
end subroutine getwallbcnode 

!Subroutine to find distance to Cell Edge in Current Direction------------------------
subroutine find_length(x,y,z,theta, phi, pl) 
!This subroutine computes the distance from the current point (defined by x,y,z) to the bounding cell surface in the current direction of emission (defined by input cone angle theta [0,PI] and azimuthal angle phi [0,2*PI]). Point may be on a cell face/edge/corner. More information on page 158 of Mahan's book.
        implicit none
        !Initialize Variables
        double precision, intent(in)  :: x,y,z,theta,phi
        double precision, intent(out) :: pl
        integer :: i,j,k
        integer :: ix,jy,kz,ijknew
        double precision :: x0,y0,z0,dl,dm,dn,x_err,y_err,z_err
        !double precision :: dx,dy,dz
        double precision, dimension(6) :: t = -1 
        integer :: l, cellID,ijk
        !Find Direction Cosines
        dl = sin(theta)*cos(phi)
        dm = sin(theta)*sin(phi)
        dn = cos(theta)
        !Find Cell Corner
        ix = max(2,ceiling(x/dx1)+1)
        jy = max(2,ceiling(y/dy1)+1)
        kz = max(2,ceiling(z/dz1)+1)
        if(ix.gt.imax+1) ix = imax+1
        if(jy.gt.jmax+1) jy = jmax+1
        if(kz.gt.kmax+1) kz = kmax+1

!        write(*,*) "find length", ix, jy,kz,dx1,dy1,dz1
        x0 = XEF(ix) -dx1/2.0d0
        y0 = YNF(jy) -dy1/2.0d0  
        z0 = ZTF(kz) -dz1/2.0d0 
!        !Fix Corner if Emission Coordinate is on Cell Face, Edge, or Corner
        x_err = mod(x,dx1)
        y_err = mod(y,dy1)
        z_err = mod(z,dz1) 
        if (x_err .le. eps .OR.dx1-x_err .le. eps)then
            if (dl > 0)then
                x0 = x
            else
                x0 = x - dx1 
            end if
        end if
        if (y_err .le. eps .OR.dy1-y_err .le. eps) then
            if (dm > 0)then
                y0 = y
            else
                y0 = y - dy1  
            endif
        end if
        if (z_err .le. eps .OR.dz1-z_err .le. eps)then
            if (dn > 0)then
                z0 = z
            else
                z0 = z - dz1  
            end if
        end if
        !Find All t's
        if (abs(dl) >= eps)then 
            t(1) = (x0-x)/dl
            t(2) = (x0+dx1-x)/dl
        end if
        if (abs(dm) >= eps)then
            t(3) = (y0-y)/dm
            t(4) = (y0+dy1 -y)/dm
        end if
        if (abs(dn) >= eps)then          
            t(5) = (z0-z)/dn
            t(6) = (z0+dz1-z)/dn
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
!
!    !Subroutine to find Reflecting Direction----------------------------------------------
subroutine  reflect(x_r,y_r,z_r,rand,rand2, theta, phi)
!        !This subroutine provides the new ray direction after a reflecting event occurs, relative to the global cartesian coordinates given the faceID at which the reflection occurs(-x=1, +x=2, -y=3, +y=4, -z=5, +z=6) and two random numbers seeded between 0 and 1 (rand and rand2). The reflecting surface is assumed to be diffuse. (This function can also be utilized when emmitting rays from surfaces not normal to the global coordinate system.)
        implicit none
!        !Initialize Variables
        integer :: IDface
        integer :: i, j, k 
        double precision, intent(in) :: rand, rand2
        double precision :: theta, phi, theta_p, phi_p
        double precision ::  x_p, y_p, z_p, x_p_t, y_p_t, z_p_t
        double precision ::  x_r, y_r, z_r
        double precision, dimension(6,2) :: face_angles
        double precision, dimension(3,3) :: Q
        double precision :: dx, dy, dz 
        
!       Define Angles (Theta, Phi) of Face Normal Vectors
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
        return 
end subroutine reflect

SUBROUTINE init_random_seed
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER              :: isize,idate(8)
      INTEGER,ALLOCATABLE  :: iseed(:)
!-----------------------------------------------
      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize)
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed)
      iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)

      DEALLOCATE( iseed )

END SUBROUTINE init_random_seed

end module rad_pmc_mod
