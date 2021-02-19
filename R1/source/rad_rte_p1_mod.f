!!****************************************************************
!!** module rte_p1
!!****************************************************************
module rad_rte_p1
! module to solve the P1 Radiative Transfer Equation (RTE)
use rad_param
use rad_fields
use SOLVE_LIN_EQ_MOD
implicit none
private

! ---- Interfaces ------
public :: rad_rte_p1_init
public :: rad_rte_p1_final
public :: rad_rte_p1_calc

interface rad_rte_p1_init
    module procedure init1
end interface

interface rad_rte_p1_final
    module procedure final1
end interface

interface rad_rte_p1_calc
    module procedure calc
end interface

contains

subroutine init1(config)
    use rad_config, only : configuration
    implicit none
    type(configuration), intent(in) :: config
end subroutine

subroutine final1()
    implicit none
end subroutine

subroutine calc()
use rad_util, only : bbem

    real(dp), dimension(dimension_3) :: k, beta, w, omega, s, Em
    real(dp), dimension(dimension_bc) :: emw
    integer :: iq, m

    do iq = 1, nq
        k = k_g(:, iq) + sum(k_s(:, :, iq), 2)
        k = max(k, BETAMIN)
        beta = k + sum(scat(:,:, iq),2)
        omega = sum(scat(:,:,iq),2)/beta
        omega = min(omega,1.d0)
        w = 3*(1.d0-omega)
        Em = E(:,0,iq)*k_g(:,iq)
        do m=1, smax
            Em=Em+E(:,m,iq)*k_s(:,m,iq)
        end do
        Em = Em/k*4.d0*pi
        call bbem (emw, Tw)
        emw=emw*4*pi
        call solve_p1_eqn(G(:,iq), beta, w, omega, Em, ew, emw)
    end do
end subroutine
!------------------------------------------------------------------
! subroutine to solve generic P1-like PDE with boundary conditions
subroutine solve_p1_eqn(G, beta, w, omega, S, ew, Sw)
!use param
!use param1
!use ambm
!use geometry
!use indices
!use compar
!use run
!use leqsol
!use bc
!use toleranc
!use output
!use residual
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use bc, only: bc_x_g, bc_xw_g, bc_hw_x_g, bc_c_x_g
      use bc, only: bc_x_s, bc_xw_s, bc_hw_x_s, bc_c_x_s
      use bc, only: bc_defined
      use ChiScheme, only: set_chi, unset_chi
      use compar, only: ijkstart3, ijkend3, mype, numpes
      use fldvar, only: ep_g, u_g, v_g, w_g, x_g, x_go, rop_go
      use fldvar, only: ep_s, u_s, v_s, w_s, x_s, x_so, rop_so
      use functions, only: fluid_at, zmax
      use geometry, only: ijkmax2, vol
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use mflux, only: flux_ge, flux_gse, flux_gn, flux_gsn
      use mflux, only: flux_se, flux_sse, flux_sn, flux_ssn
      use mflux, only: flux_gt, flux_gst, flux_st, flux_sst
      use mpi_utility, only: global_all_sum
      use param, only: dimension_3, dimension_m
      use param, only: dimension_n_s, dimension_n_g
      use param1, only: zero, one, undefined
      use physprop, only: nmax, dif_g, dif_s
      use physprop, only: smax
      use ps, only: point_source
      use ps, only: ps_x_g, ps_massflow_g
      use ps, only: ps_x_s, ps_massflow_s
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_x
      use run, only: species_eq, discretize, odt, added_mass, m_am
      use run, only: chi_scheme
      use rxns, only: sum_r_g, rox_gc, r_gp
      use rxns, only: sum_r_s, rox_sc, r_sp
      use toleranc, only: zero_x_gs
      use ur_facs, only: ur_fac
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: gas_species, solids_species
      use utilities, only: bound_x
    ! input and output
    real(dp), intent(out), dimension(DIMENSION_3) :: G
    ! G: incident radiation
    real(dp), intent(in), dimension(DIMENSION_3) :: beta, w, S, omega
    ! beta: extinction coefficient
    ! omega: scattering albedo
    ! w: general P1 equation proportional constant 
    ! S: source term in general P1 equation
    real(dp), intent(in), dimension(DIMENSION_BC) :: ew, Sw
    ! ew: wall emissivity
    ! Sw: wall emission source 
    
    ! local variables
    real(dp), dimension(DIMENSION_3) :: ZEROS,ONES, invBeta
    real(dp), dimension(DIMENSION_3) :: S_C, S_P
    real(dp), dimension(DIMENSION_BC) :: BC_HW, BC_C, BC_W, BC_F
    ! Indices
    INTEGER IJK, I, J, K, L, TOTAL_BC
    ! linear equation solver method and iterations
    INTEGER LEQM, LEQI
    ! temporary variables in residual computation 
    real(dp) :: res1, mres1, num_res, den_res
    INTEGER ires1 , ier
    ! initialize variables
    !G=ZERO
    ZEROS = ZERO 
    ONES = ONE 
	invBeta = one/(max(beta, BETAMIN))
    
    CALL LOCK_AMBM
    CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER)

!C-------------------------------------------------------------
!C  discretize the diffusion term
!c-------------------------------------------------------------
    ! it is assumed that extrapolation has been done 
    ! before entering this subroutine.
    ! variables need extrapolation to ghost cells 
    ! include beta, w, and omega.

    CALL CONV_DIF_G (G, invBeta, ZEROS, ZEROS, ZEROS, &
        ZEROS, ZEROS, ZEROS, 0, A_M)

!C-------------------------------------------------------------
!C  discretize boundary conditions
!c-------------------------------------------------------------
    BC_HW = ZERO  ! default values
    BC_C  = ZERO

	do L = 1,DIMENSION_BC
		if (.not.BC_DEFINED(L)) then
			BC_C(L) = UNDEFINED
			BC_HW(L) = UNDEFINED
		else
			BC_C(L) = ew(L)/(2d0-ew(L))*1.5d0*Sw(L)
			BC_HW(L) =  ew(L)/(2d0-ew(L))*1.5d0
		end if
	end do
	BC_W = ZERO
	BC_F = ZERO

    CALL BC_G (G, BC_F, BC_W, BC_HW, BC_C, 0, A_M, B_M, beta)
	
!C-------------------------------------------------------------
!C  discretize the source term
!c-------------------------------------------------------------
	S_C = w*beta*VOL*S 
	S_P = w*beta*VOL
	
    CALL SOURCE_PHI (S_P, S_C, ONES, G, 0, A_M, B_M, IER)

!C-------------------------------------------------------------
!C  solve the linear equation
!c-------------------------------------------------------------
    CALL CALC_RESID_S (G, A_M, B_M, 0, num_res, den_res, res1, &
        mres1, ires1, ZERO, IER) 
    ! P1 equation does not require under relaxation
    ! CALL UNDER_RELAX_S (G, A_M, B_M, 0, UR_FAC(9), IER) 
    CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), LEQI, LEQM, IER)
    CALL SOLVE_LIN_EQ ('G', 9, G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER) 

    CALL UNLOCK_AMBM
    return
end subroutine solve_P1_eqn

!------------------- local procedures --------------------------------//

! convection-diffusion discretization for the P1 equation.
! because the P1 equation does not have convection terms, 
! the convection discretization scheme is not effective.
! Here 'conv_dif_phi0' is used as a base subroutine


      SUBROUTINE CONV_DIF_G(PHI, DIF, UF, VF, WF, &
                               Flux_E, Flux_N, Flux_T, M, A_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: fluid_at
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE param
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Scalar
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: Uf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Vf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Wf(DIMENSION_3)
! Mass flux components
      DOUBLE PRECISION, INTENT(IN) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_T(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb
!---------------------------------------------------------------------//


!!!$omp      parallel do                                              &
!!!$omp&     private(IJK, IPJK, IJPK, IJKM, IMJK, IJMK, IJKM,         &
!!!$omp&             d_fe, df_w, df_n, df_s, df_t, df_b)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            !-- rad_changes_begin
            ! need to change the interpolation of the diffusivity (1/beta)
            ! near walls.
            CALL GET_GCELL_DIFF_TERMS(dif, d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)
            ! CALL GET_PHICELL_DIFF_TERMS(dif, d_fe, d_fw, d_fn, d_fs, &
            !   d_ft, d_fb, ijk)
            !-- rad_changes_end

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            IF (UF(IJK) >= ZERO) THEN
               A_M(IJK,east,M) = D_Fe
               A_M(IPJK,west,M) = D_Fe + FLUX_E(IJK)
            ELSE
               A_M(IJK,east,M) = D_Fe - FLUX_E(IJK)
               A_M(IPJK,west,M) = D_Fe
            ENDIF
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               IF (UF(IMJK) >= ZERO) THEN
                  A_M(IJK,west,M) = D_Fw + FLUX_E(IMJK)
               ELSE
                  A_M(IJK,west,M) = D_Fw
               ENDIF
            ENDIF


! North face (i, j+1/2, k)
            IF (VF(IJK) >= ZERO) THEN
               A_M(IJK,north,M) = D_Fn
               A_M(IJPK,south,M) = D_Fn + FLUX_N(IJK)
            ELSE
               A_M(IJK,north,M) = D_Fn - FLUX_N(IJK)
               A_M(IJPK,south,M) = D_Fn
            ENDIF
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               IF (VF(IJMK) >= ZERO) THEN
                  A_M(IJK,south,M) = D_Fs + FLUX_N(IJMK)
               ELSE
                  A_M(IJK,south,M) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)

! Top face (i, j, k+1/2)
               IF (WF(IJK) >= ZERO) THEN
                  A_M(IJK,top,M) = D_FT
                  A_M(IJKP,bottom,M) = D_Ft + FLUX_T(IJK)
               ELSE
                  A_M(IJK,top,M) = D_Ft - FLUX_T(IJK)
                  A_M(IJKP,bottom,M) = D_Ft
               ENDIF

! Bottom face (i, j, k-1/2)

               IF (.NOT.FLUID_AT(IJKM)) THEN
                  IF (WF(IJKM) >= ZERO) THEN
                     A_M(IJK,bottom,M) = D_Fb + FLUX_T(IJKM)
                  ELSE
                     A_M(IJK,bottom,M) = D_Fb
                  ENDIF
               ENDIF
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_G

!-- rad_changes_begin
! subroutine for cell diffusion coefficients
! based on 'GET_PHICELL_DIFF_TERMS' 
! the indices near the wall is changed
      SUBROUTINE GET_GCELL_DIFF_TERMS(Dif, D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK)
!-- rad_changes_end

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_treatment_at, cut_cell_at

      USE functions, only: fluid_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of
      USE functions, only: ip_of, jp_of, kp_of

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE geometry, only: odx_e, ody_n, odz_t
      USE geometry, only: do_k
      USE geometry, only: ox
      USE geometry, only: dx, dy, dz
      USE geometry, only: ayz, axz, axy

      USE indices, only: i_of, j_of, k_of
      USE indices, only: im1, jm1, km1

      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)

! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ipjk, ijpk, ijkp
      INTEGER :: i, j, k, im, jm, km
      INTEGER :: ijke, ijkw, ijkn, ijks, ijkt, ijkb

! area terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
         
      !-- rad_changes_begin
      ! the following three indices were only for cut-cell method before
      ! they are added here 
      IPJK = IP_OF(IJK)
      IJPK = JP_OF(IJK)
      IJKP = KP_OF(IJK)
      !-- rad_changes_end  

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = IM1(I)
      JM = JM1(J)
      KM = KM1(K)

      !-- rad_changes_begin
      ! these indices are not needed
      ! IJKE = EAST_OF(IJK)
      ! IJKN = NORTH_OF(IJK)
      ! IJKS = SOUTH_OF(IJK)
      ! IJKW = WEST_OF(IJK)
      !-- rad_changes_end  

      C_AE = ODX_E(I)*AYZ(IJK)
      C_AW = ODX_E(IM)*AYZ(IMJK)
      C_AN = ODY_N(J)*AXZ(IJK)
      C_AS = ODY_N(JM)*AXZ(IJMK)
      C_AT = OX(I)*ODZ_T(K)*AXY(IJK)
      C_AB = OX(I)*ODZ_T(KM)*AXY(IJKM)

      IF(CUT_TREATMENT_AT(IJK).AND.CUT_CELL_AT(IJK)) THEN
         !-- rad_changes_begin
         ! the following three indices were only for cut-cell method before
         ! they have been added above
         ! IPJK = IP_OF(IJK)
         ! IJPK = JP_OF(IJK)
         ! IJKP = KP_OF(IJK)
         !-- rad_changes_end
            
         IF (.NOT.FLUID_AT(IPJK)) C_AE = ODX_E(I)*DY(J)*DZ(K)
         IF (.NOT.FLUID_AT(IMJK)) C_AW = ODX_E(IM)*DY(J)*DZ(K)
         IF (.NOT.FLUID_AT(IJPK)) C_AN = ODY_N(J)*DX(I)*DZ(K)
         IF (.NOT.FLUID_AT(IJMK)) C_AS = ODY_N(JM)*DX(I)*DZ(K)
         IF (.NOT.FLUID_AT(IJKP)) C_AT = OX(I)*ODZ_T(K)*DX(I)*DY(J)
         IF (.NOT.FLUID_AT(IJKM)) C_AB = OX(I)*ODZ_T(KM)*DX(I)*DY(J)
      ENDIF

!-- rad_changes_begin
! the old calculations are commented out
!! East face (i+1/2, j, k)
!      D_Fe = AVG_X_H(DIF(IJK),DIF(IJKE),I)*C_AE
!
!! West face (i-1/1, j, k)
!      D_FW = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*C_AW
!
!! North face (i, j+1/2, k)
!      D_FN = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*C_AN
!
!! South face (i, j-1/2, k)
!      D_FS = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*C_AS
!      IF (DO_K) THEN
!         IJKT = TOP_OF(IJK)
!         IJKB = BOTTOM_OF(IJK)
!
!! Top face (i, j, k+1/2)
!         D_FT = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*C_AT
!! Bottom face (i, j, k-1/2)
!         D_FB = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*C_AB
!      ENDIF


! new calculations using different indices near wall
! East face (i+1/2, j, k)
      D_Fe = AVG_X_H(DIF(IJK),DIF(IPJK),I)*C_AE

! West face (i-1/1, j, k)
      D_FW = AVG_X_H(DIF(IMJK),DIF(IJK),IM)*C_AW

! North face (i, j+1/2, k)
      D_FN = AVG_Y_H(DIF(IJK),DIF(IJPK),J)*C_AN

! South face (i, j-1/2, k)
      D_FS = AVG_Y_H(DIF(IJMK),DIF(IJK),JM)*C_AS

      IF (DO_K) THEN
! Top face (i, j, k+1/2)
         D_FT = AVG_Z_H(DIF(IJK),DIF(IJKP),K)*C_AT
! Bottom face (i, j, k-1/2)
         D_FB = AVG_Z_H(DIF(IJKM),DIF(IJK),KM)*C_AB
      ENDIF
!-- rad_changes_end

      RETURN
      END SUBROUTINE GET_GCELL_DIFF_TERMS

      SUBROUTINE BC_G(VAR, BC_PHIF, BC_PHIW, BC_HW_PHI, &
                        BC_C_PHI, M, A_M, B_M, BETA)
! Boundary condition for P1 equation

! Modules
!--------------------------------------------------------------------//
      USE param
      USE param1
      USE geometry
      USE indices
      USE bc
      USE compar
      USE cutcell, only : CARTESIAN_GRID, CG_SAFE_MODE
      USE fun_avg
      USE functions
      IMPLICIT NONE

! Dummy arguments
!--------------------------------------------------------------------//
! The field variable being solved for:
!     e.g., T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(IN) :: VAR(DIMENSION_3)
! Boundary conditions specifications
! bc_phif = flow boundary value
! bc_phiw = wall boundary value
! bc_hw_phi = transfer coefficient
!      = 0 value means specified flux (neumann type)
!      = undefined value means specified wall value (dirichlet type)
!      = other value means mixed type
! bc_C_phi = transfer flux
      DOUBLE PRECISION, INTENT(IN) :: BC_phif(DIMENSION_BC), &
                                      BC_Phiw(DIMENSION_BC), &
                                      BC_hw_Phi(DIMENSION_BC), &
                                      BC_C_Phi(DIMENSION_BC)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Extinction coefficient
      DOUBLE PRECISION, INTENT(IN) :: BETA(DIMENSION_3)
! Local variables
!--------------------------------------------------------------------//
! Boundary condition index
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK, &
                 IM, JM, KM
!--------------------------------------------------------------------//
      DOUBLE PRECISION :: BETAW
! Set up the default walls (i.e., bc_type='dummy' or undefined/default
! boundaries) as non-conducting...
! ---------------------------------------------------------------->>>
      IF(.NOT.CARTESIAN_GRID) THEN
! when setting up default walls do not use cutcells to avoid conflict

      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
!!$omp    parallel do private(IJK, J1, I1)
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)

               IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
                  A_M(KP_OF(IJK),bottom,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value (set
! the boundary cell value equal to adjacent fluid cell value)
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ONE
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
!!$omp    parallel do private(IJK, J1, I1)
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
                  A_M(KM_OF(IJK),top,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! south xz plane
      J1 = 1
!!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JP_OF(IJK),south,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
!!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JM_OF(IJK),north,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! west yz plane
      I1 = imin2
!!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IP_OF(IJK),west,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ONE
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! east yz plane
      I1 = IMAX2
!!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IM_OF(IJK),east,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ONE
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      ENDIF   !(.NOT.CARTESIAN_GRID)

! End setting the default boundary conditions
! ----------------------------------------------------------------<<<


! Set user defined wall boundary conditions
! ---------------------------------------------------------------->>>
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
!-- rad_changes_begin
! the same formula for boundaries of different kind 
!            IF (BC_TYPE_ENUM(L)==NO_SLIP_WALL .OR. &
!                BC_TYPE_ENUM(L)==FREE_SLIP_WALL .OR. &
!                BC_TYPE_ENUM(L)==PAR_SLIP_WALL) THEN
             IF (.TRUE.) THEN
!-- rad_changes_end
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
!!$omp    parallel do private(IJK, K, J, I, IM, JM, KM)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IM = IM1(I)
                        JM = JM1(J)
                        KM = KM1(K)
! first set the boundary cell value equal to the known value in that
! cell
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = VAR(IJK)
! second modify the matrix equation according to the user specified
! boundary condition
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
! specified wall value (i.e., dirichlet type boundary)
                              A_M(IJK,east,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
! if bc_hw__phi=0 then specified flux boundary (i.e., neumann type
! boundary) otherwise a mixed type boundary
          
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(EAST_OF(IJK)))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BETAW+ODX_E(I))
                              A_M(IJK,east,M) = -(HALF*BC_HW_PHI(L)*BETAW-ODX_E(I))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations
!                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(I))
!                              A_M(IJK,east,M) = -(HALF*BC_HW_PHI(L)-ODX_E(I))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                             BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,west,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(WEST_OF(IJK)))
                              A_M(IJK,west,M) = -(HALF*BC_HW_PHI(L)*BETAW-ODX_E(IM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BETAW+ODX_E(IM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations 
!                              A_M(IJK,west,M) = -(HALF*BC_HW_PHI(L)-ODX_E(IM))
!                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(IM))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                            BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,north,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(NORTH_OF(IJK)))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BETAW+ODY_N(J))
                              A_M(IJK,north,M) = -(HALF*BC_HW_PHI(L)*BETAW-ODY_N(J))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations
!                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(J))
!                              A_M(IJK,north,M) = -(HALF*BC_HW_PHI(L)-ODY_N(J))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                             BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,south,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(SOUTH_OF(IJK)))
                              A_M(IJK,south,M) = -(HALF*BC_HW_PHI(L)*BETAW-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BETAW+ODY_N(JM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations
!                              A_M(IJK,south,M) = -(HALF*BC_HW_PHI(L)-ODY_N(JM))
!                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(JM))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                             BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,top,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(TOP_OF(IJK)))
                              A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)*BETAW+OX(I)*ODZ_T(K))
                              A_M(IJK,top,M)=-(HALF*BC_HW_PHI(L)*BETAW-OX(I)*ODZ_T(K))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations
!                              A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(K))
!                              A_M(IJK,top,M)=-(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(K))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                             BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,bottom,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
!-- rad_changes_begin                              
! interpolate beta at the boundary
                              BETAW = HALF*(BETA(IJK)+BETA(BOTTOM_OF(IJK)))
                              A_M(IJK,bottom,M) = -(HALF*BC_HW_PHI(L)*BETAW-&
                                               OX(I)*ODZ_T(KM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)*BETAW+&
                                               OX(I)*ODZ_T(KM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))*BETAW
! comment out old calculations
!                              A_M(IJK,bottom,M) = -(HALF*BC_HW_PHI(L)-&
!                                               OX(I)*ODZ_T(KM))
!                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+&
!                                               OX(I)*ODZ_T(KM))
!                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
!                                             BC_C_PHI(L))
!-- rad_changes_end
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
           ENDIF   ! end if (ns, fs, psw)
         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc
! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

!-- rad_changes_begin
! remove the implementatons for open boundaries
! because the boundary conditions for these boundaries can be 
! unified to the ones above.
!! Set user defined boundary conditions for non-wall cells
!! Setting p_inflow, p_outflow, mass_outflow or outflow flow boundary
!! conditions
!! ---------------------------------------------------------------->>>
!      DO L = 1, DIMENSION_BC
!         IF (BC_DEFINED(L)) THEN
!            IF (BC_TYPE_ENUM(L)==P_OUTFLOW .OR. &
!                BC_TYPE_ENUM(L)==MASS_OUTFLOW .OR. &
!                BC_TYPE_ENUM(L)==OUTFLOW) THEN
!               I1 = BC_I_W(L)
!               I2 = BC_I_E(L)
!               J1 = BC_J_S(L)
!               J2 = BC_J_N(L)
!               K1 = BC_K_B(L)
!               K2 = BC_K_T(L)
!!!$omp    parallel do private(IJK, K, J, I)
!               DO K = K1, K2
!                  DO J = J1, J2
!                     DO I = I1, I2
!                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!                       IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
!                        IJK = FUNIJK(I,J,K)
!! first set the flow boundary cell value equal to zero
!                        A_M(IJK,east,M) = ZERO
!                        A_M(IJK,west,M) = ZERO
!                        A_M(IJK,north,M) = ZERO
!                        A_M(IJK,south,M) = ZERO
!                        A_M(IJK,top,M) = ZERO
!                        A_M(IJK,bottom,M) = ZERO
!                        A_M(IJK,0,M) = -ONE
!                        B_M(IJK,M) = ZERO
!! now set the flow boundary cell value equal to the adjacent fluid
!! cell value
!                        SELECT CASE (TRIM(BC_PLANE(L)))
!                        CASE ('E')
!! fluid cell on the east side
!                           A_M(IJK,east,M) = ONE
!                        CASE ('W')
!! fluid cell on the west side
!                           A_M(IJK,west,M) = ONE
!                        CASE ('N')
!                           A_M(IJK,north,M) = ONE
!                        CASE ('S')
!                           A_M(IJK,south,M) = ONE
!                        CASE ('T')
!                           A_M(IJK,top,M) = ONE
!                        CASE ('B')
!                           A_M(IJK,bottom,M) = ONE
!                        END SELECT
!                     ENDDO
!                  ENDDO
!               ENDDO
!! end setting p_outflow, mass_outflow or outflow flow boundary
!! conditions
!! ----------------------------------------------------------------<<<
!
!            ELSEIF(BC_TYPE_ENUM(L)==P_INFLOW .OR. &
!                   BC_TYPE_ENUM(L)==MASS_INFLOW) THEN
!
!! Setting bc that are defined but not nsw, fsw, psw, p_outflow,
!! mass_outflow or outflow (at this time, this section addresses
!! p_inflow and mass_inflow type boundaries)
!! ----------------------------------------------------------------<<<
!               I1 = BC_I_W(L)
!               I2 = BC_I_E(L)
!               J1 = BC_J_S(L)
!               J2 = BC_J_N(L)
!               K1 = BC_K_B(L)
!               K2 = BC_K_T(L)
!!!$omp    parallel do private(IJK, K, J, I)
!               DO K = K1, K2
!                  DO J = J1, J2
!                     DO I = I1, I2
!                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
!                        IJK = FUNIJK(I,J,K)
!! setting the value in the boundary cell equal to what is known
!                        A_M(IJK,east,M) = ZERO
!                        A_M(IJK,west,M) = ZERO
!                        A_M(IJK,north,M) = ZERO
!                        A_M(IJK,south,M) = ZERO
!                        A_M(IJK,top,M) = ZERO
!                        A_M(IJK,bottom,M) = ZERO
!                        A_M(IJK,0,M) = -ONE
!!                        B_M(IJK,M) = -BC_PHIF(L)  !does not allow the profile to be changed, e.g., from usr1
!                        B_M(IJK,M) = -VAR(IJK)
!                     ENDDO
!                  ENDDO
!               ENDDO
!            ENDIF   ! end if/else (bc_type)
!! end setting of p_inflow or mass_inflow boundary conditions
!! ----------------------------------------------------------------<<<
!
!         ENDIF   ! end if (bc_defined)
!      ENDDO   ! end L do loop over dimension_bc
!-- rad_changes_end



!-- rad_changes_begin
! comment out cartesian grid modifications
! warning: there may be inconsistency
!! modifications for cartesian grid implementation
!      IF(CARTESIAN_GRID .AND. .NOT.(CG_SAFE_MODE(1)==1)) &
!         CALL BC_PHI_CG(VAR, BC_PHIF, BC_PHIW, BC_HW_PHI, &
!                        BC_C_PHI, M, A_M, B_M)
!-- rad_changes_begin

      RETURN
      END SUBROUTINE BC_G

end module 


