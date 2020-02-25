!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE compar
      USE constant
      USE energy
      USE fldvar
      USE fun_avg
      USE functions
      USE funits
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE sendrecv
      USE toleranc
      USE usr

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

!-----------------------------------------------
      INCLUDE 'species.inc'
!      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Index used for looping over all reactions.
      INTEGER L

! Phase index for reacting solids phase :: Coal
      INTEGER, parameter :: mCoal = 1
      INTEGER, parameter :: mRycl = 2

! Proximate Analysis:
      DOUBLE PRECISION :: PAC  ! Char
      DOUBLE PRECISION :: PAV  ! Volatiles
      DOUBLE PRECISION :: PAM  ! Moisture
      DOUBLE PRECISION :: PAA  ! Ash

! Specific gas consant for Oxygen (cm^3.atm/g.K)
      DOUBLE PRECISION, parameter :: R_O2 = 2.564322d0
      DOUBLE PRECISION, parameter :: R_gas = 1.987d0      

! Bounded Phase temperatues (K)
      DOUBLE PRECISION, parameter :: MAX_TEMP = 2.5d3
      DOUBLE PRECISION :: xTg   ! Gas
      DOUBLE PRECISION :: xTs   ! Solids
      DOUBLE PRECISION :: xTs_2 ! Recycled
      DOUBLE PRECISION :: xTgs  ! Average
      DOUBLE PRECISION :: xTgs_2  ! Average      

! Diffusion coefficient for O2 in N2. (cm^2/sec)
      DOUBLE PRECISION :: Diff
! Gas pressure:
      DOUBLE PRECISION :: Pg_atm     ! (atm)
      DOUBLE PRECISION :: PG_atmXMW  ! times mixture moleculre weight

! Gas phase molar concentrations: (g-mole/cm^3)
      DOUBLE PRECISION :: c_O2   ! Oxygen
      DOUBLE PRECISION :: c_CO   ! Carbon monoxide
      DOUBLE PRECISION :: c_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: c_CH4  ! Methane
      DOUBLE PRECISION :: c_H2   ! Hydrogen
      DOUBLE PRECISION :: c_H2O  ! Water Vapor
      DOUBLE PRECISION :: c_Tar  ! Tar
      DOUBLE PRECISION :: c_mix_g  ! gas mixture
      
! Gas phase molar fractions      
      DOUBLE PRECISION :: cx_CO   ! Carbon monoxide
      DOUBLE PRECISION :: cx_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: cx_H2   ! Hydrogen
      DOUBLE PRECISION :: cx_H2O  ! Water Vapor

! Gas phase partial pressures: (atm)
      DOUBLE PRECISION :: p_O2   ! Oxygen
      DOUBLE PRECISION :: p_CO   ! Carbon monoxide
      DOUBLE PRECISION :: p_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: p_CH4  ! Methane
      DOUBLE PRECISION :: p_H2   ! Hydrogen
      DOUBLE PRECISION :: p_H2O  ! Water Vapor
      DOUBLE PRECISION :: p_Tar  ! Tar

! Gas phase mole fractions: (g-mol/g-mol Mix)
      DOUBLE PRECISION :: y_H2   ! Hydrogen

! Ratio of unreacted core to particle diameter (-).
      DOUBLE PRECISION :: rd, rd_2
! Solids phase surface area per unit volume (1/cm)
      DOUBLE PRECISION :: Sa, Sa_2
! Coal phase molar concentrations.  (g-mole/cm^3)
      DOUBLE PRECISION :: c_Moisture   ! Moisture
      DOUBLE PRECISION :: c_Volatiles  ! Volatile Matter
      DOUBLE PRECISION :: mc_Char       ! Char
      DOUBLE PRECISION :: c_Ash        ! Ash  << g/cm^3 >>
! Recycled material
      DOUBLE PRECISION :: c_Char_2
      DOUBLE PRECISION :: c_Ash_2

! Logicals used to 'order' heterogeneous reactions.
! Drying --> Pyrolysis --> Gasification & Combustion
      LOGICAL :: DRIED      ! Initial moisture is gone
      LOGICAL :: PYROLYZED  ! Initial volatiles are gone.

! Mass fraction of Char - Modified for combustion
      DOUBLE PRECISION :: Xs_Char
      DOUBLE PRECISION :: Xs_Char_2
! Volume based char conversion rate
      DOUBLE PRECISION :: charConversion
      DOUBLE PRECISION :: charConversion_2      

! Char combustion rate limiting steps:
      DOUBLE PRECISION :: K_a    ! Ash layer diffusion
      DOUBLE PRECISION :: K_f    ! Film layer diffusion
      DOUBLE PRECISION :: K_r    ! Chemical kinetics
      DOUBLE PRECISION :: K_eff  ! Effective 1/(1/Ka + 1/Kf + 1/Kr)

! Intermediate values for Water gas shift rate calculation.
      DOUBLE PRECISION :: WGS  ! Common rate (forward and backward)
      DOUBLE PRECISION :: FWD  ! Concentration driving forward reaction
      DOUBLE PRECISION :: RVS  ! Concentration driving reverse reaction
      
! Intermediate values for steam gasification 	  
      DOUBLE PRECISION :: eq_SG, rate_SG      
      
! Intermediate values for CO2 gasification      
      DOUBLE PRECISION :: eq_CG, rate_CG
      
! Intermediate values for Methanation
      DOUBLE PRECISION :: eq_HG, rate_HG	

! Volume Fraction of Ash
      DOUBLE PRECISION :: EP_Ash
      DOUBLE PRECISION :: EP_Ash_2      

! These are local aliases for variables that need converted to CGS.
      DOUBLE PRECISION :: Pg                 ! Gas pressure
      DOUBLE PRECISION :: ROg                ! Gas Density
      DOUBLE PRECISION :: ROPs(DIMENSION_M)  ! Solids Bulk Density
      DOUBLE PRECISION :: Dp(DIMENSION_M)    ! Solids Diameter

! Minimum amount of species required to facilitate a reaction.
      DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-6
! Minimum amount of inital species required to stage reactions.
      DOUBLE PRECISION, parameter :: PA_Limiter = 1.0d-4
! Rate limiter to prevent depleting reactants
!      DOUBLE PRECISION:: RATE_LIMIT      

! SI units to CGS:          Conversion     !   SI       |    CGS
      Pg   = P_g(IJK)     ! *  1.0d+1       !   Pa       |   Barye
      ROg  = RO_g(IJK)    ! *  1.0d-3       !   kg/m^3   |   g/cm^3
      ROPs = ROP_s(IJK,:) ! *  1.0d-3       !   kg/m^3   |   g/cm^3
      Dp   = D_p(IJK,:)   ! *  1.0d+2       !   m        |   cm


! Gas phase quantities:
!---------------------------------------------------------------------//
! Initialize gas phase molar concentrations and partial pressures.
      c_O2  = ZERO;   p_O2  = ZERO  ! Oxygen
      c_CO  = ZERO;   p_CO  = ZERO  ! Carbon monoxide
      c_CO2 = ZERO;   p_CO2 = ZERO  ! Carbon dioxide
      c_CH4 = ZERO;   p_CH4 = ZERO  ! Methane
      c_H2  = ZERO;   p_H2  = ZERO  ! Hydrogen
      c_H2O = ZERO;   p_H2O = ZERO  ! Water Vapor
      c_Tar = ZERO;   p_Tar = ZERO  ! Tar

! Initialize gas phase mole fractions.
      y_H2  = ZERO  ! Hydrogen
! Calculate the bounded gas phase temperature (K)
      xTg  = min(MAX_TEMP, T_g(IJK))
! Gas pressure (atm)
      Pg_atm    = Pg / 101.325d4
! Gas pressure multipled by the molecular weight (atm)
      Pg_atmXmw = Pg_atm * MW_MIX_g(IJK)
! Gas mixture molar concentration
      c_mix_g = ROg /  MW_MIX_g(IJK)     
! Calculate the diffusion coefficient for O2 in N2. Field, 1967.
      DIFF = 4.26d0 * ((xTg/1.8d3)**1.75d0) / Pg_atm

! Tar
      IF(X_g(IJK,Tar) .GT. c_Limiter) THEN
         c_Tar = ROg * X_g(IJK,Tar) / MW_g(Tar)       ! (g-mole/cm^3)
!        p_Tar = Pg_atmXmw * X_g(IJK,Tar) / MW_g(Tar) ! Not used
      ENDIF
! Oxygen, O2
      IF(X_g(IJK,O2) .GT. c_Limiter) THEN
         c_O2 = ROg * X_g(IJK,O2) / MW_g(O2)          ! (g-mole/cm^3)
         p_O2 = Pg_atmXmw * X_g(IJK,O2) / MW_g(O2)    ! (atm)
      ENDIF
! Hydrogen, H2
      IF(X_g(IJK,H2) .GT. c_Limiter) THEN
         c_H2 = ROg * X_g(IJK,H2) / MW_g(H2)          ! (g-mole/cm^3)
         p_H2 = Pg_atmXmw * X_g(IJK,H2) / MW_g(H2)    ! (atm)
         y_H2 = X_g(IJK,H2) * MW_MIX_g(IJK)/MW_g(H2)  ! (mol-H2/mol-Mix)
      ENDIF
! Carbon monoxide, CO
      IF(X_g(IJK,CO) .GT. c_Limiter) THEN
         c_CO = ROg * X_g(IJK,CO) / MW_g(CO)           ! (g-mole/cm^3)
         p_CO = Pg_atmXmw * X_g(IJK,CO) / MW_g(CO)     ! (atm)
      ENDIF

! Water Vapor, H2O
      IF(X_g(IJK,H2O) .GT. c_Limiter) THEN
         c_H2O = ROg * X_g(IJK,H2O) / MW_g(H2O)        ! (g-mole/cm^3)
         p_H2O = Pg_atmXmw * X_g(IJK,H2O) / MW_g(H2O)  ! (atm)
      ENDIF

! Carbon dioxide, CO2
      IF(X_g(IJK,CO2) .GT. c_Limiter) THEN
         c_CO2 = ROg * X_g(IJK,CO2) / MW_g(CO2)        ! Not used
         p_CO2 = Pg_atmXmw * X_g(IJK,CO2) / MW_g(CO2)  ! (atm)
      ENDIF
! Methane, CH4
      IF(X_g(IJK,CH4) .GT. c_Limiter) THEN
         c_CH4 = ROg * X_g(IJK,CH4) / MW_g(CH4)        ! (g-mole/cm^3)
!        p_CH4 = Pg_atmXmw * X_g(IJK,CH4) / MW_g(CH4)  ! Not used
      ENDIF

	  cx_CO = c_CO/c_mix_g
	  cx_CO2 = c_CO2/c_mix_g
	  cx_H2 = c_H2/c_mix_g
	  cx_H2O = c_H2O/c_mix_g	  
	  
! Coal phase quantities:
!---------------------------------------------------------------------//
! Calculate the bounded solids temperature (K)
      xTs  = ZERO
      xTs_2= ZERO
! Initialize the gas/solids average temperature (K)
      xTgs = xTg
      xTgs_2 = xTg      
! Initialize ratio of unreacted core diameter to the initial
! particle diameter (cm/cm) --> (-)
      rd = ZERO
      rd_2 = ZERO
! Initialize the surface area per unit volume (1/cm)
      Sa = ZERO
      Sa_2 = ZERO
! Initialize the ash volume fraction
      EP_Ash = ZERO
      EP_Ash_2 = ZERO      
! Initialize the amount of char conversion
      charConversion = ONE

! Set the coal proximate and ultimate analysis:
      PAM = C(3)  ! Moisture
      PAV = C(2) ! Volatiles
      PAC = C(1)      ! Char
      PAA = C(4)       ! Ash

!   c(1)            = 3.0700e-01
!   c(2)            = 3.6910e-01
!   c(3)            = 1.7780e-01
!   c(4)            = 1.4610e-01


!       PAM = proxAnalysis(Moisture)  ! Moisture
!       PAV = proxAnalysis(Volatiles) ! Volatiles
!       PAC = proxAnalysis(Char)      ! Char
!       PAA = proxAnalysis(Ash)       ! Ash


! Initialize coal phase concentrations (g-mol/cm^3).
      c_Moisture  = ZERO   ! Moisture
      c_Volatiles = ZERO   ! Volatiles
      mc_Char      = ZERO   ! Char 
      c_Ash       = ZERO   ! Ash (g/cm^3)

      c_Char_2    = ZERO
      c_Ash_2     = ZERO

! If the fluid cell contains the coal phase, calculate the properties.
      IF(EP_s(IJK,mCoal) .gt. 1d-6) THEN
! Bounded coal phase temperature. (K)
         xTs = min(MAX_TEMP, T_s(IJK,mCoal))
! Surface per unit volume (1/cm)
         Sa = 6.0d0 * EP_s(IJK,mCoal) / Dp(mCoal)
! Ash volume fraction
         EP_Ash = (0.25d0 + 0.75d0*(ONE - PAA))**2.5d0
! Calculate the average gas/solids temperture (K). This value defaults 
! to the gas phase temperture if there are no solids.
         xTgs = HALF * (xTg + xTs)

! Molar concentration Volatile Matter (g-mole/cm^3)
         IF(X_s(IJK,mCoal,Volatiles) .gt. c_Limiter)                    &
            c_Volatiles = ROPs(mCoal) *                                 &
               X_s(IJK,mCoal,Volatiles) / MW_s(mCoal,Volatiles)
! Molar concentration Char (g-mole/cm^3)
         IF(X_s(IJK,mCoal,Char) .gt. c_Limiter)                         &
            mc_Char = ROPs(mCoal) *                                      &
               X_s(IJK,mCoal,Char) / MW_s(mCoal,Char)

! Molar concentration Moisture (g-mole/cm^3)
         IF(X_s(IJK,mCoal,Moisture) .gt. c_Limiter)                     &
            c_Moisture = ROPs(mCoal) *                                  &
               X_s(IJK,mCoal,Moisture) / MW_s(mCoal,Moisture)

         IF(X_s(IJK,mCoal,Ash) .gt. c_Limiter) THEN
! Mass concentration Ash (g/cm^3)
            c_Ash = ROPs(mCoal) * X_s(IJK,mCoal,Ash)
! Char conversion (-)
            charConversion = (X_s(IJK,mCoal,Char) * PAA) / &
               (X_s(IJK,mCoal,Ash) * PAC)
! Ratio of unreacted core diameter to particle diameter (-)
            rd = min(ONE, charConversion**(1.0d0/3.0d0))
         ELSE
            rd = ZERO
         ENDIF

! Logical set true if the majority of the intial moisture was driven
! from the solids. Used to suppress pyrolysis. The secondary logical 
! check is included for initially dry coal.
!!         DRIED = ((PAA * X_s(IJK,mCoal,Moisture)) .LE.                  &
!!            (PAM * X_s(IJK,mCoal,Ash) * 1.0d-3))  .OR.                  &
!!            (PAM .LT. PA_Limiter)

! Logical set true if the majority of initial volatiles were driven
! from the solids. Used to suppress gasification and combustion.
!!         PYROLYZED = DRIED .AND. ((PAA * X_s(IJK,mCoal,Volatiles)) .LE. &
!!            (PAV * X_s(IJK,mCoal,Ash) * 1.0d-3) .OR.                    &
!!             PAV .LT. PA_Limiter)         
      ENDIF         

! If the fluid cell contains the recyled solids phase, calculate the properties.
      IF(EP_s(IJK,mRycl) .gt. 1d-6) THEN
! Bounded coal phase temperature. (K)
         xTs_2 = min(MAX_TEMP, T_s(IJK,mRycl))
! Surface per unit volume (1/cm)
         Sa_2 = 6.0d0 * EP_s(IJK,mRycl) / Dp(mRycl)
! Ash volume fraction
         EP_Ash_2 = (0.25d0 + 0.75d0*(ONE - PAA))**2.5d0
! Calculate the average gas/solids temperture (K). This value defaults 
! to the gas phase temperture if there are no solids.
         xTgs_2 = HALF * (xTg + xTs_2)

! Molar concentration Char (g-mole/cm^3)
         IF(X_s(IJK,mRycl,Char_2) .gt. c_Limiter)                         &
            c_Char_2 = ROPs(mRycl) *                                      &
               X_s(IJK,mRycl,Char_2) / MW_s(mRycl,Char_2)

         IF(X_s(IJK,mRycl,Ash_2) .gt. c_Limiter) THEN
! Mass concentration Ash (g/cm^3)
            c_Ash_2 = ROPs(mRycl) * X_s(IJK,mRycl,Ash_2)
! Char conversion (-) Recyled material is from the coal
            charConversion_2 = (X_s(IJK,mRycl,Char_2) * PAA) / &
               (X_s(IJK,mRycl,Ash_2) * PAC)
! Ratio of unreacted core diameter to particle diameter (-)
            rd_2 = min(ONE, charConversion_2**(1.0d0/3.0d0))
         ELSE
            rd_2 = ZERO
         ENDIF         
     
      ENDIF

!**********************************************************************!
!                                                                      !
!              Heterogeneous and Catalytic Reaction Rates              !
!                                                                      !
!**********************************************************************!

! Water gas shift: CO + H2O <--> CO2 + H2
!---------------------------------------------------------------------//

! Ref: Wen, Chen, and Onozaki (1982)
! Common rate expression for forward and reverse reaction:
! It is catalytic reaction, hence the film temperature is used
! For the current case, the temperature of the recyled solid material is used
! as the system has more.	
! 0.014 is for sub-bitumous coal, it should be higher for lignite ash.
         WGS = 1.4d-1 * 2.877d5 * exp(-1.397d4/xTgs_2)
         WGS = WGS * (Pg_atm**(HALF - Pg_atm/2.5d2))
         WGS = WGS * EP_g(IJK) * (c_Ash + c_Ash_2)             &  !c_ash is mass concentration
            * exp(-8.91d0 + 5.553d3/xTgs_2)
	 WGS = WGS / c_mix_g / c_mix_g
	 
!  Homogeneous reaction from ??????	
!        K = 2780 exp(-1510/Tg)   kmol^-1 m^3 s^-1
!        K* = 0.0265 exp(3958/Tg) -
!        R = k * (c_co*c_h2o - c_co2*c_h2/K*)    c_co: kmol/m3
!	
	 WGS = Max(2.78d0 * exp(-1.510d3/xTg),WGS)
! Concentration driving forward reaction:
         FWD = c_CO * c_H2O
! Concentration driving reverse reaction:
         RVS = c_CO2 * c_H2 /(2.65d-2*exp(3958d0/xTg))
! Assign the net reaction rate.
         IF(FWD > RVS)THEN
            RATES(WaterGasShift_F) = (FWD - RVS)*WGS 
            RATES(WaterGasShift_F) = RATES(WaterGasShift_F)*RATE_LIMIT(c_co)*RATE_LIMIT(c_h2o)
         ELSE
            RATES(WaterGasShift_R) = (RVS - FWD)*WGS 
            RATES(WaterGasShift_R) = RATES(WaterGasShift_R)*RATE_LIMIT(c_co2)*RATE_LIMIT(c_h2)
         ENDIF
	
! Heterogeneous reactions are suppress when the solids phase reactant
! has a volume fraction less than 1.0d-6.
      IF(EP_s(IJK,mCoal) .GT. 1.0d-6) THEN

! Moisture Release: Moisture --> H2O
!---------------------------------------------------------------------//
! MGAS   
         if(xTs .ge. 398.d0)then
         RATES(Drying) = 1.1d5*exp(-2.12d4/R_gas/xTs)*c_Moisture*MW_s(mCoal,Moisture)
         RATES(Drying) = RATES(Drying)*RATE_LIMIT(c_moisture)	
         else
         RATES(Drying) = zero
         endif

! Pyrolysis:  Volatiles --> 
!---------------------------------------------------------------------//
! NOTE: Pyrolysis is suppressed if the majority of moisture has not
! been driven off.    NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
!!	if(DRIED)then
        RATES(Pyrolysis) = &
            57.21d0 * exp(-5560d0/R_gas/xTs)*c_Volatiles
        RATES(Pyrolysis) = RATES(Pyrolysis)*RATE_LIMIT(c_volatiles)
!!        endif

! Steam gasification and char combustion are suppressed if the majority
! of moisture and volatile matter have not been driven off.
!!         IF(PYROLYZED) THEN
         	 
! Gasification kinetics from PCCL based on 100% H2O, CO2 and H2 respectively
! Steam Gasification: C + H2O --> CO + H2
!---------------------------------------------------------------------//
! PCCL     NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
	rate_SG = k_h2o(cx_h2) * exp(e_h2o(cx_h2)/R_gas/xTgs)*(P_H2O**napp_h2o(cx_h2o))/(1.d0+kh2(cx_co)*P_H2)*mc_char
        RATES(SteamGasification) = rate_SG *RATE_LIMIT(mc_char)*RATE_LIMIT(c_h2o)

! CO2 Gasification: Char + CO2 --> 2CO
!---------------------------------------------------------------------//
! PCCL      NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
	rate_CG = k_co2(cx_h2) * exp(e_co2(cx_h2)/R_gas/xTgs)*(P_CO2**napp_co2(cx_co2))/(1.d0+kco(cx_h2)*P_CO)*mc_char	
        RATES(CO2_Gasification) = rate_SG *RATE_LIMIT(mc_char)*RATE_LIMIT(c_co2)	    
            
! Methanation: Char + 2H2 --> CH4
!---------------------------------------------------------------------//
! PCCL            NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014 Not true
        rate_HG = 3.3d-4 * exp(-2.0d4/R_gas/xTgs)*(P_H2**1.d0)*mc_char
	RATES(Methanation) = rate_HG *RATE_LIMIT(mc_char)*RATE_LIMIT(c_h2) 

! Char Combustion: 2C + O2 --> 2CO
!---------------------------------------------------------------------//
! Film Diffusion:  
!	Film temperature is used instead of gas temperature
            K_f = DIFF*N_Sh(IJK,mCoal) / (Dp(mCoal)*R_O2*xTgs)
! Chemical kinetics:
            K_r = 8.71d3 * exp(-1.7965d4/xTs) * rd * rd

! If rd is (near) zero, then there is no char available for combustion.
            IF(rd .gt. SMALL_NUMBER) THEN
               IF(rd .ge. ONE) THEN
! No ash layer resistance.
                  K_eff = ONE/K_f + ONE/K_r
               ELSE
! Ash Diffusion resistance.
                  K_a = 2.0d0 * DIFF * EP_Ash * rd /                   &
                     (Dp(mCoal)*(ONE - rd)*R_O2*xTs)
                  K_eff = ONE/K_f + ONE/K_r + ONE/K_a
               ENDIF
! 'Soft land' the reaction as the mass fraction of char approaches zero.
               Xs_Char = X_s(IJK,mCoal,Char) / &
                  (X_s(IJK,mCoal,Char) + 1.0d-6)
! Calculate the final rate.
               RATES(CharCombustion) = (Sa * p_O2  * Xs_Char) /        &
                  (K_eff * MW_g(O2))
            ENDIF 
!!         ENDIF ! PYROLYZED
      ENDIF ! Eps > 1.0E-6

!  REACTION FOR RECYCLED MATERIAL
!**********************************************************************!
!                                                                      !
!              Heterogeneous and Catalytic Reaction Rates              !
!                                                                      !
!**********************************************************************!
! Heterogeneous reactions are suppress when the solids phase reactant
! has a volume fraction less than 1.0d-6.
      IF(EP_s(IJK,mRycl) .GT. 1.0d-6) THEN

! Gasification kinetics from PCCL based on 100% H2O, CO2 and H2 respectively      	      
! Steam Gasification: C + H2O --> CO + H2
!---------------------------------------------------------------------//
! PCCL   NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
	rate_SG = k_h2o(cx_h2) * exp(e_h2o(cx_h2)/R_gas/xTgs_2)*(P_H2O**napp_h2o(cx_h2o))/(1.d0+kh2(cx_co)*P_H2)*c_char_2
        RATES(SteamGasification_2) = rate_SG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_h2o)	

! CO2 Gasification: Char + CO2 --> 2CO
!---------------------------------------------------------------------//
! PCCL  NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
	rate_CG = k_co2(cx_h2) * exp(e_co2(cx_h2)/R_gas/xTgs_2)*(P_CO2**napp_co2(cx_co2))/(1.d0+kco(cx_h2)*P_CO)*c_char_2	
        RATES(CO2_Gasification_2) = rate_SG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_co2)	          

! Methanation: Char + 2H2 --> CH4
!---------------------------------------------------------------------//        
! PCCL  NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014  Not true        
        rate_HG = 3.3d-4 * exp(-2.0d4/R_gas/xTgs_2)*(P_H2**1.d0)*c_char_2
	RATES(Methanation_2) = rate_HG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_h2)
            
! Char Combustion: 2C + O2 --> 2CO
!---------------------------------------------------------------------//
! Film Diffusion:
            K_f = DIFF*N_Sh(IJK,mRycl) / (Dp(mRycl)*R_O2*xTgs_2)
! Chemical kinetics:
            K_r = 8.71d3 * exp(-1.7965d4/xTs_2) * rd_2 * rd_2

! If rd is (near) zero, then there is no char available for combustion.
            IF(rd_2 .gt. SMALL_NUMBER) THEN
               IF(rd_2 .ge. ONE) THEN
! No ash layer resistance.
                  K_eff = ONE/K_f + ONE/K_r
               ELSE
! Ash Diffusion resistance.
                  K_a = 2.0d0 * DIFF * EP_Ash_2 * rd_2 /                   &
                     (Dp(mRycl)*(ONE - rd_2)*R_O2*xTs_2)
                  K_eff = ONE/K_f + ONE/K_r + ONE/K_a
               ENDIF
! 'Soft land' the reaction as the mass fraction of char approaches zero.
               Xs_Char_2 = X_s(IJK,mRycl,Char_2) / &
                  (X_s(IJK,mRycl,Char_2) + 1.0d-6)
! Calculate the final rate.
               RATES(CharCombustion_2) = (Sa_2 * p_O2  * Xs_Char_2) /        &
                  (K_eff * MW_g(O2))
            ENDIF 

      ENDIF ! Eps > 1.0E-6
	
      

!**********************************************************************!
!                                                                      !
!                 Homogeneous gas phase Reaction Rates                 !
!                                                                      !
!**********************************************************************!


! Tar-Cracking: Tar --> 1.65201*CO  + 0.452635*CO2 + ...
!---------------------------------------------------------------------//
! Ref: PCCL           YYYYYYYYYYYYYYYYY
      RATES(TarCracking) = 350.d0 * exp(-12600d0/R_gas/xTg) * EP_g(IJK) * c_Tar
      RATES(TarCracking) = RATES(TarCracking)*RATE_LIMIT(c_tar)

! Tar-Combustion: Tar + 29.8O2 --> 20.8CO2 + 24.4H2O + 0.1SO3 + 0.8N2
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(TarCombustion) =  3.8d11 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
         (c_O2**1.5d0) * (c_Tar**0.25d0)
      RATES(TarCombustion) = RATES(TarCombustion)*RATE_LIMIT(c_tar)*RATE_LIMIT(c_o2)   

! CO Combustion: 2CO + O2 --> 2CO2
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(CO_Combustion) = 3.98d14 * exp(-2.013d4/xTg) * EP_g(IJK) * &
         (c_O2**0.25d0) * c_CO * (c_H2O**0.5d0)
      RATES(CO_Combustion) = RATES(CO_Combustion)*RATE_LIMIT(c_co)*RATE_LIMIT(c_o2) 

! CH4 Combustion: CH4 + 2O2 --> CO2 + 2H2O
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(CH4_Combustion) = 6.7d12 * exp(-2.436d4/xTg) * EP_g(IJK) * &
         (c_O2**1.3d0) * (c_CH4**0.2d0)
      RATES(CH4_Combustion) = RATES(CH4_Combustion)*RATE_LIMIT(c_ch4)*RATE_LIMIT(c_o2)  
      

! H2 Combustion: 2H2 + O2 --> 2H2O
!---------------------------------------------------------------------//
! Ref: Peters (1979)
      RATES(H2_Combustion) = 1.08d16 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
         c_O2 * c_H2
      RATES(H2_Combustion) = RATES(H2_Combustion)*RATE_LIMIT(c_h2)*RATE_LIMIT(c_o2)           

!**********************************************************************!

!**********************************************************************!
!                         Save the Rates                               !
!**********************************************************************!
      DO L=1, min(nRR,NO_OF_RXNS)
         ReactionRates(IJK,L) = RATES(L)
      ENDDO

      RETURN  
	  CONTAINS	  
!----------------------------------------------------------------------------
!
!   RATE LIMITER	
!
!----------------------------------------------------------------------------	
      double precision function rate_limit( x )
      double precision x
      rate_limit = x/(x + 1.0d-5)
      return
      end function rate_limit 	  
!-------------------------
!
!  Gasification rate parameters
!
!-------------------------	  
	  double precision function k_h2o ( x )   !c_h2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  k_h2o = EXP(-10748*x**4 + 6422.4*x**3 - 1497.7*x**2 + 185.47*x + 4.1828)
	  return
	  end function k_h2o
	  
	  double precision function k_co2 ( x )   !c_h2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  k_co2 = EXP(-10600*x**4 + 6356*x**3 - 1483.6*x**2 + 178.62*x + 10.142)
	  return
	  end function k_co2
	  
	  double precision function e_h2o ( x )   !c_h2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  e_h2o = -(182719*x**5 - 121066*x**4 + 33580*x**3 - 5179.6*x**2 + 516.55*x + 18.52)*1000.d0
	  return
	  end function e_h2o
	  
	  double precision function e_co2 ( x )   !c_h2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  e_co2 = -(122746*x**5 - 91677*x**4 + 28920*x**3 - 5005.3*x**2 + 535.39*x + 35.555)*1000.d0
	  return
	  end function e_co2

	  double precision function napp_h2o ( x )  !c_h2o
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  napp_h2o = -0.5328*x + 0.9927
	  return
	  end function napp_h2o
	  
	  double precision function napp_co2 ( x )  !c_co2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  napp_co2 = -0.6936*x + 0.9972
	  return
	  end function napp_co2
	  
	  double precision function kh2 ( x )   !c_co
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  kh2 = 11.751*x**2 - 5.5042*x + 1.3411
	  return
	  end function kh2
	  
	  double precision function kco ( x )   !c_h2
	  double precision x
	  if(x .gt. 0.2d0)then
	     x = 0.2d0
	  endif
	  kco = -7850.8*x**5 + 4641.3*x**4 - 1061*x**3 + 121.11*x**2 - 7.8113*x + 0.3419
	  return
	  end function kco
	  
      END SUBROUTINE USR_RATES
      
  

	  































! !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! !                                                                      !
! !  Module name: USR_RATES                                              !
! !                                                                      !
! !  Purpose: Hook for user defined reaction rates.                      !
! !                                                                      !
! !  Author: J.Musser                                   Date: 10-Oct-12  !
! !                                                                      !
! !  Comments: Write reaction rates in units of moles/sec.cm^3 (cgs) or  !
! !  modles/sec.m^3 (SI). Units should match those specified in the data !
! !  file.
! !                                                                      !
! !  Example reaction: Methane combustion                                !
! !                                                                      !
! !  mfix.dat input:                                                     !
! !``````````````````````````````````````````````````````````````````````!
! !    @(RXNS)                                                           !
! !      CH4_Comb { chem_eq = "CH4 + 2.0*O2 --> CO2 + 2.0*H2O" }         !
! !    @(END)                                                            !
! !``````````````````````````````````````````````````````````````````````!
! !                                                                      !
! !  usr_rates.f input:                                                  !
! !``````````````````````````````````````````````````````````````````````!
! !    c_O2  = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))                          !
! !    c_CH4 = (RO_g(IJK)*X_g(IJK,CH4)/MW_g(CH4))                        !
! !    RATES(CH4_Comb) = 2.0d5 * EP_g(IJK) * c_O2 * c_CH4                !
! !``````````````````````````````````````````````````````````````````````!
! !  * Species alias and reaction names given in the data file can be    !
! !    used in reference to the reaction index in RATES and a species    !
! !    index in gas/solids phase variables.                              !
! !                                                                      !
! !  * Additional information is provided in section 4.11 of the code    !
! !    Readme.                                                           !
! !                                                                      !
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!       SUBROUTINE USR_RATES(IJK, RATES)
! 
! 
! ! Gas phase global variables:
! !`````````````````````````````````````````````````````````````````````//
! ! Volume fraction.
!       use fldvar, only: EP_g
! ! Species molecular weights
!       use physprop, only: MW_g
! ! Mixture molecular weight
!       use physprop, only: MW_MIX_g
! ! Pressure.
!       use fldvar, only: P_g
! ! Density.
!       use fldvar, only: RO_g
! ! Temperature.
!       use fldvar, only: T_g
! ! Species mass fractions.
!       use fldvar, only: X_g
! 
! 
! ! Solids phase global variables:
! !`````````````````````````````````````````````````````````````````````//
! ! Particle diameter.
!       use fldvar, only: D_p
! ! Species molecular weight.
!       use physprop, only: MW_s
! ! Solid density
!       use fldvar, only: RO_s
! ! Bulk density.
!       use fldvar, only: ROP_s
! ! Material density.
!       use physprop, only: RO_s0
! ! Temperature.
!       use fldvar, only: T_s
! ! Species mass fractions.
!       use fldvar, only: X_s
! ! Number of solids phases
!       use param, only: DIMENSION_M
! 
! 
! ! Reaction global variables:
! !`````````````````````````````````````````````````````````````````````//
! ! Number of reactions.
!       use rxns, only: NO_OF_RXNS
! ! Size of Reaction rates arrays.
!       use rxns, only: nRR
! ! Array for storing reaction rates for output
!       use rxns, only: ReactionRates
! 
! 
! ! Global parameters:
! !`````````````````````````````````````````````````````````````````````//
! ! Gas constant (cal/mol.K)
!       use constant, only: RGAS => GAS_CONST_cal
! ! Double percision aliases:
!       use param1, only: ZERO
!       use param1, only: SMALL_NUMBER
!       use param1, only: HALF
!       use param1, only: ONE
! 
! ! Global user variables:
! !`````````````````````````````````````````````````````````````````````//
! ! Sherwood number :: calculated in usr1.f
!       use usr, only: N_sh
! ! Dummy variable required by usrnlst.inc (UNUSED)
!       use usr, only: DUMMY_DP
! ! Proximate Analysis:
!       use constant, only: proxAnalysis => C
! 
! 
!       implicit none
! 
! 
! ! Passed arguments:
! !`````````````````````````````````````````````````````````````````````//
! ! Fluid cell index
!       INTEGER, INTENT(IN) :: IJK
! ! Array for storing reaction rates.
!       DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES
! 
! 
! ! Included parameters and functions:
! !`````````````````````````````````````````````````````````````````````//
!       include 'species.inc'
!       include 'ep_s1.inc'
!       include 'ep_s2.inc'
!       include 'usrnlst.inc'
! 
! 
! ! Reaction specific variables:
! !`````````````````````````````````````````````````````````````````````//
! ! Index used for looping over all reactions.
!       INTEGER L
! 
! ! Phase index for reacting solids phase :: Coal
!       INTEGER, parameter :: mCoal = 1
!       INTEGER, parameter :: mRycl = 2
! 
! ! Proximate Analysis:
!       DOUBLE PRECISION :: PAC  ! Char
!       DOUBLE PRECISION :: PAV  ! Volatiles
!       DOUBLE PRECISION :: PAM  ! Moisture
!       DOUBLE PRECISION :: PAA  ! Ash
! 
! ! Specific gas consant for Oxygen (cm^3.atm/g.K)
!       DOUBLE PRECISION, parameter :: R_O2 = 2.564322d0
!       DOUBLE PRECISION, parameter :: R_gas = 1.987d0      
! 
! ! Bounded Phase temperatues (K)
!       DOUBLE PRECISION, parameter :: MAX_TEMP = 2.5d3
!       DOUBLE PRECISION :: xTg   ! Gas
!       DOUBLE PRECISION :: xTs   ! Solids
!       DOUBLE PRECISION :: xTs_2 ! Recycled
!       DOUBLE PRECISION :: xTgs  ! Average
!       DOUBLE PRECISION :: xTgs_2  ! Average      
! 
! ! Diffusion coefficient for O2 in N2. (cm^2/sec)
!       DOUBLE PRECISION :: Diff
! ! Gas pressure:
!       DOUBLE PRECISION :: Pg_atm     ! (atm)
!       DOUBLE PRECISION :: PG_atmXMW  ! times mixture moleculre weight
! 
! ! Gas phase molar concentrations: (g-mole/cm^3)
!       DOUBLE PRECISION :: c_O2   ! Oxygen
!       DOUBLE PRECISION :: c_CO   ! Carbon monoxide
!       DOUBLE PRECISION :: c_CO2  ! Carbon dioxide
!       DOUBLE PRECISION :: c_CH4  ! Methane
!       DOUBLE PRECISION :: c_H2   ! Hydrogen
!       DOUBLE PRECISION :: c_H2O  ! Water Vapor
!       DOUBLE PRECISION :: c_Tar  ! Tar
!       DOUBLE PRECISION :: c_mix_g  ! gas mixture
!       
! ! Gas phase molar fractions      
!       DOUBLE PRECISION :: cx_CO   ! Carbon monoxide
!       DOUBLE PRECISION :: cx_CO2  ! Carbon dioxide
!       DOUBLE PRECISION :: cx_H2   ! Hydrogen
!       DOUBLE PRECISION :: cx_H2O  ! Water Vapor
! 
! ! Gas phase partial pressures: (atm)
!       DOUBLE PRECISION :: p_O2   ! Oxygen
!       DOUBLE PRECISION :: p_CO   ! Carbon monoxide
!       DOUBLE PRECISION :: p_CO2  ! Carbon dioxide
!       DOUBLE PRECISION :: p_CH4  ! Methane
!       DOUBLE PRECISION :: p_H2   ! Hydrogen
!       DOUBLE PRECISION :: p_H2O  ! Water Vapor
!       DOUBLE PRECISION :: p_Tar  ! Tar
! 
! ! Gas phase mole fractions: (g-mol/g-mol Mix)
!       DOUBLE PRECISION :: y_H2   ! Hydrogen
! 
! ! Ratio of unreacted core to particle diameter (-).
!       DOUBLE PRECISION :: rd, rd_2
! ! Solids phase surface area per unit volume (1/cm)
!       DOUBLE PRECISION :: Sa, Sa_2
! ! Coal phase molar concentrations.  (g-mole/cm^3)
!       DOUBLE PRECISION :: c_Moisture   ! Moisture
!       DOUBLE PRECISION :: c_Volatiles  ! Volatile Matter
!       DOUBLE PRECISION :: c_Char       ! Char
!       DOUBLE PRECISION :: c_Ash        ! Ash  << g/cm^3 >>
! ! Recycled material
!       DOUBLE PRECISION :: c_Char_2
!       DOUBLE PRECISION :: c_Ash_2
! 
! ! Logicals used to 'order' heterogeneous reactions.
! ! Drying --> Pyrolysis --> Gasification & Combustion
!       LOGICAL :: DRIED      ! Initial moisture is gone
!       LOGICAL :: PYROLYZED  ! Initial volatiles are gone.
! 
! ! Mass fraction of Char - Modified for combustion
!       DOUBLE PRECISION :: Xs_Char
!       DOUBLE PRECISION :: Xs_Char_2
! ! Volume based char conversion rate
!       DOUBLE PRECISION :: charConversion
!       DOUBLE PRECISION :: charConversion_2      
! 
! ! Char combustion rate limiting steps:
!       DOUBLE PRECISION :: K_a    ! Ash layer diffusion
!       DOUBLE PRECISION :: K_f    ! Film layer diffusion
!       DOUBLE PRECISION :: K_r    ! Chemical kinetics
!       DOUBLE PRECISION :: K_eff  ! Effective 1/(1/Ka + 1/Kf + 1/Kr)
! 
! ! Intermediate values for Water gas shift rate calculation.
!       DOUBLE PRECISION :: WGS  ! Common rate (forward and backward)
!       DOUBLE PRECISION :: FWD  ! Concentration driving forward reaction
!       DOUBLE PRECISION :: RVS  ! Concentration driving reverse reaction
!       
! ! Intermediate values for steam gasification 	  
!       DOUBLE PRECISION :: eq_SG, rate_SG      
!       
! ! Intermediate values for CO2 gasification      
!       DOUBLE PRECISION :: eq_CG, rate_CG
!       
! ! Intermediate values for Methanation
!       DOUBLE PRECISION :: eq_HG, rate_HG	
! 
! ! Volume Fraction of Ash
!       DOUBLE PRECISION :: EP_Ash
!       DOUBLE PRECISION :: EP_Ash_2      
! 
! ! These are local aliases for variables that need converted to CGS.
!       DOUBLE PRECISION :: Pg                 ! Gas pressure
!       DOUBLE PRECISION :: ROg                ! Gas Density
!       DOUBLE PRECISION :: ROPs(DIMENSION_M)  ! Solids Bulk Density
!       DOUBLE PRECISION :: Dp(DIMENSION_M)    ! Solids Diameter
! 
! ! Minimum amount of species required to facilitate a reaction.
!       DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-6
! ! Minimum amount of inital species required to stage reactions.
!       DOUBLE PRECISION, parameter :: PA_Limiter = 1.0d-4
! ! Rate limiter to prevent depleting reactants
! !      DOUBLE PRECISION:: RATE_LIMIT      
! 
! ! SI units to CGS:          Conversion     !   SI       |    CGS
!       Pg   = P_g(IJK)     ! *  1.0d+1       !   Pa       |   Barye
!       ROg  = RO_g(IJK)    ! *  1.0d-3       !   kg/m^3   |   g/cm^3
!       ROPs = ROP_s(IJK,:) ! *  1.0d-3       !   kg/m^3   |   g/cm^3
!       Dp   = D_p(IJK,:)   ! *  1.0d+2       !   m        |   cm
! 
! 
! ! Gas phase quantities:
! !---------------------------------------------------------------------//
! ! Initialize gas phase molar concentrations and partial pressures.
!       c_O2  = ZERO;   p_O2  = ZERO  ! Oxygen
!       c_CO  = ZERO;   p_CO  = ZERO  ! Carbon monoxide
!       c_CO2 = ZERO;   p_CO2 = ZERO  ! Carbon dioxide
!       c_CH4 = ZERO;   p_CH4 = ZERO  ! Methane
!       c_H2  = ZERO;   p_H2  = ZERO  ! Hydrogen
!       c_H2O = ZERO;   p_H2O = ZERO  ! Water Vapor
!       c_Tar = ZERO;   p_Tar = ZERO  ! Tar
! 
! ! Initialize gas phase mole fractions.
!       y_H2  = ZERO  ! Hydrogen
! ! Calculate the bounded gas phase temperature (K)
!       xTg  = min(MAX_TEMP, T_g(IJK))
! ! Gas pressure (atm)
!       Pg_atm    = Pg / 101.325d4
! ! Gas pressure multipled by the molecular weight (atm)
!       Pg_atmXmw = Pg_atm * MW_MIX_g(IJK)
! ! Gas mixture molar concentration
!       c_mix_g = ROg /  MW_MIX_g(IJK)     
! ! Calculate the diffusion coefficient for O2 in N2. Field, 1967.
!       DIFF = 4.26d0 * ((xTg/1.8d3)**1.75d0) / Pg_atm
! 
! ! Tar
!       IF(X_g(IJK,Tar) .GT. c_Limiter) THEN
!          c_Tar = ROg * X_g(IJK,Tar) / MW_g(Tar)       ! (g-mole/cm^3)
! !        p_Tar = Pg_atmXmw * X_g(IJK,Tar) / MW_g(Tar) ! Not used
!       ENDIF
! ! Oxygen, O2
!       IF(X_g(IJK,O2) .GT. c_Limiter) THEN
!          c_O2 = ROg * X_g(IJK,O2) / MW_g(O2)          ! (g-mole/cm^3)
!          p_O2 = Pg_atmXmw * X_g(IJK,O2) / MW_g(O2)    ! (atm)
!       ENDIF
! ! Hydrogen, H2
!       IF(X_g(IJK,H2) .GT. c_Limiter) THEN
!          c_H2 = ROg * X_g(IJK,H2) / MW_g(H2)          ! (g-mole/cm^3)
!          p_H2 = Pg_atmXmw * X_g(IJK,H2) / MW_g(H2)    ! (atm)
!          y_H2 = X_g(IJK,H2) * MW_MIX_g(IJK)/MW_g(H2)  ! (mol-H2/mol-Mix)
!       ENDIF
! ! Carbon monoxide, CO
!       IF(X_g(IJK,CO) .GT. c_Limiter) THEN
!          c_CO = ROg * X_g(IJK,CO) / MW_g(CO)           ! (g-mole/cm^3)
!          p_CO = Pg_atmXmw * X_g(IJK,CO) / MW_g(CO)     ! (atm)
!       ENDIF
! 
! ! Water Vapor, H2O
!       IF(X_g(IJK,H2O) .GT. c_Limiter) THEN
!          c_H2O = ROg * X_g(IJK,H2O) / MW_g(H2O)        ! (g-mole/cm^3)
!          p_H2O = Pg_atmXmw * X_g(IJK,H2O) / MW_g(H2O)  ! (atm)
!       ENDIF
! 
! ! Carbon dioxide, CO2
!       IF(X_g(IJK,CO2) .GT. c_Limiter) THEN
!          c_CO2 = ROg * X_g(IJK,CO2) / MW_g(CO2)        ! Not used
!          p_CO2 = Pg_atmXmw * X_g(IJK,CO2) / MW_g(CO2)  ! (atm)
!       ENDIF
! ! Methane, CH4
!       IF(X_g(IJK,CH4) .GT. c_Limiter) THEN
!          c_CH4 = ROg * X_g(IJK,CH4) / MW_g(CH4)        ! (g-mole/cm^3)
! !        p_CH4 = Pg_atmXmw * X_g(IJK,CH4) / MW_g(CH4)  ! Not used
!       ENDIF
! 
! 	  cx_CO = c_CO/c_mix_g
! 	  cx_CO2 = c_CO2/c_mix_g
! 	  cx_H2 = c_H2/c_mix_g
! 	  cx_H2O = c_H2O/c_mix_g	  
! 	  
! ! Coal phase quantities:
! !---------------------------------------------------------------------//
! ! Calculate the bounded solids temperature (K)
!       xTs  = ZERO
!       xTs_2= ZERO
! ! Initialize the gas/solids average temperature (K)
!       xTgs = xTg
!       xTgs_2 = xTg      
! ! Initialize ratio of unreacted core diameter to the initial
! ! particle diameter (cm/cm) --> (-)
!       rd = ZERO
!       rd_2 = ZERO
! ! Initialize the surface area per unit volume (1/cm)
!       Sa = ZERO
!       Sa_2 = ZERO
! ! Initialize the ash volume fraction
!       EP_Ash = ZERO
!       EP_Ash_2 = ZERO      
! ! Initialize the amount of char conversion
!       charConversion = ONE
! 
! ! Set the coal proximate and ultimate analysis:
!       PAM = proxAnalysis(Moisture)  ! Moisture
!       PAV = proxAnalysis(Volatiles) ! Volatiles
!       PAC = proxAnalysis(Char)      ! Char
!       PAA = proxAnalysis(Ash)       ! Ash
! 
! 
! ! Initialize coal phase concentrations (g-mol/cm^3).
!       c_Moisture  = ZERO   ! Moisture
!       c_Volatiles = ZERO   ! Volatiles
!       c_Char      = ZERO   ! Char 
!       c_Ash       = ZERO   ! Ash (g/cm^3)
! 
!       c_Char_2    = ZERO
!       c_Ash_2     = ZERO
! 
! ! If the fluid cell contains the coal phase, calculate the properties.
!       IF(EP_s(IJK,mCoal) .gt. 1d-6) THEN
! ! Bounded coal phase temperature. (K)
!          xTs = min(MAX_TEMP, T_s(IJK,mCoal))
! ! Surface per unit volume (1/cm)
!          Sa = 6.0d0 * EP_s(IJK,mCoal) / Dp(mCoal)
! ! Ash volume fraction
!          EP_Ash = (0.25d0 + 0.75d0*(ONE - PAA))**2.5d0
! ! Calculate the average gas/solids temperture (K). This value defaults 
! ! to the gas phase temperture if there are no solids.
!          xTgs = HALF * (xTg + xTs)
! 
! ! Molar concentration Volatile Matter (g-mole/cm^3)
!          IF(X_s(IJK,mCoal,Volatiles) .gt. c_Limiter)                    &
!             c_Volatiles = ROPs(mCoal) *                                 &
!                X_s(IJK,mCoal,Volatiles) / MW_s(mCoal,Volatiles)
! ! Molar concentration Char (g-mole/cm^3)
!          IF(X_s(IJK,mCoal,Char) .gt. c_Limiter)                         &
!             c_Char = ROPs(mCoal) *                                      &
!                X_s(IJK,mCoal,Char) / MW_s(mCoal,Char)
! 
! ! Molar concentration Moisture (g-mole/cm^3)
!          IF(X_s(IJK,mCoal,Moisture) .gt. c_Limiter)                     &
!             c_Moisture = ROPs(mCoal) *                                  &
!                X_s(IJK,mCoal,Moisture) / MW_s(mCoal,Moisture)
! 
!          IF(X_s(IJK,mCoal,Ash) .gt. c_Limiter) THEN
! ! Mass concentration Ash (g/cm^3)
!             c_Ash = ROPs(mCoal) * X_s(IJK,mCoal,Ash)
! ! Char conversion (-)
!             charConversion = (X_s(IJK,mCoal,Char) * PAA) / &
!                (X_s(IJK,mCoal,Ash) * PAC)
! ! Ratio of unreacted core diameter to particle diameter (-)
!             rd = min(ONE, charConversion**(1.0d0/3.0d0))
!          ELSE
!             rd = ZERO
!          ENDIF
! 
! ! Logical set true if the majority of the intial moisture was driven
! ! from the solids. Used to suppress pyrolysis. The secondary logical 
! ! check is included for initially dry coal.
! !!         DRIED = ((PAA * X_s(IJK,mCoal,Moisture)) .LE.                  &
! !!            (PAM * X_s(IJK,mCoal,Ash) * 1.0d-3))  .OR.                  &
! !!            (PAM .LT. PA_Limiter)
! 
! ! Logical set true if the majority of initial volatiles were driven
! ! from the solids. Used to suppress gasification and combustion.
! !!         PYROLYZED = DRIED .AND. ((PAA * X_s(IJK,mCoal,Volatiles)) .LE. &
! !!            (PAV * X_s(IJK,mCoal,Ash) * 1.0d-3) .OR.                    &
! !!             PAV .LT. PA_Limiter)         
!       ENDIF         
! 
! ! If the fluid cell contains the recyled solids phase, calculate the properties.
!       IF(EP_s(IJK,mRycl) .gt. 1d-6) THEN
! ! Bounded coal phase temperature. (K)
!          xTs_2 = min(MAX_TEMP, T_s(IJK,mRycl))
! ! Surface per unit volume (1/cm)
!          Sa_2 = 6.0d0 * EP_s(IJK,mRycl) / Dp(mRycl)
! ! Ash volume fraction
!          EP_Ash_2 = (0.25d0 + 0.75d0*(ONE - PAA))**2.5d0
! ! Calculate the average gas/solids temperture (K). This value defaults 
! ! to the gas phase temperture if there are no solids.
!          xTgs_2 = HALF * (xTg + xTs_2)
! 
! ! Molar concentration Char (g-mole/cm^3)
!          IF(X_s(IJK,mRycl,Char_2) .gt. c_Limiter)                         &
!             c_Char_2 = ROPs(mRycl) *                                      &
!                X_s(IJK,mRycl,Char_2) / MW_s(mRycl,Char_2)
! 
!          IF(X_s(IJK,mRycl,Ash_2) .gt. c_Limiter) THEN
! ! Mass concentration Ash (g/cm^3)
!             c_Ash_2 = ROPs(mRycl) * X_s(IJK,mRycl,Ash_2)
! ! Char conversion (-) Recyled material is from the coal
!             charConversion_2 = (X_s(IJK,mRycl,Char_2) * PAA) / &
!                (X_s(IJK,mRycl,Ash_2) * PAC)
! ! Ratio of unreacted core diameter to particle diameter (-)
!             rd_2 = min(ONE, charConversion_2**(1.0d0/3.0d0))
!          ELSE
!             rd_2 = ZERO
!          ENDIF         
!      
!       ENDIF
! 
! !**********************************************************************!
! !                                                                      !
! !              Heterogeneous and Catalytic Reaction Rates              !
! !                                                                      !
! !**********************************************************************!
! 
! ! Water gas shift: CO + H2O <--> CO2 + H2
! !---------------------------------------------------------------------//
! 
! ! Ref: Wen, Chen, and Onozaki (1982)
! ! Common rate expression for forward and reverse reaction:
! ! It is catalytic reaction, hence the film temperature is used
! ! For the current case, the temperature of the recyled solid material is used
! ! as the system has more.	
! ! 0.014 is for sub-bitumous coal, it should be higher for lignite ash.
!          WGS = 1.4d-1 * 2.877d5 * exp(-1.397d4/xTgs_2)
!          WGS = WGS * (Pg_atm**(HALF - Pg_atm/2.5d2))
!          WGS = WGS * EP_g(IJK) * (c_Ash + c_Ash_2)             &  !c_ash is mass concentration
!             * exp(-8.91d0 + 5.553d3/xTgs_2)
! 	 WGS = WGS / c_mix_g / c_mix_g
! 	 
! !  Homogeneous reaction from ??????	
! !        K = 2780 exp(-1510/Tg)   kmol^-1 m^3 s^-1
! !        K* = 0.0265 exp(3958/Tg) -
! !        R = k * (c_co*c_h2o - c_co2*c_h2/K*)    c_co: kmol/m3
! !	
! 	 WGS = Max(2.78d0 * exp(-1.510d3/xTg),WGS)
! ! Concentration driving forward reaction:
!          FWD = c_CO * c_H2O
! ! Concentration driving reverse reaction:
!          RVS = c_CO2 * c_H2 /(2.65d-2*exp(3958d0/xTg))
! ! Assign the net reaction rate.
!          IF(FWD > RVS)THEN
!             RATES(WaterGasShift_F) = (FWD - RVS)*WGS 
!             RATES(WaterGasShift_F) = RATES(WaterGasShift_F)*RATE_LIMIT(c_co)*RATE_LIMIT(c_h2o)
!          ELSE
!             RATES(WaterGasShift_R) = (RVS - FWD)*WGS 
!             RATES(WaterGasShift_R) = RATES(WaterGasShift_R)*RATE_LIMIT(c_co2)*RATE_LIMIT(c_h2)
!          ENDIF
! 	
! ! Heterogeneous reactions are suppress when the solids phase reactant
! ! has a volume fraction less than 1.0d-6.
!       IF(EP_s(IJK,mCoal) .GT. 1.0d-6) THEN
! 
! ! Moisture Release: Moisture --> H2O
! !---------------------------------------------------------------------//
! ! MGAS   
!          if(xTs .ge. 398.d0)then
!          RATES(Drying) = 1.1d5*exp(-2.12d4/R_gas/xTs)*c_Moisture*MW_s(mCoal,Moisture)
!          RATES(Drying) = RATES(Drying)*RATE_LIMIT(c_moisture)	
!          else
!          RATES(Drying) = zero
!          endif
! 
! ! Pyrolysis:  Volatiles --> 
! !---------------------------------------------------------------------//
! ! NOTE: Pyrolysis is suppressed if the majority of moisture has not
! ! been driven off.    NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
! !!	if(DRIED)then
!         RATES(Pyrolysis) = &
!             57.21d0 * exp(-5560d0/R_gas/xTs)*c_Volatiles
!         RATES(Pyrolysis) = RATES(Pyrolysis)*RATE_LIMIT(c_volatiles)
! !!        endif
! 
! ! Steam gasification and char combustion are suppressed if the majority
! ! of moisture and volatile matter have not been driven off.
! !!         IF(PYROLYZED) THEN
!          	 
! ! Gasification kinetics from PCCL based on 100% H2O, CO2 and H2 respectively
! ! Steam Gasification: C + H2O --> CO + H2
! !---------------------------------------------------------------------//
! ! PCCL     NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
! 	rate_SG = k_h2o(cx_h2) * exp(e_h2o(cx_h2)/R_gas/xTgs)*(P_H2O**napp_h2o(cx_h2o))/(1.d0+kh2(cx_co)*P_H2)*c_char
!         RATES(SteamGasification) = rate_SG *RATE_LIMIT(c_char)*RATE_LIMIT(c_h2o)
! 
! ! CO2 Gasification: Char + CO2 --> 2CO
! !---------------------------------------------------------------------//
! ! PCCL      NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
! 	rate_CG = k_co2(cx_h2) * exp(e_co2(cx_h2)/R_gas/xTgs)*(P_CO2**napp_co2(cx_co2))/(1.d0+kco(cx_h2)*P_CO)*c_char	
!         RATES(CO2_Gasification) = rate_SG *RATE_LIMIT(c_char)*RATE_LIMIT(c_co2)	    
!             
! ! Methanation: Char + 2H2 --> CH4
! !---------------------------------------------------------------------//
! ! PCCL            NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014 Not true
!         rate_HG = 3.3d-4 * exp(-2.0d4/R_gas/xTgs)*(P_H2**1.d0)*c_char
! 	RATES(Methanation) = rate_HG *RATE_LIMIT(c_char)*RATE_LIMIT(c_h2) 
! 
! ! Char Combustion: 2C + O2 --> 2CO
! !---------------------------------------------------------------------//
! ! Film Diffusion:  
! !	Film temperature is used instead of gas temperature
!             K_f = DIFF*N_Sh(IJK,mCoal) / (Dp(mCoal)*R_O2*xTgs)
! ! Chemical kinetics:
!             K_r = 8.71d3 * exp(-1.7965d4/xTs) * rd * rd
! 
! ! If rd is (near) zero, then there is no char available for combustion.
!             IF(rd .gt. SMALL_NUMBER) THEN
!                IF(rd .ge. ONE) THEN
! ! No ash layer resistance.
!                   K_eff = ONE/K_f + ONE/K_r
!                ELSE
! ! Ash Diffusion resistance.
!                   K_a = 2.0d0 * DIFF * EP_Ash * rd /                   &
!                      (Dp(mCoal)*(ONE - rd)*R_O2*xTs)
!                   K_eff = ONE/K_f + ONE/K_r + ONE/K_a
!                ENDIF
! ! 'Soft land' the reaction as the mass fraction of char approaches zero.
!                Xs_Char = X_s(IJK,mCoal,Char) / &
!                   (X_s(IJK,mCoal,Char) + 1.0d-6)
! ! Calculate the final rate.
!                RATES(CharCombustion) = (Sa * p_O2  * Xs_Char) /        &
!                   (K_eff * MW_g(O2))
!             ENDIF 
! !!         ENDIF ! PYROLYZED
!       ENDIF ! Eps > 1.0E-6
! 
! !  REACTION FOR RECYCLED MATERIAL
! !**********************************************************************!
! !                                                                      !
! !              Heterogeneous and Catalytic Reaction Rates              !
! !                                                                      !
! !**********************************************************************!
! ! Heterogeneous reactions are suppress when the solids phase reactant
! ! has a volume fraction less than 1.0d-6.
!       IF(EP_s(IJK,mRycl) .GT. 1.0d-6) THEN
! 
! ! Gasification kinetics from PCCL based on 100% H2O, CO2 and H2 respectively      	      
! ! Steam Gasification: C + H2O --> CO + H2
! !---------------------------------------------------------------------//
! ! PCCL   NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
! 	rate_SG = k_h2o(cx_h2) * exp(e_h2o(cx_h2)/R_gas/xTgs_2)*(P_H2O**napp_h2o(cx_h2o))/(1.d0+kh2(cx_co)*P_H2)*c_char_2
!         RATES(SteamGasification_2) = rate_SG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_h2o)	
! 
! ! CO2 Gasification: Char + CO2 --> 2CO
! !---------------------------------------------------------------------//
! ! PCCL  NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014
! 	rate_CG = k_co2(cx_h2) * exp(e_co2(cx_h2)/R_gas/xTgs_2)*(P_CO2**napp_co2(cx_co2))/(1.d0+kco(cx_h2)*P_CO)*c_char_2	
!         RATES(CO2_Gasification_2) = rate_SG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_co2)	          
! 
! ! Methanation: Char + 2H2 --> CH4
! !---------------------------------------------------------------------//        
! ! PCCL  NNNNNNNNNNNNNNNNNNNNNNNNFeb 5, 2014  Not true        
!         rate_HG = 3.3d-4 * exp(-2.0d4/R_gas/xTgs_2)*(P_H2**1.d0)*c_char_2
! 	RATES(Methanation_2) = rate_HG *RATE_LIMIT(c_char_2)*RATE_LIMIT(c_h2)
!             
! ! Char Combustion: 2C + O2 --> 2CO
! !---------------------------------------------------------------------//
! ! Film Diffusion:
!             K_f = DIFF*N_Sh(IJK,mRycl) / (Dp(mRycl)*R_O2*xTgs_2)
! ! Chemical kinetics:
!             K_r = 8.71d3 * exp(-1.7965d4/xTs_2) * rd_2 * rd_2
! 
! ! If rd is (near) zero, then there is no char available for combustion.
!             IF(rd_2 .gt. SMALL_NUMBER) THEN
!                IF(rd_2 .ge. ONE) THEN
! ! No ash layer resistance.
!                   K_eff = ONE/K_f + ONE/K_r
!                ELSE
! ! Ash Diffusion resistance.
!                   K_a = 2.0d0 * DIFF * EP_Ash_2 * rd_2 /                   &
!                      (Dp(mRycl)*(ONE - rd_2)*R_O2*xTs_2)
!                   K_eff = ONE/K_f + ONE/K_r + ONE/K_a
!                ENDIF
! ! 'Soft land' the reaction as the mass fraction of char approaches zero.
!                Xs_Char_2 = X_s(IJK,mRycl,Char_2) / &
!                   (X_s(IJK,mRycl,Char_2) + 1.0d-6)
! ! Calculate the final rate.
!                RATES(CharCombustion_2) = (Sa_2 * p_O2  * Xs_Char_2) /        &
!                   (K_eff * MW_g(O2))
!             ENDIF 
! 
!       ENDIF ! Eps > 1.0E-6
! 	
!       
! 
! !**********************************************************************!
! !                                                                      !
! !                 Homogeneous gas phase Reaction Rates                 !
! !                                                                      !
! !**********************************************************************!
! 
! 
! ! Tar-Cracking: Tar --> 1.65201*CO  + 0.452635*CO2 + ...
! !---------------------------------------------------------------------//
! ! Ref: PCCL           YYYYYYYYYYYYYYYYY
!       RATES(TarCracking) = 350.d0 * exp(-12600d0/R_gas/xTg) * EP_g(IJK) * c_Tar
!       RATES(TarCracking) = RATES(TarCracking)*RATE_LIMIT(c_tar)
! 
! ! Tar-Combustion: Tar + 29.8O2 --> 20.8CO2 + 24.4H2O + 0.1SO3 + 0.8N2
! !---------------------------------------------------------------------//
! ! Ref: Westbrook and Dryer (1981)
!       RATES(TarCombustion) =  3.8d11 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
!          (c_O2**1.5d0) * (c_Tar**0.25d0)
!       RATES(TarCombustion) = RATES(TarCombustion)*RATE_LIMIT(c_tar)*RATE_LIMIT(c_o2)   
! 
! ! CO Combustion: 2CO + O2 --> 2CO2
! !---------------------------------------------------------------------//
! ! Ref: Westbrook and Dryer (1981)
!       RATES(CO_Combustion) = 3.98d14 * exp(-2.013d4/xTg) * EP_g(IJK) * &
!          (c_O2**0.25d0) * c_CO * (c_H2O**0.5d0)
!       RATES(CO_Combustion) = RATES(CO_Combustion)*RATE_LIMIT(c_co)*RATE_LIMIT(c_o2) 
! 
! ! CH4 Combustion: CH4 + 2O2 --> CO2 + 2H2O
! !---------------------------------------------------------------------//
! ! Ref: Westbrook and Dryer (1981)
!       RATES(CH4_Combustion) = 6.7d12 * exp(-2.436d4/xTg) * EP_g(IJK) * &
!          (c_O2**1.3d0) * (c_CH4**0.2d0)
!       RATES(CH4_Combustion) = RATES(CH4_Combustion)*RATE_LIMIT(c_ch4)*RATE_LIMIT(c_o2)  
!       
! 
! ! H2 Combustion: 2H2 + O2 --> 2H2O
! !---------------------------------------------------------------------//
! ! Ref: Peters (1979)
!       RATES(H2_Combustion) = 1.08d16 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
!          c_O2 * c_H2
!       RATES(H2_Combustion) = RATES(H2_Combustion)*RATE_LIMIT(c_h2)*RATE_LIMIT(c_o2)           
! 
! !**********************************************************************!
! 
! !**********************************************************************!
! !                         Save the Rates                               !
! !**********************************************************************!
!       DO L=1, min(nRR,NO_OF_RXNS)
!          ReactionRates(IJK,L) = RATES(L)
!       ENDDO
! 
!       RETURN  
! 	  CONTAINS	  
! !----------------------------------------------------------------------------
! !
! !   RATE LIMITER	
! !
! !----------------------------------------------------------------------------	
!       double precision function rate_limit( x )
!       double precision x
!       rate_limit = x/(x + 1.0d-5)
!       return
!       end function rate_limit 	  
! !-------------------------
! !
! !  Gasification rate parameters
! !
! !-------------------------	  
! 	  double precision function k_h2o ( x )   !c_h2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  k_h2o = EXP(-10748*x**4 + 6422.4*x**3 - 1497.7*x**2 + 185.47*x + 4.1828)
! 	  return
! 	  end function k_h2o
! 	  
! 	  double precision function k_co2 ( x )   !c_h2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  k_co2 = EXP(-10600*x**4 + 6356*x**3 - 1483.6*x**2 + 178.62*x + 10.142)
! 	  return
! 	  end function k_co2
! 	  
! 	  double precision function e_h2o ( x )   !c_h2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  e_h2o = -(182719*x**5 - 121066*x**4 + 33580*x**3 - 5179.6*x**2 + 516.55*x + 18.52)*1000.d0
! 	  return
! 	  end function e_h2o
! 	  
! 	  double precision function e_co2 ( x )   !c_h2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  e_co2 = -(122746*x**5 - 91677*x**4 + 28920*x**3 - 5005.3*x**2 + 535.39*x + 35.555)*1000.d0
! 	  return
! 	  end function e_co2
! 
! 	  double precision function napp_h2o ( x )  !c_h2o
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  napp_h2o = -0.5328*x + 0.9927
! 	  return
! 	  end function napp_h2o
! 	  
! 	  double precision function napp_co2 ( x )  !c_co2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  napp_co2 = -0.6936*x + 0.9972
! 	  return
! 	  end function napp_co2
! 	  
! 	  double precision function kh2 ( x )   !c_co
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  kh2 = 11.751*x**2 - 5.5042*x + 1.3411
! 	  return
! 	  end function kh2
! 	  
! 	  double precision function kco ( x )   !c_h2
! 	  double precision x
! 	  if(x .gt. 0.2d0)then
! 	     x = 0.2d0
! 	  endif
! 	  kco = -7850.8*x**5 + 4641.3*x**4 - 1061*x**3 + 121.11*x**2 - 7.8113*x + 0.3419
! 	  return
! 	  end function kco
! 	  
!       END SUBROUTINE USR_RATES
!       
!   
! 
! 	  
! 
! 
