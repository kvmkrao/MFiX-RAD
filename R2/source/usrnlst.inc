!
!  First, define the user-defined namelist variables as DOUBLE PRECISION,
!  INTEGER, LOGICAL, or CHARACTER.  Second, include the variable in the
!  appropriate common.  Third, add the variable name in the NAMELIST section.
!  A dummy variable is defined, so that the code can be compiled with the
!  template include file.
!
!
!                      Proximate Analysis
      DOUBLE PRECISION PAFC, PAA

!  Radiation Keywords
      LOGICAL RAD_ON
      DOUBLE PRECISION RAD_EMIS_W(DIMENSION_BC), RAD_T_W(DIMENSION_BC)
      INTEGER RAD_NQUAD, RAD_SKIP, RAD_NRR
      CHARACTER(30) RAD_RTE, RAD_SPECTRAL
      DOUBLE PRECISION const_sfp, const_emisp,const_absg
!
!     Double Precision variables
      COMMON / USR_DATA_DP / &
             PAFC, const_absg, const_sfp, const_emisp, PAA, &
             RAD_EMIS_W, RAD_T_W
!
!     Integer variables
      COMMON / USR_DATA_I / &
		RAD_NQUAD, RAD_SKIP, RAD_NRR
!
!     Logical variables
      COMMON / USR_DATA_L / &
	     RAD_ON
!
!     Character variables
      COMMON / USR_DATA_C / &
         RAD_RTE, RAD_SPECTRAL
!
      NAMELIST / USR_INPUT_DATA / &
            PAFC, const_sfp, const_emisp, const_absg,PAA, & 
	    RAD_ON, RAD_EMIS_W, RAD_T_W, RAD_NQUAD, &
            RAD_SKIP, RAD_NRR, RAD_RTE, RAD_SPECTRAL
