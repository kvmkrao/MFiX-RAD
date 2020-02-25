!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR2                                                   C
!  Purpose: This routine is called from the outer iteration loop and   C
!           is user-definable.  The user may insert code in this       C
!           routine or call appropriate user defined subroutines. This C
!           can be used for updating quantities every iteration.  Use  C
!           this routine sparingly considering the large computational C
!           over head. This routine is not called from an              C
!           IJK loop, hence all indices are undefined.                 C               C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR2
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      Use usr
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc
      USE constant
      USE compar
      USE funits
      USE fun_avg
      USE usr
      USE functions
      use rad
!      use rad_fields, only: Srad
      use rad_fields, only: k_g, k_s, G, E
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
      INTEGER IJK, M, I, IMJK, IJMK, IJKM,k,j
      DOUBLE PRECISION DIFF, EP_g2
      DOUBLE PRECISION Sc1o3, UGC, VGC, WGC, USCM, VSCM, WSCM, VREL, Re
      DOUBLE PRECISION Tg_MFIX,xt
      INCLUDE 'usrnlst.inc'
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!-- rad_changes_begin
!      call rad_calc()
!      call rad_write_src()

      RETURN
      END SUBROUTINE USR2
