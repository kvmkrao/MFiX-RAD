SUBROUTINE USR_READ_NAMELIST(STRING, LINE_STRING, LINE_LEN, IOS)

#ifdef USR_NAMELIST
   use usr
#endif

   IMPLICIT NONE

! Holds one line in the input file
   CHARACTER(LEN=512), INTENT(INOUT) :: LINE_STRING
! Length of noncomment string
   INTEGER, INTENT(INOUT) :: LINE_LEN
   CHARACTER(len=256), INTENT(INOUT) :: STRING
   INTEGER, INTENT(INOUT) :: IOS

#ifdef USR_NAMELIST
   INCLUDE  'usrnlst.inc'

! User defined input parameters.
   STRING=''; STRING = '&USR_INPUT_DATA '//&
      trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
   READ(STRING, NML=USR_INPUT_DATA, IOSTAT=IOS)
#endif

END SUBROUTINE USR_READ_NAMELIST
