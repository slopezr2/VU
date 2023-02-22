!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!#######################################################################
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_Param

  implicit none


  ! --- in/out ---------------------------

  public


  ! --- const ------------------------------

  ! maximum number of output sets:
  integer, parameter  ::  MAORI_MAX_SET = 32

  ! maximum length for names etc:
  integer, parameter  ::  MAORI_LEN_NAME     =  32
  integer, parameter  ::  MAORI_LEN_LONGNAME = 128
  integer, parameter  ::  MAORI_LEN_UNIT     =  32

  ! maximum expected line lengths in rcfiles:
  integer, parameter  ::  MAORI_LEN_LINE = 512

  ! maximum length of file names:
  integer, parameter  ::  MAORI_LEN_FILE = 512

  ! data types:
  integer, parameter  ::  MAORI_INT       =  1
  integer, parameter  ::  MAORI_REAL      =  2
  integer, parameter  ::  MAORI_CHAR      =  3
  integer, parameter  ::  MAORI_SAMPLE    =  4
  integer, parameter  ::  MAORI_PROFILE   =  5
  integer, parameter  ::  MAORI_FIELD     =  6
  integer, parameter  ::  MAORI_SATELLITE =  7
  !
  ! valid range:
  integer, parameter  ::  MAORI_TYPE_MIN = 1
  integer, parameter  ::  MAORI_TYPE_MAX = 7
  ! type names:
  character(len=MAORI_LEN_NAME), parameter  ::  &
      MAORI_TYPE_NAME(MAORI_TYPE_MIN:MAORI_TYPE_MAX) = &
        (/ 'int      ', &
           'real     ', &
           'char     ', &
           'sample   ', &
           'profile  ', &
           'field    ', &
           'satellite' /)
           
  ! maximum number of variables:
  integer, parameter  ::  MAORI_MAX_VAR = 32
  
  ! value used for no-data :
  real, parameter     ::  MAORI_NAN = -999.9

  !
  ! Assimilation flags:
  !
  !   000
  !   001    ibset(0,0)   ! position 0, least significant bit, 'right'
  !   010    ibset(0,1)   ! position 1
  !   011
  !   100    ibset(0,2)   ! position 2, most significant bit, 'left'
  !   101
  !   110
  !   111
  !
  ! use bit-by-bit functions using positions:
  !
  !   integer, parameter  ::  bpos_outside     = 0
  !   integer, parameter  ::  bpos_nodata      = 2
  !   integer, parameter  ::  bpos_validation  = 3
  !
  ! Set and clear the bits using:
  !   status = 0
  !   status = IBSet( status, bpos_nodata )
  !   status = IBClr( status, bpos_nodata )
  ! Test if a bit is set:
  !   if ( BTest( status, bpos_nodata ) ) ...
  !
  integer, parameter            ::  MAORI_ASTAT_OUTSIDE    = 0  ! +1
  integer, parameter            ::  MAORI_ASTAT_NODATA     = 1  ! +2
  integer, parameter            ::  MAORI_ASTAT_VALIDATION = 2  ! +4
  integer, parameter            ::  MAORI_ASTAT_SCREENED   = 3  ! +8
  integer, parameter            ::  MAORI_ASTAT_ANALYSED   = 4  ! +16
  !
  character(len=*), parameter   ::  MAORI_ASTAT_DESCRIPTION = &
      'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation, +8=screened, +16=analysed'


end module MAORI_Param
