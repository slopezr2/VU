!#################################################################
!
! Fortran module around UDUnits2 .
!
! USAGE
!
!   use UDUnits2, only : UDUNITS2_NOERR
!   use UDUnits2, only : UDUnits2_StrError
!   use UDUnits2, only : UDUnits2_Init
!   use UDUnits2, only : UDUnits2_Done
!
!   integer                     ::  status
!   character(len=64)           ::  spec, spec2
!   real(8)                     ::  factor
!
!   ! * module initialisation
!
!   call UDUnits2_Init( status )
!   if (status/=UDUNITS2_NOERR) then
!     print *, trim(UDUnits2_StrError(status))
!     stop
!   end if
!
!   ! * high levell routines
!  
!   spec = 'kg/s'
!   call UDUnits2_Standard( spec, spec2, status )
!   write (*,'("standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!  
!   spec = 'gram/cm3' ; spec2 = 'kg/m3'
!   call UDUnits2_ConversionFactor( spec, spec2, slope, status )
!   write (*,'("conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!   ! * done with module
!
!   call UDUnits2_Done( status )
!  
! HISTORY
!   2015-12, Arjo Segers, TNO
!
!
!### macro's #####################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
!
#define IF_UDUNITS2_NOT_OK_RETURN(action) call UDUnits2_Get_Status(status); if (status/=0) then; action; return; end if
!
!#################################################################

module UDUnits2

  use udunits2_f, only : UT_SYSTEM_PTR

  implicit none
  

  ! --- in/out -----------------------------------
  
  private
  
  public    ::  UDUNITS2_NOERR, UDUnits2_StrError

  public    ::  UDUnits2_Init, UDUnits2_Done
  
  public    ::  UDUnits2_Standard
  public    ::  UDUnits2_ConversionFactor


  ! --- const --------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'UDUnits2'

  ! no error:
  integer, parameter  ::  UDUNITS2_NOERR = 0
  

  ! --- var --------------------------------

  ! unit system:
  type(UT_SYSTEM_PTR)     ::  unitSystem

  ! storage for latest error:
  integer, parameter  ::  error_status  = -100
  character(len=256)  ::  error_message = ''
  

contains


  ! ====================================================================
  ! ===
  ! === module routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits2_Init( status )
  
    use udunits2_f, only : f_ut_read_xml
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! read default unit system:
    unitSystem = f_ut_read_xml('')
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine UDUnits2_Init
  
  
  ! ***


  subroutine UDUnits2_Done( status )

    use udunits2_f, only : f_ut_free_system
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_Done'
    
    ! --- begin ----------------------------------
    
    ! done with system:
    call f_ut_free_system( unitSystem )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine UDUnits2_Done
  
  
  ! ====================================================================
  ! ===
  ! === error messages
  ! ===
  ! ====================================================================
  
  subroutine UDUnits2_Get_Status( status )
  
    use udunits2_f, only : f_ut_get_status
    use udunits2_f, only : UT_SUCCESS
    use udunits2_f, only : UT_BAD_ARG
    use udunits2_f, only : UT_EXISTS
    use udunits2_f, only : UT_NO_UNIT
    use udunits2_f, only : UT_OS
    use udunits2_f, only : UT_NOT_SAME_SYSTEM
    use udunits2_f, only : UT_MEANINGLESS
    use udunits2_f, only : UT_NO_SECOND
    use udunits2_f, only : UT_VISIT_ERROR
    use udunits2_f, only : UT_CANT_FORMAT
    use udunits2_f, only : UT_SYNTAX
    use udunits2_f, only : UT_UNKNOWN
    use udunits2_f, only : UT_OPEN_ARG
    use udunits2_f, only : UT_OPEN_ENV
    use udunits2_f, only : UT_OPEN_DEFAULT
    use udunits2_f, only : UT_PARSE
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_Get_Status'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! get status of latest command:
    status = f_ut_get_status()
    ! error ?
    if ( status /= UT_SUCCESS ) then
      ! fill mesage:
      select case ( status )
        ! supported:
        case ( UT_SUCCESS         ) ; error_message = 'Success'
        case ( UT_BAD_ARG         ) ; error_message = 'An argument violates the the function''s contract (e.g., it''s NULL).'
        case ( UT_EXISTS          ) ; error_message = 'Unit, prefix, or identifier already exists'
        case ( UT_NO_UNIT         ) ; error_message = 'No such unit exists'
        case ( UT_OS              ) ; error_message = 'Operating-system error. See errno for the reason.'
        case ( UT_NOT_SAME_SYSTEM ) ; error_message = 'The units belong to different unit-systems'
        case ( UT_MEANINGLESS     ) ; error_message = 'The operation on the unit or units is meaningless'
        case ( UT_NO_SECOND       ) ; error_message = 'The unit-system doesn''t have a unit named "second"'
        case ( UT_VISIT_ERROR     ) ; error_message = 'An error occurred while visiting a unit'
        case ( UT_CANT_FORMAT     ) ; error_message = 'A unit can''t be formatted in the desired manner'
        case ( UT_SYNTAX          ) ; error_message = 'String unit representation contains syntax error'
        case ( UT_UNKNOWN         ) ; error_message = 'String unit representation contains unknown word'
        case ( UT_OPEN_ARG        ) ; error_message = 'Can''t open argument-specified unit database'
        case ( UT_OPEN_ENV        ) ; error_message = 'Can''t open environment-specified unit database'
        case ( UT_OPEN_DEFAULT    ) ; error_message = 'Can''t open installed, default, unit database'
        case ( UT_PARSE           ) ; error_message = 'Error parsing unit database'
        ! unknown:
        case default
          write (error_message,'("Unknown error status from UDUnits routine : ",i6)') status
      end select
      ! break:
      status=error_status; return
    end if ! error status
    
    ! ok
    status = 0
    
  end subroutine UDUnits2_Get_Status
  
  ! *
  
  function UDUnits2_StrError( status )

    ! --- in/out ---------------------------------
    
    character(len=256)          ::  UDUnits2_StrError
    integer, intent(inout)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_StrError'
    
    ! --- begin ----------------------------------
    
    ! copy latest error message:
    UDUnits2_StrError = trim(error_message)
    
  end function UDUnits2_StrError
  

  ! ====================================================================
  ! ===
  ! === high level routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits2_Standard( spec_from, spec_to, status )
  
    use udunits2_f, only : UT_UNIT_PTR
    use udunits2_f, only : UT_ASCII
    use udunits2_f, only : UT_DEFINITION
    use udunits2_f, only : f_ut_parse
    use udunits2_f, only : f_ut_format
    use udunits2_f, only : f_ut_free
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(out)     ::  spec_to
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_Standard'
    
    ! --- local ----------------------------------
    
    type(UT_UNIT_PTR)      ::  unit_from
    integer                ::  length

    ! --- begin ----------------------------------
    
    ! fill with secification:
    unit_from = f_ut_parse( unitSystem, spec_from, UT_ASCII )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    
    ! extract standard name:
    length = f_ut_format( unit_from, spec_to, UT_ASCII+UT_DEFINITION )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)

    ! clear:
    call f_ut_free( unit_from )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine UDUnits2_Standard
  
  
  ! ***
  
  
  subroutine UDUnits2_ConversionFactor( spec_from, spec_to, factor, status )
  
    use udunits2_f, only : UT_UNIT_PTR
    use udunits2_f, only : UT_ASCII
    use udunits2_f, only : f_ut_parse
    use udunits2_f, only : f_ut_free
    use udunits2_f, only : CV_CONVERTER_PTR
    use udunits2_f, only : f_ut_get_converter
    use udunits2_f, only : f_cv_free
    use udunits2_f, only : f_cv_convert_double
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(in)      ::  spec_to
    real(8), intent(out)              ::  factor
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits2_ConversionFactor'
    
    ! --- local ----------------------------------
    
    type(UT_UNIT_PTR)      ::  unit_from
    type(UT_UNIT_PTR)      ::  unit_to
    type(CV_CONVERTER_PTR) ::  convertor
    real(8)                ::  x1, x2
    real(8)                ::  y1, y2
    real(8)                ::  offset

    ! --- begin ----------------------------------

    ! input unit:
    unit_from = f_ut_parse( unitSystem, spec_from, UT_ASCII )
    call UDUnits2_Get_Status( status )
    if ( status /= 0 ) then
      write (error_message,'(a,"; could not parse: ",a)') trim(error_message), trim(spec_from)
      status=error_status; return
    end if

    ! output unit:
    unit_to = f_ut_parse( unitSystem, spec_to, UT_ASCII )
    call UDUnits2_Get_Status( status )
    if ( status /= 0 ) then
      write (error_message,'(a,"; could not parse: ",a)') trim(error_message), trim(spec_to)
      status=error_status; return
    end if
    
    ! create convertor:
    convertor = f_ut_get_converter( unit_from, unit_to )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)

    ! original values:
    x1 = 1.0d0
    x2 = 2.0d0
    ! convert unit numbers:
    y1 = f_cv_convert_double( convertor, x1 )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    y2 = f_cv_convert_double( convertor, x2 )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    ! assume linear conversion:
    !   y = offset + factor * x
    factor = ( y2 - y1 ) / ( x2 - x1 )
    offset = y1 - factor * x1

    ! clear:
    call f_cv_free( convertor )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)

    ! clear:
    call f_ut_free( unit_to )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    ! clear:
    call f_ut_free( unit_from )
    IF_UDUNITS2_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( abs(offset) > 1.0e-6 ) then
      write (error_message,*) 'found conversion offset unequal to zero : ', offset
      status=error_status; return
    end if
    
    ! ok
    status = 0
    
  end subroutine UDUnits2_ConversionFactor
    

end module UDUnits2


!! ######################################################################
!! ###
!! ### test
!! ###
!! ######################################################################
!!
!#define IF_UDUNITS2_NOT_OK_STOP if (status/=0) then; write (*,'(a)') trim(UDUnits2_StrError(status)); TRACEBACK; stop; end if
!!
!! f90 -o test.x udunits2.F90 \
!!     -I${UDUNITS2_HOME}/include \
!!     -L${UDUNITS2_HOME}/lib -ludunits2  &&  ./test.x
!!
!
!program Test_UDUnits2
!
!  use UDUnits2
!
!  implicit none
!
!  character(len=*), parameter ::  rname = 'Test_UDUnits2'
!
!  integer                     ::  status
!  character(len=64)           ::  spec, spec2
!  real(8)                     ::  slope, offset
!
!  ! info ...
!  write (*,'(a,": ")') rname
!  write (*,'(a,": ** begin")') rname
!  write (*,'(a,": ")') rname
!  
!  ! * module initialisation
!
!  call UDUnits2_Init( status )
!  IF_UDUNITS2_NOT_OK_STOP
!
!  ! * high levell routines
!
!  spec = 'kg/s'
!  call UDUnits2_Standard( spec, spec2, status )
!  IF_UDUNITS2_NOT_OK_STOP
!  write (*,'(a,": standard name of `",a,"` is `",a,"`")') &
!                      rname, trim(spec), trim(spec2)
!
!  spec = 'gram/cm3' ; spec2 = 'kg/m3'
!  call UDUnits2_ConversionFactor( spec, spec2, slope, status )
!  IF_UDUNITS2_NOT_OK_STOP
!  write (*,'(a,": conversion factor from `",a,"` to `",a,"` is ",f12.4)') &
!                      rname, trim(spec), trim(spec2), slope
!
!  ! * done with module
!
!  call UDUnits2_Done( status )
!  IF_UDUNITS2_NOT_OK_STOP
!
!  ! *
!
!  ! info ...
!  write (*,'(a,": ")') rname
!  write (*,'(a,": ** end")') rname
!  write (*,'(a,": ")') rname
!
!end program Test_UDUnits2
