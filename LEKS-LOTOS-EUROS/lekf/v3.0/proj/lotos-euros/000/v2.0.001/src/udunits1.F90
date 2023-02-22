!#################################################################
!
! Fortran module around UDUnits .
!
! USAGE
!
!   use UDUnits1
!
!   integer                     ::  status
!   integer(UD_POINTER_KIND)    ::  unit, unit2
!   character(len=64)           ::  spec, spec2
!   real(8)                     ::  slope, offset
!
!   ! * module initialisation
!
!   call UDUnits1_Init( status )
!   if (status/=UDUNITS_NOERR) then
!     print *, trim(UDUnits1_StrError(status))
!     stop
!   end if
!
!   ! * high levell routines
!  
!   spec = 'kg/s'
!   call UDUnits1_Standard( spec, spec2, status )
!   write (*,'("standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!  
!   spec = 'gram/cm3' ; spec2 = 'kg/m3'
!   call UDUnits1_ConversionFactor( spec, spec2, slope, status )
!   write (*,'("conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!   ! * low level routines
!
!   call UDUnits1_Make( unit, status )
!   call UDUnits1_Make( unit2, status )
!
!   call UDUnits1_Decode( 'kg', unit, status )
!   call UDUnits1_Encode( unit, spec, status )
!  
!   call UDUnits1_Convert( unit, unit2, slope, offset, status )
!  
!   ! * done with module
!
!   call UDUnits1_Done( status )
!  
! HISTORY
!   2010 feb, Arjo Segers, JRC
!
!
!#################################################################

module UDUnits1

  use UDUnits1_Inc, only : UD_POINTER_KIND

  implicit none
  

  ! --- in/out -----------------------------------
  
  private
  
  !public    ::  UD_POINTER_KIND
  public    ::  UDUNITS1_NOERR, UDUnits1_StrError

  public    ::  UDUnits1_Init, UDUnits1_Done

  !public    ::  UDUnits1_Make
  !public    ::  UDUnits1_Decode, UDUnits1_Encode
  !public    ::  UDUnits1_Convert
  
  public    ::  UDUnits1_Standard
  public    ::  UDUnits1_ConversionFactor


  ! --- const --------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'UDUnits1'

  ! name of environment variable with path to data file:
  character(len=*), parameter ::  env_var = 'UDUNITS_PATH'

!  ! unit should be of type :  integer(UD_POINTER_KIND)
!  integer, parameter  ::  UD_POINTER_KIND = 4
  
  ! no error:
  integer, parameter  ::  UDUNITS1_NOERR = 0
  
!  ! error codes:
!  integer, parameter  ::  UT_EOF      =   1
!  integer, parameter  ::  UT_ENOFILE  =  -1
!  integer, parameter  ::  UT_ESYNTAX  =  -2
!  integer, parameter  ::  UT_EUNKNOWN =  -3
!  integer, parameter  ::  UT_EIO      =  -4
!  integer, parameter  ::  UT_EINVALID =  -5
!  integer, parameter  ::  UT_ENOINIT  =  -6
!  integer, parameter  ::  UT_ECONVERT =  -7
!  integer, parameter  ::  UT_EALLOC   =  -8
!  integer, parameter  ::  UT_ENOROOM  =  -9
!  integer, parameter  ::  UT_ENOTTIME = -10
!  
  !integer, parameter  ::  UT_MAXNUM_BASE_QUANTITIES = 10
  
  ! storage for latest error:
  integer, parameter  ::  error_status  = -100
  character(len=256)  ::  error_message = ''
  
  ! maximum length of specifications:
  integer, parameter  ::  spec_len = 64



contains


  ! ====================================================================
  ! ===
  ! === module routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits1_Init( status )
  
    !use UDUnits1_Inc, only : udunits_inc_test
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits1_Init'
    
    ! --- external -------------------------------
    
    ! Initialize the units package:
    integer, external  :: UTOpen

    ! --- local ----------------------------------
    
    character(len=256)        ::  UDUNITS_PATH
    integer                   ::  length
    
    ! --- begin ----------------------------------
    
    !call udunits_inc_test()
    
    ! following the manuals, the path to the udunits data file is
    ! taken from the environment variable UDUNITS_PATH if not specified;
    ! this does not seem to work properly however, and therefore
    ! the path is explicitly taken from the environment:
    call Get_Environment_Variable( env_var, UDUNITS_PATH, length, status )
    if (status/=0) then
      write (error_message,'("could not get environment variable `",a,"`")') trim(env_var)
      status=error_status; return
    end if
    
    ! Initialize the units package:
    status = UTOpen( trim(UDUNITS_PATH) )
    if (status/=0) write (error_message,'("could not initialize from data file `",a,"`")') trim(UDUNITS_PATH)
    
  end subroutine UDUnits1_Init
  
  
  ! ***


  subroutine UDUnits1_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- external -------------------------------
    
    ! --- begin ----------------------------------
    
    ! function UTFree not available in Fortran interface ...
    
    ! ok
    status = 0
    
  end subroutine UDUnits1_Done
  
  
  ! ====================================================================
  ! ===
  ! === error messages
  ! ===
  ! ====================================================================
  
  
  function UDUnits1_StrError( status )
  
    use UDUnits1_Inc, only : UT_EOF, UT_ENOFILE, UT_ESYNTAX, UT_EUNKNOWN, &
                            UT_EIO, UT_EINVALID, UT_ENOINIT, UT_ECONVERT, &
                            UT_EALLOC, UT_ENOROOM, UT_ENOTTIME

    ! --- in/out ---------------------------------
    
    character(len=256)          ::  UDUnits1_StrError
    integer, intent(inout)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits1_StrError'
    
    ! --- begin ----------------------------------
    
    ! no error?
    if ( status == UDUNITS1_NOERR ) then
      ! empty message:
      UDUnits1_StrError = ''
    else
      ! fill message:
      select case ( status )
        ! supported:
        case ( UT_EOF        ) ; UDUnits1_StrError = 'End of file'
        case ( UT_ENOFILE    ) ; UDUnits1_StrError = 'Units file does not exist'
        case ( UT_ESYNTAX    ) ; UDUnits1_StrError = 'Syntax error'
        case ( UT_EUNKNOWN   ) ; UDUnits1_StrError = 'Unknown unit specification'
        case ( UT_EIO        ) ; UDUnits1_StrError = 'I/O error while accessing the units file'
        case ( UT_EINVALID   ) ; UDUnits1_StrError = 'Invalid value'
        case ( UT_ENOINIT    ) ; UDUnits1_StrError = 'Package has not be initialized'
        case ( UT_ECONVERT   ) ; UDUnits1_StrError = 'Conversion error'
        case ( UT_EALLOC     ) ; UDUnits1_StrError = 'Memory allocation failure'
        case ( UT_ENOROOM    ) ; UDUnits1_StrError = 'No room for result'
        case ( UT_ENOTTIME   ) ; UDUnits1_StrError = 'No time value'
        ! other ...
        case ( error_status ) ; UDUnits1_StrError = ''
        ! unknown:
        case default
          write (UDUnits1_StrError,'("Unknown error status from UDUnits routine : ",i6)') status
      end select
      ! add error buffer:
      if ( status /= 0 ) then
        if ( len_trim(UDUnits1_StrError) == 0 ) then
          UDUnits1_StrError = trim(error_message)
        else
          UDUnits1_StrError = trim(UDUnits1_StrError)//'; '//trim(error_message)
        end if
      end if
    end if  ! error status
    
  end function UDUnits1_StrError


  ! ====================================================================
  ! ===
  ! === low level routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits1_Make( unit, status )
  
    ! --- in/out ---------------------------------
    
    integer(UD_POINTER_KIND), intent(out)   ::  unit
    integer, intent(out)            ::  status
    
    ! --- external -------------------------------
    
    ! set return status:
    integer(UD_POINTER_KIND), external  :: UTMake

    ! --- begin ----------------------------------
    
    ! Create a new unit:
    unit = UTMake()
    
    ! set return status:
    status = 0
    if ( unit < 0 ) status = int(unit)
    
  end subroutine UDUnits1_Make
  
  
  ! ***


  subroutine UDUnits1_Decode( spec, unit, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec
    integer(UD_POINTER_KIND), intent(in)      ::  unit
    integer, intent(out)              ::  status
    
    ! --- external -------------------------------
    
    ! Decode a formatted specification into a unit:
    integer, external  :: UTDec

    ! --- begin ----------------------------------
    
    ! Decode a formatted specification into a unit:
    status = UTDec( spec, unit )
    if (status/=0) write (error_message,'("could not decode `",a,"`")') trim(spec)
    
  end subroutine UDUnits1_Decode
  
  
  ! ***


  subroutine UDUnits1_Encode( unit, spec, status )
  
    ! --- in/out ---------------------------------
    
    integer(UD_POINTER_KIND), intent(in)      ::  unit
    character(len=*), intent(out)     ::  spec
    integer, intent(out)              ::  status
    
    ! --- external -------------------------------
    
    ! Encode a unit into a formatted specification:
    integer, external  :: UTEnc

    ! --- begin ----------------------------------
    
    ! Encode a unit into a formatted specification:
    status = UTEnc( unit, spec )
    if (status/=0) write (error_message,'("could not encode from unit into formatted specification")')
    
  end subroutine UDUnits1_Encode
  
  
  ! ***


  subroutine UDUnits1_Convert( unit_from, unit_to, slope, intercept, status )
  
    ! --- in/out ---------------------------------
    
    integer(UD_POINTER_KIND), intent(in)      ::  unit_from
    integer(UD_POINTER_KIND), intent(in)      ::  unit_to
    real(8), intent(out)              ::  slope, intercept
    integer, intent(out)              ::  status
    
    ! --- external -------------------------------
    
    ! Convert from one unit to another:
    integer, external  :: UTCvt
    
    ! --- local ----------------------------------
    
    character(len=spec_len)   ::  spec_from, spec_to

    ! --- begin ----------------------------------
    
    ! Convert from one unit to another:
    status = UTCvt( unit_from, unit_to, slope, intercept )
    if (status/=0) then
      call UDUnits1_Encode( unit_from, spec_from, status )
      if (status/=0) then
        write (error_message,'("could not convert units; failed to convert unit_from to specification")')
        status = error_status; return
      end if
      call UDUnits1_Encode( unit_to, spec_to, status )
      if (status/=0) then
        write (error_message,'("could not convert from `",a,"`; failed to convert unit_to to specification")') trim(spec_from)
        status = error_status; return
      end if
      write (error_message,'("could not convert from `",a,"` to `",a,"`")') trim(spec_from), trim(spec_to)
      status = error_status; return
    end if
    
  end subroutine UDUnits1_Convert
  
  
  ! ====================================================================
  ! ===
  ! === high level routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits1_Standard( spec_from, spec_to, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(out)     ::  spec_to
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits1_Standard'
    
    ! --- local ----------------------------------
    
    integer(UD_POINTER_KIND)      ::  unit_from

    ! --- begin ----------------------------------
    
    ! setup unit:
    call UDUnits1_Make( unit_from, status )
    if (status/=0) return
    ! fill with secification:
    call UDUnits1_Decode( spec_from, unit_from, status )
    if (status/=0) return
    ! extract standard name:
    call UDUnits1_Encode( unit_from, spec_to, status )
    if (status/=0) return
    
    ! ok
    status = 0
    
  end subroutine UDUnits1_Standard
  
  
  ! ***
  
  
  subroutine UDUnits1_ConversionFactor( spec_from, spec_to, factor, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(in)      ::  spec_to
    real(8), intent(out)              ::  factor
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits1_ConversionFactor'
    
    ! --- local ----------------------------------
    
    integer(UD_POINTER_KIND)      ::  unit_from
    integer(UD_POINTER_KIND)      ::  unit_to
    real(8)                       ::  offset

    ! --- begin ----------------------------------
    
    ! input unit:
    call UDUnits1_Make( unit_from, status )
    if (status/=0) return
    call UDUnits1_Decode( spec_from, unit_from, status )
    if (status/=0) return
    
    ! output unit:
    call UDUnits1_Make( unit_to, status )
    if (status/=0) return
    call UDUnits1_Decode( spec_to, unit_to, status )
    if (status/=0) return
    
    ! Convert from one unit to another:
    call UDUnits1_Convert( unit_from, unit_to, factor, offset, status )
    if (status/=0) return
    
    ! check ...
    if ( offset /= 0.0d0 ) then
      write (error_message,*) 'found conversion offset unequal to zero : ', offset
      status=error_status; return
    end if
    
    ! ok
    status = 0
    
  end subroutine UDUnits1_ConversionFactor
    

end module UDUnits1


!! ######################################################################
!! ###
!! ### test
!! ###
!! ######################################################################
!
!! f90 -o test.x udunits_inc.F udunits.F90 -I${UDUNITS_HOME}/include -L${UDUNITS_HOME}/lib -ludunits  &&  ./test.x
!
!program test_udunits1
!
!  use UDUnits
!
!  implicit none
!
!  integer                     ::  status
!  integer(UD_POINTER_KIND)    ::  unit, unit2
!  character(len=64)           ::  spec, spec2
!  real(8)                     ::  slope, offset
!
!  write (*,'("begin test_udunits")')
!  
!  write (*,'("   UD_POINTER_KIND : ",i4)') UD_POINTER_KIND
!
!  ! * module initialisation
!
!  call UDUnits1_Init( status )
!  if (status/=UDUNITS_NOERR) then
!    print *, trim(UDUnits1_StrError(status))
!    stop
!  end if
!
!  ! * high levell routines
!
!  spec = 'kg/s'
!  call UDUnits1_Standard( spec, spec2, status )
!  write (*,'("  standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!
!  spec = 'gram/cm3' ; spec2 = 'kg/m3'
!  call UDUnits1_ConversionFactor( spec, spec2, slope, status )
!  write (*,'("  conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!  ! * low level routines
!
!  call UDUnits1_Make( unit, status )
!  call UDUnits1_Make( unit2, status )
!
!  call UDUnits1_Decode( 'kg', unit, status )
!  call UDUnits1_Encode( unit, spec, status )
!
!  call UDUnits1_Convert( unit, unit2, slope, offset, status )
!
!  ! * done with module
!
!  call UDUnits1_Done( status )
!
!  ! *
!
!  write (*,'("end test_udunits")')
!
!end program test_udunits1
