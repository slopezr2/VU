!#################################################################
!
! Fortran module around UDUnits libraries (version 1 and 2) .
!
! VERSION WRAPPER
!
!   This module is a wrapper to either udunits1 or UDUnits.
!   Which interface is selected depends on macro definitions:
!
!     with_udunits1     : version 1
!     with_udunits2     : version 2
!
!   The macro definition should be placed in the 
!   macro include file "udunits_version.inc":
!
!     #define with_udunits2
!
!   If none of these macro's is defined, the code can still be compiled,
!   and the module can be initialized and finalized. However, calls
!   to the work routines will result in a return with a error status.
!
!   If both macro's are defined, the UDUnits2 library has priority.
!
!
! USAGE
!
!   use UDUnits, only : UDUnits_NOERR
!   use UDUnits, only : UDUnits_StrError
!   use UDUnits, only : UDUnits_Init
!   use UDUnits, only : UDUnits_Done
!
!   integer                     ::  status
!   character(len=64)           ::  spec, spec2
!   real(8)                     ::  factor
!
!   ! * module initialisation
!
!   call UDUnits_Init( status )
!   if (status/=UDUnits_NOERR) then
!     print *, trim(UDUnits_StrError(status))
!     stop
!   end if
!
!   ! * high levell routines
!  
!   ! current units:
!   spec = 'kg/s'
!   ! convert to standard units:
!   call UDUnits_Standard( spec, spec2, status )
!   if (status/=UDUnits_NOERR) stop 'error from UDUnits'
!   ! info ...
!   write (*,'("standard name of `",a,"` is `",a,"`")') trim(spec), trim(spec2)
!  
!   ! input and target units:
!   spec = 'gram/cm3' ; spec2 = 'kg/m3'
!   ! obtain conversion factor:
!   call UDUnits_ConversionFactor( spec, spec2, slope, status )
!   if (status/=UDUnits_NOERR) stop 'error from UDUnits'
!   ! info:
!   write (*,'("conversion factor from `",a,"` to `",a,"` is ",f12.4)') trim(spec), trim(spec2), slope
!
!   ! * done with module
!
!   call UDUnits_Done( status )
!   if (status/=UDUnits_NOERR) stop 'error from UDUnits'
!
!
! CONTENT
!
!   The module consists of the following files:
!
!     udunits.F90          : main module
!     udunits_version.inc  : macro include file for version selection
!
!       udunits1.F90         : high level interface to udunits v1 library
!         udunits1_inc.F       : interface to library include file
!
!       udunits2.F90         : high level interface to udunits v2 library
!         udunits2_f.f90       : fortran-to-c wrapper for udunits2
!         udunits2_f.inc       :  "
!
!
! HISTORY
!   2015-12, Arjo Segers, TNO
!     Original version.
!
!### macro's #####################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; action; return; end if
!
#include "udunits_version.inc"
!
!#################################################################


module UDUnits

  implicit none
  

  ! --- in/out -----------------------------------
  
  private
  
  public    ::  UDUNITS_NOERR, UDUnits_StrError

  public    ::  UDUnits_Init, UDUnits_Done
  
  public    ::  UDUnits_Standard
  public    ::  UDUnits_ConversionFactor


  ! --- const --------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'UDUnits'

  ! error states:
  integer, parameter  ::  UDUNITS_NOERR = 0
  integer, parameter  ::  UDUNITS_NOT_INITIALIZED     = 1001
  integer, parameter  ::  UDUNITS_ALREADY_INITIALIZED = 1002
  integer, parameter  ::  UDUNITS_NO_DEFS             = 1003
  

  ! --- var --------------------------------
  
  ! initialized flag:
  logical             ::  initialized = .false.

  ! storage for latest error:
  character(len=256)  ::  error_message = ''
  

contains


  ! ====================================================================
  ! ===
  ! === module routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits_Init( status )

#ifdef with_udunits2
    use UDUnits2, only : UDUnits2_Init, UDUNITS2_NOERR
#else
#ifdef with_udunits1
    use UDUnits1, only : UDUnits1_Init, UDUNITS1_NOERR
#endif
#endif
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! check ..
    if ( initialized ) then; status=UDUNITS_ALREADY_INITIALIZED; return; end if
    
#ifdef with_udunits2
    call UDUnits2_Init( status )
    if (status/=UDUNITS2_NOERR) return
#else
#ifdef with_udunits1
    call UDUnits1_Init( status )
    if (status/=UDUNITS1_NOERR) return
#endif
#endif

    ! reset flag:
    initialized = .true.
    
    ! ok
    status = 0
    
  end subroutine UDUnits_Init
  
  
  ! ***


  subroutine UDUnits_Done( status )

#ifdef with_udunits2
    use UDUnits2, only : UDUnits2_Done, UDUNITS2_NOERR
#else
#ifdef with_udunits1
    use UDUnits1, only : UDUnits1_Done, UDUNITS1_NOERR
#endif
#endif
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_Done'
    
    ! --- begin ----------------------------------

    ! check ...
    if ( .not. initialized ) then; status=UDUNITS_NOT_INITIALIZED; return; end if

#ifdef with_udunits2
    call UDUnits2_Done( status )
    if (status/=UDUNITS2_NOERR) return
#else
#ifdef with_udunits1
    call UDUnits1_Done( status )
    if (status/=UDUNITS1_NOERR) return
#endif
#endif

    ! reset flag:
    initialized = .false.
    
    ! ok
    status = 0
    
  end subroutine UDUnits_Done
  

  ! ====================================================================
  ! ===
  ! === error messages
  ! ===
  ! ====================================================================
  

  function UDUnits_StrError( status )

#ifdef with_udunits2
    use UDUnits2, only : UDUnits2_StrError
#else
#ifdef with_udunits1
    use UDUnits1, only : UDUnits1_StrError
#endif
#endif

    ! --- in/out ---------------------------------
    
    character(len=256)          ::  UDUnits_StrError
    integer, intent(inout)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_StrError'
    
    ! --- begin ----------------------------------

    ! switch:
    select case ( status )
      !
      case ( UDUNITS_NO_DEFS )
        error_message = 'Procedure in UDUnits module requires native udunits library; define either with_udunits1 or with_udunits2 macro'
      !
      case ( UDUNITS_NOT_INITIALIZED )
        error_message = 'UDUnits module not initialized yet'
      !
      case ( UDUNITS_ALREADY_INITIALIZED )
        error_message = 'UDUnits module already initialized'
      !
      case default
#ifdef with_udunits2
        error_message = trim(UDUnits2_StrError(status))
#else
#ifdef with_udunits1
        error_message = trim(UDUnits1_StrError(status))
#else
        status=UDUNITS_NO_DEFS; return
#endif
#endif
    end select
    
    ! copy to function result:
    UDUnits_StrError = trim(error_message)
    
  end function UDUnits_StrError
  

  ! ====================================================================
  ! ===
  ! === high level routines
  ! ===
  ! ====================================================================
  
  
  subroutine UDUnits_Standard( spec_from, spec_to, status )

#ifdef with_udunits2
    use UDUnits2, only : UDUnits2_Standard, UDUNITS2_NOERR
#else
#ifdef with_udunits1
    use UDUnits1, only : UDUnits1_Standard, UDUNITS1_NOERR
#endif
#endif
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(out)     ::  spec_to
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_Standard'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( .not. initialized ) then; status=UDUNITS_NOT_INITIALIZED; return; end if
    
    ! dummy assignment to avoid errors on un-filled output ...
    spec_to = trim(spec_from)

#ifdef with_udunits2
    call UDUnits2_Standard( spec_from, spec_to, status )
    if (status/=UDUNITS2_NOERR) return
#else
#ifdef with_udunits1
    call UDUnits1_Standard( spec_from, spec_to, status )
    if (status/=UDUNITS1_NOERR) return
#else
    status=UDUNITS_NO_DEFS; return
#endif
#endif
    
    ! ok
    status = 0
    
  end subroutine UDUnits_Standard
  
  
  ! ***
  
  
  subroutine UDUnits_ConversionFactor( spec_from, spec_to, factor, status )

#ifdef with_udunits2
    use UDUnits2, only : UDUnits2_ConversionFactor, UDUNITS2_NOERR
#else
#ifdef with_udunits1
    use UDUnits1, only : UDUnits1_ConversionFactor, UDUNITS1_NOERR
#endif
#endif
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  spec_from
    character(len=*), intent(in)      ::  spec_to
    real, intent(out)                 ::  factor
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/UDUnits_ConversionFactor'
    
    ! --- local ----------------------------------
    
    real(8)   ::  factor_r8

    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. initialized ) then; status=UDUNITS_NOT_INITIALIZED; return; end if

    ! specifications are the same ?
    if ( trim(spec_from) == trim(spec_to) ) then
    
      ! no conversion needed:
      factor = 1.0
    
    ! adhoc ..
    else if ( ( (trim(spec_from) == '(0 - 1)'                    ) .and. (trim(spec_to) == '1'            ) ) .or. &
              ( (trim(spec_from) == 'degrees'                    ) .and. (trim(spec_to) == 'degrees_east' ) ) .or. &
              ( (trim(spec_from) == 'degrees'                    ) .and. (trim(spec_to) == 'degrees_north') ) .or. &
              ( (trim(spec_from) == 'kg kg**-1'                  ) .and. (trim(spec_to) == 'kg/kg'        ) ) .or. &
              ( (trim(spec_from) == 'm s**-1'                    ) .and. (trim(spec_to) == 'm/s'          ) ) .or. &
              ( (trim(spec_from) == 'm**3 m**-3'                 ) .and. (trim(spec_to) == 'm3/m3'        ) ) .or. &
              ( (trim(spec_from) == 'J m**-2 s**-1'              ) .and. (trim(spec_to) == 'J/m2/s'       ) ) .or. &
              ( (trim(spec_from) == 'm of water equivalent'      ) .and. (trim(spec_to) == 'm'            ) ) .or. &
              ( (trim(spec_from) == 'm of water equivalent s**-1') .and. (trim(spec_to) == 'm/s'          ) )      ) then

      ! no conversion needed:
      factor = 1.0

    else

      ! call routine from version specific module:
#ifdef with_udunits1
      call UDUnits1_ConversionFactor( spec_from, spec_to, factor_r8, status )
      if (status/=UDUNITS1_NOERR) return
#else
#ifdef with_udunits2
      call UDUnits2_ConversionFactor( spec_from, spec_to, factor_r8, status )
      if (status/=UDUNITS2_NOERR) return
#else
      status=UDUNITS_NO_DEFS; return
#endif
#endif
      ! cast to parameter kind:
      factor = real(factor_r8)

    end if  ! different specifications
    
    ! ok
    status = 0
    
  end subroutine UDUnits_ConversionFactor
  
  
end module UDUnits


!! ######################################################################
!! ###
!! ### test
!! ###
!! ######################################################################
!!
!!#define IF_UDUNITS_NOT_OK_STOP if (status/=UDUNITS_NOERR) then; write (*,'(a)') trim(UDUnits_StrError(status)); TRACEBACK; stop; end if
!!
!!! f90 -o test.x udunits.F90 \
!!     -I${UDUNITS_HOME}/include \
!!     -L${UDUNITS_HOME}/lib -ludunits
!!
!
!program Test_UDUnits
!
!  use UDUnits
!
!  implicit none
!
!  character(len=*), parameter ::  rname = 'Test_UDUnits'
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
!  call UDUnits_Init( status )
!  IF_UDUNITS_NOT_OK_STOP
!
!  ! * high levell routines
!
!  spec = 'kg/s'
!  call UDUnits_Standard( spec, spec2, status )
!  IF_UDUNITS_NOT_OK_STOP
!  write (*,'(a,": standard name of `",a,"` is `",a,"`")') &
!                      rname, trim(spec), trim(spec2)
!
!  spec = 'gram/cm3' ; spec2 = 'kg/m3'
!  call UDUnits_ConversionFactor( spec, spec2, slope, status )
!  IF_UDUNITS_NOT_OK_STOP
!  write (*,'(a,": conversion factor from `",a,"` to `",a,"` is ",f12.4)') &
!                      rname, trim(spec), trim(spec2), slope
!
!  ! * done with module
!
!  call UDUnits_Done( status )
!  IF_UDUNITS_NOT_OK_STOP
!
!  ! *
!
!  ! info ...
!  write (*,'(a,": ")') rname
!  write (*,'(a,": ** end")') rname
!  write (*,'(a,": ")') rname
!
!end program Test_UDUnits
