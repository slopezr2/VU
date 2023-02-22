!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! Tools to convert gribfiles to NetCDF using GRIB API.
!
!######################################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; stop 1; end if
!
!######################################################################

program C3PO_Grib_to_NetCDF

  implicit none
  
  ! --- const --------------------------------------------------
  
  character(len=*), parameter  ::  rname = 'C3PO_Grib_to_NetCDF'
  
  
  ! --- local ---------------------------------------------------
  
  integer                         ::  status
  integer                         ::  narg, iarg
  integer                         ::  arglength
  character(len=:), allocatable   ::  argvalue
  character(len=:), allocatable   ::  gribfile
  character(len=:), allocatable   ::  ncfile
  
  ! --- begin ---------------------------------------------------
  
  ! info ..
  write (*,'("")')
  write (*,'("C3PO - Grib to NetCDF convertor")')
  write (*,'("")')
  
  !
  ! arguments
  !
  
  ! count:
  narg = Command_Argument_Count()

  ! loop:
  do iarg = 1, narg
    ! get length:
    call Get_Command_Argument( iarg, length=arglength, status=status )
    if ( status /= 0 ) then
      write (*,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength
      TRACEBACK; stop 1
    end if
    ! storage:
    allocate( character(len=arglength) :: argvalue )
    ! get value: 
    call Get_Command_Argument( iarg, value=argvalue, status=status )
    if ( status /= 0 ) then
      write (*,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength
      TRACEBACK; stop 1
    end if
    ! flag ?
    if ( argvalue(1:1) == '-' ) then
      ! switch:
      select case ( trim(argvalue) )
        case ( '-h' )
          write (*,'("Usage:")')
          write (*,'("  c3po_grib_to_netcdf.x <gribfile> <ncfile>")')
          write (*,'("  c3po_grib_to_netcdf.x -h")')
        case default
          write (*,'("unsupported argument `",a,"`")') trim(argvalue)
          TRACEBACK; stop 1
      end select
    else if ( .not. allocated(gribfile) ) then
      allocate( character(len=arglength) :: gribfile )
      gribfile = trim(argvalue)
    else if ( .not. allocated(ncfile) ) then
      allocate( character(len=arglength) :: ncfile )
      ncfile = trim(argvalue)
    else
      write (*,'("unsupported argument `",a,"`")') trim(argvalue)
      TRACEBACK; stop 1
    end if
    ! clear:
    deallocate( argvalue )
  end do   ! arguments
  
  ! check ...
  if ( .not. allocated(gribfile) ) then
    write (*,'("no gribfile argument specified")')
    TRACEBACK; stop 1
  end if
  ! check ...
  if ( .not. allocated(ncfile) ) then
    write (*,'("no ncfile argument specified")')
    TRACEBACK; stop 1
  end if
      
  ! info ...
  write (*,'("grib file : ",a)') trim(gribfile)
  write (*,'("nc   file : ",a)') trim(ncfile)

  
  !
  ! end
  !
  
  ! info ...
  write (*,'("")')
  write (*,'("End")')
  write (*,'("")')

end program C3PO_Grib_to_NetCDF

