!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Timeprof_Hourly

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Timeprof_Hourly
  public  ::  LE_Emis_Timeprof_Hourly_Init, LE_Emis_Timeprof_Hourly_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Timeprof_Hourly'
  
  
  ! --- types --------------------------------
  
  ! storage for emission time profile:
  type T_Emis_Timeprof_Hourly
    
    ! dimensions:
    integer                   ::  year
    integer                   ::  nhour  ! (8760 for a non leap year)
    ! factor assigned to hour:
    real, allocatable         ::  factor(:)  ! (hour)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Timeprof_Hourly_Init( tprof_hourly, indir, fname, status )

    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var
    use MDF, only : MDF_NETCDF, MDF_READ
    
    ! --- in/out ------------------------------
    
    type(T_Emis_Timeprof_Hourly), intent(out) ::  tprof_hourly
    character(len=*), intent(in)              ::  indir
    character(len=*), intent(in)              ::  fname
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Timeprof_Hourly_Init'
    
    ! --- local -------------------------------
    
    character(len=512)             ::  filename        
    logical                        ::  exist    
    integer                        ::  hid
    integer                        ::  dimid
    integer                        ::  varid
    
    ! --- begin -------------------------------
    
    ! store:
    !tprof_hourly%year         = year
    
    !complete file name
    write( filename, '(a,"/",a)' ) trim(indir), trim(fname)             
    
    ! file should be present:
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! open file:
    call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, hid, status )
    IF_NOTOK_RETURN(status=1)

    ! number of hours in time profile file
    call MDF_Inq_DimID( hid, 'hour', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=tprof_hourly%nhour )
    IF_NOTOK_RETURN(status=1)                    
    
    ! Init:
    allocate( tprof_hourly%factor(tprof_hourly%nhour), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! read profile
    call MDF_Inq_VarID( hid, 'factors', varid, status )
    IF_NOTOK_RETURN(status=1)           
    call MDF_Get_Var( hid, varid, tprof_hourly%factor, status )          
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( abs(sum(tprof_hourly%factor(:))/tprof_hourly%nhour -1 )  >= 1e-2)  then
      write(gol, '("Time profile do not have average 1 ")'); call goErr            
      write(gol, '("Filename: ", a)') trim(filename); call goErr
      write(gol, '("Average = ", f10.5)' ) sum(tprof_hourly%factor(:))/tprof_hourly%nhour   ; call goErr
      TRACEBACK; status=1; return
    end if 

    ! close file:
    call MDF_close( hid, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Timeprof_Hourly_Init
  
  
  ! ***
  

  subroutine LE_Emis_Timeprof_Hourly_Done( tprof_hourly, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Timeprof_Hourly), intent(inout) ::  tprof_hourly
    integer, intent(out)                        ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate(tprof_hourly%factor, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Timeprof_Hourly_Done


end module LE_Emis_Timeprof_Hourly

