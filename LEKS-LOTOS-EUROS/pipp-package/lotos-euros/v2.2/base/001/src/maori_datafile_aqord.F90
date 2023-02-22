!#######################################################################
!
! NAME
!
!   MAORI_DataFile_AQORD
!
! DESCRIPTION
!
!   Read ground obs from nc files (converted from text).
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_DataFile_AQORD

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  MAORI_DataFile_AQORD_ReadRecord


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_DataFile_AQORD'


  ! --- types ------------------------------


contains


  ! ====================================================================


  subroutine MAORI_DataFile_AQORD_ReadRecord( fname, tr, value, unit, status, concname )
  
    use GO, only : TDate, TIncrDate, Get, wrtgol
    use GO, only : Compare_Date_Num
    use GO, only : operator(-), rTotal
    use MDF
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)          ::  fname
    type(TDate), intent(in)               ::  tr(2)
    real, intent(out)                     ::  value
    character(len=*), intent(out)         ::  unit
    integer, intent(out)                  ::  status
    character(len=*), optional            ::  concname

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_DataFile_AQORD_ReadRecord'

    
    ! --- local ----------------------------------
    
    logical               ::  exist
    integer               ::  hid, dimid, varid
    integer               ::  nrec
    integer               ::  irec
    integer               ::  i
    real                  ::  conc(1)
    character(len=64)     ::  units
    character(len=64)     ::  calendar
    type(TDate)           ::  tref
    type(TIncrDate)       ::  tstep
    real, allocatable     ::  time1(:)
    real                  ::  time2(1)
    real                  ::  time_1, time_2
    
    ! --- begin ----------------------------------
    
    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! open file:
    call MDF_Open( trim(fname), MDF_NETCDF, MDF_READ, hid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of time records:
    call MDF_Inq_DimID( hid, 'time', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=nrec )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( time1(nrec) )
    
    ! start times:
    call MDF_Inq_VarID( hid, 'time1', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, time1, status )
    IF_NOTOK_RETURN(status=1)
    
    ! units:
    call MDF_Get_Att( hid, varid, 'units', units, status )
    IF_NOTOK_RETURN(status=1)
    !! calendar:
    !call MDF_Get_Att( hid, varid, 'calendar', calendar, status )
    !IF_NOTOK_RETURN(status=1)
    calendar = 'gregorian'
    
    ! required start:
    !call Get( tr(1), time6=time6_1 )
    ! search ..
    irec = -1
    do i = 1, nrec
      ! compare numeric value with start of timerange;
      ! return status -1 if not matched, 0 if ok:
      call Compare_Date_Num( tr(1), time1(i), units, calendar, status )
      IF_ERROR_RETURN(status=1)
      ! match?
      if ( status == 0 ) then
        irec = i
        exit
      end if
    end do
    ! check ..
    if ( irec < 0 ) then
      call wrtgol( 'start time not found : ', tr(1) ); call goErr
      write (gol,'("  file                   : ",a)') trim(fname); call goErr
      call wrtgol( '  requested interval     : ', tr ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! clear:
    deallocate( time1, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! end time:
    call MDF_Inq_VarID( hid, 'time2', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, time2, status, start=(/irec/), count=(/1/) )
    IF_NOTOK_RETURN(status=1)
    ! required end time:
    ! compare numeric value with start of timerange;
    ! return status -1 if not matched, 0 if ok:
    call Compare_Date_Num( tr(2), time2(1), units, calendar, status )
    IF_ERROR_RETURN(status=1)
    ! no match?
    if ( status < 0 ) then
      write (gol,'("end times do not match:")'); call goErr
      write (gol,'("  file                : ",a)') trim(fname); call goErr
      write (gol,'("  record              : ",i6)') irec; call goErr
      write (gol,'("  record end time     : ",f12.4)') time2; call goErr
      call wrtgol( '  requested           : ', tr ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! concentration:
    if ( present(concname) ) then
      call MDF_Inq_VarID( hid, trim(concname), varid, status )
      IF_NOTOK_RETURN(status=1)
    else
      call MDF_Inq_VarID( hid, 'conc', varid, status )
      IF_NOTOK_RETURN(status=1)
    end if
    call MDF_Get_Var( hid, varid, conc, status, start=(/irec/), count=(/1/) )
    IF_NOTOK_RETURN(status=1)
    
    ! store output:
    value = conc(1)
    
    ! get unit:
    if ( present(concname) ) then
      call MDF_Get_Att( hid, varid, 'units',  unit, status )
      IF_NOTOK_RETURN(status=1)
    else
      call MDF_Get_Att( hid, MDF_GLOBAL, 'unit',  unit, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! close:
    call MDF_Close( hid, status )
    IF_NOTOK_RETURN(status=1)

    
    ! ok
    status = 0
    
  end subroutine MAORI_DataFile_AQORD_ReadRecord


end module MAORI_DataFile_AQORD
