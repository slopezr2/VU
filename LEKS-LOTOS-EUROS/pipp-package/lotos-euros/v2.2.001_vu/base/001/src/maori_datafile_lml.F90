!#######################################################################
!
! NAME
!
!   MAORI_DataFile_LML
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

module MAORI_DataFile_LML

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  MAORI_DataFile_LML_ReadRecord


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_DataFile_LML'


  ! --- types ------------------------------


contains


  ! ====================================================================


  subroutine MAORI_DataFile_LML_ReadRecord( fname, tr, stationcode, value, unit, status )
  
    use GO, only : TDate, TIncrDate, Get, wrtgol
!    use GO, only : Extract_Ref_and_Step
!    use GO, only : operator(-), rTotal
    use MDF
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)          ::  fname
    type(TDate), intent(in)               ::  tr(2)
    real, intent(out)                     ::  value
    character(len=*), intent(in)          ::  stationcode
    character(len=*), intent(out)         ::  unit
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_DataFile_LML_ReadRecord'

    
    ! --- local ----------------------------------
    
    logical               ::  exist
    integer               ::  hid, dimid, varid
    integer               ::  nrec
    integer               ::  irec
    integer               ::  i
!    integer, allocatable  ::  start_time(:,:)
!    integer               ::  end_time(6)
    real                  ::  conc(1)
!    integer               ::  time6_1(6)
!    integer               ::  time6_2(6)
    character(len=64)     ::  units
!    character(len=64)     ::  calendar
    type(TDate)           ::  tref
    type(TIncrDate)       ::  tstep
    real, allocatable     ::  start_time(:,:)
!    real, allocatable     ::  time1(:)
!    real                  ::  time2(1)
!    real                  ::  time_1, time_2
    integer              ::  time6(6)
    
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
    
!    ! storage:
    allocate( start_time(6,nrec) )
!    
!    ! start times:
    call MDF_Inq_VarID( hid, 'start_time', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, start_time, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
!    allocate( time1(nrec) )
    
    ! start times:
!    call MDF_Inq_VarID( hid, 'start_time', varid, status )
!    IF_NOTOK_RETURN(status=1)
!    call MDF_Get_Var( hid, varid, time1, status )
!    IF_NOTOK_RETURN(status=1)
    
    ! units:
!    call MDF_Get_Att( hid, varid, 'units', units, status )
!    IF_NOTOK_RETURN(status=1)
    !! calendar:
    !call MDF_Get_Att( hid, varid, 'calendar', calendar, status )
    !IF_NOTOK_RETURN(status=1)
!    calendar = 'gregorian'
    ! extract offset and step from "hours since 2008-01-01 00:00:00"
!    call Extract_Ref_and_Step( units, calendar, tref, tstep, status )
!    IF_NOTOK_RETURN(status=1)
    
    ! required start:
    !call Get( tr(1), time6=time6_1 )
    !time_1 = rTotal( tr(1) - tref, 'sec' ) / rTotal(tstep,'sec')
    call Get( tr(1), time6=time6 )    
! search ..
    irec = -1
    do i = 1, nrec
!      if ( all( start_time(:,i) == time6_1 ) ) then
      if ( all( time6 == start_time(:,i) ) ) then
        irec = i
        exit
      end if
    end do
    ! check ..
    if ( irec < 0 ) then
      call wrtgol( 'start time not found : ', tr(1) ); call goErr
      write (gol,'("  file                   : ",a)') trim(fname); call goErr
      call wrtgol( '  requested interval     : ', tr(1), ' , ', tr(2) ); call goErr
      TRACEBACK; status=1; return
    end if
    
!    ! clear:
    deallocate( start_time )
    
    ! clear:
!    deallocate( time1 )
    
!     ! end time:
! !    call MDF_Inq_VarID( hid, 'end_time', varid, status )
!     call MDF_Inq_VarID( hid, 'time2', varid, status )
!     IF_NOTOK_RETURN(status=1)
! !    call MDF_Get_Var( hid, varid, end_time, start=(/1,irec/), count=(/6,1/), status )
!     call MDF_Get_Var( hid, varid, time2, status, start=(/irec/), count=(/1/) )
!     IF_NOTOK_RETURN(status=1)
!     ! required end time:
! !    call Get( tr(2), time6=time6_2 )
!     time_2 = rTotal( tr(2) - tref, 'sec' ) / rTotal(tstep,'sec')
!     ! check ...
! !    if ( any( end_time /= time6_2 ) ) then
!     if ( abs( time_2 - time2(1) ) > 0.1 ) then
!       write (gol,'("end times do not match:")'); call goErr
!       write (gol,'("  file                : ",a)') trim(fname); call goErr
!       write (gol,'("  record              : ",i6)') irec; call goErr
!       write (gol,'("  record end time     : ",f12.4)') time2; call goErr
!       call wrtgol( '  requested           : ', tr(1), ' , ', tr(2) ); call goErr
! !      write (gol,'("  record end time  : ",i4,5i3)') end_time; call goErr
!       write (gol,'("  requested end time  : ",f12.4)') time_2; call goErr
!       TRACEBACK; status=1; return
!     end if
    
    ! concentration:
    call MDF_Inq_VarID( hid, trim(stationcode)//'_conc', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, conc, status, start=(/irec/), count=(/1/) )
    IF_NOTOK_RETURN(status=1)
    
    ! store output:
    value = conc(1)
    
    ! get unit:
    call MDF_Get_Att( hid, MDF_GLOBAL, 'unit',  unit, status )
    IF_NOTOK_RETURN(status=1)
    
    ! close:
    call MDF_Close( hid, status )
    IF_NOTOK_RETURN(status=1)

    
    ! ok
    status = 0
    
  end subroutine MAORI_DataFile_LML_ReadRecord


end module MAORI_DataFile_LML
