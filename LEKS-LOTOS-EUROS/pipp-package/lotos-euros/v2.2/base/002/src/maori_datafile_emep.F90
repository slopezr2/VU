!#################################################################
!
! NAME
!
!   file_emep  -  read daily PM10 observations
!
! DATA FILE
!
!   Example of data file:
!
!    ---[Stations_DailyTimeSeries_PM10_2006.csv]--------------------------
!    compound,date_time,AT02,AT05,CH01,CH02,CH03,CH04,CH05,...
!    PM10,1-1-2006,34.14,10.05,-9999.99,3.5,5.8,1.46,1.2,...
!    :
!    ---------------------------------------------------------------------
!
! USAGE
!
!   use MAORI_DataFile_EMEP
!
!   type(T_MAORI_DataFile_EMEP)        ::  emep_file
!
!   ! open file:
!   call Init( emep_file, 'Stations_DailyTimeSeries_PM10_2006.csv', status )
!   if (status/=0) stop
!
!   ! read next record:
!   call ReadRecord( emep, status )
!   if (status/=0) stop
!
!   ! .. or search record for a certain time interval:
!   call FindRecord( emep, t1, t2, status )
!   if (status/=0) stop
!
!   ! extract record stuff:
!   !   nobs           : integer; number of observations
!   !   t              : TDate ; time in UCT
!   !   nodata         : real ; value used to describe no data
!   call Get( emep, status, nobs=nobs, t, nodata=nodata )
!   if (status/=0) stop
!
!   ! return observation data from current record
!   call GetObservation( emep, iobs, status, &
!                                 station_code, comp, unit, value )
!  
!   ! close:
!   call Done( emep_file, status )
!   if (status/=0) stop
!  
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#################################################################

module MAORI_DataFile_EMEP

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_DataFile_EMEP
  public  ::  Init, Done
  public  ::  FindRecord
  public  ::  ReadRecord
  public  ::  Get
  public  ::  GetObservation
  public  ::  SearchStation

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'emep_data_file'
  
  ! value for no data ...
  real, parameter   ::  emep_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_DataFile_EMEP
    ! file name:
    character(len=512)          ::  fname
    ! seperation character:
    character(len=1)            ::  sep
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
    ! measured component:
    character(len=8)            ::  compname
    character(len=8)            ::  compunit
    ! number of observation:
    integer                     ::  nobs
    ! info per observation:
    character(len=4), pointer   ::  station_code(:)
    ! measured values:
    real, pointer               ::  value(:)
    logical                     ::  filled
    type(TDate)                 ::  t
  end type T_MAORI_DataFile_EMEP


  ! --- interfaces -------------------------
  
  interface Init
    module procedure emep_Init
  end interface
  
  interface Done
    module procedure emep_Done
  end interface
    
  interface ReadRecord
    module procedure emep_ReadRecord
  end interface
    
  interface FindRecord
    module procedure emep_FindRecord
  end interface
  
  interface Get
    module procedure emep_Get
  end interface
  
  interface GetObservation
    module procedure emep_GetObservation
  end interface

  interface SearchStation
    module procedure emep_SearchStation
  end interface


contains


  ! ======================================================================


  subroutine emep_Init( emep, fname, status )
  
    use GO, only : goGetFU
    use GO, only : goReadFromLine, goSplitString
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(out)    ::  emep
    character(len=*), intent(in)          ::  fname
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=4000)   ::  line
    character(len=16)     ::  key
    !integer               ::  iobs
    integer               ::  i
    
    ! --- begin ----------------------------------
    
    ! seperation character:
    emep%sep = ','
    
    ! store file name:
    emep%fname = trim(fname)
    
    ! check ...
    inquire( file=trim(emep%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("EMEP data file not found:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( emep%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( emep%fu, file=trim(emep%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EMEP data file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no observations counted yet:
    emep%nobs = 0
    
    ! 
    ! compound,date_time,AT02,AT05,CH01,CH02,CH03,CH04,CH05,...
    ! PM10,1-1-2006,34.14,10.05,-9999.99,3.5,5.8,1.46,1.2,...
    ! :
    ! 

    ! loop over header lines:
    emep%iline = 0
    do
    
      ! next line ...
      emep%iline = emep%iline + 1
      
      ! read line:
      read (emep%fu,'(a)',iostat=status) line
      if ( status/=0 ) then
        write (gol,'("reading line from EMEP data file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(emep%fname); call goErr
        write (gol,'("  line   : ",i6)') emep%iline; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! skip blank lines:
      if ( len_trim(line) == 0 ) cycle
      
      ! read column header:
      call goReadFromLine( line, key, status, sep=emep%sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( trim(key) /= 'compound' ) then
        write (gol,'("first column should be `compound`, not `",a,"`")') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
      ! read column header:
      call goReadFromLine( line, key, status, sep=emep%sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( trim(key) /= 'date_time' ) then
        write (gol,'("first column should be `date_time`, not `",a,"`")') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
      ! number of fields is equal to number of seperation characters plus 1;
      !   AT02,AT05,CH01,CH02,CH03,CH04,CH05
      emep%nobs = 1
      do i = 1, len_trim(line)
        if ( line(i:i) == emep%sep ) emep%nobs = emep%nobs + 1
      end do
      ! storage:
      allocate( emep%station_code(emep%nobs) )
      allocate( emep%value       (emep%nobs) )
      ! fill station names:
      call goSplitString( line, emep%nobs, emep%station_code, status, sep=emep%sep )
      IF_NOTOK_RETURN(status=1)
      ! leave:
      exit

    end do   ! lines
        
    ! no data values read yet:
    emep%filled = .false.
       
    ! ok
    status = 0
    
  end subroutine emep_Init
  
  
  ! ***
  
  
  subroutine emep_Done( emep, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(inout)       ::  emep
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( emep%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EMEP data file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! clear:
    deallocate( emep%station_code )
    deallocate( emep%value        )
   
    ! ok
    status = 0
    
  end subroutine emep_Done
  

  ! ***
  
  
  subroutine emep_ReadRecord( emep, status )

    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine

    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(inout)       ::  emep
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=4000)   ::  line
    character(len=13)     ::  stime
    integer               ::  day, month, year
    integer               ::  irec
    
    ! --- begin ----------------------------------
    
    ! next line ...
    emep%iline = emep%iline + 1
      
    ! read line:
    read (emep%fu,'(a)',iostat=status) line
    if (status<0) then
      status=-1; return   ! eof
    end if
    if ( status/=0 ) then
      write (gol,'("reading line from EMEP data file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(emep%fname); call goErr
      write (gol,'("  line   : ",i6)') emep%iline; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! 
    ! compound,date_time,AT02,AT05,CH01,CH02,CH03,CH04,CH05,...
    ! PM10,1-1-2006,34.14,10.05,-9999.99,3.5,5.8,1.46,1.2,...
    ! :
    ! 

    ! extract component name:
    call goReadFromLine( line, emep%compname, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! set expected unit ...
    select case ( trim(emep%compname) )
      case ( 'PM10' ) ; emep%compunit = 'ug/m3'
      case default
        write (gol,'("could not guess unit for component `",a,"` ...")') trim(emep%compname); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! extract time field:
    call goReadFromLine( line, stime, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)
    ! split into fields:
    call goReadFromLine( stime,   day, status, sep='-' )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( stime, month, status, sep='-' )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( stime,  year, status, sep='-' )
    IF_NOTOK_RETURN(status=1)
    ! store time:
    emep%t = NewDate( year=year, month=month, day=day )
    
    ! read values:
    do irec = 1, emep%nobs
      ! extract value:
      call goReadFromLine( line, emep%value(irec), status, sep=emep%sep )
      IF_NOTOK_RETURN(status=1)
    end do

    ! record read now:
    emep%filled = .true.

    ! ok
    status = 0
    
  end subroutine emep_ReadRecord
  

  ! ***
  

  ! read records until time is valid

  subroutine emep_FindRecord( emep, t1, t2, status )

    use GO, only : TDate, NewDate, IncrDate, Get
    use GO, only : operator(+), operator(>), operator(==), wrtgol

    ! --- in/out ----------------------------------

    type(T_MAORI_DataFile_EMEP), intent(inout)     ::  emep
    type(TDate), intent(in)                  ::  t1, t2
    integer, intent(out)                     ::  status
 
    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/emep_FindRecord'

    ! --- local ------------------------------------
    
    integer         ::  year, month, day
    type(TDate)     ::  tday

    ! --- begin ------------------------------------
    
    ! daily records, set time of day:
    call Get( t1, year=year, month=month, day=day )
    tday = NewDate( year=year, month=month, day=day )
    
    ! check ...
    if ( t2 > tday+IncrDate(day=1) ) then
      call wrtgol( 'time interval should be within a day : ', t1, ' - ', t2 ); call goErr
      TRACEBACK; status=1; return
    end if

    ! loop until time is found:
    do
      ! something filled ?
      if ( emep%filled ) then
        ! match ? then leave:
        if ( emep%t == tday ) exit
        ! after requested time ? problem ...
        if ( emep%t > tday ) then
          write (gol,'("time of current record exceeds requested time:")'); call goErr
          call wrtgol( '  record t         : ', emep%t ); call goErr
          call wrtgol( '  requested tday   : ', tday ); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! read first or next record:
      call ReadRecord( emep, status )
      if ( status < 0 ) then
        write (gol,'("eof reached but requested time not found ...")'); call goErr
        write (gol,'("  file          : ",a)') trim(emep%fname); call goErr
        call wrtgol( '  requested t   : ', tday ); call goErr
        TRACEBACK; status=1; return
      else if ( status > 0 ) then
        TRACEBACK; status=1; return
      end if
    end do
    
    ! ok
    status = 0
    
  end subroutine emep_FindRecord
          

  ! ***
  
  
  subroutine emep_Get( emep, status, nobs, t, nodata )
  
    use GO, only : TDate
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(inout)  ::  emep
    integer, intent(out)                  ::  status
    integer, intent(out), optional        ::  nobs
    type(TDate), intent(out), optional    ::  t
    real, intent(out), optional           ::  nodata

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! return number of observations ?
    if ( present(nobs) ) nobs = emep%nobs
    
    ! return value used for no data ?
    if ( present(nodata) ) nodata = emep_nodata
    
    ! time of value:
    if ( present(t) ) then
      ! check ...
      if ( .not. emep%filled ) then
        write (gol,'("no values read yet ...")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      t = emep%t
    end if
   
    ! ok
    status = 0
    
  end subroutine emep_Get
  

  ! ***
  
  
  subroutine emep_GetObservation( emep, iobs, status, &
                                 station_code, comp, unit, value )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(inout)        ::  emep
    integer, intent(in)                         ::  iobs
    integer, intent(out)                        ::  status
    
    character(len=*), intent(out), optional     ::  station_code
    character(len=*), intent(out), optional     ::  comp
    character(len=*), intent(out), optional     ::  unit
    real, intent(out), optional                 ::  value

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_GetObservation'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (iobs <= 0) .or. (iobs > emep%nobs) ) then
      write (gol,'("strange observation number:")'); call goErr
      write (gol,'("  iobs   : ",i6)') iobs; call goErr
      write (gol,'("  nobs   : ",i6)') emep%nobs; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return observation property:
    if ( present(station_code) ) station_code = emep%station_code(iobs)
    if ( present(comp        ) ) comp         = emep%compname
    if ( present(unit        ) ) unit         = emep%compunit
    
    ! return observerd value:
    if ( present(value) ) then
      ! check ...
      if ( .not. emep%filled ) then
        write (gol,'("no values read yet ...")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      value = emep%value(iobs)
    end if
   
    ! ok
    status = 0
    
  end subroutine emep_GetObservation
  

  ! ***
  
  
  subroutine emep_SearchStation( emep, station_code, iobs, status )
  
    use GO, only : goMatchValue
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_DataFile_EMEP), intent(inout)        ::  emep
    character(len=*), intent(in)                ::  station_code
    integer, intent(out)                        ::  iobs
    integer, intent(out)                        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_SearchStation'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! search for station name in list, return index:
    call goMatchValue( trim(station_code), emep%station_code, iobs, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine emep_SearchStation
  


end module MAORI_DataFile_EMEP

