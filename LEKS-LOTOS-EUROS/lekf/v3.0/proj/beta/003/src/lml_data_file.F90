!#################################################################
!
! NAME
!
!   file_lml  -  read observation file from lml network
!
! DATA FILE
!
!   Example of data file:
!
!    ---[o3ned12m2006o.spr]---------------------------------------------
!
!    Component: O3               
!
!     unit        ;ug/m3     ;ug/m3     ;...
!     latitude (y);        51.12;        51.54;...
!     longitude(x);         6.04;         5.85;...
!     station numb;          107;          131;...
!     Height/type ;observation;observation;...
!     Origin      ;RILPLUS   ;RILPLUS   ;...
!    Hour\Station ;Posterholt-Vlodropperweg                ;Vredepeel-Vredeweg                      ;...
!    -------------;     ----------;     ----------;...
!    01/01/2006 01; .4815E+02; .5583E+02; ...
!    01/01/2006 02; .5089E+02; .5480E+02; ...
!    01/01/2006 03; .5189E+02; .5378E+02; ...
!    :
!    ---------------------------------------------------------------------
!
!   NOTE: time in CET, not UTC
!
! USAGE
!
!   use file_lml
!
!   type(T_LML_Data_File)        ::  lml_file
!
!   ! open file:
!   call Init( lml_file, 'o3ned12m2006o.spr', status )
!   if (status/=0) stop
!
!   ! read next record:
!   call ReadRecord( lml, status )
!   if (status/=0) stop
!
!   ! extract record stuff:
!   !   nobs           : integer; number of observations
!   !   t_cet, t_uct   : TDate ; time in CET (native) or UCT
!   !   nodata         : real ; value used to describe no data
!   call Get( lml, status, nobs=n, t_cet=t, t_utc=t, nodata=nodata )
!   if (status/=0) stop
!
!   ! return observation data from current record
!   call GetObservation( lml, iobs, status, &
!                                 station, unit, longitude, latitude, &
!                                 station_numb, value )
!  
!   ! close:
!   call Done( lml_file, status )
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
!#################################################################

module lml_data_file

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_LML_Data_File
  public  ::  Init, Done
  public  ::  FindRecord
  public  ::  ReadRecord
  public  ::  Get
  public  ::  GetObservation
  public  ::  SearchStation

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'lml_data_file'
  
  ! value for no data ...
  real, parameter   ::  lml_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_LML_Data_File
    ! file name:
    character(len=512)          ::  fname
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
    ! measured component:
    character(len=8)            ::  compname
    ! number of observation:
    integer                     ::  nobs
    ! info per observation:
    character(len=40), pointer  ::  station(:)
    character(len=6), pointer   ::  unit(:)
    real, pointer               ::  longitude(:)
    real, pointer               ::  latitude(:)
    integer, pointer            ::  station_numb(:)
    ! measured values:
    real, pointer               ::  value(:)
    logical                     ::  filled
    type(TDate)                 ::  t_cet, t_utc
  end type T_LML_Data_File


  ! --- interfaces -------------------------
  
  interface Init
    module procedure lml_Init
  end interface
  
  interface Done
    module procedure lml_Done
  end interface
    
  interface ReadRecord
    module procedure lml_ReadRecord
  end interface
    
  interface FindRecord
    module procedure lml_FindRecord
  end interface
  
  interface Get
    module procedure lml_Get
  end interface
  
  interface GetObservation
    module procedure lml_GetObservation
  end interface

  interface SearchStation
    module procedure lml_SearchStation
  end interface


contains


  ! ======================================================================


  subroutine lml_Init( lml, fname, status )
  
    use GO, only : goGetFU
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(out)    ::  lml
    character(len=*), intent(in)          ::  fname
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Init'
    
    ! required headers:
    integer, parameter    ::  ihead_unit    = 1
    integer, parameter    ::  ihead_lat     = 2
    integer, parameter    ::  ihead_lon     = 3
    integer, parameter    ::  ihead_numb    = 4
    integer, parameter    ::  ihead_station = 5
    integer, parameter    ::  nhead = ihead_station
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=4000)   ::  line
    character(len=16)     ::  key
    integer               ::  iobs
    logical               ::  heads(nhead)
    
    ! --- begin ----------------------------------
    
    ! store file name:
    lml%fname = trim(fname)
    
    ! check ...
    inquire( file=trim(lml%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("LML data file not found:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( lml%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( lml%fu, file=trim(lml%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening LML data file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no observations counted yet:
    lml%nobs = 0
    
    ! no header values found yet:
    heads = .false.
    ! loop over header lines:
    lml%iline = 0
    do
    
      ! next line ...
      lml%iline = lml%iline + 1
      
      ! read line:
      read (lml%fu,'(a)',iostat=status) line
      if ( status/=0 ) then
        write (gol,'("reading line from LML data file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
        write (gol,'("  line   : ",i6)') lml%iline; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! skip blank lines:
      if ( len_trim(line) == 0 ) cycle
    
      ! end of header ? then leave:
      if ( line(1:10) == '----------' ) exit
      
      ! component line ?
      if ( line(1:10) == 'Component:' ) then
        ! extract first part:
        call goReadFromLine( line, key, status, sep=':' )
        IF_NOTOK_RETURN(status=1)
        ! remainder is component name:
        lml%compname = trim(line)
        ! next line ...
        cycle
      end if
        
      ! count observations if not done yet:
      if ( lml%nobs == 0 ) then
        ! number of fields is equal to number of ';' characters:
        !   unit        ;ug/m3     ;ug/m3     ;...
        do iobs = 1, len_trim(line)
          if ( line(iobs:iobs) == ';' ) lml%nobs = lml%nobs + 1
        end do
        ! check ...
        if ( lml%nobs == 0 ) then
          write (gol,'("could not detect number of observations from header line:")'); call goErr
          write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
          write (gol,'("  line   : ",i6)') lml%iline; call goErr
          TRACEBACK; status=1; return
        end if
        ! storage:
        allocate( lml%unit        (lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lml%longitude   (lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lml%latitude    (lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lml%station_numb(lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lml%station     (lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lml%value       (lml%nobs), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
      
      ! extract key:
      call goReadFromLine( line, key, status, sep=';' )
      IF_NOTOK_RETURN(status=1)

      ! extract values:
      select case ( trim(key) )
        case ( 'unit' )
          do iobs = 1, lml%nobs
            call goReadFromLine( line, lml%unit(iobs), status, sep=';' )
            IF_NOTOK_RETURN(status=1)
          end do
          heads(ihead_unit) = .true.
        case ( 'latitude (y)' )
          do iobs = 1, lml%nobs
            call goReadFromLine( line, lml%latitude(iobs), status, sep=';' )
            IF_NOTOK_RETURN(status=1)
          end do
          heads(ihead_lat) = .true.
        case ( 'longitude(x)' )
          do iobs = 1, lml%nobs
            call goReadFromLine( line, lml%longitude(iobs), status, sep=';' )
            IF_NOTOK_RETURN(status=1)
          end do
          heads(ihead_lon) = .true.
        case ( 'station numb' )
          do iobs = 1, lml%nobs
            call goReadFromLine( line, lml%station_numb(iobs), status, sep=';' )
            IF_NOTOK_RETURN(status=1)
          end do
          heads(ihead_numb) = .true.
        case ( 'Height/type' )
          ! skip ...
        case ( 'Origin' )
          ! skip ...
        case ( 'Hour\Station' )
          do iobs = 1, lml%nobs
            call goReadFromLine( line, lml%station(iobs), status, sep=';' )
            IF_NOTOK_RETURN(status=1)
          end do
          heads(ihead_station) = .true.
        case default
          !write (gol,'("unsupported key in header line:")'); call goErr
          !write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
          !write (gol,'("  line   : ",i6)') lml%iline; call goErr
          !write (gol,'("  key    : ",a)') trim(key); call goErr
          !TRACEBACK; status=1; return
          ! continue ...
          write (gol,'("WARNING - unsupported LML header line: ",a)') trim(key); call goPr
      end select
      
    end do  ! header lines
    
    ! check ...
    if ( .not. all(heads) ) then
      write (gol,'("not all required header lines found : ")'); call goErr
      write (gol,'("  unit      : ",l2)') heads(ihead_unit   ); call goErr
      write (gol,'("  lat       : ",l2)') heads(ihead_lat    ); call goErr
      write (gol,'("  lon       : ",l2)') heads(ihead_lon    ); call goErr
      write (gol,'("  numb      : ",l2)') heads(ihead_numb   ); call goErr
      write (gol,'("  station   : ",l2)') heads(ihead_station); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no data values read yet:
    lml%filled = .false.
       
    ! ok
    status = 0
    
  end subroutine lml_Init
  
  
  ! ***
  
  
  subroutine lml_Done( lml, status )
  
    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(inout)       ::  lml
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( lml%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing LML data file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! clear:
    if ( lml%nobs > 0 ) then
      deallocate( lml%unit         )
      deallocate( lml%longitude    )
      deallocate( lml%latitude     )
      deallocate( lml%station_numb )
      deallocate( lml%station      )
      deallocate( lml%value        )
    end if
   
    ! ok
    status = 0
    
  end subroutine lml_Done
  

  ! ***
  
  
  subroutine lml_ReadRecord( lml, status )

    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine

    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(inout)       ::  lml
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=4000)   ::  line
    character(len=13)     ::  stime
    character(len=10)     ::  sval
    integer               ::  day, month, year, hour
    integer               ::  iobs
    logical               ::  warned_negative
    
    ! --- begin ----------------------------------
    
    ! next line ...
    lml%iline = lml%iline + 1
      
    ! read line:
    read (lml%fu,'(a)',iostat=status) line
    if (status<0) then
      status=-1; return   ! eof
    end if
    if ( status/=0 ) then
      write (gol,'("reading line from LML data file:")'); call goErr
      write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
      write (gol,'("  line   : ",i6)') lml%iline; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract time field:
    !   01/01/2006 01
    call goReadFromLine( line, stime, status, sep=';' )
    IF_NOTOK_RETURN(status=1)
    ! read values:
    read (stime,'(i2,x,i2,x,i4,x,i2)',iostat=status) day, month, year, hour
    if ( status/=0 ) then
      write (gol,'("reading time from LML record:")'); call goErr
      write (gol,'("  file         : ",a)') trim(lml%fname); call goErr
      write (gol,'("  line number  : ",i6)') lml%iline; call goErr
      write (gol,'("  time string  : ",a)') trim(stime); call goErr
      write (gol,'("  time values  : ",4i5)') day, month, year, hour; call goErr
      TRACEBACK; status=1; return
    end if
    ! store:
    lml%t_cet = NewDate( year=year, month=month, day=day, hour=hour )
    lml%t_utc = lml%t_cet - IncrDate(hour=1)
    
    ! no warnings yet:
    warned_negative = .false.
    
    ! extract values:
    do iobs = 1, lml%nobs
      ! read field into character variable:
      call goReadFromLine( line, sval, status, sep=';' )
      IF_NOTOK_RETURN(status=1)
      ! empty ?
      if ( len_trim(sval) == 0 ) then
        ! no data ...
        lml%value(iobs) = lml_nodata
      else
        ! read from character variable:
        read (sval,*,iostat=status) lml%value(iobs)
        if (status/=0) then
          write (gol,'("reading real value from character variable:")'); call goErr
          write (gol,'("  sval  : `",a,"`")') trim(sval); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ..
        if ( lml%value(iobs) < 0.0 ) then
          ! warning ?
          if ( .not. warned_negative ) then
            !write (gol,'("WARNING - negative values in line ",i6," of ",a)') &
            !         lml%iline, trim(lml%fname); call goErr
            warned_negative = .true.
          end if
          ! set to no-data:
          lml%value(iobs) = lml_nodata
        end if
      end if
    end do
    
    ! record read now:
    lml%filled = .true.

    ! ok
    status = 0
    
  end subroutine lml_ReadRecord
  

  ! ***
  

  ! read records until time is valid

  subroutine lml_FindRecord( lml, t_utc, status )

    use GO, only : TDate, operator(>), operator(==), wrtgol

    ! --- in/out ----------------------------------

    type(T_LML_Data_File), intent(inout)     ::  lml
    type(TDate), intent(in)                  ::  t_utc
    integer, intent(out)                     ::  status
 
    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/lml_FindRecord'

    ! --- local ------------------------------------

    ! --- begin ------------------------------------

    ! loop until time is found:
    do
      ! something filled ?
      if ( lml%filled ) then
        ! match ? then leave:
        if ( lml%t_utc == t_utc ) exit
        ! after requested time ? problem ...
        if ( lml%t_utc > t_utc ) then
          write (gol,'("time of current record exceeds requested time:")'); call goErr
          call wrtgol( '  record t_utc      : ', lml%t_utc ); call goErr
          call wrtgol( '  requested t_utc   : ', t_utc ); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! read first or next record:
      call ReadRecord( lml, status )
      if ( status < 0 ) then
        write (gol,'("eof reached but requested time not found ...")'); call goErr
        write (gol,'("  file              : ",a)') trim(lml%fname); call goErr
        call wrtgol( '  requested t_utc   : ', t_utc ); call goErr
        TRACEBACK; status=1; return
      else if ( status > 0 ) then
        TRACEBACK; status=1; return
      end if
    end do
    
    ! ok
    status = 0
    
  end subroutine lml_FindRecord
          

  ! ***
  
  
  subroutine lml_Get( lml, status, nobs, t_cet, t_utc, nodata )
  
    use GO, only : TDate
    
    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(inout)  ::  lml
    integer, intent(out)                  ::  status
    integer, intent(out), optional        ::  nobs
    type(TDate), intent(out), optional    ::  t_cet
    type(TDate), intent(out), optional    ::  t_utc
    real, intent(out), optional           ::  nodata

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! return number of observations ?
    if ( present(nobs) ) nobs = lml%nobs
    
    ! return value used for no data ?
    if ( present(nodata) ) nodata = lml_nodata
    
    ! time of values:
    if ( present(t_cet) .or. present(t_utc) ) then
      ! check ...
      if ( .not. lml%filled ) then
        write (gol,'("no values read yet ...")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      if ( present(t_cet) ) t_cet = lml%t_cet
      if ( present(t_utc) ) t_utc = lml%t_utc
    end if
   
    ! ok
    status = 0
    
  end subroutine lml_Get
  

  ! ***
  
  
  subroutine lml_GetObservation( lml, iobs, status, &
                                 station, comp, unit, longitude, latitude, &
                                 station_numb, value )
  
    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(inout)        ::  lml
    integer, intent(in)                         ::  iobs
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  station
    character(len=*), intent(out), optional     ::  comp
    character(len=*), intent(out), optional     ::  unit
    real, intent(out), optional                 ::  longitude
    real, intent(out), optional                 ::  latitude
    integer, intent(out), optional              ::  station_numb
    real, intent(out), optional                 ::  value

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_GetObservation'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (iobs <= 0) .or. (iobs > lml%nobs) ) then
      write (gol,'("strange observation number:")'); call goErr
      write (gol,'("  iobs   : ",i6)') iobs; call goErr
      write (gol,'("  nobs   : ",i6)') lml%nobs; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return observation property:
    if ( present(station     ) ) station      = lml%station     (iobs)
    if ( present(comp        ) ) comp         = lml%compname
    if ( present(unit        ) ) unit         = lml%unit        (iobs)
    if ( present(longitude   ) ) longitude    = lml%longitude   (iobs)
    if ( present(latitude    ) ) latitude     = lml%latitude    (iobs)
    if ( present(station_numb) ) station_numb = lml%station_numb(iobs)
    
    ! return observerd value:
    if ( present(value) ) then
      ! check ...
      if ( .not. lml%filled ) then
        write (gol,'("no values read yet ...")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      value = lml%value(iobs)
    end if
   
    ! ok
    status = 0
    
  end subroutine lml_GetObservation
  

  ! ***
  
  
  subroutine lml_SearchStation( lml, station, iobs, status )
  
    use GO, only : goMatchValue
    
    ! --- in/out --------------------------------
    
    type(T_LML_Data_File), intent(inout)        ::  lml
    character(len=*), intent(in)                ::  station
    integer, intent(out)                        ::  iobs
    integer, intent(out)                        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_SearchStation'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! search for station name in list, return index:
    call goMatchValue( trim(station), lml%station, iobs, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine lml_SearchStation
  


end module lml_data_file

