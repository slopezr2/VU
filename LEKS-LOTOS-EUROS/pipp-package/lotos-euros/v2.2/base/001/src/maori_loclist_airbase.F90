!#######################################################################
!
! NAME
!
!   MAORI_LocList_AirBase  -  read list of AirBase station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    oI_code;station name;lat;long;alt;type of station;type of area;characteristics of zone;city;local code;network code;time reference;AirBase code
!    CZ0TNUJ;Navsi u Jablunkova;49.5933;18.7439;380;Background;rural;natural;;TNUJ;CZ001A;UTC;CZ0135A
!    DENW063;Horn-Bad Meinberg Egge;51.8267;8.9506;430;Background;rural;unknown;VELDROM;DENW063;DE004A;UTC;DE0680A
!     :
!
! USAGE
!
!   use MAORI_LocList_AirBase_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_AirBase)            ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'airbase-loclist.txt', nrec, status )
!   if (status/=0) stop
!
!   ! loop over records:
!   do irec = 1, nrec
!     ! read location, fill meta data with station name etc:
!     call ReadRecord( file, lon, lat, alt, meta, status )
!     if (status<0) exit   ! eof
!     if (status/=0) stop
!   end do
!
!   call Done( file, status )
!   if (status/=0) stop
!  
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_LocList_AirBase

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_AirBase
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_AirBase'
  
  ! value for no data ...
  real, parameter   ::  airbase_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_AirBase
    ! file name:
    character(len=512)          ::  fname
    ! comment character:
    character(len=1)            ::  comment
    ! field seperation character:
    character(len=1)            ::  sep
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
    ! header line:
    character(len=1024)         ::  headers
  end type T_MAORI_LocList_AirBase


  ! --- interfaces -------------------------
  
  interface Init
    module procedure airbase_Init
  end interface
  
  interface Done
    module procedure airbase_Done
  end interface

  interface ReadRecord
    module procedure airbase_ReadRecord
  end interface  
  


contains


  ! ======================================================================


  subroutine airbase_Init( airbase, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AirBase), intent(out)     ::  airbase
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/airbase_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    airbase%fname = trim(fname)
    
    ! set special characers:
    airbase%comment = '#'
    airbase%sep     = ';'
    
    ! check ...
    inquire( file=trim(airbase%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("AirBase loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(airbase%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( airbase%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( airbase%fu, file=trim(airbase%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AirBase loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(airbase%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    airbase%iline = 0
    
    ! record counter:
    nrec = 0
    do
      ! next line:
      airbase%iline = airbase%iline + 1
      ! read line:
      read (airbase%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from AirBase loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(airbase%fname); call goErr
        write (gol,'("  line   : ",i6)') airbase%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == airbase%comment ) cycle
      ! increase counter:
      nrec = nrec + 1
    end do
    ! do not count header line:
    nrec = nrec - 1
    
    ! close file:
    close( airbase%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AirBase loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(airbase%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( airbase%fu, file=trim(airbase%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AirBase loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(airbase%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    airbase%iline = 0
    ! loop until header line is read:
    do
      ! next line:
      airbase%iline = airbase%iline + 1
      ! read line:
      read (airbase%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from AirBase loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(airbase%fname); call goErr
        write (gol,'("  line   : ",i6)') airbase%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == airbase%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    ! store:
    airbase%headers = trim(line)

    ! ok
    status = 0
    
  end subroutine airbase_Init
  
  
  ! ***
  
  
  subroutine airbase_Done( airbase, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AirBase), intent(inout)    ::  airbase
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/airbase_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( airbase%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AirBase loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(airbase%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine airbase_Done
  

  ! ***
  
  
  subroutine airbase_ReadRecord( airbase, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AirBase), intent(inout)      ::  airbase
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/airbase_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=1024)    ::  line
    character(len=64)      ::  field
    character(len=1024)    ::  headers
    character(len=128)     ::  header
    
    ! --- begin ----------------------------------
    
    ! loop until first non-empty line or end-of-file:
    do
      ! next line:
      airbase%iline = airbase%iline + 1
      ! read line:
      read (airbase%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from AirBase loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(airbase%fname); call goErr
        write (gol,'("  line   : ",i6)') airbase%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == airbase%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! set dummy values:
    name = ''
    code = ''
    lon  = -999
    lat  = -999
    alt  = -999

    ! copy headers, parts will be chopped:
    headers = trim(airbase%headers)   
    ! extract all parts:
    do
      ! extract header:
      call goReadFromLine( headers, header, status, sep=airbase%sep )
      IF_NOTOK_RETURN(status=1)
      ! switch:
      select case ( trim(header) )
        case ( 'EoI_code' )
          ! extract station code:
          call goReadFromLine( line, code, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'lat' )
          ! extract latitude:
          call goReadFromLine( line, lat, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'long' )
          ! extract longitude:
          call goReadFromLine( line, lon, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'alt' )
          ! extract altitude:
          call goReadFromLine( line, alt, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
          ! trap dummy values ...
          if ( alt < -900 ) alt = 0.0
        case ( 'station name' )
          ! extract station name:
          call goReadFromLine( line, name, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
        case default
          ! just extract this field, don't use it ...
          call goReadFromLine( line, field, status, sep=airbase%sep )
          IF_NOTOK_RETURN(status=1)
      end select
      ! empty ? then leave:
      if ( len_trim(line) == 0 ) exit
    end do
    
    ! check ...
    if ( len_trim(name) == 0 ) then
      write (gol,'("could not find value for `name` in headers : ",a)') trim(airbase%headers); call goErr
      write (gol,'("  file : ",a)') trim(airbase%fname); call goErr
      write (gol,'("  line : ",i6)') airbase%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( len_trim(code) == 0 ) then
      write (gol,'("could not find value for `code` in headers : ",a)') trim(airbase%headers); call goErr
      write (gol,'("  file : ",a)') trim(airbase%fname); call goErr
      write (gol,'("  line : ",i6)') airbase%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( lon < -900 ) then
      write (gol,'("could not find value for `lon` in headers : ",a)') trim(airbase%headers); call goErr
      write (gol,'("  file : ",a)') trim(airbase%fname); call goErr
      write (gol,'("  line : ",i6)') airbase%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( lat < -900 ) then
      write (gol,'("could not find value for `lat` in headers : ",a)') trim(airbase%headers); call goErr
      write (gol,'("  file : ",a)') trim(airbase%fname); call goErr
      write (gol,'("  line : ",i6)') airbase%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( alt < -900 ) then
      write (gol,'("could not find value for `alt` in headers : ",a)') trim(airbase%headers); call goErr
      write (gol,'("  file : ",a)') trim(airbase%fname); call goErr
      write (gol,'("  line : ",i6)') airbase%iline; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no meta data:
    meta = ''
    
    ! ok
    status = 0
    
  end subroutine airbase_ReadRecord
  

end module MAORI_LocList_AirBase

