!#######################################################################
!
! NAME
!
!   MAORI_LocList_AQORD  -  read list of AQORD station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    oI_code;station name;lat;long;alt;type of station;type of area;characteristics of zone;city;local code;network code;time reference;AQORD code
!    CZ0TNUJ;Navsi u Jablunkova;49.5933;18.7439;380;Background;rural;natural;;TNUJ;CZ001A;UTC;CZ0135A
!    DENW063;Horn-Bad Meinberg Egge;51.8267;8.9506;430;Background;rural;unknown;VELDROM;DENW063;DE004A;UTC;DE0680A
!     :
!
! USAGE
!
!   use MAORI_LocList_AQORD_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_AQORD)            ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'aqord-loclist.txt', nrec, status )
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

module MAORI_LocList_AQORD

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_AQORD
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_AQORD'
  
  ! value for no data ...
  real, parameter   ::  aqord_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_AQORD
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
  end type T_MAORI_LocList_AQORD


  ! --- interfaces -------------------------
  
  interface Init
    module procedure aqord_Init
  end interface
  
  interface Done
    module procedure aqord_Done
  end interface

  interface ReadRecord
    module procedure aqord_ReadRecord
  end interface  
  


contains


  ! ======================================================================


  subroutine aqord_Init( aqord, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AQORD), intent(out)     ::  aqord
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aqord_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    aqord%fname = trim(fname)
    
    ! set special characers:
    aqord%comment = '#'
    aqord%sep     = ','
    
    ! check ...
    inquire( file=trim(aqord%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("AQORD loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(aqord%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( aqord%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( aqord%fu, file=trim(aqord%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AQORD loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aqord%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    aqord%iline = 0
    
    ! record counter:
    nrec = 0
    do
      ! next line:
      aqord%iline = aqord%iline + 1
      ! read line:
      read (aqord%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from AQORD loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(aqord%fname); call goErr
        write (gol,'("  line   : ",i6)') aqord%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == aqord%comment ) cycle
      ! no comment (incl header), increase counter:
      nrec = nrec + 1
    end do

    ! NOT needed, header is comment
    !! do not count header line:
    !nrec = nrec - 1
    
    ! close file:
    close( aqord%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AQORD loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aqord%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( aqord%fu, file=trim(aqord%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AQORD loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aqord%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    aqord%iline = 0
    ! loop until header line is read:
    do
      ! next line:
      aqord%iline = aqord%iline + 1
      ! read line:
      read (aqord%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from AQORD loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(aqord%fname); call goErr
        write (gol,'("  line   : ",i6)') aqord%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! headers ?
      if ( line(1:9) == '#headers:' ) then
      	aqord%headers = trim(line)
        exit
      end if
      ! skip comment:
      if ( line(1:1) == aqord%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    !! store:
    !aqord%headers = trim(line)

    ! ok
    status = 0
    
  end subroutine aqord_Init
  
  
  ! ***
  
  
  subroutine aqord_Done( aqord, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AQORD), intent(inout)    ::  aqord
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aqord_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( aqord%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AQORD loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aqord%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine aqord_Done
  

  ! ***
  
  
  subroutine aqord_ReadRecord( aqord, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AQORD), intent(inout)      ::  aqord
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aqord_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=1024)    ::  line
    character(len=64)      ::  field
    character(len=1024)    ::  headers
    character(len=128)     ::  header
    
    ! --- begin ----------------------------------
    
    ! loop until first non-empty line or end-of-file:
    do
      ! next line:
      aqord%iline = aqord%iline + 1
      ! read line:
      read (aqord%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from AQORD loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(aqord%fname); call goErr
        write (gol,'("  line   : ",i6)') aqord%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == aqord%comment ) cycle
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
    headers = trim(aqord%headers)   
    ! extract all parts:
    do
      ! extract header:
      call goReadFromLine( headers, header, status, sep=aqord%sep )
      IF_NOTOK_RETURN(status=1)
      ! switch:
      select case ( trim(header) )
        case ( 'code' )
          ! extract station code:
          call goReadFromLine( line, code, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'lat' )
          ! extract latitude:
          call goReadFromLine( line, lat, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'lon' )
          ! extract longitude:
          call goReadFromLine( line, lon, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
        case ( 'height' )
          ! extract altitude:
          call goReadFromLine( line, alt, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
          ! trap dummy values ...
          if ( alt < -900 ) alt = 0.0
        case ( 'name' )
          ! extract station name:
          call goReadFromLine( line, name, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
        case default
          ! just extract this field, don't use it ...
          call goReadFromLine( line, field, status, sep=aqord%sep )
          IF_NOTOK_RETURN(status=1)
      end select
      ! empty ? then leave:
      if ( len_trim(line) == 0 ) exit
    end do
    
    ! check ...
    if ( len_trim(name) == 0 ) then
      write (gol,'("could not find value for `name` in headers : ",a)') trim(aqord%headers); call goErr
      write (gol,'("  file : ",a)') trim(aqord%fname); call goErr
      write (gol,'("  line : ",i6)') aqord%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( len_trim(code) == 0 ) then
      write (gol,'("could not find value for `code` in headers : ",a)') trim(aqord%headers); call goErr
      write (gol,'("  file : ",a)') trim(aqord%fname); call goErr
      write (gol,'("  line : ",i6)') aqord%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( lon < -900 ) then
      write (gol,'("could not find value for `lon` in headers : ",a)') trim(aqord%headers); call goErr
      write (gol,'("  file : ",a)') trim(aqord%fname); call goErr
      write (gol,'("  line : ",i6)') aqord%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( lat < -900 ) then
      write (gol,'("could not find value for `lat` in headers : ",a)') trim(aqord%headers); call goErr
      write (gol,'("  file : ",a)') trim(aqord%fname); call goErr
      write (gol,'("  line : ",i6)') aqord%iline; call goErr
      TRACEBACK; status=1; return
    end if
    if ( alt < -900 ) then
      write (gol,'("could not find value for `alt` in headers : ",a)') trim(aqord%headers); call goErr
      write (gol,'("  file : ",a)') trim(aqord%fname); call goErr
      write (gol,'("  line : ",i6)') aqord%iline; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no meta data:
    meta = ''
    
    ! ok
    status = 0
    
  end subroutine aqord_ReadRecord
  

end module MAORI_LocList_AQORD

