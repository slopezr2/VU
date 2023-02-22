!#######################################################################
!
! NAME
!
!   MAORI_LocList_EARLINET  -  read list of EARLINET station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    ---[stations.txt]---------------------------------------------
!    #
!    # Locations of EARLINET stations.
!    #
!    # Extracted from : EARLINET_stations.doc
!    #
!    # Records:
!    #   name
!    #   country code
!    #   station code
!    #   latitude north/south
!    #   longitude east/west
!    #   altitude (above sea level)
!    #   Raman equipped station ?
!    #
!    Aberystwith           , UK, ab, 52.4    N,  4.06   W,           , raman
!    Athens                , GR, at, 37.9716 N, 23.7875 E,           , raman
!    Barcelona             , ES, ba, 41.393  N,  2.120  E,           ,
!     :
!    ---------------------------------------------------------------------
!
! USAGE
!
!   use MAORI_LocList_EARLINET_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_EARLINET)        ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'earlinet-loclist.txt', nrec, status )
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

module MAORI_LocList_EARLINET

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_EARLINET
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_EARLINET'
  
  ! value for no data ...
  real, parameter   ::  earlinet_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_EARLINET
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
  end type T_MAORI_LocList_EARLINET


  ! --- interfaces -------------------------
  
  interface Init
    module procedure earlinet_Init
  end interface
  
  interface Done
    module procedure earlinet_Done
  end interface

  interface ReadRecord
    module procedure earlinet_ReadRecord
  end interface  
  


contains


  ! ======================================================================


  subroutine earlinet_Init( earlinet, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EARLINET), intent(out)     ::  earlinet
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/earlinet_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    earlinet%fname = trim(fname)
    
    ! set special characers:
    earlinet%comment = '#'
    earlinet%sep     = ','
    
    ! check ...
    inquire( file=trim(earlinet%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("EARLINET loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(earlinet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( earlinet%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( earlinet%fu, file=trim(earlinet%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EARLINET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(earlinet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    earlinet%iline = 0
    
    ! count records:
    nrec = 0
    do
      ! next line:
      earlinet%iline = earlinet%iline + 1
      ! read line:
      read (earlinet%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from EARLINET loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(earlinet%fname); call goErr
        write (gol,'("  line   : ",i6)') earlinet%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == earlinet%comment ) cycle
      ! increase counter:
      nrec = nrec + 1
    end do
    
    ! close file:
    close( earlinet%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EARLINET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(earlinet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( earlinet%fu, file=trim(earlinet%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EARLINET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(earlinet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    earlinet%iline = 0

    ! ok
    status = 0
    
  end subroutine earlinet_Init
  
  
  ! ***
  
  
  subroutine earlinet_Done( earlinet, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EARLINET), intent(inout)    ::  earlinet
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/earlinet_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( earlinet%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EARLINET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(earlinet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine earlinet_Done
  

  ! ***
  
  
  subroutine earlinet_ReadRecord( earlinet, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EARLINET), intent(inout)      ::  earlinet
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/earlinet_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=128)    ::  line
    character(len=128)    ::  field
    
    ! --- begin ----------------------------------
    
    ! next line ...
    earlinet%iline = earlinet%iline + 1

    ! loop until first non-empty line or end-of-file:
    do
      ! read line:
      read (earlinet%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from EARLINET loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(earlinet%fname); call goErr
        write (gol,'("  line   : ",i6)') earlinet%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == earlinet%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! extract station name:
    call goReadFromLine( line, name, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! extract country code:
    call goReadFromLine( line, field, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    ! add to meta data:
    write (meta,'("country=",a)') trim(field)
    
    ! extract station code:
    call goReadFromLine( line, code, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! extract latitude:
    call goReadFromLine( line, field, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    ! first part is number:
    call goReadFromLine( field, lat, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! remainder is sign:
    select case ( trim(field) )
      case ( 'N', 'n' )
      case ( 'S', 's' ) ; lat = -1.0 * lat
      case default
        write (gol,'("unsupported latitude direction: `",a,"`")') trim(field); call goErr
        write (gol,'("  file    : ",a)') trim(earlinet%fname); call goErr
        write (gol,'("  record  : ",i6)') earlinet%iline; call goErr
        TRACEBACK; status=1; return
    end select

    ! extract longitude:
    call goReadFromLine( line, field, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    ! first part is number:
    call goReadFromLine( field, lon, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! remainder is sign:
    select case ( trim(field) )
      case ( 'E', 'e' )
      case ( 'W', 'w' ) ; lon = -1.0 * lon
      case default
        write (gol,'("unsupported longitude direction: `",a,"`")') trim(field); call goErr
        TRACEBACK; status=1; return
    end select

    ! extract altitude:
    call goReadFromLine( line, field, status, sep=earlinet%sep )
    IF_NOTOK_RETURN(status=1)
    ! empty ?
    if ( len_trim(field) == 0 ) then
      ! set default:
      alt = 0.0
    else
      ! first part is number:
      call goReadFromLine( field, alt, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! remainder is unit:
      select case ( trim(field) )
        case ( 'm asl' )
        case default
          write (gol,'("unsupported altitude unit: `",a,"`")') trim(field); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! add rest to meta data:
    meta = trim(meta)//'; '//trim(line)

    ! ok
    status = 0
    
  end subroutine earlinet_ReadRecord
  

end module MAORI_LocList_EARLINET

