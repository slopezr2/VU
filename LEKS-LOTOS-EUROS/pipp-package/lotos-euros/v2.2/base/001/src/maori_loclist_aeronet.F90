!#######################################################################
!
! NAME
!
!   MAORI_LocList_AERONET  -  read list of AERONET station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    ---[stations.txt]---------------------------------------------
!    #
!    # Locations of AERONET stations.
!    #
!    # Extracted from : AERONET_stations.doc
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
!  The AERONET sites measure AOT at different wavelengts;
!  the model simulates at 550nm :
!
!        wavelengths [nm]
!        AERONET  LOTOS-EUROS
!         340
!         380
!         440
!                  450  (commented)
!         500
!         532
!         535
!                  550
!                  650  (commented)
!         670
!         870
!        1020
!        1640
!
! USAGE
!
!   use MAORI_LocList_AERONET_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_AERONET)        ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'aeronet-loclist.txt', nrec, status )
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

module MAORI_LocList_AERONET

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_AERONET
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_AERONET'
  
  ! value for no data ...
  real, parameter   ::  aeronet_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_AERONET
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
  end type T_MAORI_LocList_AERONET


  ! --- interfaces -------------------------
  
  interface Init
    module procedure aeronet_Init
  end interface
  
  interface Done
    module procedure aeronet_Done
  end interface

  interface ReadRecord
    module procedure aeronet_ReadRecord
  end interface  
  


contains


  ! ======================================================================


  subroutine aeronet_Init( aeronet, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AERONET), intent(out)     ::  aeronet
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aeronet_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    aeronet%fname = trim(fname)
    
    ! set special characers:
    aeronet%comment = '#'
    aeronet%sep     = ' '
    
    ! check ...
    inquire( file=trim(aeronet%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("AERONET loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(aeronet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( aeronet%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( aeronet%fu, file=trim(aeronet%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AERONET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aeronet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    aeronet%iline = 0
    
    ! count records:
    nrec = 0
    do
      ! next line:
      aeronet%iline = aeronet%iline + 1
      ! read line:
      read (aeronet%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from AERONET loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(aeronet%fname); call goErr
        write (gol,'("  line   : ",i6)') aeronet%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == aeronet%comment ) cycle
      ! increase counter:
      nrec = nrec + 1
    end do
    
    ! close file:
    close( aeronet%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AERONET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aeronet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( aeronet%fu, file=trim(aeronet%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening AERONET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aeronet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    aeronet%iline = 0

    ! ok
    status = 0
    
  end subroutine aeronet_Init
  
  
  ! ***
  
  
  subroutine aeronet_Done( aeronet, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AERONET), intent(inout)    ::  aeronet
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aeronet_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( aeronet%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing AERONET loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(aeronet%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine aeronet_Done
  

  ! ***
  
  
  subroutine aeronet_ReadRecord( aeronet, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_AERONET), intent(inout)      ::  aeronet
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/aeronet_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=128)    ::  line
    
    ! --- begin ----------------------------------
    
    ! next line ...
    aeronet%iline = aeronet%iline + 1

    ! loop until first non-empty line or end-of-file:
    do
      ! read line:
      read (aeronet%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from AERONET loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(aeronet%fname); call goErr
        write (gol,'("  line   : ",i6)') aeronet%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == aeronet%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do

    ! extract station name:
    call goReadFromLine( line, name, status, sep=aeronet%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! no station code ..
    code = '-'
    
    ! extract longitude:
    call goReadFromLine( line, lon, status, sep=aeronet%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract latitude:
    call goReadFromLine( line, lat, status, sep=aeronet%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract altitude:
    call goReadFromLine( line, alt, status, sep=aeronet%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! no meta data:
    meta = ''

    ! ok
    status = 0
    
  end subroutine aeronet_ReadRecord
  

end module MAORI_LocList_AERONET

