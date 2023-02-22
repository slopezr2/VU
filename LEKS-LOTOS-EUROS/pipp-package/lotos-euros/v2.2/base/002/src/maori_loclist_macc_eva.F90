!#######################################################################
!
! NAME
!
!   MAORI_LocList_MACC_EVA  -  read list of MACC_EVA station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    ---[sO3_Background_rural_loclist_assimilation-set1_v6]---------------
!    AT0ENK1;13.671114;48.391670;525;Background;rural;Enzenkirchen im Sauwald
!    AT0ILL1;16.766403;47.770000;117;Background;rural;Illmitz
!    ...
!
!    ---------------------------------------------------------------------
!
!
! USAGE
!
!   use MAORI_LocList_loclist_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_MACC_EVA)  ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'MACC_EVA-loclist.txt', nrec, status )
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

module MAORI_LocList_MACC_EVA

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_MACC_EVA
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_MACC_EVA'
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_MACC_EVA
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
  end type T_MAORI_LocList_MACC_EVA


  ! --- interfaces -------------------------
  
  interface Init
    module procedure loclist_Init
  end interface
  
  interface Done
    module procedure loclist_Done
  end interface

  interface ReadRecord
    module procedure loclist_ReadRecord
  end interface 



contains


  ! ======================================================================


  subroutine loclist_Init( loclist, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_MACC_EVA), intent(out)     ::  loclist
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/loclist_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    loclist%fname = trim(fname)
    
    ! set special characers:
    loclist%comment = '#'
    loclist%sep     = ';'
    
    ! check ...
    inquire( file=trim(loclist%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("MACC_EVA loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(loclist%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( loclist%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( loclist%fu, file=trim(loclist%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening MACC_EVA loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(loclist%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    loclist%iline = 0
    
    ! count records:
    nrec = 0
    do
      ! next line:
      loclist%iline = loclist%iline + 1
      ! read line:
      read (loclist%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(loclist%fname); call goErr
        write (gol,'("  line   : ",i6)') loclist%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == loclist%comment ) cycle
      ! increase counter:
      nrec = nrec + 1
    end do
    
    ! close file:
    close( loclist%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(loclist%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( loclist%fu, file=trim(loclist%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(loclist%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    loclist%iline = 0

    ! ok
    status = 0
    
  end subroutine loclist_Init
  
  
  ! ***
  
  
  subroutine loclist_Done( loclist, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_MACC_EVA), intent(inout)    ::  loclist
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/loclist_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( loclist%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing MACC_EVA loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(loclist%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine loclist_Done
  

  ! ***
  
  
  subroutine loclist_ReadRecord( loclist, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_MACC_EVA), intent(inout)      ::  loclist
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/loclist_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=128)    ::  line
    character(len=32)     ::  cat, subcat
    
    ! --- begin ----------------------------------
    
    ! next line ...
    loclist%iline = loclist%iline + 1

    ! loop until first non-empty line or end-of-file:
    do
      ! read line:
      read (loclist%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from MACC_EVA loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(loclist%fname); call goErr
        write (gol,'("  line   : ",i6)') loclist%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == loclist%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do

    ! extract station code:
    call goReadFromLine( line, code, status, sep=loclist%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! extract longitude:
    call goReadFromLine( line, lon, status, sep=loclist%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract latitude:
    call goReadFromLine( line, lat, status, sep=loclist%sep )
    IF_NOTOK_RETURN(status=1)

    !! extract altitude:
    !call goReadFromLine( line, alt, status, sep=loclist%sep )
    !IF_NOTOK_RETURN(status=1)
    ! dummy ....
    alt = 0.0

    ! extract category:
    call goReadFromLine( line, cat, status, sep=loclist%sep )
    IF_NOTOK_RETURN(status=1)
    ! extract subcategory:
    call goReadFromLine( line, subcat, status, sep=loclist%sep )
    IF_NOTOK_RETURN(status=1)

    !! extract station name:
    !call goReadFromLine( line, name, status, sep=loclist%sep )
    !IF_NOTOK_RETURN(status=1)
    ! dummy ...
    name = trim(code)
        
    ! no meta data:
    meta = ''

    ! ok
    status = 0
    
  end subroutine loclist_ReadRecord
  

end module MAORI_LocList_MACC_EVA

