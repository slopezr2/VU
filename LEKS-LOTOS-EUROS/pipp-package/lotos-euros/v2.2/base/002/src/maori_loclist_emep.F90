!#######################################################################
!
! NAME
!
!   MAORI_LocList_EMEP  -  read list of EMEP station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    ---[stations.txt]---------------------------------------------
!    AT02   47.77   16.77  117.00 Illmitz
!    AT03   47.55   11.72  960.00 Achenkirch
!    AT04   47.65   13.20  851.00 St.-Koloman
!     :
!    ---------------------------------------------------------------------
!
! USAGE
!
!   use MAORI_LocList_EMEP_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_EMEP)            ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'emep-loclist.txt', nrec, status )
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

module MAORI_LocList_EMEP

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_EMEP
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_EMEP'
  
  ! value for no data ...
  real, parameter   ::  emep_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_EMEP
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
  end type T_MAORI_LocList_EMEP


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
  


contains


  ! ======================================================================


  subroutine emep_Init( emep, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EMEP), intent(out)     ::  emep
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    emep%fname = trim(fname)
    
    ! set special characers:
    emep%comment = '#'
    emep%sep     = ' '
    
    ! check ...
    inquire( file=trim(emep%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("EMEP loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( emep%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( emep%fu, file=trim(emep%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    emep%iline = 0
    
    ! count records:
    nrec = 0
    do
      ! next line:
      emep%iline = emep%iline + 1
      ! read line:
      read (emep%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
      if ( status/=0 ) then
        write (gol,'("reading line from EMEP loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(emep%fname); call goErr
        write (gol,'("  line   : ",i6)') emep%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == emep%comment ) cycle
      ! increase counter:
      nrec = nrec + 1
    end do
    
    ! close file:
    close( emep%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( emep%fu, file=trim(emep%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    emep%iline = 0

    ! ok
    status = 0
    
  end subroutine emep_Init
  
  
  ! ***
  
  
  subroutine emep_Done( emep, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EMEP), intent(inout)    ::  emep
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( emep%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing EMEP loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(emep%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine emep_Done
  

  ! ***
  
  
  subroutine emep_ReadRecord( emep, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_EMEP), intent(inout)      ::  emep
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/emep_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=128)    ::  line
    
    ! --- begin ----------------------------------
    
    ! next line ...
    emep%iline = emep%iline + 1

    ! loop until first non-empty line or end-of-file:
    do
      ! read line:
      read (emep%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from EMEP loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(emep%fname); call goErr
        write (gol,'("  line   : ",i6)') emep%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == emep%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! extract station code:
    call goReadFromLine( line, code, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! extract latitude:
    call goReadFromLine( line, lat, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract longitude:
    call goReadFromLine( line, lon, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract altitude:
    call goReadFromLine( line, alt, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract station name:
    call goReadFromLine( line, name, status, sep=emep%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! no meta data:
    meta = ''
    
    ! ok
    status = 0
    
  end subroutine emep_ReadRecord
  

end module MAORI_LocList_EMEP

