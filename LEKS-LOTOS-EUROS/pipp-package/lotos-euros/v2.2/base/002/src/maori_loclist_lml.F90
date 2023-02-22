!#######################################################################
!
! NAME
!
!   MAORI_LocList_LML  -  read list of LML station locations and meta data
!
! DATA FILE
!
!   Example of loclist file:
!
!    ---[stations.txt]---------------------------------------------
!    22
!    lml_107   6.04306   51.12028   0   'Posterholt-Vlodropperweg     ' 
!    lml_131   5.85361   51.54111   0   'Vredepeel-Vredeweg           '
!    ...
!
!    ---------------------------------------------------------------------
!
!
! USAGE
!
!   use MAORI_LocList_LML_file
!
!   integer                         ::  status
!   type(T_MAORI_LocList_LML)        ::  file
!   integer                         ::  nrec
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! open file, count records:
!   call Init( file, 'lml-loclist.txt', nrec, status )
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

module MAORI_LocList_LML

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_LML
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_LML'
  
  ! value for no data ...
  real, parameter   ::  lml_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_LML
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
  end type T_MAORI_LocList_LML


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
  


contains


  ! ======================================================================


  subroutine lml_Init( lml, fname, nrec, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_LML), intent(out)     ::  lml
    character(len=*), intent(in)         ::  fname
    integer, intent(out)                 ::  nrec
    integer, intent(out)                 ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=512)    ::  line
    
    ! --- begin ----------------------------------
    
    ! store file name:
    lml%fname = trim(fname)
    
    ! set special characers:
    lml%comment = '#'
    lml%sep     = ' '
    
    ! check ...
    inquire( file=trim(lml%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("LML loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( lml%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( lml%fu, file=trim(lml%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening LML loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    lml%iline = 0
    
    ! read no of records:
    read (lml%fu,*,iostat=status) nrec
    if ( status/=0 ) then
      write (gol,'("reading nr of records from LML loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! reset counter:
    lml%iline = lml%iline + 1

    ! ok
    status = 0
    
  end subroutine lml_Init
  
  
  ! ***
  
  
  subroutine lml_Done( lml, status )
  
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_LML), intent(inout)    ::  lml
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( lml%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing LML loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine lml_Done
  

  ! ***
  
  
  subroutine lml_ReadRecord( lml, name, code, lon, lat, alt, meta, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_MAORI_LocList_LML), intent(inout)      ::  lml
    character(len=*), intent(out)            ::  name
    character(len=*), intent(out)            ::  code
    real, intent(out)                        ::  lon, lat   ! [degree]
    real, intent(out)                        ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)            ::  meta
    integer, intent(out)                     ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=128)    ::  line
    
    ! --- begin ----------------------------------
    
    ! next line ...
    lml%iline = lml%iline + 1

    ! loop until first non-empty line or end-of-file:
    do
      ! read line:
      read (lml%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from LML loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
        write (gol,'("  line   : ",i6)') lml%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == lml%comment ) cycle
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do

    ! extract station code:
    call goReadFromLine( line, code, status, sep=lml%sep )
    IF_NOTOK_RETURN(status=1)
    
    ! extract longitude:
    call goReadFromLine( line, lon, status, sep=lml%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract latitude:
    call goReadFromLine( line, lat, status, sep=lml%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract altitude:
    call goReadFromLine( line, alt, status, sep=lml%sep )
    IF_NOTOK_RETURN(status=1)

    ! extract station name:
    call goReadFromLine( line, name, status, sep=lml%sep )
    IF_NOTOK_RETURN(status=1)
        
    ! no meta data:
    meta = ''

    ! ok
    status = 0
    
  end subroutine lml_ReadRecord
  

end module MAORI_LocList_LML

