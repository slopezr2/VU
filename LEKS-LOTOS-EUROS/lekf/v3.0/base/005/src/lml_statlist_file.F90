!#################################################################
!
! NAME
!
!   lml_statlist_file  -  read list of LML station names and meta data
!
! DATA FILE
!
!   Example of statlist file:
!
!    ---[lml-statlist.txt]---------------------------------------------
!    Posterholt-Vlodropperweg                
!    Vredepeel-Vredeweg                      
!    Wijnandsrade-Opfergeltstraat            
!     :
!    ---------------------------------------------------------------------
!
! USAGE
!
!   use lml_statlist_file
!
!   type(T_LML_Statlist_File)        ::  file
!
!   call Init( file, 'lml-statlist.txt', status )
!   if (status/=0) stop
!
!   call ReadRecord( file, station_name, status )
!   if (status/=0) stop
!
!   call Done( file, status )
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

module lml_statlist_file

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_LML_Statlist_File
  public  ::  Init, Done
  public  ::  ReadRecord
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'lml_statlist_file'
  
  ! value for no data ...
  real, parameter   ::  lml_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_LML_Statlist_File
    ! file name:
    character(len=512)          ::  fname
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
  end type T_LML_Statlist_File


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


  subroutine lml_Init( lml, fname, status )
  
    use GO, only : goGetFU
    
    ! --- in/out --------------------------------
    
    type(T_LML_Statlist_File), intent(out)         ::  lml
    character(len=*), intent(in)          ::  fname
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    
    ! --- begin ----------------------------------
    
    ! store file name:
    lml%fname = trim(fname)
    
    ! check ...
    inquire( file=trim(lml%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("LML statlist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( lml%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( lml%fu, file=trim(lml%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening LML statlist file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no line read yet:
    lml%iline = 0

    ! ok
    status = 0
    
  end subroutine lml_Init
  
  
  ! ***
  
  
  subroutine lml_Done( lml, status )
  
    ! --- in/out --------------------------------
    
    type(T_LML_Statlist_File), intent(inout)       ::  lml
    integer, intent(out)                  ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/lml_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( lml%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing LML statlist file:")'); call goErr
      write (gol,'("  ",a)') trim(lml%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine lml_Done
  

  ! ***
  
  
  subroutine lml_ReadRecord( lml, station, status )
  
    use GO, only : TDate, NewDate, IncrDate, operator(-)
    use GO, only : goReadFromLine
    
    ! --- in/out --------------------------------
    
    type(T_LML_Statlist_File), intent(inout)       ::  lml
    character(len=*), intent(out)         ::  station
    integer, intent(out)                  ::  status

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
        write (gol,'("reading line from LML statlist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(lml%fname); call goErr
        write (gol,'("  line   : ",i6)') lml%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! extract station name:
    call goReadFromLine( line, station, status, sep=';' )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine lml_ReadRecord
  

end module lml_statlist_file

