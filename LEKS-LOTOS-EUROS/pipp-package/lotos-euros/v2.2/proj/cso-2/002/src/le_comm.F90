!###############################################################################
!
! NAME
!   LE_Comm - communication between domains I/O task
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Comm

  use GO, only : gol, goPr, goErr
  use GO, only : TgoComm

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Comm_Init, LE_Comm_Done
  public  ::  goc_dom
  public  ::  goc_io


  ! --- const ------------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Comm'
  

  ! --- var --------------------------------------

  ! i/o tasks:
  type(TgoComm)                     ::  goc_io
  ! domain tasks:
  type(TgoComm)                     ::  goc_dom


contains


  ! ====================================================================
  
  
  subroutine LE_Comm_Init( rcF, status )
  
    use GO     , only : TrcFile
    use GO     , only : goc
      
    ! --- in/out ---------------------------------
    
    integer, intent(out)                      ::  status
    type(TrcFile), intent(in)                 ::  rcF
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Comm_Init'
    
    ! --- local ----------------------------------

    logical                           ::  with_io_proc
    integer                           ::  i

    ! --- begin ----------------------------------
    
    ! use seperate i/o task?
    call rcF%Get( 'le.io.proc', with_io_proc, status )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'(a,": global group:")') rname; call goPr
    write (gol,'(a,":   number of processes       : ",i0)') rname, goc%npes; call goPr
    write (gol,'(a,":   process rank (0-based)    : ",i0)') rname, goc%id; call goPr

    ! ~~ i/o group

    ! serial, or no designated i/o processor?
    if ( (goc%npes == 1) .or. (.not. with_io_proc) ) then
      ! root performs all tasks:
      call goc_io%InitGroup( goc, (/goc%root_id/), status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! i/o tasks are performed by last process, controlled by root:
      call goc_io%InitGroup( goc, (/goc%root_id,goc%npes-1/), status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! info ...
    write (gol,'(a,": i/o group:")') rname; call goPr
    write (gol,'(a,":   included in this process  : ",l1)') rname, goc_io%enabled; call goPr
    if ( goc_io%enabled ) then
      write (gol,'(a,":   number of processes       : ",i0)') rname, goc_io%npes; call goPr
      write (gol,'(a,":   process rank (0-based)    : ",i0)') rname, goc_io%id; call goPr
    end if

    ! ~~ dom group

    ! serial?
    if ( goc%npes == 1 ) then
      ! domain group has just the root:
      call goc_dom%InitGroup( goc, (/goc%root_id/), status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! with i/o task?
      if ( with_io_proc ) then
        ! domain group has all-except-last processors:
        call goc_dom%InitGroup( goc, (/(i,i=0,goc%npes-2)/), status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! domain group has all processors:
        call goc_dom%InitGroup( goc, (/(i,i=0,goc%npes-1)/), status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end if

    ! info ...
    write (gol,'(a,": domains group:")') rname; call goPr
    write (gol,'(a,":   included in this process  : ",l1)') rname, goc_dom%enabled; call goPr
    if ( goc_dom%enabled ) then
      write (gol,'(a,":   number of processes       : ",i0)') rname, goc_dom%npes; call goPr
      write (gol,'(a,":   process rank (0-based)    : ",i0)') rname, goc_dom%id; call goPr
    end if

    ! ok
    status = 0

  end subroutine LE_Comm_Init
  
  ! ***
  
  
  subroutine LE_Comm_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Comm_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    !call goc_io%Done( status )
    !IF_NOT_OK_RETURN(status=1)
    call goc_dom%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Comm_Done



end module LE_Comm


