!###############################################################################
!
! NAME
!   LE_Config - general configurations used in multiple modules
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Config

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Config_Init, LE_Config_Done
  
  public  ::  rcF, rcFP
  public  ::  outputdir


  ! --- const ------------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Config'
  

  ! --- var --------------------------------------

  ! settings:
  type(TRcFile), target     ::  rcF
  ! pointer to settings, passed to routines for usage at run-time:
  type(TrcFile), pointer    ::  rcFP
  
  ! output directory, might include processor number:
  character(len=1024)       ::  outputdir


contains


  ! ====================================================================
  
  
  subroutine LE_Config_Init( status, rcfile )
  
    use GO     , only : TrcFile
    use GO     , only : pathsep
    use GO     , only : goSystem
    use GO     , only : CheckDir
    !use LE_Comm, only : goc_dom
      
    ! --- in/out ---------------------------------
    
    integer, intent(out)                      ::  status
    character(len=*), intent(in), optional    ::  rcfile
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Config_Init'
    
    ! --- local ----------------------------------

    character(len=1024)       ::  rcfile0
    character(len=1024)       ::  basedir
    logical                   ::  exist

    ! --- begin ----------------------------------


    !
    ! ** settings file
    !
    
    ! argument?
    if ( present(rcfile) ) then
      ! copy:
      rcfile0 = rcfile
    else
      ! extract argument, return name of rcfile:
      call LE_Arguments( rcfile0, status )
      if (status/=0) then
        call LE_Print_Usage( status )
        status=1; return
      end if
    end if

    ! info ...
    write (gol,'("Name of LOTOS-EUROS rc file : ",a)') trim(rcfile0); call goPr

    ! open rcfile to read other settings:
    call rcF%Init( rcfile0, status )
    IF_NOTOK_RETURN(status=1)
    
    ! assign pointer:
    rcFP => rcF


    !
    ! ** output directory
    !
    
    ! output directory:
    call rcF%Get( 'le.output.outdir', outputdir, status )
    IF_NOTOK_RETURN(status=1)
    
    !! full path including domain number if not the root domain:
    !if ( goc%id > 0 ) then
    !  write (outputdir,'(a,a,"d",i2.2)') trim(basedir), pathsep, goc%id
    !else
    !  outputdir = trim(basedir)
    !end if
    
    ! create if necessary:
    call CheckDir( trim(outputdir)//pathsep//'.', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Config_Init
  
  ! ***
  
  
  subroutine LE_Config_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Config_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! done with settings:
    call rcF%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Config_Done


  ! ===================================================================
  ! ===
  ! === arguments
  ! ===
  ! ===================================================================


  subroutine LE_Arguments( rcfile, status )

    use GO, only : goArgCount, goGetArg

    ! --- in/out ----------------------------------

    character(len=*), intent(out)     ::  rcfile
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Arguments'

    ! --- local -----------------------------------

    integer               ::  narg
    integer               ::  iarg
    character(len=1024)   ::  line

    ! --- begin -----------------------------------

    ! on root only, since some mpirun version do not parse
    ! all arguments to each executable:

    ! number of arguments:
    call goArgCount( narg, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( narg == 0 ) then
      write (gol,'("no arguments found ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! defaults:
    rcfile = 'None'

    ! loop over arguments:
    iarg = 0
    do
      ! next:
      iarg = iarg + 1
      ! get argument:
      call goGetArg( iarg, line, status )
      IF_NOTOK_RETURN(status=1)
      ! not filled yet ?
      if ( trim(rcfile) == 'None' ) then
        rcfile = trim(line)
      else
        write (gol,'("unsupported argument : ",a)') trim(line); call goErr
        TRACEBACK; status=1; return
      end if
      ! last one is processed now ?
      if ( iarg == narg ) exit
    end do

    ! ok
    status = 0

  end subroutine LE_Arguments


  ! ***


  subroutine LE_Print_Usage( status )

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- begin ----------------------------------

    ! display usage line:
    write (gol,'("Usage: lotos-euros.x <rcfile>")'); call goPr

    ! ok
    status = 0

  end subroutine LE_Print_Usage


end module LE_Config


