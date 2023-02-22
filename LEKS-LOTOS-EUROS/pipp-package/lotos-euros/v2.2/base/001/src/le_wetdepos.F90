!###############################################################################
!
! Wet deposition.
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

module LE_WetDepos

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private
  
  public  ::  LE_WetDepos_Init, LE_WetDepos_Done
  public  ::  LE_WetDepos_Apply


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'LE_WetDepos'

  ! ---
  
  character(len=16)    ::  wet_depos_type
  

contains


  ! ====================================================================


  subroutine LE_WetDepos_Init( rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use LE_WetDepos_Emep, only : LE_WetDepos_Emep_Init
    use LE_WetDepos_Camx, only : LE_WetDepos_Camx_Init


    ! --- in/out ------------------------------
    
    type(TrcFile), intent(in)             ::  rcF
    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! data file for wet scavenging parameters
    call ReadRc( rcF, 'le.wet_depos.type', wet_depos_type, status )
    IF_NOTOK_RETURN(status=1)
    
    ! which wet deposition type
    select case (trim( wet_depos_type ) )
      
      !~ emep
      case ('EMEP' )
        
        call LE_WetDepos_Emep_Init( status )
        IF_NOTOK_RETURN(status=1)
      
      !~ CAMx  
      case ( 'CAMx' )
        
        call LE_WetDepos_Camx_Init( rcF, status )
        IF_NOTOK_RETURN(status=1)
      
      case default
        write (gol,'("Wrong wet deposition type chosen, must be EMEP or CAMx : ",a)') trim(wet_depos_type); call goErr
        TRACEBACK; status=1; return

    end select    
    
    ! ok
    status = 0

  end subroutine LE_WetDepos_Init


  ! ***


  subroutine LE_WetDepos_Done( status )

    use LE_WetDepos_Emep, only : LE_WetDepos_Emep_Done
    use LE_WetDepos_Camx, only : LE_WetDepos_Camx_Done

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! which wet deposition type
    select case (trim( wet_depos_type ) )
      
      !~ emep
      case ('EMEP' )
        
        call LE_WetDepos_Emep_Done( status )
        IF_NOTOK_RETURN(status=1)
      
      !~ CAMx  
      case ( 'CAMx' )
        
        call LE_WetDepos_Camx_Done( status )
        IF_NOTOK_RETURN(status=1)
      
      case default
        write (gol,'("Wrong wet deposition type chosen, must be EMEP or CAMx : ",a)') trim(wet_depos_type); call goErr
        TRACEBACK; status=1; return

    end select    
    
    ! ok
    status = 0
    
  end subroutine LE_WetDepos_Done


  ! ***


  subroutine LE_WetDepos_Apply( c, cwet, update_cwet, deltat, i_step, status )
    
    use LE_WetDepos_Emep, only : LE_WetDepos_Emep_Apply
    use LE_WetDepos_Camx, only : LE_WetDepos_Camx_Apply
    use dims,     only  : nx, ny, nz, nspec
    
    ! --- in/out ------------------------------
  
    real, intent(inout)   ::  c(nx,ny,nz,nspec)
    real, intent(inout)   ::  cwet(nx,ny,nz,nspec)
    logical, intent(in)   ::  update_cwet
    real, intent(in)      ::  deltat    
    integer, intent(in)   ::  i_step
    integer, intent(out)  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Apply'

    ! which wet deposition type
    select case (trim( wet_depos_type ) )
      
      !~ emep
      case ('EMEP' )
        
        call LE_WetDepos_Emep_Apply( c, cwet, update_cwet, deltat, status )
        IF_NOTOK_RETURN(status=1)
      
      !~ CAMx  
      case ( 'CAMx' )
        
        call LE_WetDepos_Camx_Apply( c, cwet, update_cwet, deltat, i_step, status )
        IF_NOTOK_RETURN(status=1)
      
      case default
        write (gol,'("Wrong wet deposition type chosen, must be EMEP or CAMx : ",a)') trim(wet_depos_type); call goErr
        TRACEBACK; status=1; return

    end select    
    
    ! ok
    status = 0

  end subroutine LE_WetDepos_Apply


end module LE_WetDepos
