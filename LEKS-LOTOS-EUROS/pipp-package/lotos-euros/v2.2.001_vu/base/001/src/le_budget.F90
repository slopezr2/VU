!###############################################################################
!
! LE_Budget - budget stuff
!
! NOTE: to have ensemble of budgets, the budget type
!   should become part of the state vector
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

module LE_Budget

  use GO, only : gol, goPr, goErr

  use LE_Budget_DryDepos, only : T_DryDepos_Budget
  use LE_Budget_WetDepos, only : T_WetDepos_Budget

  implicit none


  ! --- in/out --------------------------------
  
  private

  public  ::  T_Budget
  
  public  ::  Budget_Init, Budget_Done
  public  ::  Budget_Update
  public  ::  Budget_Update_O3max
  public  ::  Budget_Reset_Hour
  public  ::  Budget_Reset_Day
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Budget'


  ! --- types --------------------------------------
  
  type T_Budget
    ! collect:
    type(T_DryDepos_Budget)   ::  drydepos
    type(T_WetDepos_Budget)   ::  wetdepos
    ! ozone maximum
    real, allocatable         ::  o3max(:,:)
  end type T_Budget


contains



  ! ========================================================================
  

  subroutine Budget_Init( bud, status )
  
    use Dims              , only : nx, ny
    use Indices           , only : i_SOxa, i_NOya, i_NHxa
    use LE_Budget_Param   , only : nbud, ispec_bud
    use LE_Budget_DryDepos, only : DryDepos_Budget_Init
    use LE_Budget_WetDepos, only : WetDepos_Budget_Init
  
    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(out)             ::  bud
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Budget_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------

    ! set accumulated tracers corresponding to budgets:
    !                        SOx     NOy    NHx
    ispec_bud(1:nbud) = (/ i_SOxa, i_NOya, i_NHxa /)
    
    ! init budgets:
    call DryDepos_Budget_Init( bud%drydepos, status )
    IF_NOTOK_RETURN(status=1)
    call WetDepos_Budget_Init( bud%wetdepos, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for ozone maximum; init to zero:
    allocate( bud%o3max(nx,ny) ) ; bud%o3max = 0.0

    ! ok
    status = 0
    
  end subroutine Budget_Init
  
  
  ! ***
  
  
  subroutine Budget_Done( bud, status )
  
    use LE_Budget_DryDepos, only : DryDepos_Budget_Done
    use LE_Budget_WetDepos, only : WetDepos_Budget_Done
  
    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(inout)           ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Budget_Done'
    
    ! --- begin ----------------------------------
    
    ! done with budgets:
    call DryDepos_Budget_Done( bud%drydepos, status )
    IF_NOTOK_RETURN(status=1)
    call WetDepos_Budget_Done( bud%wetdepos, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( bud%o3max )

    ! ok
    status = 0
    
  end subroutine Budget_Done


  ! ***
  
  
  subroutine Budget_Update( bud, t, status ) 
  
    use GO                , only : TDate
    use LE_Budget_DryDepos, only : DryDepos_Budget_Update
    use LE_Budget_WetDepos, only : WetDepos_Budget_Update

    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(inout)           ::  bud
    type(TDate), intent(in)                 ::  t
    integer, intent(inout)                  ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Budget_Update'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! Update with budgets:
    call DryDepos_Budget_Update( bud%drydepos, t, status )
    IF_NOTOK_RETURN(status=1)
    call WetDepos_Budget_Update( bud%wetdepos, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Budget_Update


  ! ***
  
  
  subroutine Budget_Update_O3max( bud, cg_o3, status ) 
  
    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(inout)           ::  bud
    real, intent(in)                        ::  cg_o3(:,:)  ! ozone at ground level
    integer, intent(inout)                  ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Budget_Update_O3max'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! Update daily ozone maximum.
    ! Not fully correct, in fact le_output should be called first,
    ! but since the daily maximum is somewhere at 16h it does not matter  
    where ( cg_o3 > bud%o3max )
      bud%o3max = cg_o3
    end where

    ! ok
    status = 0
    
  end subroutine Budget_Update_O3max


  ! ***
  
  
  subroutine Budget_Reset_Hour( bud, status )
  
    use LE_Budget_DryDepos, only : DryDepos_Budget_Reset_Hour
    use LE_Budget_WetDepos, only : WetDepos_Budget_Reset_Hour
  
    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(inout)           ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Budget_Reset_Hour'
    
    ! --- begin ----------------------------------
    
    ! Reset with budgets:
    call DryDepos_Budget_Reset_Hour( bud%drydepos, status )
    IF_NOTOK_RETURN(status=1)
    call WetDepos_Budget_Reset_Hour( bud%wetdepos, status )
    IF_NOTOK_RETURN(status=1)
    
    ! reset maximum:
    bud%o3max = 0.0

    ! ok
    status = 0
    
  end subroutine Budget_Reset_Hour


  ! ***
  
  
  subroutine Budget_Reset_Day( bud, status )
  
    use LE_Budget_DryDepos, only : DryDepos_Budget_Reset_Day
    use LE_Budget_WetDepos, only : WetDepos_Budget_Reset_Day
  
    ! --- in/out ---------------------------------
    
    type(T_Budget), intent(inout)           ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Budget_Reset_Day'
    
    ! --- begin ----------------------------------
    
    ! Reset with budgets:
    call DryDepos_Budget_Reset_Day( bud%drydepos, status )
    IF_NOTOK_RETURN(status=1)
    call WetDepos_Budget_Reset_Day( bud%wetdepos, status )
    IF_NOTOK_RETURN(status=1)
    
    ! reset maximum:
    bud%o3max = 0.0

    ! ok
    status = 0
    
  end subroutine Budget_Reset_Day
  

end module LE_Budget

