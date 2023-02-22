!#################################################################
!
! LE_Budget_WetDepos  -   budget stuff
!
! There are budgets for different types of exchange processes; 
! nex     = 1: wet deposition
!
! There are three time intervals for which budgets are kept;
!  bud%ex_day (nx,ny,nbud ,nex) : daily budget for groups
!  bud%ex_hour(nx,ny,nspec,nex) : hourly budget for each species
!  bud%ex_now (nx,ny,nspec,nex) : budget for current time step for each species
! 
! For wet deposition, separate budgets are computed per land use class for each species
!  bud%ex_day_lu (nx,ny,nspec,nlu_depac,nex_wet) : daily budget per land use class
!  bud%ex_hour_lu(nx,ny,nspec,nlu_depac,nex_wet): hourly budget per land use class
!  bud%ex_now_lu (nx,ny,nspec,nlu_depac) : budget for current time step per land use class
!
! The last index is reserved for the type of exchange process:
!  ex_wetdepo = 1
!  bud%ex...(...,ex_wetdepo): wet deposition
!
! The "now" budgets are only present in subroutine depos; 
! hourly budgets are updated in subroutine depos;
! daily budgets are updated in WetDepos_Budget_update.
!
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

module LE_Budget_WetDepos

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out --------------------------------
  
  private

  public  ::  ex_wetdepo

  public  ::  T_WetDepos_Budget
  
  public  ::  WetDepos_Budget_Init, WetDepos_Budget_Done
  public  ::  WetDepos_Budget_Update
  public  ::  WetDepos_Budget_Reset_Hour
  public  ::  WetDepos_Budget_Reset_Day
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Budget_WetDepos'


  ! --- const ------------------------------------
  
  ! Set types of echange processes:
  ! wet deposition, Wet deposition, 'Wet emission' (in case of a compensation point)
  integer, parameter :: nex = 1
  integer, parameter :: ex_wetdepo = 1


  ! --- types --------------------------------------
  
  type T_WetDepos_Budget
    ! hourly budget (mass per grid cell area)  (ppb.m ug/m2), integrated in time
    real, allocatable   ::  ex_hour   (:,:,:,:,:)  ! (nx,ny,nz,nspec,nex)
    real, allocatable   ::  ex_day    (:,:,:,:)  ! (nx,ny,nbud ,nex)
  end type T_WetDepos_Budget


contains



  ! ========================================================================
  

  subroutine WetDepos_Budget_Init( bud, status )
  
    use Dims           , only : nx, ny, nz, nspec
    use LE_Budget_Param, only : nbud
    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    type(T_WetDepos_Budget), intent(out)    ::  bud
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/WetDepos_Budget_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------

    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( bud%ex_hour(nx,ny,nz,nspec,nex) )
    allocate( bud%ex_day (nx,ny,nbud ,nex) )
    ! init to zero:
    bud%ex_hour = 0.0
    bud%ex_day  = 0.0
    
    ! ok
    status = 0
    
  end subroutine WetDepos_Budget_Init
  
  
  ! ***
  
  
  subroutine WetDepos_Budget_Done( bud, status )
  
    ! --- in/out ---------------------------------
    
    type(T_WetDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/WetDepos_Budget_Done'
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( bud%ex_hour )
    deallocate( bud%ex_day  ) 

    ! ok
    status = 0
    
  end subroutine WetDepos_Budget_Done


  ! ***
  
  
  ! Update mass budgets

  subroutine WetDepos_Budget_Update( bud, status ) 

    use Binas          , only : xm_air
    use Indices        , only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
    use LE_Budget_Param, only : nbud, ispec_bud
    use Dims           , only : nx, ny, nz

    use LE_Data               , only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    type(T_WetDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/WetDepos_Budget_Update'
    
    ! --- local ----------------------------------
    
    integer             ::  ibud
    integer             ::  ispec
    integer             ::  icomp
    integer             ::  iex
    real, allocatable   ::  convfact(:,:,:)

    ! meteo data:
    real, pointer       ::  dens(:,:,:)   ! (lon,lat,alt)            

    ! --- begin ----------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')            
    IF_NOTOK_RETURN(status=1)

    ! Update daily budgets.
    ! Note that budgets are computed as concentration*layer_height in
    ! ppb m (gas) or ug/m2 (aerosol). The multiplication factor 'accum_ww'
    ! is either 1.0 (aerosols in ug/m3), or converts from ppb to ug/m3,
    ! which ensures that all budgets are expressed in ug/m2.
    
    ! storage:
    allocate( convfact(nx,ny,nz) )

    ! loop over all budget species:
    do ibud = 1, nbud
      ! corresponding accumulated tracer:
      ispec = ispec_bud(ibud)
      ! enabled ?
      if ( ispec > 0 ) then
        ! loop over contributing components:
        do icomp = 1, accum_n(ispec)
          ! conversion needed ?
          convfact = 1.0
          if ( accum_ppb_to_ugm3(ispec,icomp) ) convfact = dens(1:nx,1:ny,1:nz)/xm_air
          ! loop over exchange processes:
          do iex = 1, nex
            ! add contribution; 
            ! factors in accum_ww include conversion from ppb to ug/m3 if necessary:
            bud%ex_day(:,:,ibud,iex) = bud%ex_day(:,:,ibud,iex) &
                   + sum(bud%ex_hour(:,:,:,accum_ii(ispec,icomp),iex) * convfact * accum_ww(ispec,icomp),dim=3)
          end do
        end do
      end if
    end do
    !! reset hourly budgets to zero:
    !bud%ex_hour  = 0.0 
    
    ! clear:
    deallocate( convfact )
        
    ! ok
    status = 0
    
  end subroutine WetDepos_Budget_Update
  
  
  ! ***
  
  
  subroutine WetDepos_Budget_Reset_Hour( bud, status )
  
    ! --- in/out ---------------------------------
    
    type(T_WetDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/WetDepos_Budget_Reset_Hour'
    
    ! --- begin ----------------------------------

    ! reset all daily budgets to zero:
    bud%ex_hour  = 0.0 



    ! ok
    status = 0
    
  end subroutine WetDepos_Budget_Reset_Hour
  
  
  ! ***
  
  
  subroutine WetDepos_Budget_Reset_Day( bud, status )
  
    ! --- in/out ---------------------------------
    
    type(T_WetDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/WetDepos_Budget_Reset_Day'
    
    ! --- begin ----------------------------------

    ! reset all daily budgets to zero:
    bud%ex_day  = 0.0 

    ! ok
    status = 0
    
  end subroutine WetDepos_Budget_Reset_Day
    
  
  
end module LE_Budget_WetDepos

