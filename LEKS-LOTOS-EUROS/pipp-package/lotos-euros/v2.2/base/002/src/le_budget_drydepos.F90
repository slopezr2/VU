!###############################################################################
!
! LE_Budget_DryDepos  -   budget stuff
!
! There are budgets for different types of exchange processes; 
! nex     = 2: dry deposition, dry emission (in case of compensation point)
!
! There are three time intervals for which budgets are kept;
!  bud%ex_day(nx,ny,nbud,nex)  : daily budget for groups
!  bud%ex_hour(nx,ny,nspec,nex): hourly budget for each species
!  bud%ex_now(nx,ny,nspec,nex) : budget for current time step for each species
! 
! For dry deposition, separate budgets are computed per land use class for each species
!  bud%ex_day_lu(nx,ny,nspec,nlu_depac,nex_dry) : daily budget per land use class
!  bud%ex_hour_lu(nx,ny,nspec,nlu_depac,nex_dry): hourly budget per land use class
!  bud%ex_now_lu(nx,ny,nspec,nlu_depac) : budget for current time step per land use class
!
! The last index is reserved for the type of exchange process:
!  ex_drydepo = 1
!  ex_dryemis = 2
!  bud%ex...(...,ex_drydepo): dry deposition
!  bud%ex...(...,ex_dryemis): dry emission, in case there is a compensation point
!
! The "now" budgets are only present in subroutine depos; 
! hourly budgets are updated in subroutine depos;
! daily budgets are updated in DryDepos_Budget_update.

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

module LE_Budget_DryDepos

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out --------------------------------
  
  private

  public  ::  nex, ex_drydepo, ex_dryemis
  
  public  ::  T_DryDepos_Budget
  
  public  ::  DryDepos_Budget_Init, DryDepos_Budget_Done
  public  ::  DryDepos_Budget_Update
  public  ::  DryDepos_Budget_Reset_Hour
  public  ::  DryDepos_Budget_Reset_Day
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Budget_DryDepos'


  ! --- const ------------------------------------
  
  ! Set types of echange processes:
  ! wet deposition, dry deposition, 'dry emission' (in case of a compensation point)
  integer, parameter :: nex = 2
  integer, parameter :: ex_drydepo = 1
  integer, parameter :: ex_dryemis = 2


  ! --- types --------------------------------------
  
  type T_DryDepos_Budget
    ! deposition budgets:
    real, allocatable   ::  ex_hour   (:,:,:  ,:)  ! (nx,ny,nspec          ,nex)  [conc*m]
    real, allocatable   ::  ex_day    (:,:,:  ,:)  ! (nx,ny,nbud           ,nex)  [ug/m2]
    ! deposition budgets per landuse class:
    real, allocatable   ::  ex_hour_lu(:,:,:,:,:)  ! (nx,ny,nspec,nlu_depac,nex)  [conc*m]
    real, allocatable   ::  ex_day_lu (:,:,:,:,:)  ! (nx,ny,nspec,nlu_depac,nex)  [ug/m2]
    ! ozone flux through stomata:
    real, allocatable   ::  o3flx_stomata_hour_lu(:,:,:)    !  (nx,ny,nlu_depac) [nmol.m-2 PLA]
    real, allocatable   ::  o3flx_stomata_day_lu (:,:,:)    !  (nx,ny,nlu_depac) [umol.m-2 PLA]
    ! ozone flux through stomata exceedance:
    real, allocatable   ::  o3exc_stomata_hour_lu(:,:,:)    !  (nx,ny,nlu_depac) [nmol.m-2 PLA]
    real, allocatable   ::  o3exc_stomata_day_lu (:,:,:)    !  (nx,ny,nlu_depac) [umol.m-2 PLA]
    !
    ! * compensation point for NH3
    !
    ! summed concentration up till current time step (ppb)
    real, allocatable         :: cnh3_sum(:,:)
    ! number of time steps for which cnh3_sum is accumulated.
    integer                   :: cnh3_nsum
    ! average NH3 concentration over previous month (ppb)
    real, allocatable         :: cnh3_ave_prev(:,:)
    !
    ! * accumulated so2 concentration (co-deposition)
    !
    ! summed concentration up till current time step (ppb)
    real, allocatable         :: cso2_sum(:,:)
    ! number of time steps for which cnh3_sum is accumulated.
    integer                   :: cso2_nsum
    ! average SO2 concentration over previous month (ppb)
    real, allocatable         :: cso2_ave_prev(:,:)
  end type T_DryDepos_Budget


contains



  ! ========================================================================
  

  subroutine DryDepos_Budget_Init( bud, status )
  
    use Dims           , only : nx, ny, nspec
    use Indices        , only : i_NH3, i_SO2
    use LE_Landuse_Data, only : nlu
    use LE_Budget_Param, only : nbud
    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    type(T_DryDepos_Budget), intent(out)    ::  bud
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/DryDepos_Budget_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------

    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! exchange:
    allocate( bud%ex_hour   (nx,ny,nspec          ,nex) ) ; bud%ex_hour = 0.0
    allocate( bud%ex_day    (nx,ny,nbud           ,nex) ) ; bud%ex_day  = 0.0
    ! exchange per landuse class: 
    allocate( bud%ex_hour_lu(nx,ny,nspec,nlu,nex) ) ; bud%ex_hour_lu = 0.0
    allocate( bud%ex_day_lu (nx,ny,nspec,nlu,nex) ) ; bud%ex_day_lu  = 0.0
    
    ! O3 flux through stomata:
    allocate( bud%o3flx_stomata_hour_lu(nx,ny,nlu) ) ; bud%o3flx_stomata_hour_lu = 0.0
    allocate( bud%o3flx_stomata_day_lu (nx,ny,nlu) ) ; bud%o3flx_stomata_day_lu  = 0.0
    ! O3 flux through stomata, exceedance over threshold:
    allocate( bud%o3exc_stomata_hour_lu(nx,ny,nlu) ) ; bud%o3exc_stomata_hour_lu = 0.0
    allocate( bud%o3exc_stomata_day_lu (nx,ny,nlu) ) ; bud%o3exc_stomata_day_lu  = 0.0

    ! * NH3 compensation point
    
    ! enabled ?
    if ( i_nh3 > 0  ) then
      ! accumulation array for NH3 concentrations:    
      allocate( bud%cnh3_sum(nx,ny) )
      ! nothing added yet; init sum and counter to zero:
      bud%cnh3_sum      = 0.0
      bud%cnh3_nsum     = 0
      ! average of previous month:
      allocate( bud%cnh3_ave_prev(nx,ny) )
      ! init to very small number; not really necessary,
      ! but a world without any NH3 molecule is not very
      ! realistic either ...
      bud%cnh3_ave_prev = 1.0e-4   ! ppb?
    end if ! nh3 enabled
    
    if ( i_so2 > 0 ) then
      ! accumulation array for SO2 concentrations:    
      allocate( bud%cso2_sum(nx,ny) )
      ! nothing added yet; init sum and counter to zero:
      bud%cso2_sum      = 0.0
      bud%cso2_nsum     = 0
      ! average of previous month:
      allocate( bud%cso2_ave_prev(nx,ny) )
      ! init to very small number; not really necessary,
      ! but a world without any NH3 molecule is not very
      ! realistic either ...
      bud%cso2_ave_prev = 1.0e-4   ! ppb?
    end if ! so2 enabled

    
    ! *

    ! ok
    status = 0
    
  end subroutine DryDepos_Budget_Init
  
  
  ! ***
  
  
  subroutine DryDepos_Budget_Done( bud, status )
  
    use Indices        , only : i_NH3, i_SO2

    ! --- in/out ---------------------------------
    
    type(T_DryDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/DryDepos_Budget_Done'
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( bud%ex_hour )
    deallocate( bud%ex_day  ) 
    deallocate( bud%ex_hour_lu )
    deallocate( bud%ex_day_lu  ) 

    ! clear:
    deallocate( bud%o3flx_stomata_hour_lu )
    deallocate( bud%o3flx_stomata_day_lu  )
    deallocate( bud%o3exc_stomata_hour_lu )
    deallocate( bud%o3exc_stomata_day_lu  )

    ! * NH3 compensation point
    
    ! enabled ?
    if ( i_nh3 > 0 ) then
      deallocate( bud%cnh3_sum )
      deallocate( bud%cnh3_ave_prev )
    end if   ! NH3 enabled
    
    if (i_so2 > 0 ) then
      deallocate( bud%cso2_sum )
      deallocate( bud%cso2_ave_prev )
    end if ! SO2 enabled
    ! *
    
    ! ok
    status = 0
    
  end subroutine DryDepos_Budget_Done


  ! ***
  
  
  ! Update mass budgets

  subroutine DryDepos_Budget_Update( bud, t, status ) 

    use Binas          , only : xm_O3, xm_air
    use GO             , only : TDate, MidNight
    use Indices        , only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
    use Indices        , only : i_NH3, i_SO2
    use LE_Budget_Param, only : nbud, ispec_bud
    use Dims           , only : nx, ny

    use LE_Data        , only : LE_Data_GetPointer

    ! --- in/out ---------------------------------
    
    type(T_DryDepos_Budget), intent(inout)  ::  bud
    type(TDate), intent(in)                 ::  t
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/DryDepos_Budget_Update'
    
    ! --- local ----------------------------------
    
    integer             ::  ibud
    integer             ::  ispec
    integer             ::  icomp
    integer             ::  iex
    real, allocatable   ::  convfact(:,:)

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
    allocate( convfact(nx,ny) )  

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
          if ( accum_ppb_to_ugm3(ispec,icomp) ) convfact = dens(1:nx,1:ny,1)/xm_air
          ! add contribution:
          ! factors in accum_ww include conversion from ppb to ug/m3 if necessary:
          do iex = 1, nex
            bud%ex_day(:,:,ibud,iex) = bud%ex_day(:,:,ibud,iex) &
                     + bud%ex_hour(:,:,accum_ii(ispec,icomp),iex) * convfact * accum_ww(ispec,icomp)
          end do
        end do
      end if
    end do

    ! Update daily dry budgets per landuse category:
    bud%ex_day_lu = bud%ex_day_lu + bud%ex_hour_lu
    
    ! *

    ! calculate ozone fluxes and budget
    ! land-use dependent ozone stomata flux budget, 
    ! units in O3 mole (from [nmol m-2 PLA] to [umol O3.m-2 PLA])
    !                             umol O3.m-2 PLA           nmol.m-2 PLA               nmol/umol
    bud%o3flx_stomata_day_lu = bud%o3flx_stomata_day_lu + bud%o3flx_stomata_hour_lu/1000

    ! land-use dependent ozone stomata flux exceedence accumulated per day, 
    !  units in [umol m-2 PLA]
    bud%o3exc_stomata_day_lu = bud%o3exc_stomata_day_lu + bud%o3exc_stomata_hour_lu/1000
    
    ! * NH3 budgets

    ! Compute averaged concentration of previous month and reset; only for NH3.
    ! Is needed for the computation of the compensation point in DEPAC.

    ! enabled ?
    if ( i_nh3 > 0 .and. i_so2 > 0) then

      ! end of current day ?
      if ( MidNight(t) ) then
        ! new month ?
        if ( t%day == 1 ) then

          !! xxx debug printing
          !outF%debug_print = .true.
          !outF%debug_ix = 6
          !outF%debug_iy = 14

          ! Compute average over previous month:
          bud%cnh3_ave_prev = bud%cnh3_sum / bud%cnh3_nsum   ! ppb
          bud%cso2_ave_prev = bud%cso2_sum / bud%cso2_nsum   ! ppb

          !if (outF%debug_print) then
          !  ix = outF%debug_ix
          !  iy = outF%debug_iy
          !  write(*,'(a,1x,a6,4(1x,i4),2(1x,e12.5))') &
          !  'drydep init_month; spec,ix,iy,lu,cnh3_nsum,cnh3_sum,cnh3_ave_prev: ', &
          !  'nh3 ',ix,iy,0,cnh3_nsum,cnh3_sum(ix,iy),cnh3_ave_prev(ix,iy)
          !endif

          ! Reset:
          bud%cnh3_sum  = 0.0
          bud%cnh3_nsum = 0
          bud%cso2_sum  = 0.0
          bud%cso2_nsum = 0

          !! xxx debug printing
          !outF%debug_print = .false.

        endif  ! first day of month
      end if  ! midnight

    end if  ! NH3 enabled
    
    ! *
    
    ! clear:
    deallocate( convfact, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine DryDepos_Budget_Update
  
  
  ! ***
  
  
  subroutine DryDepos_Budget_Reset_Hour( bud, status )
  
   ! --- in/out ---------------------------------
    
    type(T_DryDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/DryDepos_Budget_Reset_Hour'
    
    ! --- begin ----------------------------------

    ! reset all daily budgets to zero:
    bud%ex_hour  = 0.0 
    bud%ex_hour_lu  = 0.0 
    bud%o3flx_stomata_hour_lu  = 0.0
    bud%o3exc_stomata_hour_lu  = 0.0

    ! ok
    status = 0
    
  end subroutine DryDepos_Budget_Reset_Hour
  
  
  ! ***
  
  
  subroutine DryDepos_Budget_Reset_Day( bud, status )
  
    ! --- in/out ---------------------------------
    
    type(T_DryDepos_Budget), intent(inout)  ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/DryDepos_Budget_Reset_Day'
    
    ! --- begin ----------------------------------

    ! reset all daily budgets to zero:
    bud%ex_day  = 0.0 
    bud%ex_day_lu  = 0.0 
    bud%o3flx_stomata_day_lu  = 0.0
    bud%o3exc_stomata_day_lu  = 0.0

    ! ok
    status = 0
    
  end subroutine DryDepos_Budget_Reset_Day
  

end module LE_Budget_DryDepos

