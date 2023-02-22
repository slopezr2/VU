!###############################################################################
!
! Meteo Data
!
! (To be moved to LE_Data structures)
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


module LE_Meteo_Data

  use GO, only : gol, goPr, goErr
  use grid   , only : TllGridInfo
  
  implicit none


  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_Meteo_Data_Init, LE_Meteo_Data_Done
  public  ::  LE_Meteo_Data_Setup

  public  ::  volume, ovolume
  public  ::  afluxx, afluxy, afluxz
  public  ::  cflx, cfly, cflz
  
  public  ::  airmass, oairmass
  public  ::  hpres, ohpres

  public  ::  freepathlen
  public  ::  viscos
  
  public  ::  day_mean_tsurf


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_LE_Meteo_Data'


  ! --- var --------------------------------------
  
  ! grid cell volumes, incl. single cell halo:
  real, allocatable     ::   volume(:,:,:)    ! (0:nx+1,0:ny+1,1:nz+1) [km3] (zucht)
  real, allocatable     ::  ovolume(:,:,:)    ! (0:nx+1,0:ny+1,1:nz+1) [km3] (zucht)
  
  ! copy of old volume in m3,
  ! and current and original airmasses;
  ! needed for unit conversions in 'adjust'
  real, allocatable, target     ::   airmass(:,:,:)  ! kg
  real, allocatable, target     ::  oairmass(:,:,:)  ! kg
  real, allocatable, target     ::   hpres(:,:,:)    ! Pa
  real, allocatable, target     ::  ohpres(:,:,:)    ! Pa

  ! advective fluxes:
  real, allocatable     ::  afluxx(:,:,:)    ! km3/min (zucht)
  real, allocatable     ::  afluxy(:,:,:)    ! km3/min (zucht)
  real, allocatable     ::  afluxz(:,:,:)    ! km3/min (zucht)

  ! courant numbers:
  real, allocatable     ::  cflx(:,:,:)    ! 1
  real, allocatable     ::  cfly(:,:,:)    ! 1
  real, allocatable     ::  cflz(:,:,:)    ! 1
  
  ! viscosity profiles:
  real, allocatable     ::  viscos(:,:,:)       ! (nx,ny,nz)
  ! free path length
  real, allocatable     ::  freepathlen(:,:,:)  ! (nx,ny,nz)

  ! daily mean surface temperature 
  real, allocatable     ::  day_mean_tsurf(:,:) ! deduced meteo field

  
contains



  ! ========================================================================
  

  subroutine LE_Meteo_Data_Init( rcF, status )
  
    use GO     , only : TrcFile
    use Dims   , only : nx, ny, nz
    use LE_Data, only : LE_Data_Enable
    use LE_Data, only : levels_type
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Meteo_Data_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'p', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'hp', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'airm', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'uflux', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'wflux', status )
    IF_NOTOK_RETURN(status=1)
   
    ! switch:
    select case ( trim(levels_type) )
      !~
      case ( 'mixlayer' )
        ! enable fields required for advection:
        ! ~ current volume (mid of timestep):
        call LE_Data_Enable( 'vol', status )
        IF_NOTOK_RETURN(status=1)
      !~
      case ( 'hyblevel', 'metlevel' )
        ! enable fields required for advection:
        ! ~ previous volume (at start of timestep):
        call LE_Data_Enable( 'pvol', status )
        IF_NOTOK_RETURN(status=1)
        ! ~ current volume (mid of timestep):
        call LE_Data_Enable( 'vol', status )
        IF_NOTOK_RETURN(status=1)
        ! ~ next volume (at end of timestep):
        call LE_Data_Enable( 'nvol', status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported level type `",a,"`")') trim(levels_type); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage for volumes:   
    allocate( volume(0:nx+1,0:ny+1,1:nz+1), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate( ovolume(0:nx+1,0:ny+1,1:nz+1), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)

    ! storage for (old) airmass:
    allocate(airmass(nx,ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate(oairmass(nx,ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)

    ! storage for (old) half-level ppressure:
    allocate( ohpres(nx,ny,0:nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)
    allocate( hpres(nx,ny,0:nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)
    
    ! fluxes through cell edges:
    allocate( afluxx(0:nx,ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate( afluxy(nx,0:ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate( afluxz(nx,ny,0:nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)

    ! storage for cfl criteria (diagnostic output):  
    allocate( cflx(nx,ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate( cfly(nx,ny,nz), stat=status, source=0.0 ) 
    IF_NOTOK_RETURN(status=1)
    allocate( cflz(nx,ny,nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)

    ! storage:    
    allocate( viscos(nx,ny,nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)
    allocate( freepathlen(nx,ny,nz), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)

    ! storage:    
    allocate( day_mean_tsurf(nx,ny), stat=status, source=0.0 )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Meteo_Data_Init
  
  
  ! ***
  
  
  subroutine LE_Meteo_Data_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Meteo_Data_Done'
    
    ! --- begin ----------------------------------
    
    ! clear:    
    deallocate( volume, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( ovolume, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( airmass, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( oairmass, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( hpres, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( ohpres, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( afluxx, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( afluxy, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( afluxz, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! clear:    
    deallocate( cflx, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( cfly, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( cflz, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( viscos, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( freepathlen, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( day_mean_tsurf, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Meteo_Data_Done
  
  
  ! ***
  
  
  ! update fields derived from meteo variables
  
  subroutine LE_Meteo_Data_Setup( t, dt, status )

    use GO                , only : TDate, Midnight
    use GO                , only : TIncrDate, rTotal
    
    use JAQL_DryDeposition, only : Free_Path_Length
    use JAQL_DryDeposition, only : Viscosity

    use LE_Data           , only : LE_Data_GetPointer    
    use Dims              , only : nx, ny, nz
    use Dims              , only : runF
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt
    integer, intent(out)          ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Meteo_Data_Setup'
    
    ! --- local ----------------------------------
    
    real, pointer     ::  temp(:,:,:)   ! (lon,lat,lev)
    real, pointer     ::  pres(:,:,:)   ! (lon,lat,lev)
    real, pointer     ::  tsurf(:,:,:)   ! (lon,lat,1)
    
    ! --- begin ----------------------------------

    ! copy volume to original array,
    ! fill original mass flux arrays:
    call MassFlux( runF%first, status )
    IF_NOTOK_RETURN(status=1)

    ! pointers to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units ='K', &
                                check_lbo=(/ 1, 1, 1/), &
                                check_ubo=(/nx,ny,nz/) )
    IF_NOTOK_RETURN(status=1)  
    call LE_Data_GetPointer( 'p', pres, status, check_units ='Pa', &
                                check_lbo=(/ 1, 1, 1/), &
                                check_ubo=(/nx,ny,nz/) )
    IF_NOTOK_RETURN(status=1)  
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K', &
                                check_lbo=(/ 1, 1,1/), &
                                check_ubo=(/nx,ny,1/) )
    IF_NOTOK_RETURN(status=1)

    ! update meteo dependent parameters:
    freepathlen = Free_Path_Length( temp(1:nx,1:ny,1:nz), pres(1:nx,1:ny,1:nz) )

    ! fill viscosity:
    viscos = Viscosity( temp(1:nx,1:ny,1:nz) )

    ! Build up daily mean temperature per day;
    ! reset at start of day:
    if ( MidNight(t) ) day_mean_tsurf = 0.0
    ! add contribution for this time step:
    day_mean_tsurf(1:nx,1:ny) = day_mean_tsurf(1:nx,1:ny) + tsurf(1:nx,1:ny,1) * rTotal(dt,'hour')/24.0
    
    ! ok
    status = 0
    
  end subroutine LE_Meteo_Data_Setup


  ! ***
  
  
  subroutine MassFlux( lfirst, status )

    use LE_Data, only : levels_type
    use LE_Data, only : LE_Data_GetPointer
    use Dims   , only : nx, ny, nz
    
    ! --- in/out --------------------------------------
    
    logical, intent(in)     ::  lfirst
    integer, intent(out)    ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MassFlux'
    
    ! conversion factors:
    real, parameter  ::  km3_per_m3 = 1.0e-9   ! km3/m3
    real, parameter  ::  s_per_min  = 60.0     ! s/min
    
    ! --- local ---------------------------------------

    real, pointer   ::  pdata(:,:,:)
    
    ! --- begin ---------------------------------------

    !
    ! copy volume and volume fluxes computed in LE_Data
    !

    ! switch:
    select case ( trim(levels_type) )
      !~
      case ( 'mixlayer' )
        ! if not first step: copy volume to ovolume
        if ( .not. lfirst ) then
          ovolume  = volume   ! km3 (zucht)
          oairmass = airmass
          ohpres    = hpres  
        end if

        ! get pointer to volume valid for mid of timestep:
        call LE_Data_GetPointer( 'vol', pdata, status, check_units='m3', &
                                    check_lbo=(/   0,   0,   1/), &
                                    check_ubo=(/nx+1,ny+1,nz+1/) )
        IF_NOTOK_RETURN(status=1)
        ! copy data, convert:
        !                m3               km3/m3
        volume = pdata(0:nx+1,0:ny+1,1:nz+1) * km3_per_m3   ! km3 (zucht)

        ! get pointer to airmass valid for mid of timestep:
        call LE_Data_GetPointer( 'airm', pdata, status, check_units='kg', &
                                    check_lbo=(/   0,   0,   1/), &
                                    check_ubo=(/nx+1,ny+1,nz+1/) )
        IF_NOTOK_RETURN(status=1)
        ! copy data:
        airmass = pdata(1:nx,1:ny,1:nz)  ! kg
    
        ! get pointer to half level pressure for mid of timestep
        call LE_Data_GetPointer( 'hp', pdata, status, check_units='Pa', &
                                    check_lbo=(/   0,   0,   0/), &
                                    check_ubo=(/nx+1,ny+1,nz+1/) )
        IF_NOTOK_RETURN(status=1)
        ! Copy data
        hpres = pdata(1:nx,1:ny,0:nz) ! Pa
    
        ! if first step, assum the old volume was equal to the new volume
        if ( lfirst ) then
          ovolume  = volume   ! km3 (zucht)
          oairmass = airmass  ! kg
          ohpres   = hpres
        end if
    
      !~
      case ( 'hyblevel', 'metlevel' )
      
        ! get pointer to volume valid for start of timestep:
        call LE_Data_GetPointer( 'pvol', pdata, status, check_units='m3', &
                                    check_lbo=(/   0,   0,   1/), &
                                    check_ubo=(/nx+1,ny+1,nz+1/) )
        IF_NOTOK_RETURN(status=1)
        ! copy data, convert:
        !           m3            km3/m3
        volume = pdata(0:nx+1,0:ny+1,1:nz+1) * km3_per_m3   ! km3 (zucht)
        
        ! volume will be changed during advection;
        ! store into old volume (not used anymore?)
        ovolume = volume   ! km3
        
      !~
      case default
        write (gol,'("unsupported level type `",a,"`")') trim(levels_type); call goErr
        TRACEBACK; status=1; return
    end select

    !! also keep version in decent units ...
    !ovolume_m3 = ovolume / km3_per_m3  ! m3
    
    ! get pointer:
    call LE_Data_GetPointer( 'uflux', pdata, status, check_units='m3/s', &
                                check_lbo=(/ 0, 1, 1/), &
                                check_ubo=(/nx,ny,nz/) )
    IF_NOTOK_RETURN(status=1)
    ! copy data, convert:
    !               m3/s               km3/m3       s/min
    afluxx = pdata(0:nx,1:ny,1:nz) * km3_per_m3 * s_per_min   ! km3/min (zucht)
    
    ! get pointer:
    call LE_Data_GetPointer( 'vflux', pdata, status, check_units='m3/s', &
                                check_lbo=(/ 1, 0, 1/), &
                                check_ubo=(/nx,ny,nz/) )
    IF_NOTOK_RETURN(status=1)
    ! copy data, convert:
    !               m3/s               km3/m3       s/min
    afluxy = pdata(1:nx,0:ny,1:nz) * km3_per_m3 * s_per_min   ! km3/min (zucht)
    
    ! get pointer:
    call LE_Data_GetPointer( 'wflux', pdata, status, check_units='m3/s', &
                                check_lbo=(/ 1, 1, 0/), &
                                check_ubo=(/nx,ny,nz/) )
    IF_NOTOK_RETURN(status=1)
    ! copy data, convert:
    !               m3/s               km3/m3       s/min
    afluxz = pdata(1:nx,1:ny,0:nz) * km3_per_m3 * s_per_min   ! km3/min (zucht)
    
    ! ok
    status = 0
    
  end subroutine MassFlux


end module LE_Meteo_Data

