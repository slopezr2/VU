!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################

module Pollen_Data_Grass

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  use GO, only : TDate

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Pollen_Data_Grass
  public  ::  Pollen_Data_Grass_Init, Pollen_Data_Grass_Done
  public  ::  Pollen_Data_Grass_Get

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'Pollen_Data_Grass'
  

  ! --- local --------------------------------

  ! version number for source function:
  integer                             :: polg_ver

  ! emission data base
  type T_Pollen_Data_Grass
    ! label assigned to this emission:
    character(len=32)                 ::  label
    ! settings:
    character(len=64)                 ::  rckey
    ! input file template:
    character(len=1024)               ::  fname_template
    ! current month:
    integer                           ::  year
    ! emisison field:
    real, allocatable                 ::  emis(:,:)  ! (nx,ny)
    real, allocatable                 ::  ripened_left(:,:) ! (nx,ny)
    real, allocatable                 ::  max_avail_grains(:,:)  ! (nx,ny)
    real, allocatable                 ::  rest_avail_grains(:,:)  ! (nx,ny)
    real, allocatable                 ::  intensity(:,:) ! (nx,ny)
    character(len=32)                 ::  emis_units
    character(len=32)                 ::  tracer
    character(len=32)                 ::  source
    integer                           ::  itr_glob  ! global tracer index

    character(len=512)                ::  file
    character(len=512)                ::  parameterfile
    integer                           ::  nlon, nlat
    
    real, allocatable                 ::  lons(:)
    real, allocatable                 ::  lats(:)
    
    real, allocatable                 ::  landuse_fracs(:,:) ! (nx,ny)
        
    real, allocatable                 ::  startday_em(:,:)  ! (nx,ny) startday of emitting pollen (for grass, ..)
    real, allocatable                 ::  endday_em(:,:)   ! (nx,ny) end day of emitting pollen (for grass, ..)
    real, allocatable                 ::  nday_em(:,:)     ! (nx,ny) length of emitting pollen season (for grass, ..)
    
    real, allocatable                 ::  gamma_field(:,:)    
    real                              ::  low_hum_th, high_hum_th, prec_th  ! humidity and preciptation threshold
    character(len=8)                  ::  low_hum_th_unit, high_hum_th_unit, prec_th_unit

    real                              ::  wspd_sat_lev, wspd_max_impact     ! windspeed constraints
    character(len=8)                  ::  wspd_sat_lev_unit
    real                              ::  max_avail_grains_m2    ! available grains per m2
    
    real                              ::  unc_max_avail_grains
    character(len=8)                  ::  unc_max_avail_grains_unit
    real                              ::  unc_startday
    
    logical                           ::  allow_zero
    
  end type T_Pollen_Data_Grass



contains


  ! ===============================================================


  subroutine Pollen_Data_Grass_Init( emg, rcF, rckey, t, status )
    
    use GO , only : ReadRc
    use GO , only : Set
    use GO , only : goSplitString, goReplace, operator(<), operator(<=), DayNumber
    use GO , only : TDate, Get
    use GO , only : goVarValue

    use LE_Grid         , only : ugg
    use C3PO            , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common  , only : Grid_Convertors
    use LE_Data         , only : LE_Data_Enable
    use LE_Restart_Tools, only : LE_Restart_Restore

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Grass), intent(inout)  ::  emg
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rckey
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Grass_Init'

    ! --- local ------------------------------

    character(len=128), allocatable :: descriptionlist(:) ! list of variables needed per pollen tracer
    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  values_in(:,:)
    character(len=64)     ::  units_in
    real                  ::  missing_value
    
    integer               ::  polg_ver
    integer               ::  idesc
    
    character(len=256)    ::  restart_path, restart_key
    logical               ::  do_restart
    
    ! --- begin -------------------------------

    ! enable data:
    call LE_Data_Enable( 'rain', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'wspd_surf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'rh', status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.version', polg_ver, status )
    IF_NOTOK_RETURN(status=1)

    ! read the correct input file for grass pollen data (mask, stard/end-days etc.)
    select case ( polg_ver )
      case ( 2 )  
        call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.file', emg%file, status )
        IF_NOTOK_RETURN(status=1)
      case ( 3 )
        call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.file3', emg%file, status )
        IF_NOTOK_RETURN(status=1)
      case default 
        write(gol, '("unsupported version for grass pollen data and code: ")' )  ; call goErr
        TRACEBACK;status=1;return
    end select
 
    ! year of emissions
    emg%year = t%year

    ! find landuse fractions (for birch the lon.lat naming is different compared with grass, olive
    if (polg_ver == 2 ) then 

      ! start/length of emitting period (grid is different)
      allocate( descriptionlist(3) )
      descriptionlist = (/ 'standard_name=StartCD         ', &
                           'standard_name=lenCD           ', &
                           'standard_name=Grass Percentage' /)
    else if ( polg_ver == 3 ) then
      allocate( descriptionlist(3) )
      descriptionlist = (/ 'var_name=start_calday_th_POLLEN_GRASS_m32', &
                           'var_name=end_calday_th_POLLEN_GRASS_m32  ', &
                           'var_name=emis_mask_POLLEN_GRASS_m32      ' /)
    end if
    
    ! loop over input variables
    do idesc = 1, 3
      description = descriptionlist(idesc)

      call file_in%Open( trim(emg%file), status )
      IF_NOTOK_RETURN(status=1)

      ! variable id:
      call file_in%Inq_VarID( trim(description), varid, status )
      IF_NOTOK_RETURN(status=1)
      ! init grid definition
      call file_in%Get_Grid( varid, grid_in, status )
      IF_NOTOK_RETURN(status=1)

      ! storage:
      allocate( values_in( grid_in%nlon, grid_in%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! read: 
      call file_in%Get_Var( trim(description), values_in, units_in, status )
      IF_NOTOK_RETURN(status=1)
      ! define target array
      if ( trim(description) == 'standard_name=StartCD' .or. trim(description) == 'var_name=start_calday_th_POLLEN_GRASS_m32') then
        ! check units
        if ( trim(units_in) == 'Julian day' .or. trim(units_in) == 'day' ) then
          allocate( emg%startday_em( ugg%nlon, ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emg%startday_em, status )
          IF_NOTOK_RETURN(status=1)    
        else 
          write( gol, '("Unsupported unit: ", a " for startCD field in file: ", a )' ) trim(units_in), trim(emg%file) ; call GoErr
          TRACEBACK;status=1;return
        end if
      else if ( trim(description) == 'standard_name=lenCD'  ) then
        ! check units
        if ( trim(units_in) == 'days' ) then
          allocate( emg%nday_em( ugg%nlon, ugg%nlat),  stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emg%nday_em, status )
          IF_NOTOK_RETURN(status=1)    
        else 
          write( gol, '("Unsupported unit: ", a " for lenCD field in file: ", a )' ) trim(units_in), trim(emg%file) ; call GoErr
          TRACEBACK;status=1;return
        end if
      else if ( trim(description) == 'var_name=end_calday_th_POLLEN_GRASS_m32'  ) then
        ! check units
        if ( trim(units_in) == 'day' ) then
          allocate( emg%endday_em( ugg%nlon, ugg%nlat),  stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emg%endday_em, status )
          IF_NOTOK_RETURN(status=1)    
        else 
          write( gol, '("Unsupported unit: ", a " for endCD field in file: ", a )' ) trim(units_in), trim(emg%file) ; call GoErr
          TRACEBACK;status=1;return
        end if
      else if ( trim(description) == 'standard_name=Grass Percentage' .or. trim(description) == 'var_name=emis_mask_POLLEN_GRASS_m32' ) then
        ! check units
        if ( trim(units_in) == '%' ) then
          ! convert to fractions:
          values_in = values_in * 0.01
          allocate( emg%landuse_fracs( ugg%nlon, ugg%nlat),  stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emg%landuse_fracs, status )
          IF_NOTOK_RETURN(status=1)
          emg%landuse_fracs = max( emg%landuse_fracs, 0.0)    
        else if ( trim(units_in) == 'eff.fr.' ) then
          allocate( emg%landuse_fracs( ugg%nlon, ugg%nlat),  stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emg%landuse_fracs, status )
          IF_NOTOK_RETURN(status=1)
          emg%landuse_fracs = max( emg%landuse_fracs, 0.0) 
        
        else    
          ! unsupported units
          write( gol, '("Unsupported unit: ", a " for landuse fraction field in file: ", a )' ) trim(units_in), trim(emg%file) ; call GoErr
          TRACEBACK;status=1;return
        end if
      else 
        write( gol, '("Unsupprted description: ", a, " for file: " )' ) trim(description), trim(emg%file) ; call GoErr
        TRACEBACK;status=1;return
      end if

      ! clear: 
      deallocate( values_in )
      ! Done with pollen input grid
      call grid_in%Done(status)
      IF_NOTOK_RETURN(status=1)

      ! close:
      call file_in%Close( status )
      IF_NOTOK_RETURN(status=1)

    end do ! descriptions
    
    if ( polg_ver == 2 ) then
      ! fill in sarray with end days of emitting period
      allocate( emg%endday_em(ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      emg%endday_em = emg%startday_em + emg%nday_em 
    else if ( polg_ver == 3 ) then
      allocate( emg%nday_em( ugg%nlon, ugg%nlat ), stat=status )
      IF_NOTOK_RETURN(status=1)
      emg%nday_em = emg%endday_em - emg%startday_em
    endif !polg_ver?

    ! clear:
    deallocate( descriptionlist )
    
    !       
    allocate( emg%gamma_field(ugg%nlon,ugg%nlat) )
    emg%gamma_field = 0.0

    ! fill in intensity array with only values 1        
    allocate( emg%intensity(ugg%nlon,ugg%nlat) )
    emg%intensity = 1.0
        
    ! initialize emissions
    allocate(emg%emis(ugg%nlon,ugg%nlat) )
    emg%emis = 0.0   
    
    allocate(emg%ripened_left(ugg%nlon,ugg%nlat) )
    emg%ripened_left = 0.0
    
    ! Read in other parameters
        
    ! humiduty tresholds, below low_th:full emission, above high_th: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.low_hum_th', emg%low_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.low_hum_th_unit', emg%low_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.high_hum_th', emg%high_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.high_hum_th_unit', emg%high_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check humidty units
    if ( emg%low_hum_th_unit == '%' ) then
      emg%low_hum_th = emg%low_hum_th   ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emg%high_hum_th_unit
      TRACEBACK;status=1;return
    endif      
    if ( emg%high_hum_th_unit == '%' ) then
      emg%high_hum_th = emg%high_hum_th ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emg%high_hum_th_unit
      TRACEBACK;status=1;return
    endif
    
    ! preciptation threshold, above: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.prec_th', emg%prec_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.prec_th_unit', emg%prec_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! check units
    if ( emg%prec_th_unit == 'mm/hr' ) then
      ! LE unit is in m/s;
      !    m/s          mm/hr   /  mm/m  / (s/hr)
      emg%prec_th = emg%prec_th / 1000.0 / 3600.0
    else
      write(gol, '( "wrong unit for preciptation threshold: ", a )' ) emg%prec_th_unit
      TRACEBACK;status=1;return
    endif     
    
    ! windspeed constraints
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.wspd_sat_lev', emg%wspd_sat_lev, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.wspd_sat_lev_unit', emg%wspd_sat_lev_unit, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.wspd_max_impact', emg%wspd_max_impact, status )
    IF_ERROR_RETURN(status=1)    
    
    ! check units
    if ( emg%wspd_sat_lev_unit == 'm/s' ) then
      emg%wspd_sat_lev = emg%wspd_sat_lev  ! wsurf (u10) LE unit is in m/s  !!!! note that u and v are in km/min in LE
    else
      write(gol, '( "wrong unit for wind saturation level threshold: ", a )' ) emg%wspd_sat_lev_unit
      TRACEBACK;status=1;return
    endif    
    
    ! maximum available grains per m2
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.max_avail_grns', emg%max_avail_grains_m2, status )
    IF_ERROR_RETURN(status=1)    
    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.unc_max_avail_grns', emg%unc_max_avail_grains, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.unc_max_avail_grns_unit', emg%unc_max_avail_grains_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! uncertainties in start of growing season and nr ov available grains
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.unc_startday', emg%unc_startday, status )    
    IF_ERROR_RETURN(status=1)
    
    if ( emg%unc_max_avail_grains_unit == '%' ) then
      emg%unc_max_avail_grains = 0.01 * emg%unc_max_avail_grains ! uncertainty in fraction
    else
      write(gol, '( "wrong unit for uncertainty of max available grains: ", a )' ) emg%unc_max_avail_grains_unit ; call goErr
      TRACEBACK;status=1;return
    endif       
                
    ! available grains per grid cell
    allocate(emg%max_avail_grains(ugg%nlon,ugg%nlat) )
    emg%max_avail_grains = emg%max_avail_grains_m2 * ugg%area * emg%landuse_fracs * emg%intensity
    
    allocate(emg%rest_avail_grains(ugg%nlon,ugg%nlat) )
    emg%rest_avail_grains = emg%max_avail_grains_m2 * ugg%area * emg%landuse_fracs * emg%intensity
        
    call ReadRc( rcF, 'le.restart', do_restart, status )
    IF_NOTOK_RETURN(status=1)
        
    ! allow pollen zero for the rest of the year?
    call ReadRc( rcF, trim(rckey)//'.'//trim(emg%tracer)//'.allow_zero', emg%allow_zero, status)
    IF_NOTOK_RETURN(status=1)
    
    if ( all( DayNumber(t) < emg%startday_em(:,:) - emg%unc_startday .and. emg%startday_em > 0.0 ) .and. (.not.do_restart) ) then    
      ! do nothing
      ! normal calculation of Pollen season
      !emg%rest_avail_grains = 0.0
      !emg%ripened_left = 0.0
    else 
      if ( emg%allow_zero ) then      
        write( gol, '(" Allow rest amount of Pollen to be zero for: ", a )' ) trim(emg%tracer) ; call GoPr
        emg%rest_avail_grains = 0.0
        emg%ripened_left = 0.0

      else if ( any( DayNumber(t) > emg%startday_em(:,:) - emg%unc_startday .and. emg%startday_em > 0.0 ) .and. (.not.do_restart) ) then
        write(gol, '(" Restart option is not enabled")' ) ; call GoErr
        write(gol, '(" model run starts after begin of growing season for Pollen")' )  ; call GoErr
        write(gol, '(" so information of currently released amount of pollen must be in restart file")' ) ; call GoErr
        write(gol, '(" Start day of pollen season: ", i4)' ) nint( minval( emg%startday_em ) ) ; call goErr
        write(gol, '(" Uncertainty in start time is: ", i3)' ) nint(emg%unc_startday) ; call goErr
        write(gol, '(" Currentmodel time: ", i4, i2, i2)' ) t%year, t%month, t%day ; call goErr
        TRACEBACK; status=1; return

      else if ( do_restart ) then

        call ReadRc( rcF, 'le.restart.path', restart_path, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'le.restart.key', restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of pollen left ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'amt_polg_left', emg%rest_avail_grains, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of ripened pollen left ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'ripened_polg_left', emg%ripened_left, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

      end if
    end if  
       
    ! ok
    status = 0

  end subroutine Pollen_Data_Grass_Init

  
  ! ***


  subroutine Pollen_Data_Grass_Done( emg, status )

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Grass), intent(inout)    ::  emg
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Grass_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! clear ?
    deallocate( emg%emis              )
    deallocate( emg%ripened_left      ) 
    deallocate( emg%landuse_fracs     ) 
    deallocate( emg%startday_em       )
    deallocate( emg%endday_em         )
    deallocate( emg%nday_em           )
    deallocate( emg%rest_avail_grains )

    !  ok
    status = 0

  end subroutine Pollen_Data_Grass_Done


  ! ***

  subroutine Pollen_Data_Grass_Get( emg, t1, t2, status )
    
    use GO, only : TDate, TIncrDate, operator(-), operator(<), operator(>)
    use GO, only : operator(<=), Midnight, DayNumber
    use dims, only : nx, ny
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : day_mean_tsurf
    
    ! --- in/out ---------------------------

    type(T_Pollen_Data_Grass), intent(inout)    ::  emg
    type(TDate), intent(in)               ::  t1, t2   
    integer, intent(out)                  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Grass_Get'
    
    ! --- local ----------------------------
    
    integer                     ::  ix, iy
    type(TIncrDate)             ::  dt
    real                        ::  nmin
    real                        ::  season_ratio, gamma_f, nsec, delta
    real                        ::  ripened_grains, constr
    real                        ::  fu_fade_in, fu_fade_out, fu_fade_out_days

    real, pointer               ::  rain (:,:,:)   ! (nlon,nlat,1)
    real, pointer               ::  wsurf(:,:,:)   ! (nlon,nlat,1)
    real, pointer               ::  rh   (:,:,:)   ! (nlon,nlat,nz)
    
    ! --- begin ----------------------------
    ! meteo data:

    call LE_Data_GetPointer( 'rain'     , rain , status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'wspd_surf', wsurf, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh'       , rh   , status, check_units ='%' )
    IF_NOTOK_RETURN(status=1)
    
    dt = t2-t1

    ! number of minutes and seconds in dt
    nmin = dt%day*24*60 + dt%hour*60 + dt%min + ( dt%sec / 60. )
    nsec = dt%day*24*60*60 + dt%hour*60*60 + dt%min*60 + dt%sec

    polg_ver = 3

    ! update cumulative pollen release (R(t) in paper: Agricultural and Forest meteorology)
    do ix = 1, nx
    do iy = 1, ny 

      if ( emg%startday_em(ix,iy) < 0.0 ) then
        ! no information about start/end emitting season
        cycle
      end if

      select case ( polg_ver )  
        case (2) !Linear fade-in/out

          if ( DayNumber(t2) < emg%startday_em(ix,iy) - emg%unc_startday ) then
            ! pollen not yet released
            cycle
          end if

          ! uncertainty of start and end days    
          call fu_fade_in_linear( DayNumber(t2)/emg%startday_em(ix,iy), &
                                    emg%unc_startday/emg%startday_em(ix,iy), &
                                  fu_fade_in, &
                                        status )
          IF_NOTOK_RETURN(status=1)

          call fu_fade_out_linear( DayNumber(t2)/ (emg%startday_em(ix,iy)+emg%nday_em(ix,iy)) , &
                                    emg%unc_startday/(emg%startday_em(ix,iy)+emg%nday_em(ix,iy)), &
                                      fu_fade_out_days, &
                                        status )
          IF_NOTOK_RETURN(status=1)

        case ( 3 ) !Gamma with tails
        
          ! [-]    s       s/day     nday
          delta = nsec / 86400. / emg%nday_em(ix,iy) ! dt relative to season length
          season_ratio = (DayNumber(t2) - emg%startday_em(ix,iy)) / emg%nday_em(ix,iy)
          call fu_mk_new_pollen_gamma_w_tails(season_ratio, delta, gamma_f, status)
          IF_NOTOK_RETURN(status=1)

          ! fill gamma field
          emg%gamma_field(ix,iy) = gamma_f
      end select

      if ( emg%rest_avail_grains(ix,iy) > 0.0 ) then

          ! calculate other meteorological constraints, should hold for all pollen code versions
          constr = &
                 max(0., (1. - rain(ix,iy,1)/ emg%prec_th ) ) * &  ! preciptation threshold 
                 min( max( ( emg%high_hum_th - rh(ix,iy,1) ) / (emg%high_hum_th - emg%low_hum_th), 0.), 1.) * &  ! humidity thresholds
                 (emg%wspd_max_impact - exp( -1.0*(wsurf(ix,iy,1)+0.0) / emg%wspd_sat_lev ) )   ! windspeed constraint>


          select case ( polg_ver )  
            case ( 2 ) !Linear fade-in/out
              ! calculate ripened grains
              ripened_grains = & 
                     emg%max_avail_grains(ix,iy) / &     ! total annual grains
                      ( emg%nday_em(ix,iy) ) / 24.0  ! divide by number emitting days and hours per day

              ! uncertainty of total grains
              call fu_fade_out_linear( 1. - (emg%rest_avail_grains(ix,iy)/emg%max_avail_grains(ix,iy)), emg%unc_max_avail_grains , fu_fade_out, status)
              IF_NOTOK_RETURN(status=1)           

              emg%emis(ix,iy) = &
                min( emg%rest_avail_grains(ix,iy), &  ! resting available grains.            
                   ( ripened_grains + emg%ripened_left(ix,iy) )* & 
                   constr * &
                   fu_fade_in * &  ! Startday fade in
                   fu_fade_out_days * & ! Endday fade out
                   fu_fade_out &  ! total pollen fade out
                   ) 

            case ( 3 ) !Gamma with tails

              ! calculate ripened grains
              ripened_grains = & 
                     emg%max_avail_grains(ix,iy) * &     ! total annual grains
                      gamma_f !gamma distribution with tails
              emg%emis(ix,iy) = & !emission per hour
                min( emg%rest_avail_grains(ix,iy), (ripened_grains+emg%ripened_left(ix,iy))*constr )

          end select

          ! ripened grains left in catkins due to meteorological circumstances
          emg%ripened_left(ix,iy) = max(0.,  emg%ripened_left(ix,iy) + ripened_grains - emg%emis(ix,iy) )

      else
        emg%emis(ix,iy) = 0.0
      end if
    end do
    end do

    ! emg%emis is the number of grains emitted in the period t1 --> t2 (dt)
    ! so rest of the grains after t2 will be:
    emg%rest_avail_grains = emg%rest_avail_grains - emg%emis

    ! emitted grains per minute will be :   
    emg%emis = emg%emis/nmin

    ! ok
    status = 0

  end subroutine Pollen_Data_Grass_Get


  ! ***


  subroutine fu_fade_in_linear( value, uncer, fade_in, status )
    
    ! --- in/out -----
    real, intent(in)      ::  value
    real, intent(in)      ::  uncer
    real, intent(out)     ::  fade_in
    integer, intent(out)  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/fu_fade_in_linear'
    
    ! --- local -----
    real :: uncer_local
    
    if ( uncer < 1.e-5 ) then
      uncer_local = 1.e-5
    else
      uncer_local = uncer
    end if
    
    if ( value > 1. + uncer_local ) then
      fade_in = 1.
    else if ( value < 1. - uncer_local ) then
      fade_in = 0.
    else
      fade_in = ( value - 1. + uncer_local) / ( 2. * uncer_local)
    endif          
    
    ! ok
    status = 0
  
  end subroutine fu_fade_in_linear 
  
  
  ! ***
  
  
  subroutine fu_fade_out_linear( value, uncer, fade_out, status )
    
    ! --- in/out -----
    real, intent(in)      ::  value
    real, intent(in)      ::  uncer
    real, intent(out)     ::  fade_out
    integer, intent(out)  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/fu_fade_in_linear'
    
    ! --- local -----
    real :: uncer_local
    
    if ( uncer < 1.e-5 ) then
      uncer_local = 1.e-5
    else
      uncer_local = uncer
    end if
    
    if ( value > 1. + uncer_local ) then
      fade_out = 0.
    else if ( value < 1. - uncer_local ) then
      fade_out = 1.
    else
      fade_out = ( 1. + uncer_local - value) / ( 2. * uncer_local)
    endif          
    
    ! ok
    status = 0
  
  end subroutine fu_fade_out_linear 
   
      
  subroutine fu_mk_new_pollen_gamma_w_tails(fNowRel, delta_t, gamma_factor, status )

    ! FMI commentary:
    ! Returns the pollen prepared for release assuming the modified "taily" Gamma distribution
    ! of the season. Tails are the reason for many parameters: have to describe the main peak
    ! via gamma-type distribution, and both elevated tails via add-on corrections.
    ! formula: rate(x)=exp(-a_1/beta)* sum(scale_i * a_i^power_i), i=1:3
    ! where a_i = max(x-timesRel_i,0)
    !   
    ! Fitting for grass showed that:
    ! beta = 0.155, 
    ! nTerms = 3
    ! timesRel_1=0.164, scale1=13.1, power1=1.16, ! main season shape
    ! timesRel_2=-0.26, scale2=12.6, power2=2.8,   ! correction 1
    ! timesRel_3=-0.7, scale3=0.25, power3=1.3     ! correction 2
    !   
    implicit none
  
    ! Parameters
    integer, parameter               :: nTerms = 3 
    real,parameter                   :: beta = 0.155
    real,parameter,dimension(nTerms) :: timesRel = (/ 0.164, -0.26, -0.7 /)
    real,parameter,dimension(nTerms) :: scales = (/ 13.1, 12.6, 0.25 /)
    real,parameter,dimension(nTerms) :: powers = (/ 1.16, 2.8, 1.3 /)
    real, intent(in)                 :: fNowRel, delta_t ! normalised time, total scaling, timestep, beta
    real, intent(out)                :: gamma_factor
    integer, intent(out)             :: status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/fu_mk_new_pollen_gamma_w_tails'
         
  ! Local variables
    real :: a1, a, sum, time
    integer :: iTmp, jTmp

    a1 = max(0.0,fNowRel-timesRel(1))

    if(a1 > beta*10.)then
      gamma_factor = 0.0  ! too far from the season peak (as decided by timesRel 1)
      return
    else
      ! Be careful: the rise of the function can be quite steep.  First-order integration
      !
      sum = 0.0
      time = fNowRel
      do jTmp = 1,2
        do iTmp = nTerms, 1, -1
          if(time > timesRel(iTmp)) sum = sum + scales(iTmp) * (time - timesRel(iTmp))**powers(iTmp)
        end do
        sum = sum * exp(-a1 / beta)
        time = time + delta_t
      end do ! trapezoid integration cycle

      gamma_factor = 0.5 * sum * delta_t

    endif  ! if after the season

    ! ok
    status = 0

  end subroutine fu_mk_new_pollen_gamma_w_tails

end module Pollen_Data_Grass
