!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module Pollen_Data_Olive

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  use GO, only : TDate

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Pollen_Data_Olive

  public  ::  Pollen_Data_Olive_Init, Pollen_Data_Olive_Done
  public  ::  Pollen_Data_Olive_Get

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'Pollen_Olive_Data'
  

  ! --- local --------------------------------

  ! emission data base
  type T_Pollen_Data_Olive
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
    real, allocatable                 ::  heatsum_start(:,:) ! (nx,ny)
    real, allocatable                 ::  heatsum_end(:,:) ! (nx,ny)
    real, allocatable                 ::  heatsum(:,:)      ! (nx,ny)
    real, allocatable                 ::  heatsum_diff(:,:) ! (nx,ny)
    real, allocatable                 ::  heatsum_grow(:,:) ! (nx,ny)
        
    type(TDate)                       ::  startday_hs   ! start day of evaluating heat sum
    type(TDate)                       ::  endday_hs   ! end day of evaluating heat sum 
            
    real                              ::  temp_th  ! minimum temperature to get grow of heatsum
    character(len=8)                  ::  temp_th_unit

    real                              ::  low_hum_th, high_hum_th, prec_th  ! humidity and preciptation threshold
    character(len=8)                  ::  low_hum_th_unit, high_hum_th_unit, prec_th_unit

    real                              ::  wspd_sat_lev, wspd_max_impact     ! windspeed constraints
    character(len=8)                  ::  wspd_sat_lev_unit
    real                              ::  max_avail_grains_m2    ! available grains per m2
    
    real                              ::  unc_heatsum_start, unc_max_avail_grains
    character(len=8)                  ::  unc_heatsum_start_unit, unc_max_avail_grains_unit
    
    logical                           ::  allow_zero

  end type T_Pollen_Data_Olive



contains


  ! ===============================================================


  subroutine Pollen_Data_Olive_Init( emo, rcF, rckey, t, status )
    
    use GO , only : ReadRc
    use GO , only : Set
    use GO , only : goSplitString, goReplace, operator(<), operator(<=)
    use GO , only : TDate, Get
    use GO , only : goVarValue

    use LE_Grid         , only : ugg
    use C3PO            , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common  , only : Grid_Convertors
    use LE_Data         , only : LE_Data_Enable
    use LE_Restart_Tools, only : LE_Restart_Restore

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Olive), intent(inout)  ::  emo
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rckey
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Olive_Init'
    integer, parameter          ::  max_field = 100

    ! --- local ------------------------------

    character(len=128), allocatable :: descriptionlist(:) ! list of variables needed per pollen tracer
    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  values_in(:,:)
    character(len=64)     ::  units_in
    real                  ::  missing_value   

    integer               ::  idesc

    real                  ::  heatsum_diff_olives
    character(len=32)     ::  heatsum_diff_olives_unit
    
    integer               ::  nfield
    character(len=128)    ::  field
    integer               ::  vals(max_field)
    
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
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)

    ! read in input files for specific emissions
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.file', emo%file, status )
    IF_NOTOK_RETURN(status=1)
    
    ! year of emissions
    emo%year = t%year
    
    ! input descriptions
    allocate( descriptionlist(2) )
    descriptionlist = (/ 'standard_name=Olive Percentage', &
                         'standard_name=Heatsum_thresh  ' /)
    
    do idesc = 1, 2
      
      ! which description
      description = descriptionlist(idesc)
      
      ! open file:
      call file_in%Open( trim(emo%file), status )
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
      if ( trim(description) == 'standard_name=Olive Percentage' ) then
        ! check units:
        if ( trim(units_in) == '%' ) then
          values_in = values_in * 0.01
          allocate( emo%landuse_fracs(ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emo%landuse_fracs, status )
          IF_NOTOK_RETURN(status=1)    
          emo%landuse_fracs = max( emo%landuse_fracs, 0.0)
        else
          write(gol, '("unsupported unit: ", a, " for landuse percentage, file: ", a)' ) trim(units_in), trim(emo%file) ; call goErr
          TRACEBACK;status=1;return
        end if
      else if ( trim(description) == 'standard_name=Heatsum_thresh' ) then
        ! check units:
        if ( trim(units_in) == 'degree days' ) then
          allocate( emo%heatsum_start(ugg%nlon,ugg%nlat), stat=status ) 
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emo%heatsum_start, status )
          IF_NOTOK_RETURN(status=1)   
        else
          write(gol, '("unsupported unit: ", a, " for heatsum threshold, file: ", a)' ) trim(units_in), trim(emo%file) ; call goErr
          TRACEBACK;status=1;return
        end if 
      else 
        write( gol, '("Unsupprted description: ", a, " for file: " )' ) trim(description), trim(emo%file) ; call GoErr
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
    
    ! clear:
    deallocate( descriptionlist )
    
    ! Read constant heatsum difference for olive pollen
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.heatsum_diff', heatsum_diff_olives, status)
    IF_NOTOK_RETURN(status=1)       
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.heatsum_diff_unit', heatsum_diff_olives_unit, status)
    IF_NOTOK_RETURN(status=1)

    if ( .not. (heatsum_diff_olives_unit == 'degree_days') ) then
      write(gol, '("unsupported unit for heatsum diff pollen: ", a)' ) trim(heatsum_diff_olives_unit) ; call goErr
      TRACEBACK;status=1;return
    end if
          
    ! storage:
    allocate( emo%heatsum_diff(ugg%nlon,ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! Olive pollen has constant heatsum diff
    emo%heatsum_diff = heatsum_diff_olives

    ! storage:
    allocate( emo%heatsum_end(ugg%nlon,ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    emo%heatsum_end = emo%heatsum_start + emo%heatsum_end

    ! storage:
    allocate( emo%heatsum_grow(ugg%nlon,ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    emo%heatsum_grow = 0.0

    allocate( emo%intensity(ugg%nlon,ugg%nlat), stat=status  )
    IF_NOTOK_RETURN(status=1)
    ! fill in intensity array with only values 1        
    emo%intensity = 1.0

       
    ! initialize emissions
    allocate( emo%emis(ugg%nlon,ugg%nlat), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    allocate( emo%ripened_left(ugg%nlon,ugg%nlat), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! Read in other parameters
    
    ! temperature threshold to get grow of heatsum
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.temp_th', emo%temp_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.temp_th_unit', emo%temp_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check units  (LE unit is Kelvin)
    if ( emo%temp_th_unit == 'K' .or. emo%temp_th_unit == 'Kelvin' ) then
      emo%temp_th = emo%temp_th 
    else if ( emo%temp_th_unit == 'C' .or. emo%temp_th_unit == 'Celsius' ) then
      emo%temp_th = emo%temp_th + 273.15
    else
      write(gol, '( "wrong unit for temperature threshold: ", a )' ) emo%temp_th_unit
      TRACEBACK;status=1;return
    endif
    
    ! humiduty tresholds, below low_th:full emission, above high_th: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.low_hum_th', emo%low_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.low_hum_th_unit', emo%low_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.high_hum_th', emo%high_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.high_hum_th_unit', emo%high_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check humidty units
    if ( emo%low_hum_th_unit == '%' ) then
      emo%low_hum_th = emo%low_hum_th   ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emo%high_hum_th_unit
      TRACEBACK;status=1;return
    endif      
    if ( emo%high_hum_th_unit == '%' ) then
      emo%high_hum_th = emo%high_hum_th ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emo%high_hum_th_unit
      TRACEBACK;status=1;return
    endif
    
    ! preciptation threshold, above: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.prec_th', emo%prec_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.prec_th_unit', emo%prec_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! check units
    if ( emo%prec_th_unit == 'mm/hr' ) then
      ! LE unit is in m/s;
      !    m/s          mm/hr   /  mm/m  / (s/hr)
      emo%prec_th = emo%prec_th / 1000.0 / 3600.0
    else
      write(gol, '( "wrong unit for preciptation threshold: ", a )' ) emo%prec_th_unit
      TRACEBACK;status=1;return
    endif     
    
    ! windspeed constraints
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.wspd_sat_lev', emo%wspd_sat_lev, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.wspd_sat_lev_unit', emo%wspd_sat_lev_unit, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.wspd_max_impact', emo%wspd_max_impact, status )
    IF_ERROR_RETURN(status=1)    
    
    ! check units
    if ( emo%wspd_sat_lev_unit == 'm/s' ) then
      emo%wspd_sat_lev = emo%wspd_sat_lev  ! wsurf (u10) LE unit is in m/s  !!!! note that u and v are in km/min in LE
    else
      write(gol, '( "wrong unit for wind saturation level threshold: ", a )' ) emo%wspd_sat_lev_unit
      TRACEBACK;status=1;return
    endif    
    
    ! maximum available grains per m2
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.max_avail_grns', emo%max_avail_grains_m2, status, default=-999.0 )
    IF_ERROR_RETURN(status=1)    
    
    ! uncertainties in start of growing season and nr ov available grains
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.unc_heatsum_start', emo%unc_heatsum_start, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.unc_heatsum_start_unit', emo%unc_heatsum_start_unit, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.unc_max_avail_grns', emo%unc_max_avail_grains, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.unc_max_avail_grns_unit', emo%unc_max_avail_grains_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check uncertainty units
    if ( emo%unc_heatsum_start_unit == '%' ) then
      emo%unc_heatsum_start = 0.01 * emo%unc_heatsum_start   ! uncertainty in fraction
    else
      write(gol, '( "wrong unit for uncertainty of heatsum start: ", a )' ) emo%unc_heatsum_start_unit
      TRACEBACK;status=1;return
    endif      
    if ( emo%unc_max_avail_grains_unit == '%' ) then
      emo%unc_max_avail_grains = 0.01 * emo%unc_max_avail_grains ! uncertainty in fraction
    else
      write(gol, '( "wrong unit for uncertainty of max available grains: ", a )' ) emo%unc_max_avail_grains_unit
      TRACEBACK;status=1;return
    endif       
    
    
    ! Start and end of the accumulation period for the grow of the heatsum
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.start_accum_period', field, status )
    IF_ERROR_RETURN(status=1)    
    call goReplace( field, '%{year}', '(i4.4)', t%year, status )
    IF_NOTOK_RETURN(status=1)    
    call goSplitString( field, nfield, vals, status, sep='-')
    IF_NOTOK_RETURN(status=1)    
    call Set( emo%startday_hs, year=vals(1), month=vals(2), day=vals(3), calendar='gregorian'  )   
    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.end_accum_period', field, status )
    IF_ERROR_RETURN(status=1)    
    call goReplace( field, '%{year}', '(i4.4)', t%year, status )
    IF_NOTOK_RETURN(status=1)    
    call goSplitString( field, nfield, vals, status, sep='-')
    IF_NOTOK_RETURN(status=1)    
    call Set( emo%endday_hs  , year=vals(1), month=vals(2), day=vals(3), calendar='gregorian' )
        
    ! available grains per grid cell
    allocate(emo%max_avail_grains(ugg%nlon,ugg%nlat) )
    emo%max_avail_grains = emo%max_avail_grains_m2 * ugg%area * emo%landuse_fracs * emo%intensity
    
    allocate(emo%rest_avail_grains(ugg%nlon,ugg%nlat) )
    emo%rest_avail_grains = emo%max_avail_grains_m2 * ugg%area * emo%landuse_fracs * emo%intensity
    
    ! initialize heatsum
    allocate( emo%heatsum(ugg%nlon,ugg%nlat), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! check if t < start day, then heatsum = 0.0 else it should be read from restart file    
    call ReadRc( rcF, 'le.restart', do_restart, status )
    IF_NOTOK_RETURN(status=1)

    ! allow pollen zero for the rest of the year?
    call ReadRc( rcF, trim(rckey)//'.'//trim(emo%tracer)//'.allow_zero', emo%allow_zero, status)
    IF_NOTOK_RETURN(status=1)
    
    if ( t <= emo%startday_hs .and. (.not.do_restart) ) then
      ! info ..
      write (gol,'(" set heatsum for olive to zero ...")'); call goPr
      ! default:
      emo%heatsum = 0.0
    else
      
      ! No pollen in simulation for the rest of this year
      if ( emo%allow_zero ) then      
        write( gol, '(" Allow rest amount of Pollen to be zero for: ", a )' ) trim(emo%tracer) ; call GoPr
        emo%rest_avail_grains = 0.0
        emo%ripened_left      = 0.0
        
      else if ( .not. do_restart ) then
        write(gol, '(" Restart option is not enabled")' ) ; call GoErr
        write(gol, '(" model run starts after begin of growing season for Pollen")' )  ; call GoErr
        write(gol, '(" so information of current heatsum must be in restart file")' ) ; call GoErr
        write(gol, '(" Start day of growing season: ", i4, i2, i2)' ) emo%startday_hs%year, emo%startday_hs%month, emo%startday_hs%day ; call goErr
        write(gol, '(" Currentmodel time: ", i4, i2, i2)' ) t%year, t%month, t%day ; call goErr
       TRACEBACK; status=1; return
      
      else

        call ReadRc( rcF, 'le.restart.path', restart_path, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'le.restart.key', restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore heatsum for pollen emissions ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'heatsum_polo', emo%heatsum, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of pollen left for pollen emissions ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'amt_polo_left', emo%rest_avail_grains, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of ripened pollen left ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'ripened_polo_left', emo%ripened_left, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

      end if
    end if  
       
    ! ok
    status = 0

  end subroutine Pollen_Data_Olive_Init

  
  ! ***


  subroutine Pollen_Data_Olive_Done( emo, status )

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Olive), intent(inout)    ::  emo
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Olive_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! clear ?
    deallocate( emo%emis              )
    deallocate( emo%ripened_left      ) 
    deallocate( emo%landuse_fracs     ) 
    deallocate( emo%heatsum_start     ) 
    deallocate( emo%heatsum_end       ) 
    deallocate( emo%heatsum           ) 
    deallocate( emo%heatsum_diff      ) 
    deallocate( emo%heatsum_grow      )
    deallocate( emo%rest_avail_grains )

    !  ok
    status = 0

  end subroutine Pollen_Data_Olive_Done


  ! ***

  
  subroutine Pollen_Data_Olive_Get( emo, t1, t2, status )
    
    use GO, only : TDate, TIncrDate, operator(-), operator(<), operator(>)
    use GO, only : operator(<=), Midnight
    use dims, only : nx, ny
    use LE_Meteo_Data, only : day_mean_tsurf
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------

    type(T_Pollen_Data_Olive), intent(inout)    ::  emo
    type(TDate), intent(in)               ::  t1, t2   
    integer, intent(out)                  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Olive_Get'
    
    ! --- local ----------------------------
    
    integer                     ::  ix, iy
    type(TIncrDate)             ::  dt
    real                        ::  nmin
    real                        ::  ripened_grains, constr
    real                        ::  fu_fade_in, fu_fade_out, fu_fade_out_days
    
    ! meteo data:
    real, pointer               ::  rain (:,:,:)   ! (nlon,nlat,1)
    real, pointer               ::  wsurf(:,:,:)   ! (nlon,nlat,1)
    real, pointer               ::  rh   (:,:,:)   ! (nlon,nlat,nz)
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    
    ! --- begin ----------------------------
    
    ! meteo data:
    call LE_Data_GetPointer( 'rain'     , rain , status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'wspd_surf', wsurf, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh'       , rh   , status, check_units ='%' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'tsurf'    , tsurf, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    
    dt = t2-t1
    
    if ( t2 > emo%startday_hs .and. Midnight(t2) ) then

      do ix = 1, nx
      do iy = 1, ny
        if (day_mean_tsurf(ix,iy) - emo%temp_th > 0.0 ) then
            emo%heatsum(ix,iy) = emo%heatsum(ix,iy) + day_mean_tsurf(ix,iy) - emo%temp_th
        end if
        !! testing ...
        !if ( ix==1 .and. iy==3 ) then
        !  write (gol,*) 'xxx0 ', ix, iy, day_mean_tsurf(ix,iy), emo%temp_th, emo%heatsum(ix,iy), day_mean_tsurf(ix,iy); call goPr
        !  write (gol,*) '  x0 ', ix, iy, emo%heatsum(ix,iy); call goPr
        !end if
      end do
      end do

    else if ( t2 <= emo%startday_hs ) then
      emo%heatsum = 0.0
    endif

    ! update cumulative pollen release (R(t) in paper: Agricultural and Forest meteorology)
    do ix = 1, nx
    do iy = 1, ny 
      !
      !! testing ...
      !write (gol,*) 'xxx ', ix, iy, emo%heatsum(ix,iy), emo%heatsum_start(ix,iy), emo%unc_heatsum_start, &
      !                 emo%heatsum_diff(ix,iy), emo%rest_avail_grains(ix,iy); call goPr
      !
      if ( emo%heatsum(ix,iy) > emo%heatsum_start(ix,iy) * (1-emo%unc_heatsum_start) &
          .and. emo%heatsum_diff(ix,iy) > 0.0 &
          .and. emo%rest_avail_grains(ix,iy) > 0.0 ) then

          ! grow in heatsum per hour is
          emo%heatsum_grow(ix,iy) = max( 0.0, &
            (dt%day + ( 1/24.)*dt%hour + (1/(24.*60.))*dt%min + (1/(24.*60.*60.))*dt%sec ) * &
              ( tsurf(ix,iy,1) - emo%temp_th ) )

          ! calculate the fade in and fade out to deal with uncertainties of heatsum start and max avail grains          
          call fu_fade_in_linear( emo%heatsum(ix,iy)/emo%heatsum_start(ix,iy) , emo%unc_heatsum_start, fu_fade_in, status )
          IF_NOTOK_RETURN(status=1)

          call fu_fade_out_linear( 1. - (emo%rest_avail_grains(ix,iy)/emo%max_avail_grains(ix,iy)), emo%unc_max_avail_grains , fu_fade_out, status)
          IF_NOTOK_RETURN(status=1)           

          ! calculate heatsum dependent ripened grains
          ripened_grains = & 
                 emo%max_avail_grains(ix,iy) * &     ! total annual grains
                 emo%heatsum_grow(ix,iy) / ( emo%heatsum_diff(ix,iy) )

          ! calculate other meteorological constraints
          constr = &
                 max(0., (1. - rain(ix,iy,1)/ emo%prec_th ) ) * &  ! preciptation threshold 
                 min( max( ( emo%high_hum_th - rh(ix,iy,1) ) / (emo%high_hum_th - emo%low_hum_th), 0.), 1.) * &  ! humidity thresholds
                 (emo%wspd_max_impact - exp( -1.0*(wsurf(ix,iy,1)+0.0) / emo%wspd_sat_lev ) )   ! windspeed constraints

          emo%emis(ix,iy) = &
            min( emo%rest_avail_grains(ix,iy), &  ! resting available grains.            
                 ( ripened_grains + emo%ripened_left(ix,iy) )* & 
                 constr * &
                 fu_fade_in * &  ! Heatsum fade in
                 fu_fade_out &  ! total pollen fade out
                 ) 
          ! ripened grains left in catkins due to meteorological circumstances
          emo%ripened_left(ix,iy) = max(0.,  emo%ripened_left(ix,iy) + ripened_grains - emo%emis(ix,iy) )

      else
        emo%emis(ix,iy) = 0.0
      end if
    end do
    end do

    ! emo%emis is the number of grains emitted in the period t1 --> t2 (dt)
    ! so rest of the grains after t2 will be:
    emo%rest_avail_grains = emo%rest_avail_grains - emo%emis

    ! emitted grains per minute will be :   
    nmin = dt%day*24*60 + dt%hour*60 + dt%min + ( dt%sec / 60. )

    emo%emis = emo%emis/nmin

                                      
    ! ok
    status = 0

  end subroutine Pollen_Data_Olive_Get
  
  
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
  
  
end module Pollen_Data_Olive
