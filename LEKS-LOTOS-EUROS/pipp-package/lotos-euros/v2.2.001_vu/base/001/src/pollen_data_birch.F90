!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module Pollen_Data_Birch

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  use GO, only : TDate

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Pollen_Data_Birch

  public  ::  Pollen_Data_Birch_Init, Pollen_Data_Birch_Done
  public  ::  Pollen_Data_Birch_Get

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'Pollen_Data_Birch'
  

  ! --- local --------------------------------

  ! emission data base
  type T_Pollen_Data_Birch
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
    real, allocatable                 ::  scale_factor(:,:) ! (nx,ny)
    character(len=32)                 ::  emis_units
    character(len=32)                 ::  tracer
    character(len=32)                 ::  source
    integer                           ::  itr_glob  ! global tracer index

    character(len=512)                ::  file
    character(len=512)                ::  scalefile
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
  end type T_Pollen_Data_Birch



contains


  ! ===============================================================


  subroutine Pollen_Data_Birch_Init( emb, rcF, rckey, t, status )
    
    use GO , only : ReadRc
    use GO , only : Set
    use GO , only : goSplitString, goReplace, operator(<), operator(<=), operator(==)
    use GO , only : TDate, Get, NewDate
    use GO , only : goVarValue
    
    use LE_Grid         , only : ugg
    use LE_Grid         , only : dom
    use C3PO            , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common  , only : Grid_Convertors
    use LE_Data         , only : LE_Data_Enable
    use LE_Restart_Tools, only : LE_Restart_Restore

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Birch), intent(inout)  ::  emb
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rckey
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Birch_Init'
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

    logical               ::  sf_available
    integer               ::  nfield
    character(len=128)    ::  field
    integer               ::  vals(max_field)
    
    character(len=256)    ::  restart_path, restart_key
    
    ! restart vars
    logical               ::  do_restart
    integer               ::  glbo(2)
    
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
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.file', emb%file, status )
    IF_NOTOK_RETURN(status=1)
    
    ! in parallel, read in input files for scaling factors file
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.scale-factor.available', sf_available, status )
    IF_NOTOK_RETURN(status=1)

    ! year of emissions
    emb%year = t%year
    
    ! input descriptions
    allocate( descriptionlist(4), stat=status )
    IF_NOTOK_RETURN(status=1)
    descriptionlist = (/ 'standard_name=Birch Percentage   ', &
                         'standard_name=Heatsum_start      ', &
                         'standard_name=Heatsum_end        ', &
                         'standard_name=Intensity of Pollen' /)

    ! loop over input variables
    do idesc = 1, 4
      description = descriptionlist(idesc)

      call file_in%Open( trim(emb%file), status )
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
      if ( description == 'standard_name=Birch Percentage' ) then
        if ( trim(units_in) == '%' ) then
          values_in = values_in * 0.01
          allocate( emb%landuse_fracs(ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emb%landuse_fracs, status )
          IF_NOTOK_RETURN(status=1)
        else
          write(gol, '("unsupported unit: ", a, ", in Birch percentage, file: ", a)' ) trim(units_in), trim(emb%file) ; call goErr
          TRACEBACK;status=1;return
        end if
      else if ( description == 'standard_name=Heatsum_start' ) then
        if ( trim(units_in) == '1' ) then
          allocate( emb%heatsum_start(ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emb%heatsum_start, status )
          IF_NOTOK_RETURN(status=1)
        else
          write(gol, '("unsupported unit: ", a, ", in heatsum start, file: ", a)' ) trim(units_in), trim(emb%file) ; call goErr
          TRACEBACK;status=1;return
        end if
      else if ( description == 'standard_name=Heatsum_end' ) then
        if ( trim(units_in) == '1' ) then
          allocate( emb%heatsum_end(ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emb%heatsum_end, status )
          IF_NOTOK_RETURN(status=1)
        else
          write(gol, '("unsupported unit: ", a, ", in heatsum end, file: ", a)' ) trim(units_in), trim(emb%file) ; call goErr
          TRACEBACK;status=1;return
        end if
      else if ( description == 'standard_name=Intensity of Pollen' ) then
        if ( trim(units_in) == '1' ) then
          allocate( emb%intensity(ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, emb%intensity, status )
          IF_NOTOK_RETURN(status=1)
        else
          write(gol, '("unsupported unit: ", a, ", in intensity, file: ", a)' ) trim(units_in), trim(emb%file) ; call goErr
          TRACEBACK;status=1;return
        end if
      else 
        write( gol, '("Unsupprted description: ", a, " for file: " )' ) trim(description), trim(emb%file) ; call GoErr
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
    
    ! If scaling factor file exists
    if ( sf_available ) then 
      !info 
      write(gol, '(" Using birch pollen scaling factors...")' ) ; call goPr

      call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.scale-factor.file', emb%scalefile, status )
      IF_NOTOK_RETURN(status=1)

      call file_in%Open( trim(emb%scalefile), status )
      IF_NOTOK_RETURN(status=1)
      
      description='var_name=scale_factor;units=1'
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

      ! check units
      if ( trim(units_in) == '1' ) then
        allocate( emb%scale_factor(ugg%nlon,ugg%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
        call Grid_Convertors%Ugg_Interpol( grid_in, values_in, ugg, emb%scale_factor, status )
        IF_NOTOK_RETURN(status=1)  
      else
        write(gol, '("unsupported unit for pollen scaling factor, file: ", a)' ) trim(emb%scalefile) ; call goErr
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
      
    endif ! end with scale factors file
        
    ! Other heatsum arrays
    allocate( emb%heatsum_diff(ugg%nlon,ugg%nlat) )
    emb%heatsum_diff = emb%heatsum_end - emb%heatsum_start

    allocate( emb%heatsum_grow(ugg%nlon,ugg%nlat) )
    emb%heatsum_grow = 0.0
   
    ! initialize emissions
    allocate(emb%emis(ugg%nlon,ugg%nlat) )
    emb%emis = 0.0   
    
    allocate(emb%ripened_left(ugg%nlon,ugg%nlat) )
    emb%ripened_left = 0.0
    
    ! Read in other parameters
    
    ! temperature threshold to get grow of heatsum
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.temp_th', emb%temp_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.temp_th_unit', emb%temp_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check units  (LE unit is Kelvin)
    if ( emb%temp_th_unit == 'K' .or. emb%temp_th_unit == 'Kelvin' ) then
      emb%temp_th = emb%temp_th 
    else if ( emb%temp_th_unit == 'C' .or. emb%temp_th_unit == 'Celsius' ) then
      emb%temp_th = emb%temp_th + 273.15
    else
      write(gol, '( "wrong unit for temperature threshold: ", a )' ) emb%temp_th_unit
      TRACEBACK;status=1;return
    endif
    
    ! humiduty tresholds, below low_th:full emission, above high_th: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.low_hum_th', emb%low_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.low_hum_th_unit', emb%low_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.high_hum_th', emb%high_hum_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.high_hum_th_unit', emb%high_hum_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check humidty units
    if ( emb%low_hum_th_unit == '%' ) then
      emb%low_hum_th = emb%low_hum_th   ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emb%high_hum_th_unit
      TRACEBACK;status=1;return
    endif      
    if ( emb%high_hum_th_unit == '%' ) then
      emb%high_hum_th = emb%high_hum_th ! LE unit is %
    else
      write(gol, '( "wrong unit for humidity threshold: ", a )' ) emb%high_hum_th_unit
      TRACEBACK;status=1;return
    endif
    
    ! preciptation threshold, above: no emission
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.prec_th', emb%prec_th, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.prec_th_unit', emb%prec_th_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! check units
    if ( emb%prec_th_unit == 'mm/hr' ) then
      ! LE unit is in m/s;
      !    m/s          mm/hr   /  mm/m  / (s/hr)
      emb%prec_th = emb%prec_th / 1000.0 / 3600.0
    else
      write(gol, '( "wrong unit for preciptation threshold: ", a )' ) emb%prec_th_unit
      TRACEBACK;status=1;return
    endif     
    
    ! windspeed constraints
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.wspd_sat_lev', emb%wspd_sat_lev, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.wspd_sat_lev_unit', emb%wspd_sat_lev_unit, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.wspd_max_impact', emb%wspd_max_impact, status )
    IF_ERROR_RETURN(status=1)    
    
    ! check units
    if ( emb%wspd_sat_lev_unit == 'm/s' ) then
      emb%wspd_sat_lev = emb%wspd_sat_lev  ! wsurf (u10) LE unit is in m/s  !!!! note that u and v are in km/min in LE
    else
      write(gol, '( "wrong unit for wind saturation level threshold: ", a )' ) emb%wspd_sat_lev_unit
      TRACEBACK;status=1;return
    endif    
    
    ! maximum available grains per m2
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.max_avail_grns', emb%max_avail_grains_m2, status )
    IF_ERROR_RETURN(status=1)    
    
    ! uncertainties in start of growing season and nr ov available grains
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.unc_heatsum_start', emb%unc_heatsum_start, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.unc_heatsum_start_unit', emb%unc_heatsum_start_unit, status )
    IF_ERROR_RETURN(status=1)    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.unc_max_avail_grns', emb%unc_max_avail_grains, status )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.unc_max_avail_grns_unit', emb%unc_max_avail_grains_unit, status )
    IF_ERROR_RETURN(status=1)
    
    ! Check uncertainty units
    if ( emb%unc_heatsum_start_unit == '%' ) then
      emb%unc_heatsum_start = 0.01 * emb%unc_heatsum_start   ! uncertainty in fraction
    else
      write(gol, '( "wrong unit for uncertainty of heatsum start: ", a )' ) emb%unc_heatsum_start_unit
      TRACEBACK;status=1;return
    endif      
    if ( emb%unc_max_avail_grains_unit == '%' ) then
      emb%unc_max_avail_grains = 0.01 * emb%unc_max_avail_grains ! uncertainty in fraction
    else
      write(gol, '( "wrong unit for uncertainty of max available grains: ", a )' ) emb%unc_max_avail_grains_unit
      TRACEBACK;status=1;return
    endif       
    
    
    ! Start and end of the accumulation period for the grow of the heatsum
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.start_accum_period', field, status )
    IF_ERROR_RETURN(status=1)    
    call goReplace( field, '%{year}', '(i4.4)', t%year, status )
    IF_NOTOK_RETURN(status=1)    
    call goSplitString( field, nfield, vals, status, sep='-')
    IF_NOTOK_RETURN(status=1)    
    call Set( emb%startday_hs, year=vals(1), month=vals(2), day=vals(3), calendar='gregorian'  )   
    
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.end_accum_period', field, status )
    IF_ERROR_RETURN(status=1)    
    call goReplace( field, '%{year}', '(i4.4)', t%year, status )
    IF_NOTOK_RETURN(status=1)    
    call goSplitString( field, nfield, vals, status, sep='-')
    IF_NOTOK_RETURN(status=1)    
    call Set( emb%endday_hs  , year=vals(1), month=vals(2), day=vals(3), calendar='gregorian' )
        
    ! available grains per grid cell
    allocate(emb%max_avail_grains(ugg%nlon,ugg%nlat) )
    if ( sf_available ) then 
      emb%max_avail_grains = emb%max_avail_grains_m2 * ugg%area * emb%landuse_fracs * emb%intensity &
                           * emb%scale_factor 
      deallocate( emb%scale_factor     ) 
    else
      emb%max_avail_grains = emb%max_avail_grains_m2 * ugg%area * emb%landuse_fracs * emb%intensity
    endif
    
    allocate(emb%rest_avail_grains(ugg%nlon,ugg%nlat) )
    emb%rest_avail_grains = emb%max_avail_grains_m2 * ugg%area * emb%landuse_fracs * emb%intensity
    
    ! initialize heatsum
    allocate(emb%heatsum(ugg%nlon,ugg%nlat) )    
    
    ! check if t < start day, then heatsum = 0.0 else it should be read from restart file    
    call ReadRc( rcF, 'le.restart', do_restart, status )
    IF_NOTOK_RETURN(status=1)
    
    ! allow pollen zero for the rest of the year?
    call ReadRc( rcF, trim(rckey)//'.'//trim(emb%tracer)//'.allow_zero', emb%allow_zero, status)
    IF_NOTOK_RETURN(status=1)
    
    if ( t <= emb%startday_hs .and. (.not.do_restart) ) then
      emb%heatsum = 0.0

    ! Before/at the start of the season when restart file is used
    ! Set maximum value for available pollen and reset to zero for ripened pollen
    else if ( t <= emb%startday_hs .and. do_restart ) then
      !emb%heatsum = 0.0 
      emb%rest_avail_grains = emb%max_avail_grains_m2 * ugg%area * emb%landuse_fracs * emb%intensity
      emb%ripened_left = 0.0
      write( gol, '(" Run starts before or at start of birch pollen season with a &
      restart file: Set max value for avail. pollen and reset to 0.0 &
      for ripened pollen.")') ; call GoPr

    !After the start of the season with or without restart file:
    else
      
      ! No pollen in simulation for the rest of this year
      ! regardless if restart file is used
      if ( emb%allow_zero ) then      
        write( gol, '(" Allow rest amount of Pollen to be zero for: ", a )' ) trim(emb%tracer) ; call GoPr
        emb%rest_avail_grains = 0.0
        emb%ripened_left = 0.0
      !~ a restart file is required
      else if ( .not. do_restart ) then
        write(gol, '(" Restart option is not enabled")' ) ; call GoErr
        write(gol, '(" model run starts after begin of growing season for Pollen")' )  ; call GoErr
        write(gol, '(" so information of current heatsum must be in restart file")' ) ; call GoErr
        write(gol, '(" Start day of growing season: ", i4, i2, i2)' ) emb%startday_hs%year, emb%startday_hs%month, emb%startday_hs%day ; call goErr
        write(gol, '(" Currentmodel time: ", i4, i2, i2)' ) t%year, t%month, t%day ; call goErr
        TRACEBACK; status=1; return
      
      ! Restart file exist, read heatsum, available and ripened pollen from
      ! there
      else 
      
        ! domain start in global index space:
        call dom%Get( status, glbo=glbo )
        IF_NOTOK_RETURN(status=1)

        call ReadRc( rcF, 'le.restart.path', restart_path, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'le.restart.key', restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore heatsum for pollen emissions ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'heatsum_polb', emb%heatsum, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of pollen left ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'amt_polb_left', emb%rest_avail_grains, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'("LE:       restore amount of ripened pollen left ...")'); call goPr
        ! read from restart file:
        call LE_Restart_Restore( 'ripened_polb_left', emb%ripened_left, &
                                   t, restart_path, restart_key, status )
        IF_NOTOK_RETURN(status=1)

      end if
    end if  
       
    ! ok
    status = 0

  end subroutine Pollen_Data_Birch_Init

  
  ! ***


  subroutine Pollen_Data_Birch_Done( emb, status )

    ! --- in/out ------------------------------

    type(T_Pollen_Data_Birch), intent(inout)    ::  emb
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! clear ?
    deallocate( emb%emis              )
    deallocate( emb%ripened_left      ) 
    deallocate( emb%landuse_fracs     ) 
    deallocate( emb%heatsum_start     ) 
    deallocate( emb%heatsum_end       ) 
    deallocate( emb%heatsum           ) 
    deallocate( emb%heatsum_diff      ) 
    deallocate( emb%heatsum_grow      )
    deallocate( emb%rest_avail_grains )

    !  ok
    status = 0

  end subroutine Pollen_Data_Birch_Done


  ! ***

  
  subroutine Pollen_Data_Birch_Get( emb, t1, t2, status )
    
    use GO, only : TDate, TIncrDate, operator(-), operator(<), operator(>)
    use GO, only : operator(<=), Midnight, DayNumber
    use dims, only : nx, ny
    use LE_Meteo_Data, only : day_mean_tsurf
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------

    type(T_Pollen_Data_Birch), intent(inout)    ::  emb
    type(TDate), intent(in)               ::  t1, t2   
    integer, intent(out)                  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Birch_Get'
    
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
    
    ! timestep:
    dt = t2-t1
    
    ! update heatsum
    if ( t2 > emb%startday_hs .and. Midnight(t2) ) then

      do ix = 1, nx
      do iy = 1, ny
        if (day_mean_tsurf(ix,iy) - emb%temp_th > 0.0 ) then
            emb%heatsum(ix,iy) = emb%heatsum(ix,iy) + day_mean_tsurf(ix,iy) - emb%temp_th
        end if
      end do
      end do

    else if ( t2 <= emb%startday_hs ) then
      emb%heatsum = 0.0
    endif

    ! update cumulative pollen release (R(t) in paper: Agricultural and Forest meteorology)
    do ix = 1, nx
    do iy = 1, ny 
      if ( emb%heatsum(ix,iy) > emb%heatsum_start(ix,iy) * (1-emb%unc_heatsum_start) &
          .and. emb%heatsum_diff(ix,iy) > 0.0 &
          .and. emb%rest_avail_grains(ix,iy) > 0.0 ) then

          ! grow in heatsum per hour is
          emb%heatsum_grow(ix,iy) = max( 0.0, &
            (dt%day + ( 1/24.)*dt%hour + (1/(24.*60.))*dt%min + (1/(24.*60.*60.))*dt%sec ) * &
              ( tsurf(ix,iy,1) - emb%temp_th ) )

          ! calculate the fade in and fade out to deal with uncertainties of heatsum start and max avail grains          
          call fu_fade_in_linear( emb%heatsum(ix,iy)/emb%heatsum_start(ix,iy) , emb%unc_heatsum_start, fu_fade_in, status )
          IF_NOTOK_RETURN(status=1)

          call fu_fade_out_linear( 1. - (emb%rest_avail_grains(ix,iy)/emb%max_avail_grains(ix,iy)), emb%unc_max_avail_grains , fu_fade_out, status)
          IF_NOTOK_RETURN(status=1)           

          ! calculate heatsum dependent ripened grains
          ripened_grains = & 
                 emb%max_avail_grains(ix,iy) * &     ! total annual grains
                 emb%heatsum_grow(ix,iy) / ( emb%heatsum_diff(ix,iy) )

          ! calculate other meteorological constraints
          constr = &
                 max(0., (1. - rain(ix,iy,1)/ emb%prec_th ) ) * &  ! preciptation threshold 
                 min( max( ( emb%high_hum_th - rh(ix,iy,1) ) / (emb%high_hum_th - emb%low_hum_th), 0.), 1.) * &  ! humidity thresholds
                 (emb%wspd_max_impact - exp( -1.0*(wsurf(ix,iy,1)+0.0) / emb%wspd_sat_lev ) )   ! windspeed constraints

          emb%emis(ix,iy) = &
            min( emb%rest_avail_grains(ix,iy), &  ! resting available grains.            
                 ( ripened_grains + emb%ripened_left(ix,iy) )* & 
                 constr * &
                 fu_fade_in * &  ! Heatsum fade in
                 fu_fade_out &  ! total pollen fade out
                 ) 
          ! ripened grains left in catkins due to meteorological circumstances
          emb%ripened_left(ix,iy) = max(0.,  emb%ripened_left(ix,iy) + ripened_grains - emb%emis(ix,iy) )

      else
        emb%emis(ix,iy) = 0.0
      end if
    end do
    end do

    ! emb%emis is the number of grains emitted in the period t1 --> t2 (dt)
    ! so rest of the grains after t2 will be:
    emb%rest_avail_grains = emb%rest_avail_grains - emb%emis

    ! emitted grains per minute will be :   
    nmin = dt%day*24*60 + dt%hour*60 + dt%min + ( dt%sec / 60. )

    emb%emis = emb%emis/nmin
      
    ! ok
    status = 0

  end subroutine Pollen_Data_Birch_Get
  
  
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
  
  
end module Pollen_Data_Birch
