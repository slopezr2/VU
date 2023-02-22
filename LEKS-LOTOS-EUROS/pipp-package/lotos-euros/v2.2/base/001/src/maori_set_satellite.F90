!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_Set_Satellite

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate
!  use GO, only : T_ISTA

  use Grid, only : TllGridInfo

  use MAORI_Param, only : MAORI_LEN_NAME, MAORI_LEN_FILE, MAORI_LEN_LINE
!  use MAORI_LocList, only : T_MAORI_LocList
  use MAORI_VarList, only : T_MAORI_VarList

!  use EMEP_DataFile, only : T_EMEP_DataFile
  
#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  implicit none

  ! --- in/out ---------------------------

  private

  public  ::  T_MAORI_Data_Satellite
  public  ::  MAORI_Data_Satellite_Init, MAORI_Data_Satellite_Done
  public  ::  MAORI_Data_Satellite_Start, MAORI_Data_Satellite_Setup
  public  ::  MAORI_Data_Satellite_Inq
  public  ::  MAORI_Data_Satellite_Param_Inq, MAORI_Data_Satellite_Param_Put
  public  ::  MAORI_Data_Satellite_Obs_Get, MAORI_Data_Satellite_Obs_Put
  public  ::  MAORI_Data_Satellite_Var_Inq
  public  ::  MAORI_Data_Satellite_Loc_Inq

  public  ::  T_MAORI_State_Satellite
  public  ::  MAORI_State_Satellite_Init, MAORI_State_Satellite_Done
  public  ::  MAORI_State_Satellite_Start, MAORI_State_Satellite_Setup
  public  ::  MAORI_State_Satellite_Put
  public  ::  MAORI_State_Satellite_Post
  public  ::  MAORI_State_Satellite_Values_Get, MAORI_State_Satellite_Values_Put
  public  ::  MAORI_State_Satellite_LocValues_Get, MAORI_State_Satellite_LocValues_Put
  public  ::  MAORI_State_Satellite_Obs_Get

  public  ::  T_MAORI_Output_Satellite
  public  ::  MAORI_Output_Satellite_Init, MAORI_Output_Satellite_Done
  public  ::  MAORI_Output_Satellite_Start
  public  ::  MAORI_Output_Satellite_Write


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_Set_Satellite'

  ! parameter indices:
  integer, parameter  ::  ipar_bound_west      = 1
  integer, parameter  ::  ipar_bound_east      = 2
  integer, parameter  ::  ipar_bound_south     = 3
  integer, parameter  ::  ipar_bound_north     = 4
  integer, parameter  ::  ipar_grid_nlon       = 5
  integer, parameter  ::  ipar_grid_nlat       = 6
  integer, parameter  ::  ipar_nlayer          = 7
  integer, parameter  ::  ipar_bias_correction = 8
  integer, parameter  ::  npar = 8
  
  ! maximum number of satellite platforms:
  integer, parameter    ::  maxplatform = 2


  ! --- types ------------------------------

  type T_MAORI_Data_Satellite
    ! name:
    character(len=MAORI_LEN_NAME)         ::  name
    ! flags:
    logical                               ::  started
    ! instrument:
    character(len=MAORI_LEN_NAME)         ::  instrument
    ! archive description:
    character(len=MAORI_LEN_FILE)         ::  query
    !
    ! parameters retrieved from model
    !
    ! parameter flags:
    logical                               ::  param_set(npar)
    ! model parameters:
    real                                  ::  bound_west, bound_east, bound_south, bound_north
    integer                               ::  grid_nlon, grid_nlat
    type(TllGridInfo)                     ::  grid_lli
    character(len=MAORI_LEN_LINE)         ::  bias_correction
    ! time resolution:
    character(len=5)                      ::  tres_period    ! 'day'  (in future: 'month', ...)
    !character(len=4)                      ::  tres_type      ! 'aver' (in future: 'inst', ..)
    !character(len=6)                      ::  tres_unit       ! 'hour'
    !real                                  ::  tres_start
    !real                                  ::  tres_step
    !type(TIncrDate)                       ::  tres_dt
    !logical                               ::  tres_hist
    ! time step stuff:
    logical                               ::  the_end
    type(TDate)                           ::  setup_tr(2)
    type(TDate)                           ::  period_tr(2)
    logical                               ::  end_of_period
    logical                               ::  new_period
    !integer                               ::  nstep
    !integer                               ::  istep
    ! ~ history
    !integer                               ::  nhist
    !integer                               ::  ihist
    !! ~ averaging
    !type(T_ISTA)                          ::  ista
    !real                                  ::  w_aver, w_prev, w_curr
    !logical                               ::  finished_aver
    ! horizontal grid:
    integer                               ::  grid_nx, grid_ny
    ! vertical resolution:
    logical                               ::  profile
    integer                               ::  nlayer   ! in the model
    integer                               ::  maxlev   ! 1 or nlayer+1
    ! records (either pixels or 1)
    integer                               ::  pix_max
    integer                               ::  pix_n
    ! list with variables:
    type(T_MAORI_VarList)                 ::  varlist
    integer                               ::  nvar
    integer, pointer                      ::  nlev(:)  ! nvar ; 1, nlayer, or nlayer+1
    !
    ! observations
    !
    ! list with input variables:
    type(T_MAORI_VarList)                 ::  obs_varlist
    integer                               ::  obs_nvar
    integer                               ::  obs_nlev
    !character(len=MAORI_LEN_LINE)         ::  obs_query
    !character(len=MAORI_LEN_FILE)         ::  obs_file
    !logical                               ::  obs_opened
    !type(T_EMEP_DataFile), pointer        ::  obs_emep
    real, pointer                         ::  obs_pix_value (:,:,:) ! obs_nlev,maxpix,obs_nvar
    real, pointer                         ::  obs_pix_errstd(:,:,:)
    integer, pointer                      ::  obs_pix_astat(:,:,:)
    real, pointer                         ::  obs_grid_value (:,:,:,:) ! nlon,nlat,obs_nlev,obs_nvar
    real, pointer                         ::  obs_grid_errstd(:,:,:,:)
    integer, pointer                      ::  obs_grid_astat(:,:,:,:)
    !real                                  ::  obs_err_frac
    !integer, pointer                      ::  obs_iloc2iobs(:)
    integer, pointer                      ::  obs_ivar(:)
    !type(TDate)                           ::  obs_t
    !
    ! assimilation stuff
    !
    ! error std.dev. fraction:
    real                                  ::  sigma_offset
    real                                  ::  sigma_scale
    !
    ! analyse this set ?
    logical                               ::  assim_analyse
    ! localisation length:
    real                                  ::  assim_rho
    ! screening factor: measurements are rejected
    ! if square of departure exceeds factor times variance:
    real                                  ::  assim_alfa
    ! number of 'state' values:
    integer                               ::  nvalue
    integer                               ::  nlocvalue
    !
    ! storage
    !
    ! how stored ?
    logical                               ::  store_pix
    logical                               ::  store_grid
    !
    !integer                               ::  pix_ipl(maxp)
    !integer                               ::  pix_jday(maxp)
    real, pointer                         ::  pix_time(:)
    integer, pointer                      ::  pix_date(:,:)
    real, pointer                         ::  pix_longitude(:)
    real, pointer                         ::  pix_latitude (:)
    real, pointer                         ::  pix_proba    (:)
    integer, pointer                      ::  pix_slots(:)
    real, pointer                         ::  pix_grid_i(:)
    real, pointer                         ::  pix_grid_j(:)
    !real                                  ::  pix_Cloud_Fraction_Land(maxp)
    !real                                  ::  pix_Cloud_Fraction_Ocean(maxp)
    !real                                  ::  pix_Optical_Depth_Land_And_Ocean(maxp)
    !real                                  ::  pix_Optical_Depth_Land_And_Ocean_sigma(maxp)
    !real, pointer                         ::  pix_obs(:)
    !real, pointer                         ::  pix_obs_sigma(:)
    !character(len=len_units)              ::  pix_Cloud_Fraction_Land__units
    !character(len=len_units)              ::  pix_Cloud_Fraction_Ocean__units
    !character(len=len_units)              ::  pix_Optical_Depth_Land_And_Ocean__units
    !character(len=len_long_name)          ::  pix_Cloud_Fraction_Land__long_name
    !character(len=len_long_name)          ::  pix_Cloud_Fraction_Ocean__long_name
    !character(len=len_long_name)          ::  pix_Optical_Depth_Land_And_Ocean__long_name
    !! cell indices in model grid:
    !integer                               ::  pix_ix(maxp)
    !integer                               ::  pix_iy(maxp)
    !! assimilation flag:
    !integer                               ::  pix_assim_status(maxp)
    !
    ! stored grid:
    ! ~ pixels are averaged over some hours (differs per pixel),
    !   and also over a grid cell; 
    !   pixels per grid cell, and weights of instantaneous hourly concentration:
    integer, pointer                      ::  grid_npix(:,:)    ! nx,ny
    real, pointer                         ::  grid_w24(:,:,:)  ! nx,ny,0:24
!    !
!    ! o modis
!    ! 
!    character(len=MAORI_LEN_FILE)         ::  modis_path
!    integer                               ::  modis_nplatform
!    character(len=MAORI_LEN_NAME)         ::  modis_platforms(maxplatform)
!    character(len=MAORI_LEN_NAME)         ::  modis_products (maxplatform)
!    character(len=MAORI_LEN_FILE)         ::  modis_list_file   (maxplatform)
!    logical                               ::  modis_list_opened (maxplatform)
!    integer                               ::  modis_list_fu     (maxplatform)
!    integer                               ::  modis_list_iline  (maxplatform)
!    character(len=MAORI_LEN_FILE)         ::  modis_list_dfile  (maxplatform)
!    type(TDate)                           ::  modis_list_t      (maxplatform)
!    integer                               ::  modis_list_jday   (maxplatform)
!    !
!    ! o seviri
!    !
!    character(len=MAORI_LEN_FILE)         ::  seviri_path
!    character(len=MAORI_LEN_NAME)         ::  seviri_prefix
!    character(len=MAORI_LEN_FILE)         ::  seviri_basename
  end type T_MAORI_Data_Satellite
  

  ! simulation from a state:
  type T_MAORI_State_Satellite
    ! id
    character(len=MAORI_LEN_NAME)         ::  name
    ! flags:
    logical                               ::  started
    ! simulated values:
    logical, pointer                      ::  pix_set (:,:,:)    ! maxlev,maxpix,nvar
    real, pointer                         ::  pix_curr(:,:,:)    ! maxlev,maxpix,nvar
    ! simulated values:
    logical, pointer                      ::  grid_set (:,:,:,:)    ! grid_nx,grid_ny,maxlev,nvar
    real, pointer                         ::  grid_curr(:,:,:,:)    ! grid_nx,grid_ny,maxlev,nvar
    real, pointer                         ::  grid_prev(:,:,:,:)    ! grid_nx,grid_ny,maxlev,nvar
    real, pointer                         ::  grid_aver(:,:,:,:)    ! grid_nx,grid_ny,maxlev,nvar
    ! adhoc: for testing only ..
    real, pointer                         ::  grid_all(:,:,:,:,:)    ! grid_nx,grid_ny,maxlev,0:24,nvar
  end type T_MAORI_State_Satellite
  

  ! output: file units etc
  type T_MAORI_Output_Satellite
    ! id
    character(len=MAORI_LEN_NAME)         ::  name
    ! file opened ?
    logical                               ::  opened
    ! record counters:
    integer                               ::  pix_irec
#ifdef with_netcdf
    ! file handle:
    integer                               ::  ncid
    ! dimension handles:
    !integer                               ::  dimid_namelen
    !integer                               ::  dimid_codelen
    integer                               ::  dimid_grid_nx
    integer                               ::  dimid_grid_ny
    integer                               ::  dimid_column
    integer                               ::  dimid_pix
    integer                               ::  dimid_hhour
    !integer                               ::  dimid_slev
    !integer                               ::  dimid_flev
    !integer                               ::  dimid_hlev
    !integer                               ::  dimid_time
    integer                               ::  dimid_datelen
    ! dimension variables:
    integer                               ::  varid_pix_time !, varid_pix_time2
    integer                               ::  varid_pix_date !, varid_pix_date2
    integer                               ::  varid_pix_lon
    integer                               ::  varid_pix_lat
    integer                               ::  varid_pix_proba
    integer                               ::  varid_pix_slots
    integer                               ::  varid_pix_grid_i
    integer                               ::  varid_pix_grid_j
    !! station variables:
    !integer                               ::  varid_stn_name
    !integer                               ::  varid_stn_code
    integer                               ::  varid_grid_lon
    integer                               ::  varid_grid_lat
    !integer                               ::  varid_stn_alt
    integer                               ::  varid_grid_n
    integer                               ::  varid_grid_w
    ! simulated variables:
    integer, pointer                      ::  varid_pix (:)
    integer, pointer                      ::  varid_grid(:)
    integer, pointer                      ::  varid_grid_all(:)
    ! observations:
    integer, pointer                      ::  varid_pix_y(:)
    integer, pointer                      ::  varid_pix_r(:)
    integer, pointer                      ::  varid_pix_astat(:)
    integer, pointer                      ::  varid_grid_y(:)
    integer, pointer                      ::  varid_grid_r(:)
    integer, pointer                      ::  varid_grid_astat(:)
#endif
  end type T_MAORI_Output_Satellite


  ! --- interfaces -------------------------------
  
  interface MAORI_Data_Satellite_Param_Put
    module procedure MAORI_Data_Satellite_Param_Put_i
    module procedure MAORI_Data_Satellite_Param_Put_r
    module procedure MAORI_Data_Satellite_Param_Put_s
  end interface


contains


  ! ********************************************************************
  ! ***
  ! *** data
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Satellite_Init( mdata, name, rcF, rckey, t0, status )

    use GO, only : TRcFile, ReadRc
    use GO, only : goVarValue, goSplitString
    use GO, only : TDate, NewDate, AnyDate !, IncrDate, operator(*)
    use GO, only : goGetFU
    !use GO, only : ISTA_Init

    use MAORI_Param  , only : MAORI_LEN_LINE
    use MAORI_VarList, only : MAORI_VarList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(out)     ::  mdata
    character(len=*), intent(in)                  ::  name
    type(TRcFile), intent(in)                     ::  rcF
    character(len=*), intent(in)                  ::  rckey
    type(TDate), intent(in)                       ::  t0
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Init'

    ! --- local -----------------------------------

    character(len=MAORI_LEN_LINE)       ::  key
    integer                             ::  n
    integer                             ::  ipl
    logical                             ::  exist

    ! --- begin -----------------------------------

    ! store name:
    mdata%name = trim(name)
    
    ! by default, store nothing:
    mdata%store_pix  = .false.
    mdata%store_grid = .false.

    ! read key that describes archive:
    call ReadRc( rcF, trim(rckey)//'.query', mdata%query, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read key that describes instrument:
    call ReadRc( rcF, trim(rckey)//'.instrument', mdata%instrument, status )
    IF_NOTOK_RETURN(status=1)
    
    ! instrument dependend stuff:
    select case ( trim(mdata%instrument) )
    
!      ! ~ modis:
!      case ( 'modis' )
!      
!        ! store pixels:
!        mdata%store_pix = .true.
!
!        ! path to archive:
!        mdata%modis_path = 'no-path'
!          call goVarValue( mdata%query, ';', 'path', '=', mdata%modis_path, status )
!          IF_ERROR_RETURN(status=1)
!
!        ! space seperated list with platform names:
!        key = 'no-platforms'
!          call goVarValue( mdata%query, ';', 'platforms', '=', key, status )
!          IF_ERROR_RETURN(status=1)
!        ! split into platforms:
!        call goSplitString( trim(key), mdata%modis_nplatform, mdata%modis_platforms, status, sep=',' )
!        IF_NOTOK_RETURN(status=1)
!
!        ! space seperated list with product names (should be one per platform):
!        key = 'no-products'
!          call goVarValue( mdata%query, ';', 'products', '=', key, status )
!          IF_ERROR_RETURN(status=1)
!        ! split into platforms:
!        call goSplitString( trim(key), n, mdata%modis_products, status, sep=',' )
!        IF_NOTOK_RETURN(status=1)
!        ! check ..
!        if ( n /= mdata%modis_nplatform ) then
!          write (gol,'("numper of products defined in query is ",i4," while number of platforms is ",i4)') n, mdata%modis_nplatform; call goErr
!          TRACEBACK; status=1; return
!        end if
!    
!        ! list files; loop over platforms:
!        do ipl = 1, mdata%modis_nplatform
!          ! full path:
!          write (mdata%modis_list_file(ipl),'(a,"/",a,"/",a,".list")') &
!                   trim(mdata%modis_path), trim(mdata%modis_platforms(ipl)), trim(mdata%modis_products(ipl))
!          ! check ...
!          inquire( file=trim(mdata%modis_list_file(ipl)), exist=exist )
!          if ( .not. exist ) then
!            write (gol,'("list file not found: ",a)') trim(mdata%modis_list_file(ipl)); call goErr
!            write (gol,'("this should be a text file with on every line the name of a data file")'); call goErr
!            TRACEBACK; status=1; return
!          end if
!          ! free file unit:
!          call goGetFU( mdata%modis_list_fu(ipl), status )
!          IF_NOTOK_RETURN(status=1)
!          ! open :
!          open( mdata%modis_list_fu(ipl), file=trim(mdata%modis_list_file(ipl)), &
!                  status='old', form='formatted', iostat=status )
!          if (status/=0) then
!            write (gol,'("opening modis list file:")'); call goErr
!            write (gol,'("  platform    : ",a," (",i2,")")') trim(mdata%modis_platforms(ipl)), ipl; call goErr
!            write (gol,'("  list file   : ",a)') trim(mdata%modis_list_file(ipl)); call goErr
!            TRACEBACK; status=1; return
!          end if
!          ! set flag:
!          mdata%modis_list_opened(ipl) = .true.
!          ! nothing read yet:
!          mdata%modis_list_iline(ipl) = 0
!          mdata%modis_list_t(ipl) = AnyDate()
!        end do
!        
!      ! ~ seviri:
!      case ( 'seviri' )
!      
!        ! path to archive:
!        mdata%seviri_path = 'no-path'
!          call goVarValue( mdata%query, ';', 'path', '=', mdata%seviri_path, status )
!          IF_ERROR_RETURN(status=1)
!
!        ! prefix for filenames:
!        mdata%seviri_prefix = 'no-prefix'
!          call goVarValue( mdata%query, ';', 'prefix', '=', mdata%seviri_prefix, status )
!          IF_ERROR_RETURN(status=1)

      ! ~ unknown ...
      case default

        write (gol,'("unsupported instrument : ",a)') trim(mdata%instrument); call goErr
        TRACEBACK; status=1; return
    
    end select
    
    ! maximum number of pixels:
    mdata%pix_max = 1000
      call goVarValue( mdata%query, ';', 'maxpix', '=', mdata%pix_max, status )
      IF_ERROR_RETURN(status=1)

    ! time resolution key:
    call ReadRc( rcF, trim(rckey)//'.tres', key, status )
    IF_NOTOK_RETURN(status=1)
    ! extract values:
    mdata%tres_period = 'day'
      call goVarValue( key, ';', 'period', '=', mdata%tres_period, status )
      IF_ERROR_RETURN(status=1)
    !mdata%tres_type = 'aver'
    !  call goVarValue( key, ';', 'type', '=', mdata%tres_type, status )
    !  IF_ERROR_RETURN(status=1)
    !mdata%tres_unit = 'hour'
    !  call goVarValue( key, ';', 'unit', '=', mdata%tres_unit, status )
    !  IF_ERROR_RETURN(status=1)
    !mdata%tres_step = 1.0
    !  call goVarValue( key, ';', 'step', '=', mdata%tres_step, status )
    !  IF_ERROR_RETURN(status=1)
    !mdata%tres_hist = .false.
    !  call goVarValue( key, ';', 'hist', '=', mdata%tres_hist, status )
    !  IF_ERROR_RETURN(status=1)
    !! set time resoultion start and increment:
    !select case ( trim(mdata%tres_unit) )
    !  case ( 'hour' )
    !    ! maximum number of steps in a period:
    !    select case ( trim(mdata%tres_period) )
    !      case ( 'day' ) ; maxstep = 24
    !      case default
    !        write (gol,'("unsupported period to determine nstep : ",a)') trim(mdata%tres_period); call goErr
    !        TRACEBACK; status=1; return
    !    end select
    !    ! start hour:
    !    mdata%tres_start = t0%hour
    !      call goVarValue( key, ';', 'start', '=', mdata%tres_start, status )
    !      IF_ERROR_RETURN(status=1)
    !    ! fill start time:
    !    mdata%tres_t0 = NewDate( year=t0%year, month=t0%month, day=t0%day, hour=nint(mdata%tres_start) )
    !    ! fill time step:
    !    if ( modulo(real(maxstep),mdata%tres_step) /= 0.0 ) then
    !      write (gol,'("tres hourly only support integer number of steps within ",i4," hour:")') mdata%nstep; call goErr
    !      write (gol,'("  tres (hours)    : ",f8.2)') mdata%tres_step; call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !    mdata%tres_dt = IncrDate(hour=1) * mdata%tres_step
    !    ! number of steps within a period:
    !    mdata%nstep = nint(real(maxstep)/mdata%tres_step)
    !    ! store history ?
    !    if ( mdata%tres_hist ) then
    !      ! store history for whole file:
    !      mdata%nhist = mdata%nstep
    !    else
    !      ! no history, only one value:
    !      mdata%nhist = 1
    !    end if
    !  !case ( 'day'  )
    !  !  mdata%tres_t0 = NewData( year=t0%year, month=t0%month, day=mdata%tres_start )
    !  !  mdata%tres_dt = IncrDate( day=1)*mdata%tres_step
    !  case default
    !    write (gol,'("unsupported tres unit : ",a)') trim(mdata%tres_unit)
    !    TRACEBACK; status=1; return
    !end select
    
    !! setup sample averaging:
    !call ISTA_Init( mdata%ista, mdata%tres_t0, mdata%tres_dt, status )
    !IF_NOTOK_RETURN(status=1)
    !! init flags:
    !mdata%finished_aver = .false.
    
    ! init time range:
    mdata%setup_tr(1) = t0
    mdata%setup_tr(2) = t0
    ! init period time range:
    mdata%period_tr(1) = AnyDate()
    mdata%period_tr(2) = AnyDate()
    
    ! profile ? otherwise surface only:
    call ReadRc( rcF, trim(rckey)//'.profile', mdata%profile, status )
    IF_NOTOK_RETURN(status=1)

    ! store pixels ?
    call ReadRc( rcF, trim(rckey)//'.pixels', mdata%store_pix, status )
    IF_NOTOK_RETURN(status=1)
    ! store grid cell averages ?
    call ReadRc( rcF, trim(rckey)//'.grided', mdata%store_grid, status )
    IF_NOTOK_RETURN(status=1)

    ! read space seperated list with variables:
    call ReadRc( rcF, trim(rckey)//'.var', key, status )
    IF_NOTOK_RETURN(status=1)
    ! store list:
    call MAORI_VarList_Init( mdata%varlist, key, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read space seperated list with input variables:
    call ReadRc( rcF, trim(rckey)//'.obs.var', key, status )
    IF_NOTOK_RETURN(status=1)
    ! store list:
    call MAORI_VarList_Init( mdata%obs_varlist, key, status )
    IF_NOTOK_RETURN(status=1)
    ! store number:
    mdata%obs_nvar = mdata%obs_varlist%n
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! single level for the moment ...
      mdata%obs_nlev = 1
      ! assumed error standard deviation:
      call ReadRc( rcF, trim(rckey)//'.obs.sigma', key, status )
      IF_NOTOK_RETURN(status=1)
      ! extract values:
      mdata%sigma_offset = 0.0
        call goVarValue( key, ';', 'offset', '=', mdata%sigma_offset, status )
        IF_ERROR_RETURN(status=1)
      mdata%sigma_scale = 0.0
        call goVarValue( key, ';', 'scale', '=', mdata%sigma_scale, status )
        IF_ERROR_RETURN(status=1)
    end if

    ! analyse in assimilation ?
    call ReadRc( rcF, trim(rckey)//'.assim.analyse', mdata%assim_analyse, status )
    IF_NOTOK_RETURN(status=1)
    ! assimilation settings:
    if ( mdata%assim_analyse ) then
      ! check ...
      if ( mdata%obs_nvar < 1 ) then
        write (gol,'("could not analyse `",a,"` if no observations are available ...")') trim(name); call goErr
        TRACEBACK; status=1; return
      end if
      ! spatial correlation scale:
      call ReadRc( rcF, trim(rckey)//'.assim.rho', mdata%assim_rho, status )
      IF_NOTOK_RETURN(status=1)
      ! screening factor:
      call ReadRc( rcF, trim(rckey)//'.assim.alfa', mdata%assim_alfa, status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      mdata%assim_rho  = 0.0   ! no correlation
      mdata%assim_alfa = 999.9   ! no screening ...
    end if
    
    ! no parameters set yet:
    mdata%param_set = .false.
    
    ! not started yet:
    mdata%started = .false.

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Init


  ! ***


  subroutine MAORI_Data_Satellite_Start( mdata, status )

    !use GO, only : AnyDate
    use GO, only : goMatchValue, goVarValue
    use Grid, only : Init
    !use MAORI_LocList, only : MAORI_LocList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)  ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Start'
    
    ! --- local ----------------------------------
    
    integer                         ::  ivar
    real                            ::  dlon, dlat
    !character(len=MAORI_LEN_NAME)   ::  obs_type

    ! --- begin -----------------------------------
    
    ! all parameters set ?
    if ( .not. all(mdata%param_set) ) then
      write (gol,'("not all parameters set: ")'); call goErr
      write (gol,'("  ",1000l2)') mdata%param_set
      TRACEBACK; status=1; return
    end if

    ! grid resolution:
    dlon = ( mdata%bound_east  - mdata%bound_west  )/mdata%grid_nlon
    dlat = ( mdata%bound_north - mdata%bound_south )/mdata%grid_nlat
    ! setup grid description:
    call Init( mdata%grid_lli, mdata%bound_west +dlon/2.0, dlon, mdata%grid_nlon, &
                               mdata%bound_south+dlat/2.0, dlat, mdata%grid_nlat, status )
    IF_NOTOK_RETURN(status=1)

    ! storage for pixels:
    if ( mdata%store_pix ) then
      !mdata%pix_max = mdata%pix_max  ! already set ...
      mdata%pix_n = 0  ! empty
      ! storage:
      allocate( mdata%pix_date   (6,mdata%pix_max) )
      allocate( mdata%pix_time     (mdata%pix_max) )
      allocate( mdata%pix_longitude(mdata%pix_max) )
      allocate( mdata%pix_latitude (mdata%pix_max) )
      allocate( mdata%pix_proba    (mdata%pix_max) )
      allocate( mdata%pix_slots    (mdata%pix_max) )
      ! mapping ..
      if ( mdata%store_grid ) then
        allocate( mdata%pix_grid_i(mdata%pix_max) )
        allocate( mdata%pix_grid_j(mdata%pix_max) )
      end if
    else
      ! empty ...
      mdata%pix_n = 0
    end if
    
    ! storage for grid cell averages:
    if ( mdata%store_grid ) then
      ! dimensions:
      mdata%grid_nx = mdata%grid_nlon
      mdata%grid_ny = mdata%grid_nlat
      ! storage:
      allocate( mdata%grid_npix(mdata%grid_nx,mdata%grid_ny) )
      allocate( mdata%grid_w24 (mdata%grid_nx,mdata%grid_ny,0:24) )
    else
      ! dimensions:
      mdata%grid_nx = 0
      mdata%grid_ny = 0
    end if

    ! extract dimensions:
    mdata%nvar = mdata%varlist%n
    
    ! setup storage for number of output levels:
    allocate( mdata%nlev(mdata%nvar) )
    ! set maximum and actual number of output levels:
    do ivar = 1, mdata%nvar
      ! profile output or surface only ?
      if ( mdata%profile ) then
        ! all layers plus 1 extra for half levels
        mdata%maxlev = mdata%nlayer+1
        ! half level field ?
        if ( mdata%varlist%var(ivar)%halflevel ) then
          mdata%nlev(ivar) = mdata%nlayer + 1
        else
          mdata%nlev(ivar) = mdata%nlayer
        end if
      else
        ! single level output only:
        mdata%maxlev     = 1
        mdata%nlev(ivar) = 1
      end if
    end do
    
    ! number of state values on a single location:
    mdata%nlocvalue = 0
    do ivar = 1, mdata%nvar
      ! add contribution of this variable (times 3 for 'curr', 'prev', and 'aver'):
      mdata%nlocvalue = mdata%nlocvalue + mdata%nlev(ivar) * 3
    end do
    
    ! current total number of state values;
    ! will be computed online in MAORI_Data_Satellite_Setup:
    mdata%nvalue = 0
    
    ! info ...
    write (gol,'("  MAORI model variables:")'); call goPr
    write (gol,'("    -------------------------")'); call goPr
    write (gol,'("      nr name                ")'); call goPr
    write (gol,'("    ---- --------------------")'); call goPr
    do ivar = 1, mdata%nvar
      write (gol,'("    ",i4," ",a20)') ivar, trim(mdata%varlist%var(ivar)%name); call goPr
    end do
    write (gol,'("    -------------------------")'); call goPr
    
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then

      ! storage:
      if ( mdata%store_pix ) then
        allocate( mdata%obs_pix_value (mdata%obs_nlev,mdata%pix_max,mdata%obs_nvar) )
        allocate( mdata%obs_pix_errstd(mdata%obs_nlev,mdata%pix_max,mdata%obs_nvar) )
        allocate( mdata%obs_pix_astat (mdata%obs_nlev,mdata%pix_max,mdata%obs_nvar) )
      end if
      if ( mdata%store_grid ) then
        allocate( mdata%obs_grid_value (mdata%grid_nx,mdata%grid_ny,mdata%obs_nlev,mdata%obs_nvar) )
        allocate( mdata%obs_grid_errstd(mdata%grid_nx,mdata%grid_ny,mdata%obs_nlev,mdata%obs_nvar) )
        allocate( mdata%obs_grid_astat (mdata%grid_nx,mdata%grid_ny,mdata%obs_nlev,mdata%obs_nvar) )
      end if

      ! output variable corresponding to observed variable:    
      allocate( mdata%obs_ivar(mdata%obs_nvar) )
      ! loop over observed variables:
      do ivar = 1, mdata%obs_nvar
        ! search in list of input variables:
        call goMatchValue( trim(mdata%obs_varlist%var(ivar)%name), &
                mdata%varlist%var(1:mdata%nvar)%name, mdata%obs_ivar(ivar), status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! info ...
      write (gol,'("  MAORI observed variables:")'); call goPr
      write (gol,'("    ----------------------------------------------")'); call goPr
      write (gol,'("      nr name                 model variable      ")'); call goPr
      write (gol,'("    ---- -------------------- --------------------")'); call goPr
      do ivar = 1, mdata%obs_nvar
        write (gol,'("    ",i4," ",a20," ",a20)') ivar, trim(mdata%obs_varlist%var(ivar)%name), trim(mdata%varlist%var(mdata%obs_ivar(ivar))%name); call goPr
      end do
      write (gol,'("    ----------------------------------------------")'); call goPr

    end if  ! obs_nvar > 0

    ! reset flag:
    mdata%started = .true.

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Start


  ! ***


  subroutine MAORI_Data_Satellite_Setup( mdata, t1, t2, the_end, status )

    use GO            , only : NewDate, IncrDate, AnyDate, IsAnyDate, Get
    use GO            , only : operator(>), operator(<=), operator(+), operator(-), rTotal
    use GO            , only : operator(==), operator(/=), operator(<=), operator(>=)
    use GO            , only : wrtgol
    use GO            , only : goReadFromLine
    !use GO           , only : goReplace, goVarValue
    !use GO           , only : ISTA_Add, ISTA_Next, ISTA_Shift
    use Grid   , only : GetLocation

!    use MODIS  , only : T_MODIS, Init, Done
!    use SEVIRI , only : T_SEVIRI, Init, Done

    !use EMEP_DataFile, only : Init, Done, FindRecord, Get, GetObservation

    use MAORI_Common , only : com_t0
    use MAORI_Param  , only : MAORI_LEN_FILE, MAORI_NAN
    use MAORI_Param  , only : MAORI_ASTAT_VALIDATION
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)  ::  mdata
    type(TDate), intent(in)                   ::  t1, t2
    logical, intent(in)                       ::  the_end
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Setup'

    ! --- local -----------------------------------

    ! time settings stuff:
    integer                         ::  date1(6)
    !real                            ::  w_aver, w_prev, w_curr
    !real                            ::  frac

    ! ~ common:
    integer                             ::  year, month, day, hour, minu
    real                                ::  value, errstd
    
!    ! ~ modis specials
!    character(len=32)                   ::  prompt
!    integer                             ::  ipl
!    integer                             ::  i, j
!    character(len=MAORI_LEN_FILE)       ::  line
!    character(len=32)                   ::  key
!    integer                             ::  jday
!    type(T_MODIS)                       ::  tmodis
!    character(len=MAORI_LEN_FILE)       ::  fname
!    logical                             ::  indomain
!    logical                             ::  isfirst
!    real(8)                             ::  sst_first
!    type(TDate)                         ::  t_first
!    type(TDate)                         ::  pix_t0, pix_t
!    integer                             ::  time6(6)
!    integer                             ::  ipix
!    
!    ! ~ severi specials
!    type(T_SEVIRI)                  ::  seviri
!    integer                         ::  np, ip
!    real                            ::  plon, plat
!    integer                         ::  ilon, ilat
!    integer                         ::  islot
!    logical                         ::  slots(32)
!    integer                         ::  ih
!    real                            ::  w
!    integer                         ::  used_ncell, used_nslot
    
    !character(len=MAORI_LEN_NAME)   ::  obs_type
    !character(len=MAORI_LEN_FILE)   ::  fname
    !integer                         ::  nobs, iobs
    !character(len=MAORI_LEN_NAME)   ::  vname, aname
    !character(len=MAORI_LEN_NAME)   ::  station_code
    !integer                         ::  iloc
    integer                         ::  ivar
    !integer                         ::  k
    !type(TDate)                     ::  tcur
    !character(len=MAORI_LEN_NAME)   ::  comp, unit
    !real                            ::  val
    !real                            ::  unitfac
    
    ! --- begin -----------------------------------
    
    !
    ! time step settings
    !
    
    ! setup at begin of time step over [t1,t2] :
    !  o a sample will be received at t2
    
    ! info ...
    write (gol,'("MAORI/SAT Setup:  setup for ",a)') trim(mdata%name); call goPr
    !call wrtgol( 'MAORI/SAT Setup:    time step       : ', t1, ' - ', t2 ); call goPr
    
    ! check ...
    if ( t1 /= mdata%setup_tr(2) ) then
      write (gol,'("setup interval not connected with previous interval:")'); call goErr
      call wrtgol( '  previous setup interval : ', mdata%setup_tr(1), ' - ', mdata%setup_tr(2) ); call goErr
      call wrtgol( '  t1, t2                  : ', t1, ' - ', t2 ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract date arrays:
    call Get( t1, time6=date1 )
    
    ! new period ?
    if ( IsAnyDate(mdata%period_tr(1)) ) then
      mdata%new_period = .true.
    else
      mdata%new_period =  t2 > mdata%period_tr(2)
    end if
    ! info ...
    write (gol,'("MAORI Setup:    new period      : ",l2)') mdata%new_period; call goPr
    
    ! set time period including t1:
    select case ( trim(mdata%tres_period) )
      case ( 'day' )
        mdata%period_tr(1) = NewDate( year=date1(1), month=date1(2), day=date1(3), hour=00 )
        mdata%period_tr(2) = mdata%period_tr(1) + IncrDate(hour=24)
      case default
        write (gol,'("could not determine time range period for : ",a)') trim(mdata%tres_period); call goErr
        TRACEBACK; status=1; return
    end select
    ! check ...
    if ( t2 > mdata%period_tr(2) ) then
      write (gol,'("setup period partly outside periode defined by t1:")'); call goErr
      call wrtgol( '  t1, t2     : ', t1, ' - ', t2 ); call goErr
      call wrtgol( '  period     : ', mdata%period_tr(1), ' - ', mdata%period_tr(2) ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! end of period ?
    mdata%end_of_period = t2 == mdata%period_tr(2)
    
    !! info ...
    !call wrtgol( 'MAORI Setup:    current period  : ', mdata%period_tr(1), ' - ', mdata%period_tr(2) ); call goPr
    !write (gol,'("MAORI Setup:    end of period   : ",l2)') mdata%end_of_period; call goPr

    !! previously finished an averaging interval ? then step to new inteval:
    !if ( mdata%finished_aver ) then
    !  !! info ...
    !  !write (gol,'("MAORI Setup:    next aver interval ...")'); call goPr
    !  ! init new interval; return weights to build new average:
    !  !    aver :=   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
    !  call ISTA_Next( mdata%ista, w_prev, w_curr, status )
    !  IF_NOTOK_RETURN(status=1)
    !  ! store:
    !  mdata%w_aver = 0.0
    !  mdata%w_prev = w_prev
    !  mdata%w_curr = w_curr
    !  ! reset flag:
    !  mdata%finished_aver = .false.
    !  ! increase step counter:
    !  mdata%istep = mdata%istep + 1
    !  if ( mdata%istep > mdata%nstep ) mdata%istep = 1
    !  !! info ...
    !  !write (gol,'("MAORI Setup:      weights       : ",3f6.2)') 0.0, w_prev, w_curr; call goPr
    !  !write (gol,'("MAORI Setup:      istep         : ",i6)') mdata%istep; call goPr
    !else
    !  ! reset weights:
    !  !    aver :=   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
    !  mdata%w_aver = 1.0
    !  mdata%w_prev = 0.0
    !  mdata%w_curr = 0.0
    !end if
    !
    !! obtain weights and flags for samples added at t2:
    !call ISTA_Add( mdata%ista, t2, w_aver, w_prev, w_curr, mdata%finished_aver, status )
    !IF_NOTOK_RETURN(status=1)
    !! add weights:
    !!    aver :=   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
    !!    aver :=   aver * w_aver2 + prev * w_prev2 + curr * w_curr2
    !!          =   aver * (w_aver1 * w_aver2          ) + 
    !!              prev * (w_prev1 * w_aver2 + w_prev2) + 
    !!              curr * (w_curr1 * w_aver2 + w_curr2)
    !mdata%w_aver = mdata%w_aver * w_aver
    !mdata%w_prev = mdata%w_prev * w_aver + w_prev
    !mdata%w_curr = mdata%w_curr * w_aver + w_curr
    !!! info ...
    !!write (gol,'("MAORI Setup:    add weights     : ",3f6.2)') w_aver, w_prev, w_curr; call goPr
    !!write (gol,'("MAORI Setup:    tot weights     : ",3f6.2)') mdata%w_aver, mdata%w_prev, mdata%w_curr; call goPr
    !! check ...
    !if ( abs( mdata%w_aver + mdata%w_prev + mdata%w_curr - 1.0 ) > 0.01 ) then
    !  write (gol,'("weights do not sum to 1 : ")'); call goErr
    !  write (gol,'("  w_aver, w_prev, w_curr  : ",3f8.4)') mdata%w_aver, mdata%w_prev, mdata%w_curr; call goErr
    !  write (gol,'("  sum                     : ",3f8.4)') mdata%w_aver+mdata%w_prev+mdata%w_curr; call goErr
    !  TRACEBACK; status=1; return
    !end if
    !
    !! determine step within period; trap initial setup for instant time:
    !frac = rTotal( t2 - mdata%period_tr(1), mdata%tres_unit ) / mdata%tres_step
    !mdata%istep = ceiling(frac)
    !if ( (t1 == t2) .and. (frac == ceiling(frac)) ) mdata%istep = mdata%istep + 1
    !! check ...
    !if ( (mdata%istep < 1) .or. (mdata%istep > mdata%nstep) ) then
    !  write (gol,'("step index ",i6," out of range 1 .. ",i6)') mdata%istep, mdata%nstep; call goErr
    !  TRACEBACK; status=1; return
    !end if
    !
    !! history:
    !if ( mdata%tres_hist ) then
    !  mdata%ihist = mdata%istep
    !else
    !  mdata%ihist = 1
    !end if
    !
    !!! info ...
    !!write (gol,'("MAORI Setup:    current step    : ",i6," / ",i6)') mdata%istep, mdata%nstep; call goPr
    !!write (gol,'("MAORI Setup:    current hist    : ",i6," / ",i6)') mdata%ihist, mdata%nhist; call goPr
    !
    !! shift current into previous:
    !call ISTA_Shift( mdata%ista, status )
    !IF_NOTOK_RETURN(status=1)
    
    ! store:
    mdata%setup_tr(1) = t1
    mdata%setup_tr(2) = t2
    mdata%the_end     = the_end
    
    !
    ! observed data
    !
        
    !! clear buffer:
    !mdata%nloc = 0
    
    ! instrument dependend stuff:
    select case ( trim(mdata%instrument) )
    
!      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      case ( 'modis' )
!      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        
!        ! clear buffer:
!        mdata%pix_n = 0
!    
!        ! loop over platforms:
!        do ipl = 1, mdata%modis_nplatform
!
!          ! set prompt:
!          prompt = 'MAORI: MODIS/'//trim(mdata%modis_platforms(ipl))//': '
!          ! info:
!          call wrtgol( trim(prompt)//'setup for ', t1, ' - ', t2 ); call goPr
!
!          ! if list file is closed, no new data is available anymore ...
!          if ( .not. mdata%modis_list_opened(ipl) ) cycle
!
!          ! loop until time fits:
!          do
!
!            ! item available ?
!            if ( .not. IsAnyDate(mdata%modis_list_t(ipl)) ) then
!
!              ! too new ? then leave:
!              if ( mdata%modis_list_t(ipl) > t2 ) then
!                call wrtgol( trim(prompt)//'  keep orbit in buffer; valid for : ', mdata%modis_list_t(ipl) ); call goPr
!                exit
!              end if
!
!              ! in interval (t1,t2] ?
!              if ( (mdata%modis_list_t(ipl) > t1) .and. (mdata%modis_list_t(ipl) <= t2) ) then
!                ! info ...
!                write (gol,'(a,"  read orbit ...")') trim(prompt); call goPr
!                ! create full file name: <archive>/<platform>/<orbitfile>
!                write (fname,'(a,"/",a,"/",a)')  trim(mdata%modis_path), trim(mdata%modis_platforms(ipl)), trim(mdata%modis_list_dfile(ipl))
!                !! info ...
!                !write (gol,'(a,"    filename     : ",a)') trim(prompt), trim(fname); call goPr
!                ! read input file:
!                call Init( tmodis, trim(fname), status )
!                IF_NOTOK_RETURN(status=1)
!                ! handle time:
!                if ( trim(tmodis%Scan_Start_Time__units) == 'Seconds since 1993-1-1 00:00:00.0 0' ) then
!                  pix_t0 = NewDate( 1993, 01, 01, 00, 00, 00 )
!                else
!                  write (gol,'("unsupported MODIS time unit : ",a)') trim(tmodis%Scan_Start_Time__units); call goErr
!                  TRACEBACK; status=1; return
!                end if
!                !! copy units:
!                !mdata%pix_Cloud_Fraction_Land__units              = modis%Cloud_Fraction_Land__units
!                !mdata%pix_Cloud_Fraction_Ocean__units             = modis%Cloud_Fraction_Ocean__units
!                !mdata%pix_Optical_Depth_Land_And_Ocean__units     = modis%Optical_Depth_Land_And_Ocean__units
!                !! copy long_name:
!                !mdata%pix_Cloud_Fraction_Land__long_name          = modis%Cloud_Fraction_Land__long_name
!                !mdata%pix_Cloud_Fraction_Ocean__long_name         = modis%Cloud_Fraction_Ocean__long_name
!                !mdata%pix_Optical_Depth_Land_And_Ocean__long_name = modis%Optical_Depth_Land_And_Ocean__long_name
!                !! info ...
!                !write (gol,'(a,"    checking new pixels : ",i6)') trim(prompt), modis%nx*modis%ny; call goPr
!                ! loop over pixels:
!                isfirst = .true.
!                do j = 1, tmodis%ny
!                  do i = 1, tmodis%nx
!                    ! no-data ? then skip:
!                    if ( tmodis%Optical_Depth_Land_And_Ocean(i,j) < 0.0 ) cycle
!                    ! check location:
!                    indomain = (tmodis%Longitude(i,j) >= mdata%bound_west ) .and. (tmodis%Longitude(i,j) <= mdata%bound_east ) .and. &
!                               (tmodis%Latitude (i,j) >= mdata%bound_south) .and. (tmodis%Latitude (i,j) <= mdata%bound_north)
!                    ! not in domain ? then skip:
!                    if ( .not. indomain ) cycle
!                    ! full ?
!                    if ( mdata%pix_n == mdata%pix_max ) then
!                      write (gol,'("reached maximum number of pixels : ",i6)') mdata%pix_max; call goErr
!                      write (gol,'("please increase value of `maxpix` in query ...")'); call goErr
!                      TRACEBACK; status=1; return
!                    end if
!                    ! pixel counter:
!                    ipix = mdata%pix_n + 1
!                    !! store platform index:
!                    !mdata%pix_ipl (ipix) = ipl
!                    !! store julian day:
!                    !mdata%pix_jday(ipix) = mdata%modis_list_jday(ipl)
!                    ! pre-compute expensive time stuf:
!                    if ( isfirst ) then
!                      ! store:
!                      sst_first = tmodis%Scan_Start_Time(i,j)
!                      ! full time:
!                      t_first = pix_t0 + IncrDate( sec=nint(sst_first) )
!                      ! reset flag:
!                      isfirst = .false.
!                    end if
!                    ! store ?
!                    if ( mdata%store_pix ) then
!                      ! convert from 'seconds since ...' to year,month,etc
!                      pix_t = t_first + IncrDate( sec=nint(tmodis%Scan_Start_Time(i,j)-sst_first) )
!                      call Get( pix_t, time6=time6 )
!                      mdata%pix_date(:,ipix) = time6
!                      ! convert to 'days since t0' :
!                      mdata%pix_time(ipix) = rTotal( pix_t - com_t0, 'day' )
!                      ! copy data:
!                      mdata%pix_longitude(ipix)                    = tmodis%Longitude(i,j)
!                      mdata%pix_latitude (ipix)                    = tmodis%Latitude (i,j)
!                      !mdata%pix_Cloud_Fraction_Land (ipix)         = modis%Cloud_Fraction_Land (i,j)                    
!                      !mdata%pix_Cloud_Fraction_Ocean(ipix)         = modis%Cloud_Fraction_Ocean(i,j)                    
!                      !mdata%pix_Optical_Depth_Land_And_Ocean(ipix) = modis%Optical_Depth_Land_And_Ocean(i,j)     
!                      ! no slots ...
!                      mdata%pix_slots(ipix) = 0
!                    end if
!
!                    ! observations ?
!                    if ( mdata%obs_nvar > 0 ) then
!                      ! loop over observed variables:
!                      do ivar = 1, mdata%obs_nvar
!                        ! copy:
!                        select case ( trim(mdata%obs_varlist%var(ivar)%name) )
!                          case ( 'aod' )
!                            ! copy value:
!                            mdata%obs_pix_value(1,ipix,ivar) = tmodis%Optical_Depth_Land_And_Ocean(i,j)
!                            ! assumed error std.dev. is a base value + fraction of the observed value:
!                            mdata%obs_pix_errstd(1,ipix,ivar) = &
!                                    mdata%sigma_offset + mdata%sigma_scale * mdata%obs_pix_value(1,ipix,ivar)
!                          case default
!                            write (gol,'("unsupported obs var : ",a)') trim(mdata%obs_varlist%var(ivar)%name); call goErr
!                            TRACEBACK; status=1; return
!                        end select
!                      end do
!                    end if
!                    !mdata%pix_obs(ipix)                         = modis%Optical_Depth_Land_And_Ocean(i,j)
!                    !! assumed error std.dev. is a base value + fraction of the observed value:
!                    !mdata%pix_obs_sigma(ipix) = &
!                    !      mdata%sigma_offset + mdata%sigma_scale * mdata%pix_obs(ipix) ! m
!                    !! grid cell with pixel center:
!                    !call GetLocation( lli, mdata%pix_Longitude(ipix), mdata%pix_Latitude(ipix), &
!                    !                       mdata%pix_ix(ipix), mdata%pix_iy(ipix), status )
!                    !! set flags:
!                    !!mdata%pix_assim_status(ipix) = status_default
!                    !mdata%pix_assim_status(ipix) = 0
!
!                    ! store counter:
!                    mdata%pix_n = ipix
!                  end do
!                end do
!                ! info ...
!                write (gol,'(a,"    last pixel   : ",i6)') trim(prompt), mdata%pix_n; call goPr
!                ! close:
!                call Done( tmodis, status )
!                IF_NOTOK_RETURN(status=1)
!              else
!                ! file too old ...
!                call wrtgol( trim(prompt)//'  orbit in buffer too old : ', mdata%modis_list_t(ipl) ); call goPr
!              end if
!
!            end if
!
!            ! increase counter:
!            mdata%modis_list_iline(ipl) = mdata%modis_list_iline(ipl) + 1
!
!            ! read next line:
!            read (mdata%modis_list_fu(ipl),'(a)',iostat=status) line
!            ! eof ?
!            if ( status < 0 ) then
!              ! close:
!              close( unit=mdata%modis_list_fu(ipl), iostat=status )
!              IF_NOTOK_RETURN(status=1)
!              ! reset flag:
!              mdata%modis_list_opened(ipl) = .false.
!              ! dummy time ..
!              mdata%modis_list_t(ipl) = AnyDate()
!              ! leave list loop:
!              exit
!            end if
!            ! error ...
!            if (status/=0) then
!              write (gol,'("reading line from list file : ")'); call goErr
!              write (gol,'("  file   : ",a)') trim(mdata%modis_list_file(ipl)); call goErr
!              write (gol,'("  line   : ",i6)') trim(line); call goErr
!              TRACEBACK; status=1; return
!            end if
!
!            ! store for later usage:
!            mdata%modis_list_dfile(ipl) = trim(line)
!
!            ! example file name:
!            !   MYD04_L2.A2006001.0850.005.2007128013650.hdf
!
!            ! extract fields:
!            ! MYD04_L2
!            call goReadFromLine( line, key, status, sep='.' )
!            IF_NOTOK_RETURN(status=1)
!            ! A2006001
!            call goReadFromLine( line, key, status, sep='.' )
!            IF_NOTOK_RETURN(status=1)
!            read (key(2:8),'(i4,i3)',iostat=status) year, jday
!            if (status/=0 ) then
!              write (gol,'("reading year and jday from line : ",a)') trim(line); call goErr
!              TRACEBACK; status=1; return
!            end if
!            ! 0850
!            call goReadFromLine( line, key, status, sep='.' )
!            IF_NOTOK_RETURN(status=1)
!            read (key,'(2i2)',iostat=status) hour, minu
!            if (status/=0 ) then
!              write (gol,'("reading hour and minu from line : ",a)') trim(line); call goErr
!              TRACEBACK; status=1; return
!            end if
!
!            ! set time; jday 1..365 should be converted to correct month/day :
!            mdata%modis_list_jday(ipl) = jday
!            mdata%modis_list_t   (ipl) = NewDate( year=year, month=1, day=jday, hour=hour, min=minu )
!
!          end do  ! lines in file list
!
!        end do  ! platforms
!        
!      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      case ( 'seviri' )
!      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!        ! new period started ? then read:
!        if ( mdata%new_period ) then
!        
!          ! day files; extract day from t1
!          call Get( t1, year=year, month=month, day=day )
!          ! store current input time range:
!          mdata%period_tr(1) = NewDate( year=year, month=month, day=day, hour=00, min=00 )
!          mdata%period_tr(2) = mdata%period_tr(1) + IncrDate(day=1)
!
!          ! filename:
!          write (mdata%seviri_basename,'(a,"_",i4.4,2i2.2)') &
!                              trim(mdata%seviri_prefix), year, month, day
!
!          ! read all data for this day:
!          call Init( seviri, trim(mdata%seviri_path), trim(mdata%seviri_basename), status )
!          IF_NOTOK_RETURN(status=1)
!
!          ! number of pixels:
!          np = seviri%np
!
!          ! info ...
!          write (gol,'("MAORI/SAT:   total number of pixels   : ",i8)') np; call goPr
!
!          ! no pixels yet:
!          if ( mdata%store_pix ) then
!            mdata%pix_n = 0
!          end if
!          ! no pixels averaged yet:
!          if ( mdata%store_grid ) then
!            ! no average:
!            mdata%grid_npix = 0
!            mdata%grid_w24 = 0.0
!            ! for info ...
!            used_ncell = 0
!            ipix  = 0
!            used_nslot = 0
!          end if
!
!          ! read observations ?
!          if ( mdata%obs_nvar > 0 ) then
!            ! no observations averaged yet ...
!            if ( mdata%store_grid ) then
!              mdata%obs_grid_value  = 0.0
!              mdata%obs_grid_errstd = 0.0
!            end if
!          end if
!          
!          ! current accepted pixel:
!          ipix = 0
!
!          ! loop over pixels:
!          do ip = 1, np
!
!            ! get location:
!            plon = seviri%Longitude(ip)
!            plat = seviri%Latitude (ip)
!            
!            ! get corrresponding grid cell:
!            call GetLocation( mdata%grid_lli, plon, plat, ilon, ilat, status )
!            if ( status < 0 ) cycle   ! not in domain, skip ...
!            IF_NOTOK_RETURN(status=1)
!
!            ! full ?
!            if ( mdata%pix_n == mdata%pix_max ) then
!              write (gol,'("reached maximum number of pixels : ",i6)') mdata%pix_max; call goErr
!              write (gol,'("please value of `maxpix` in query ...")'); call goErr
!              TRACEBACK; status=1; return
!            end if
!            ! increase pixel counter:
!            ipix = ipix + 1
!
!            ! info ...
!            !if ( ipix < 100 ) then; write (gol,'("MAORI/SAT:   add pixel    : ",i6)') ip; call goPr; end if
!
!            ! extract slots:
!            islot = seviri%UISlotSeries(ip)
!            ! convert to logicals:
!            do ih = 1, 32
!              slots(ih) = BTest( islot, ih-1 )
!            end do
!!! testing: sample last hour only:
!!slots = .false.
!!slots(24) = .true.
!            ! check ...
!            if ( (.not. any(slots(1:24))) .or. any(slots(25:32)) ) then
!              write (gol,'("no slots or after hour 24 in seviri pixel:")'); call goErr
!              write (gol,'("  file         : ",a)') trim(mdata%seviri_basename); call goErr
!              write (gol,'("  pixel        : ",i8)') ip; call goErr
!              write (gol,'("  slot number  : ",i8)') islot; call goErr
!              write (gol,'("  slots        : ",32l2)') slots; call goErr
!              TRACEBACK; status=1; return
!            end if
!
!            ! store pixels ?
!            if ( mdata%store_pix ) then
!              ! store:
!              mdata%pix_longitude(ipix) = plon
!              mdata%pix_latitude (ipix) = plat
!              mdata%pix_proba    (ipix) = seviri%Proba(ip)
!              mdata%pix_slots    (ipix) = islot
!              ! grid mapping ?
!              if ( mdata%store_grid ) then
!                mdata%pix_grid_i(ipix) = ilon
!                mdata%pix_grid_j(ipix) = ilat
!              end if
!            end if
!            
!            ! store grids ?
!            if ( mdata%store_grid ) then
!              ! assign pixel to grid cell:
!              mdata%grid_npix(ilon,ilat) = mdata%grid_npix(ilon,ilat) + 1
!              ! weight per slot:
!              w = 1.0 / real(count(slots))
!              ! loop over hours:
!              do ih = 1, 24
!                ! slot used ? than add weight:
!                if ( slots(ih) ) then
!                  ! contribution is average of instantaneous fields at begin and end of hour:
!                  mdata%grid_w24(ilon,ilat,ih-1) = mdata%grid_w24(ilon,ilat,ih-1) + 0.5 * w
!                  mdata%grid_w24(ilon,ilat,ih  ) = mdata%grid_w24(ilon,ilat,ih  ) + 0.5 * w
!                  ! info ...
!                  used_nslot = used_nslot + 1
!                  !if ( sum(mdata%grid_npix) < 100 ) then; write (gol,'("MAORI/SAT:     slot at    : ",i8)') ih; call goPr; end if
!                end if ! slot used
!              end do  ! hours
!            end if  ! store grid
!
!            ! observations ?
!            if ( mdata%obs_nvar > 0 ) then
!              ! loop over observed variables:
!              do ivar = 1, mdata%obs_nvar
!                ! copy:
!                select case ( trim(mdata%obs_varlist%var(ivar)%name) )
!                  case ( 'aod' )
!                    ! set values:
!                    value  = seviri%AOD(ip)
!                    errstd = seviri%AOD_Err(ip)
!                    ! store pixel:
!                    if ( mdata%store_pix ) then
!                      mdata%obs_pix_value (1,ipix,ivar) = value
!                      mdata%obs_pix_errstd(1,ipix,ivar) = errstd
!                      mdata%pix_n = ipix
!                    end if
!                    ! add values to sum:
!                    if ( mdata%store_grid ) then
!                      mdata%obs_grid_value (ilon,ilat,1,ivar) = mdata%obs_grid_value (ilon,ilat,1,ivar) + value
!                      !mdata%obs_grid_errstd(ilon,ilat,1,ivar) = mdata%obs_grid_errstd(ilon,ilat,1,ivar) + errstd**2  ! variance
!                    end if
!                  case default
!                    write (gol,'("unsupported obs var : ",a)') trim(mdata%obs_varlist%var(ivar)%name); call goErr
!                    TRACEBACK; status=1; return
!                end select
!              end do
!            end if  ! observerd values
!
!          end do  ! raw pixels
!          
!          ! post processing:
!          ! ~ weights
!          if ( mdata%store_grid ) then
!            ! loop over sampled hours:
!            do ih = 0, 24
!              ! compute grid cell averages:
!              where ( mdata%grid_npix > 0 )
!                mdata%grid_w24(:,:,ih) = mdata%grid_w24(:,:,ih) / real(mdata%grid_npix)
!              end where
!            end do  ! instant hours
!          end if  ! store grid
!          ! ~ observations ?
!          if ( mdata%obs_nvar > 0 ) then
!            ! loop over observed variables:
!            do ivar = 1, mdata%obs_nvar
!              !
!              if ( mdata%store_grid ) then
!                ! set to average if pixels are found, otherwise to nan:
!                where ( mdata%grid_npix > 0 )
!                  ! compute grid cell averages of observed values:
!                  mdata%obs_grid_value(:,:,1,ivar) = mdata%obs_grid_value(:,:,1,ivar) / real(mdata%grid_npix)
!                  ! compute average variance, convert to standard deviation
!                  !mdata%obs_grid_errstd(:,:,1,ivar) = sqrt( mdata%obs_grid_errstd(:,:,1,ivar) / real(mdata%grid_npix) )
!                  ! use offset and scaling:
!                  mdata%obs_grid_errstd(:,:,1,ivar) =  &
!                          mdata%sigma_offset + mdata%sigma_scale * mdata%obs_grid_value(:,:,1,ivar)
!                elsewhere
!                  ! no data ...
!                  mdata%obs_grid_value (:,:,1,ivar) = MAORI_NAN
!                  mdata%obs_grid_errstd(:,:,1,ivar) = MAORI_NAN
!                end where
!              end if
!            end do
!          end if  ! observerd values
!
!          ! info ...
!          used_ncell = count( mdata%grid_npix > 0 )
!          write (gol,'("MAORI/SAT:     selected slots         : ",i8)') used_nslot; call goPr
!          write (gol,'("MAORI/SAT:     selected pixels        : ",i8)') ipix ; call goPr
!          write (gol,'("MAORI/SAT:     selected cells         : ",i8)') used_ncell; call goPr
!
!          ! destroy:
!          call Done( seviri, status )
!          IF_NOTOK_RETURN(status=1)
!
!          !! set flag:
!          !mdata%inp_opened = .true.
!
!        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported instrument : ",a)') trim(mdata%instrument); call goErr
        TRACEBACK; status=1; return
    
    end select

    ! initialise assimilation status:
    if ( mdata%store_pix ) then
      ! clear all bits:
      mdata%obs_pix_astat = 0
      ! set validation bit if necessary:
      if ( .not. mdata%assim_analyse ) mdata%obs_pix_astat = IBSet( mdata%obs_pix_astat, MAORI_ASTAT_VALIDATION )
    end if
    if ( mdata%store_grid ) then
      ! clear all bits:
      mdata%obs_grid_astat = 0
      ! set validation bit if necessary:
      if ( .not. mdata%assim_analyse ) mdata%obs_grid_astat = IBSet( mdata%obs_grid_astat, MAORI_ASTAT_VALIDATION )
    end if

    ! set number of state values:
    mdata%nvalue = 0.0
    if ( mdata%store_pix  ) mdata%nvalue = mdata%nvalue + mdata%nlocvalue * mdata%pix_n
    if ( mdata%store_grid ) mdata%nvalue = mdata%nvalue + mdata%nlocvalue * mdata%grid_nx * mdata%grid_ny
   
    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Setup


  ! ***


  subroutine MAORI_Data_Satellite_Done( mdata, status )

!    use GO, only : goVarValue
!    use GO, only : ISTA_Done
    use Grid, only : Done
!    
!    use MAORI_LocList, only : MAORI_LocList_Done
    use MAORI_VarList, only : MAORI_VarList_Done
!#ifdef with_netcdf
!    use NetCDF       , only : NF90_Close
!#endif

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)       ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Done'
    
    ! --- local ----------------------------------
    
!    character(len=MAORI_LEN_NAME)   ::  obs_type
    integer                   ::  ipl

    ! --- begin -----------------------------------

    ! instrument dependend stuff:
    select case ( trim(mdata%instrument) )
    
!      ! ~ modis:
!      case ( 'modis' )
!    
!        ! loop over list files:
!        do ipl = 1, mdata%modis_nplatform
!          ! list file open ?
!          if ( mdata%modis_list_opened(ipl) ) then
!            ! close ...
!            close( mdata%modis_list_fu(ipl), iostat=status )
!            IF_NOTOK_RETURN(status=1)
!            ! reset flag:
!            mdata%modis_list_opened(ipl) = .false.
!          end if
!        end do
!        
!      ! ~ seviri:
!      case ( 'seviri' )
      
        ! nothing to be done ...
      
      ! unknown ...
      case default
        write (gol,'("unsupported instrument : ",a)') trim(mdata%instrument); call goErr
        TRACEBACK; status=1; return
    
    end select
    
    ! pixels stored ?
    if ( mdata%store_pix ) then
      ! clear:
      deallocate( mdata%pix_date      )
      deallocate( mdata%pix_time      )
      deallocate( mdata%pix_longitude )
      deallocate( mdata%pix_latitude  )
      deallocate( mdata%pix_proba     )
      !deallocate( mdata%pix_obs       )
      !deallocate( mdata%pix_obs_sigma )
      ! mapping ..
      if ( mdata%store_grid ) then
        deallocate( mdata%pix_grid_i )
        deallocate( mdata%pix_grid_j )
      end if
    end if

    ! grid stored ?
    if ( mdata%store_grid ) then
      ! clear:
      deallocate( mdata%grid_npix )
      deallocate( mdata%grid_w24 )
    end if

    ! done with grid:
    call Done( mdata%grid_lli, status )
    IF_NOTOK_RETURN(status=1)

!    ! done with sample averaging:
!    call ISTA_Done( mdata%ista, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! done with locations:
!    call MAORI_LocList_Done( mdata%loclist, status )
!    IF_NOTOK_RETURN(status=1)

    ! done with variables:
    call MAORI_VarList_Done( mdata%varlist, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear levels:
    deallocate( mdata%nlev )

    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
!      ! extract from query:
!      obs_type='none'
!        call goVarValue( mdata%obs_query, ';', 'type', '=', obs_type, status )
!        IF_ERROR_RETURN(status=1)
!      ! file input ?
!      select case ( trim(obs_type) )
!        case ( 'dummy' )
!          ! nothing extra required ...
!        case ( 'emep-daily-csv' )
!          ! setup access to observation data:
!          deallocate( mdata%obs_emep )
!        case default
!          write (gol,'("unsupported input type : ",a)') trim(obs_type); call goErr
!          TRACEBACK; status=1; return
!      end select
      ! clear:
      !deallocate( mdata%obs_iloc2iobs )
      if ( mdata%store_pix ) then
        deallocate( mdata%obs_pix_value  )
        deallocate( mdata%obs_pix_errstd )
        deallocate( mdata%obs_pix_astat  )
      end if
      if ( mdata%store_grid ) then
        deallocate( mdata%obs_grid_value  )
        deallocate( mdata%obs_grid_errstd )
        deallocate( mdata%obs_grid_astat  )
      end if
      deallocate( mdata%obs_ivar )
    end if

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Done


  ! ***


  subroutine MAORI_Data_Satellite_Inq( mdata, status, &
                              nparam, nlon, nlat, nloc, nvar, obs_nvar, &
                              assim_analyse, assim_analyse_now, assim_rho, &
                              nvalue, nlocvalue )
                              
    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(in)  ::  mdata
    integer, intent(out)                      ::  status

    integer, intent(out), optional            ::  nparam
    integer, intent(out), optional            ::  nlon, nlat, nloc
    integer, intent(out), optional            ::  nvar
    integer, intent(out), optional            ::  obs_nvar
    logical, intent(out), optional            ::  assim_analyse
    logical, intent(out), optional            ::  assim_analyse_now
    real, intent(out), optional               ::  assim_rho
    integer, intent(out), optional            ::  nvalue
    integer, intent(out), optional            ::  nlocvalue

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Inq'

    ! --- begin -----------------------------------

    ! number of required parameters:
    if ( present(nparam) ) nparam = npar

    ! return number of locations:
    if ( present(nlon) ) nlon = mdata%grid_nx
    if ( present(nlat) ) nlat = mdata%grid_ny
    if ( present(nloc) ) nloc = mdata%pix_n

    ! return number of variables to be put out:
    if ( present(nvar) ) nvar = mdata%nvar

    ! return number of observation variables:
    if ( present(obs_nvar) ) obs_nvar = mdata%obs_nvar

    ! analyse during assimilation ?
    if ( present(assim_analyse) ) assim_analyse = mdata%assim_analyse
    
    ! something to analyse now ?
    if ( present(assim_analyse_now) ) then
      ! true if this set is to be analysed, and if this is the end of the averaging period:
      assim_analyse_now = mdata%assim_analyse .and. mdata%end_of_period
    end if

    ! assimilation correlation scale:
    if ( present(assim_rho) ) assim_rho = mdata%assim_rho

    ! return number of elements in the state:
    if ( present(nvalue) ) nvalue = mdata%nvalue

    ! return number of elements in the state for a single location:
    if ( present(nlocvalue) ) nlocvalue = mdata%nlocvalue

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Inq


  ! ********************************************************************
  ! ***
  ! *** parameters
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Satellite_Param_Inq( mdata, ipar, type, name, unit, status )

    use MAORI_Param, only : MAORI_INT, MAORI_REAL, MAORI_CHAR
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(out)                      ::  type
    character(len=*), intent(out)             ::  name
    character(len=*), intent(out)             ::  unit
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Param_Inq'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bound_west      ) ; type = MAORI_REAL ; name = 'bound_west'      ; unit = 'degree_east'
      case ( ipar_bound_east      ) ; type = MAORI_REAL ; name = 'bound_east'      ; unit = 'degree_east'
      case ( ipar_bound_south     ) ; type = MAORI_REAL ; name = 'bound_south'     ; unit = 'degree_north'
      case ( ipar_bound_north     ) ; type = MAORI_REAL ; name = 'bound_north'     ; unit = 'degree_north'
      case ( ipar_grid_nlon       ) ; type = MAORI_INT  ; name = 'grid_nlon'       ; unit = '1'
      case ( ipar_grid_nlat       ) ; type = MAORI_INT  ; name = 'grid_nlat'       ; unit = '1'
      case ( ipar_nlayer          ) ; type = MAORI_INT  ; name = 'nlayer'          ; unit = '1'
      case ( ipar_bias_correction ) ; type = MAORI_CHAR ; name = 'bias_correction' ; unit = 'free'
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        write (gol,'("  set sample name : ",a)') trim(mdata%name); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Param_Inq


  ! ***
  

  subroutine MAORI_Data_Satellite_Param_Put_i( mdata, ipar, ival, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(in)                       ::  ival
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Param_Put_i'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_grid_nlon ) ; mdata%grid_nlon = ival
      case ( ipar_grid_nlat ) ; mdata%grid_nlat = ival
      case ( ipar_nlayer ) ; mdata%nlayer = ival
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Param_Put_i


  ! ***
  

  subroutine MAORI_Data_Satellite_Param_Put_r( mdata, ipar, rval, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    real, intent(in)                          ::  rval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Param_Put_r'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bound_west  ) ; mdata%bound_west  = rval
      case ( ipar_bound_east  ) ; mdata%bound_east  = rval
      case ( ipar_bound_south ) ; mdata%bound_south = rval
      case ( ipar_bound_north ) ; mdata%bound_north = rval
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Param_Put_r


  ! ***
  

  subroutine MAORI_Data_Satellite_Param_Put_s( mdata, ipar, sval, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)        ::  mdata
    integer, intent(in)                             ::  ipar
    character(len=*), intent(in)                    ::  sval
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Param_Put_s'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bias_correction ) ; mdata%bias_correction = sval
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Param_Put_s


  ! ********************************************************************
  ! ***
  ! *** locations
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Satellite_Loc_Inq( mdata, iloc, status, &
                                             lon, lat )

    use MAORI_Param  , only : MAORI_SAMPLE
    use MAORI_Param  , only : MAORI_TYPE_NAME
    use MAORI_LocList, only : MAORI_Loc_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(in)    ::  mdata
    integer, intent(in)                         ::  iloc
    integer, intent(out)                        ::  status

    real, intent(out), optional                 ::  lon
    real, intent(out), optional                 ::  lat

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Loc_Inq'

    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. mdata%store_pix ) then
      write (gol,'("could not inquire location if pixels are not stored")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract value:
    if ( present(lon) ) lon = mdata%pix_longitude(iloc)
    if ( present(lat) ) lat = mdata%pix_latitude (iloc)

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Loc_Inq


  ! ********************************************************************
  ! ***
  ! *** variables
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Satellite_Var_Inq( mdata, ivar, status, &
                                             name, unit, nlev )

    use MAORI_VarList, only : MAORI_Var_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(in)            ::  mdata
    integer, intent(in)                         ::  ivar
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  unit
    integer, intent(out), optional              ::  nlev

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Satellite_Var_Inq'

    ! --- begin -----------------------------------

    ! inquire output set:
    call MAORI_Var_Inq( mdata%varlist, ivar, status, &
                             name=name, unit=unit )
    IF_NOTOK_RETURN(status=1)

    ! return number of levels to be put out:
    if ( present(nlev) ) nlev = mdata%nlev(ivar)

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Var_Inq


  ! ********************************************************************
  ! ***
  ! *** observed values
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_Data_Satellite_Obs_Get( mdata, obs_ivar, status, iloc, ilon, ilat, y, r, alfa )

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(in)        ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status
    
    integer, intent(in), optional                   ::  iloc
    integer, intent(in), optional                   ::  ilon, ilat
    real, intent(out), optional                     ::  y       ! measured value
    real, intent(out), optional                     ::  r       ! error std.dev.
    real, intent(out), optional                     ::  alfa    ! screening factor

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Obs_Get'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( mdata%obs_nlev /= 1  ) then
      write (gol,'("getobs not implemented for obs_nlev ",i4)') mdata%obs_nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( present(iloc) ) then
      if ( any((/present(ilon),present(ilat)/)) ) then
        write (gol,'("either iloc or both ilon and ilat should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( .not. mdata%store_pix ) then
        write (gol,'("iloc present while pixels not stored")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( iloc > mdata%pix_n ) then
        write (gol,'("iloc ",i10," exceeds number of pixels ",i10)') iloc, mdata%pix_n; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(ilon) .or. present(ilat) ) then
      if ( .not. all((/present(ilon),present(ilat)/)) ) then
        write (gol,'("both ilon and ilat should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( .not. mdata%store_grid ) then
        write (gol,'("ilon,ilat present while grid not stored")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( (ilon > mdata%grid_nx) .or. (ilat > mdata%grid_ny) ) then
        write (gol,'("ilon,ilat ",2i10," exceeds grid size ",2i10)') ilon,ilat, mdata%grid_nx,mdata%grid_ny; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! copy measured value:
    if ( present(y) ) then
      ! fill:
      if ( present(iloc) ) y = mdata%obs_pix_value (1,iloc     ,obs_ivar)
      if ( present(ilon) ) y = mdata%obs_grid_value(ilon,ilat,1,obs_ivar)
    end if

    ! copy error std.dev.
    if ( present(r) ) then
      ! fill:
      if ( present(iloc) ) r = mdata%obs_pix_errstd (1,iloc     ,obs_ivar)
      if ( present(ilon) ) r = mdata%obs_grid_errstd(ilon,ilat,1,obs_ivar)
    end if

    ! screening factor:
    if ( present(alfa) ) then
      ! fill:
      alfa = mdata%assim_alfa
    end if

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Obs_Get


  ! ***


  subroutine MAORI_Data_Satellite_Obs_Put( mdata, obs_ivar, status, iloc, ilon, ilat, astat_ibset )

    ! --- in/out ---------------------------------

    type(T_MAORI_Data_Satellite), intent(inout)     ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status
    
    integer, intent(in), optional                   ::  iloc
    integer, intent(in), optional                   ::  ilon, ilat
    integer, intent(in), optional                   ::  astat_ibset  ! MAORI_ASTAT_*

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Satellite_Obs_Get'

    ! --- local -----------------------------------

    ! --- begin ----------------------------------
    
    ! check ...
    if ( present(iloc) ) then
      if ( any((/present(ilon),present(ilat)/)) ) then
        write (gol,'("either iloc or both ilon and ilat should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( .not. mdata%store_pix ) then
        write (gol,'("iloc present while pixels not stored")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( iloc > mdata%pix_n ) then
        write (gol,'("iloc ",i10," exceeds number of pixels ",i10)') iloc, mdata%pix_n; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(ilon) .or. present(ilat) ) then
      if ( .not. all((/present(ilon),present(ilat)/)) ) then
        write (gol,'("both ilon and ilat should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( .not. mdata%store_grid ) then
        write (gol,'("ilon,ilat present while grid not stored")'); call goErr
        TRACEBACK; status=1; return
      end if
      if ( (ilon > mdata%grid_nx) .or. (ilat > mdata%grid_ny) ) then
        write (gol,'("ilon,ilat ",2i10," exceeds grid size ",2i10)') ilon,ilat, mdata%grid_nx,mdata%grid_ny; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! set assimilation flag:
    if ( present(astat_ibset) ) then
      ! fill:
      if ( present(iloc) ) mdata%obs_pix_astat      (:,iloc,obs_ivar) = IBSet( mdata%obs_pix_astat      (:,iloc,obs_ivar), astat_ibset )
      if ( present(ilon) ) mdata%obs_grid_astat(ilon,ilat,:,obs_ivar) = IBSet( mdata%obs_grid_astat(ilon,ilat,:,obs_ivar), astat_ibset )
    end if

    ! ok
    status = 0

  end subroutine MAORI_Data_Satellite_Obs_Put
  

  ! ********************************************************************
  ! ***
  ! *** simulated values
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_State_Satellite_Init( mstate, mdata, name, status )

!    use GO  , only : AnyDate
!
!    use MAORI_Param, only : MAORI_LEN_LINE

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(out)       ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)         ::  mdata
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Init'

    ! --- local -----------------------------------
    
!    logical     ::  ldum

    ! --- begin -----------------------------------
    
    !write (gol,'("MAORI:     init sampling ...")'); call goPr
    
!    ! dummy to avoid warnings about unused variables:
!    ldum = mdata%started

    ! store name:
    mstate%name = trim(name)
    
    ! no values stored yet:
    if ( mdata%store_pix ) then
      nullify( mstate%pix_set  )
      nullify( mstate%pix_curr )
    end if
    if ( mdata%store_grid ) then
      nullify( mstate%grid_set  )
      nullify( mstate%grid_curr )
      nullify( mstate%grid_prev )
      nullify( mstate%grid_aver )
      nullify( mstate%grid_all )
    end if
    
    ! not started yet:
    mstate%started = .false.

    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Init


  ! ***


  subroutine MAORI_State_Satellite_Done( mstate, mdata, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)       ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Done'

    ! --- local -----------------------------------
    
!    logical     ::  ldum

    ! --- begin ----------------------------------
    
!    ! dummy to avoid warnings about unused variables:
!    ldum = mdata%started

    ! clear values:
    if ( mdata%store_pix ) then
      deallocate( mstate%pix_set  )
      deallocate( mstate%pix_curr )
    end if
    if ( mdata%store_grid ) then
      deallocate( mstate%grid_set  )
      deallocate( mstate%grid_curr )
      deallocate( mstate%grid_prev )
      deallocate( mstate%grid_aver )
      deallocate( mstate%grid_all )
    end if
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Done


  ! ***


  subroutine MAORI_State_Satellite_Start( mstate, mdata, status )

!    use GO, only : AnyDate
!    use MAORI_LocList, only : MAORI_LocList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)       ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Start'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. mdata%started ) then
      write (gol,'("could not start state while sample not started ...")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! setup storage for simulated values:
    if ( mdata%store_pix ) then
      ! storage:
      allocate( mstate%pix_set (mdata%maxlev,mdata%pix_max,mdata%nvar) )
      allocate( mstate%pix_curr(mdata%maxlev,mdata%pix_max,mdata%nvar) )
      ! initial value:
      mstate%pix_curr = 0.0
    end if
    if ( mdata%store_grid ) then
      ! storage:
      allocate( mstate%grid_set (mdata%grid_nx,mdata%grid_ny,mdata%maxlev,mdata%nvar) )
      allocate( mstate%grid_curr(mdata%grid_nx,mdata%grid_ny,mdata%maxlev,mdata%nvar) )
      allocate( mstate%grid_prev(mdata%grid_nx,mdata%grid_ny,mdata%maxlev,mdata%nvar) )
      allocate( mstate%grid_aver(mdata%grid_nx,mdata%grid_ny,mdata%maxlev,mdata%nvar) )
      allocate( mstate%grid_all (mdata%grid_nx,mdata%grid_ny,mdata%maxlev,0:24,mdata%nvar) )
      ! initial average:
      mstate%grid_curr = 0.0
      mstate%grid_prev = 0.0
      mstate%grid_aver = 0.0
      mstate%grid_all  = 0.0
    end if

    ! reset flag:
    mstate%started = .true.
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Start


  ! ***


  subroutine MAORI_State_Satellite_Setup( mstate, mdata, status )
    
    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)       ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Setup'

    ! --- local ----------------------------------
    
!    logical   ::  ldum
    
    ! --- begin ----------------------------------
    
!    ! dummy to avoid warnings about unused variables:
!    ldum = mdata%started
    
    ! store grids ?
    if ( mdata%store_grid ) then
      ! shift current values into previous:
      mstate%grid_prev = mstate%grid_curr
    end if

    ! reset current simulated values:
    if ( mdata%store_pix ) then
      mstate%pix_curr = 0.0
      mstate%pix_set  = .false.
    end if
    if ( mdata%store_grid ) then
      mstate%grid_curr = 0.0
      mstate%grid_set  = .false.
    end if

    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Setup


  ! ***


  subroutine MAORI_State_Satellite_Put( mstate, mdata, ivar, status, ilon, ilat, iloc, values )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)    ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)        ::  mdata
    integer, intent(in)                             ::  ivar
    integer, intent(out)                            ::  status

    integer, intent(in), optional                   ::  ilon, ilat, iloc
    real, intent(in), optional                      ::  values(:)
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Satellite_Put'

    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. present(values) ) then
      write (gol,'("optional arguments `values` should be provided ...")'); call goErr
      TRACEBACK; status=1; return 
    end if
    ! check ...
    if ( size(values) /= mdata%nlev(ivar) ) then
      write (gol,'("put values with ",i4," levels, while expected ",i4)') &
               size(values), mdata%nlev(ivar); call goErr
      TRACEBACK; status=1; return
    end if

    ! store value in appropriate levels:
    if ( any( (/present(iloc)/) ) ) then
      ! check ...
      if ( any( (/present(ilon),present(ilat)/) ) ) then
        write (gol,'("optional arguments `ilon` and `ilat` should not be provived with `iloc` ...")'); call goErr
        TRACEBACK; status=1; return 
      end if
      if ( .not. mdata%store_pix ) then
        write (gol,'("values on `iloc` provided while no storage on pixels required")'); call goErr
        TRACEBACK; status=1; return 
      end if
      ! store:
      mstate%pix_curr(1:mdata%nlev(ivar),iloc,ivar) = values
      ! set flag at all levels:
      mstate%pix_set(:,iloc,ivar) = .true.
    end if

    ! store value in appropriate levels:
    if ( any( (/present(ilon),present(ilat)/) ) ) then
      ! check ...
      if ( .not. all( (/present(ilon),present(ilat)/) ) ) then
        write (gol,'("optional arguments `ilon` and `ilat` should be both provided ...")'); call goErr
        TRACEBACK; status=1; return 
      end if
      if ( any( (/present(iloc)/) ) ) then
        write (gol,'("optional argument `iloc` should not be provived with `ilon` and `ilat` ...")'); call goErr
        TRACEBACK; status=1; return 
      end if
      if ( .not. mdata%store_grid ) then
        write (gol,'("values on `ilon` and `ilat` provided while no storage on grid required")'); call goErr
        TRACEBACK; status=1; return 
      end if
      ! store:
      mstate%grid_curr(ilon,ilat,1:mdata%nlev(ivar),ivar) = values
      ! set flag at all levels:
      mstate%grid_set(ilon,ilat,:,ivar) = .true.
    end if

    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Put
  

  ! ***
  
  
  subroutine MAORI_State_Satellite_Post( mstate, mdata, status )

    use GO, only : Get, wrtgol, Midnight, Precisely
    use MAORI_Param, only : MAORI_NAN

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)    ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)        ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Satellite_Post'

    ! --- local ----------------------------------
    
    integer       ::  ih
    integer       ::  ilev, ivar

    ! --- begin -----------------------------------
    
    ! store pixels ?
    if ( mdata%store_pix ) then
    
      ! any pixels yet ?
      if ( mdata%pix_n > 0 ) then

        ! check ...
        if ( .not. all(mstate%pix_set(:,1:mdata%pix_n,:)) ) then
          write (gol,'("not all simulated pixels seem to be set")'); call goErr
          TRACEBACK; status=1; return
        end if
        
      end if  ! pixels available
      
    end if
    
    ! store grids ?
    if ( mdata%store_grid ) then
    
      ! info ...
      call wrtgol('MAORI: post process `'//trim(mdata%name)//'` samples for : ',mdata%setup_tr(2)); call goPr

      ! clean history ?
      if ( mdata%new_period ) then
        ! info ...
        write (gol,'("MAORI:   clean history ...")'); call goPr
        !! check ...
        !if ( mdata%w_aver /= 0.0 ) then
        !  write (gol,'("new period but non-zero w_aver : ",f12.6)') mdata%w_aver; call goPr
        !  TRACEBACK; status=1; return
        !end if
        ! clear data:
        mstate%grid_aver = 0.0
        ! initial hour:
        ih = 0
        ! apply weights for first sample:
        do ivar = 1, mdata%nvar
          do ilev = 1, mdata%nlev(ivar)
            mstate%grid_aver(:,:,ilev,ivar) = mstate%grid_aver(:,:,ilev,ivar) &
                       + mstate%grid_prev(:,:,ilev,ivar) * mdata%grid_w24(:,:,ih)
          end do
        end do
      end if  ! new period

      ! hour in {0,..,24} :
      !  check:
      if ( .not. Precisely(mdata%setup_tr(2),1.0,'hour') ) then
        write (gol,'("Severi requires simulations on hourly times")'); call goErr
        TRACEBACK; status=1; return
      end if
      !  first guess:
      call Get( mdata%setup_tr(2), hour=ih )
      !  midnight at end of period ? then this is hour 24 :
      if ( Midnight(mdata%setup_tr(2)) .and. mdata%end_of_period ) ih = 24

      ! loop over variables:
      do ivar = 1, mdata%nvar
        ! loop over levels:
        do ilev = 1, mdata%nlev(ivar)
          ! apply weights for this hour:
          mstate%grid_aver(:,:,ilev,ivar) = mstate%grid_aver(:,:,ilev,ivar) &
                     + mstate%grid_curr(:,:,ilev,ivar) * mdata%grid_w24(:,:,ih)
          ! store all:
          mstate%grid_all(:,:,ilev,ih,ivar) = mstate%grid_curr(:,:,ilev,ivar)
        end do
      end do
      
    end if  ! store grid

    ! ok
    status = 0
    
  end subroutine MAORI_State_Satellite_Post
  

  ! ***


  subroutine MAORI_State_Satellite_Obs_Get( mstate, mdata, obs_ivar, status, iloc, ilon, ilat, value )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(in)       ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)        ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status

    integer, intent(in), optional                   ::  iloc
    integer, intent(in), optional                   ::  ilon, ilat
    real, intent(out), optional                     ::  value

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Satellite_Obs_Get'

    ! --- begin -----------------------------------
    
    ! return simulated value:
    if ( present(value) ) then
      ! check ...
      if ( mdata%obs_nlev /= 1  ) then
        write (gol,'("getobs not implemented for obs_nlev ",i4)') mdata%obs_nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      if ( present(iloc) ) then
        value = mstate%pix_curr(1,iloc,mdata%obs_ivar(obs_ivar))
      else if ( all((/present(ilon),present(ilat)/)) ) then
        value = mstate%grid_curr(ilon,ilat,1,mdata%obs_ivar(obs_ivar))
      else
        write (gol,'("either `iloc` or both `ilon` and `ilat` arguments requried")'); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Obs_Get


  ! ***


  subroutine MAORI_State_Satellite_Values_Get( mstate, mdata, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(in)           ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)            ::  mdata
    real, intent(out)                                   ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Values_Get'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nvalue ) then
      write (gol,'("size of output array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy pixel values:
      if ( mdata%store_pix .and. (mdata%pix_n > 0) ) then
        n = mdata%nlev(ivar) * mdata%pix_n
        values(i+1:i+n) = reshape(mstate%pix_curr(1:mdata%nlev(ivar),1:mdata%pix_n,ivar),(/n/)) ; i = i + n
      end if
      ! copy grid values:
      if ( mdata%store_grid ) then
        n = mdata%grid_nx * mdata%grid_ny * mdata%nlev(ivar)
        values(i+1:i+n) = reshape(mstate%grid_curr(:,:,1:mdata%nlev(ivar),ivar),(/n/)) ; i = i + n
        values(i+1:i+n) = reshape(mstate%grid_prev(:,:,1:mdata%nlev(ivar),ivar),(/n/)) ; i = i + n
        values(i+1:i+n) = reshape(mstate%grid_aver(:,:,1:mdata%nlev(ivar),ivar),(/n/)) ; i = i + n
      end if
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Values_Get


  ! ***


  subroutine MAORI_State_Satellite_Values_Put( mstate, mdata, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)        ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)            ::  mdata
    real, intent(in)                                    ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_Values_Put'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nvalue ) then
      write (gol,'("size of input array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy pixel values:
      if ( mdata%store_pix .and. (mdata%pix_n > 0) ) then
        n = mdata%nlev(ivar) * mdata%pix_n
        mstate%pix_curr(1:mdata%nlev(ivar),1:mdata%pix_n,ivar) = reshape(values(i+1:i+n),(/mdata%nlev(ivar),mdata%pix_n/)) ; i = i + n
      end if
      ! copy grid values:
      if ( mdata%store_grid ) then
        n = mdata%grid_nx * mdata%grid_ny * mdata%nlev(ivar)
        mstate%grid_curr(:,:,1:mdata%nlev(ivar),ivar) = reshape(values(i+1:i+n),(/mdata%grid_nx,mdata%grid_ny,mdata%nlev(ivar)/)) ; i = i + n
        mstate%grid_prev(:,:,1:mdata%nlev(ivar),ivar) = reshape(values(i+1:i+n),(/mdata%grid_nx,mdata%grid_ny,mdata%nlev(ivar)/)) ; i = i + n
        mstate%grid_aver(:,:,1:mdata%nlev(ivar),ivar) = reshape(values(i+1:i+n),(/mdata%grid_nx,mdata%grid_ny,mdata%nlev(ivar)/)) ; i = i + n
      end if
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_Values_Put


  ! ***


  subroutine MAORI_State_Satellite_LocValues_Get( mstate, mdata, ilon, ilat, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(in)           ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)            ::  mdata
    integer, intent(in)                                 ::  ilon, ilat
    real, intent(out)                                   ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_LocValues_Get'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nlocvalue ) then
      write (gol,'("size of output array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nlocvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( .not. mdata%store_grid ) then
      write (gol,'("assim only implemented for gridded fields yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nlev(ivar)
      values(i+1:i+n) = mstate%grid_curr(ilon,ilat,1:mdata%nlev(ivar),ivar) ; i = i + n
      values(i+1:i+n) = mstate%grid_prev(ilon,ilat,1:mdata%nlev(ivar),ivar) ; i = i + n
      values(i+1:i+n) = mstate%grid_aver(ilon,ilat,1:mdata%nlev(ivar),ivar) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_LocValues_Get


  ! ***


  subroutine MAORI_State_Satellite_LocValues_Put( mstate, mdata, ilon, ilat, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_State_Satellite), intent(inout)        ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)            ::  mdata
    integer, intent(in)                                 ::  ilon, ilat
    real, intent(in)                                    ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Satellite_LocValues_Put'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nlocvalue ) then
      write (gol,'("size of input array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nlocvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( .not. mdata%store_grid ) then
      write (gol,'("assim only implemented for gridded fields yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nlev(ivar)
      mstate%grid_curr(ilon,ilat,1:mdata%nlev(ivar),ivar) = values(i+1:i+n) ; i = i + n
      mstate%grid_prev(ilon,ilat,1:mdata%nlev(ivar),ivar) = values(i+1:i+n) ; i = i + n
      mstate%grid_aver(ilon,ilat,1:mdata%nlev(ivar),ivar) = values(i+1:i+n) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_State_Satellite_LocValues_Put


  ! ********************************************************************
  ! ***
  ! *** output
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_Output_Satellite_Init( moutput, mdata, name, status )

!    use GO, only : AnyDate

    ! --- in/out ---------------------------------

    type(T_MAORI_Output_Satellite), intent(out)       ::  moutput
    type(T_MAORI_Data_Satellite), intent(in)          ::  mdata
    character(len=*), intent(in)                      ::  name
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Satellite_Init'

    ! --- local -----------------------------------
    
    logical     ::  ldum
    
    ! --- begin -----------------------------------
    
    ! dummy to avoid warnings about unused variables:
    ldum = mdata%started
    
    ! store:
    moutput%name = trim(name)
    
    ! no file opened yet:
    moutput%opened = .false.

    ! ok
    status = 0

  end subroutine MAORI_Output_Satellite_Init
  

  ! ***
  
  
  subroutine MAORI_Output_Satellite_Done( moutput, mdata, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Output_Satellite), intent(inout)     ::  moutput
    type(T_MAORI_Data_Satellite), intent(in)          ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Satellite_Done'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! clear storage for netcdf id's :
    if ( mdata%store_pix ) then
      deallocate( moutput%varid_pix )
    end if
    if ( mdata%store_grid ) then
      deallocate( moutput%varid_grid )
    end if
    
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      if ( mdata%store_pix ) then
        deallocate( moutput%varid_pix_y )
        deallocate( moutput%varid_pix_r )
        deallocate( moutput%varid_pix_astat )
      end if
      if ( mdata%store_grid ) then
        deallocate( moutput%varid_grid_y )
        deallocate( moutput%varid_grid_r )
        deallocate( moutput%varid_grid_astat )
      end if
    end if
    
    ! ok
    status = 0

  end subroutine MAORI_Output_Satellite_Done
  

  ! ***
  
  
  subroutine MAORI_Output_Satellite_Start( moutput, mdata, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Output_Satellite), intent(inout)        ::  moutput
    type(T_MAORI_Data_Satellite), intent(in)             ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Satellite_Start'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! setup storage for netcdf id's :
    if ( mdata%store_pix ) then
      allocate( moutput%varid_pix(mdata%nvar) )
    end if
    if ( mdata%store_grid ) then
      allocate( moutput%varid_grid(mdata%nvar) )
      allocate( moutput%varid_grid_all(mdata%nvar) )
    end if

    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! setup storage for netcdf id's :
      if ( mdata%store_pix ) then
        allocate( moutput%varid_pix_y(mdata%obs_nvar) )
        allocate( moutput%varid_pix_r(mdata%obs_nvar) )
        allocate( moutput%varid_pix_astat(mdata%obs_nvar) )
      end if
      if ( mdata%store_grid ) then
        allocate( moutput%varid_grid_y(mdata%obs_nvar) )
        allocate( moutput%varid_grid_r(mdata%obs_nvar) )
        allocate( moutput%varid_grid_astat(mdata%obs_nvar) )
      end if
    end if

    ! no file opened yet:
    moutput%opened = .false.

    ! ok
    status = 0

  end subroutine MAORI_Output_Satellite_Start
  

  ! ***
  
  
  subroutine MAORI_Output_Satellite_Create( moutput, mdata, t, status )

    use GO, only : Get, Pretty
!    use GO, only : TDate, NewDate, IncrDate, AnyDate
!    use GO, only : operator(>), operator(+), operator(-)
    
    use MAORI_Param , only : MAORI_ASTAT_DESCRIPTION
    use MAORI_Common, only : com_path, com_model_name, com_exp_id, com_t0
    use MAORI_Common, only : MAORI_Common_NF90_GlobalAttributes
    
#ifdef with_netcdf
    use NetCDF, only : NF90_Create, NF90_Close
    use NetCDF, only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_NOCLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF, only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
#endif
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Output_Satellite), intent(inout)        ::  moutput
    type(T_MAORI_Data_Satellite), intent(in)             ::  mdata
    type(TDate), intent(in)                           ::  t
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Satellite_Create'

    ! --- local -----------------------------------
    
    integer                         ::  time6(6)
    character(len=MAORI_LEN_FILE)   ::  fname
    integer                         ::  cmode
    character(len=MAORI_LEN_LINE)   ::  units
    integer                         ::  ivar
    integer                         ::  dimid_z
    integer                         ::  varid
    
    ! --- begin -----------------------------------
    
    ! extract time fields from time assigned to current record:
    call Get( t, time6=time6 )

    ! new file name:
    write (fname,'(a,"_",a,"_",a,"_",i4.4,2i2.2)') &
              trim(com_model_name), trim(com_exp_id), &
              trim(mdata%name), time6(1:3)
    if ( len_trim(moutput%name) > 0 ) fname = trim(fname)//'_'//trim(moutput%name)
    fname = trim(fname)//'.nc'

    ! info ...
    write (gol,'("MAORI:     file name : ",a)') trim(fname); call goPr
    
    ! add path:
    fname = trim(com_path)//'/'//trim(fname)

#ifdef with_netcdf
    ! set creation mode flag:
    cmode = NF90_NOCLOBBER     ! do not overwrite existing files

    ! create file:
    status = NF90_Create( fname, cmode, moutput%ncid )
    if ( status /= 0 ) then
      write (gol,'("creating file :")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
#endif

    ! reset flag:
    moutput%opened = .true.

    ! write global attributes:
    call MAORI_Common_NF90_GlobalAttributes( moutput%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! extra global attributes:
    status = NF90_Put_Att( moutput%ncid, NF90_GLOBAL, 'bias_correction', trim(mdata%bias_correction) )
    IF_NF90_NOTOK_RETURN(status=1)

#ifdef with_netcdf
    ! define dimensions:

    ! gridded fields:
    if ( mdata%store_pix ) then
      status = NF90_Def_Dim( moutput%ncid, 'pixel', NF90_UNLIMITED, moutput%dimid_pix )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    ! gridded fields:
    if ( mdata%store_grid ) then
      status = NF90_Def_Dim( moutput%ncid, 'grid_nx', mdata%grid_nx, moutput%dimid_grid_nx )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( moutput%ncid, 'grid_ny', mdata%grid_ny, moutput%dimid_grid_ny )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( moutput%ncid, 'hhour', 1+24, moutput%dimid_hhour )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    ! output levels:
    status = NF90_Def_Dim( moutput%ncid, 'column', 1, moutput%dimid_column )
    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Def_Dim( moutput%ncid, 'surf_level', 1, moutput%dimid_slev )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Def_Dim( moutput%ncid, 'full_level', mdata%nlayer, moutput%dimid_flev )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Def_Dim( moutput%ncid, 'half_level', mdata%nlayer+1, moutput%dimid_hlev )
!    IF_NF90_NOTOK_RETURN(status=1)
!    ! timing:
!    status = NF90_Def_Dim( moutput%ncid, 'time', NF90_UNLIMITED, moutput%dimid_time )
!    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'datelen', 6, moutput%dimid_datelen )
    IF_NF90_NOTOK_RETURN(status=1)

    ! define pixel variables:

    if ( mdata%store_pix ) then

      status = NF90_Def_Var( moutput%ncid, 'pix_time', NF90_REAL, (/moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'time' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'time' )
      IF_NF90_NOTOK_RETURN(status=1)
      write (units,'("days since ",a)') trim(Pretty(com_t0))
      status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(units) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'calendar', 'standard' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_time = varid

      status = NF90_Def_Var( moutput%ncid, 'pix_date', NF90_INT, (/moutput%dimid_datelen,moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'date and time' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'year, month, day' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_date = varid

      status = NF90_Def_Var( moutput%ncid, 'pix_longitude', NF90_REAL, (/moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_east' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_lon = varid

      status = NF90_Def_Var( moutput%ncid, 'pix_latitude', NF90_REAL, (/moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_north' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_lat = varid

      status = NF90_Def_Var( moutput%ncid, 'pix_proba', NF90_REAL, (/moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_proba = varid

      status = NF90_Def_Var( moutput%ncid, 'pix_slots', NF90_INT, (/moutput%dimid_pix/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'bit flags for hours 1..24' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_pix_slots = varid

      if ( mdata%store_grid ) then
      
        status = NF90_Def_Var( moutput%ncid, 'pix_grid_i', NF90_INT, (/moutput%dimid_pix/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'date and time' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'units', 'year, month, day' )
        IF_NF90_NOTOK_RETURN(status=1)
        moutput%varid_pix_grid_i = varid
      
        status = NF90_Def_Var( moutput%ncid, 'pix_grid_j', NF90_INT, (/moutput%dimid_pix/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'date and time' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'units', 'year, month, day' )
        IF_NF90_NOTOK_RETURN(status=1)
        moutput%varid_pix_grid_j = varid
        
      end if  ! store grids
      
    end if  ! store pixels

    ! define grid variables:

    if ( mdata%store_grid ) then

      status = NF90_Def_Var( moutput%ncid, 'grid_longitude', NF90_REAL, (/moutput%dimid_grid_nx/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_east' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_grid_lon = varid

      status = NF90_Def_Var( moutput%ncid, 'grid_latitude', NF90_REAL, (/moutput%dimid_grid_ny/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_north' )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_grid_lat = varid

      status = NF90_Def_Var( moutput%ncid, 'grid_n', NF90_INT, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_grid_n = varid

      status = NF90_Def_Var( moutput%ncid, 'grid_w', NF90_REAL, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,moutput%dimid_hhour/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      moutput%varid_grid_w = varid
      
    end if  ! store grid

    ! loop over variables:
    do ivar = 1, mdata%nvar

      ! select level id:
      if ( mdata%nlev(ivar) == 1 ) then
        dimid_z = moutput%dimid_column
!      else if ( mdata%nlev(ivar) == mdata%nlayer ) then
!        dimid_z = moutput%dimid_flev
!      else if ( mdata%nlev(ivar) == mdata%nlayer+1 ) then
!        dimid_z = moutput%dimid_hlev
      else
        write (gol,'("unsupported number of levels : ")'); call goErr
        write (gol,'("  variable    : ",i6," (",a,")")') ivar, trim(mdata%varlist%var(ivar)%name); call goErr
        write (gol,'("  nlev        : ",i6)') mdata%nlev(ivar); call goErr
        write (gol,'("  nlayer      : ",i6)') mdata%nlayer; call goErr
        TRACEBACK; status=1; return
      end if

      ! store pixels ?
      if ( mdata%store_pix ) then
        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'pix_'//trim(mdata%varlist%var(ivar)%name), &
                     NF90_REAL, (/dimid_z,moutput%dimid_pix/), &
                     varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'model simulation' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_pix(ivar) = varid
      end if
      ! store grid ?
      if ( mdata%store_grid ) then
        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'grid_'//trim(mdata%varlist%var(ivar)%name), &
                     NF90_REAL, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,dimid_z/), &
                     varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'model simulation' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_grid(ivar) = varid
        
        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'grid_'//trim(mdata%varlist%var(ivar)%name)//'_all', &
                     NF90_REAL, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,dimid_z,moutput%dimid_hhour/), &
                     varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'model simulation' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_grid_all(ivar) = varid
      end if

    end do

    ! loop over observerd variables:
    do ivar = 1, mdata%obs_nvar

      ! select level id:
      if ( mdata%obs_nlev == 1 ) then
        dimid_z = moutput%dimid_column
      else
        write (gol,'("unsupported number of obs levels : ",i6)') mdata%obs_nlev; call goErr
        TRACEBACK; status=1; return
      end if

      ! store pixels ?
      if ( mdata%store_pix ) then

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'pix_'//trim(mdata%obs_varlist%var(ivar)%name)//'_y', &
                        NF90_REAL, (/dimid_z,moutput%dimid_pix/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_pix_y(ivar) = varid

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'pix_'//trim(mdata%obs_varlist%var(ivar)%name)//'_r', &
                        NF90_REAL, (/dimid_z,moutput%dimid_pix/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation error std.dev.' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_pix_r(ivar) = varid

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'pix_'//trim(mdata%obs_varlist%var(ivar)%name)//'_astat', &
                        NF90_INT, (/dimid_z,moutput%dimid_pix/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'assimilation status' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(MAORI_ASTAT_DESCRIPTION) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_pix_astat(ivar) = varid

      end if

      ! store pixels ?
      if ( mdata%store_grid ) then

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'grid_'//trim(mdata%obs_varlist%var(ivar)%name)//'_y', &
                        NF90_REAL, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,dimid_z/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_grid_y(ivar) = varid

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'grid_'//trim(mdata%obs_varlist%var(ivar)%name)//'_r', &
                        NF90_REAL, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,dimid_z/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation error std.dev.' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_grid_r(ivar) = varid

        ! define veriable:
        status = NF90_Def_Var( moutput%ncid, 'grid_'//trim(mdata%obs_varlist%var(ivar)%name)//'_astat', &
                        NF90_INT, (/moutput%dimid_grid_nx,moutput%dimid_grid_ny,dimid_z/), &
                        varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add attributes:
        status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'assimilation status' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(MAORI_ASTAT_DESCRIPTION) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        moutput%varid_grid_astat(ivar) = varid
      end if

      
    end do

    ! end defintion mode:

    status = NF90_EndDef( moutput%ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    ! no records written yet:
    if ( mdata%store_pix ) then
      moutput%pix_irec = 0
    end if

    ! ok
    status = 0

  end subroutine MAORI_Output_Satellite_Create
  

  ! ***
  
  
  subroutine MAORI_Output_Satellite_Write( moutput, mstate, mdata, status )

!    use GO, only : TDate, NewDate, AnyDate, Get
!    use GO, only : operator(>), operator(-)
!    use GO, only : rTotal
    
    use MAORI_Common, only : com_t0

#ifdef with_netcdf
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Close
#endif
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Output_Satellite), intent(inout)     ::  moutput
    type(T_MAORI_State_Satellite), intent(in)         ::  mstate
    type(T_MAORI_Data_Satellite), intent(in)          ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Satellite_Write'

    ! --- local -----------------------------------
    
!    integer                         ::  i, j
!    character(len=MAORI_LEN_NAME)   ::  sval
    integer                         ::  ivar
!    type(TDate)                     ::  t
!    integer                         ::  time6(6)
!    real                            ::  time
!    integer                         ::  ihist
    
    ! --- begin -----------------------------------
    
    ! info ...
    write (gol,'("MAORI: write set ",a," `",a,"` ...")') trim(mdata%name), trim(moutput%name); call goPr
      
    ! write if:
    !  o end of collection period is reached;
    !  o end of the run is reached:
    if ( mdata%end_of_period .or. mdata%the_end ) then
      ! info ...
      write (gol,'("MAORI:   write now ...")'); call goPr
    else
      ! info ...
      write (gol,'("MAORI:   do not write yet: end_of_period ",l2)') mdata%end_of_period; call goPr
      ! ok
      status=0; return
    end if

    ! new file ?
    if ( .not. moutput%opened ) then
      ! info ...
      write (gol,'("MAORI:   create new file ...")'); call goPr
      ! new file:
      call MAORI_Output_Satellite_Create( moutput, mdata, mdata%period_tr(1), status )
      IF_NOTOK_RETURN(status=1)
    end if
    
!    ! info ...
!    write (gol,'("MAORI:   write records; ihist is ",i4," of ",i4)') mdata%ihist, mdata%nhist; call goPr

#ifndef with_netcdf
    ! check libraries ...
    write (gol,'("not compiled with netcdf support")'); call goErr
    TRACEBACK; status=1; return
#endif

!    ! loop over historical records:
!    do ihist = 1, mdata%nhist
!
!      ! next time record:
!      moutput%pix_irec = moutput%pix_irec + 1
!
!      ! write station info only once ...    
!      if ( moutput%pix_irec == 1 ) then
!#ifdef with_netcdf
!        ! write station info; 
!        ! for some reason, character arrays need to be written one-by-one:
!        do i = 1, mdata%nloc
!          sval = trim(mdata%loclist%loc(i)%name)
!          status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_name, trim(sval), &
!                                   start=(/1,i/), count=(/len_trim(sval),1/) )
!          IF_NF90_NOTOK_RETURN(status=1)
!        end do
!        do i = 1, mdata%nloc
!          sval = trim(mdata%loclist%loc(i)%code)
!          status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_code, trim(sval), &
!                                   start=(/1,i/), count=(/len_trim(sval),1/) )
!          IF_NF90_NOTOK_RETURN(status=1)
!        end do
!        ! write non-character fields:
!        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_lon, mdata%loclist%loc(1:mdata%nloc)%lon )
!        IF_NF90_NOTOK_RETURN(status=1)
!        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_lat, mdata%loclist%loc(1:mdata%nloc)%lat )
!        IF_NF90_NOTOK_RETURN(status=1)
!        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_alt, mdata%loclist%loc(1:mdata%nloc)%alt )
!        IF_NF90_NOTOK_RETURN(status=1)
!#endif
!      end if  ! first record

    ! info ...
    write (gol,'("MAORI:   add pixels ",i10," .. ",i10)') moutput%pix_irec+1, moutput%pix_irec+mdata%pix_n; call goPr

#ifdef with_netcdf
    ! pixels ?
    if ( mdata%store_pix ) then
      ! any pixels ?
      if ( mdata%pix_n > 0 ) then
        ! add pixels:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_lon, &
                       mdata%pix_longitude(1:mdata%pix_n), &
                       start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add pixels:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_lat, &
                       mdata%pix_latitude(1:mdata%pix_n), &
                       start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add pixels for proba:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_proba, &
                       mdata%pix_proba(1:mdata%pix_n), &
                       start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add pixels:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_slots, &
                       mdata%pix_slots(1:mdata%pix_n), &
                       start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! also grids ?
        if ( mdata%store_grid ) then
          ! add pixels:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_grid_i, &
                         mdata%pix_grid_i(1:mdata%pix_n), &
                         start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! add pixels:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_grid_j, &
                         mdata%pix_grid_j(1:mdata%pix_n), &
                         start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        ! add pixels:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_time, &
                       mdata%pix_time(1:mdata%pix_n), &
                       start=(/moutput%pix_irec+1/), count=(/mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add pixels:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_date, &
                       mdata%pix_date(:,1:mdata%pix_n), &
                       start=(/1,moutput%pix_irec+1/), count=(/6,mdata%pix_n/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if  ! any pixels
    end if  ! store pixels
    ! gridded ?
    if ( mdata%store_grid ) then
      ! add ax:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_lon, mdata%grid_lli%lon_deg )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add ax:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_lat, mdata%grid_lli%lat_deg )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add field:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_n, mdata%grid_npix )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add field:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_w, mdata%grid_w24 )
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! store pixels ?
      if ( mdata%store_pix ) then
        ! any pixels ?
        if ( mdata%pix_n > 0 ) then
          ! write simulations; order should be (lev,pixel):
          status = NF90_Put_Var( moutput%ncid, moutput%varid_pix(ivar), &
                         mstate%pix_curr(1:mdata%nlev(ivar),1:mdata%pix_n,ivar), &
                         start=(/1,moutput%pix_irec+1/), &
                         count=(/mdata%nlev(ivar),mdata%pix_n/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if  ! any pixels
      end if
      ! store grid ?
      if ( mdata%store_grid ) then
        ! write simulations:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_grid(ivar), &
                       mstate%grid_aver(:,:,1:mdata%nlev(ivar),ivar) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! write simulations:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_all(ivar), &
                       mstate%grid_all(:,:,1:mdata%nlev(ivar),:,ivar) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
    end do

    ! write observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! loop over observerd variables:
      do ivar = 1, mdata%obs_nvar
        ! store pixels ?
        if ( mdata%store_pix ) then
          ! any pixels ?
          if ( mdata%pix_n > 0 ) then
            ! write observations; order should be (lev,pixel):
            status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_y(ivar), &
                           mdata%obs_pix_value (1:mdata%obs_nlev,1:mdata%pix_n,ivar), &
                           start=(/1,moutput%pix_irec+1/), &
                           count=(/mdata%nlev(ivar),mdata%pix_n/) )
            IF_NF90_NOTOK_RETURN(status=1)
            ! write error standard deviations:
            status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_r(ivar), &
                           mdata%obs_pix_errstd(1:mdata%obs_nlev,1:mdata%pix_n,ivar), &
                           start=(/1,moutput%pix_irec+1/), &
                           count=(/mdata%nlev(ivar),mdata%pix_n/) )
            IF_NF90_NOTOK_RETURN(status=1)
            ! write assimilation status:
            status = NF90_Put_Var( moutput%ncid, moutput%varid_pix_astat(ivar), &
                           mdata%obs_pix_astat(1:mdata%obs_nlev,1:mdata%pix_n,ivar), &
                           start=(/1,moutput%pix_irec+1/), &
                           count=(/mdata%nlev(ivar),mdata%pix_n/) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if  ! any pixels
        end if  ! store pixels
        ! store grid ?
        if ( mdata%store_grid ) then
          ! write observations; order should be (lev,pixel):
          status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_y(ivar), &
                         mdata%obs_grid_value (1:mdata%grid_nx,1:mdata%grid_ny,1:mdata%obs_nlev,ivar) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write error standard deviations:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_r(ivar), &
                         mdata%obs_grid_errstd(1:mdata%grid_nx,1:mdata%grid_ny,1:mdata%obs_nlev,ivar) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write assimilation status:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_grid_astat(ivar), &
                         mdata%obs_grid_astat(1:mdata%grid_nx,1:mdata%grid_ny,1:mdata%obs_nlev,ivar) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    ! increase record counters:
    ! pixels records ?
    if ( mdata%store_pix ) then
      ! any pixels ?
      if ( mdata%pix_n > 0 ) then
        moutput%pix_irec = moutput%pix_irec + mdata%pix_n
      end if  ! any pixels
    end if  ! store pixels
#endif

!    end do  ! hist

    ! close if end of period or end of run is reached:
    if ( mdata%end_of_period .or. mdata%the_end ) then
      ! info ...
      write (gol,'("MAORI:   close file ...")'); call goPr
      ! close:
#ifdef with_netcdf
      status = NF90_Close( moutput%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif
      ! reset flag:
      moutput%opened = .false.
    end if

    ! ok
    status = 0

  end subroutine MAORI_Output_Satellite_Write



end module MAORI_Set_Satellite
