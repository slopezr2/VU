!#################################################################
!
! NAME
!   LEKF_Meas_OMI_TRC  -  interface to OMI TRC data
!      latest adjustments:
!         7 june 2010, Henk Eskes, KNMI
!
! Interface:
!   call InitMeas_OMI_TRC( rcfile, status )
!       Set relevant parameters (based on rc file) and initialise
!   call DoneMeas_OMI_TRC( status )
! Finish and release storage if needed
!   call GetMeas_OMI_TRC( tint1, tint2, status )
! Make measurements between time "tint1" and "tint2" available to assimilation (read the OMI file when needed)
!   call CalcMeas_OMI_TRC( le_vmr_grid, le_vmr_aloft, le2omi_vcd_trop, status )
! Simulate OMI TRC obs "le2omi_vcd_trop" given the model concentrations "le_vmr_grid" and "le_vmr_aloft"
!       (only for observations in time range set with GetMeas_OMI_TRC)
!   call MeasUpdate_OMI_TRC( t1, t2, nmodes, status )
!       Analyse observations between times "t1" and "t2" for "nmodes" modes
!
! Output routines:
!   call Init( kfmo, rcfile, rckey, status )
! Initialise OMI output
!   call Done( kfmo, status )
! Finalise
!   call PutOut( kfmo, key, t, status )
! Write netcdf file with OMI data and model equivalent values
! Routine should be called every timestep
! Output written only once per day
! Note: PutOut will call "Calcmeas_OMI_TRC" !
!
! OMI TRC DATA
!
!  Production chain:
!   1. Global OMI files (hdf4/5), (very large)
!   2. are converted to smaller files with
!      selected pixels and fields for Europe, hdf format (Henk Eskes )
!   3. these are automatically converted to netcdf (Henk Eskes)
!
!  Example netcdf file:
!
!    netcdf omi_trc_EU_20060701 {
!    dimensions:
!            fakeDim0 = 17814 ;   nmeas
!            fakeDim1 = 6     ;   year,month,day,hour,minu,sec
!            fakeDim2 = 17814 ;   nmeas
!            fakeDim3 = 17814 ;   nmeas
!            fakeDim4 = 17814 ;   nmeas
!            fakeDim5 = 4     ;   corners
!            fakeDim6 = 17814 ;   nmeas
!            fakeDim7 = 4     ;   corners
!            fakeDim8 = 17814 ;   nmeas
!            fakeDim9 = 17814 ;   nmeas
!            fakeDim10 = 17814 ;   nmeas
!            fakeDim11 = 17814 ;   nmeas
!            fakeDim12 = 30    ;   levels
!            fakeDim13 = 17814 ;   nmeas
!            fakeDim14 = 31    ;   half levels
!            fakeDim15 = 17814 ;   nmeas
!            fakeDim16 = 17814 ;   nmeas
!            fakeDim17 = 17814 ;   nmeas
!            fakeDim18 = 17814 ;   nmeas
!            fakeDim19 = 17814 ;   nmeas
!    variables:
!            int date_time(fakeDim0, fakeDim1) ;
!            float longitude(fakeDim2) ;
!            float latitude(fakeDim3) ;
!            float corner_longitudes(fakeDim4, fakeDim5) ;
!            float corner_latitudes(fakeDim6, fakeDim7) ;
!            float vcd_trop(fakeDim8)                     ; 1e15 (mlc TRC)/cm2
!            float sigma_vcd_trop(fakeDim9) ;
!            float sigma_vcd_trop_ak(fakeDim10) ;
!            float kernel(fakeDim11, fakeDim12) ;
!            float pressure_levels(fakeDim13, fakeDim14) ;
!            float cloud_top_pessure(fakeDim15) ;
!            float cloud_radiance_fraction(fakeDim16) ;
!            float cloud_radiance_fraction(fakeDim17) ;
!            float pixel_number(fakeDim18) ;
!            float image_number(fakeDim19) ;
!
!    // global attributes:
!                    :Author = "Henk Eskes" ;
!                    :Affiliation = "KNMI (Royal Netherlands Meteorological Institute)" ;
!                    :Email = "eskes@knmi.nl" ;
!                    :Instrument = "OMI" ;
!                    :Species = "TRC" ;
!                    :OMI_retrieval_software = "tm4trca_omi: version 0.9.3.5, 15 April 2007" ;
!                    :Unit_of_TRC_column = "1e15 molecules/cm2" ;
!                    :Unit_of_pressure = "Pascal" ;
!                    :Number_of_measurements = 17814 ;
!                    :Number_of_pressure_levels = 30 ;
!                    :Total_number_of_measurements_on_this_day = 985900 ;
!                    :Window_longitude_minimum = -15.f ;
!                    :Window_longitude_maximum = 35.f ;
!                    :Window_latitude_minimum = 35.f ;
!                    :Window_latitude_maximum = 70.f ;
!                    :Extra_filter = "Cloud radiance fraction < 0.5" ;
!                    :Track_id_list = "   10436" ;
!                    :Track_count = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
!                    :Track_orbit_ident = "0607010050" ;
!                    :Track_input_pointer = "OMI-Aura_L2-OMTRC_2006m0701t0023-o10423_v002-2006m0704t202304.he5" ;
!                    :Track_stripe_corr = "Polynomial" ;
!                    :Track_orbit_number = 10423, 10424, 10425, 10426, 10427, 10428, 10429, 10430, 10431, 10432, 10433, 10434, 10435, 10436 ;
!                    :Track_start_date = "20060701005047...";
!                    :Track_end_date = "20060701013625...";
!    }
!
!
!  KERNEL
!
!    yr =  ya +  A(H'x-ya) + v   ,   v  ~ N(o, R  )
!   Cyr = Cya + CA(H'x-ya) + vc  ,   vc ~ N(o,CRC')   total column, C=[1,..,1]
!         ya = 0 in our case (linear approximation) and is not provided in data file
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
!###############################################################################

module LEKF_Meas_OMI_TRC

  use GO,         only : gol, goPr, goErr
  use GO,         only : TDate

  use NetCDF,     only : NF90_StrError, NF90_NOERR

  use Dims            , only : nx, ny
  use LE_Output_Common, only : T_LE_Output_Common


  implicit none


  ! --- in/out ----------------------------

  private

  public  ::  Initmeas_OMI_TRC, Donemeas_OMI_TRC
  public  ::  Getmeas_OMI_TRC
  public  ::  CalcMeas_OMI_TRC

  public  ::  Meas_Output_OMI_TRC

  public  ::  omi_trc
  public  ::  omi_trc_itracer


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_OMI_TRC'

  ! extra layers on top of LE:
  !integer, parameter     ::  nhz = 2   ! free trop, strato
  integer, parameter     ::  nhz = 0   ! profile from boundary conditions

  !! define two additional layers on top of the LE model
  !! first layer, between top LE (3.5 km) and topFreeTropLayer will contain the free tropospheric column
  !real, parameter  ::  topFreeTropLayer = 7.0e3     ! meter
  !real, parameter  ::  topAtmosphere    = 200.0e3   ! meter

  ! which tracer entities simulated from state should be put out ?
  integer, parameter            ::  nent = 4
  integer, parameter            ::  ient_le2omi_vcd_trop = 1
  integer, parameter            ::  ient_le2omi_vcd      = 2
  integer, parameter            ::  ient_le_vcd          = 3
  integer, parameter            ::  ient_le_vmr          = 4
  character(len=15), parameter   ::  entname(nent) = &
       (/'le2omi_vcd_trop','le2omi_vcd     ','le_vcd         ','le_vmr         '/)

  ! filter states to be written:
  !integer, parameter            ::  nfs = 2
  !character(len=2), parameter   ::  fsname(nfs) = (/'xb','x '/)
  integer, parameter            ::  nfs = 3
  character(len=2), parameter   ::  fsname(nfs) = (/'xb','x ','s '/)

  ! filter moments: 1=analysis, 2=forecast
  integer, parameter            ::  nfm = 2
  integer, parameter            ::  ifm_forecast = 1
  integer, parameter            ::  ifm_analysis = 2
  character(len=*), parameter   ::  fm_descr(nfm) = &
       (/ 'after filter forecast', &
       'after filter analysis' /)

  ! Status number is sum of one or more numbers 2^n,
  ! for example   1+4+8
  ! To test if it is composed of a certain number, use:
  !   if ( iand(status,4) /= 0 )  ..
  !
  integer, parameter            ::  status_default           = 0
  ! not used, we skip pixels outside the domain immediately ...
  !integer, parameter            ::  status_outside           = 1
  integer, parameter            ::  status_nodata            = 2
  integer, parameter            ::  status_validation        = 4
  character(len=*), parameter   ::  status_description = &
       'status flag, 0=default, +2=no-data, +4=validation'
       !'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation'

  ! --- types --------------------------------

  ! output stuff

  type Meas_Output_OMI_TRC
     ! common stuff:
     type(T_LE_Output_Common)    ::  com
     ! file opened ?
     logical                     ::  opened
     ! current time range:
     type(TDate)                 ::  tr(2)
     ! latest time:
     type(TDate)                 ::  tprev
     ! time record counter:
     integer                     ::  itrec
     ! file name:
     character(len=512)          ::  fname
     ! replace existing ?
     logical                     ::  replace
     ! file handle:
     integer                     ::  ncid
     ! dimension handles:
     integer                     ::  dimid_omi_obs
     integer                     ::  dimid_omi_lev
     integer                     ::  dimid_omi_hlev
     integer                     ::  dimid_le_lev
     integer                     ::  dimid_le_hlev
     integer                     ::  dimid_time
     integer                     ::  dimid_datelen
     ! variable handles:
     integer                     ::  varid_pixel_id
     integer                     ::  varid_lon
     integer                     ::  varid_lat
!     integer                     ::  varid_time
!     integer                     ::  varid_date
     integer                     ::  varid_omi_time
     integer                     ::  varid_omi_date
     integer                     ::  varid_omi_pressure_levels
     integer                     ::  varid_omi_kernel
     integer                     ::  varid_omi_cloud_radiance_fraction
     integer                     ::  varid_omi_status
     integer                     ::  varid_omi_analysed
     integer                     ::  varid_omi_screened
     integer                     ::  varid_le_altitude_levels
     integer                     ::  varid_le_pressure_levels
     ! tracer variables:
     !real                        ::  unitconv
     integer                     ::  varid_y
     integer                     ::  varid_r
     integer                     ::  varid(nent,nfs,nfm)
  contains
    procedure   ::  Init        =>  meas_output_omi_Init
    procedure   ::  Done        =>  meas_output_omi_Done
    procedure   ::  PutOut      =>  meas_output_omi_PutOut
  end type Meas_Output_OMI_TRC


  ! --- var ----------------------------------

  ! archive description:
  character(len=256)        ::  arch_filenames
  character(len=256)        ::  arch_file

  ! Storage to read the file with 1 day of measurements
  type OMI_TRC_type
    ! tracer name:
    character(len=32)    ::  tracer_name
    character(len=32)    ::  top_ph_entry             ! input parameter for pressure (to merge bc from top )
    character(len=32)    ::  top_trc_entry            ! input parameter for tracers (to merge bc from top )
    ! pixels:
    integer              ::  nmeas                    ! nr of OMI obs on this day
    integer              ::  nmeas_intimestep         ! nr of valid OMI obs in current assim timestep
    Integer              ::  nlev                     ! nr of pressure/kernel levels
    logical              ::  nodataleft               ! Are there OMI obs after this time step ?
    integer, pointer     ::  date(:,:)                ! date of OMI obs
    real, pointer        ::  time(:)                  ! time of OMI obs ("days since ...")
    real, pointer        ::  longitude(:)
    real, pointer        ::  latitude(:)
    real, pointer        ::  corner_longitudes(:,:)
    real, pointer        ::  corner_latitudes (:,:)
    real, pointer        ::  vcd_trop(:)              ! OMI tropospheric column
    real, pointer        ::  sigma_vcd_trop(:)        ! and uncertainty
    real, pointer        ::  kernel(:,:)              ! Averaging kernel vector
    real, pointer        ::  pressure_levels(:,:)     ! Pressure levels on which kernel is defined
    real, pointer        ::  cloud_radiance_fraction(:)        ! Cloud fraction of OMI obs
    !
    ! pixel number within swath:
    integer, pointer     ::  pixel_number(:)          ! Pixel number of OMI obs.
    integer              ::  pixel_range(2)           ! start end
    !
    ! unique pixel number, used in domain decomposition:
    !integer, allocatable ::  pixel_id(:)
    ! according to JNt, "allocatable" gives strange errors on KNMI/HPC ...
    integer, pointer     ::  pixel_id(:)
    !
    integer, pointer     ::  ix(:)                    ! index of LE grid cell containing OMI obs
    integer, pointer     ::  iy(:)                    ! index of LE grid cell containing OMI obs
    integer, pointer     ::  status(:)                ! flag for OMI obs - see above
    logical, pointer     ::  analysed(:)              ! Is this obs analysed ?
    logical, pointer     ::  screened(:)              ! Is this obs rejected during the analysis ?
    logical, pointer     ::  intimestep(:)            ! Is this obs part of current Kalman time step ?
    ! le_date ?
    real, pointer        ::  le_altitude_levels(:,:)  ! LE altitudes @ OMI obs
    real, pointer        ::  le_pressure_levels(:,:)  ! LE pressures @ OMI obs
    real, pointer        ::  le_trc_vmr_prof(:,:)     ! LE vmr profile @ OMI obs (on LE levels)
    real, pointer        ::  le_trc_vcd_prof(:,:)     ! LE vcd profile @ OMI obs (on LE levels)
    real, pointer        ::  le2omi_trc_vcd_prof(:,:) ! LE vcd profile @ OMI obs (on OMI levels)
    real, pointer        ::  le2omi_trc_vcd_trop(:,:,:)   ! LE estimate of OMI obs (tropospheric column)
    !
    ! analyse this set ?
    logical              ::  analyse
    !
    !! TESTING: analysis pixel order:
    !character(len=32)    ::  analyse_order
    !
    ! how to define LE-OMI repr.error. ?
    !  'data'   : from file, eventually with extra scaling
    !  'frac'   : as fraction of observed value
    character(len=4)        ::  r_type
    ! the scaling factor for the OMI observation error
    real                    ::  r_data_scaling
    ! EnerGEO special:  r = frac * y, and bounded to range
    real                    ::  r_frac_factor
    real                    ::  r_frac_min
    real                    ::  r_frac_max
    !
    ! screening factor: measurements are rejected
    ! if square of departure exceeds factor times variance:
    real                 ::  screening_factor
    !
    ! the radius of influence for TRC values
    ! "wild guess" to be sophisticated in future
    real                 ::  rho
    !
  end type OMI_TRC_type

  ! Subset of observations for Kalman time step
  type(OMI_TRC_type)      ::  omi_trc

  !! Storage to read the file with 1 day of measurements
  !type TRC_col_freetrop_type
  !   integer                       ::  year, month, day       ! day of simulation
  !   real, pointer, dimension(:,:) ::  trc_col_freetrop       ! estimate/simulation of free troposphere TRC column @ OMI overpass
  !end type TRC_col_freetrop_type

  !! Free troposphere column estimate
  !type(TRC_col_freetrop_type)  ::  nft

  ! species analysis weights:
  logical, allocatable    ::  analyse_spec(:)  ! (nnoise)
  logical                 ::  analyse_sia

  ! noise analysis weights:
  logical, allocatable    ::  analyse_noise(:)  ! (nnoise)

  ! time range:
  type(TDate)             ::  omi_trc_tr(2)  ! time range with valid data (24 hour)

  ! top of troposphere to be used for kernels:
  real                    ::  omi_trc_ptropo

  ! LE species index:
  integer                 ::  omi_trc_itracer


contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================


  ! opens the measurement file (if present)
  ! end sets some parameters

  subroutine InitMeas_OMI_TRC( rcF, status )

    use GO             , only : TrcFile
    use GO             , only : TDate, NewDate, Get, AnyDate
    use GO             , only : goMatchValue
    use GO             , only : goReadFromLine
    use Indices        , only : nspec_all
    use Indices        , only : specname
    use LEKF_Meas_Tools, only : GetSpecApply
    use LEKF_Noise     , only : nnoise, GetNoiseApply

    ! --- in/out --------------------------------

    type(TRcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/initmeas_OMI_TRC'

    ! --- local ----------------------------------

    character(len=1024)   ::  key
    character(len=1024)   ::  line

    ! --- begin ----------------------------------

    write (gol,'("KF:     OMI TRC - setup measurements ...")'); call goPr

    ! ~~ settings from .rc file

    ! tracer:
    call rcF%Get( 'kf.meas.omi_trc.tracer', omi_trc%tracer_name, status )
    IF_NOTOK_RETURN(status=1)
    ! search index:
    call goMatchValue( trim(omi_trc%tracer_name), specname, omi_trc_itracer, status )
    IF_NOTOK_RETURN(status=1)

    ! filename template:
    call rcF%Get( 'kf.meas.omi_trc.filenames', arch_filenames, status )
    IF_NOTOK_RETURN(status=1)

    call rcF%Get( 'kf.meas.omi_trc.analyse', omi_trc%analyse, status )
    IF_NOTOK_RETURN(status=1)

    !! testing ....
    !call rcF%Get( 'kf.meas.omi_trc.analyse.order', omi_trc%analyse_order, status )
    !IF_NOTOK_RETURN(status=1)

    ! localization length scale:
    call rcF%Get( 'kf.meas.omi_trc.rho', omi_trc%rho, status )
    IF_NOTOK_RETURN(status=1)

    ! how to define representation error?
    call rcF%Get( 'kf.meas.omi_trc.r.type', omi_trc%r_type, status )
    IF_NOTOK_RETURN(status=1)
    ! switch:
    select case ( trim(omi_trc%r_type) )
      !~ from data?
      case ( 'data' )
        ! extra scaling factor:
        call rcF%Get( 'kf.meas.omi_trc.r.data.scaling', omi_trc%r_data_scaling, status, default=1.0 )
        IF_NOTOK_RETURN(status=1)
      !~ fraction of data:
      case ( 'frac' )
        ! scaling and bounds:
        call rcF%Get( 'kf.meas.omi_trc.r.frac.factor', omi_trc%r_frac_factor, status )
        IF_NOTOK_RETURN(status=1)
        call rcF%Get( 'kf.meas.omi_trc.r.frac.min'   , omi_trc%r_frac_min   , status )
        IF_NOTOK_RETURN(status=1)
        call rcF%Get( 'kf.meas.omi_trc.r.frac.max'   , omi_trc%r_frac_max   , status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported OMI TRC r type `",a,"`")') trim(omi_trc%r_type); call goErr
        TRACEBACK; status=1; return
    end select

    ! screening factor:
    call rcF%Get( 'kf.meas.omi_trc.alfa', omi_trc%screening_factor, status )
    IF_NOTOK_RETURN(status=1)

    ! top of tropo sphere to be used with kernel (Pa):
    call rcF%Get( 'kf.meas.omi_trc.ptropo', omi_trc_ptropo, status )
    IF_NOTOK_RETURN(status=1)

    ! boundary condition variables used above model top:
    call rcF%Get( 'kf.meas.omi_trc.top_hp', omi_trc%top_ph_entry, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( 'kf.meas.omi_trc.top_trc', omi_trc%top_trc_entry, status )
    IF_NOTOK_RETURN(status=1)

    ! pixel range:
    call rcF%Get( 'kf.meas.omi_trc.pixels', line, status )
    IF_NOTOK_RETURN(status=1)
    ! defined?
    if ( len_trim(line) > 0 ) then
      call goReadFromLine( line, omi_trc%pixel_range(1), status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( line, omi_trc%pixel_range(2), status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
    else
      omi_trc%pixel_range(1) = -999
      omi_trc%pixel_range(2) = -999
    end if

    ! analysis settings neeeded ?
    if ( omi_trc%analyse ) then

      ! read description line for correlations with species
      call rcF%Get( 'kf.meas.omi_trc.spec', key, status )
      IF_NOTOK_RETURN(status=1)
      ! storage for noise correlation weights:
      allocate( analyse_spec(nspec_all), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! fill:
      call GetSpecApply( key, analyse_spec, analyse_sia, status )
      IF_NOTOK_RETURN(status=1)

      ! read description line for correlations with noise factors:
      call rcF%Get( 'kf.meas.omi_trc.noise', key, status )
      IF_NOTOK_RETURN(status=1)
      ! storage for noise correlation weights:
      allocate( analyse_noise(nnoise), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! fill:
      call GetNoiseApply( key, analyse_noise, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! nothing stored yet ...
    omi_trc%nmeas = 0

    ! no data read yet:
    omi_trc_tr(1) = AnyDate()
    omi_trc_tr(2) = AnyDate()

    ! ok
    status = 0

  end subroutine InitMeas_OMI_TRC


  ! ***


  ! read station locations etc

  subroutine DoneMeas_OMI_TRC( status )

    ! --- in/out --------------------------------

    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/DoneMeas_OMI_TRC'

    ! --- begin ------------------------------

    ! remove storage for OMI data
    call Dealloc_OMI_TRC( omi_trc, status )
    IF_NOTOK_RETURN(status=1)

    ! analysis settings ?
    if ( omi_trc%analyse ) then
      ! clear storage for species correlation weights:
      deallocate( analyse_spec, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! clear storage for noise correlation weights:
      deallocate( analyse_noise, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    !! clear:
    !deallocate ( nft%trc_col_freetrop )

    ! ok
    status = 0

  end subroutine DoneMeas_OMI_TRC


  ! ***


  subroutine Alloc_OMI_TRC( omitrc, nrmeas, nrlev, status )

    use dims, only : nz

    ! --- in/out -----------------------------

    type(OMI_TRC_type), intent(inout) ::  omitrc
    integer, intent(in)               ::  nrmeas, nrlev
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/Alloc_OMI_TRC'

    ! --- local -----------------------------

    if ( omitrc%nmeas == 0 ) then
       omitrc%nmeas = nrmeas
       omitrc%nlev  = nrlev
       omitrc%nmeas_intimestep = 0
       allocate( omitrc%time(nrmeas)                         )
       allocate( omitrc%date(6,nrmeas)                       )
       allocate( omitrc%longitude(nrmeas)                    )
       allocate( omitrc%latitude (nrmeas)                    )
       allocate( omitrc%corner_longitudes(4,nrmeas)          )
       allocate( omitrc%corner_latitudes (4,nrmeas)          )
       allocate( omitrc%vcd_trop(nrmeas)                     )
       allocate( omitrc%sigma_vcd_trop(nrmeas)               )
       allocate( omitrc%kernel(nrlev,nrmeas)                 )
       allocate( omitrc%pressure_levels(0:nrlev,nrmeas)      )
       allocate( omitrc%cloud_radiance_fraction(nrmeas)      )
       allocate( omitrc%pixel_number(nrmeas)                 )
       allocate( omitrc%pixel_id(nrmeas)                     )
       allocate( omitrc%ix(nrmeas)                           )
       allocate( omitrc%iy(nrmeas)                           )
       allocate( omitrc%status(nrmeas)                       )
       allocate( omitrc%analysed(nrmeas)                     )
       allocate( omitrc%screened(nrmeas)                     )
       allocate( omitrc%intimestep(nrmeas)                   )
       allocate( omitrc%le_altitude_levels(0:nz+nhz,nrmeas)    )
       allocate( omitrc%le_pressure_levels(0:nz+nhz,nrmeas)    )
       allocate( omitrc%le_trc_vmr_prof(1:nz+nhz,nrmeas)       )
       allocate( omitrc%le_trc_vcd_prof(1:nz+nhz,nrmeas)       )
       allocate( omitrc%le2omi_trc_vcd_prof(nrlev,nrmeas)    )
       allocate( omitrc%le2omi_trc_vcd_trop(nrmeas,nfs,nfm)  )
    else
       write (gol,'("KF:     OMI TRC - ERROR arrays already allocated ?!?")') ; call goErr
       TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine Alloc_OMI_TRC


  ! ***


  subroutine Dealloc_OMI_TRC( omitrc, status )

    ! --- in/out -----------------------------

    type(OMI_TRC_type), intent(inout) ::  omitrc
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/ClearMeas_OMI_TRC'

    ! --- local -----------------------------

    if ( omitrc%nmeas > 0 ) then
       deallocate( omitrc%time                )
       deallocate( omitrc%date                )
       deallocate( omitrc%longitude           )
       deallocate( omitrc%latitude            )
       deallocate( omitrc%corner_longitudes   )
       deallocate( omitrc%corner_latitudes    )
       deallocate( omitrc%vcd_trop            )
       deallocate( omitrc%sigma_vcd_trop      )
       deallocate( omitrc%kernel              )
       deallocate( omitrc%pressure_levels     )
       deallocate( omitrc%cloud_radiance_fraction      )
       deallocate( omitrc%pixel_number        )
       deallocate( omitrc%pixel_id            )
       deallocate( omitrc%ix                  )
       deallocate( omitrc%iy                  )
       deallocate( omitrc%status              )
       deallocate( omitrc%analysed            )
       deallocate( omitrc%screened            )
       deallocate( omitrc%intimestep          )
       deallocate( omitrc%le_altitude_levels  )
       deallocate( omitrc%le_pressure_levels  )
       deallocate( omitrc%le_trc_vmr_prof     )
       deallocate( omitrc%le_trc_vcd_prof     )
       deallocate( omitrc%le2omi_trc_vcd_prof )
       deallocate( omitrc%le2omi_trc_vcd_trop )
       omitrc%nmeas = 0
       omitrc%nlev  = 0
       omitrc%nmeas_intimestep = 0
    else
       write (gol,'("KF:     OMI TRC - WARNING Dealloc_OMI_TRC: Arrays are not allocated ?!?")') ; call goPr
    end if

    ! ok
    status = 0

  end subroutine Dealloc_OMI_TRC


  ! ***


  subroutine GetMeas_OMI_TRC( tint1, tint2, status )

    use GO,      only : TDate, NewDate, IncrDate, Get, wrtgol
    use GO,      only : operator(<), operator(<=), operator(+), operator(-), operator(/)
    use GO,      only : goReplace, goUpCase, goLoCase

    use dims,    only : runF
    use LE_Grid, only : ugg

    !use LEKF_dims, only : maxomi

    use NetCDF,  only : NF90_NOWRITE
    use NetCDF,  only : NF90_Open, NF90_Close
    use NetCDF,  only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF,  only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var

    ! --- in/out ----------------------------

    type(TDate), intent(in)   ::  tint1, tint2
    integer, intent(out)      ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/Getmeas_OMI_TRC'

    ! --- local -----------------------------

    integer             ::  year, month, day
    type(TDate)         ::  omi_trc_t
    character(len=256)  ::  fname
    logical             ::  exist

    integer             ::  ncid, varid
    integer             ::  dimid
    integer             ::  dimids(2)

    integer             ::  nrmeas, nrlev, imeas, i, iobs, count_timestep
    integer             ::  nrmeas_in_domain
    character(len=16)   ::  trname
    !character(len=14)   ::  sint1, sint2

    ! --- begin -----------------------------

    !write (sint1,'(i4,2i2.2,1x,i2.2,":",i2.2)') tint1%year,tint1%month,tint1%day,tint1%hour,tint1%min
    !write (sint2,'(i4,2i2.2,1x,i2.2,":",i2.2)') tint2%year,tint2%month,tint2%day,tint2%hour,tint2%min
    !write (gol,*) "KF:     Getmeas_OMI_TRC - Get data for time from ", sint1, " to ", sint2 ;     call goPr
    call wrtgol( rname//' get OMI data for timerange ', (/tint1,tint2/) ); call goPr

    ! no data for this time yet ?
    if ( (tint2 < omi_trc_tr(1)) .or. (omi_trc_tr(2) < tint2) ) then

       ! already something loaded ?
       if ( omi_trc%nmeas > 0 ) then
          ! clear OMI storage:
          call Dealloc_OMI_TRC( omi_trc, status )
          IF_NOTOK_RETURN(status=1)
       end if

       ! extract time values:
       call Get( tint2, year=year, month=month, day=day )

       ! new time range:
       omi_trc_tr(1) = NewDate( year=year, month=month, day=day )
       omi_trc_tr(2) = omi_trc_tr(1) + IncrDate(day=1)

       ! info ...
       write (gol,'(a," - open file for ",i4,2i2.2)') rname, year, month, day; call goPr

       ! start filename with template:
       fname = trim(arch_filenames)
       ! replace some values if necessary:
       call goReplace( fname, '%{year}', '(i4.4)', year, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{yyyy}', '(i4.4)', year, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{mm}', '(i2.2)', month, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{yyyymm}', '(i6.6)', year*100+month, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{yyyymmdd}', '(i8.8)', year*10000+month*100+day, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{TRACER}', goUpCase(trim(omi_trc%tracer_name)), status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{tracer}', goLoCase(trim(omi_trc%tracer_name)), status )
       IF_ERROR_RETURN(status=1)

       ! info ...
       write (gol,'(a," - read from file ",a)') rname, fname; call goPr

       ! available ?
       inquire( file=trim(fname), exist=exist )
       if ( .not. exist ) then
         write (gol,'("WARNING - file not found : ",a)') trim(fname); call goPr
         return
         !write (gol,'("KF:     OMI TRC - file not found : ",a)') trim(fname); call goErr
         !TRACEBACK; status=1; return
       end if

       ! open file:
       status = NF90_Open( fname, NF90_NOWRITE, ncid )
       IF_NF90_NOTOK_RETURN(status=1)

       ! pixel dimension:
       status = NF90_Inq_DimID( ncid, 'pixel', dimid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! number of pixels:
       status = NF90_Inquire_Dimension( ncid, dimid, len=nrmeas )
       IF_NF90_NOTOK_RETURN(status=1)

       ! info ...
       write (gol,'(a," - nr of OMI obs       : ",i6)') rname, nrmeas; call goPr

       ! any pixels?
       if ( nrmeas > 0 ) then

         ! get variable id:
         status = NF90_Inq_Varid( ncid, 'pressure_levels', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         ! get dimension id's:
         status = NF90_Inquire_Variable( ncid, varid, dimids=dimids )
         IF_NF90_NOTOK_RETURN(status=1)
         ! number of (half)levels:
         status = NF90_Inquire_Dimension( ncid, dimids(1), len=nrlev )
         IF_NF90_NOTOK_RETURN(status=1)
         ! convert to full levels:
         nrlev = nrlev - 1

         ! info ...
         write (gol,'(a," - nr of kernel levels : ",i6)') rname, nrlev; call goPr

         ! create storage
         call Alloc_OMI_TRC( omi_trc, nrmeas, nrlev, status )
         IF_NOTOK_RETURN(status=1)

         ! default field values: status etc
         omi_trc%status(:)   = 0
         omi_trc%analysed(:) = .false.
         omi_trc%screened(:) = .false.

         ! fill:
         status = NF90_Inq_Varid( ncid, 'date_time', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%date)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'longitude', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%longitude)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'latitude', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%latitude)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'corner_longitudes', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%corner_longitudes )
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'corner_latitudes', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%corner_latitudes )
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'vcd_trop', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%vcd_trop)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'sigma_vcd_trop', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%sigma_vcd_trop)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'kernel', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%kernel)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'pressure_levels', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%pressure_levels)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'cloud_radiance_fraction', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%cloud_radiance_fraction)
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill:
         status = NF90_Inq_Varid( ncid, 'pixel_number', varid )
         IF_NF90_NOTOK_RETURN(status=1)
         status = NF90_Get_Var( ncid, varid, omi_trc%pixel_number )
         IF_NF90_NOTOK_RETURN(status=1)

         ! fill record numbers:
         do i = 1, nrmeas
           omi_trc%pixel_id(i) = i
         end do

       end if  ! any pixels

       ! close:
       status = NF90_Close( ncid )
       IF_NF90_NOTOK_RETURN(status=1)

       ! info ...
       write (gol,'(a," - file closed ")') rname; call goPr

       ! info ...
       if ( (trim(omi_trc%r_type) == 'data') .and. (abs(omi_trc%r_data_scaling-1.0) > 1.0e-5) ) then
         write (gol,'(a," - obs error scaling   : ",f6.3)') rname, omi_trc%r_data_scaling; call goPr
       end if

       ! init counters:
       omi_trc%nmeas = 0
       nrmeas_in_domain = 0

       ! interpolation from LE grid cells:
       do imeas = 1, nrmeas

         ! * check time:
         omi_trc_t = NewDate( time6=omi_trc%date(:,imeas) )
         if ( ( omi_trc_t < omi_trc_tr(1) ) .or. ( omi_trc_tr(2) < omi_trc_t ) ) then
           write (gol,'("observation time outside expected bounds:")'); call goErr
           write (gol,'("  file name    : ",a)') trim(fname); call goErr
           write (gol,'("  observation  : ",i6)') imeas; call goErr
           call wrtgol( '  time         : ', omi_trc_t ); call goErr
           call wrtgol( ', expected     : ', omi_trc_tr(1), ' - ', omi_trc_tr(2) ); call goErr
           !TRACEBACK; status=1; return
           write (gol,'("  WARNING Just skip and continue  please check   : ",a)') trim(fname); call goErr
           cycle
         end if

         ! * check if meas is in domain;
         ! return cell indices (ix,iy) or status=-1 if not in domain
         call ugg%GetLocation( omi_trc%longitude(imeas), omi_trc%latitude(imeas), &
                                omi_trc%ix(imeas), omi_trc%iy(imeas), status, &
                                quiet=.true. )
         IF_ERROR_RETURN(status=1)
         ! not in domain ? skip immediately:
         if ( status < 0 ) cycle

         ! check validity of pixel number:
         if ( omi_trc%pixel_number(imeas) <= 0 .or. omi_trc%pixel_number(imeas) > 60 ) then
           write( gol, '("Pixel number outside possible bounds [1-60], number is: ",i0)' ) omi_trc%pixel_number(imeas) ; call goErr
           TRACEBACK;status=1;return
         end if

         ! limitted pixel range?
         if ( omi_trc%pixel_range(1) > 0 ) then
           if ( omi_trc%pixel_number(imeas) < omi_trc%pixel_range(1) ) cycle
         end if
         if ( omi_trc%pixel_range(2) > 0 ) then
           if ( omi_trc%pixel_number(imeas) > omi_trc%pixel_range(2) ) cycle
         end if

         !! TESTING: pixel numbers [26,54] are not valid in 2014; cancel out of assimilation
         !if ( omi_trc%pixel_number(imeas) >= 26 .and. omi_trc%pixel_number(imeas) <= 54 ) cycle

         ! increase counter:
         nrmeas_in_domain = nrmeas_in_domain + 1

         ! switch:
         select case ( trim(omi_trc%r_type) )
           !~ from data?
           case ( 'data' )
             ! Apply a scaling to the OMI observation error
             omi_trc%sigma_vcd_trop(imeas) = omi_trc%sigma_vcd_trop(imeas) * omi_trc%r_data_scaling
           !~ fraction of data:
           case ( 'frac' )
             ! fraction of observed value, bounded:
             omi_trc%sigma_vcd_trop(imeas) = min( max( omi_trc%r_frac_min, omi_trc%vcd_trop(imeas) * omi_trc%r_frac_factor ), omi_trc%r_frac_max )
           !~
           case default
             write (gol,'("unsupported OMI TRC r type `",a,"`")') trim(omi_trc%r_type); call goErr
             TRACEBACK; status=1; return
         end select

         ! increase counter:
         omi_trc%nmeas = omi_trc%nmeas + 1
         ! restore:
         omi_trc%pixel_id(omi_trc%nmeas)          = omi_trc%pixel_id(imeas)
         omi_trc%date(:,omi_trc%nmeas)            = omi_trc%date(:,imeas)
         omi_trc%time(omi_trc%nmeas)              = omi_trc%time(imeas)
         omi_trc%longitude(omi_trc%nmeas)         = omi_trc%longitude(imeas)
         omi_trc%latitude(omi_trc%nmeas)          = omi_trc%latitude(imeas)
         omi_trc%vcd_trop(omi_trc%nmeas)          = omi_trc%vcd_trop(imeas)
         omi_trc%sigma_vcd_trop(omi_trc%nmeas)    = omi_trc%sigma_vcd_trop(imeas)
         omi_trc%kernel(:,omi_trc%nmeas)          = omi_trc%kernel(:,imeas)
         omi_trc%pressure_levels(:,omi_trc%nmeas) = omi_trc%pressure_levels(:,imeas)
         omi_trc%cloud_radiance_fraction(omi_trc%nmeas)    = omi_trc%cloud_radiance_fraction(imeas)
         omi_trc%ix(omi_trc%nmeas)                = omi_trc%ix(imeas)
         omi_trc%iy(omi_trc%nmeas)                = omi_trc%iy(imeas)

         ! * check value:
         ! HE: Do not exclude negative values!
         ! HE: In very clean areas there are negative values that should be included
         ! HE: since otherwise a bias towards positive values is introduced ....
         ! if ( omi_trc%vcd_trop(imeas) < 0.0 ) omi_trc%status(imeas) = omi_trc%status(imeas) + status_nodata

       end do  ! pixels

       write (gol,'(a," - nr of observations  = ",i0)') rname, omi_trc%nmeas; call goPr
       write (gol,'(a," - nr of obs in domain = ",i0)') rname, nrmeas_in_domain; call goPr
       !if ( omi_trc%nmeas > 0 ) then
       !  write (gol,*) "KF:   OMI TRC - last OMI observation read : ";     call goPr
       !  write (gol,*) "KF:     time ",omi_trc%time(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     lon  ",omi_trc%longitude(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     lat  ",omi_trc%latitude(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     vcd  ",omi_trc%vcd_trop(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     sig  ",omi_trc%sigma_vcd_trop(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     kern ",omi_trc%kernel(:,omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     pres ",omi_trc%pressure_levels(:,omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     clfr ",omi_trc%cloud_radiance_fraction(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     ix   ",omi_trc%ix(omi_trc%nmeas);     call goPr
       !  write (gol,*) "KF:     iy   ",omi_trc%iy(omi_trc%nmeas);     call goPr
       !end if

    else

       write (gol,'(a," - data already available for this day")') rname; call goPr

    end if

    ! data from file has been read for the current day and is stored at this point
    if ( omi_trc%nmeas == 0 ) then
      ! info ...
      write (gol,'(a," - no observations for today ... ")') rname; call goPr
      ! no data at all, so nothing left either ...
      omi_trc%nodataleft = .true.

    else

       write (gol,*) "KF:     OMI TRC - first time ",omi_trc%date(:,1);     call goPr
       write (gol,*) "KF:     OMI TRC - last  time ",omi_trc%date(:,omi_trc%nmeas);     call goPr

      ! -----------------------------------------------
      ! now fill "omi_trc%intimestep", which fags only those measurement
      ! of relevance for this Kalman time step
      ! -----------------------------------------------

      omi_trc%intimestep(:) = .false.
      omi_trc%nmeas_intimestep = 0
      omi_trc%nodataleft = .true.

      do i = 1, omi_trc%nmeas

        ! get OMI time
        omi_trc_t = NewDate( time6=omi_trc%date(:,i) )
        ! check if there is still valid data after this time interval, with times > tint2:
        if ( (tint2 < omi_trc_t) .and. (omi_trc%status(i) == status_default) ) omi_trc%nodataleft = .false.
        ! skip if outside (tint1,tint2] :
        if ( (omi_trc_t <= tint1) .or. (tint2 < omi_trc_t) ) then
               cycle
        end if

        ! skip if not accepted (validation, no data?):
        if ( omi_trc%status(i) /= status_default ) then
          cycle
        end if

        ! flag to signal that this observation needs to be processed in this assimilation time step
        omi_trc%intimestep(i) = .true.

        ! count number of relevant obs
        omi_trc%nmeas_intimestep = omi_trc%nmeas_intimestep + 1

      end do

    end if  ! obs avail ?

    write (gol,'("KF:   OMI TRC - nr of valid obs in timestep : ",i6)') omi_trc%nmeas_intimestep ; call goPr

    ! ok
    status = 0

  end subroutine GetMeas_OMI_TRC


  ! ***


  !
  ! Simulate omi trc given the model concentrations
  !
  ! On return the allocatable array 'le2omi_vcd_trop' is allocated to 'npix'
  ! which the number of pixels; if necessary, the array is dealocated first.
  !

  subroutine CalcMeas_OMI_TRC( le_vmr_grid, le_vmr_aloft, &
                                 npix, le2omi_vcd_trop, status )

    use Binas       , only : Avog, xm_air, grav
    use JAQL        , only : PotentialPressure
    use Num         , only : IntervalSum

    use dims        , only : nx, ny, nz, nspec
    use LE_Bound    , only : caloft
    use LE_Data     , only : LE_Data_GetPointer
    use indices     , only : specname, specmolm

    use LE_Bound_Top, only : LE_Bound_Top_MergeProfiles
    use LE_Data     , only : LE_Data_GetPointer

    ! --- in/out -----------------------------

    real, intent(in)                  ::  le_vmr_grid(nx,ny,nz)   ! volume ppb
    real, intent(in)                  ::  le_vmr_aloft(nx,ny)     ! volume ppb
    integer, intent(out)              ::  npix
    real, allocatable, intent(out)    ::  le2omi_vcd_trop(:)      ! unit 1e15 molec/cm2
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/CalcMeas_OMI_TRC'

    ! --- local --------------------------------

    real, pointer                   ::  temp(:,:,:)
    real, pointer                   ::  hm(:,:,:)
    real, pointer                   ::  oro(:,:,:)

    integer, parameter :: nlevmax = 100
    integer     ::  imeas
    integer     ::  ix, iy, iz
    integer     ::  ilev
    real        ::  omi_ph(0:nlevmax)
    real        ::  le2omi_vcd(1:nlevmax)
    real, allocatable  ::  le_hh(:)
    real, allocatable  ::  le_ph(:)
    real, allocatable  ::  le_temp(:)
    real, allocatable  ::  le_rh(:)
    real, allocatable  ::  le_airm(:)
    real, allocatable  ::  le_vmr(:)
    real, allocatable  ::  le_vcd(:)
    !real, allocatable  ::  x(:)   ! a shifted copy of le_ph

    integer     ::  ilast, nlevx, nmeasx, nmeas_processed
    real        ::  le_ph_dum

    real                ::  du_to_kgm2
    real, pointer       ::  mer_phlev(:)
    real, pointer       ::  mer_vcd(:)


    ! --- local --------------------------------

    ! meteo:
    call LE_Data_GetPointer( 't', temp, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', hm, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'oro', oro, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)

    ! allocate
    allocate( le_hh(0:nz+nhz) )
    allocate( le_ph(0:nz+nhz) )
    allocate( le_temp(1:nz+nhz) )
    allocate( le_rh(1:nz+nhz) )
    allocate( le_airm(1:nz+nhz) )
    allocate( le_vmr(1:nz+nhz) )
    allocate( le_vcd(1:nz+nhz) )
    !allocate( x(1:nz+3) )         ! a shifted copy of le_ph

    ! clear output:
    if ( allocated(le2omi_vcd_trop) ) then
      deallocate( le2omi_vcd_trop, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! set length of output:
    npix = omi_trc%nmeas
    ! any values ?
    if ( npix > 0 ) then
      ! storage:
      allocate( le2omi_vcd_trop(npix) )
      ! init values for safety:
      le2omi_vcd_trop = 0.0
    end if

    ! no data ?
    if ( omi_trc%nmeas < 1 .or. omi_trc%nmeas_intimestep < 1 ) then
       write (gol,'(a," - skip, no observations")') rname; call goPr
       status = 0; return
    end if

    write (gol,'(a," - compute LE predicted OMI TRC, ",i6," obs ...")') rname, omi_trc%nmeas_intimestep; call goPr

    nmeasx = omi_trc%nmeas
    nlevx = omi_trc%nlev
    nmeas_processed = 0

    ! init top model profiles:
    nullify( mer_phlev )
    nullify( mer_vcd )

    ! loop over omi measurements:
    do imeas = 1, nmeasx

       ! first check if this obs needs processing
       if ( .not. omi_trc%intimestep(imeas) ) cycle

       nmeas_processed = nmeas_processed + 1

       ! omi half level pressures:
       omi_ph(0:nlevx) = omi_trc%pressure_levels(:,imeas)

       ! model grid cell:
       ix = omi_trc%ix(imeas)
       iy = omi_trc%iy(imeas)

       ! model half level heights:
       le_hh(0)    = oro(ix,iy,1)  ! m
       le_hh(1:nz) = le_hh(0) + hm(ix,iy,1:nz) ! m
       !! the "aloft" layer with free troposphere TRC is distributed between 3.5 and 7 km ..
       !le_hh(nz+1) = topFreeTropLayer
       !! the aloft-aloft layer has TRC=0 and ranges between 7 and 200 km
       !le_hh(nz+2) = topAtmosphere

       ! model temperature:
       le_temp(1:nz) = temp(ix,iy,1:nz)
       !le_temp(nz+1) = le_temp(nz)
       !le_temp(nz+2) = le_temp(nz)

       ! model relative humidity;
       ! assume zero since LE rh is probably not what you expect ...
       le_rh = 0.0

       ! model pressure profile:
       le_ph(0) = omi_ph(0)   ! use same surface pressure as omi
       do iz = 1, nz+nhz
          ! pressure decrease given layer heigt and temperature,
          ! assume zero relative humidity since LE rh is probably not
          ! what you expect ...
          call PotentialPressure( le_ph_dum, le_ph(iz-1), le_temp(iz), le_rh(iz), le_hh(iz)-le_hh(iz-1) )
          le_ph(iz) = le_ph_dum
       end do
       !! set pressure top at zero, otherwise problems with OMI kernels that reach to this level too:
       !le_ph(nz+2) = 0.0

       !! check ...
       !if ( le_ph(nz+2) > omi_trc%pressure_levels(nlevx,imeas) ) then
       !   write (gol,'("top of LE layers should be at or over OMI TRC layers:")'); call goErr
       !   write (gol,'("  OMI pressure levels : ",100f12.2)') omi_trc%pressure_levels(:,imeas); call goErr
       !   write (gol,'("  LE pressure levels : ",100f12.2)') le_ph; call goErr
       !   write (gol,*) '  LE top  : ', le_ph(nz+2); call goErr
       !   write (gol,*) '  OMI top : ', omi_trc%pressure_levels(nlevx,imeas); call goErr
       !   write (gol,'("    (pixel ",i6,", level ",i6,")")') imeas, nlevx; call goErr
       !   TRACEBACK; status=1; return
       !end if

       ! LE air mass:
       do iz = 1, nz+nhz
          le_airm(iz) = ( le_ph(iz-1) - le_ph(iz) )/grav  ! (kg air)/m2
       end do

       ! le trc vmr:
       le_vmr(1:nz) = le_vmr_grid (ix,iy,:)      ! trc volume ppb
       !! le_vmr(nz+1) = le_vmr_aloft(ix,iy)        ! trc volume ppb
       !! Henk Eskes: set initial aloft concentrations equal to 0
       !le_vmr(nz+1) = 0.0                        ! trc volume ppb
       !le_vmr(nz+2) = 0.0                        ! trc volume ppb

       ! Compute partial TRC columns in (1e15 mlc trc)/cm2 :
       !    c   : LE volume mixing ratios (ppb)
       !
       !  1e15 mlc trc    mol trc    mol air  kg air  m2
       !  ---- -------  -----------  -------  ------  ---  / 1e15  = (1e15 mlc trc)/cm2
       !   1   mol trc  1e9 mol air  kg air     m2    cm2
       !
       do iz = 1, nz
          le_vcd(iz) = 1.0e-15 * Avog * le_vmr(iz) / 1.0e9 / xm_air * le_airm(iz) * 1.0e-4
          ! unit (1e15 mlc trc)/cm2
       end do
       !le_vcd(nz+1) = 0.0
       !le_vcd(nz+2) = 0.0
       !! Henk Eskes: the LE aloft is replaced by the free tropospheric column read from file
       !!             unit: 1e15 mlc trc / cm2
       !if ( omi_trc_addfreetrop ) then
       !   le_vcd(nz+1) = nft%trc_col_freetrop(ix,iy)
       !   le_vmr(nz+1) = le_vcd(nz+1) / 1.0e-15 / Avog * 1.0e9 * xm_air / le_airm(nz+1) / 1.0e-4
       !end if

       !... combined LE and parent model profile ...
       ! convert from "(1e15 mlc)/cm2" to "kg/m2" :
       ! (1e15 mlc)/cm2 * 1e15 / (mlc/mol) *       (kg/mol)            * cm2/m2
       du_to_kgm2       = 1e15 /   Avog    * specmolm(omi_trc_itracer) * 1.0e4
       ! get parent model profile; return status:
       !    0 : ok, profile extracted
       !   -1 : no top field stored for requested tracer
       !   -2 : top field defined, but not filled
       call LE_Bound_Top_MergeProfiles( ix, iy, omi_trc_itracer, &
                                         omi_trc%top_ph_entry, omi_trc%top_trc_entry, &
                                         le_ph(0:nz), le_vcd(1:nz)*du_to_kgm2, 'kg/m2', &
                                         mer_phlev, mer_vcd, status )
       IF_NOTOK_RETURN(status=1)
       ! covnert to units of original model profile:
       mer_vcd = mer_vcd / du_to_kgm2

       !! testing ...
       !print *, 'xxx1 ', ix, iy, omi_trc_itracer
       !print *, '  x1 phlev = ', mer_phlev
       !print *, '  x1 vcd   = ', mer_vcd
       !stop

       !...

       ! LE trc profile at OMI layers:
       ilast = 1
       !x(1:nz+3) = -1.0*le_ph(0:nz+2)  ! an increasing array starting with index=1
       ! loop over omi pressure layers:
       do ilev = 1, nlevx

          ! only layers completely below ptropo:
          if ( min( omi_ph(ilev-1), omi_ph(ilev) ) >= omi_trc_ptropo ) then

            ! fill each omi layer with sum of one or more fractions of LE layers;
            ! use pressure axis to have mass-conservation; negate to have increasing axis:
            ! old: call IntervalSum( -1.0*le_ph, l
            !call IntervalSum( x, le_vcd, &
            ! use merged axes:
            call IntervalSum( -1.0*mer_phlev, mer_vcd, &
                               -1.0*omi_ph(ilev-1), -1.0*omi_ph(ilev), &
                               le2omi_vcd(ilev), ilast, status )
            if ( status /= 0 ) then
              write (gol,*) 'from IntervalSum:'; call goErr
              write (gol,*) '  imeas   = ', imeas, ' / ', nmeasx; call goErr
              write (gol,*) '  loc     = ', omi_trc%longitude(imeas), omi_trc%latitude(imeas); call goErr
              write (gol,*) '  cell    = ', omi_trc%ix(imeas), omi_trc%iy(imeas); call goErr
              write (gol,*) '  ilev    = ', ilev; call goErr
              write (gol,*) '  le_hh   = ', le_hh; call goErr
              write (gol,*) '  le_ph   = ', le_ph; call goErr
              write (gol,*) '  le_arim = ', le_airm; call goErr
              write (gol,*) '  le_vmr  = ', le_vmr; call goErr
              write (gol,*) '  le_vcd  = ', le_vcd; call goErr
              write (gol,*) '  nlevx   = ', nlevx; call goErr
              write (gol,*) '  omi_ph  = ', omi_ph(0:nlevx); call goErr
              !write (gol,*) '  x       = ', x; call goErr
              write (gol,*) '  dp      = ', -1.0*omi_ph(ilev-1), -1.0*omi_ph(ilev); call goErr
              TRACEBACK; status=1; return
            end if

          else

            ! not in troposphere ...
            le2omi_vcd(ilev) = 0.0

          end if

       end do  ! omi profile layers

       ! apply kernel (Eskian apriori is zero ..)
       le2omi_vcd_trop(imeas) = 0.0
       do ilev = 1, nlevx
          le2omi_vcd_trop(imeas) = le2omi_vcd_trop(imeas) + &
               omi_trc%kernel(ilev,imeas) * le2omi_vcd(ilev)
       end do

       ! store LE diagnostic values
       omi_trc%le_altitude_levels  (:,imeas)  = le_hh(:)
       omi_trc%le_pressure_levels  (:,imeas)  = le_ph(:)
       omi_trc%le_trc_vmr_prof     (:,imeas)  = le_vmr(:)
       omi_trc%le_trc_vcd_prof     (:,imeas)  = le_vcd(:)
       omi_trc%le2omi_trc_vcd_prof (:,imeas)  = 0.0
       omi_trc%le2omi_trc_vcd_prof (1:nlevx,imeas)  = le2omi_vcd(1:nlevx)

       ! Eskes print statement
       ! write(gol,'("CalcMeas_OMI_TRC: imeas, date, le2omi_vcd_trop = ",i5,f8.4)') imeas, omi_trc%date(:,imeas),le2omi_vcd_trop(imeas) ; call goPr
       ! write(gol,'("CalcMeas_OMI_TRC: imeas, date, le2omi_vcd_trop = ",i5,f8.4,f8.4)') imeas, omi_trc%date(:,imeas),le2omi_vcd_trop(imeas) ; call goPr
       ! write (gol,*) "KF:   OMI TRC - Nr of observations  = ", omi_trc%nmeas;     call goPr
       ! write (gol,*) "CalcMeas_OMI_TRC: imeas, date, le2omi_vcd_trop =  ", imeas, omi_trc%date(:,imeas),le2omi_vcd_trop(imeas) ;   call goPr

       !! testing ..
       !stop

    end do  ! omi meas

    if ( omi_trc%nmeas_intimestep /= nmeas_processed ) then
       write (gol,'(a," - WARNING CalcMeas: In total ",i6," obs processed")') rname, nmeas_processed ; call goPr
    end if

    ! deallocate
    deallocate( le_hh )
    deallocate( le_ph )
    deallocate( le_temp )
    deallocate( le_rh )
    deallocate( le_airm )
    deallocate( le_vmr )
    deallocate( le_vcd )
    !deallocate( x )

    ! clear top model profiles if ncessary:
    if ( associated(mer_phlev) ) deallocate( mer_phlev )
    if ( associated(mer_vcd  ) ) deallocate( mer_vcd   )

    ! ok
    status = 0

  end subroutine CalcMeas_OMI_TRC


  ! ======================================================================
  ! ===
  ! === output
  ! ===
  ! ======================================================================

  subroutine meas_output_omi_Init( self, rcF, rckey, status )

    use GO              , only : AnyDate
    use GO              , only : TrcFile
    use LE_Output_Common, only : Init

    ! --- in/out --------------------------------

    class(Meas_Output_OMI_TRC), intent(out) ::  self
    type(TRcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/meas_output_omi_Init'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! init common stuff:
    call Init( self%com, rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files?
    call rcF%Get( trim(rckey)//'.replace', self%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! files not open yet:
    self%opened = .false.

    ! no time range set yet:
    self%tr(1) = AnyDate()
    self%tr(2) = AnyDate()

    ! ok
    status = 0

  end subroutine meas_output_omi_Init


  ! ***


  subroutine meas_output_omi_Done( self, status )

    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    class(Meas_Output_OMI_TRC), intent(inout)     ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/meas_output_omi_Done'

    ! --- begin ---------------------------------

    ! file opened ?
    if ( self%opened ) then
       ! close:
       status = NF90_Close( self%ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! reset flag:
       self%opened = .true.
    end if

    ! done with common stuff ...
    call Done( self%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine meas_output_omi_Done


  ! ***

  ! Output is generate only once per day, after all OMI data has been processed.
  ! The subroutine will accumulate the relevant information to print until
  ! all data has been analysed

  subroutine meas_output_omi_PutOut( self, key, t, status, last )

    use GO, only : TDate, IncrDate, Get, NewDate, AnyDate
    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal
    use GO, only : Precisely, MidNight, wrtgol
    use GO, only : goc

#ifdef _MPI
    use NetCDF , only : NF90_Create_Par
    use NetCDF , only : NF90_NETCDF4
#else
    use NetCDF , only : NF90_Create
#endif
    use NetCDF , only : NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE

    use Dims            , only : nz
    use LE_Bound        , only : caloft
    use LE_Output_Common, only : PutOut_GlobalAttributes

    use LEKF_State, only : kf_with_xb, kf_with_xm
    use LEKF_State, only : xb, x, sigma

    ! --- in/out --------------------------------

    class(Meas_Output_OMI_TRC), intent(inout) ::  self
    character(len=*), intent(in)              ::  key
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    logical, intent(in), optional  ::  last

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/kfmo_meas_PutOut'

    ! --- local ---------------------------------

    logical               ::  islast

    logical               ::  newfile
    logical               ::  writetofile
    type(TDate)           ::  tday
    type(TDate)           ::  omi_trc_t
    integer               ::  time6(6)
    real                  ::  time

    integer               ::  cmode
    integer               ::  varid
    type(TDate)           ::  t0

    integer               ::  ifs
    character(len=32)     ::  varname
    character(len=512)    ::  descr
    integer               ::  imeas
    integer               ::  ifm
    integer               ::  ient
    integer               ::  nmeas, nlev, i
    integer               ::  nmeas_tot, istart
    integer               ::  nmeas_intimestep_tot

    integer               ::  n_omi_trc
    real, allocatable     ::  le2omi_trc_vcd_trop(:)
    integer, allocatable  ::  dummy(:)

    ! --- begin ---------------------------------

    ! flag:
    islast = .false.
    if ( present(last) ) islast = last

    ! dims:
    nmeas = omi_trc%nmeas
    nlev  = omi_trc%nlev

    ! total number over all processors, and start index in output:
    call goc%ParInfo( nmeas, nmeas_tot, istart, status )
    IF_NOTOK_RETURN(status=1)

    ! no measurements ? then leave
    if ( nmeas_tot < 1 ) then
       ! info ...
       write (gol,'(a," - no trc obs for today")') rname; call goPr
       return
    end if

    ! necessary?
    nmeas_intimestep_tot = omi_trc%nmeas_intimestep
    call goc%AllReduce( 'sum',  nmeas_intimestep_tot, status )
    IF_NOTOK_RETURN(status=1)
    ! no measurements within timestep ? then leave
    if ( nmeas_intimestep_tot < 1 ) then
       ! info ...
       call wrtgol( rname//' - no trc obs for time : ',t); call goPr
       return
    end if

    ! which filter moment ?
    select case ( key )
    case ( 'forecast' ) ; ifm = 1
    case ( 'analysis' ) ; ifm = 2
    case default
       write (gol,'("unsupported key : ",a)') trim(key); call goErr
       TRACEBACK; status=1; return
    end select

    ! Compute LE predictions of the OMI measurements
    ! for all the filter states and current filter moment
    do ifs = 1, nfs

       ! simulate omi stuff:
       select case ( trim(fsname(ifs)) )
       case ( 'xb' )
          ! skip ?
          if ( .not. kf_with_xb ) cycle
          ! simulate omi stuff from le fields, allocate result if necessary:
          call CalcMeas_OMI_TRC( xb%c(:,:,:,omi_trc_itracer), caloft(:,:,omi_trc_itracer), &
                                   n_omi_trc, le2omi_trc_vcd_trop, status )
          IF_NOTOK_RETURN(status=1)
       case ( 'x' )
          ! skip ?
          if ( .not. kf_with_xm ) cycle
          ! simulate omi stuff from le fields, allocate result if necessary:
          call CalcMeas_OMI_TRC( x%c(:,:,:,omi_trc_itracer), caloft(:,:,omi_trc_itracer), &
                                   n_omi_trc, le2omi_trc_vcd_trop, status )
          IF_NOTOK_RETURN(status=1)
       case ( 's' )
          ! skip ?
          if ( .not. kf_with_xm ) cycle
          ! simulate omi stuff from le fields, allocate result if necessary:
          call CalcMeas_OMI_TRC( sigma%c(:,:,:,omi_trc_itracer), caloft(:,:,omi_trc_itracer), &
                                   n_omi_trc, le2omi_trc_vcd_trop, status )
          IF_NOTOK_RETURN(status=1)
       case default
          write (gol,'("unsupported filter state `",a,"`")') trim(fsname(ifs)); call goErr
          TRACEBACK; status=1; return
       end select

       ! Store the fields for all filter states
       do i = 1, nmeas
          ! First check if this obs is relevant for this time step
          if ( omi_trc%intimestep(i) ) then
             omi_trc%le2omi_trc_vcd_trop(i,ifs,ifm) = le2omi_trc_vcd_trop(i)
          end if
       end do

    end do  ! ifs = 1, nfs

    ! storag allocated?
    if ( allocated(le2omi_trc_vcd_trop) ) then
      ! clear:
      deallocate( le2omi_trc_vcd_trop, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! check if all valid data has been analysed
    writetofile = ( omi_trc%nodataleft .and. (key == 'analysis') ) .or. islast
    ! info ..
    write (gol,'(a," - write flag ",l1,": no data left (",l1,"), key analysis (",a,"), is last (",l1,")")') &
                          rname, writetofile, omi_trc%nodataleft, trim(key), islast; call goPr

    ! all processes should agree on writing ...
    call goc%AllReduce( 'and', writetofile, status )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( .not. writetofile ) then
      write (gol,'(a," - no output yet; not all processes agree to write ...")') rname; call goPr
      return
    end if

    ! ------------------------------------------------------------------
    ! Now write to file because all status=0 OMI data has been analysed
    ! ------------------------------------------------------------------

    ! close old file if necessary (should not be necessary)
    if (  self%opened ) then
       ! close:
       status = NF90_Close( self%ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! reset flag:
       self%opened = .false.
    end if

    ! info ...
    call wrtgol( rname//' - put out all data for this day at : ',t); call goPr

    ! day is defined for (00,24]
    tday = t
    if ( MidNight(t) ) tday = tday - IncrDate(day=1)

    ! extract time fields:
    call Get( tday, time6=time6 )

    ! time range for this file is (00,24]
    self%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
    self%tr(2) = self%tr(1) + IncrDate( day=1 )

    ! new file name:
    write (self%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
         trim(self%com%outdir), trim(self%com%model), trim(self%com%expid), &
         'meas-omi-trc', time6(1:3)

    ! set creation mode flag:
    if ( self%replace ) then
      cmode = NF90_CLOBBER       ! overwrite existing files
    else
      cmode = NF90_NOCLOBBER     ! do not overwrite existing files
    end if

    ! create file:
#ifdef _MPI
    ! create netcdf4 file for parallel output:
    cmode = ior( cmode, NF90_NETCDF4 )
    ! create:
    status = NF90_Create_Par( self%fname, cmode, &
                              goc%comm%mpi_val, goc%info%mpi_val, &
                              self%ncid )
    if (  status /= NF90_NOERR ) then
      gol=nf90_strerror(status); call goErr
      write (gol,'("creating file :")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
#else
    status = NF90_Create( self%fname, cmode, self%ncid )
    if ( status /= 0 ) then
       write (gol,'("creating file :")'); call goErr
       write (gol,'("  ",a)') trim(self%fname); call goErr
       TRACEBACK; status=1; return
    end if
#endif

    ! reset flag:
    self%opened = .true.

    ! write global attributes:
    call PutOut_GlobalAttributes( self%com, self%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! define dimensions:
    status = NF90_Def_Dim( self%ncid, 'omi_obs' , nmeas_tot, self%dimid_omi_obs )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( self%ncid, 'omi_lev' , omi_trc%nlev, self%dimid_omi_lev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( self%ncid, 'omi_hlev', omi_trc%nlev+1, self%dimid_omi_hlev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( self%ncid, 'le_lev'  , nz+nhz, self%dimid_le_lev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( self%ncid, 'le_hlev' , (nz+nhz)+1, self%dimid_le_hlev )
    IF_NF90_NOTOK_RETURN(status=1)
    !       status = NF90_Def_Dim( self%ncid, 'time', NF90_UNLIMITED, self%dimid_time )
    !       IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( self%ncid, 'datelen', 6, self%dimid_datelen )
    IF_NF90_NOTOK_RETURN(status=1)

    ! define variable:
    status = NF90_Def_Var( self%ncid, 'pixel_id', NF90_INT, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'pixel number in input file' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_pixel_id = varid

    ! define daily constant variables:
    status = NF90_Def_Var( self%ncid, 'lon', NF90_REAL, self%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_east' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_lon = varid

    status = NF90_Def_Var( self%ncid, 'lat', NF90_REAL, self%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_lat = varid

    status = NF90_Def_Var( self%ncid, 'omi_time', NF90_REAL, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'days since 2000-01-01 00:00:00' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'calender', 'gregorian' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_time = varid

    status = NF90_Def_Var( self%ncid, 'omi_date', NF90_INT, (/self%dimid_datelen,self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI date and time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'year, month, day, hour, minute, second' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_date = varid

    status = NF90_Def_Var( self%ncid, 'omi_pressure_levels', NF90_REAL, (/self%dimid_omi_hlev,self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI pressure halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_pressure_levels = varid

    status = NF90_Def_Var( self%ncid, 'omi_kernel', NF90_REAL, (/self%dimid_omi_lev,self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI kernel' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_kernel = varid

    status = NF90_Def_Var( self%ncid, 'omi_cloud_radiance_fraction', NF90_REAL, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI cloud fraction' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '[0,1]' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_cloud_radiance_fraction = varid

    status = NF90_Def_Var( self%ncid, 'omi_status', NF90_INT, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI filter status' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(status_description) )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_status = varid

    status = NF90_Def_Var( self%ncid, 'omi_analysed', NF90_INT, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI analysed flag' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'flag: 0 = false, 1 = true' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_analysed = varid

    status = NF90_Def_Var( self%ncid, 'omi_screened', NF90_INT, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI screening flag' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'flag: 0 = false, 1 = true' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_omi_screened = varid

    status = NF90_Def_Var( self%ncid, 'omi_vcd_trop_y', NF90_REAL, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI TRC tropospheric column density' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '(1e15 mlc '//trim(omi_trc%tracer_name)//')/cm2' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_y = varid

    status = NF90_Def_Var( self%ncid, 'omi_vcd_trop_r', NF90_REAL, (/self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'OMI TRC tropospheric column density error std.dev.' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '(1e15 mlc '//trim(omi_trc%tracer_name)//')/cm2' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_r = varid

    ! define temporal variables:

!     status = NF90_Def_Var( self%ncid, 'time', NF90_REAL, (/self%dimid_omi_obs/), varid )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'time' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'long_name', 'time' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'units', 'days since 2000-01-01 00:00:00' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'calender', 'gregorian' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     self%varid_time = varid

!     status = NF90_Def_Var( self%ncid, 'date', NF90_INT, (/self%dimid_datelen,self%dimid_omi_obs/), varid )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'long_name', 'date and time' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     status = NF90_Put_Att( self%ncid, varid, 'units', 'year, month, day, hour, minute, second' )
!     IF_NF90_NOTOK_RETURN(status=1)
!     self%varid_date = varid

    status = NF90_Def_Var( self%ncid, 'le_altitude_levels', NF90_REAL, (/self%dimid_le_hlev,self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS altitude halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'm' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_le_altitude_levels = varid

    status = NF90_Def_Var( self%ncid, 'le_pressure_levels', NF90_REAL, (/self%dimid_le_hlev,self%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS pressure halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOTOK_RETURN(status=1)
    self%varid_le_pressure_levels = varid

    ! loop over filter moments (forecast, analysis)
    do ifm = 1, nfm

       ! loop over filter states:
       do ifs = 1, nfs

          ! filter state description, skip if possible
          select case ( trim(fsname(ifs)) )
          case ( 'xb' )
             ! skip ?
             if ( .not. kf_with_xb ) cycle
             if ( ifm == ifm_forecast ) cycle
             ! description:
             write (descr,'("model background run")')
          case ( 'x' )
             ! skip ?
             if ( .not. kf_with_xm ) cycle
             ! description:
             write (descr,'("ensemble mean; ",a)') trim(fm_descr(ifm))
          case ( 's' )
             ! skip ?
             if ( .not. kf_with_xm ) cycle
             ! description:
             write (descr,'("ensemble std.dev.; ",a)') trim(fm_descr(ifm))
          case default
             write (gol,'("unsupported filter state : ",a)') trim(fsname(ifs)); call goErr
             TRACEBACK; status=1; return
          end select

          ! variable name:
          write (varname,'("conc_",a)') trim(fsname(ifs))

          ! extend for forecast moment:
          if ( ifm == ifm_forecast ) varname = trim(varname)//'_f'

          ! loop over tracer entities to be written:
          do ient = 1, nent

             ! variable name:
             write (varname,'(a,"_",a)') trim(entname(ient)), trim(fsname(ifs))

             ! extend for forecast moment:
             if ( ifm == ifm_forecast ) varname = trim(varname)//'_f'

             ! define variables for each entity:
             select case ( ient )

             case ( ient_le2omi_vcd_trop )

                status = NF90_Def_Var( self%ncid, trim(varname), NF90_REAL, (/self%dimid_omi_obs/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS '//trim(omi_trc%tracer_name)//' tropospheric column density' )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'units', '(1e15 mlc '//trim(omi_trc%tracer_name)//')/cm2' )
                IF_NF90_NOTOK_RETURN(status=1)

             case ( ient_le2omi_vcd )

                status = NF90_Def_Var( self%ncid, trim(varname), NF90_REAL, (/self%dimid_omi_lev,self%dimid_omi_obs/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS '//trim(omi_trc%tracer_name)//' vertical column density at OMI layers' )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'units', '(1e15 mlc '//trim(omi_trc%tracer_name)//')/cm2' )
                IF_NF90_NOTOK_RETURN(status=1)

             case ( ient_le_vcd )

                status = NF90_Def_Var( self%ncid, trim(varname), NF90_REAL, (/self%dimid_le_lev,self%dimid_omi_obs/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS '//trim(omi_trc%tracer_name)//' vertical column density' )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'units', '(1e15 mlc '//trim(omi_trc%tracer_name)//')/cm2' )
                IF_NF90_NOTOK_RETURN(status=1)

             case ( ient_le_vmr )

                status = NF90_Def_Var( self%ncid, trim(varname), NF90_REAL, (/self%dimid_le_lev,self%dimid_omi_obs/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'long_name', 'LOTOS-EUROS '//trim(omi_trc%tracer_name)//' volume-mixing-ratio' )
                IF_NF90_NOTOK_RETURN(status=1)
                status = NF90_Put_Att( self%ncid, varid, 'units', 'ppb' )
                IF_NF90_NOTOK_RETURN(status=1)

             case default

                write (gol,'("unsupported entity : ",i6)') ient; call goErr
                TRACEBACK; status=1; return

             end select

             ! standard attributes:
             status = nf90_put_att( self%ncid, varid, 'description', trim(descr) )
             IF_NF90_NOTOK_RETURN(status=1)

             ! store variable id:
             self%varid(ient,ifs,ifm) = varid

          end do  ! output entities

       end do  ! filter moments

    end do  ! filter states

    ! end defintion mode:

    status = NF90_EndDef( self%ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! no records written yet:
    self%itrec = 0
    self%tprev = AnyDate()

    ! next time record ?
    if ( t > self%tprev ) self%itrec = self%itrec + 1

    ! store this time:
    self%tprev = t

    ! time since 2000-1-1 00:00
    t0 = NewDate( time6=(/2000,01,01,00,00,00/) )
    time = rTotal( t - t0, 'day' )

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_pixel_id, omi_trc%pixel_id(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_lon, omi_trc%longitude(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_lat, omi_trc%latitude(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! date in year,month,etc
    status = NF90_Put_Var( self%ncid, self%varid_omi_date, omi_trc%date(:,1:nmeas), &
                                  start=(/1,istart/), count=(/6,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! convert date to days since ...
    do imeas = 1, nmeas
       omi_trc_t = NewDate( time6=omi_trc%date(:,imeas) )
       omi_trc%time(imeas) = rTotal( omi_trc_t - t0, 'day' )
    end do
    status = NF90_Put_Var( self%ncid, self%varid_omi_time, omi_trc%time(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_omi_pressure_levels, omi_trc%pressure_levels(:,1:nmeas), &
                                  start=(/1,istart/), count=(/omi_trc%nlev+1,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_omi_kernel         , omi_trc%kernel(:,1:nmeas)         , &
                                  start=(/1,istart/), count=(/omi_trc%nlev,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_omi_cloud_radiance_fraction , omi_trc%cloud_radiance_fraction(1:nmeas) , &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_y, omi_trc%vcd_trop(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_r, omi_trc%sigma_vcd_trop(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( self%ncid, self%varid_omi_status, omi_trc%status(1:nmeas), &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    allocate ( dummy(nmeas) )
    dummy(:) = 0
    do imeas = 1, nmeas
       if ( omi_trc%analysed(imeas) ) dummy(imeas) = 1
    end do
    status = NF90_Put_Var( self%ncid, self%varid_omi_analysed, dummy, &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    deallocate ( dummy )

    allocate ( dummy(nmeas) )
    dummy(:) = 0
    do imeas = 1, nmeas
       if ( omi_trc%screened(imeas) ) dummy(imeas) = 1
    end do
    status = NF90_Put_Var( self%ncid, self%varid_omi_screened, dummy, &
                                  start=(/istart/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    deallocate ( dummy )

    ! current time up to seconds:
    call Get( t, time6=time6 )

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_le_altitude_levels, omi_trc%le_altitude_levels(:,1:nmeas), &
         start=(/1,istart,self%itrec/), count=(/nz+nhz+1,nmeas,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( self%ncid, self%varid_le_pressure_levels, omi_trc%le_pressure_levels(:,1:nmeas), &
         start=(/1,istart,self%itrec/), count=(/nz+nhz+1,nmeas,1/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! loop over filter moments
    do ifm = 1, nfm

       ! loop over filter states:
       do ifs = 1, nfs

          ! filter state description, skip if possible
          select case ( trim(fsname(ifs)) )
          case ( 'xb' )
             ! skip ?
             if ( .not. kf_with_xb ) cycle
             if ( ifm == ifm_forecast ) cycle
          case ( 'x' )
             ! skip ?
             if ( .not. kf_with_xm ) cycle
          case ( 's' )
             ! skip ?
             if ( .not. kf_with_xm ) cycle
          case default
             write (gol,'("unsupported filter state : ",a)') trim(fsname(ifs)); call goErr
             TRACEBACK; status=1; return
          end select

          ! write:
          allocate( le2omi_trc_vcd_trop(1:nmeas), stat=status )
          IF_NF90_NOTOK_RETURN(status=1)
          le2omi_trc_vcd_trop(:) = omi_trc%le2omi_trc_vcd_trop(1:nmeas,ifs,ifm)
          status = NF90_Put_Var( self%ncid, self%varid(ient_le2omi_vcd_trop,ifs,ifm), le2omi_trc_vcd_trop(1:nmeas), &
               start=(/istart,self%itrec/), count=(/nmeas,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
          deallocate( le2omi_trc_vcd_trop, stat=status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Var( self%ncid, self%varid(ient_le2omi_vcd,ifs,ifm), omi_trc%le2omi_trc_vcd_prof(:,1:nmeas), &
               start=(/1,istart,self%itrec/), count=(/nlev,nmeas,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Var( self%ncid, self%varid(ient_le_vcd,ifs,ifm), omi_trc%le_trc_vcd_prof(:,1:nmeas), &
               start=(/1,istart,self%itrec/), count=(/nz+nhz,nmeas,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Var( self%ncid, self%varid(ient_le_vmr,ifs,ifm), omi_trc%le_trc_vmr_prof(:,1:nmeas), &
               start=(/1,istart,self%itrec/), count=(/nz+nhz,nmeas,1/) )
          IF_NF90_NOTOK_RETURN(status=1)

       end do  ! filter states

    end do  ! filter moments

    ! ok
    status = 0

  end subroutine meas_output_omi_PutOut


end module LEKF_Meas_OMI_TRC
