!###############################################################################
!
! NAME
!
!   le_output_omi_trc  -  interface to OMI TRC data
!
! Interface:
!
! o call LE_Output_OMI_TRC_Init( rcfile, status )
!       Set relevant parameters (based on rc file) and initialise
!   call LE_Output_OMI_TRC_Done( status )
!
!    	Finish and release storage if needed
!
! o call LE_Output_OMI_TRC_GetMeas( tint1, tint2, status )
!
!     Make measurements between time "tint1" and "tint2" available to assimilation (read the OMI file when needed)
!
! o call CalcMeas_OMI_TRC( le_vmr_grid, le_vmr_aloft, le2omi_vcd_trop, status )
!	    Simulate OMI TRC obs "le2omi_vcd_trop" given the model concentrations "le_vmr_grid" and "le_vmr_aloft" 
!     (only for observations in time range set with LE_Output_OMI_TRC_GetMeas)
!
! o call MeasUpdate_OMI_TRC( t1, t2, nmodes, status )
!     Analyse observations between times "t1" and "t2" for "nmodes" modes
!
! Output routines:
!   call Init( leo, rcfile, rckey, status )
!	Initialise OMI output
!   call Done( leo, status )
!	Finalise
!   call PutOut( leo, key, t, status )
!	Write netcdf file with OMI data and model equivalent values
!	Routine should be called every timestep
!	Output written only once per day
!	Note: PutOut will call "Calcmeas_OMI_TRC" !
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
!            float cloud_fraction(fakeDim16) ;
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
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################
  
module le_output_omi_trc

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use NetCDF, only : NF90_StrError, NF90_NOERR
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  use dims,      only : nx, ny

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  LE_Output_OMI_TRC_Init, LE_Output_OMI_TRC_Done
  public  ::  LE_Output_OMI_TRC_GetMeas
  public  ::  CalcMeas_OMI_TRC
  
  public  ::  T_LE_Output_OMI_TRC
  public  ::  Init, Done, PutOut
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'le_output_omi_trc'

  ! maximum number of omi_trc measurements:
  integer, parameter :: maxomi = 25000

  ! which tracer entities simulated from state should be put out ?
  integer, parameter            ::  nent = 4
  integer, parameter            ::  ient_le2omi_vcd_trop = 1
  integer, parameter            ::  ient_le2omi_vcd      = 2
  integer, parameter            ::  ient_le_vcd          = 3
  integer, parameter            ::  ient_le_vmr          = 4
  character(len=15), parameter   ::  entname(nent) = &
       (/'le2omi_vcd_trop','le2omi_vcd     ','le_vcd         ','le_vmr         '/)
  
  ! Status number is sum of one or more numbers 2^n,
  ! for example   1+4+8
  ! To test if it is composed of a certain number, use:
  !   if ( iand(status,4) /= 0 )  ..
  !   
  integer, parameter            ::  status_default           = 0
  integer, parameter            ::  status_outside           = 1
  integer, parameter            ::  status_nodata            = 2
  integer, parameter            ::  status_validation        = 4
  character(len=*), parameter   ::  status_description = &
       'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation'


  ! --- types --------------------------------
    
  ! output stuff
  
  type T_LE_Output_OMI_TRC
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! state name:
    character(len=16)           ::  state
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! replace existing files
    logical                     ::  replace
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! file name:
    character(len=512)          ::  fname
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
    ! meta variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_ilon
    integer                     ::  varid_ilat
    integer                     ::  varid_omi_time
    integer                     ::  varid_omi_date
    integer                     ::  varid_omi_pressure_levels
    integer                     ::  varid_omi_kernel
    integer                     ::  varid_omi_cloud_fraction
    integer                     ::  varid_omi_status
    integer                     ::  varid_omi_analysed
    integer                     ::  varid_omi_screened
    integer                     ::  varid_le_altitude_levels
    integer                     ::  varid_le_pressure_levels
    ! tracer variables:
    !real                        ::  unitconv
    integer                     ::  varid_y
    integer                     ::  varid_r
    integer                     ::  varid(nent)

  end type T_LE_Output_OMI_TRC

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
    real, pointer        ::  time(:)                  ! time of OMI obs
    real, pointer        ::  longitude(:) 
    real, pointer        ::  latitude(:)
    real, pointer        ::  vcd_trop(:)              ! OMI tropospheric column
    real, pointer        ::  sigma_vcd_trop(:)        ! and uncertainty
    real, pointer        ::  kernel(:,:)              ! Averaging kernel vector 
    real, pointer        ::  pressure_levels(:,:)     ! Pressure levels on which kernel is defined
    real, pointer        ::  cloud_fraction(:)        ! Cloud fraction of OMI obs
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
    real, pointer        ::  le2omi_trc_vcd_trop(:)   ! LE estimate of OMI obs (tropospheric column)
  end type OMI_TRC_type


  ! --- interfaces -------------------------
  
  interface Init
     module procedure leo_Init
  end interface
  
  interface Done
     module procedure leo_Done
  end interface
  
  interface PutOut
     module procedure leo_PutOut
  end interface
  
  ! --- var ----------------------------------

  ! archive description:
  character(len=256)        ::  arch_filenames
  character(len=256)        ::  arch_file
  
  ! Subset of observations for Kalman time step
  type(OMI_TRC_type)        ::  omi_trc
  
  ! time range:
  type(TDate)               ::  omi_trc_tr(2)  ! time range with valid data (24 hour)

  ! top of troposphere to be used for kernels:
  real                      ::  omi_trc_ptropo

  ! LE species index:
  integer                   ::  omi_trc_itracer

  
contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================
  
  
  ! opens the measurement file (if present)
  ! end sets some parameters
  
  subroutine LE_Output_OMI_TRC_Init( rcF, status )

    use GO          , only : TrcFile
    use GO          , only : TDate, NewDate, Get, AnyDate
    use GO          , only : goMatchValue
    use Indices     , only : specname
    use LE_Data, only : LE_Data_Enable

    ! --- in/out --------------------------------

    type(TRcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_OMI_TRC_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    write (gol,'("    OMI/TRC: initialise ...")') ; call goPr

    ! ~~ settings from .rc file

    ! tracer:
    call rcF%Get( 'le.output.omi_trc.tracer', omi_trc%tracer_name, status )
    IF_NOTOK_RETURN(status=1)
    ! search index:
    call goMatchValue( trim(omi_trc%tracer_name), specname, omi_trc_itracer, status )
    IF_NOTOK_RETURN(status=1)
    
    ! filename template:
    call rcF%Get( 'le.output.omi_trc.filenames', arch_filenames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! top of tropo sphere to be used with kernel (Pa):
    call rcF%Get( 'le.output.omi_trc.ptropo', omi_trc_ptropo, status )
    IF_NOTOK_RETURN(status=1)
    
    call rcF%Get( 'le.output.omi_trc.top_hp', omi_trc%top_ph_entry, status )
    IF_NOTOK_RETURN(status=1)
    
    call rcF%Get( 'le.output.omi_trc.top_trc', omi_trc%top_trc_entry, status )
    IF_NOTOK_RETURN(status=1)
    
    ! nothing stored yet ...
    omi_trc%nmeas = 0

    ! no data read yet:
    omi_trc_tr(1) = AnyDate()
    omi_trc_tr(2) = AnyDate()

    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'oro', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_OMI_TRC_Init


  ! ***


  ! read station locations etc

  subroutine LE_Output_OMI_TRC_Done( status )

    !use kf_meas, only : Done

    ! --- in/out --------------------------------

    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_OMI_TRC_Done'

    ! --- begin ------------------------------

    ! remove storage for OMI data
    call Dealloc_OMI_TRC( omi_trc, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Output_OMI_TRC_Done


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
       allocate( omitrc%latitude(nrmeas)                     )
       allocate( omitrc%vcd_trop(nrmeas)                     )
       allocate( omitrc%sigma_vcd_trop(nrmeas)               )
       allocate( omitrc%kernel(nrlev,nrmeas)                 )
       allocate( omitrc%pressure_levels(0:nrlev,nrmeas)      )
       allocate( omitrc%cloud_fraction(nrmeas)               )
       allocate( omitrc%ix(nrmeas)                           )
       allocate( omitrc%iy(nrmeas)                           )
       allocate( omitrc%status(nrmeas)                       )
       allocate( omitrc%analysed(nrmeas)                     )
       allocate( omitrc%screened(nrmeas)                     )
       allocate( omitrc%intimestep(nrmeas)                   )
       allocate( omitrc%le_altitude_levels(0:nz,nrmeas)    )
       allocate( omitrc%le_pressure_levels(0:nz,nrmeas)    )
       allocate( omitrc%le_trc_vmr_prof(1:nz,nrmeas)       )
       allocate( omitrc%le_trc_vcd_prof(1:nz,nrmeas)       )
       allocate( omitrc%le2omi_trc_vcd_prof(nrlev,nrmeas)    )
       allocate( omitrc%le2omi_trc_vcd_trop(nrmeas)          )

       ! reset arrays:
       omi_trc%le2omi_trc_vcd_trop = -999.9
       omi_trc%le2omi_trc_vcd_prof = -999.9
       omi_trc%le_trc_vcd_prof     = -999.9
       omi_trc%le_trc_vmr_prof     = -999.9

    else
       write (gol,'("    OMI/TRC: ERROR: arrays already allocated ?!?")') ; call goErr
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
       deallocate( omitrc%vcd_trop            )
       deallocate( omitrc%sigma_vcd_trop      )
       deallocate( omitrc%kernel              )
       deallocate( omitrc%pressure_levels     )
       deallocate( omitrc%cloud_fraction      )
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
       write (gol,'("    OMI/TRC: WARNING Dealloc_OMI_TRC: Arrays are not allocated ?!?")') ; call goPr       
    end if

    ! ok
    status = 0

  end subroutine Dealloc_OMI_TRC


  ! ***

  subroutine LE_Output_OMI_TRC_GetMeas( tint1, tint2, status )

    use GO,      only : TDate, NewDate, IncrDate, Get, wrtgol
    use GO,      only : operator(<), operator(<=), operator(+), operator(-), operator(/)
    use GO,      only : goReplace, goUpCase, goLoCase

    
    use LE_Grid, only : ugg
    !use kf_dims, only : maxomi

    use NetCDF,  only : NF90_NOWRITE
    use NetCDF,  only : NF90_Open, NF90_Close
    use NetCDF,  only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF,  only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var

    ! --- in/out ----------------------------

    type(TDate), intent(in)   ::  tint1, tint2
    integer, intent(out)      ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_OMI_TRC_GetMeas'

    ! --- local -----------------------------

    integer             ::  year, month, day
    type(TDate)         ::  omi_trc_t  
    character(len=256)  ::  fname
    logical             ::  exist

    integer             ::  ncid, varid
    integer             ::  dimids(2)

    integer             ::  nrmeas, nrlev, imeas, i, nrmeas_in_domain
    !integer             ::  iobs, count_timestep
    Logical             ::  indomain

    ! --- begin -----------------------------

    call wrtgol( '    OMI/TRC: LE_Output_OMI_TRC_GetMeas - Get data for time from ',tint1,' to ',tint2); call goPr 

    ! no data for this time yet ?
    if ( (tint1 < omi_trc_tr(1)) .or. (omi_trc_tr(2) < tint1) ) then

       ! already something loaded ?
       if ( omi_trc%nmeas > 0 ) then
          ! clear OMI storage:
          call Dealloc_OMI_TRC( omi_trc, status )
          IF_NOTOK_RETURN(status=1)
       end if

       ! extract time values:
       call Get( tint1, year=year, month=month, day=day )

       ! new time range:
       omi_trc_tr(1) = NewDate( year=year, month=month, day=day )
       omi_trc_tr(2) = omi_trc_tr(1) + IncrDate(day=1)

       ! info ...
       write (gol,'("    OMI/TRC: open file for ",i4,2i2.2)') year, month, day; call goPr

       ! start filename with template:
       fname = trim(arch_filenames)
       ! replace some values if necessary:
       call goReplace( fname, '%{year}', '(i4.4)', year, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{month}', '(i2.2)', month, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{yyyy}', '(i4.4)', year, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{mm}', '(i2.2)', month, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{yyyymmdd}', '(i8.8)', year*10000+month*100+day, status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{TRACER}', goUpCase(trim(omi_trc%tracer_name)), status )
       IF_ERROR_RETURN(status=1)
       call goReplace( fname, '%{tracer}', goLoCase(trim(omi_trc%tracer_name)), status )
       IF_ERROR_RETURN(status=1)

       ! info ...
       write (gol,'("    OMI/TRC: read from file ",a)') fname; call goPr

       ! available ?
       inquire( file=trim(fname), exist=exist )
       if ( .not. exist ) then
          !write (gol,'("    OMI/TRC: file not found : ",a)') trim(fname); call goErr
          !TRACEBACK; status=1; return
          ! info ...
          write (gol,'("    OMI/TRC: WARNING - no OMI TRC data file ...")'); call goPr
          ! set to empty:
          omi_trc%nmeas = 0
          ! done:
          status=0; return
       end if

       ! open file:
       status = NF90_Open( fname, NF90_NOWRITE, ncid )
       IF_NF90_NOTOK_RETURN(status=1)

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
       ! number of pixels:
       status = NF90_Inquire_Dimension( ncid, dimids(2), len=nrmeas )
       IF_NF90_NOTOK_RETURN(status=1)
       if ( nrmeas > maxomi ) then
          write (gol,'("number of pixels in omi trc file exceeds maximum : ")'); call goErr
          write (gol,'("  file    : ",i6)') nrmeas; call goErr
          write (gol,'("  maxomi  : ",i6)') maxomi; call goErr
          TRACEBACK; status=1; return
       end if

       ! info ...
       !write (gol,'("    OMI/TRC: nr of kernel levels : ",i6)') nrlev; call goPr
       !write (gol,'("    OMI/TRC: nr of OMI obs       : ",i6)') nrmeas; call goPr

       ! create storage
       call Alloc_OMI_TRC( omi_trc, nrmeas, nrlev, status )

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
       status = NF90_Inq_Varid( ncid, 'cloud_fraction', varid )
       if ( status /= 0 ) then
         status = NF90_Inq_Varid( ncid, 'cloud_radiance_fraction', varid)
         IF_NF90_NOTOK_RETURN(status=1)
       end if
       status = NF90_Get_Var( ncid, varid, omi_trc%cloud_fraction)
       IF_NF90_NOTOK_RETURN(status=1)

       ! close:
       status = NF90_Close( ncid )
       IF_NF90_NOTOK_RETURN(status=1)

       ! info ...
       !write (gol,'("    OMI/TRC: file closed ")') ; call goPr

       ! default field values: status etc
       omi_trc%status(:)   = 0
       omi_trc%analysed(:) = .false.
       omi_trc%screened(:) = .false.

       nrmeas_in_domain = 0

       ! interpolation from LE grid cells:
       do imeas = 1, nrmeas

          ! * check time:
          omi_trc_t = NewDate( time6=omi_trc%date(:,imeas) )
          if ( ( omi_trc_t < omi_trc_tr(1) ) .or. ( omi_trc_tr(2) < omi_trc_t ) ) then
             ! write (gol,'("observation time outside expected bounds:")'); call goErr
             ! write (gol,'("  file name    : ",a)') trim(fname); call goErr
             ! write (gol,'("  observation  : ",i6)') imeas; call goErr
             ! call wrtgol( '  time         : ', omi_trc_t ); call goErr
             ! call wrtgol( ', expected     : ', omi_trc_tr(1), ' - ', omi_trc_tr(2) ); call goErr
             ! TRACEBACK; status=1; return
             !
             ! Henk Eskes: Note that in Summer months there can be some OMI observations from an other date,
             !    related to night-time orbits 
             !    for safety these obs are not used
             write (gol,'("WARNING - observation time outside expected bounds:")'); call goPr
             write (gol,'("  file name    : ",a)') trim(fname); call goPr
             write (gol,'("  observation  : ",i6)') imeas; call goPr
             call wrtgol( '  time         : ', omi_trc_t ); call goPr
             call wrtgol( ', expected     : ', omi_trc_tr(1), ' - ', omi_trc_tr(2) ); call goPr
             omi_trc%status(imeas) = omi_trc%status(imeas) + status_nodata
          end if

          ! * check if meas is in domain;
          ! return cell indices (ix,iy) and flag:
          call ugg%InDomain( omi_trc%longitude(imeas), omi_trc%latitude(imeas), indomain, status )
          IF_NOTOK_RETURN(status=1)
          ! not in domain ?
          if ( .not. indomain ) omi_trc%status(imeas) = omi_trc%status(imeas) + status_outside
          if ( indomain )  then
            nrmeas_in_domain = nrmeas_in_domain + 1
            call ugg%GetLocation( omi_trc%longitude(imeas), omi_trc%latitude(imeas), &
                                  omi_trc%ix(imeas), omi_trc%iy(imeas), status )
            IF_NOTOK_RETURN(status=1)
          end if

          ! * check value:
          ! HE: Do not exclude negative values!
          ! HE: In very clean areas there are negative values that should be included
          ! HE: since otherwise a bias towards positive values is introduced ....
          ! if ( omi_trc%vcd_trop(imeas) < 0.0 ) omi_trc%status(imeas) = omi_trc%status(imeas) + status_nodata
       end do

       ! info ...
       !write (gol,*) "    OMI/TRC: Nr of observations  = ", omi_trc%nmeas;     call goPr
       !write (gol,*) "    OMI/TRC: Nr of obs in domain = ", nrmeas_in_domain;     call goPr
       !write (gol,*) "    OMI/TRC: last OMI observation read : ";     call goPr
       !write (gol,*) "    OMI/TRC:   time ",omi_trc%time(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   lon  ",omi_trc%longitude(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   lat  ",omi_trc%latitude(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   vcd  ",omi_trc%vcd_trop(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   sig  ",omi_trc%sigma_vcd_trop(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   kern ",omi_trc%kernel(:,omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   pres ",omi_trc%pressure_levels(:,omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   clfr ",omi_trc%cloud_fraction(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   ix   ",omi_trc%ix(omi_trc%nmeas);     call goPr
       !write (gol,*) "    OMI/TRC:   iy   ",omi_trc%iy(omi_trc%nmeas);     call goPr

    else

       write (gol,*) "    OMI/TRC: data already available for this day";     call goPr

    end if

    ! no data ?
    if ( omi_trc%nmeas < 1 ) then
       write (gol,'("    OMI/TRC: GetMeas: skip, no observations")'); call goPr
       status = 0; return
    end if

    ! data from file has been read for the current day and is stored at this point

    ! info ...
    !write (gol,*) "    OMI/TRC: first time ",omi_trc%date(:,1);     call goPr
    !write (gol,*) "    OMI/TRC: last  time ",omi_trc%date(:,omi_trc%nmeas);     call goPr

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
       if ( (omi_trc_t <= tint1) .or. (tint2 < omi_trc_t) ) cycle

       ! skip if not accepted (validation, no data?):
       if ( omi_trc%status(i) /= status_default ) cycle
 
       ! flag to signal that this observation needs to be processed in this assimilation time step
       omi_trc%intimestep(i) = .true.

       ! count number of relevant obs
       omi_trc%nmeas_intimestep = omi_trc%nmeas_intimestep + 1

    end do

    write (gol,'("    OMI/TRC: - nr of valid obs in timestep : ",i6)') omi_trc%nmeas_intimestep ; call goPr

    ! ok
    status = 0

  end subroutine LE_Output_OMI_TRC_GetMeas


  ! ***


  ! simulate omi trc given the model concentrations

  subroutine CalcMeas_OMI_TRC( le_vmr_grid, le_vmr_aloft, le2omi_vcd_trop, status )

    use Binas, only : Avog, xm_air, grav
    use JAQL , only : PotentialPressure
    use Num  , only : IntervalSum

    use dims, only : nx, ny, nz, nspec
    use indices, only : specmolm
    
    use LE_Bound_Top, only : LE_Bound_Top_MergeProfiles
    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out -----------------------------

    real, intent(in)          ::  le_vmr_grid(nx,ny,nz)   ! volume ppb
    real, intent(in)          ::  le_vmr_aloft(nx,ny)     ! volume ppb
    real, intent(out)         ::  le2omi_vcd_trop(maxomi) ! unit 1e15 molec/cm2
    integer, intent(out)      ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/CalcMeas_OMI_TRC'

    ! --- local --------------------------------

    integer, parameter :: nlevmax = 100
    integer     ::  imeas
    integer     ::  ix, iy, iz
    integer     ::  ilev
    real, allocatable  ::  omi_ph(:)
    real, allocatable  ::  le_hh(:)
    real, allocatable  ::  le_ph(:)
    real, allocatable  ::  le_temp(:)
    real, allocatable  ::  le_rh(:)
    real, allocatable  ::  le_airm(:)
    real, allocatable  ::  le_vmr(:)
    real, allocatable  ::  le_vcd(:)
    real, allocatable  ::  le2omi_vcd(:)

    integer     ::  ilast, nlevx, nmeasx, nmeas_processed
    real        ::  le_ph_dum
    
    real                ::  du_to_kgm2
    real, pointer       ::  mer_phlev(:)
    real, pointer       ::  mer_vcd(:)
    
    ! meteo data:
    real, pointer       ::  h_m(:,:,:)   ! (lon,lat,1)
    real, pointer       ::  oro(:,:,:)   ! (lon,lat,1)
    real, pointer       ::  temp(:,:,:)   ! (lon,lat,1)    

    ! --- local --------------------------------
    
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'oro', oro, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units ='K')    
    IF_NOTOK_RETURN(status=1)

    ! allocate
    allocate( omi_ph(0:nlevmax) )
    allocate( le_hh(0:nz) )
    allocate( le_ph(0:nz) )
    allocate( le_temp(1:nz) )
    allocate( le_rh(1:nz) )
    allocate( le_airm(1:nz) )
    allocate( le_vmr(1:nz) )
    allocate( le_vcd(1:nz) )
    allocate( le2omi_vcd(1:nlevmax) )

    ! init
    nullify( mer_vcd )
    nullify( mer_phlev )

    le2omi_vcd_trop = 0.0

    ! no data ?
    if ( omi_trc%nmeas < 1 .or. omi_trc%nmeas_intimestep < 1 ) then
       write (gol,'("    OMI/TRC: CalcMeas: skip, no observations")'); call goPr
       status = 0; return
    end if

    write (gol,'("    OMI/TRC: CalcMeas: compute LE predicted OMI TRC, ",i6," obs ...")') omi_trc%nmeas_intimestep; call goPr

    nmeasx = omi_trc%nmeas
    nlevx = omi_trc%nlev
    nmeas_processed = 0

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
      le_hh(1:nz) = le_hh(0) + h_m(ix,iy,:)  ! m 

      ! model temperature:
      le_temp(1:nz) = temp(ix,iy,:)

      ! model relative humidity;
      ! assume zero since LE rh is probably not what you expect ...
      le_rh = 0.0

      ! model pressure profile:
      le_ph(0) = omi_ph(0)   ! use same surface pressure as omi
      do iz = 1, nz
         ! pressure decrease given layer heigt and temperature,
         ! assume zero relative humidity since LE rh is probably not
         ! what you expect ...
         call PotentialPressure( le_ph_dum, le_ph(iz-1), le_temp(iz), le_rh(iz), le_hh(iz)-le_hh(iz-1) )
         le_ph(iz) = le_ph_dum
      end do

      ! LE air mass:
      do iz = 1, nz
         le_airm(iz) = ( le_ph(iz-1) - le_ph(iz) )/grav  ! (kg air)/m2
      end do

      ! le trc vmr:
      le_vmr(1:nz) = le_vmr_grid (ix,iy,:)      ! trc volume ppb

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

      ! convert from "(1e15 mlc)/cm2" to "kg/m2" :
      ! (1e15 mlc)/cm2 * 1e15 / (mlc/mol) *       (kg/mol)            * cm2/m2
      du_to_kgm2       = 1e15 /   Avog    * specmolm(omi_trc_itracer) * 1.0e4
      ! get parent model profile; return status:
      !    0 : ok, profile extracted
      !   -1 : no top field stored for requested tracer
      !   -2 : top field defined, but not filled
      call LE_Bound_Top_MergeProfiles(  ix, iy, omi_trc_itracer, &
                                         omi_trc%top_ph_entry, omi_trc%top_trc_entry, &
                                         le_ph(0:nz), le_vcd(1:nz)*du_to_kgm2, 'kg/m2', &
                                         mer_phlev, mer_vcd, status )
      IF_NOTOK_RETURN(status=1)
      ! convert to units of original model profile:
      mer_vcd = mer_vcd / du_to_kgm2

      ! LE trc profile at OMI layers:
      ilast = 1
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
              print *, 'ERROR - from IntervalSum:'
              print *, 'ERROR - imeas   = ', imeas, ' / ', nmeasx
              print *, 'ERROR - loc     = ', omi_trc%longitude(imeas), omi_trc%latitude(imeas)
              print *, 'ERROR - cell    = ', omi_trc%ix(imeas), omi_trc%iy(imeas)
              print *, 'ERROR - ilev    = ', ilev
              print *, 'ERROR - le_hh   = ', le_hh
              print *, 'ERROR - le_ph   = ', le_ph
              print *, 'ERROR - le_arim = ', le_airm
              print *, 'ERROR - le_vmr  = ', le_vmr
              print *, 'ERROR - le_vcd  = ', le_vcd
              print *, 'ERROR - nlevx   = ', nlevx
              print *, 'ERROR - omi_ph  = ', omi_ph(0:nlevx)
              !print *, 'ERROR - x       = ', x
              print *, 'ERROR - dp      = ', -1.0*omi_ph(ilev-1), -1.0*omi_ph(ilev)
            end if
            IF_NOTOK_RETURN(status=1)
            
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

    end do  ! omi meas

    if ( omi_trc%nmeas_intimestep /= nmeas_processed ) then
       write (gol,'("    OMI/TRC: WARNING CalcMeas: In total ",i6," obs processed")') nmeas_processed ; call goPr
    end if

    ! deallocate
    deallocate( omi_ph )
    deallocate( le_hh )
    deallocate( le_ph )
    deallocate( le_temp )
    deallocate( le_rh )
    deallocate( le_airm )
    deallocate( le_vmr )
    deallocate( le_vcd )
    deallocate( le2omi_vcd )
    
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


  subroutine leo_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : AnyDate
    use LE_Output_Common, only : Init

    ! --- in/out --------------------------------

    type(T_LE_Output_OMI_TRC), intent(out)  ::  leo
    type(TRcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rckey
    character(len=*), intent(in)            ::  typ
    character(len=*), intent(in)            ::  name
    character(len=*), intent(in)            ::  state
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Init'

    ! --- local ----------------------------------

    ! --- begin ---------------------------------

    ! store file name keys:
    leo%typ   = typ
    leo%name  = name
    ! store state name:
    leo%state = trim(state)
    
    ! init common stuff:
    call Init( leo%com, rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)
 
    ! replace existing files ?
    call rcF%Get( trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! files not open yet:
    leo%opened = .false.

    ! no time range set yet:
    leo%tr(1) = AnyDate()
    leo%tr(2) = AnyDate()
    
    ! ok
    status = 0

  end subroutine leo_Init


  ! ***


  subroutine leo_Done( leo, status )

    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_OMI_TRC), intent(inout)      ::  leo
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Done'

    ! --- begin ---------------------------------

    ! file opened ?
    if ( leo%opened ) then
       ! close:
       status = NF90_Close( leo%ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! reset flag:
       leo%opened = .false.
    end if

    ! done with common stuff ...
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine leo_Done


  ! ***


  ! Output is generate only once per day, after all OMI data has been processed.
  ! The subroutine will accumulate the relevant information to print until
  ! all data has been analysed

  subroutine leo_PutOut( leo, t, c, status, last )

    use GO, only : TDate, IncrDate, Get, NewDate, AnyDate
    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal
    use GO, only : Precisely, wrtgol

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE

    use Dims   , only : nx, ny, nz, nspec
    use LE_Bound        , only : caloft
    use LE_Output_Common, only : PutOut_GlobalAttributes

    ! --- in/out --------------------------------

    type(T_LE_Output_OMI_TRC), intent(inout)  ::  leo
    type(TDate), intent(in)                   ::  t
    real, intent(in)                          ::  c(nx,ny,nz,nspec)
    integer, intent(out)                      ::  status
    logical, intent(in), optional             ::  last

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_meas_PutOut'

    ! --- local ---------------------------------

    !logical               ::  newfile    
    logical               ::  writetofile
    type(TDate)           ::  omi_trc_t
    integer               ::  time6(6)
    real                  ::  time

    integer               ::  cmode
    integer               ::  varid
    type(TDate)           ::  t0

    character(len=32)     ::  varname
    character(len=512)    ::  descr
    integer               ::  imeas   
    integer               ::  ient
    integer               ::  nlev, i
    ! filter pixels:
    integer               ::  nmeas
    integer, allocatable  ::  ii(:)

    real, allocatable     ::  le2omi_trc_vcd_trop(:)
    integer, allocatable  ::  dummy(:)

    ! --- begin ---------------------------------

    ! no measurements ? then leave
    if ( omi_trc%nmeas < 1 .or. omi_trc%nmeas_intimestep < 1 ) then
       ! info ...
       call wrtgol('    OMI/TRC: put out return: no TRC obs for time : ',t); call goPr
       status=0; return
    end if

    ! storage:
    allocate( le2omi_trc_vcd_trop(1:maxomi) )

    ! simulate omi stuff:
    call CalcMeas_OMI_TRC( c(:,:,:,omi_trc_itracer), caloft(:,:,omi_trc_itracer), le2omi_trc_vcd_trop, status )
    IF_NOTOK_RETURN(status=1)

    ! Store the fields for all filter states
    do i = 1, omi_trc%nmeas
      ! First check if this obs is relevant for this time step
      if ( omi_trc%intimestep(i) ) then
        omi_trc%le2omi_trc_vcd_trop(i) = le2omi_trc_vcd_trop(i)
      end if
    end do

    ! storage:
    deallocate( le2omi_trc_vcd_trop )

    ! check if all valid data has been filled ?
    writetofile = omi_trc%nodataleft
  
    ! optional argument provided ?  
    if ( present(last) ) then
      ! last call to putout routines ?
      if ( last ) then
        ! info if write was originally not intended ...
        if ( .not. writetofile ) then
          write (gol,'("    OMI/TRC: this is the last output, so write although not all pixels were simulated yet ...")'); call goPr
        end if
        ! reset flag:
        writetofile = .true.
      end if  ! last putout ?
    end if  ! present last

    ! leave here ?
    if ( .not. writetofile ) then
      status=0; return
    end if

    ! ------------------------------------------------------------------
    ! Now write to file because all status=0 OMI data has been analysed
    ! ------------------------------------------------------------------
    
    ! filter:
    allocate( ii(1:omi_trc%nmeas) )
    nmeas = 0
    do i = 1, omi_trc%nmeas
      ! skip pixels outside domain:
      if ( iand(omi_trc%status(i),status_outside) /= 0 ) cycle
      ! skip pixels without data:
      if ( iand(omi_trc%status(i),status_nodata) /= 0 ) cycle
      ! add index:
      nmeas = nmeas + 1
      ii(nmeas) = i
    end do
    write (gol,'("    OMI/TRC: pixels in domain and with data : ",i6," / ",i6)') nmeas, omi_trc%nmeas; call goPr
    ! empty ?
    if ( nmeas < 1 ) then
      write (gol,'("    OMI/TRC: no pixels that survived the filter ...")'); call goPr
      deallocate( ii )
      status=0; return
    end if

    ! close old file if necessary (should not be necessary)
    if (  leo%opened ) then
       ! close:
       status = NF90_Close( leo%ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! reset flag:
       leo%opened = .false.
    end if

    ! info ...
    call wrtgol('    OMI/TRC: Now put out all data for this day at : ',t); call goPr

    ! extract time fields:
    call Get( t, time6=time6 )
    
    ! time range for this file is (00,24]
    leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
    leo%tr(2) = leo%tr(1) + IncrDate( day=1 )

    ! new file name:
    write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
         trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
         trim(leo%name), time6(1:3)

    ! set creation mode flag:
    if ( leo%replace ) then
      cmode = NF90_CLOBBER       ! overwrite existing files
    else
      cmode = NF90_NOCLOBBER     ! do not overwrite existing files
    end if

    ! create file:
    status = NF90_Create( leo%fname, cmode, leo%ncid )
    if ( status /= 0 ) then
       write (gol,'("creating file :")'); call goErr
       write (gol,'("  ",a)') trim(leo%fname); call goErr
       TRACEBACK; status=1; return
    end if

    ! reset flag:
    leo%opened = .true.

     print *, 'meas-omi name created ', trim(leo%fname)


    ! write global attributes:
    call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! define dimensions:
    status = NF90_Def_Dim( leo%ncid, 'omi_obs', nmeas, leo%dimid_omi_obs )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( leo%ncid, 'omi_lev' , omi_trc%nlev, leo%dimid_omi_lev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( leo%ncid, 'omi_hlev' , omi_trc%nlev+1, leo%dimid_omi_hlev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( leo%ncid, 'le_lev' , nz, leo%dimid_le_lev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( leo%ncid, 'le_hlev' , nz+1, leo%dimid_le_hlev )
    IF_NF90_NOTOK_RETURN(status=1)
    !       status = NF90_Def_Dim( leo%ncid, 'time', NF90_UNLIMITED, leo%dimid_time )
    !       IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( leo%ncid, 'datelen', 6, leo%dimid_datelen )
    IF_NF90_NOTOK_RETURN(status=1)

    ! define daily constant variables:
    status = NF90_Def_Var( leo%ncid, 'lon', NF90_REAL, leo%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'degrees_east' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_lon = varid

    status = NF90_Def_Var( leo%ncid, 'lat', NF90_REAL, leo%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'degrees_north' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_lat = varid

    ! define daily constant variables:
    status = NF90_Def_Var( leo%ncid, 'ilon', NF90_REAL, leo%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'index', 'model_gridcell_longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_ilon = varid

    status = NF90_Def_Var( leo%ncid, 'ilat', NF90_REAL, leo%dimid_omi_obs, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'index', 'model_gridcell_latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_ilat = varid

    status = NF90_Def_Var( leo%ncid, 'omi_time', NF90_REAL, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'days since 2000-01-01 00:00:00' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'calendar', 'standard' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_time = varid

    status = NF90_Def_Var( leo%ncid, 'omi_date', NF90_INT, (/leo%dimid_datelen,leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI date and time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'year, month, day, hour, minute, second' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_date = varid

    status = NF90_Def_Var( leo%ncid, 'omi_pressure_levels', NF90_REAL, (/leo%dimid_omi_hlev,leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI pressure halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_pressure_levels = varid

    status = NF90_Def_Var( leo%ncid, 'omi_kernel', NF90_REAL, (/leo%dimid_omi_lev,leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI kernel' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', '1' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_kernel = varid

    status = NF90_Def_Var( leo%ncid, 'omi_cloud_fraction', NF90_REAL, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI cloud fraction' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', '[0,1]' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_cloud_fraction = varid

    status = NF90_Def_Var( leo%ncid, 'omi_status', NF90_INT, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI filter status' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', trim(status_description) )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_status = varid

    status = NF90_Def_Var( leo%ncid, 'omi_analysed', NF90_INT, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI analysed flag' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'flag: 0 = false, 1 = true' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_analysed = varid

    status = NF90_Def_Var( leo%ncid, 'omi_screened', NF90_INT, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI screening flag' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'flag: 0 = false, 1 = true' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_omi_screened = varid

    status = NF90_Def_Var( leo%ncid, 'omi_vcd_trop_y', NF90_REAL, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI TRC tropospheric column density' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', '(1e15 mlc TRC)/cm2' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_y = varid

    status = NF90_Def_Var( leo%ncid, 'omi_vcd_trop_r', NF90_REAL, (/leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'OMI TRC tropospheric column density error std.dev.' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', '(1e15 mlc TRC)/cm2' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_r = varid

    ! define temporal variables:

    status = NF90_Def_Var( leo%ncid, 'le_altitude_levels', NF90_REAL, (/leo%dimid_le_hlev,leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS altitude halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'm' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_le_altitude_levels = varid

    status = NF90_Def_Var( leo%ncid, 'le_pressure_levels', NF90_REAL, (/leo%dimid_le_hlev,leo%dimid_omi_obs/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS pressure halflevels' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( leo%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOTOK_RETURN(status=1)
    leo%varid_le_pressure_levels = varid

    ! state description:
    descr = 'model background run'

    ! variable name:
    varname = 'conc'

    ! loop over tracer entities to be written:
    do ient = 1, nent

      ! variable name:
      varname = trim(entname(ient))

      ! define variables for each entity:
      select case ( ient )

        case ( ient_le2omi_vcd_trop )

          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, (/leo%dimid_omi_obs/), varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS TRC tropospheric column density' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'units', '(1e15 mlc TRC)/cm2' )
          IF_NF90_NOTOK_RETURN(status=1)

        case ( ient_le2omi_vcd )

          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, (/leo%dimid_omi_lev,leo%dimid_omi_obs/), varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS TRC vertical column density at OMI layers' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'units', '(1e15 mlc TRC)/cm2' )
          IF_NF90_NOTOK_RETURN(status=1)

        case ( ient_le_vcd )

          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, (/leo%dimid_le_lev,leo%dimid_omi_obs/), varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS TRC vertical column density' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'units', '(1e15 mlc TRC)/cm2' )
          IF_NF90_NOTOK_RETURN(status=1)

        case ( ient_le_vmr )

          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, (/leo%dimid_le_lev,leo%dimid_omi_obs/), varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'LOTOS-EUROS TRC volume-mixing-ratio' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( leo%ncid, varid, 'units', 'ppb' )
          IF_NF90_NOTOK_RETURN(status=1)

        case default

          write (gol,'("unsupported entity : ",i6)') ient; call goErr
          TRACEBACK; status=1; return

      end select

      ! standard attributes:
      status = nf90_put_att( leo%ncid, varid, 'description', trim(descr) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! store variable id:
      leo%varid(ient) = varid

   end do  ! output entities

    ! end defintion mode:

    status = NF90_EndDef( leo%ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! time since 2000-1-1 00:00
    t0 = NewDate( time6=(/2000,01,01,00,00,00/) )
    time = rTotal( t - t0, 'day' )
    
    status = NF90_Put_Var( leo%ncid, leo%varid_lon, omi_trc%longitude(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_lat, omi_trc%latitude(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_ilon, omi_trc%ix(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_ilat, omi_trc%iy(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! date in year,month,etc
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_date, omi_trc%date(:,ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! convert date to days since ...
    do imeas = 1, omi_trc%nmeas
       omi_trc_t = NewDate( time6=omi_trc%date(:,imeas) )
       omi_trc%time(imeas) = rTotal( omi_trc_t - t0, 'day' )
    end do
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_time, omi_trc%time(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_pressure_levels, omi_trc%pressure_levels(:,ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_kernel         , omi_trc%kernel(:,ii(1:nmeas))          )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_cloud_fraction , omi_trc%cloud_fraction(ii(1:nmeas))  )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_y, omi_trc%vcd_trop(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( leo%ncid, leo%varid_r, omi_trc%sigma_vcd_trop(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Put_Var( leo%ncid, leo%varid_omi_status, omi_trc%status(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)

    allocate ( dummy(omi_trc%nmeas) )
    dummy(:) = 0
    do imeas = 1, omi_trc%nmeas 
       if ( omi_trc%analysed(imeas) ) dummy(imeas) = 1
    end do
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_analysed, dummy(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    deallocate ( dummy )

    allocate ( dummy(omi_trc%nmeas) )
    dummy(:) = 0
    do imeas = 1, omi_trc%nmeas 
       if ( omi_trc%screened(imeas) ) dummy(imeas) = 1
    end do
    status = NF90_Put_Var( leo%ncid, leo%varid_omi_screened, dummy(ii(1:nmeas)) )
    IF_NF90_NOTOK_RETURN(status=1)
    deallocate ( dummy )

    ! current time up to seconds:
    call Get( t, time6=time6 )

    !nmeas = omi_trc%nmeas  ! <-- now filtered !
    nlev = omi_trc%nlev

    ! write:
    status = NF90_Put_Var( leo%ncid, leo%varid_le_altitude_levels, &
         omi_trc%le_altitude_levels(:,ii(1:nmeas)), &
         start=(/1,1/), count=(/nz+1,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( leo%ncid, leo%varid_le_pressure_levels, &
         omi_trc%le_pressure_levels(:,ii(1:nmeas)), &
         start=(/1,1/), count=(/nz+1,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! write:
    status = NF90_Put_Var( leo%ncid, leo%varid(ient_le2omi_vcd_trop), &
         omi_trc%le2omi_trc_vcd_trop(ii(1:nmeas)), &
         start=(/1/), count=(/nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( leo%ncid, leo%varid(ient_le2omi_vcd), &
         omi_trc%le2omi_trc_vcd_prof(:,ii(1:nmeas)), &
         start=(/1,1/), count=(/nlev,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( leo%ncid, leo%varid(ient_le_vcd), &
         omi_trc%le_trc_vcd_prof(:,ii(1:nmeas)), &
         start=(/1,1/), count=(/nz,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( leo%ncid, leo%varid(ient_le_vmr), &
         omi_trc%le_trc_vmr_prof(:,ii(1:nmeas)), &
         start=(/1,1/), count=(/nz,nmeas/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! reset arrays after contents has been written:
    omi_trc%le2omi_trc_vcd_trop = -999.9
    omi_trc%le2omi_trc_vcd_prof = -999.9
    omi_trc%le_trc_vcd_prof     = -999.9
    omi_trc%le_trc_vmr_prof     = -999.9

    ! clear:
    deallocate( ii )

    ! ok
    status = 0

  end subroutine leo_PutOut


end module le_output_omi_trc



