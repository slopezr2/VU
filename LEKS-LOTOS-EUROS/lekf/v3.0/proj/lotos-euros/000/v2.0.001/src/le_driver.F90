!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Driver

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  LE_Message_Init, LE_Message_Done
  public    ::  LE_Timers_Init, LE_Timers_Done
  public    ::  LE_Model_Init, LE_Model_Done
  public    ::  LE_State_Init
  public    ::  LE_TimeStep_Init, LE_TimeStep_Run, LE_TimeStep_Done


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Driver'

  ! from minutes to seconds:  dt_sec = dt_min * min2sec
  real, parameter  ::  min2sec  = 60.0


  ! --- var --------------------------------------
  
  ! settings:
  type(TRcFile), target   ::  rcF
  ! pointer to settings, passed to routines for usage at run-time:
  type(TrcFile), pointer  ::  rcFP
  
  ! read OH fields ?
  logical   ::  read_OH

  ! timers:
  integer   ::  itim_setup, itim_setup_data, itim_setup_meteo, itim_setup_emis, itim_setup_bound
  integer   ::  itim_adjust, itim_chem, itim_vdif, itim_depos, itim_sedim
#ifdef with_m7
  integer   ::  itim_m7
#endif
  integer   ::  itim_wetdep, itim_emis, itim_adv, itim_rad
  integer   ::  itim_partupd


contains


  ! ===================================================================
  ! ===
  ! === arguments
  ! ===
  ! ===================================================================


  subroutine LE_Arguments( rcfile, status )

    use GO, only : goArgCount, goGetArg

    ! --- in/out ----------------------------------

    character(len=*), intent(out)     ::  rcfile
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Arguments'

    ! --- local -----------------------------------

    integer               ::  narg
    integer               ::  iarg
    character(len=1024)   ::  line

    ! --- begin -----------------------------------

    ! on root only, since some mpirun version do not parse
    ! all arguments to each executable:

    ! number of arguments:
    call goArgCount( narg, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( narg == 0 ) then
      write (gol,'("no arguments found ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! defaults:
    rcfile = 'None'

    ! loop over arguments:
    iarg = 0
    do
      ! next:
      iarg = iarg + 1
      ! get argument:
      call goGetArg( iarg, line, status )
      IF_NOTOK_RETURN(status=1)
      ! not filled yet ?
      if ( trim(rcfile) == 'None' ) then
        rcfile = trim(line)
      else
        write (gol,'("unsupported argument : ",a)') trim(line); call goErr
        TRACEBACK; status=1; return
      end if
      ! last one is processed now ?
      if ( iarg == narg ) exit
    end do

    ! ok
    status = 0

  end subroutine LE_Arguments


  ! ***


  subroutine LE_Print_Usage( status )

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- begin ----------------------------------

    ! display usage line:
    write (gol,'("Usage: lotos-euros.x <rcfile>")'); call goPr

    ! ok
    status = 0

  end subroutine LE_Print_Usage


  ! ====================================================================
  ! ===
  ! === messages
  ! ===
  ! ====================================================================


  subroutine LE_Message_Init( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Message_Init'

    ! --- external -------------------------------

#ifdef _OPENMP
    integer, external  ::  omp_get_thread_num
    integer, external  ::  omp_get_num_threads
#endif

    ! --- begin ----------------------------------

    write (gol,'(" " )'); call goPr
    write (gol,'("***********************************************************")'); call goPr
    write (gol,'("***********************************************************")'); call goPr
    write (gol,'("***                                                     ***")'); call goPr
    write (gol,'("***                     LOTOS-EUROS                     ***")'); call goPr
    write (gol,'("***                                                     ***")'); call goPr
    write (gol,'("*** (c) TNO, RIVM, PBL, KNMI                            ***")'); call goPr
    write (gol,'("***                                                     ***")'); call goPr
    write (gol,'("***********************************************************")'); call goPr
    write (gol,'("***********************************************************")'); call goPr
    write (gol,'(" " )'); call goPr
#ifdef _OPENMP
    write (gol,'("OpenMP: enabled                                   : yes")'); call goPr
    !$OMP PARALLEL
    write (gol,'("OpenMP: thread number / number of threads         : ",i4," / ",i4)') omp_get_thread_num(), omp_get_num_threads(); call goPr
    !$OMP END PARALLEL
#else
    write (gol,'("OpenMP: enabled                                   : no")'); call goPr
#endif
    write (gol,'(" " )'); call goPr

    ! ok
    status = 0

  end subroutine LE_Message_Init


  ! ***


  subroutine LE_Message_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Message_Done'

    ! --- begin ----------------------------------

    write (gol,'(" ")'); call goPr
    write (gol,'("***                          ")'); call goPr
    write (gol,'("*** end of simulation reached")'); call goPr
    write (gol,'("***                          ")'); call goPr
    write (gol,'(" ")'); call goPr

    ! ok
    status = 0

  end subroutine LE_Message_Done


  ! ====================================================================
  ! ===
  ! === timers
  ! ===
  ! ====================================================================


  subroutine LE_Timers_Init( status )

    use GO, only : GO_Timer_Init
  
    ! --- in/out -------------------------

    integer, intent(out)      ::  status
    
    ! --- const --------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Timers_Init'
    
    ! --- local --------------------------
    
    ! --- begin --------------------------

    ! setup timing tables:
    call GO_Timer_Init( status )
    IF_NOTOK_RETURN(status=1)
    
  end subroutine LE_Timers_Init


  ! ***
  
  
  subroutine LE_Timers_Done( rcfile, status, verbose, write_profile )

    use GO, only : pathsep, GetFileNameParts
    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : GO_Timer_Post, GO_Timer_Done
  
    ! --- in/out -------------------------

    character(len=*), intent(in)              ::  rcfile
    integer, intent(out)                      ::  status
    logical, intent(in), optional             ::  verbose
    logical, intent(in), optional             ::  write_profile
    
    ! --- const --------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Timers_Done'
    
    ! --- local --------------------------
    
    logical                   ::  write_file
    type(TrcFile)             ::  rcF
    character(len=512)        ::  outdir
    character(len=512)        ::  basename
    character(len=512)        ::  filename

    ! --- begin --------------------------
    
    ! write?
    write_file = .true.
    if ( present(write_profile) ) write_file = write_profile
    
    ! write profile?
    if ( write_file ) then

      ! open rcfile:
      call Init( rcF, rcfile, status )
      IF_NOTOK_RETURN(status=1)
      ! output directory:
      call ReadRc( rcF, 'le.output.outdir', outdir, status )
      IF_NOTOK_RETURN(status=1)
      ! close:
      call Done( rcF, status )
      IF_NOTOK_RETURN(status=1)

      ! extract basename without extension:
      call GetFileNameParts( rcfile, status, basename=basename, noext=.true. )
      IF_NOTOK_RETURN(status=1)
      ! name of profile output:
      write (filename,'(a,a,a,".prf")') trim(outdir), pathsep, trim(basename)

      ! stop timing, write timing profile, eventually print it too:
      call GO_Timer_Post( status, file=trim(filename), verbose=verbose )
      IF_NOTOK_RETURN(status=1)
      
    else

      ! stop timing, eventually print profile:
      call GO_Timer_Post( status, verbose=verbose )
      IF_NOTOK_RETURN(status=1)
    
    end if

    ! stop timing:
    call GO_Timer_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok:
    status = 0

  end subroutine LE_Timers_Done


  ! ====================================================================
  ! ===
  ! === model
  ! ===
  ! ====================================================================


  subroutine LE_Model_Init( t, dt, nhour, rcfile, status )

    use GO, only : TDate, NewDate, TIncrDate, IncrDate, Get, iTotal, operator(-)
    use GO, only : Init, ReadRc
    use GO, only : GO_Timer_Def

    use UDUnits         , only : UDUnits_Init, UDUnits_StrError
    use CF_Standard     , only : CF_Standard_Init
    use MDF             , only : MDF_Init

    use dims            , only : runF
    use dims            , only : nx, ny, nz
    use dims            , only : LE_Dims_Alloc
    use LE_LandUse      , only : LE_LandUse_Init
    use LE_Logging      , only : LE_Logging_Init
    use LE_Logging      , only : u_log
    use LE_Grid         , only : LE_Grid_Init
    use LE_Time         , only : LE_Time_Init
    use LE_Data         , only : LE_Data_Init
    use LE_Meteo_Data   , only : LE_Meteo_Data_Init
    use LE_Advec        , only : LE_Advec_Init
    use LE_VDiff        , only : LE_VDiff_Init
    use LE_Stability    , only : LE_Stability_Init
    use LE_DryDepos     , only : LE_DryDepos_Init
    use LE_Particle_Data, only : LE_Particle_Data_Init
    use LE_Sedim        , only : LE_Sedim_Init
    use LE_WetDepos     , only : LE_WetDepos_Init
    use LE_Chem         , only : LE_Chem_Init
    use LE_Emis         , only : LE_Emis_Init
    use Indices         , only : nspec
    use Indices         , only : LE_Indices_Init
    use LE_Chem_OH      , only : LE_Chem_OH_Init
    use LE_Solar        , only : solardec
    use LE_Bound        , only : LE_Bound_Init
    use LE_BiasCorr     , only : LE_Biascorr_Init 
    use Indices         , only : specname
#ifdef with_m7    
    use LE_M7           , only : LE_M7_Init
#endif
#ifdef with_radiation
    use LE_Radiation    , only:  LE_Radiation_Init
#endif
#ifdef with_vbs
    use LE_VBS          , only : LE_VBS_Init
#endif
#ifdef with_labeling
    use SA_Labeling, only : SA_Labeling_Init
#endif
!#ifdef _MPI
!    use MPI_Send        , only : MPI_Send_Dims
!#endif  

    ! --- in/out ---------------------------------

    type(TDate), intent(out)        ::  t
    type(TIncrDate), intent(out)    ::  dt
    integer, intent(out)            ::  nhour
#ifdef with_openda
    character(len=*), intent(in)    ::  rcfile
#else
    character(len=*), intent(out)   ::  rcfile
#endif
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Model_Init'

    ! --- local ----------------------------------

    character(len=1024)   ::  cf_standard_table
    integer               ::  timestep_min
    integer               ::  ip

    ! --- begin ----------------------------------

    ! define timers:
    call GO_Timer_Def( itim_setup   , 'time step setup', status )
    IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_setup_data , 'setup data' , status )
      IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_setup_meteo, 'setup meteo', status )
      IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_setup_bound, 'setup bound', status )
      IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_setup_emis , 'setup emis' , status )
      IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_adjust  , 'adjust', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_adv     , 'advection', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_chem    , 'chemistry', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_vdif    , 'vertical diffusion', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_depos   , 'dry deposition', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_partupd , 'particle update', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_sedim   , 'sedimentation', status )
    IF_NOTOK_RETURN(status=1)
#ifdef with_m7    
    call GO_Timer_Def( itim_m7      , 'M7', status )
    IF_NOTOK_RETURN(status=1)
#endif
    call GO_Timer_Def( itim_wetdep  , 'wet deposition', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_emis    , 'emission', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_rad     , 'radiation', status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** rcfile name
    !

#ifdef with_openda
#else
    ! extract argument, return name of rcfile:
    call LE_Arguments( rcfile, status )
    if (status/=0) then
      call LE_Print_Usage( status )
      status=1; return
    end if
#endif

    ! info ...
    write (gol,'("Name of LOTOS-EUROS rc file : ",a)') trim(rcfile); call goPr

    ! open rcfile to read other settings:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! assign pointer:
    rcFP => rcF


    !
    ! ** HDF/NetCDF interface
    !

    call MDF_Init( status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** UDUnits interface
    !

    ! init:
    call UDUnits_Init( status )
    if ( status /= 0 ) then
      gol = trim(UDUnits_StrError(status)); call goErr
      TRACEBACK; status=1; return
    end if


    !
    ! ** CF conventions
    !

    ! read name of xml table:
    call ReadRc( rcF, 'cf.standard.table', cf_standard_table, status )
    IF_NOTOK_RETURN(status=1)
    ! init table:
    call CF_Standard_Init( cf_standard_table, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** logging
    !

    ! open log files etc:
    call LE_Logging_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! write details to the log file
    write (u_log,*) 'input specifications: see "'//trim(rcfile)//'" in this directory'
    write (u_log,*)


    !
    ! ** grid definition
    !

    ! grid dimension and position:
    call LE_Grid_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

!#ifdef _MPI  
!    call MPI_Send_Dims( 'OPS', status )    
!    IF_NOTOK_RETURN(status=1)
!#endif


    !
    ! ** data
    !

    ! init model data
    call LE_Data_Init( rcFP, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** init flags
    !

    ! first time step is true
    runF%first = .true.


    !
    ! ** bias correction
    !

    ! read bias correction parameters,
    ! should be done before accumulated tracer are defined:
    call LE_BiasCorr_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    

    !
    ! ** set tracer indices
    !

    ! initialize indices variables:
    call LE_Indices_Init( status )
    IF_NOTOK_RETURN(status=1)
    

    !
    ! ** allocate
    !

    ! NOW we know everything to initialise arrays
    call LE_Dims_Alloc()


    !
    ! ** time definition
    !

    ! time range etc:
    call LE_Time_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill start time:
    t = NewDate( year=runF%yy_s, month=runF%mm_s, day=runF%dd_s, hour=runF%hh_s )

    !
    ! ** landuse
    !

    call LE_LandUse_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** meteo data
    !

    ! init some meteo data that is not part of LE_Data yet:
    call LE_Meteo_Data_Init( status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ** restart
    !

    ! read from restart file ?
    call ReadRc( rcF, 'le.restart', runF%restart, status )
    IF_NOTOK_RETURN(status=1)

    ! read from restart file ?
    if ( runF%restart ) then
       ! where to read from:
       call ReadRc( rcF, 'le.restart.path', runF%restart_path, status )
       IF_NOTOK_RETURN(status=1)
       ! file name description:
       call ReadRc( rcF, 'le.restart.key', runF%restart_key, status )
       IF_NOTOK_RETURN(status=1)
    else
       ! dummy ...
       runF%restart_path = 'none'
       runF%restart_key  = 'none'
    end if

    ! write restart files ?
    call ReadRc( rcF, 'le.restart.save', runF%restart_save, status )
    IF_NOTOK_RETURN(status=1)

    ! write restart files ?
    if ( runF%restart_save ) then
       ! when to write:
       call ReadRc( rcF, 'le.restart.save.dhour', runF%restart_save_dhour, status, default=24.0 )
       IF_NOTOK_RETURN(status=1)
       ! where to write:
       call ReadRc( rcF, 'le.restart.save.path', runF%restart_save_path, status )
       IF_NOTOK_RETURN(status=1)
       ! file name description:
       call ReadRc( rcF, 'le.restart.save.key', runF%restart_save_key, status )
       IF_NOTOK_RETURN(status=1)
    else
       runF%restart_save_dhour = 9999.9
       runF%restart_save_path  = 'none'
       runF%restart_save_key   = 'none'
    end if

    !
    ! ** advection/diffusion
    !

    ! setup advection
    call LE_Advec_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! setup vdiff
    call LE_VDiff_Init( status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ** timing
    !

    ! main time step (between output times) in minutes:
    call ReadRc( rcF, 'timestep.output', timestep_min, status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    dt = IncrDate( min=timestep_min )

    ! hours since start of the day:
    call Get( t, hour=nhour )
    if ( runF%forecast_mode ) then
       nhour = iTotal( t - runF%t_base, 'hour' )
       write (gol,'(" forecast mode, set hour offset nhour = ",i3)') nhour ; call goPr
    end if

#ifdef with_labeling
  call SA_Labeling_Init( rcfile, 'le.output', status)
  IF_NOTOK_RETURN(status=1)
#endif
    
    !
    ! ** chemistry
    !

    call LE_Chem_Init( rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
#ifdef with_m7
    !
    ! ** M7 aerosol scheme
    !
    
    call LE_M7_Init( status )
    IF_NOTOK_RETURN(status=1)
#endif
    
#ifdef with_vbs
    !
    ! ** VBS gas-to-aerosol condensation
    !
    
    call LE_VBS_Init( rcfile, status )
    IF_NOTOK_RETURN(status=1)
#endif

    !
    ! ** stability
    !

    call LE_Stability_Init( status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ** deposition
    !

    ! before sedim init!
    call LE_Particle_Data_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! before depos init!
    call LE_Sedim_Init( status )
    IF_NOTOK_RETURN(status=1)

    call LE_DryDepos_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

    call LE_WetDepos_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_radiation
    !
    !  ** radiation properties
    !

    call LE_Radiation_Init(rcF, status)
    IF_NOTOK_RETURN(status=1)
#endif

    !
    ! ** emissions
    !

    ! read the emissions
    ! setup emissions:
    call LE_Emis_Init( rcF, t, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ** OH fields
    !    Used for sulphur- or methane-only runs.
    !
    
    ! read OH fields ?
    call ReadRc( rcF, 'le.chem.OH', read_OH, status )
    IF_NOTOK_RETURN(status=1)

    ! enabled ?
    if ( read_OH ) then
      call LE_Chem_OH_Init( rcfile, status )
      IF_NOTOK_RETURN(status=1)
    end if

    !
    ! ** sun
    !

    ! solar declination (per day)
    call solardec( t%month, t%day )

    !
    ! ** boundary conditions
    !

    ! setup boundary conditions:
    call LE_Bound_Init( rcfile, t, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! **
    !

    ! ok
    status = 0

  end subroutine LE_Model_Init


  ! ***


  subroutine LE_Model_Done( status )

    use GO              , only : Done
    use UDUnits         , only : UDUnits_Done, UDUnits_StrError
    use CF_Standard     , only : CF_Standard_Done
    use MDF             , only : MDF_Done
    use dims            , only : runF
    use LE_BiasCorr     , only : LE_BiasCorr_Done
    use Indices         , only : LE_Indices_Done
    use LE_LandUse      , only : LE_LandUse_Done
    use LE_Advec        , only : LE_Advec_Done
    use LE_VDiff        , only : LE_VDiff_Done
    use LE_Stability    , only : LE_Stability_Done
    use LE_Particle_Data, only : LE_Particle_Data_Done
    use LE_Sedim        , only : LE_Sedim_Done
    use LE_DryDepos     , only : LE_DryDepos_Done
    use LE_WetDepos     , only : LE_WetDepos_Done
    use LE_Chem         , only : LE_Chem_Done
    use LE_Emis         , only : LE_Emis_Done
    use LE_Bound        , only : LE_Bound_Done
    use LE_Grid         , only : LE_Grid_Done
    use LE_Time         , only : LE_Time_Done
    use LE_Meteo_Data   , only : LE_Meteo_Data_Done
    use LE_Data         , only : LE_Data_Done
    use LE_Logging      , only : LE_Logging_Done
    use LE_Logging      , only : u_log
    use LE_Chem_OH      , only : LE_Chem_OH_Done
#ifdef with_m7    
    use LE_M7           , only : LE_M7_Done
#endif
#ifdef with_radiation
    use LE_Radiation    , only:  LE_Radiation_Done
#endif
#ifdef with_vbs
    use LE_VBS          , only : LE_VBS_Done
#endif
#ifdef with_labeling
    use SA_Labeling     , only : SA_Labeling_Done
#endif        

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Model_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with advection
    call LE_Advec_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with vdiff
    call LE_VDiff_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with OH fields:
    if( read_OH ) then
      call LE_Chem_OH_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! done with emissions:
    call LE_Emis_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with boundary conditions:
    call LE_Bound_Done( status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_vbs
    ! done with VBS scheme:
    call LE_VBS_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_radiation
    ! done with radiation:
    call LE_Radiation_Done( status)
    IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_m7
    ! done with aerosols:
    call LE_M7_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! done with chemistry:
    call LE_Chem_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with deposition:
    call LE_DryDepos_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with sedimentation:
    call LE_Sedim_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with Partsize:
    call LE_Particle_Data_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with wet deposition:
    call LE_WetDepos_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with stability:
    call LE_Stability_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with meteo data:
    call LE_Meteo_Data_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with landuse:
    call LE_LandUse_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with time defintion:
    call LE_Time_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with grid defintion structure:
    call LE_Grid_Done( status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
    call SA_Labeling_Done( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! done with CF standard:
    call CF_Standard_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with UDUnits interface:
    call UDUnits_Done( status )
    if ( status /= 0 ) then
      gol = trim(UDUnits_StrError(status)); call goErr
      TRACEBACK; status=1; return
    end if

    ! done with tracer indices:
    call LE_Indices_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with bias correction parametes:
    call LE_BiasCorr_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with data:
    call LE_Data_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with HDF/NetCDF interface:
    call MDF_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with settings:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! done with logging:
    call LE_Logging_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! write dummy file to indicate proper end:
    call LE_Write_OkFile( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Model_Done


  ! Write dummy file 'le.ok'.
  ! Existence of this file is used by the scripts to check
  ! if a run ended properly.
  ! Checking exit status would be better, but this does
  ! not trap 'stop' statements and other obscure endings.

  subroutine LE_Write_OkFile( status )

    use GO, only : goGetFU

    ! --- in/out ----------------------------------

    integer, intent(out)          ::  status

    ! --- const ------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Write_OkFile'

    ! --- local ----------------------------------

    integer           ::  fu

    ! --- begin -----------------------------------

    ! get free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open file:
    open( unit=fu, file='le.ok', form='formatted', status='unknown', iostat=status )
    if ( status/=0 ) then
      write (gol,'("from opening okfile")'); call goErr
    else
      ! write happy message:
      write (fu,'("Program terminated normally")',iostat=status)
      if ( status/=0 ) then
        write (gol,'("from writing to okfile")'); call goErr
      else
        ! close:
        close( fu, iostat=status )
        if ( status/=0 ) then
          write (gol,'("from closing okfile")'); call goErr
        end if
      end if
    end if

    ! ok
    status = 0

  end subroutine LE_Write_OkFile


  ! ====================================================================
  ! ===
  ! === state
  ! ===
  ! ====================================================================


  subroutine LE_State_Init( c, t, status )

    use GO        , only : TDate
    use dims      , only : nx, ny, nz, nspec
    use LE_Logging , only : ident1
    use LE_Bound  , only : LE_Bound_Initial
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ---------------------------------

    real, intent(out)               ::  c(nx,ny,nz,nspec)
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_State_Init'

    ! --- begin ----------------------------------

    ! setup initial field from coarse model if possible,
    ! otherwise, interpolate between east/west boundary conditions:
    call LE_Bound_Initial( c, t, status )
    IF_NOTOK_RETURN(status=1)

    !! info ...
    !print *, ident1,'<initialised concentrations>'

    ! enable data:
    call LE_Data_Enable( 'occ', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_State_Init


  ! ====================================================================
  ! ===
  ! === timestep
  ! ===
  ! ====================================================================


  subroutine LE_TimeStep_Init( t, dt, nhour, status )

    use GO              , only : TDate, TIncrDate, DayNumber, MidNight, operator(+), wrtgol
    use GO              , only : GO_Timer_Start, GO_Timer_End
    use LE_Logging      , only : ident1, ident2
    use dims            , only : runF
    use dims            , only : coszen
    use LE_Data         , only : LE_Data_GetPointer
    use LE_Data         , only : LE_Data_Setup
    use LE_Solar        , only : solardec
    use LE_Solar        , only : zangle
    use LE_Meteo_Data   , only : LE_Meteo_Data_Setup
    use LE_Restart      , only : LE_Restart_Restore_Data
    use Indices         , only : i_no, i_no2, i_o3
    use LE_Chem         , only : calc_r1r3
    use LE_LandUse      , only : LE_Landuse_Update
    use LE_Stability    , only : monin_ustar
    use LE_Emis         , only : LE_Emis_Setup
    use LE_Bound        , only : LE_Bound_Setup

    ! --- in/out ---------------------------------

    type(TDate), intent(in)         ::  t
    type(TIncrDate), intent(in)     ::  dt
    integer, intent(in)             ::  nhour
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_TimeStep_Init'

    ! --- local ----------------------------------
    
    ! meteo data:
    real, pointer        ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  occ  (:,:,:)   ! (lon,lat,jz)

    ! --- begin ----------------------------------
    
    ! start timing:
    call GO_Timer_Start(itim_setup, status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    call wrtgol('LE: prepare timestep ',t,' to ',t+dt); call goPr

    ! now also calculate the solar declination for the new day
    if ( MidNight(t) ) then
      !print *, ident2,'<solardec>'
      call solardec( t%month, t%day )
      !print *, ident2,'  ok'
    end if
    
    ! ~
    
    ! start timing:
    call GO_Timer_Start( itim_setup_data, status )
    IF_NOTOK_RETURN(status=1)
    ! data for present time step:
    call LE_Data_Setup( runF%t_base, t, t+dt, status )
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End( itim_setup_data, status )
    IF_NOTOK_RETURN(status=1)

    ! start timing:
    call GO_Timer_Start( itim_setup_meteo, status )
    IF_NOTOK_RETURN(status=1)
    ! meteo data for present time step:
    call LE_Meteo_Data_Setup( t, dt, status )
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End( itim_setup_meteo, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ~

    ! saved parameters
    if ( runF%first ) then
      ! read saved model data ?
      if ( runF%restart ) then
        ! Read saved parameters:
        !  o old volume is used in 'adjust_c';
        !  o old exposure class is used in 'exposure' called from 'monin_ustar' .
        call LE_Restart_Restore_Data( t, trim(runF%restart_path), trim(runF%restart_key), status )
        IF_NOTOK_RETURN(status=1)
      end if  ! restore from restart file ?
    end if  ! first time step ?

    ! cosine of solar angle
    call zangle(t%hour, 30.0)

    if ( all( (/i_no2,i_no,i_o3/) > 0 ) ) then
      ! reaction ratio for photostationary equilibrium
      ! use surface temperature array, since it involves
      ! sub-grid scale concentrations close to the surface
      ! point to meteo data:
      call LE_Data_GetPointer( 'tsurf', tsurf, status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_GetPointer( 'occ', occ, status, check_units ='1' )
      IF_NOTOK_RETURN(status=1)
      ! no 3D yet, use total cloud cover (overhead cloud cover at bottom layer):
      call calc_r1r3( occ(:,:,1), coszen, tsurf, status )
      IF_NOTOK_RETURN(status=1)
    endif
    
    ! update surface roughness, on sea it depends on winspeed:
    call LE_Landuse_Update( t, runF%first .and. (.not. runF%restart), status )
    IF_NOTOK_RETURN(status=1)

    ! compute stability stuff;
    ! if first time that this routine is called and old exposure class
    ! has not been read from restart file, then set lfirst to .true. :
    call monin_ustar( runF%first .and. (.not. runF%restart), status )
    IF_NOTOK_RETURN(status=1)

    call GO_Timer_Start( itim_setup_emis, status )
    IF_NOTOK_RETURN(status=1)
    call LE_Emis_Setup( t, t+dt, status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_End( itim_setup_emis, status )
    IF_NOTOK_RETURN(status=1)

    ! fill boundary conditions:
    call GO_Timer_Start( itim_setup_bound, status )
    IF_NOTOK_RETURN(status=1)
    call LE_Bound_Setup( t, t+dt, nhour, status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_End( itim_setup_bound, status )
    IF_NOTOK_RETURN(status=1)
    
    !! info:
    !print *, ident1,'<end preparation input for model step>'

    ! end timing:
    call GO_Timer_End(itim_setup, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_TimeStep_Init


  ! ***


  subroutine LE_TimeStep_Run( t, dt, c, cg, aerh2o, bud, update_bud, status )

    use GO              , only : TDate, TIncrDate, operator(+), wrtgol, rTotal
    use GO              , only : GO_Timer_Start, GO_Timer_End
    use LE_Logging      , only : ident1, ident2
    use dims            , only : runF
    use dims            , only : nx, ny, nz, nspec
    use dims            , only : OPS_mass
    use LE_Advec        , only : LE_Advec_Get_NStep
    use LE_Advec        , only : LE_Advec_Apply
    use LE_VDiff        , only : LE_VDiff_Apply
    use LE_DryDepos     , only : LE_DryDepos_Setup_vd
    use LE_DryDepos     , only : LE_DryDepos_Apply
    use LE_DryDepos     , only : mix2ground
    use LE_Sedim        , only : LE_Sedim_Apply
#ifdef with_m7    
    use LE_M7           , only : LE_M7_Apply
#endif
    use LE_Particle_Data, only : LE_Particle_Data_Update
    use LE_Chem         , only : LE_Chem_Step
    use LE_Chem_OH      , only : LE_Chem_OH_Setup
    use LE_WetDepos     , only : LE_WetDepos_Apply
    use LE_Adjust       , only : LE_Adjust_Apply
    use LE_Emis         , only : LE_Emis_Add
#ifdef with_radiation
    use LE_Radiation    , only: LE_Radiation_Calc
#endif
    use LE_Logging        , only : PrintDebug
    use LE_Budget         , only : T_Budget
    use LE_Budget         , only : Budget_Reset_Hour, Budget_Update
    use LE_Budget_DryDepos, only : ex_drydepo
    use LE_Budget_WetDepos, only : ex_wetdepo
    use LE_BiasCorr       , only : LE_BiasCorr_Fill
#ifdef with_labeling
    use SA_Labeling,        only : SA_Synchronize_Conc, SA_frac
    use SA_Labeling,        only : SA_DryDepos_Reset, SA_WetDepos_Reset
#endif    
!#ifdef _MPI    
!    use MPI_Get           , only :  MPI_Get_OPS_emis
!#endif    
    
    ! --- in/out ---------------------------------

    type(TDate), intent(in)         ::  t
    type(TIncrDate), intent(in)     ::  dt
    real, intent(inout)             ::  c(nx,ny,nz,nspec)
    real, intent(inout)             ::  cg(nx,ny   ,nspec)
    real, intent(inout)             ::  aerh2o(nx,ny,nz)
    type(T_Budget), intent(inout)   ::  bud
    logical, intent(in)             ::  update_bud
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_TimeStep_Run'

    ! --- local ----------------------------------

    integer     ::  it
    real        ::  timestep_min  ! operator splitting half-step in minutes
    real        ::  halfstep_min  ! operator splitting half-step in minutes
    integer     ::  nadvstep
    integer     ::  nfullstep
    logical     ::  update_vd_only_if_depends_on_c
    integer     ::  i_wetdepo_step  ! step number for wet deposition

    ! --- begin ----------------------------------

    ! info ...
    call wrtgol( '    model step over ', t, ' - ', t+dt ); call goPr

    ! reset hourly budgets:
    if ( update_bud ) then
      call Budget_Reset_Hour( bud, status )
      IF_NOTOK_RETURN(status=1)
    end if

    !!>>> enable for compatibility with v1.9.000 :
    !! deposition velocities: concentration needed,
    !! since dependend on nh3/so2 ratio regime;
    !! timing included:
    !call LE_DryDepos_Setup_vd( c, bud%drydepos, .false., t, status )
    !IF_NOTOK_RETURN(status=1)
    !!<<<
    ! first call to setup vd should update everything:
    update_vd_only_if_depends_on_c = .false.

    ! adjust concentrations due to change in layer heigths ;
    ! uses the old volume from the previous time step;
    ! no idea why this is not called at the start of the run,
    ! the old volume is set to the current volume in this case
    if ( (.not. runF%first) .or. (runF%first .and. runF%restart) ) then
       call GO_Timer_Start( itim_adjust, status )
       IF_NOTOK_RETURN(status=1)
       call LE_Adjust_Apply( c, status )
       IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
       call SA_Synchronize_Conc(c, status)
       IF_NOTOK_RETURN(status=1)
#endif              
       
       call GO_Timer_End( itim_adjust, status )
       IF_NOTOK_RETURN(status=1)
       call printdebug( c, 'adjust', status )
       IF_NOTOK_RETURN(status=1)
    endif

#ifdef with_labeling
    ! reset hourly dry and wet deposition budgets
    call SA_DryDepos_Reset( status )
    IF_NOTOK_RETURN(status=1)
    call SA_WetDepos_Reset( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! read the OH concentrations for this hour for sulphur- or methane-only run
    if ( read_OH ) then
      call LE_Chem_OH_Setup( t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! interval length:
    timestep_min = rTotal( dt, 'min' )

    ! minium number of  steps required for advection 
    ! within current interval:
    call LE_Advec_Get_NStep( timestep_min, nadvstep, status )
    IF_NOTOK_RETURN(status=1)
    
    ! the two advection half-steps are performed as a full-step,
    ! therefore the number of operator splitting steps is the
    ! same as the minimum number of advection steps:
    nfullstep = nadvstep
    
    ! set time step of half-step of operator splitting sequence:
    halfstep_min = 0.5 * timestep_min/nfullstep

    ! info ...
    write (gol,'("     perform ",i2," operator sequences over ",f7.2," minutes each")') nfullstep, 2*halfstep_min; call goPr
    
    ! setup step number of wetdeposition
    i_wetdepo_step = 1
    
    ! debug ...
    call printdebug( c, 'begin', status )
    IF_NOTOK_RETURN(status=1)
        
    ! perform full operator splitting steps:
    do it = 1, nfullstep

#ifdef without_physics
      ! no physics ...
#else
    
      ! ~ chemistry (1/2)

      ! chemistry;
      ! only first halfstep, since the chemistry is applied
      ! at the end too with twice the timestep:
      ! aerosol dynamics is treated right after chemistry, to account for 
      ! h2so4 nucleation/condensation and other aerosol processes
      if (it==1) then

        ! start timing:
        call GO_Timer_Start( itim_chem, status )
        IF_NOTOK_RETURN(status=1)
        ! run:
        call LE_Chem_Step( c, aerh2o, t, halfstep_min, runF%restart, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_chem, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug( c, 'chem', status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_m7
        ! start timing:
        call GO_Timer_Start( itim_m7, status )
        IF_NOTOK_RETURN(status=1)
        ! unit conversions and call to m7
        call LE_M7_Apply( c, halfstep_min, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_m7, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug( c, 'm7', status )
        IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_labeling
        ! labeling:
        call SA_Synchronize_Conc( c, status )
        IF_NOTOK_RETURN(status=1)
#endif        
      endif
      
      ! ~ vertical diffusion (1/2)

      ! start timing:
      call GO_Timer_Start( itim_vdif, status )
      IF_NOTOK_RETURN(status=1)

      ! vertical diffusion
      call LE_VDiff_Apply( c, halfstep_min, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif        

      ! end timing:
      call GO_Timer_End( itim_vdif, status )
      IF_NOTOK_RETURN(status=1)
      
      ! info ...
      call printdebug( c, 'vdif', status )
      IF_NOTOK_RETURN(status=1)

      
      ! ~ sedimentation (1/2)
      
      ! start timing:
      call GO_Timer_Start(itim_sedim, status )
      IF_NOTOK_RETURN(status=1)

      ! compute new M7 wet radius, particle size, slip correction, etc:
      call LE_Particle_Data_Update( c, status )
      IF_NOTOK_RETURN(status=1)

      ! apply sedimentaiton:
      call LE_Sedim_Apply( c, bud%drydepos,update_bud, halfstep_min*min2sec, status )
      IF_NOTOK_RETURN(status=1)
      
#ifdef with_labeling
      ! update labels:
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif        

      ! end timing:
      call GO_Timer_End(itim_sedim, status )
      IF_NOTOK_RETURN(status=1)
      
      ! debug ...
      call printdebug( c, 'sedim', status )
      IF_NOTOK_RETURN(status=1)
      
      ! ~ dry deposition (1/2)

      ! start timing:
      call GO_Timer_Start( itim_partupd, status )
      IF_NOTOK_RETURN(status=1)
      ! compute new M7 wet radius, particle size, slip correction, etc:
      call LE_Particle_Data_Update( c, status )
      IF_NOTOK_RETURN(status=1)
      ! end timing:
      call GO_Timer_End( itim_partupd, status )
      IF_NOTOK_RETURN(status=1)

      !>>> skip for compatibility with v1.9.000 :
      ! deposition velocities: concentration needed,
      ! since dependend on nh3/so2 ratio regime;
      ! budgets needed for NH3 compensation points;
      ! timing included in routine:
      call LE_DryDepos_Setup_vd( c, bud%drydepos, update_vd_only_if_depends_on_c, t, status )
      IF_NOTOK_RETURN(status=1)
      ! from now only update vd if depends on concentrations:
      update_vd_only_if_depends_on_c = .true.
      !<<<

      ! start timing:
      call GO_Timer_Start(itim_depos, status )
      IF_NOTOK_RETURN(status=1)

      ! apply dry deposition
      call LE_DryDepos_Apply( c, halfstep_min, bud%drydepos,update_bud, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
      ! update labels:
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif       

      ! end timing:
      call GO_Timer_End(itim_depos, status )
      IF_NOTOK_RETURN(status=1)
      
      ! info ...
      call printdebug( c, 'depos', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ wet deposition (1/2)

      ! start timing:
      call GO_Timer_Start(itim_wetdep, status )
      IF_NOTOK_RETURN(status=1)

      ! wet deposition
      call LE_WetDepos_Apply( c, bud%wetdepos%ex_hour(:,:,:,:,ex_wetdepo),  update_bud, &
                          halfstep_min, i_wetdepo_step, status )
      IF_NOTOK_RETURN(status=1)
      ! add one step number for wet depositions steps
      i_wetdepo_step = i_wetdepo_step + 1
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif  
      ! end timing:
      call GO_Timer_End(itim_wetdep, status )
      IF_NOTOK_RETURN(status=1)
      
      ! info ...
      call printdebug( c, 'wetdep', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ emissions (1/2)

      ! start timing:
      call GO_Timer_Start(itim_emis, status )
      IF_NOTOK_RETURN(status=1)

      ! add emissions to concentration arrays given timestep in minutes:
      call LE_Emis_Add( c, halfstep_min, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif       

      ! end timing:
      call GO_Timer_End(itim_emis, status )
      IF_NOTOK_RETURN(status=1)
      
      ! info ...
      call printdebug( c, 'gas_emis', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ advection (1/2 and 2/2 at once)

      ! start timing:
      call GO_Timer_Start(itim_adv, status )
      IF_NOTOK_RETURN(status=1)

      ! advection:
      call LE_Advec_Apply( c, 2.0*halfstep_min, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif 
      ! end timing:
      call GO_Timer_End(itim_adv, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      call printdebug(c,'advec', status )
      IF_NOTOK_RETURN(status=1)

#endif

      ! set first time step flag false
      runF%first = .false.

#ifdef without_physics
      ! no physics ...
#else

      ! ~ emissions (2/2)

      ! start timing:
      call GO_Timer_Start(itim_emis, status )
      IF_NOTOK_RETURN(status=1)

      ! add emissions to concentration arrays given timestep in minutes:
      call LE_Emis_Add( c, halfstep_min, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif 
      ! end timing:
      call GO_Timer_End(itim_emis, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      call printdebug(c,'gas_emis', status )
      IF_NOTOK_RETURN(status=1)
      

      ! ~ wet deposition (2/2)

      ! start timing:
      call GO_Timer_Start(itim_wetdep, status )
      IF_NOTOK_RETURN(status=1)
      
      ! wet deposition:
      call LE_WetDepos_Apply( c, bud%wetdepos%ex_hour(:,:,:,:,ex_wetdepo),  update_bud, &
                                   halfstep_min, i_wetdepo_step, status )
      IF_NOTOK_RETURN(status=1)
      ! add one step number for wet depositions steps
      i_wetdepo_step = i_wetdepo_step + 1
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif
      ! end timing:
      call GO_Timer_End(itim_wetdep, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      call printdebug(c,'wetdep', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ dry deposition (2/2)

      ! start timing:
      call GO_Timer_Start( itim_partupd, status )
      IF_NOTOK_RETURN(status=1)
      ! compute new M7 wet radius, particle size, slip correction, etc:
      call LE_Particle_Data_Update( c, status )
      IF_NOTOK_RETURN(status=1)
      ! end timing:
      call GO_Timer_End( itim_partupd, status )
      IF_NOTOK_RETURN(status=1)

      !>>> skip for compatibility with v1.9.000 :
      ! deposition velocities: concentration needed, 
      ! since dependend on nh3/so2 ratio regime and particle sizes;
      ! budgets needed for NH3 compensation point
      ! timing included in routine:
      call LE_DryDepos_Setup_vd( c, bud%drydepos, update_vd_only_if_depends_on_c, t, status )
      IF_NOTOK_RETURN(status=1)
      !<<<

      ! start timing:
      call GO_Timer_Start(itim_depos, status )
      IF_NOTOK_RETURN(status=1)

      ! apply dry deposition:
      call LE_DryDepos_Apply( c, halfstep_min, bud%drydepos, update_bud, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
      ! update labels:
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif

      ! end timing:
      call GO_Timer_End(itim_depos, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      call printdebug(c,'depos', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ sedimentation (2/2)

      ! start timing:
      call GO_Timer_Start(itim_sedim, status )
      IF_NOTOK_RETURN(status=1)

      ! compute new M7 wet radius, particle size, slip correction, etc:
      call LE_Particle_Data_Update( c, status )
      IF_NOTOK_RETURN(status=1)

      ! sedimentation:
      call LE_Sedim_Apply( c, bud%drydepos, update_bud, halfstep_min*min2sec, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
      ! update labels:
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif

      ! end timing:
      call GO_Timer_End(itim_sedim, status )
      IF_NOTOK_RETURN(status=1)

      ! debug ...
      call printdebug(c,'sedim', status )
      IF_NOTOK_RETURN(status=1)
      
      ! ~ vertical diffusion (2/2)

      ! start timing:
      call GO_Timer_Start(itim_vdif, status )
      IF_NOTOK_RETURN(status=1)

      ! vertical diffusion
      call LE_VDiff_Apply( c, halfstep_min, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
#endif        

      ! end timing:
      call GO_Timer_End(itim_vdif, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      call printdebug(c,'vdif', status )
      IF_NOTOK_RETURN(status=1)

      ! ~ chemistry (2/2) combined with next (1/2) if possible
      
      ! last halfstep ?
      if ( it == nfullstep ) then

        ! start timing:
        call GO_Timer_Start( itim_chem, status )
        IF_NOTOK_RETURN(status=1)
        ! apply chemistry over half time step:
        call LE_Chem_Step( c, aerh2o, t,     halfstep_min, runF%restart, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_chem, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug(c,'chem', status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_m7
        ! start timing:
        call GO_Timer_Start( itim_m7, status )
        IF_NOTOK_RETURN(status=1)
        ! apply:
        call LE_M7_Apply( c, halfstep_min, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_m7, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug(c,'m7', status )
        IF_NOTOK_RETURN(status=1)
#endif

      else

        ! start timing:
        call GO_Timer_Start( itim_chem, status )
        IF_NOTOK_RETURN(status=1)
        ! apply chemistry over double half time step:
        call LE_Chem_Step( c, aerh2o, t, 2.0*halfstep_min, runF%restart, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_chem, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug(c,'chem', status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_m7
        ! start timing:
        call GO_Timer_Start( itim_m7, status )
        IF_NOTOK_RETURN(status=1)
        ! apply:
        call LE_M7_Apply( c, 2*halfstep_min, status )
        IF_NOTOK_RETURN(status=1)
        ! end timing:
        call GO_Timer_End( itim_m7, status )
        IF_NOTOK_RETURN(status=1)

        ! debug ...
        call printdebug(c,'m7', status )
        IF_NOTOK_RETURN(status=1)
#endif
      end if

#ifdef with_labeling
      ! start timing:
      call GO_Timer_Start( itim_chem, status )
      IF_NOTOK_RETURN(status=1)
      ! update labels:
      call SA_Synchronize_Conc(c, status)
      IF_NOTOK_RETURN(status=1)
      ! end timing:
      call GO_Timer_End( itim_chem, status )
      IF_NOTOK_RETURN(status=1)
#endif
      
      ! ~

#endif

    enddo ! loop within one hour
    
    ! compute new M7 wet radius, particle size, slip correction, etc:
    call LE_Particle_Data_Update( c, status )
    IF_NOTOK_RETURN(status=1)

    ! update budgets ?
    if ( update_bud ) then
      ! update budgets:
      call Budget_Update( bud, t+dt, status )
      IF_NOTOK_RETURN(status=1)
    end if
     
    ! re-compute deposition velocities if these depend on concentrations;
    ! budgets needed for NH3 compensation point:
    call LE_DryDepos_Setup_vd( c, bud%drydepos, .true., t, status )
    IF_NOTOK_RETURN(status=1)
    ! surface concentrations:
    call mix2ground( c, cg, status )
    IF_NOTOK_RETURN(status=1)

    ! fill bias corrected tracers:
    call LE_BiasCorr_Fill( c, cg, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_radiation
    ! start timing:
    call GO_Timer_Start( itim_rad, status )
    IF_NOTOK_RETURN(status=1)
    ! radiation properties:
    call LE_Radiation_Calc( c, status)
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End( itim_rad, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! print message
    write (gol,'(a,"<end of model step>")') ident1; call goPr

    !! info for OpenMP tests:
    !do is = 1, nspec
    !  write (*,'("XXX ",a10, 3f16.8)') trim(specname(is)), &
    !        minval(c(:,:,:,is)), maxval(c(:,:,:,is)), &
    !        sum(c(:,:,:,is))/size(c(:,:,:,is))
    !end do

    ! ok
    status = 0

  end subroutine LE_TimeStep_Run


  ! ***


  subroutine LE_TimeStep_Done( t, dt, nhour, status )

    use GO        , only : TDate, TIncrDate, MidNight, operator(+)
    use dims      , only : runF
    use LE_Logging, only : ident1, ident2

    ! --- in/out ---------------------------------

    type(TDate), intent(inout)      ::  t
    type(TIncrDate), intent(in)     ::  dt
    integer, intent(inout)          ::  nhour
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_TimeStep_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! info ...
    write (gol,'(a,"<postprocessing  model step>")') ident2; call goPr

    ! update actual model time
    t = t + dt
    nhour = nhour+1
    ! end of current day ?
    if ( MidNight(t) ) then
      ! reset hour counter if no forecast run:
      if ( .not. runF%forecast_mode ) nhour = 0
    end if

    ! info ...
    write (gol,'(a,"<end postprocessing  model step>")') ident2; call goPr

    ! ok
    status = 0

  end subroutine LE_TimeStep_Done


end module LE_Driver
