!#################################################################
!
!  KALMAN - Kalman Filter around LOTOS-EUROS
!
! Source tree:
!
!  lotos-euros-kf
!
!    ! modules per measurement type, selected using #ifdef's
!    meas_maori         ! ground observations using MAORI interface
!    meas_omi_no2       ! OMI NO2 tropospheric columns
!    meas_modis         ! MODIS AOD columns
!
!      ! shared routines for analysis etc:
!      meas
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; call goc%Abort(status); end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; call goc%Abort(status); end if

!
#include "lekf.inc"
!
!###############################################################################

program LOTOS_EUROS_KF

  use GO, only : GO_Init, GO_Done
  use GO, only : gol, goPr, goErr
  use GO, only : goExit
  use GO, only : goc
  use GO, only : lgr

  use GO, only : TrcFile
  use GO, only : TDate, TIncrDate
  use GO, only : IncrDate
  use GO, only : operator(==), operator(+), operator(-), Midnight, rTotal
  use GO, only : wrtgol, Precisely
  
  use GO, only : GO_Timer_Def, GO_Timer_Start, GO_Timer_End, GO_Timer_Switch

  use dims , only : nx, ny, nz, nspec
  use dims , only : runF

  use LE_Config , only : LE_Config_Init, LE_Config_Done

  use LE_Driver , only : LE_Timers_Init, LE_Timers_Done
  use LE_Driver , only : LE_Message_Init, LE_Message_Done
  use LE_Driver , only : LE_Model_Init, LE_Model_Done
  use LE_Driver , only : LE_TimeStep_Init, LE_TimeStep_Done

  use LEKF_state, only : LEKF_State_Init, LEKF_State_Done
  use LEKF_state, only : nmodes
  use LEKF_state, only : xb, x, sigma
  use LEKF_state, only : Ens
!  use LEKF_state, only : imodes
  use LEKF_state, only : kf_with_xb, kf_with_xm
  use LEKF_state, only : Ens_Mean, Ens_Mean_and_Sigma
  use LEKF_state, only : Ens_Setup
  use LEKF_state, only : Save_State, Restore_State
  !use LEKF_state, only : SendBudgetData

  use LEKF_noise, only : nnoise, nhist
  use LEKF_Noise, only : LEKF_Noise_Init, LEKF_Noise_Done
  use LEKF_Noise, only : SetDC, ResetDC, SetDC_Means

  use LEKF_Driver, only : LEKF_Arguments  
  use LEKF_Driver, only : LEKF_Driver_Init_State
  use LEKF_Driver, only : LEKF_TimeStep_Run
  use LEKF_Driver, only : LEKF_CalcMeas
  use LEKF_Driver, only : TruncateStates
  !use LEKF_Driver, only : biascorr_kf_surface_ozone
  
  use LEKF_Meas_Update, only : T_Meas_Update
  use LEKF_Meas_Update, only : LEKF_Meas_Update_Timers

#ifdef with_kf_meas_maori
  use LE_MAORI       , only : LE_MAORI_Data_Init, LE_MAORI_Data_Done, LE_MAORI_Data_Setup
  use LEKF_Data      , only : mad
#endif
#ifdef with_kf_meas_sat
  use LEKF_Meas_Sat    , only : meas_sat
#endif
#ifdef with_kf_meas_modis
  use LEKF_Meas_MODIS, only : InitMeas_MODIS, DoneMeas_MODIS
  use LEKF_Meas_MODIS, only : GetMeas_MODIS, MeasUpdate_MODIS
#endif

  use LEKF_Output   , only : LEKF_Output_Init, LEKF_Output_Done
  use LEKF_Output   , only : LEKF_Output_Setup, LEKF_Output_PutOut
  use LEKF_Output_DC, only : ReadDC
  
  implicit none

  ! --- const ------------------------------

  character(len=*), parameter   ::  rname = 'LEKF'

  ! --- local ------------------------------

  integer             ::  i, j
!  integer             ::  imode
  logical             ::  isfirst
  logical             ::  xb_first

  ! rcfile:
  character(len=1024)   ::  rcfile
  type(TRcFile)         ::  rcF
  
  ! logging:
  !character(len=512)    ::  outpath
  character(len=512)    ::  logfile


  ! current time and time step:
  type(TDate)           ::  t
  type(TIncrDate)       ::  dt
  logical               ::  the_end

  ! the time in hours from January 1th
  integer               ::  nhour
  
  ! return status:
  integer               ::  status

  ! timers:
  integer   ::  itim_init, itim_first, itim_run, itim_done
  integer   ::  itim_init_model, itim_init_maori, itim_init_filter, &
                 itim_init_states, itim_init_output, itim_init_meas
  integer   ::  itim_noise
  integer   ::  itim_prop_init, itim_prop_run, itim_prop_done
  integer   ::  itim_meas_load, itim_meas_sim
  integer   ::  itim_analysis
  integer   ::  itim_output_setup, itim_output_step
  integer   ::  itim_restart

  ! propagate ensemble ? analyse measurements ?
  logical             ::  kf_with_propagation
  !logical             ::  kf_with_analysis

!  ! time correlation in noise:
!  real                ::  kf_noise_sig
!  real                ::  kf_noise_tau_days
  ! refresh noise:
  real                ::  kf_noise_dt_hours

  ! restart stuff
  character(len=512)  ::  filekey
  logical             ::  restart_xb_from_mean

  ! specials to do something special with dc values
  logical             ::  restart_as_forecast
  character(len=16)   ::  restart_xb_dc_option
  integer             ::  restart_xb_dc_starthour
  integer             ::  restart_xb_dc_endhour
  integer             ::  restart_xb_dc_days
  integer, parameter  ::  dc_nt = 24
  integer, parameter  ::  dc_lag = 1
  real, allocatable   ::  dc_day(:,:,:,:)
  character(len=512)  ::  restart_xb_dc_dir
  character(len=32)   ::  restart_xb_dc_model
  character(len=32)   ::  restart_xb_dc_expid
  type(TDate)         ::  t_dc
  integer             ::  iday
  real                ::  dc_weight

  ! udpate data:
  type(T_Meas_Update) ::  Meas_Update
  
  ! --- begin --------------------------------
  
  ! init GO routines:
  call GO_Init( status )
  IF_NOTOK_RETURN(status=1)

  !! init logging, include process number in messages
  !! (not really necessary if 'mpiexec' is configured 
  !! to write stdout to seperate files per task):
  !call GO_Print_Init( status, prompt_pe=.true. )
  !IF_NOTOK_RETURN(status=1)
  
  !! write messages to individual log files per task;
  !! not needed if mpiexec is taking care of this:
  !write (logfile,'("lekf-",i2.2,".log")') goc%myid
  !! switch:
  !call GO_Print_Logfile( status, file=trim(logfile) )
  !IF_NOTOK_RETURN(status=1) 

  ! extract argument, return name of rcfile:
  call LEKF_Arguments( rcfile, status )
  IF_NOTOK_RETURN(status=1)

  ! open settings:
  call rcF%Init( trim(rcfile), status )
  IF_NOTOK_RETURN(status=1)
  
  ! info ...
  call LE_Message_Init( status )
  IF_NOTOK_RETURN(status=1)

  ! *

  ! start timing:
  call LE_Timers_Init( status )
  IF_NOTOK_RETURN(status=1)

  ! define timers:
  call GO_Timer_Def( itim_init        , 'lekf init'             , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_model  , 'lekf init model'       , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_maori  , 'lekf init maori'       , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_filter , 'lekf init filter'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_states , 'lekf init states'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_output , 'lekf init output'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_init_meas   , 'lekf init meas'        , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_first       , 'lekf first'            , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_run         , 'lekf run'              , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_done        , 'lekf done'             , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_noise       , 'lekf noise'            , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_prop_init   , 'lekf propagation init' , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_prop_run    , 'lekf propagation run'  , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_prop_done   , 'lekf propagation done' , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_meas_load   , 'lekf meas load'        , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_meas_sim    , 'lekf meas sim'         , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_analysis    , 'lekf analysis'         , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_output_setup, 'lekf output setup'     , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_output_step , 'lekf output step'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_restart     , 'lekf restart'          , status )
  IF_NOTOK_RETURN(status=1)

  ! start timing:
  call LEKF_Meas_Update_Timers( status )
  IF_NOTOK_RETURN(status=1)

  !
  ! setup
  !

  ! start timing:
  call GO_Timer_Start(  itim_init, status )
  IF_NOTOK_RETURN(status=1)
  
  ! ~~ model

  ! switch timing:
  call GO_Timer_Start( itim_init_model, status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'(a,": init LOTOS-EUROS ...")') rname; call goPr
  
  ! setup model configuration,
  ! pass the rcfile read from lekf arguments:
  call LE_Config_Init( status, rcfile=rcfile )
  IF_NOTOK_RETURN(status=1)

  ! init model, return:
  !  o start time, time step
  !  o hours since start
  call LE_Model_Init( t, dt, nhour, status )
  IF_NOTOK_RETURN(status=1)

  
  ! ~~ maori

  ! switch timing:
  call GO_Timer_Switch( itim_init_model, itim_init_maori, status )
  IF_NOTOK_RETURN(status=1)
  
#ifdef with_kf_meas_maori
  ! init MAORI stuff:
  call LE_MAORI_Data_Init( mad, trim(rcfile), t, status )
  IF_NOTOK_RETURN(status=1)
#endif

  
  ! ~~ filter

  ! switch timing:
  call GO_Timer_Switch( itim_init_maori, itim_init_filter, status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'(a,": setup filter ...")') rname; call goPr

  ! info ...
  write (gol,'(a,":   read flags to enable filter elements ...")') rname; call goPr

  ! with background run?
  call rcF%Get( 'kf.with.xb', kf_with_xb, status )
  IF_NOTOK_RETURN(status=1)
  !! root only ...
  !kf_with_xb = kf_with_xb .and. goc%root

  ! with modes?
  call rcF%Get( 'kf.with.xm', kf_with_xm, status )
  IF_NOTOK_RETURN(status=1)

  ! for testing, model propagation might be skipped:
  call rcF%Get( 'kf.with_propagation', kf_with_propagation, status )
  IF_NOTOK_RETURN(status=1)
  !call rcF%Get( 'kf.with_analysis', kf_with_analysis, status )
  !IF_NOTOK_RETURN(status=1)

  ! setup noise:
  call LEKF_Noise_Init( rcF, status )
  IF_NOTOK_RETURN(status=1)

  call rcF%Get( 'kf.noise.dt_hours', kf_noise_dt_hours, status )
  IF_NOTOK_RETURN(status=1)

  call rcF%Get( 'kf.restart', runF%restart, status )
  IF_NOTOK_RETURN(status=1)

  ! to allow reading of dc into xb ...
  call rcF%Get( 'kf.restart.as_forecast', restart_as_forecast, status )
  IF_NOTOK_RETURN(status=1)

  if ( runF%restart ) then
    ! where to read from:
    call rcF%Get( 'kf.restart.path', runF%restart_path, status )
    IF_NOTOK_RETURN(status=1)
    ! file name description:
    call rcF%Get( 'kf.restart.key', runF%restart_key, status )
    IF_NOTOK_RETURN(status=1)
    ! ... adhoc: add a standard 'name' to specify the file from which
    !     old data should be read (expcls, volume) :
    write (runF%restart_key,'(a,";name=",a)') trim(runF%restart_key), 'model'
    ! forecast mode ?
    if ( restart_as_forecast ) then
      ! load mean into model state (for forecast ?)
      call rcF%Get( 'kf.restart.xb.from_mean', restart_xb_from_mean, status )
      IF_NOTOK_RETURN(status=1)
      ! ... adhoc: add a standard 'name' to specify the file from which
      !     old data should be read (expcls, volume) :
      if ( restart_xb_from_mean ) then
         write (runF%restart_key,'(a,";name=",a)') trim(runF%restart_key), 'mean'
      end if 
      ! special options for how to use assimilated dc values:
      call rcF%Get( 'kf.restart.xb.dc_option'   , restart_xb_dc_option   , status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.starthour', restart_xb_dc_starthour, status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.endhour'  , restart_xb_dc_endhour  , status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.days'     , restart_xb_dc_days     , status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.dir'      , restart_xb_dc_dir      , status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.model'    , restart_xb_dc_model    , status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( 'kf.restart.xb.dc.expid'    , restart_xb_dc_expid    , status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      restart_xb_from_mean = .false.
    end if ! forecast mode
  else
    ! dummy ...
    runF%restart_path    = 'none'
    runF%restart_key     = 'none'
    restart_xb_from_mean = .false.
  end if

  ! write restart files ?
  call rcF%Get( 'kf.restart.save', runF%restart_save, status )
  IF_NOTOK_RETURN(status=1)

  ! write restart files ?
  if ( runF%restart_save ) then
     ! when to write:
     call rcF%Get( 'kf.restart.save.dhour', runF%restart_save_dhour, status, default=24.0 )
     IF_NOTOK_RETURN(status=1)
     ! where to write:
     call rcF%Get( 'kf.restart.save.path', runF%restart_save_path, status )
     IF_NOTOK_RETURN(status=1)
     ! file name description:
     call rcF%Get( 'kf.restart.save.key', runF%restart_save_key, status )
     IF_NOTOK_RETURN(status=1)
  else
     ! dummy ...
     runF%restart_save_dhour = 9999.9
     runF%restart_save_path  = 'none'
     runF%restart_save_key   = 'none'
  end if

  !! ozone bias correction:
  !call rcF%Get( 'kf.biascorr.surface_ozone', &
  !                 biascorr_kf_surface_ozone, status ) !, default='none' )
  !IF_NOTOK_RETURN(status=1)


  ! first time step is true
  runF%first = .true.

  ! switch timing:
  call GO_Timer_Switch( itim_init_filter, itim_init_output, status )
  IF_NOTOK_RETURN(status=1)
  
  ! info ...
  write (gol,'(a,":   setup output ...")') rname; call goPr
  ! init output stuff
  call LEKF_Output_Init( rcF, t, status )
  IF_NOTOK_RETURN(status=1)

  ! switch timing:
  call GO_Timer_Switch( itim_init_output, itim_init_states, status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'(a,":   setup states ...")') rname; call goPr
  ! setup state module:
  call LEKF_State_Init( rcF, status )
  IF_NOTOK_RETURN(status=1)

  ! switch timing:
  call GO_Timer_Switch( itim_init_states, itim_init_meas, status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'(a,":   setup measurements ...")') rname; call goPr
  ! initialize the measurements
#ifdef with_kf_meas_sat
  call meas_sat%Init( rcF, 'kf.meas.sat', status )
  IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_modis
  call InitMeas_MODIS( status )
  IF_NOTOK_RETURN(status=1)
#endif
  ! setup measurement update:
  call Meas_Update%Init( rcF, status )
  IF_NOTOK_RETURN(status=1)
  
  ! done with settings:
  call rcF%Done( status )
  IF_NOTOK_RETURN(status=1)

  ! end timing:
  call GO_Timer_End( itim_init_meas, status )
  IF_NOTOK_RETURN(status=1)

  ! end timing:
  call GO_Timer_End( itim_init, status )
  IF_NOTOK_RETURN(status=1)


  !
  ! main
  !

  ! start timing:
  call GO_Timer_Start( itim_run, status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'(a,": start time loop ...")') rname; call goPr

  ! set flag:
  isfirst = .true.

  ! the time loop 
  do

    ! info ...
    write (gol,'(a,": ")') rname; call goPr
    call wrtgol( rname//': >>> ', t ); call goPr
    write (gol,'(a,": ")') rname; call goPr

    !
    ! init time step
    !

    ! start timing:
    call GO_Timer_Start( itim_prop_init, status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    call wrtgol( rname//': prepare timestep ', t, ' to ', t+dt ); call goPr

    ! setup data:
    call LE_TimeStep_Init( t, dt, nhour, status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_prop_init, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! first init
    !

    ! init states ?
    if ( isfirst ) then
      write (gol,'(a,": setup for start ...")') rname; call goPr

      ! start timing:
      call GO_Timer_Start( itim_first, status )
      IF_NOTOK_RETURN(status=1)

      ! ensemble present?
      if ( kf_with_xm ) then
        !
        ! init ensemble members:
        do j = 1, nmodes
          !! global index:
          !imode = imodes(j)
          ! info ...
          write (gol,'(a,":   initialize mode : ",i3)') rname, j; call goPr
          ! new or restore ?
          if ( runF%restart ) then
            ! name of mode:
            write (filekey,'(a,";name=ensm",i3.3)') trim(runF%restart_key), j
            ! read state from file:
            call Restore_State( Ens(j), t, trim(runF%restart_path), trim(filekey), status )
            IF_NOTOK_RETURN(status=1)
          else
            ! init ensemble member: 
            call LEKF_Driver_Init_State( Ens(j), t, status )
            IF_NOTOK_RETURN(status=1)
            ! set mean values:
            call SetDC_Means( Ens(j)%dc(:,:,:,1), status )
            IF_NOTOK_RETURN(status=1)
            ! set new noise
            do i = 1, nnoise
              ! random factor around mean; truncate negative values:
              call SetDC( Ens(j)%dc, i, Ens(j)%rnd, status )
              IF_NOTOK_RETURN(status=1)
            enddo
          end if
        end do
        !
        ! new ensemble mean:
        call Ens_Mean( status )
        IF_NOTOK_RETURN(status=1)
        !
      end if ! with xm

      ! set base run (root only):
      if ( kf_with_xb ) then
        ! new or restore ?
        if ( runF%restart ) then
          if ( restart_xb_from_mean ) then
            write (gol,'("LEKF: restore state from ensemble mean")'); call goPr
            ! name of state:
            write (filekey,'(a,";name=mean")') trim(runF%restart_key)
            ! read state from file:
            call Restore_State( xb, t, trim(runF%restart_path), trim(filekey), status )
            IF_NOTOK_RETURN(status=1)
          else
            write (gol,'("LEKF: restore state from model")'); call goPr
            ! name of state:
            write (filekey,'(a,";name=model")') trim(runF%restart_key)
            ! read state from file:
            call Restore_State( xb, t, trim(runF%restart_path), trim(filekey), status )
            IF_NOTOK_RETURN(status=1)
          end if
        else
          ! initializes the state vector ! no 24 measurements...
          call LEKF_Driver_Init_State( xb, t, status )
          IF_NOTOK_RETURN(status=1)
          ! set current dc values to mean value:
          call SetDC_Means( xb%dc(:,:,:,1), status )
          IF_NOTOK_RETURN(status=1)
        end if
      end if  ! with xb

      ! info ...
      write (gol,'("LEKF: setup output ...")'); call goPr

      ! setup output for current interval, this is not the end yet;
      ! incl. load of some satelite data valid for this 'interval':
      call LEKF_Output_Setup( t, t, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_maori
      ! setup maori stuff for current time interval if necessary,
      ! this is not the end yet:
      call LE_MAORI_Data_Setup( mad, t, t, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif
      ! setup observations in augmented states of ensemble:
      call Ens_Setup( status )
      IF_NOTOK_RETURN(status=1)
       
      ! ~ load measurements (initial)

      write (gol,'("LEKF: load measurements ...")'); call goPr
      
      ! MAORI data already loaded above by "LE_MAORI_Data_Setup" ...
#ifdef with_kf_meas_modis
      write (gol,'("LEKF:   load MODIS ...")'); call goPr
      call GetMeas_MODIS( status )
      IF_NOTOK_RETURN(status=1)
#endif

      ! ~ fill simulated measurements in state augmentations (initial)

      ! background available ?
      if ( kf_with_xb ) then
        call LEKF_CalcMeas( xb, t, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! mean/modes available ?
      if ( kf_with_xm ) then
        ! update the present modes
        do i = 1, nmodes
          call LEKF_CalcMeas( Ens(i), t, status )
          IF_NOTOK_RETURN(status=1)
        enddo
        ! recompute ensemble mean/stdv:
        call Ens_Mean_and_Sigma( status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! ~  put out (initial)

      call wrtgol( 'LEKF: put out initial for ', t ); call goPr

      ! try to put out for start time; both forecast and analysis ...
      call LEKF_Output_PutOut( 'forecast', t,status )
      IF_NOTOK_RETURN(status=1)
      call LEKF_Output_PutOut( 'analysis', t,status )
      IF_NOTOK_RETURN(status=1)

      ! ~ finished initialization

      ! end timing:
      call GO_Timer_End( itim_first, status )
      IF_NOTOK_RETURN(status=1)

    end if  ! first


    !
    ! new noise
    !

    !
    ! colored noise with zero mean, std.dev. s, and time correlation tau :
    !
    !   g(t)           ~ N(0,s**2)
    !
    !   E[g(t)g(t+dt)] = exp( - |dt|/tau )
    !
    ! implemented with:
    !
    !   g(t+dt) = a(dt) g(t) + sqrt(1-a(dt)**2) s w(t)
    !
    ! where:
    !
    !   a(dt) = exp( - |dt|/tau )   : time correlation parameter
    !

    ! start timing:
    call GO_Timer_Start( itim_noise, status )
    IF_NOTOK_RETURN(status=1)

    if ( Precisely(t,kf_noise_dt_hours,'hour') ) then

       ! set noise for filter ?
       ! for compatability with OpenDA version, 
       ! skip if this is the initial time and SetDC is already called:
       if ( kf_with_xm .and. (.not. isfirst) ) then
          ! shift history:
          do j = 1, nmodes
             call ResetDC( Ens(j)%dc )
          enddo
          ! set new noise
          do j = 1, nmodes
             do i = 1, nnoise
                ! new noise ~ N(mu,sigma^2) with time correlation:
                call SetDC( Ens(j)%dc, i, Ens(j)%rnd, status, &
                              dt_days=kf_noise_dt_hours/24.0 )
                IF_NOTOK_RETURN(status=1)
             end do ! noise numbers
          end do ! modes
          !
          ! new ensemble mean:
          call Ens_Mean( status )
          IF_NOTOK_RETURN(status=1)
          !
       end if  ! with_xm

       ! set noise for background run ?
       if ( kf_with_xb ) then
         ! initialised from a restart ?
         if ( runF%restart .and. restart_as_forecast ) then
           ! optionally read some other dc value ...
           select case ( trim(restart_xb_dc_option) )

             case ( 'default' )

               ! info
               write (gol,'("    --> xb%dc to mean values ...")'); call goPr
               ! set to mean values:
               call SetDC_Means( xb%dc(:,:,:,1), status )
               IF_NOTOK_RETURN(status=1)

             case ( 'latest' )

               ! use value from restart, keep constant;
               ! nothing extra to be done ...

               ! info
               write (gol,'("    --> xb%dc as is, for example from restart ...")'); call goPr

             case ( 'hourly' )

               ! gedeelte dagverloop zie rc file voor tstart en tend
               if ( (t%hour >= restart_xb_dc_starthour) .and. (t%hour <= restart_xb_dc_endhour) ) then
                 write (gol,'("      read stored dc ...")'); call goPr
                 ! storage:
                 allocate( dc_day(nx,ny,nnoise,dc_nt) )
                 ! reset to zero:
                 xb%dc(:,:,:,1) = 0.0
                 ! loop over days:
                 do iday = 1, restart_xb_dc_days
                   ! dc available for day prior to start (assume we start at 00:00)
                   t_dc = runF%t_start - IncrDate(day=1)
                   !t_dc = runF%t_base - IncrDate(day=iday-1)
                   ! weight in sum:
                   dc_weight = 1.0/real(restart_xb_dc_days)
                   ! info
                   write (gol,'("    --> xb%dc from hour ",i2," (",i4,2("-",i2.2),") weight ",f5.2)') &
                                   t%hour,t_dc%year,t_dc%month,t_dc%day, dc_weight ; call goPr
                   ! read dc from file: volle dagverloop
                   call ReadDC( restart_xb_dc_dir, restart_xb_dc_model, restart_xb_dc_expid,&
                               t_dc%year,t_dc%month,t_dc%day, dc_lag, dc_nt, dc_day, status )
                   IF_NOTOK_RETURN(status=1)
                   write (gol,*) '      dc_day range : ', minval(dc_day(:,:,:,t%hour+1)), maxval(dc_day(:,:,:,t%hour+1)); call goPr
                   ! extract for current hour; [00,01) is hour 1 :
                   xb%dc(:,:,:,1) = xb%dc(:,:,:,1) + dc_day(:,:,:,t%hour+1)*dc_weight
                 end do
                 ! info ...
                 write (gol,*) '      dc range : ', minval(xb%dc(:,:,:,1)), maxval(xb%dc(:,:,:,1)); call goPr
                 ! clear:
                 deallocate( dc_day )
               else
                 ! info
                 write (gol,'("    --> xb%dc to mean value (hour not in range) ...")'); call goPr
                 ! set to mean values:
                 call SetDC_Means( xb%dc(:,:,:,1), status )
                 IF_NOTOK_RETURN(status=1)
               end if

             case ( 'average' )

               ! dc available for day prior to start (assume we start at 00:00)
               !t_dc = runF%t_start - IncrDate(day=1)
               t_dc = runF%t_base
               ! info
               write (gol,'("    --> xb%dc average over hours ",i2,"-",i2," (",i4,2("-",i2.2),")")') &
                                  restart_xb_dc_starthour, restart_xb_dc_endhour, &
                                  t_dc%year,t_dc%month,t_dc%day; call goPr
               ! storage:
               allocate( dc_day(nx,ny,nnoise,dc_nt) )
               ! read dc from file:
               call ReadDC( restart_xb_dc_dir, restart_xb_dc_model, restart_xb_dc_expid,&
                             t_dc%year,t_dc%month,t_dc%day, dc_lag, dc_nt, dc_day, status )
               IF_NOTOK_RETURN(status=1)
               ! part of the day: tstart- tend
               xb%dc(:,:,:,1) = sum( dc_day(:,:,:,restart_xb_dc_starthour:restart_xb_dc_endhour), 4 ) / &
                         real(restart_xb_dc_endhour-restart_xb_dc_starthour+1)
               ! clear:
               deallocate( dc_day )

             case ( 'hour24' )

               ! dc available for day prior to start (assume we start at 00:00)
               !t_dc = runF%t_start - IncrDate(day=1)
               t_dc = runF%t_base
               ! info
               write (gol,'("    --> xb%dc from hour ",i2," (",i4,2("-",i2.2),")")') &
                                  dc_nt, t_dc%year,t_dc%month,t_dc%day; call goPr
               ! storage:
               allocate( dc_day(nx,ny,nnoise,dc_nt) )
               ! read dc from file:
               call ReadDC( restart_xb_dc_dir, restart_xb_dc_model, restart_xb_dc_expid,&
                             t_dc%year,t_dc%month,t_dc%day, dc_lag, dc_nt, dc_day, status )
               if (status<0) then
                 write (gol,'("could not read dc")'); call goErr
                 TRACEBACK; status=1; call goExit(status)
               end if
               IF_NOTOK_RETURN(status=1)
               ! latest value:
               xb%dc(:,:,:,1) = dc_day(:,:,:,dc_nt)
               ! clear:
               deallocate( dc_day )

             case default

               ! info
               write (gol,'("unsupported xb dc option : ",a)') trim(restart_xb_dc_option); call goErr
               TRACEBACK; status=1; call goExit(status)

           end select
         end if  ! restart and forecast mode
       end if ! with_xb

    end if ! new noise

    ! end timing:
    call GO_Timer_End( itim_noise, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! forecast step
    !

    ! start timing:
    call GO_Timer_Start( itim_prop_run, status )
    IF_NOTOK_RETURN(status=1)

    if ( kf_with_propagation ) then
    
       ! store initial flag:
       xb_first = runF%first

       !! NOTE: In the serial v1, the budgets were used from the just propagated xb,
       !!     but that would mean that all modes have to wait for xb ...
       !!     Now send around the budget data first, this will give a small difference!
       !! Send some budget data from root (where xb is propagated) to all pe;
       !! this is used to compute deposition velocities:
       !call SendBudgetData( status )
       !IF_NOTOK_RETURN(status=1)

       ! do time step for base run
       if ( kf_with_xb ) then 
          ! info ...
          write (gol,'("LEKF:   propagate background run ...")'); call goPr
          ! run, flag indicates that his is the model run
          ! (if enabled, then the budgets will be updated too):
          call LEKF_TimeStep_Run( t, dt, xb, .true., status )
          IF_NOTOK_RETURN(status=1)
       end if

       ! timestep for mean/modes
       if ( kf_with_xm ) then
          ! update the present modes
          do i = 1, nmodes
             ! restore flag:
             runF%first = xb_first
             ! info ...
             write (gol,'("LEKF:   propagate mode ",i3," ...")') i; call goPr
             ! run for ensemble member; 
             ! might require xb for additive noise
             call LEKF_TimeStep_Run( t, dt, Ens(i), .false., status )!, xb=xb )
             IF_NOTOK_RETURN(status=1)
          enddo
          ! ensemble mean:
          call Ens_Mean( status )
          IF_NOTOK_RETURN(status=1)
       end if

    else

       write (gol,'("WARNING - skip propagation of ensemble ...")'); call goPr

    end if

    ! first time step is no longer true
    runF%first = .false.

    ! end timing:
    call GO_Timer_End( itim_prop_run, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! forecast done
    !

    ! start timing:
    call GO_Timer_Start( itim_prop_done, status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'("LEKF: finish timestep ...")'); call goPr

    ! end of time step;
    !   at return, t := t+dt
    call LE_TimeStep_Done( t, dt, nhour, status )
    IF_NOTOK_RETURN(status=1)

    ! end time ?
    the_end = t == runF%t_end
 
    ! end timing:
    call GO_Timer_End( itim_prop_done, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ setup output
    !

    ! start timing:
    call GO_Timer_Start( itim_output_setup, status )
    IF_NOTOK_RETURN(status=1)
          
    ! setup output for current interval;
    ! incl. load of some satelite data valid for this 'interval':
    call LEKF_Output_Setup( t-dt, t, status )
    IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_maori
    ! setup maori stuff for current time interval if necessary:
    call LE_MAORI_Data_Setup( mad, t-dt, t, the_end, status )
    IF_NOTOK_RETURN(status=1)
#endif
    ! setup observations in augmented states of ensemble:
    call Ens_Setup( status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_output_setup, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ~~ load measurements
    !

    ! start timing:
    call GO_Timer_Start( itim_meas_load, status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'("LEKF: load measurements ...")'); call goPr

    ! MAORI data already loaded above by "LE_MAORI_Data_Setup" ...
#ifdef with_kf_meas_modis
    write (gol,'("LEKF:   load MODIS ...")'); call goPr
    call GetMeas_MODIS( status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! end timing:
    call GO_Timer_End( itim_meas_load, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! fill simulated measurements in state augmentations
    !

    ! start timing:
    call GO_Timer_Start( itim_meas_sim, status )
    IF_NOTOK_RETURN(status=1)

    ! background available ?
    if ( kf_with_xb ) then
      call LEKF_CalcMeas( xb, t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! mean/modes available ?
    if ( kf_with_xm ) then
      ! update the present modes
      do i = 1, nmodes
        call LEKF_CalcMeas( Ens(i), t, status )
        IF_NOTOK_RETURN(status=1)
      enddo
      ! recompute ensemble mean:
      call Ens_Mean_and_Sigma( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! end timing:
    call GO_Timer_End( itim_meas_sim, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ output (forecast)
    !

    ! start timing:
    call GO_Timer_Start( itim_output_step, status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    call wrtgol( 'LEKF: put out forecast for ', t ); call goPr

    ! try to put out for this time; forecast before assim only !
    call LEKF_Output_PutOut( 'forecast', t, status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_output_step, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ analysis
    !

    ! start timing:
    call GO_Timer_Start( itim_analysis, status )
    IF_NOTOK_RETURN(status=1)

    write (gol,'("LEKF: analysis ...")'); call goPr

    ! filter application ?
    if ( kf_with_xm ) then
   
      ! measurement update:
      call Meas_Update%Apply( t-dt, t, status )
      IF_NOTOK_RETURN(status=1)

!#ifdef with_kf_meas_maori
!       ! analyse Sample measurents
!       call MeasUpdate_MAORI( nmodes, status )
!       IF_NOTOK_RETURN(status=1)
!#endif
!#ifdef with_kf_meas_modis
!       ! analyse MODIS measurents
!       call MeasUpdate_MODIS( nmodes, status )
!       IF_NOTOK_RETURN(status=1)
!#endif

      ! remove negative concentrations etc;
      ! also updated: measurements x%m, ensemble mean 'x', and stdv 'sigma':
      call TruncateStates( status )
      IF_NOTOK_RETURN(status=1)

      ! update measurement simulations in modes:
      do i = 1, nmodes
        call LEKF_CalcMeas( Ens(i), t, status )
        IF_NOTOK_RETURN(status=1)
      enddo
      ! recompute ensemble mean:
      call Ens_Mean_and_Sigma( status )
      IF_NOTOK_RETURN(status=1)

    end if  ! mean, modes

    ! end timing:
    call GO_Timer_End( itim_analysis, status )
    IF_NOTOK_RETURN(status=1)


    !
    ! ~~ output (analysis)
    !

    ! start timing:
    call GO_Timer_Start( itim_output_step, status )
    IF_NOTOK_RETURN(status=1)

    call wrtgol( 'LEKF: put out analysis for ', t ); call goPr

    ! try to put out for this time:
    call LEKF_Output_PutOut( 'analysis', t, status, last=the_end )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_output_step, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ restart stuff
    !

    ! start timing:
    call GO_Timer_Start( itim_restart, status )
    IF_NOTOK_RETURN(status=1)

    ! save for restart ?
    if ( runF%restart_save ) then
      ! every 24 hour ...
      if ( Precisely(t,runF%restart_save_dhour,'hour') ) then
        ! ensemble present?
        if ( kf_with_xm ) then
          ! loop over ensemble members
          do j = 1, nmodes
            !! global index:
            !imode = imodes(j)
            ! name of ensemble member:
            write (filekey,'(a,";name=ensm",i3.3)') trim(runF%restart_save_key), j
            ! read state from file:
            call Save_State( Ens(j), t, trim(runF%restart_save_path), trim(filekey), status )
            IF_NOTOK_RETURN(status=1)
          end do
          ! ensure that mean/sigma is updated:
          call Ens_Mean_and_Sigma( status )
          IF_NOTOK_RETURN(status=1)
          ! target filename:
          write (filekey,'(a,";name=mean")') trim(runF%restart_save_key)
          ! write state to file:
          call Save_State( x, t, trim(runF%restart_save_path), trim(filekey), status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! model run present?
        if ( kf_with_xb ) then
          ! save background run:
          write (filekey,'(a,";name=model")') trim(runF%restart_save_key)
          ! write state to file:
          call Save_State( xb, t, trim(runF%restart_save_path), trim(filekey), status )
          IF_NOTOK_RETURN(status=1)
        end if
      end if
    end if

    ! end timing:
    call GO_Timer_End( itim_restart, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ end time ?
    !

        !! testing ..
        !write (gol,'("WARNING - break in time loop!")'); call goPr
        !the_end = .true.

    ! check if end of simulation is reached
    if ( the_end ) exit

    ! reset flag:l
    isfirst = .false.

  end do ! time loop

  ! info ...
  write (gol,'("LEKF: end time loop")'); call goPr

  ! switch timing:
  call GO_Timer_Switch( itim_run, itim_done, status )
  IF_NOTOK_RETURN(status=1)

  write (gol,'("LEKF: done ...")'); call goPr

  ! info ..
  write (gol,'("LEKF:   done with measurements ...")'); call goPr
#ifdef with_kf_meas_sat
  call meas_sat%Done( status )
  IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_modis
  call DoneMeas_MODIS( status )
  IF_NOTOK_RETURN(status=1)
#endif
  ! done with measurement update:
  call Meas_Update%Done( status )
  IF_NOTOK_RETURN(status=1)


  ! info ...
  write (gol,'("LEKF:   done with filter output ...")'); call goPr

  ! done with output stuff
  call LEKF_Output_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'("LEKF:   done with states ...")'); call goPr
  ! done with states:
  call LEKF_State_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! info ...
  write (gol,'("LEKF:   done with noise ...")'); call goPr
  ! done with noise:
  call LEKF_Noise_Done( status )
  IF_NOTOK_RETURN(status=1)
  
! moved to lekf.F90      
#ifdef with_kf_meas_maori
  ! info ...
  write (gol,'("LEKF:   done with maori ...")'); call goPr
  ! done with MAORI stuff:
  call LE_MAORI_Data_Done( mad, status )
  IF_NOTOK_RETURN(status=1)
#endif

  ! info ...
  write (gol,'("LEKF:   done with model ...")'); call goPr
  ! done with model:
  call LE_Model_Done( status )
  IF_NOTOK_RETURN(status=1)
  
  ! done with model configuration:
  call LE_Config_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! end timing:
  call GO_Timer_End( itim_done, status )
  IF_NOTOK_RETURN(status=1)

  ! *

  ! stop timing, print timing profile:
  call LE_Timers_Done( status, write_profile=goc%root )
  IF_NOTOK_RETURN(status=1)

  ! *

        !! dummy ..
        !open( unit=123, file='le.ok' )
        !write (123,'("Program terminated normally")')
        !close( unit=123 )

  ! info ...
  call LE_Message_Done( status )
  IF_NOTOK_RETURN(status=1)
  
  !! done with logging:
  !call GO_Print_Done( status )
  !IF_NOTOK_RETURN(status=1)
  
  ! done with GO routines:
  call GO_Done( status )
  IF_NOTOK_RETURN(status=1)

end program LOTOS_EUROS_KF
