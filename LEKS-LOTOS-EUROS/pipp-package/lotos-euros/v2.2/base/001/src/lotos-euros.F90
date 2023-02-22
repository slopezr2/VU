!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; call goc%Abort(status); end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; call goc%Abort(status); end if
!
#include "le.inc"
!
!#################################################################

program LOTOS_EUROS

  use GO, only : GO_Init, GO_Done
  use GO, only : gol, goPr, goErr
  use GO, only : goc

  use GO, only : ReadRc
  use GO, only : TDate, TIncrDate
  use GO, only : operator(==), operator(-), operator(+), wrtgol, Precisely

  use GO, only : GO_Timer_Def, GO_Timer_Start, GO_Timer_End
  
  use dims       , only : nx,  ny, nz, nspec
  use dims       , only : runF

  use LE_Config  , only : LE_Config_Init, LE_Config_Done, rcF
  use LE_Driver  , only : LE_Timers_Init, LE_Timers_Done
  use LE_Driver  , only : LE_Message_Init, LE_Message_Done
  use LE_Driver  , only : LE_Model_Init, LE_Model_Done
  use LE_Driver  , only : LE_State_Init
  use LE_Driver  , only : LE_TimeStep_Init, LE_TimeStep_Run, LE_TimeStep_Done

  use LE_Particle_Data, only : LE_Particle_Data_Update

  use LE_Budget  , only : T_Budget, Budget_Init, Budget_Done

  use LE_DryDepos, only : LE_DryDepos_Setup_vd
  use LE_DryDepos, only : mix2ground
  
  use LE_BiasCorr, only : LE_BiasCorr_Fill


  use LE_Output  , only : T_LE_Output, T_LE_OutputState, Init, Done, Setup, PutOut, EndPutOut

  use LE_Restart , only : LE_Restart_Save, LE_Restart_Restore_State


  implicit none


  ! --- const ------------------------------

  character(len=*), parameter   ::  rname = 'LOTOS-EUROS'

  ! --- local ------------------------------

  ! the concentration vector and budget arrays
  real, allocatable         ::  c(:,:,:,:)     ! (nx,ny,nz,nspec)
  real, allocatable         ::  cg(:,:,:)      ! (nx,ny,nspec)
  real, allocatable         ::  aerh2o(:,:,:)  ! (nx,ny,nz)
  type(T_Budget)            ::  bud

  ! timers:
  integer                   ::  itim_first
  integer                   ::  itim_init, itim_timeloop, itim_done
  integer                   ::  itim_output, itim_save

  ! output
  type(T_LE_Output)         ::  leo
  type(T_LE_OutputState)    ::  leos

  ! current time and time step:
  type(TDate)               ::  t
  type(TIncrDate)           ::  dt

  ! the time in hours from start of run:
  integer                   ::  nhour

  ! return status:
  integer                   ::  status

  ! --- begin --------------------------------------

  ! init GO routines:
  call GO_Init( status )
  IF_NOTOK_RETURN(status=1)
  
  ! general config:
  call LE_Config_Init( status )
  IF_NOTOK_RETURN(status=1)

  ! prepare info ...
  call LE_Message_Init( status )
  IF_NOTOK_RETURN(status=1)

  ! start timing:
  call LE_Timers_Init( status )
  IF_NOTOK_RETURN(status=1)

  ! define timers:
  call GO_Timer_Def( itim_init    , 'model init'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_timeloop, 'model time loop' , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_done    , 'model done'      , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_first   , 'time step first' , status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_output  , 'time step output', status )
  IF_NOTOK_RETURN(status=1)
  call GO_Timer_Def( itim_save    , 'time step save'  , status )
  IF_NOTOK_RETURN(status=1)

  ! *

  call GO_Timer_Start(itim_init, status )
  IF_NOTOK_RETURN(status=1)

  ! init model, return:
  !  o start time, time step
  !  o hours since start
  call LE_Model_Init( t, dt, nhour, status )
  IF_NOTOK_RETURN(status=1)

  ! init output stuff:
  call Init( leo, rcF, 'le.output', status )
  IF_NOTOK_RETURN(status=1)
  call Init( leos, leo, rcF, '', status )
  IF_NOTOK_RETURN(status=1)
    
  ! state:
  allocate( c     (nx,ny,nz,nspec), source=0.0, stat=status )
  IF_NOTOK_RETURN(status=1)
  allocate( cg    (nx,ny   ,nspec), source=0.0, stat=status )
  IF_NOTOK_RETURN(status=1)
  allocate( aerh2o(nx,ny,nz      ), source=0.0, stat=status )
  IF_NOTOK_RETURN(status=1)
  
  ! budgets:
  call Budget_Init( bud, status )
  IF_NOTOK_RETURN(status=1)

  ! end timing:
  call GO_Timer_End(itim_init, status )
  IF_NOTOK_RETURN(status=1)

  ! *

  ! set up time loop with time step of one hour
  ! below, this is subdivided into smaller time steps (Strang operator splitting)
  do

    ! info ...
    write (gol,'("LE: ")'); call goPr
    call wrtgol( 'LE: >>> ', t ); call goPr
    write (gol,'("LE: ")'); call goPr
    
    ! start timing:
    call GO_Timer_Start(itim_timeloop, status )
    IF_NOTOK_RETURN(status=1)

    ! setup data:
    call LE_TimeStep_Init( t, dt, nhour, status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End(itim_timeloop, status )
    IF_NOTOK_RETURN(status=1)

    ! initial concentrations:
    if ( runF%first ) then

      ! start timing:
      call GO_Timer_Start(itim_first, status )
      IF_NOTOK_RETURN(status=1)

      ! setup concentrations:
      call LE_State_Init( c, t, status )
      IF_NOTOK_RETURN(status=1)

      ! overwrite with restart file ?
      if ( runF%restart ) then
        ! read state from file:
        call LE_Restart_Restore_State( c, cg, aerh2o, bud, t, &
                                        trim(runF%restart_path), &
                                        trim(runF%restart_key), &
                                        status )
        IF_NOTOK_RETURN(status=1)
      end if
    
      ! compute new M7 wet radius, particle size, slip correction, etc:
      call LE_Particle_Data_Update( c, status )
      IF_NOTOK_RETURN(status=1)

      ! re-compute deposition velocities,
      ! not only those that depend on the chemical regime:
      call LE_DryDepos_Setup_vd( c, bud%drydepos, .false., t, status )
      IF_NOTOK_RETURN(status=1)

      ! surface concentrations:
      if ( .not. runF%restart ) then  
        ! In case of restart, 
        ! --Ra/vd/vs of previous timestep not available in restart option
        ! --surface concentrations are read from restart file
        call mix2ground( c, cg, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! fill bias corrected tracers:
      call LE_BiasCorr_Fill( c, cg, status )
      IF_NOTOK_RETURN(status=1)


      ! setup output for current time 'interval' if necessary:
      call Setup( leo, t, t, status )
      IF_NOTOK_RETURN(status=1)
      
      ! put out data if necessary:
      call PutOut( leo, t, status )
      IF_NOTOK_RETURN(status=1)
      ! put out concentrations if necessary:
      call PutOut( leos, leo, t, c, cg, bud, status )
      IF_NOTOK_RETURN(status=1)
      ! post process if necessary:
      call EndPutOut( leo, t, status )
      IF_NOTOK_RETURN(status=1)

      ! end timing:
      call GO_Timer_End(itim_first, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! start timing:
    call GO_Timer_Start(itim_timeloop, status )
    IF_NOTOK_RETURN(status=1)

    ! time step, update budgets:
    call LE_TimeStep_Run( t, dt, c, cg, aerh2o, bud, .true., status )
    IF_NOTOK_RETURN(status=1)

    ! end of time step;
    !   at return, t := t+dt
    call LE_TimeStep_Done( t, dt, nhour, status )
    IF_NOTOK_RETURN(status=1)

    ! save for restart ?
    if ( runF%restart_save ) then
      if ( Precisely(t,runF%restart_save_dhour,'hour') ) then
        ! start timing:
        call GO_Timer_Start(itim_save, status )
        IF_NOTOK_RETURN(status=1)
        ! write state to file:
        call LE_Restart_Save( c, cg, aerh2o, &
                              bud, t, &
                              trim(runF%restart_save_path), &
                              trim(runF%restart_save_key), &
                              status )
        IF_NOTOK_RETURN(status=1)
        
        ! end timing:
        call GO_Timer_End(itim_save, status )
        IF_NOTOK_RETURN(status=1)
      end if
    end if

    ! start timing:
    call GO_Timer_Start(itim_output, status )
    IF_NOTOK_RETURN(status=1)
    ! setup output for current time interval if necessary:
    call Setup( leo, t-dt, t, status )
    IF_NOTOK_RETURN(status=1)
    ! put out data if necessary:
    ! (note: meteo is valid for time halfway the time step)
    call PutOut( leo, t, status )
    IF_NOTOK_RETURN(status=1)
    ! put out concentrations if necessary:
    call PutOut( leos, leo, t, c, cg, bud, status )
    IF_NOTOK_RETURN(status=1)
    ! post process if necessary:
    call EndPutOut( leo, t, status )
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End(itim_output, status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End(itim_timeloop, status )
    IF_NOTOK_RETURN(status=1)

    ! end time ? then leave:
    if ( t == runF%t_end ) then
      ! info ...
      write (gol,'("LE: ")'); call goPr
      write (gol,'("LE: <<< end time reached")'); call goPr
      write (gol,'("LE: ")'); call goPr
      ! leave time loop:
      exit
    end if

  end do ! time loop

  ! *

  ! start timing:
  call GO_Timer_Start(itim_done, status )
  IF_NOTOK_RETURN(status=1)

  ! done with output stuff
  call Done( leos, leo, status )
  IF_NOTOK_RETURN(status=1)
  call Done( leo, status )
  IF_NOTOK_RETURN(status=1)

  ! done with budgets:
  call Budget_Done( bud, status )
  IF_NOTOK_RETURN(status=1)
    
  ! clear:
  deallocate( c, stat=status )
  IF_NOTOK_RETURN(status=1)
  deallocate( cg, stat=status )
  IF_NOTOK_RETURN(status=1)
  deallocate( aerh2o, stat=status )
  IF_NOTOK_RETURN(status=1)
  
  ! done with model:
  call LE_Model_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! end timing:
  call GO_Timer_End(itim_done, status )
  IF_NOTOK_RETURN(status=1)

  ! done with timing, write profile:
  call LE_Timers_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! done with info ...
  call LE_Message_Done( status )
  IF_NOTOK_RETURN(status=1)
  
  ! done with config:
  call LE_Config_Done( status )
  IF_NOTOK_RETURN(status=1)

  ! done with GO routines:
  call GO_Done( status )
  IF_NOTOK_RETURN(status=1)

end program LOTOS_EUROS

