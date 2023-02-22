!###############################################################################
!
! NAME
!
!   KF_OUTPUT  -  KF LOTOS-EUROS output
!
! DESCRIPTION
!
!   Write output for all states in Kalman filter.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################


module LEKF_output

  use GO, only : gol, goPr, goErr
  
  use LE_Output       , only : T_LE_OutputState
#ifdef with_kf_meas_maori
  use LE_MAORI        , only : T_MAORI_Output
#endif
  use LEKF_Output_DC  , only : T_LEKF_Output_DC_State

  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  LEKF_Output_Init, LEKF_Output_Done
  public  ::  LEKF_Output_Setup, LEKF_Output_PutOut

  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Output'

  
  ! --- types --------------------------------
  
  type T_LEKF_Output_State
    type(T_LE_OutputState)                ::  leos
    type(T_MAORI_Output)                  ::  mao
    type(T_LEKF_Output_DC_State)          ::  dcos
  contains
    procedure ::  Init        => LEKF_Output_State_Init
    procedure ::  Done        => LEKF_Output_State_Done
    procedure ::  PutOut      => LEKF_Output_State_PutOut
  end type T_LEKF_Output_State


  ! --- var --------------------------------

  ! output state of model run:
  type(T_LEKF_Output_State)                 ::  lekfo_xb
  ! output of ensemble mean:
  type(T_LEKF_Output_State)                 ::  lekfo_xf
  type(T_LEKF_Output_State)                 ::  lekfo_xa
  ! output of ensemble stdv:
  type(T_LEKF_Output_State)                 ::  lekfo_sf
  type(T_LEKF_Output_State)                 ::  lekfo_sa
  ! output of ensemble members:
  type(T_LEKF_Output_State), allocatable    ::  lekfo_xif(:)
  type(T_LEKF_Output_State), allocatable    ::  lekfo_xia(:)

  ! timers:
  integer   ::  itim_init_output_model, itim_init_output_state, itim_init_output_maori
  integer   ::  itim_init_output_meas, itim_init_output_dc
  
    
contains


  ! ====================================================================
  ! ===
  ! === LEKF Output State
  ! ===
  ! ====================================================================
  
  
  subroutine LEKF_Output_State_Init( self, rcF, name, status )

    use GO       , only : TrcFile
    use LE_Output, only : Init
    use LEKF_Data, only : leo
#ifdef with_kf_meas_maori
    use LEKF_Data, only : mad
    use LE_MAORI , only : LE_MAORI_Output_Init
#endif
    
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_State), intent(out)       ::  self
    type(TrcFile), intent(in)                     ::  rcF
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_State_Init'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------

    ! setup standard output of base run concentrations:
    call Init( self%leos, leo, rcF, name, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_maori
    ! init MAORI output stuff:
    call LE_MAORI_Output_Init( self%mao, mad, name, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! init dc output:
    call self%dcos%Init( rcF, 'le.output', name, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LEKF_Output_State_Init
  
  
  ! ***
  

  subroutine LEKF_Output_State_Done( self, status )

    use LEKF_Data, only : leo
    use LE_Output, only : Done
#ifdef with_kf_meas_maori
    use LEKF_Data, only : mad
    use LE_MAORI , only : LE_MAORI_Output_Done
#endif
  
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_State), intent(inout)       ::  self
    integer, intent(out)                            ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_State_Done'
    
    ! --- begin ---------------------------------

    ! done with output of base run:
    call Done( self%leos, leo, status )
    IF_NOTOK_RETURN(status=1)
      
#ifdef with_kf_meas_maori
    ! done with MAORI output stuff:
    call LE_MAORI_Output_Done( self%mao, mad, status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! done with dc output:
    call self%dcos%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LEKF_Output_State_Done
  
  
  ! ***
  
  
  subroutine LEKF_Output_State_PutOut( self, x, t, status, &
                                         without_data )

    use GO        , only : TDate
    use LE_Output , only : PutOut
    use LEKF_State, only : TState
    use LEKF_State, only : bud0
    use LEKF_Data , only : leo
#ifdef with_kf_meas_maori
    use LE_MAORI  , only : LE_MAORI_Output_Write
    use LEKF_Data , only : mad
#endif
 
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_State), intent(inout)       ::  self
    type(TState), intent(inout)                     ::  x
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status
    
    logical, intent(in), optional                   ::  without_data
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_State_PutOut'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------

    ! put out model fields;
    ! only one budget array, filled by xb run:
    call PutOut( self%leos, leo, t, x%c, x%cg, bud0, status, without_data=without_data )
    !print *,'Santiago xxx ', sum(x%c)/size(x%c)
    !print *,'Santiago xxxxx', x%name
    IF_NOTOK_RETURN(status=1)
    
#ifdef with_kf_meas_maori
    ! write via MAORI:
    call LE_MAORI_Output_Write( self%mao, x%mas, mad, status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! put out dc fields:
    call self%dcos%PutOut( x, t, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Output_State_PutOut


  ! ====================================================================
  ! ===
  ! === LEKF Output
  ! ===
  ! ====================================================================
  
  
  subroutine LEKF_Output_Init( rcF, t0, status )

    use GO              , only : goc
    use GO              , only : TDate
    use GO              , only : TrcFile
    !use GO              , only : GO_Timer_Def, GO_Timer_Start, GO_Timer_End, GO_Timer_Switch

    use LEKF_Data       , only : leo
    use LE_Output       , only : Init
    !use LE_Output_Tools , only : itim_init_output_conc, itim_init_output_conc_com, itim_init_output_conc_rc

    use LEKF_Data       , only : lekfo_replace
    use LEKF_Data       , only : lekfo_with_xi
!    use LEKF_Data       , only : leos_xb, leos_xf, leos_sf, leos_xif, leos_xa, leos_sa, leos_xia
!    use LEKF_Output_DC  , only : Init
    use LEKF_State      , only : kf_with_xb, kf_with_xm
    use LEKF_State      , only : nmodes
  
!#ifdef with_kf_meas_maori
!    use LEKF_Data, only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xif, mao_xa, mao_sa, mao_xia
!    !use LE_MAORI , only : LE_MAORI_Data_Init
!    use LE_MAORI , only : LE_MAORI_Output_Init
!#endif

    ! --- in/out --------------------------------
    
    type(TrcFile), intent(in)             ::  rcF
    type(TDate), intent(in)               ::  t0
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_Init'
    
    ! --- local ---------------------------------
    
    integer             ::  j
    character(len=32)   ::  name
    
    ! --- begin ---------------------------------
    
    !! root only ..
    !if ( .not. goc%root) then
    !  status=0; return
    !end if

    !! define timers:
    !call GO_Timer_Def( itim_init_output_model, 'lekf init output model', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_state, 'lekf init output state', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_maori, 'lekf init output maori', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_meas , 'lekf init output meas' , status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_dc   , 'lekf init output dc'   , status )
    !IF_NOTOK_RETURN(status=1)

    !! define timers from model:
    !call GO_Timer_Def( itim_init_output_conc     , 'le init output conc'    , status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_conc_com , 'le init output conc com', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_init_output_conc_rc  , 'le init output conc rc' , status )
    !IF_NOTOK_RETURN(status=1)
    
    ! read total number of modes;
    ! usually done in "LEKF_State_Init", but already required here ...
    call rcF%Get( 'kf.nmodes', nmodes, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files ?
    call rcf%Get( 'kf.output.replace', lekfo_replace, status )
    IF_NOTOK_RETURN(status=1)
    
    ! output for all ensemble members?
    call rcf%Get( 'kf.output.xi', lekfo_with_xi, status )
    IF_NOTOK_RETURN(status=1)

    !! start timing:
    !call GO_Timer_Start( itim_init_output_model, status )
    !IF_NOTOK_RETURN(status=1)
  
    ! setup standard output of model data:
    call Init( leo, rcF, 'le.output', status )
    IF_NOTOK_RETURN(status=1)

    !! end timing:
    !call GO_Timer_End( itim_init_output_model, status )
    !IF_NOTOK_RETURN(status=1)

!! now in lekf.f90 ...    
!#ifdef with_kf_meas_maori
!    ! init MAORI stuff:
!    call LE_MAORI_Data_Init( mad, trim(rcfile), t0, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
    
    ! with background run ?
    if ( kf_with_xb ) then

      !! start timing:
      !call GO_Timer_Start( itim_init_output_state, status )
      !IF_NOTOK_RETURN(status=1)

!      ! setup standard output of base run concentrations:
!      call Init( leos_xb, leo, rcF, 'xb', status )
!      IF_NOTOK_RETURN(status=1)
!
!      !! switch timing:
!      !call GO_Timer_Switch( itim_init_output_state, itim_init_output_maori, status )
!      !IF_NOTOK_RETURN(status=1)
!
!#ifdef with_kf_meas_maori
!      ! init MAORI output stuff:
!      call LE_MAORI_Output_Init( mao_xb, mad, 'xb', status )
!      IF_NOTOK_RETURN(status=1)
!#endif

      ! init output state:
      call lekfo_xb%Init( rcF, 'xb', status )
      IF_NOTOK_RETURN(status=1)

      !! end timing:
      !call GO_Timer_End( itim_init_output_maori, status )
      !IF_NOTOK_RETURN(status=1)

    end if

    ! with mean and modes ?
    if ( kf_with_xm ) then

      !! start timing:
      !call GO_Timer_Start( itim_init_output_state, status )
      !IF_NOTOK_RETURN(status=1)
    
!      ! setup standard output of mean state concentrations:
!      call Init( leos_xf, leo, rcF, 'xf', status )
!      IF_NOTOK_RETURN(status=1)
!      call Init( leos_xa, leo, rcF, 'xa', status )
!      IF_NOTOK_RETURN(status=1)
!
!      ! setup standard output of standard deviation concentrations:
!      call Init( leos_sf, leo, rcF, 'sf', status )
!      IF_NOTOK_RETURN(status=1)
!      call Init( leos_sa, leo, rcF, 'sa', status )
!      IF_NOTOK_RETURN(status=1)
!      
!      ! all members?
!      if ( lekfo_xi ) then
!        ! storage:
!        allocate( leos_xif(nmodes), stat=status )
!        IF_NOTOK_RETURN(status=1)
!        allocate( leos_xia(nmodes), stat=status )
!        IF_NOTOK_RETURN(status=1)
!        ! loop over members:
!        do j = 1, nmodes
!          ! output key:
!          write (name,'("xi",i2.2)') j
!          ! setup output:
!          call Init( leos_xif(j), leo, rcF, trim(name)//'f', status )
!          IF_NOTOK_RETURN(status=1)
!          ! setup output:
!          call Init( leos_xia(j), leo, rcF, trim(name)//'a', status )
!          IF_NOTOK_RETURN(status=1)
!        end do ! ens
!      end if ! put out members
!
!      !! switch timing:
!      !call GO_Timer_Switch( itim_init_output_state, itim_init_output_maori, status )
!      !IF_NOTOK_RETURN(status=1)
!
!#ifdef with_kf_meas_maori
!      ! init MAORI output stuff:
!      call LE_MAORI_Output_Init( mao_xf, mad, 'xf', status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Init( mao_xa, mad, 'xa', status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Init( mao_sf, mad, 'sf', status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Init( mao_sa, mad, 'sa', status )
!      IF_NOTOK_RETURN(status=1)
!      ! all members?
!      if ( lekfo_xi ) then
!        ! storage:
!        allocate( mao_xif(nmodes), stat=status )
!        IF_NOTOK_RETURN(status=1)
!        allocate( mao_xia(nmodes), stat=status )
!        IF_NOTOK_RETURN(status=1)
!        ! loop over members:
!        do j = 1, nmodes
!          ! output key:
!          write (name,'("xi",i2.2)') j
!          ! setup output:
!          call LE_MAORI_Output_Init( mao_xif(j), mad, trim(name)//'f', status )
!          IF_NOTOK_RETURN(status=1)
!          ! setup output:
!          call LE_MAORI_Output_Init( mao_xia(j), mad, trim(name)//'a', status )
!          IF_NOTOK_RETURN(status=1)
!        end do ! ens
!      end if ! put out members
!#endif

      ! init output for ensemble mean/stdv:
      call lekfo_xf%Init( rcF, 'xf', status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_xa%Init( rcF, 'xa', status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_sf%Init( rcF, 'sf', status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_sa%Init( rcF, 'sa', status )
      IF_NOTOK_RETURN(status=1)
      ! all members?
      if ( lekfo_with_xi ) then
        ! storage:
        allocate( lekfo_xif(nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lekfo_xia(nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! loop over members:
        do j = 1, nmodes
          ! output key:
          write (name,'("xi",i2.2)') j
          ! setup output:
          call lekfo_xif(j)%Init( rcF, trim(name)//'f', status )
          IF_NOTOK_RETURN(status=1)
          ! setup output:
          call lekfo_xia(j)%Init( rcF, trim(name)//'a', status )
          IF_NOTOK_RETURN(status=1)
        end do ! ens
      end if ! put out members

      !! end timing:
      !call GO_Timer_End( itim_init_output_maori, status )
      !IF_NOTOK_RETURN(status=1)

    end if

    !! start timing:
    !call GO_Timer_Start( itim_init_output_meas, status )
    !IF_NOTOK_RETURN(status=1)

!    ! setup output of measurements:    
!    call kfom%Init( rcF, 'le.output', status )
!    IF_NOTOK_RETURN(status=1)
!
!    !! switch timing:
!    !call GO_Timer_Switch( itim_init_output_meas, itim_init_output_dc, status )
!    !IF_NOTOK_RETURN(status=1)
!
!    ! setup output of dc values:    
!    call Init( kfodc, rcF, 'le.output', status )
!    IF_NOTOK_RETURN(status=1)

    !! end timing:
    !call GO_Timer_End( itim_init_output_dc, status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Output_Init


  ! ***
  

  subroutine LEKF_Output_Done( status )

    use GO              , only : goc
    use LE_Output       , only : Done
    use LEKF_Data       , only : leo
!    use LEKF_Data       , only : leos_xb, leos_xf, leos_sf, leos_xif, leos_xa, leos_sa, leos_xia
    use LEKF_Data       , only : lekfo_with_xi
!    use LEKF_Output_DC  , only : Done
    use LEKF_State      , only : kf_with_xb, kf_with_xm
    use LEKF_State      , only : nmodes
  
!#ifdef with_kf_meas_maori
!    use LEKF_Data, only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xif, mao_xa, mao_sa, mao_xia
!!    use LE_MAORI , only : LE_MAORI_Data_Done
!    use LE_MAORI , only : LE_MAORI_Output_Done
!#endif
  
    ! --- in/out --------------------------------
    
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_Done'
    
    ! --- local ---------------------------------
    
    integer             ::  j
    
    ! --- begin ---------------------------------
    
    !! root only ..
    !if ( .not. goc%root) then
    !  status=0; return
    !end if
    
    if ( kf_with_xb ) then
!      ! done with output of base run:
!      call Done( leos_xb, leo, status )
!      IF_NOTOK_RETURN(status=1)
!      
!#ifdef with_kf_meas_maori
!      ! done with MAORI output stuff:
!      call LE_MAORI_Output_Done( mao_xb, mad, status )
!      IF_NOTOK_RETURN(status=1)
!#endif
      ! done with output state:
      call lekfo_xb%Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    if ( kf_with_xm ) then
!      ! done with output of mean state:
!      call Done( leos_xf, leo, status )
!      IF_NOTOK_RETURN(status=1)
!      call Done( leos_xa, leo, status )
!      IF_NOTOK_RETURN(status=1)
! 
!      ! done with output of standard deviation:
!      call Done( leos_sf, leo, status )
!      IF_NOTOK_RETURN(status=1)
!      call Done( leos_sa, leo, status )
!      IF_NOTOK_RETURN(status=1)
!
!      ! all members?
!      if ( lekfo_xi ) then
!        ! loop over members:
!        do j = 1, nmodes
!          ! done with output:
!          call Done( leos_xif(j), leo, status )
!          IF_NOTOK_RETURN(status=1)
!          ! done with output:
!          call Done( leos_xia(j), leo, status )
!          IF_NOTOK_RETURN(status=1)
!        end do
!        ! clear:
!        deallocate( leos_xif, stat=status )
!        IF_NOTOK_RETURN(status=1)
!        deallocate( leos_xia, stat=status )
!        IF_NOTOK_RETURN(status=1)
!      end if ! put out members
!      
!#ifdef with_kf_meas_maori
!      ! done with MAORI output stuff:
!      call LE_MAORI_Output_Done( mao_xf, mad, status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Done( mao_xa, mad, status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Done( mao_sf, mad, status )
!      IF_NOTOK_RETURN(status=1)
!      call LE_MAORI_Output_Done( mao_sa, mad, status )
!      IF_NOTOK_RETURN(status=1)
!      ! all members?
!      if ( lekfo_xi ) then
!        ! loop over members:
!        do j = 1, nmodes
!          ! done with output:
!          call LE_MAORI_Output_Done( mao_xif(j), mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! done with output:
!          call LE_MAORI_Output_Done( mao_xia(j), mad, status )
!          IF_NOTOK_RETURN(status=1)
!        end do
!        ! clear:
!        deallocate( mao_xif, stat=status )
!        IF_NOTOK_RETURN(status=1)
!        deallocate( mao_xia, stat=status )
!        IF_NOTOK_RETURN(status=1)
!      end if ! put out members
!#endif

      ! done with output states:
      call lekfo_xf%Done( status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_xa%Done( status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_sf%Done( status )
      IF_NOTOK_RETURN(status=1)
      call lekfo_sa%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! all members?
      if ( lekfo_with_xi ) then
        ! loop over members:
        do j = 1, nmodes
          ! done with output:
          call lekfo_xif(j)%Done( status )
          IF_NOTOK_RETURN(status=1)
          ! done with output:
          call lekfo_xia(j)%Done( status )
          IF_NOTOK_RETURN(status=1)
        end do
        ! clear:
        deallocate( lekfo_xif, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( lekfo_xia, stat=status )
        IF_NOTOK_RETURN(status=1)
      end if ! put out members

    end if  ! ensemble present

    ! done with output of model data:
    call Done( leo, status )
    IF_NOTOK_RETURN(status=1)

!! moved to lekf.F90      
!#ifdef with_kf_meas_maori
!    ! done with MAORI stuff:
!    call LE_MAORI_Data_Done( mad, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
!      
!    ! done with measurement output:
!    call kfom%Done( status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! done with dc output:
!    call Done( kfodc, status )
!    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Output_Done


  ! ***
  

  !subroutine LEKF_Output_Setup( t1, t2, the_end, status )
  subroutine LEKF_Output_Setup( t1, t2, status )

    use GO        , only : goc
    use GO        , only : TDate
    use LE_Output , only : Setup
!#ifdef with_kf_meas_maori
!    use LE_MAORI  , only : LE_MAORI_Data_Setup
!    use LEKF_Data , only : mad
!#endif
    use LEKF_Data , only : leo
!    use LEKF_State, only : Ens_Setup
  
    ! --- in/out --------------------------------
    
    type(TDate), intent(in)       ::  t1, t2
    !logical                       ::  the_end
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_Setup'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
    !! root only ..
    !if ( .not. goc%root) then
    !  status=0; return
    !end if
    
    ! setup output for current time interval if necessary:
    call Setup( leo, t1, t2, status )
    IF_NOTOK_RETURN(status=1)

!! now in lekf.F90
!#ifdef with_kf_meas_maori
!    ! setup maori stuff for current time interval if necessary:
!    call LE_MAORI_Data_Setup( mad, t1, t2, the_end, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
!    
!    ! setup ensemble:
!    call Ens_Setup( status )
!    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LEKF_Output_Setup


  ! ***
  

  subroutine LEKF_Output_PutOut( key, t, status, last )

    use GO              , only : goc
    use GO              , only : TDate
    use LEKF_State      , only : xb, x, sigma, Ens
    !use LEKF_State      , only : Ens_Mean_and_Sigma
    use LEKF_State      , only : kf_with_xb, kf_with_xm
    use LEKF_State      , only : nmodes

!    use LEKF_State      , only : bud0
!    use LEKF_Output_DC  , only : PutOut
    use LE_Output       , only : PutOut
    use LEKF_Data       , only : leo
!    use LEKF_Data       , only : leos_xb, leos_xf, leos_xif, leos_sf, leos_xa, leos_sa, leos_xia
    use LEKF_Data       , only : lekfo_with_xi
!#ifdef with_kf_meas_maori
!    use LE_MAORI        , only : LE_MAORI_Output_Write
!    use LEKF_Data       , only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xif, mao_xa, mao_sa, mao_xia
!#endif

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  key   ! 'forecast', 'analysis'
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status
    logical, intent(in), optional         ::  last
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_PutOut'
    
    ! --- local ---------------------------------
    
    integer     ::  j
    
    ! --- begin ---------------------------------
    
    !! root only ..
    !if ( .not. goc%root) then
    !  status=0; return
    !end if
    
    ! ** data

    ! at this moment ?
    select case ( key )
      case ( 'forecast' )
        ! not yet ..
      case ( 'analysis' )
        ! put out LE data:
        call PutOut( leo, t, status )
        IF_NOTOK_RETURN(status=1)
      case default
        write (gol,'("unuspported key : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
    end select

    ! ** states

    ! background run available ?
    if ( kf_with_xb ) then

      ! at this moment ?
      select case ( key )
        !~ forecast
        case ( 'forecast' )
          ! not yet ..
        !~ analysis
        case ( 'analysis' )
!          ! put out base run;
!          ! also measurements and assimilation flags might be put out ;
!          ! only one budget array, filled by xb run:
!          call PutOut( leos_xb, leo, t, xb%c, xb%cg, bud0, status )
!          IF_NOTOK_RETURN(status=1)
!#ifdef with_kf_meas_maori
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_xb, xb%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!#endif
          ! put out state:
          call lekfo_xb%PutOut( xb, t, status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("unuspported key : ",a)') trim(key); call goErr
          TRACEBACK; status=1; return
      end select

    end if

    ! mean and modes available ?
    ! not necessary to put out data (lon,lat,observations,etc) again ...
    if ( kf_with_xm ) then

      ! ... already done outside (requires collective call!)
      !! fill ensemble mean in 'x' and ensemble standard deviation in 'sigma' :
      !call Ens_Mean_and_Sigma( status )
      !IF_NOTOK_RETURN(status=1)
      ! ...

      ! put out to different files depending on key:
      select case ( key )
        case ( 'forecast' )
!          ! put out mean state; only one budget array, filled by xb run:
!          call PutOut( leos_xf, leo, t, x%c, x%cg, bud0, status, without_data=.true. )
!          IF_NOTOK_RETURN(status=1)
!          ! put out standard deviation:
!          call PutOut( leos_sf, leo, t, sigma%c, sigma%cg, bud0, status, without_data=.true. )
!          IF_NOTOK_RETURN(status=1)
!          ! all members?
!          if ( lekfo_xi ) then
!            ! put out per member:
!            do j = 1, nmodes
!              call PutOut( leos_xif(j), leo, t, Ens(j)%c, Ens(j)%cg, bud0, status, without_data=.true. )
!              IF_NOTOK_RETURN(status=1)
!            end do
!          end if ! put out members
!#ifdef with_kf_meas_maori
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_xf, x%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_sf, sigma%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! all members?
!          if ( lekfo_xi ) then
!            ! put out per member:
!            do j = 1, nmodes
!              call LE_MAORI_Output_Write( mao_xif(j), Ens(j)%mas, mad, status )
!              IF_NOTOK_RETURN(status=1)
!            end do
!          end if ! put out members
!#endif
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_sf, sigma%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! all members?
!          if ( lekfo_xi ) then
!            ! put out per member:
!            do j = 1, nmodes
!              call LE_MAORI_Output_Write( mao_xif(j), Ens(j)%mas, mad, status )
!              IF_NOTOK_RETURN(status=1)
!            end do
!          end if ! put out members
          !
          ! put out ensemble mean:
          call lekfo_xf%PutOut( x, t, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! put out ensemble stdv:
          call lekfo_sf%PutOut( sigma, t, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! all members?
          if ( lekfo_with_xi ) then
            ! put out per member:
            do j = 1, nmodes
              ! put out state:
              call lekfo_xif(j)%PutOut( Ens(j), t, status, without_data=.true. )
              IF_NOTOK_RETURN(status=1)
            end do
          end if ! put out members
        !
        case ( 'analysis' )
!          ! put out mean state; only one budget array, filled by xb run:
!          call PutOut( leos_xa, leo, t, x%c, x%cg, bud0, status, without_data=.true. )
!          IF_NOTOK_RETURN(status=1)
!          ! put out standard deviation:
!          call PutOut( leos_sa, leo, t, sigma%c, sigma%cg, bud0, status, without_data=.true. )
!          IF_NOTOK_RETURN(status=1)
!          ! all members?
!          if ( lekfo_xi ) then
!            ! put out per member:
!            do j = 1, nmodes
!              call PutOut( leos_xia(j), leo, t, Ens(j)%c, Ens(j)%cg, bud0, status, without_data=.true. )
!              IF_NOTOK_RETURN(status=1)
!            end do
!          end if ! put out members
!#ifdef with_kf_meas_maori
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_xa, x%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! write via MAORI:
!          call LE_MAORI_Output_Write( mao_sa, sigma%mas, mad, status )
!          IF_NOTOK_RETURN(status=1)
!          ! all members?
!          if ( lekfo_xi ) then
!            ! put out per member:
!            do j = 1, nmodes
!              call LE_MAORI_Output_Write( mao_xia(j), Ens(j)%mas, mad, status )
!              IF_NOTOK_RETURN(status=1)
!            end do
!          end if ! put out members
!#endif
          !
          ! put out ensemble mean:
          call lekfo_xa%PutOut( x, t, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! put out ensemble stdv:
          call lekfo_sa%PutOut( sigma, t, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! all members?
          if ( lekfo_with_xi ) then
            ! put out per member:
            do j = 1, nmodes
              ! put out state:
              call lekfo_xia(j)%PutOut( Ens(j), t, status, without_data=.true. )
              IF_NOTOK_RETURN(status=1)
            end do
          end if ! put out members
        !
        case default
          write (gol,'("unuspported key : ",a)') trim(key); call goErr
          TRACEBACK; status=1; return
      end select

    end if
    
    ! NOTE: standard deviation sigma has been re-computed now ...
    
!    ! ** noise
!
!    ! put out dc values:
!    call PutOut( kfodc, key, t, status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ** measurements
!
!    ! put out measurements:
!    call kfom%PutOut( key, t, status, last=last )
!    IF_NOTOK_RETURN(status=1)

    ! **

    ! ok
    status = 0
    
  end subroutine LEKF_Output_PutOut


end module LEKF_output
