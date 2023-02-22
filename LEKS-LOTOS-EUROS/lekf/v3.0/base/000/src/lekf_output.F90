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
  
  use LEKF_Output_Meas, only : T_LEKF_Output_Meas
  use LEKF_Output_DC  , only : T_LEKF_Output_DC
  
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  LEKF_Output_Init, LEKF_Output_Done
  public  ::  LEKF_Output_Setup, LEKF_Output_PutOut

  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Output'
  
  
  ! --- var --------------------------------
  
  type(T_LEKF_Output_Meas)   ::  kfom
  type(T_LEKF_Output_DC)     ::  kfodc
  
    
contains


  ! ====================================================
  
  
  subroutine LEKF_Output_Init( rcfile, t0, status )

    use GO              , only : goc
    use GO              , only : TDate
    use GO              , only : TrcFile, Init, Done, ReadRc
    use LE_Output       , only : Init
    use LEKF_Data       , only : lekfo_replace
    use LEKF_Data       , only : leo, leos_xb, leos_xf, leos_sf, leos_xa, leos_sa
    use LEKF_Output_DC  , only : Init
    use LEKF_State      , only : kf_with_xb, kf_with_xm
  
#ifdef with_kf_meas_maori
    use LEKF_Data, only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xa, mao_sa
    !use LE_MAORI , only : LE_MAORI_Data_Init
    use LE_MAORI , only : LE_MAORI_Output_Init
#endif

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  rcfile
    type(TDate), intent(in)               ::  t0
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_Init'
    
    ! --- local ---------------------------------
    
    type(TrcFile)   ::  rcF
    
    ! --- begin ---------------------------------
    
    ! root only ..
    if ( .not. goc%root) then
      status=0; return
    end if
    
    ! settings:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    ! replace existing files ?
    call ReadRc( rcF, 'kf.output.replace', lekfo_replace, status )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
  
    ! setup standard output of model data:
    call Init( leo, rcfile, 'le.output', status )
    IF_NOTOK_RETURN(status=1)

!! now in lekf.f90 ...    
!#ifdef with_kf_meas_maori
!    ! init MAORI stuff:
!    call LE_MAORI_Data_Init( mad, trim(rcfile), t0, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
    
    ! with background run ?
    if ( kf_with_xb ) then

      ! setup standard output of base run concentrations:
      call Init( leos_xb, leo, 'xb', status )
      IF_NOTOK_RETURN(status=1)

      ! setup standard output of base run concentrations:
      call Init( leos_xb, leo, 'xb', status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_maori
      ! init MAORI output stuff:
      call LE_MAORI_Output_Init( mao_xb, mad, 'xb', status )
      IF_NOTOK_RETURN(status=1)
#endif

    end if

    ! with mean and modes ?
    if ( kf_with_xm ) then
    
      ! setup standard output of mean state concentrations:
      call Init( leos_xf, leo, 'xf', status )
      IF_NOTOK_RETURN(status=1)
      call Init( leos_xa, leo, 'xa', status )
      IF_NOTOK_RETURN(status=1)

      ! setup standard output of standard deviation concentrations:
      call Init( leos_sf, leo, 'sf', status )
      IF_NOTOK_RETURN(status=1)
      call Init( leos_sa, leo, 'sa', status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_maori
      ! init MAORI output stuff:
      call LE_MAORI_Output_Init( mao_xf, mad, 'xf', status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Init( mao_xa, mad, 'xa', status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Init( mao_sf, mad, 'sf', status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Init( mao_sa, mad, 'sa', status )
      IF_NOTOK_RETURN(status=1)
#endif

    end if

    ! setup output of measurements:    
    call kfom%Init( rcfile, 'le.output', status )
    IF_NOTOK_RETURN(status=1)

    ! setup output of dc values:    
    call Init( kfodc, rcfile, 'le.output', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Output_Init


  ! ***
  

  subroutine LEKF_Output_Done( status )

    use GO              , only : goc
    use LE_Output       , only : Done
    use LEKF_Data       , only : leo, leos_xb, leos_xf, leos_sf, leos_xa, leos_sa
    use LEKF_Output_DC  , only : Done
    use LEKF_State      , only : kf_with_xb, kf_with_xm
  
#ifdef with_kf_meas_maori
    use LEKF_Data, only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xa, mao_sa
!    use LE_MAORI , only : LE_MAORI_Data_Done
    use LE_MAORI , only : LE_MAORI_Output_Done
#endif
  
    ! --- in/out --------------------------------
    
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_Done'
    
    ! --- begin ---------------------------------
    
    ! root only ..
    if ( .not. goc%root) then
      status=0; return
    end if
    
    if ( kf_with_xb ) then
      ! done with output of base run:
      call Done( leos_xb, leo, status )
      IF_NOTOK_RETURN(status=1)
      
#ifdef with_kf_meas_maori
      ! done with MAORI output stuff:
      call LE_MAORI_Output_Done( mao_xb, mad, status )
      IF_NOTOK_RETURN(status=1)
#endif
    end if

    if ( kf_with_xm ) then
      ! done with output of mean state:
      call Done( leos_xf, leo, status )
      IF_NOTOK_RETURN(status=1)
      call Done( leos_xa, leo, status )
      IF_NOTOK_RETURN(status=1)
 
      ! done with output of standard deviation:
      call Done( leos_sf, leo, status )
      IF_NOTOK_RETURN(status=1)
      call Done( leos_sa, leo, status )
      IF_NOTOK_RETURN(status=1)
      
#ifdef with_kf_meas_maori
      ! done with MAORI output stuff:
      call LE_MAORI_Output_Done( mao_xf, mad, status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Done( mao_xa, mad, status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Done( mao_sf, mad, status )
      IF_NOTOK_RETURN(status=1)
      call LE_MAORI_Output_Done( mao_sa, mad, status )
      IF_NOTOK_RETURN(status=1)
#endif
    end if

    ! done with output of model data:
    call Done( leo, status )
    IF_NOTOK_RETURN(status=1)

!! moved to lekf.F90      
!#ifdef with_kf_meas_maori
!    ! done with MAORI stuff:
!    call LE_MAORI_Data_Done( mad, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
      
    ! done with measurement output:
    call kfom%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with dc output:
    call Done( kfodc, status )
    IF_NOTOK_RETURN(status=1)

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
    
    ! root only ..
    if ( .not. goc%root) then
      status=0; return
    end if
    
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
    use LEKF_State      , only : xb, x, sigma
    !use LEKF_State      , only : Ens_Mean_and_Sigma
    use LEKF_State      , only : kf_with_xb, kf_with_xm

    use LEKF_State      , only : bud0
    use LEKF_Output_DC  , only : PutOut
    use LE_Output       , only : PutOut
    use LEKF_Data       , only : leo, leos_xb, leos_xf, leos_sf, leos_xa, leos_sa
#ifdef with_kf_meas_maori
    use LE_MAORI        , only : LE_MAORI_Output_Write
    use LEKF_Data       , only : mad, mao_xb, mao_xb, mao_xf, mao_sf, mao_xa, mao_sa
#endif

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  key
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status
    logical, intent(in), optional         ::  last
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_PutOut'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
    ! root only ..
    if ( .not. goc%root) then
      status=0; return
    end if
    
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
          ! put out base run;
          ! also measurements and assimilation flags might be put out ;
          ! only one budget array, filled by xb run:
          call PutOut( leos_xb, leo, t, xb%c, xb%cg, bud0, status )
          IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_maori
          ! write via MAORI:
          call LE_MAORI_Output_Write( mao_xb, xb%mas, mad, status )
          IF_NOTOK_RETURN(status=1)
#endif
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
          ! put out mean state; only one budget array, filled by xb run:
          call PutOut( leos_xf, leo, t, x%c, x%cg, bud0, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! put out standard deviation:
          call PutOut( leos_sf, leo, t, sigma%c, sigma%cg, bud0, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_maori
          ! write via MAORI:
          call LE_MAORI_Output_Write( mao_xf, x%mas, mad, status )
          IF_NOTOK_RETURN(status=1)
          ! write via MAORI:
          call LE_MAORI_Output_Write( mao_sf, sigma%mas, mad, status )
          IF_NOTOK_RETURN(status=1)
#endif
        case ( 'analysis' )
          ! put out mean state; only one budget array, filled by xb run:
          call PutOut( leos_xa, leo, t, x%c, x%cg, bud0, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
          ! put out standard deviation:
          call PutOut( leos_sa, leo, t, sigma%c, sigma%cg, bud0, status, without_data=.true. )
          IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_maori
          ! write via MAORI:
          call LE_MAORI_Output_Write( mao_xa, x%mas, mad, status )
          IF_NOTOK_RETURN(status=1)
          ! write via MAORI:
          call LE_MAORI_Output_Write( mao_sa, sigma%mas, mad, status )
          IF_NOTOK_RETURN(status=1)
#endif
        case default
          write (gol,'("unuspported key : ",a)') trim(key); call goErr
          TRACEBACK; status=1; return
      end select

    end if
    
    ! NOTE: standard deviation sigma has been re-computed now ...
    
    ! ** noise
    
    ! put out dc values:
    call PutOut( kfodc, key, t, status )
    IF_NOTOK_RETURN(status=1)
      
    ! ** measurements
    
    ! put out measurements:
    call kfom%PutOut( key, t, status, last=last )
    IF_NOTOK_RETURN(status=1)
      
    ! **

    ! ok
    status = 0
    
  end subroutine LEKF_Output_PutOut


end module LEKF_output
