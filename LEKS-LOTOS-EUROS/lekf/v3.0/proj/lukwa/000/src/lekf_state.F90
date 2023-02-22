!#################################################################
!
! describes the state vector :
!
! c  (nx,ny,nz,nspec)
! AOD(nx,ny,nz_aod)
! dc (nx,ny,nnoise,nhist)
! m  (maxmeas)
! m24(nx,ny,maxcomp24)
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

module LEKF_State

  use GO, only : gol, goPr, goErr

  use Num, only : T_Random

  use NetCDF

  use dims, only : nx, ny, nz, nspec
  use LE_Budget, only : T_Budget

  !use LEKF_dims, only : nnoise, nhist
  use LEKF_noise, only : nnoise, nhist
#ifdef with_kf_meas_modis
  use LEKF_dims, only : modis_maxpix
#endif
#ifdef with_kf_meas_maori
  use LE_MAORI, only : T_MAORI_State
#endif

  implicit none


  ! --- in/out ----------------------------

  private

  public    ::  TState
  public    ::  substate_c
  public    ::  substate_cg
  public    ::  substate_aerh2o
  public    ::  substate_dc
#ifdef with_kf_meas_omi_trc
  public    ::  substate_omi_trc
#endif
#ifdef with_kf_meas_modis
  public    ::  substate_modis
#endif
#ifdef with_kf_meas_maori
  public    ::  substate_maori
#endif
  public    ::  substates

  public    ::  kf_with_xb
  public    ::  xb

  public    ::  kf_with_xm
  public    ::  nmodes_all, nmodes_loc
  public    ::  x, sigma
  public    ::  Ens
  public    ::  imodes, imode_pe

  public    ::  bud0
  public    ::  SendBudgetData

  public    ::  LEKF_State_Init, LEKF_State_Done

  public    ::  Ens_Mean_and_Sigma, Ens_Mean
  public    ::  Ens_Setup

  public    ::  Save_State, Restore_State


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LEKF_State'

  ! parts of the state:
  integer, parameter    ::  substate_c       = 1
  integer, parameter    ::  substate_cg      = 2
  integer, parameter    ::  substate_aerh2o  = 3
  integer, parameter    ::  substate_dc      = 4
#ifdef with_kf_meas_omi_trc
  integer, parameter    ::  substate_omi_trc = 5
#endif
#ifdef with_kf_meas_modis
  integer, parameter    ::  substate_modis   = 6
#endif
#ifdef with_kf_meas_maori
  integer, parameter    ::  substate_maori   = 7
#endif
  ! max:
  integer, parameter    ::  nsubstate = 7
  ! names:
  character(len=8), parameter  ::  substates(nsubstate) = &
                                        (/ 'c       ', &
                                           'cg      ', &
                                           'aerh2o  ', &
                                           'dc      ', &
                                           'omi_trc ', &
                                           'modis   ', &
                                           'maori   ' /)
                                       


  ! --- types --------------------------

  type TState
    !
    ! ~ concentrations
    !
    real, allocatable         ::  c  (:,:,:,:)
    real, allocatable         ::  cg (:,:  ,:)
    real, allocatable         ::  aerh2o(:,:,:)
    !
    ! ~ uncertainty fields
    !
    real, allocatable         ::  dc (:,:,:,:)
    !
    ! ~ simulated observations
    !
#ifdef with_kf_meas_omi_trc
    integer                   ::  n_omi_trc
    real, allocatable         ::  omi_trc(:)
#endif
#ifdef with_kf_meas_modis
    ! o modis aod
    integer                   ::  nmodis
    real, allocatable         ::  modis(:)
#endif
#ifdef with_kf_meas_maori
    ! o maori simulations
    type(T_MAORI_State)       ::  mas
#endif
    !
    ! random number generator:
    type(T_Random)            ::  rnd
    !
  contains
    procedure ::  Init        => state_Init
    procedure ::  Done        => state_Done
    procedure ::  GetValue    => state_GetValue
  end type TState



  ! --- var -----------------------------

  ! flag to enable background run:
  logical                      ::  kf_with_xb
  ! state for base run (root only):
  type(TState), allocatable    ::  xb

  ! single version used for all ensemble members:
  type(T_Budget)               ::  bud0

  ! flag to enable ensemble:
  logical                      ::  kf_with_xm   ! mean, modes
  ! the actual number of modes
  integer                      ::  nmodes_all   ! total number
  integer                      ::  nmodes_loc   ! on this pe
  ! ensemble modes:
  type(TState), allocatable    ::  Ens(:)       ! (nmodes_loc)
  ! global mode number:
  integer, allocatable         ::  imodes(:)    ! (nmodes_loc)
  ! pe holding mode:
  integer, allocatable         ::  imode_pe(:)  ! (nmodes_all)
  ! ensemble mean and stdv:
  type(TState)                 ::  x
  type(TState)                 ::  sigma


contains


  ! ********************************************************************
  ! ***
  ! *** module
  ! ***
  ! ********************************************************************


  subroutine LEKF_State_Init( rcF, status )

    use GO, only : goc
    use GO, only : TrcFile, ReadRc
    use LE_Budget, only : Budget_Init
    use LEKF_dims, only : maxmodes

    ! --- in/out -------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_State_Init'

    ! --- local --------------------------

    integer     ::  j
    integer     ::  id

    ! --- begin ---------------------------

    ! budgets:
    call Budget_Init( bud0, status )
    IF_NOTOK_RETURN(status=1)

    ! with base run ?
    if ( kf_with_xb ) then

      ! storage:
      allocate( xb, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! initialize background state with dummy random seed:
      call xb%Init( 0, status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! with xb

    ! with mean/modes ?
    if ( kf_with_xm ) then

      ! read total number of modes:
      call ReadRc( rcF, 'kf.nmodes', nmodes_all, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      write (gol,'("total number of modes : ",i0)') nmodes_all; call goPr

      ! check ..
      if ( nmodes_all < 2 ) then
        write (gol,'("at least 2 modes needed for an ensemble, requested ",i0)') nmodes_all; call goErr
        TRACEBACK; status=1; return
      end if

      ! init local counter:
      nmodes_loc = 0
      ! storage for global mode number 1,..,nmodes_all per local mode;
      ! allocate maximum space since nmodes_loc is not known yet ...
      allocate( imodes(nmodes_all), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! storage for pe id that holds a mode:
      allocate( imode_pe(nmodes_all), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! assign modes to pe's;
      ! xb will be processed by root:
      id = 0
      ! loop over all modes:
      do j = 1, nmodes_all
        ! assign to next:
        id = modulo( id+1, goc%npes )
        ! info ..
        write (gol,'("  assign mode ",i3," to pe ",i2)') j, id; call goPr
        ! store:
        imode_pe(j) = id
        ! current pe?
        if ( id == goc%myid ) then
          ! increase counter:
          nmodes_loc = nmodes_loc + 1
          ! store global index:
          imodes(nmodes_loc) = j
        end if
      end do
      ! info ...
      write (gol,'("local number of modes : ",i0)') nmodes_loc; call goPr

      ! root has xb, other pe should have at least 1 mode:
      if ( (nmodes_loc < 1) .and. (.not. goc%root) ) then
        write (gol,'("no local number of modes, something wrong in configuration?")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ..
      if ( nmodes_loc > maxmodes ) then
        write (gol,'("nmodes_loc ",i0," exceeds maxmodes ",i0," (not sure if maxmodes is needed?)")') &
                        nmodes_loc, maxmodes; call goErr
        TRACEBACK; status=1; return
      end if

      ! storage for modes:
      allocate( Ens(nmodes_loc), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop:
      do j = 1, nmodes_loc
        ! initialize mode, set user defined seed to mode number:
        call Ens(j)%Init( imodes(j), status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! ensemble mean, no random seed:
      call x%Init( 0, status )
      IF_NOTOK_RETURN(status=1)
      ! standard deviation:
      call sigma%Init( 0, status )
      IF_NOTOK_RETURN(status=1)

    end if  ! with ensemble

    ! ok
    status = 0

  end subroutine LEKF_State_Init


  ! ***


  subroutine LEKF_State_Done( status )

    use GO       , only : goc
    use LE_Budget, only : Budget_Done

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_State_Done'

    ! --- local ----------------------------------

    integer     ::  j

    ! --- begin ----------------------------------

    ! with mean/modes ?
    if ( kf_with_xm ) then

      ! done with ensemble mean:
      call x%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! done with stdv:
      call sigma%Done( status )
      IF_NOTOK_RETURN(status=1)

      ! loop over modes:
      do j = 1, nmodes_loc
        call Ens(j)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( Ens, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( imodes, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( imode_pe, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! with base run ?
    if ( kf_with_xb ) then
      ! done:
      call xb%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( xb, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! with xb

    ! budgets:
    call Budget_Done( bud0, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LEKF_State_Done


  ! ********************************************************************
  ! ***
  ! *** state
  ! ***
  ! ********************************************************************


  subroutine state_Init( x, seed0, status )

#ifdef with_kf_meas_ground
    use LEKF_Dims, only : maxmeas
#endif
#ifdef with_kf_meas_maori
    use LE_MAORI , only : LE_MAORI_State_Init
    use LEKF_Data, only : mad
#endif

    ! --- in/out ---------------------------

    class(TState), intent(out)    ::  x
    integer, intent(in)           ::  seed0
    integer, intent(out)          ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_Init'

    ! --- begin ---------------------------
    
    ! storage:    
    allocate( x%c  (nx,ny,nz,nspec)    , stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( x%cg (nx,ny   ,nspec)    , stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( x%aerh2o(nx,ny,nz)       , stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( x%dc (nx,ny,nnoise,nhist), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! init:
    x%c   = 0.0
    x%cg  = 0.0
    x%aerh2o = 0.0
    x%dc  = 0.0
   !allocate( x%AOD(nx,ny,nz_aod)       ) ;  x%AOD = 0.0

#ifdef with_kf_meas_omi_trc
    ! dummy length, will be allocated dynamically:
    x%n_omi_trc = 0
    allocate( x%omi_trc(1)         , stat=status     )
    IF_NOTOK_RETURN(status=1)
    x%omi_trc = 0.0
#endif
#ifdef with_kf_meas_modis
    ! maximum storage; set current storage to maximum:
    allocate( x%modis(modis_maxpix), stat=status    )
    IF_NOTOK_RETURN(status=1)
    x%modis = 0.0
    x%nmodis = modis_maxpix
#endif
#ifdef with_kf_meas_maori
    ! setup MAORI state:
    call LE_MAORI_State_Init( x%mas, mad, 'lekf', status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! init random number generator:
    call x%rnd%Init(  seed0, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine state_Init


  ! ***


  subroutine state_Done( x, status )

#ifdef with_kf_meas_maori
    use LE_MAORI , only : LE_MAORI_State_Done
    use LEKF_Data, only : mad
#endif

    ! --- in/out ---------------------------

    class(TState), intent(inout)    ::  x
    integer, intent(out)            ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_Done'

    ! --- begin ---------------------------

    deallocate(  x%c , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  x%cg, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  x%aerh2o, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate(  x%dc, stat=status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_omi_trc
    if ( allocated(x%omi_trc) ) then
      deallocate(  x%omi_trc, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
#ifdef with_kf_meas_modis
    deallocate(  x%modis, stat=status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_maori
    ! done with MAORI state:
    call LE_MAORI_State_Done( x%mas, mad, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! done with random number generator:
    call x%rnd%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine state_Done


  ! ***


  subroutine State_GetValue( x, inds, value, status )

    !use Grid   , only : Interpol
    !use LE_Grid, only : lli

    ! --- in/out -----------------------------------

    class(TState), intent(in)     ::  x
    integer, intent(in)           ::  inds(0:4)
    real, intent(out)             ::  value
    integer, intent(out)          ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_GetValue'

    ! --- begin ------------------------------------

    ! select value or interpolate ?
    if ( inds(0) < 1000 ) then
      ! copy value ...
      select case ( inds(0) )
        case ( substate_c       ) ; value = x%c     (inds(1),inds(2),inds(3),inds(4))
        case ( substate_cg      ) ; value = x%cg    (inds(1),inds(2),inds(3))
        case ( substate_aerh2o  ) ; value = x%aerh2o(inds(1),inds(2),inds(3))
        case ( substate_dc      ) ; value = x%dc (inds(1),inds(2),inds(3),inds(4))
#ifdef with_kf_meas_omi_trc
        case ( substate_omi_trc )
          if ( inds(1) <= x%n_omi_trc ) then
            value = x%omi_trc(inds(1))
          else
            stop 'requested invalid element of omi_trc'
          end if
#endif
#ifdef with_kf_meas_modis
        case ( substate_modis   ) ; value = x%modis  (inds(1))
#endif
#ifdef with_kf_meas_maori
        case ( substate_maori   ) !; value = x%maori  (inds(1))
          stop 'state_GetValue not implemented for maori yet ...'
#endif
        case default
          write (gol,'("unsupported inds(0) : ",i6)') inds(0); call goErr
          TRACEBACK; status=1; return
      end select
    else
      stop 'interpol not supported by kf_state/getvalue yet'
      !! interpolate to location in mili-degrees:
      !select case ( inds(0)-1000 )
      !  ! state parts
      !  case ( substate_c   ); call Interpol( lli, x%c  (:,:,inds(3),inds(4)), inds(1)*1000.0, inds(2)*1000.0, value )
      !  !case ( substate_aod ); call Interpol( lli, x%aod(:,:,inds(3)        ), inds(1)*1000.0, inds(2)*1000.0, value )
      !  case ( substate_dc  ); call Interpol( lli, x%dc (:,:,inds(3),inds(4)), inds(1)*1000.0, inds(2)*1000.0, value )
      !  !case ( substate_m24 ); call Interpol( lli, x%m24(:,:,inds(3)        ), inds(1)*1000.0, inds(2)*1000.0, value )
      !  ! error ...
      !  case default
      !    write (gol,'("unsupported inds(0) : ",i6)') inds(0); call goErr
      !    TRACEBACK; status=1; return
      !end select
    end if


    ! ok
    status = 0

  end  subroutine State_GetValue


  ! ***


  subroutine Save_State( x, t, path, key, status )

    use GO, only : TDate
    use LE_Restart, only : T_LE_Restart_File
    use LE_Restart, only : LE_Restart_Create, LE_Restart_Write, LE_Restart_Close

    ! --- in/out ----------------------------------

    type(TState)                    ::  x
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Save_State'

    ! --- local -----------------------------------

    type(T_LE_Restart_File)   ::  F
    integer       ::  dimid_nnoise, dimid_nhist
    integer       ::  varid_dc
#ifdef with_kf_meas_omi_trc
    integer       ::  dimid_omi_trc
    integer       ::  varid_omi_trc
#endif
#ifdef with_kf_meas_modis
    integer       ::  dimid_modis
    integer       ::  varid_modis
#endif

    ! --- begin -----------------------------------

    ! create file, define standard dimensions:
    call LE_Restart_Create( F, t, path, key, status )
    IF_NOTOK_RETURN(status=1)

    ! create extra dimensions:
    status = nf90_def_dim( F%ncid, 'nnoise', nnoise, dimid_nnoise )
    IF_NF90_NOTOK_RETURN(status=1)
    status = nf90_def_dim( F%ncid, 'nhist', nhist, dimid_nhist )
    IF_NF90_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_omi_trc
    ! always write omi_trc since it's needed for restart from
    ! restart files that were produced when OMI data was missing ;
    ! use minimum length to avoid errors on to many unlimitted dimensions ..
    status = nf90_def_dim( F%ncid, 'omi_trc', max(1,x%n_omi_trc), dimid_omi_trc )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_modis
    if ( x%nmodis > 0 ) then
      status = nf90_def_dim( F%ncid, 'modis', x%nmodis, dimid_modis )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
#endif

    ! create extra variables:
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Var( F%ncid, 'dc', NF90_FLOAT, &
                           (/F%dimid_nx,F%dimid_ny,dimid_nnoise,dimid_nhist/), &
                           varid_dc )
    IF_NF90_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_omi_trc
    ! always write omi_trc since it's needed for restart from
    ! restart files that were produced when OMI data was missing:
    status = NF90_Def_Var( F%ncid, 'omi_trc', NF90_FLOAT, &
                           (/dimid_omi_trc/), &
                           varid_omi_trc )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_modis
    if ( x%nmodis > 0 ) then
      status = NF90_Def_Var( F%ncid, 'modis', NF90_FLOAT, &
                             (/dimid_modis/), &
                             varid_modis )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
#endif

    ! end definition, write standard fields:
    call LE_Restart_Write( F, x%c, x%cg, x%aerh2o, bud0, status )
    IF_NOTOK_RETURN(status=1)

    ! write data:
    status = NF90_Put_Var( F%ncid, varid_dc , x%dc  )
    IF_NF90_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_omi_trc
    if ( x%n_omi_trc > 0 ) then
      status = NF90_Put_Var( F%ncid, varid_omi_trc, x%omi_trc )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
#endif
#ifdef with_kf_meas_modis
    if ( x%nmodis > 0 ) then
      status = NF90_Put_Var( F%ncid, varid_modis, x%modis(1:x%nmodis) )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
#endif
    ! MAORI simulations should be recomputed ...

    ! close file:
    call LE_Restart_Close( F, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Save_State


  ! ***


  subroutine Restore_State( x, t, path, key, status )

    use GO, only : TDate
    use LE_Restart, only : LE_Restart_Restore_Data, LE_Restart_Restore_State
    use LE_Restart, only : LE_Restart_Restore

    ! --- in/out ----------------------------------

    type(TState), intent(inout)     ::  x
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Restore_State'

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

    ! restore data:
    call LE_Restart_Restore_Data( t, path, key, status )
    IF_NOTOK_RETURN(status=1)

    ! restore concentration field(s):
    call LE_Restart_Restore_State( x%c, x%cg, x%aerh2o, bud0, t, path, key, status )
    IF_NOTOK_RETURN(status=1)

    ! restore data:
    call LE_Restart_Restore( 'dc', x%dc, t, path, key, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_omi_trc
    ! restore data:
    call LE_Restart_Restore( 'omi_trc', x%omi_trc, t, path, key, status )
    ! adhoc fix, in new runs size is increased ..
    if ( status /= 0 ) then
      ! first fill all with zeros:
      x%omi_trc = 0.0
      ! try to read the first part from the file:
      call LE_Restart_Restore( 'omi_trc', x%omi_trc(1:20000), t, path, key, status )
      IF_NOTOK_RETURN(status=1)
    endif
#endif

#ifdef with_kf_meas_modis
    ! restore data:
    call LE_Restart_Restore( 'modis', x%modis, t, path, key, status, cnt=x%nmodis )
    if (status<0) then
      write (gol,'("    WARNING - no modis data found to restore; continue ...")'); call goPr
      x%nmodis = 0   ! no data was found ...
    else
      TRACEBACK; status=1; return
    end if
#endif

#ifdef with_kf_meas_maori
    !stop 'Restore_State not implemented for maori; recompute from state'
    write (gol,'("    WARNING - restore not implemented for MAORI yet (not necessary); rcontinue ...")'); call goPr
#endif

    ! ok
    status = 0

  end subroutine Restore_State



  ! ********************************************************************
  ! ***
  ! *** modes
  ! ***
  ! ********************************************************************




  subroutine Ens_Mean( status )

    use GO       , only : goc
#ifdef with_kf_meas_maori
    use MAORI    , only : MAORI_Data_Inquire
    use MAORI    , only : MAORI_State_Values_Get, MAORI_State_Values_Put
    use LEKF_Data, only : mad
#endif

    ! --- in/out -------------------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------------------

    character(len=*), parameter   ::  rname = mname//'/Ens_Mean'

    ! --- local --------------------------------------------

    type(TState)            ::  lsum
    integer                 ::  j

#ifdef with_kf_meas_maori
    integer                 ::  nvalue
    real, allocatable       ::  lsum_values(:)
    real, allocatable       ::  ensj_values(:)
    real, allocatable       ::     x_values(:)
#endif

    ! --- begin ---------------------------------------------

    ! ~ storage

    ! init local sum, not random seed needed:
    call lsum%Init( 0, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_maori
    ! get number of values in maori state:
    call MAORI_Data_Inquire( mad, status, nvalue=nvalue )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    if ( nvalue > 0 ) then
      ! storage:
      allocate(    x_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( ensj_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lsum_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ init sum

    ! set elements to zero:
    lsum%c       = 0.0
    lsum%cg      = 0.0
    lsum%aerh2o  = 0.0
    lsum%dc      = 0.0
#ifdef with_kf_meas_omi_trc
    if ( lsum%n_omi_trc > 0 ) then
      lsum%omi_trc = 0.0
    end if
#endif
#ifdef with_kf_meas_modis
    lsum%modis   = 0.0
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      lsum_values = 0.0
    end if
#endif

    ! ~ sum local modes

    ! loop over modes:
    do j = 1, nmodes_loc
      ! add contribution of current mode:
      lsum%c      = lsum%c        + Ens(j)%c
      lsum%cg     = lsum%cg       + Ens(j)%cg
      lsum%aerh2o = lsum%aerh2o   + Ens(j)%aerh2o
      lsum%dc     = lsum%dc       + Ens(j)%dc
#ifdef with_kf_meas_omi_trc
      if ( lsum%n_omi_trc > 0 ) then
        lsum%omi_trc = lsum%omi_trc + Ens(j)%omi_trc
      end if
#endif
#ifdef with_kf_meas_modis
      lsum%modis   = lsum%modis   + Ens(j)%modis
#endif
#ifdef with_kf_meas_maori
      if ( nvalue > 0 ) then
        ! extract values:
        call MAORI_State_Values_Get( Ens(j)%mas, mad, ensj_values, status )
        IF_NOTOK_RETURN(status=1)
        ! add contribution:
        lsum_values = lsum_values + ensj_values
      end if
#endif
    end do  ! local modes

    ! ~ global sum
    
    ! add all local sums, send result to all pe:
    call goc%AllReduce( lsum%c     , x%c     , status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%cg    , x%cg    , status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%aerh2o, x%aerh2o, status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%dc    , x%dc    , status )
    IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_omi_trc
    if ( lsum%n_omi_trc > 0 ) then
      call goc%AllReduce( lsum%omi_trc, x%omi_trc, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
#ifdef with_kf_meas_modis
    call goc%AllReduce( lsum%modis, x%modis, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      call goc%AllReduce( lsum_values, x_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ ensemble mean
    
    ! aver:
    x%c      = x%c      / nmodes_all
    x%cg     = x%cg     / nmodes_all
    x%aerh2o = x%aerh2o / nmodes_all
    x%dc     = x%dc     / nmodes_all
#ifdef with_kf_meas_omi_trc
    if ( lsum%n_omi_trc > 0 ) then
      x%omi_trc = x%omi_trc / nmodes_all
    end if
#endif
#ifdef with_kf_meas_modis
    x%modis   = x%modis / nmodes_all
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      ! average:
      x_values = x_values / nmodes_all
      ! flush array into original structure:
      call MAORI_State_Values_Put( x%mas, mad, x_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ clear

#ifdef with_kf_meas_maori
    ! defined?
    if ( nvalue > 0 ) then
      ! clear:
      deallocate(    x_values, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( ensj_values, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lsum_values, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
    
    ! clear:
    call lsum%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Ens_Mean


  ! ***


  subroutine Ens_Mean_and_Sigma( status )

    use GO       , only : goc
#ifdef with_kf_meas_maori
    use MAORI    , only : MAORI_Data_Inquire
    use MAORI    , only : MAORI_State_Values_Get, MAORI_State_Values_Put
    use LEKF_Data, only : mad
#endif

    ! --- in/out -------------------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------------------

    character(len=*), parameter   ::  rname = mname//'/Ens_Mean_and_Sigma'

    ! --- local --------------------------------------------

    type(TState)            ::  lsum
    integer                 ::  j

#ifdef with_kf_meas_maori
    integer                 ::  nvalue
    real, allocatable       ::   lsum_values(:)
    real, allocatable       ::   ensj_values(:)
    real, allocatable       ::      x_values(:)
    real, allocatable       ::  sigma_values(:)
#endif

    ! --- begin ---------------------------------------------

    ! info ..
    write (gol,*) trim(rname)//': begin'; call goPr
    write (gol,*) trim(rname)//': mean ...'; call goPr

    ! ensure that ensemble mean is updated:
    call Ens_Mean( status )
    IF_NOTOK_RETURN(status=1)

    ! info ..
    write (gol,*) trim(rname)//': sigma ...'; call goPr

    ! ~ storage

    ! info ..
    write (gol,*) trim(rname)//':   init lsum ...'; call goPr

    ! init local sum, not random seed needed:
    call lsum%Init( 0, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_maori
    ! get number of values in maori state:
    call MAORI_Data_Inquire( mad, status, nvalue=nvalue )
    IF_NOTOK_RETURN(status=1)
    ! any data ?
    if ( nvalue > 0 ) then
      ! storage:
      allocate(    x_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( ensj_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lsum_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( sigma_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! init local sum:
      lsum_values = 0.0
      ! extract from mean:
      call MAORI_State_Values_Get( x%mas, mad, x_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ init sum

    ! info ..
    write (gol,*) trim(rname)//':   fill zero ...'; call goPr

    ! init sigma to zero:
    lsum%c       = 0.0
    lsum%cg      = 0.0
    lsum%aerh2o  = 0.0
    lsum%dc      = 0.0
#ifdef with_kf_meas_omi_trc
    if ( x%n_omi_trc > 0 ) then
      lsum%omi_trc = 0.0
    end if
#endif
#ifdef with_kf_meas_modis
    lsum%modis   = 0.0
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      lsum_values = 0.0
    end if
#endif

    ! ~ sum local modes

    ! info ..
    write (gol,*) trim(rname)//':   local sums ...'; call goPr

    ! loop over modes:
    do j = 1, nmodes_loc
      !  lsum := lsum + ( Ens(j) - x )**2
      lsum%c      = lsum%c        + ( Ens(j)%c       - x%c      )**2
      lsum%cg     = lsum%cg       + ( Ens(j)%cg      - x%cg     )**2
      lsum%aerh2o = lsum%aerh2o   + ( Ens(j)%aerh2o  - x%aerh2o )**2
      lsum%dc     = lsum%dc       + ( Ens(j)%dc      - x%dc     )**2
#ifdef with_kf_meas_omi_trc
      if ( x%n_omi_trc > 0 ) then
        lsum%omi_trc = lsum%omi_trc + ( Ens(j)%omi_trc - x%omi_trc)**2
      end if
#endif
#ifdef with_kf_meas_modis
      lsum%modis   = lsum%modis   + ( Ens(j)%modis   - x%modis  )**2
#endif
#ifdef with_kf_meas_maori
      if ( nvalue > 0 ) then
        ! extract values:
        call MAORI_State_Values_Get( Ens(j)%mas, mad, ensj_values, status )
        IF_NOTOK_RETURN(status=1)
        ! add contribution:
        lsum_values   = lsum_values   + ( ensj_values   - x_values  )**2
      end if
#endif

    end do  ! local modes

    ! ~ global sum
    
    ! info ..
    write (gol,*) trim(rname)//':   global sums ...'; call goPr

    ! add all local sums, send result to all pe:
    call goc%AllReduce( lsum%c     , sigma%c     , status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%cg    , sigma%cg    , status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%aerh2o, sigma%aerh2o, status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( lsum%dc    , sigma%dc    , status )
    IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_omi_trc
    if ( lsum%n_omi_trc > 0 ) then
      call goc%AllReduce( lsum%omi_trc, sigma%omi_trc, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
#ifdef with_kf_meas_modis
    call goc%AllReduce( lsum%modis, sigma%modis, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      call goc%AllReduce( lsum_values, sigma_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ ensemble mean, post:
    
    ! aver:
    sigma%c      = sqrt( sigma%c      / (nmodes_all-1) )
    sigma%cg     = sqrt( sigma%cg     / (nmodes_all-1) )
    sigma%aerh2o = sqrt( sigma%aerh2o / (nmodes_all-1) )
    sigma%dc     = sqrt( sigma%dc     / (nmodes_all-1) )
#ifdef with_kf_meas_omi_trc
    if ( lsum%n_omi_trc > 0 ) then
      sigma%omi_trc = sqrt( sigma%omi_trc / (nmodes_all-1) )
    end if
#endif
#ifdef with_kf_meas_modis
    sigma%modis   = sqrt( sigma%modis / (nmodes_all-1) )
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      ! average:
      sigma_values = sqrt( sigma_values / (nmodes_all-1) )
      ! flush array into original structure:
      call MAORI_State_Values_Put( sigma%mas, mad, sigma_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ clear

#ifdef with_kf_meas_maori
    ! defined?
    if ( nvalue > 0 ) then
      ! clear:
      deallocate(    x_values, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( ensj_values, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lsum_values, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( sigma_values, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
    
    ! clear:
    call lsum%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! info ..
    write (gol,*) trim(rname)//': end'; call goPr

    ! ok
    status = 0

  end subroutine Ens_Mean_and_Sigma


  ! ***


  subroutine Ens_Setup( status )

#ifdef with_kf_meas_maori
    use LE_MAORI , only : LE_MAORI_State_Setup
    use LEKF_Data, only : mad
#endif

    ! --- in/out -------------------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------------------

    character(len=*), parameter   ::  rname = mname//'/Ens_Setup'

    ! --- local --------------------------------------------

    integer   ::  j

    ! --- begin ---------------------------------------------

#ifdef with_kf_meas_maori
    ! with model run?
    if ( kf_with_xb ) then
      ! setup MAORI state for current time interval:
      call LE_MAORI_State_Setup( xb%mas, mad, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! with mean/modes ?
    if ( kf_with_xm ) then
      ! loop over modes:
      do j = 1, nmodes_loc
#ifdef with_kf_meas_maori
        ! setup MAORI state for current time interval:
        call LE_MAORI_State_Setup( Ens(j)%mas, mad, status )
        IF_NOTOK_RETURN(status=1)
#endif
      end do
    end if  ! with mean,modes

    ! ok
    status = 0

  end subroutine Ens_Setup


  ! ********************************************************************
  ! ***
  ! *** budget data
  ! ***
  ! ********************************************************************
  
  
  !
  ! Scatter selected budgets from root (where xb is propagated)
  ! to all other pe:
  !
  !   bud0%drydepos%cnh3_ave_prev(:,:)
  !                 cso2_ave_prev(:,:)
  ! 
  ! The drydpos part is used in "LEKF_Driver_Init_State" in the
  ! the call that updates deposition velocities:
  !
  !   ! re-compute deposition velocities if these depend on concentrations;
  !   ! also the adhoc factor o3fac is applied if necessary ;
  !   ! use overal budget as input, used for NH3 compensation point:
  !   call LE_DryDepos_Setup_vd( x%c, bud0%drydepos, .true., t, status )
  !   IF_NOTOK_RETURN(status=1)
  !

  subroutine SendBudgetData( status )
  
    use GO, only : goc

    ! --- in/out -------------------------------------------

    integer, intent(out)    ::  status

    ! --- const --------------------------------------------

    character(len=*), parameter   ::  rname = mname//'/SendBudgetData'

    ! --- local --------------------------------------------

    ! --- begin ---------------------------------------------
    
    ! broadcast from master pe to other pe's:
    call goc%BCast( 0, bud0%drydepos%cnh3_ave_prev, status )
    IF_NOTOK_RETURN(status=1)
    call goc%BCast( 0, bud0%drydepos%cso2_ave_prev, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SendBudgetData
  
  

end module LEKF_State
