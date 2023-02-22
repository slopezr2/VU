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
#ifdef with_kf_meas_sat
  use LE_Output, only : T_LE_Output_Sat_State
#endif

  implicit none


  ! --- in/out ----------------------------

  private

  public    ::  TState
  public    ::  substate_c
  public    ::  substate_cg
  public    ::  substate_aerh2o
  public    ::  substate_dc
#ifdef with_kf_meas_sat
  public    ::  substate_sat
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
  !public    ::  nmodes_all, nmodes_loc
  public    ::  nmodes
  public    ::  x, sigma
  public    ::  Ens
  !public    ::  imodes, imode_pe

  public    ::  bud0
  !public    ::  SendBudgetData

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
#ifdef with_kf_meas_modis
  integer, parameter    ::  substate_modis   = 5
#endif
#ifdef with_kf_meas_maori
  integer, parameter    ::  substate_maori   = 6
#endif
#ifdef with_kf_meas_sat
  integer, parameter    ::  substate_sat     = 7
#endif
  ! max:
  integer, parameter    ::  nsubstate = 7
  ! names:
  character(len=8), parameter  ::  substates(nsubstate) = &
                                        (/ 'c       ', &
                                           'cg      ', &
                                           'aerh2o  ', &
                                           'dc      ', &
                                           'modis   ', &
                                           'maori   ', &
                                           'sat     ' /)
                                       


  ! --- types --------------------------

  type TState
    ! name: 'xb', 'xa', 'ensm0001', ..
    character(len=32)         ::  name
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
#ifdef with_kf_meas_sat
    ! number of satellite sets in LE output:
    integer                                     ::  nsat
    ! output state per satellite set
    ! (simulated retrieval, and eventually intermediate concentration profiles);
    ! note that LEKF output simulates retrievals from "c" array!
    type(T_LE_Output_Sat_State), allocatable    ::  sat_state(:)  ! (nsat)
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
  !integer                      ::  nmodes_all   ! total number
  !integer                      ::  nmodes_loc   ! on this pe
  integer                      ::  nmodes
  ! ensemble modes:
  type(TState), allocatable    ::  Ens(:)       ! (nmodes[_loc])
  !! global mode number:
  !integer, allocatable         ::  imodes(:)    ! (nmodes_loc)
  !! pe holding mode:
  !integer, allocatable         ::  imode_pe(:)  ! (nmodes_all)
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
    use GO, only : TrcFile
    use LE_Budget, only : Budget_Init
    use LEKF_dims, only : maxmodes

    ! --- in/out -------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_State_Init'

    ! --- local --------------------------

    integer               ::  j
    integer               ::  id
    character(len=32)     ::  name

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
      call xb%Init( 'xb', 0, status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! with xb

    ! with mean/modes ?
    if ( kf_with_xm ) then

      ! read total number of modes:
      call rcF%Get( 'kf.nmodes', nmodes, status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      write (gol,'("total number of modes : ",i0)') nmodes; call goPr

      ! check ..
      if ( nmodes < 2 ) then
        write (gol,'("at least 2 modes needed for an ensemble, requested ",i0)') nmodes; call goErr
        TRACEBACK; status=1; return
      end if

!      ! init local counter:
!      nmodes_loc = 0
!      ! storage for global mode number 1,..,nmodes_all per local mode;
!      ! allocate maximum space since nmodes_loc is not known yet ...
!      allocate( imodes(nmodes_all), stat=status )
!      IF_NOTOK_RETURN(status=1)
!      ! storage for pe id that holds a mode:
!      allocate( imode_pe(nmodes_all), stat=status )
!      IF_NOTOK_RETURN(status=1)
!      ! assign modes to pe's;
!      ! xb will be processed by root:
!      id = 0
!      ! loop over all modes:
!      do j = 1, nmodes_all
!        ! assign to next:
!        id = modulo( id+1, goc%npes )
!        ! info ..
!        write (gol,'("  assign mode ",i3," to pe ",i2)') j, id; call goPr
!        ! store:
!        imode_pe(j) = id
!        ! current pe?
!        if ( id == goc%myid ) then
!          ! increase counter:
!          nmodes_loc = nmodes_loc + 1
!          ! store global index:
!          imodes(nmodes_loc) = j
!        end if
!      end do
!      ! info ...
!      write (gol,'("local number of modes : ",i0)') nmodes_loc; call goPr
!
!      ! root has xb, other pe should have at least 1 mode:
!      if ( (nmodes_loc < 1) .and. (.not. goc%root) ) then
!        write (gol,'("no local number of modes, something wrong in configuration?")'); call goErr
!        TRACEBACK; status=1; return
!      end if
!      ! check ..
!      if ( nmodes_loc > maxmodes ) then
!        write (gol,'("nmodes_loc ",i0," exceeds maxmodes ",i0," (not sure if maxmodes is needed?)")') &
!                        nmodes_loc, maxmodes; call goErr
!        TRACEBACK; status=1; return
!      end if

      ! storage for modes:
      allocate( Ens(nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop:
      do j = 1, nmodes
        ! initialize mode, set user defined seed to mode number:
        !call Ens(j)%Init( imodes(j), status )
        write (name,'("ensm",i2.2)') j
        call Ens(j)%Init( trim(name), j, status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! ensemble mean, no random seed:
      call x%Init( 'mu', 0, status )
      IF_NOTOK_RETURN(status=1)
      ! standard deviation:
      call sigma%Init( 'sigma', 0, status )
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
      do j = 1, nmodes
        call Ens(j)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( Ens, stat=status )
      IF_NOTOK_RETURN(status=1)
      !! clear:
      !deallocate( imodes, stat=status )
      !IF_NOTOK_RETURN(status=1)
      !! clear:
      !deallocate( imode_pe, stat=status )
      !IF_NOTOK_RETURN(status=1)

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


  subroutine state_Init( self, name, seed0, status )

#ifdef with_kf_meas_ground
    use LEKF_Dims, only : maxmeas
#endif
#ifdef with_kf_meas_maori
    use LE_MAORI , only : LE_MAORI_State_Init
    use LEKF_Data, only : mad
#endif
#ifdef with_kf_meas_sat
    use LEKF_Data, only : leo
#endif

    ! --- in/out ---------------------------

    class(TState), intent(out)    ::  self
    character(len=*), intent(in)  ::  name
    integer, intent(in)           ::  seed0
    integer, intent(out)          ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_Init'

    ! --- local ---------------------------
    
    integer       ::  k

    ! --- begin ---------------------------
    
    ! store:
    self%name = trim(name)
    
    ! storage:    
    allocate( self%c  (nx,ny,nz,nspec)    , source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%cg (nx,ny   ,nspec)    , source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%aerh2o(nx,ny,nz)       , source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%dc (nx,ny,nnoise,nhist), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_sat
    ! any sat output ?
    if ( leo%nsat > 0 ) then
      ! store number:
      self%nsat = leo%nsat
      ! storage for output stuff:
      allocate( self%sat_state(self%nsat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over files to be written:
      do k = 1, self%nsat
        ! setup state (see "le_output.F90" for original):
        call self%sat_state(k)%Init( status, key=trim(self%name), &
                                       description='simulated '//trim(leo%file_sat(k))//' retrievals' )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
    else
      ! no sat output:
      leo%nsat = 0
    end if  ! any output?
#endif
#ifdef with_kf_meas_modis
    ! maximum storage; set current storage to maximum:
    allocate( self%modis(modis_maxpix), stat=status    )
    IF_NOTOK_RETURN(status=1)
    self%modis = 0.0
    self%nmodis = modis_maxpix
#endif
#ifdef with_kf_meas_maori
    ! setup MAORI state:
    call LE_MAORI_State_Init( self%mas, mad, 'lekf', status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! init random number generator:
    call self%rnd%Init(  seed0, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine state_Init


  ! ***


  subroutine state_Done( self, status )

#ifdef with_kf_meas_maori
    use LE_MAORI , only : LE_MAORI_State_Done
    use LEKF_Data, only : mad
#endif

    ! --- in/out ---------------------------

    class(TState), intent(inout)    ::  self
    integer, intent(out)            ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_Done'
    
    ! --- local ---------------------------

    integer       ::  k

    ! --- begin ---------------------------

    deallocate(  self%c , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  self%cg, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  self%aerh2o, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate(  self%dc, stat=status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_sat
    ! any sat output ?
    if ( self%nsat > 0 ) then
      ! loop over files to be written:
      do k = 1, self%nsat
        ! done with state:
        call self%sat_state(k)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
      ! clear:
      deallocate( self%sat_state, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if  ! any output?
#endif
#ifdef with_kf_meas_modis
    deallocate(  self%modis, stat=status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_maori
    ! done with MAORI state:
    call LE_MAORI_State_Done( self%mas, mad, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! done with random number generator:
    call self%rnd%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine state_Done


  ! ***


  subroutine State_GetValue( self, inds, value, status )

    !use Grid   , only : Interpol
    !use LE_Grid, only : lli

    ! --- in/out -----------------------------------

    class(TState), intent(in)     ::  self
    integer, intent(in)           ::  inds(0:4)
    real, intent(out)             ::  value
    integer, intent(out)          ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/state_GetValue'

    ! --- local ------------------------------------
    
    integer     ::  isat
    integer     ::  ipix

    ! --- begin ------------------------------------

    ! select value or interpolate ?
    if ( inds(0) < 1000 ) then
      ! copy value ...
      select case ( inds(0) )
        case ( substate_c       ) ; value = self%c     (inds(1),inds(2),inds(3),inds(4))
        case ( substate_cg      ) ; value = self%cg    (inds(1),inds(2),inds(3))
        case ( substate_aerh2o  ) ; value = self%aerh2o(inds(1),inds(2),inds(3))
        case ( substate_dc      ) ; value = self%dc (inds(1),inds(2),inds(3),inds(4))

#ifdef with_kf_meas_sat
        ! satellite simulations:
        !   inds(0:3) = (/substate_sat,isat,ipix/)
        case ( substate_sat )
          ! extract:
          isat = inds(1)
          ipix = inds(2)
          ! check ...
          if ( (isat < 1) .or. (isat > self%nsat) ) then
            write (gol,'("requested isat=",i0," while nsat=",i0)') isat, self%nsat; call goErr
            TRACEBACK; status=1; return
          end if
          ! check ...
          if ( (ipix < 1) .or. (ipix > self%sat_state(isat)%npix) ) then
            write (gol,'("requested ipix=",i0," while npix=",i0," for isat=",i0)') &
                        ipix, self%sat_state(isat)%npix, isat; call goErr
            TRACEBACK; status=1; return
          end if
          ! get value:
          call self%sat_state(isat)%GetPixel( ipix, status, y=value )
          IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_kf_meas_modis
        case ( substate_modis   ) ; value = self%modis  (inds(1))
#endif
#ifdef with_kf_meas_maori
        case ( substate_maori   ) !; value = self%maori  (inds(1))
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
      !  case ( substate_c   ); call Interpol( lli, self%c  (:,:,inds(3),inds(4)), inds(1)*1000.0, inds(2)*1000.0, value )
      !  !case ( substate_aod ); call Interpol( lli, self%aod(:,:,inds(3)        ), inds(1)*1000.0, inds(2)*1000.0, value )
      !  case ( substate_dc  ); call Interpol( lli, self%dc (:,:,inds(3),inds(4)), inds(1)*1000.0, inds(2)*1000.0, value )
      !  !case ( substate_m24 ); call Interpol( lli, self%m24(:,:,inds(3)        ), inds(1)*1000.0, inds(2)*1000.0, value )
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

    use GO        , only : TDate
    use GO        , only : goc
    use LE_Grid   , only : dom
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
#ifdef with_kf_meas_modis
    integer       ::  dimid_modis
    integer       ::  varid_modis
#endif

    ! --- begin -----------------------------------

    ! create file, define standard dimensions:
    call LE_Restart_Create( F, t, path, key, status )
    IF_NOTOK_RETURN(status=1)

    ! root only:
    if ( goc%root ) then

      ! create extra dimensions:
      status = nf90_def_dim( F%ncid, 'nnoise', nnoise, dimid_nnoise )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_def_dim( F%ncid, 'nhist', nhist, dimid_nhist )
      IF_NF90_NOTOK_RETURN(status=1)
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
      
#ifdef with_kf_meas_modis
      if ( x%nmodis > 0 ) then
        status = NF90_Def_Var( F%ncid, 'modis', NF90_FLOAT, &
                               (/dimid_modis/), &
                               varid_modis )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
#endif

    end if  ! root

    ! end definition, write standard fields:
    call LE_Restart_Write( F, x%c, x%cg, x%aerh2o, bud0, status )
    IF_NOTOK_RETURN(status=1)

    ! write data:
    call dom%Put_Var( F%ncid, varid_dc, -1, x%dc, status )
    IF_NOTOK_RETURN(status=1)

    ! no parallel write yet, just write whatever is present on root ...
    if ( goc%root ) then

#ifdef with_kf_meas_modis
      if ( x%nmodis > 0 ) then
        ! write data:
        status = NF90_Put_Var( F%ncid, varid_modis, x%modis(1:x%nmodis) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
#endif

    end if ! root

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

!    type(TState)            ::  lsum
    integer                 ::  j

#ifdef with_kf_meas_maori
    integer                 ::  nvalue
!    real, allocatable       ::  lsum_values(:)
    real, allocatable       ::  ensj_values(:)
    real, allocatable       ::     x_values(:)
#endif

    ! --- begin ---------------------------------------------

    ! ~ storage

!    ! init local sum, not random seed needed:
!    call lsum%Init( 0, status )
!    IF_NOTOK_RETURN(status=1)

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
      !allocate( lsum_values(nvalue), stat=status )
      !IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ init sum

    ! set elements to zero:
    x%c       = 0.0
    x%cg      = 0.0
    x%aerh2o  = 0.0
    x%dc      = 0.0
#ifdef with_kf_meas_modis
    x%modis   = 0.0
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      x_values = 0.0
    end if
#endif

    ! ~ sum modes

    ! loop over modes:
    do j = 1, nmodes
      ! add contribution of current mode:
      x%c      = x%c        + Ens(j)%c
      x%cg     = x%cg       + Ens(j)%cg
      x%aerh2o = x%aerh2o   + Ens(j)%aerh2o
      x%dc     = x%dc       + Ens(j)%dc
#ifdef with_kf_meas_modis
      x%modis   = x%modis   + Ens(j)%modis
#endif
#ifdef with_kf_meas_maori
      if ( nvalue > 0 ) then
        ! extract values:
        call MAORI_State_Values_Get( Ens(j)%mas, mad, ensj_values, status )
        IF_NOTOK_RETURN(status=1)
        ! add contribution:
        x_values = x_values + ensj_values
      end if
#endif
    end do  ! local modes

!    ! ~ global sum
!    
!    ! add all local sums, send result to all pe:
!    call goc%AllReduce( lsum%c     , x%c     , status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%cg    , x%cg    , status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%aerh2o, x%aerh2o, status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%dc    , x%dc    , status )
!    IF_NOTOK_RETURN(status=1)
!#ifdef with_kf_meas_modis
!    call goc%AllReduce( lsum%modis, x%modis, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
!#ifdef with_kf_meas_maori
!    if ( nvalue > 0 ) then
!      call goc%AllReduce( lsum_values, x_values, status )
!      IF_NOTOK_RETURN(status=1)
!    end if
!#endif

    ! ~ ensemble mean
    
    ! aver:
    x%c      = x%c      / nmodes
    x%cg     = x%cg     / nmodes
    x%aerh2o = x%aerh2o / nmodes
    x%dc     = x%dc     / nmodes
#ifdef with_kf_meas_modis
    x%modis   = x%modis / nmodes
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      ! average:
      x_values = x_values / nmodes
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
!      deallocate( lsum_values, stat=status )
!      IF_NOTOK_RETURN(status=1)
    end if
#endif
    
!    ! clear:
!    call lsum%Done( status )
!    IF_NOTOK_RETURN(status=1)

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

!    type(TState)            ::  lsum
    integer                 ::  j
    
          integer   ::  qi

#ifdef with_kf_meas_maori
    integer                 ::  nvalue
!    real, allocatable       ::   lsum_values(:)
    real, allocatable       ::   ensj_values(:)
    real, allocatable       ::      x_values(:)
    real, allocatable       ::  sigma_values(:)
#endif

    ! --- begin ---------------------------------------------

    !! info ..
    !write (gol,*) trim(rname)//': begin'; call goPr
    !write (gol,*) trim(rname)//': mean ...'; call goPr

    ! ensure that ensemble mean is updated:
    call Ens_Mean( status )
    IF_NOTOK_RETURN(status=1)

    !! info ..
    !write (gol,*) trim(rname)//': sigma ...'; call goPr

    ! ~ storage

    !! info ..
    !write (gol,*) trim(rname)//':   init sum ...'; call goPr

    !! init local sum, not random seed needed:
    !call lsum%Init( 0, status )
    !IF_NOTOK_RETURN(status=1)

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
      !allocate( lsum_values(nvalue), stat=status )
      !IF_NOTOK_RETURN(status=1)
      allocate( sigma_values(nvalue), stat=status )
      IF_NOTOK_RETURN(status=1)
      !! init local sum:
      !sigma_values = 0.0
      ! extract from mean:
      call MAORI_State_Values_Get( x%mas, mad, x_values, status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! ~ init sum

    !! info ..
    !write (gol,*) trim(rname)//':   fill zero ...'; call goPr

    ! init sigma to zero:
    sigma%c       = 0.0
    sigma%cg      = 0.0
    sigma%aerh2o  = 0.0
    sigma%dc      = 0.0
#ifdef with_kf_meas_modis
    sigma%modis   = 0.0
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      sigma_values = 0.0
    end if
#endif

    ! ~ sum modes

    !! info ..
    !write (gol,*) trim(rname)//':   sums ...'; call goPr

    ! loop over modes:
    do j = 1, nmodes
      !  lsum := lsum + ( Ens(j) - x )**2
      sigma%c      = sigma%c        + ( Ens(j)%c       - x%c      )**2
      sigma%cg     = sigma%cg       + ( Ens(j)%cg      - x%cg     )**2
      sigma%aerh2o = sigma%aerh2o   + ( Ens(j)%aerh2o  - x%aerh2o )**2
      sigma%dc     = sigma%dc       + ( Ens(j)%dc      - x%dc     )**2
#ifdef with_kf_meas_modis
      sigma%modis   = sigma%modis   + ( Ens(j)%modis   - x%modis  )**2
#endif
#ifdef with_kf_meas_maori
      if ( nvalue > 0 ) then
        ! extract values:
        call MAORI_State_Values_Get( Ens(j)%mas, mad, ensj_values, status )
        IF_NOTOK_RETURN(status=1)
        !! testing ...
        !do qi = 1, nvalue
        !  if (qi/=281) cycle
        !  write (gol,*) 'sss ', qi, sigma_values(qi), ensj_values(qi), x_values(qi); call goPr
        !  write (gol,*) '  s ', sigma_values(qi) + ( ensj_values(qi)   - x_values(qi)  )**2; call goPr
        !end do
        ! add contribution:
        sigma_values   = sigma_values   + ( ensj_values   - x_values  )**2
      end if
#endif

    end do  ! local modes

!    ! ~ global sum
!    
!    ! info ..
!    write (gol,*) trim(rname)//':   global sums ...'; call goPr
!
!    ! add all local sums, send result to all pe:
!    call goc%AllReduce( lsum%c     , sigma%c     , status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%cg    , sigma%cg    , status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%aerh2o, sigma%aerh2o, status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%AllReduce( lsum%dc    , sigma%dc    , status )
!    IF_NOTOK_RETURN(status=1)
!#ifdef with_kf_meas_modis
!    call goc%AllReduce( lsum%modis, sigma%modis, status )
!    IF_NOTOK_RETURN(status=1)
!#endif
!#ifdef with_kf_meas_maori
!    if ( nvalue > 0 ) then
!      call goc%AllReduce( lsum_values, sigma_values, status )
!      IF_NOTOK_RETURN(status=1)
!    end if
!#endif

    ! ~ ensemble mean, post:
    
    ! aver:
    sigma%c      = sqrt( sigma%c      / (nmodes-1) )
    sigma%cg     = sqrt( sigma%cg     / (nmodes-1) )
    sigma%aerh2o = sqrt( sigma%aerh2o / (nmodes-1) )
    sigma%dc     = sqrt( sigma%dc     / (nmodes-1) )
#ifdef with_kf_meas_modis
    sigma%modis   = sqrt( sigma%modis / (nmodes-1) )
#endif
#ifdef with_kf_meas_maori
    if ( nvalue > 0 ) then
      ! average:
      sigma_values = sqrt( sigma_values / (nmodes-1) )
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
      !deallocate( lsum_values, stat=status )
      !IF_NOTOK_RETURN(status=1)
      deallocate( sigma_values, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif
    
    !! clear:
    !call lsum%Done( status )
    !IF_NOTOK_RETURN(status=1)

    !! info ..
    !write (gol,*) trim(rname)//': end'; call goPr

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
      do j = 1, nmodes
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


!  ! ********************************************************************
!  ! ***
!  ! *** budget data
!  ! ***
!  ! ********************************************************************
!  
!  
!  !
!  ! Scatter selected budgets from root (where xb is propagated)
!  ! to all other pe:
!  !
!  !   bud0%drydepos%cnh3_ave_prev(:,:)
!  !                 cso2_ave_prev(:,:)
!  ! 
!  ! The drydpos part is used in "LEKF_Driver_Init_State" in the
!  ! the call that updates deposition velocities:
!  !
!  !   ! re-compute deposition velocities if these depend on concentrations;
!  !   ! also the adhoc factor o3fac is applied if necessary ;
!  !   ! use overal budget as input, used for NH3 compensation point:
!  !   call LE_DryDepos_Setup_vd( x%c, bud0%drydepos, .true., t, status )
!  !   IF_NOTOK_RETURN(status=1)
!  !
!
!  subroutine SendBudgetData( status )
!  
!    use GO, only : goc
!
!    ! --- in/out -------------------------------------------
!
!    integer, intent(out)    ::  status
!
!    ! --- const --------------------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/SendBudgetData'
!
!    ! --- local --------------------------------------------
!
!    ! --- begin ---------------------------------------------
!    
!    ! broadcast from master pe to other pe's:
!    call goc%BCast( 0, bud0%drydepos%cnh3_ave_prev, status )
!    IF_NOTOK_RETURN(status=1)
!    call goc%BCast( 0, bud0%drydepos%cso2_ave_prev, status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ok
!    status = 0
!
!  end subroutine SendBudgetData
  
  

end module LEKF_State
