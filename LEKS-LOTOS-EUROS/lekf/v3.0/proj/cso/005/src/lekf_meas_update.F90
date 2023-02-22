!###############################################################################
!
! Measurement updates
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################

module LEKF_Meas_Update

  use GO             , only : gol, goPr, goErr
  use NetCDF         , only : NF90_StrError, NF90_NOERR
  use Grid           , only : TllGridInfo
  use LEKF_Meas_AnDom, only : T_AnDom
  use LEKF_Meas_AnDom, only : rho_frac_zero
  use LEKF_Meas_Tools, only : T_XCorr

  implicit none
  
  
  ! --- in/out -----------------------------

  private

  public    ::  T_Meas_Update
  public    ::  LEKF_Meas_Update_Timers
  

  ! --- const --------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_Update'


  ! --- types --------------------------------
  
  type T_Meas_Update
    ! update form: "EnKF", "LETKF"
    character(len=32)                      ::  method
    ! random errors?
    logical                                ::  with_v
    ! single rho/dist/corr per analysis domain?
    logical                                ::  with_rdc
    ! number of analysis domains:
    integer                                ::  ndom
    integer                                ::  n_w_test_error,n_w_test_pass
    
    !! analysis domains:
    !type(T_AnDom), pointer                 ::  dom(:)  ! (ndom)
    ! analysis domain (subdomain and area around):
    type(T_AnDom)                          ::  andom
    ! maori info (chemical and noise correlations):
    type(T_XCorr), allocatable             ::  xcorr(:)  ! (nset)
    ! work arrays for LETKF update:    
    real, allocatable                      ::  xi(:)      ! (nmodes)
    real, allocatable                      ::  xi_a(:)    ! (nmodes)
    real, allocatable                      ::  xx(:)      ! (nmodes)
  contains
    procedure ::  Init                       =>  LEKF_Meas_Update_Init
    procedure ::  Done                       =>  LEKF_Meas_Update_Done
    procedure ::  Apply                      =>  LEKF_Meas_Update_Apply
    procedure ::  CollectObservations        =>  LEKF_Meas_Update_CollectObservations
    procedure ::  EnKF_Apply                 =>  LEKF_Meas_Update_EnKF_Apply
    procedure ::  EnKF_Apply_Domain          =>  LEKF_Meas_Update_EnKF_Apply_Domain
    procedure ::  EnKF_Apply_Domain_Substate =>  LEKF_Meas_Update_EnKF_Apply_Domain_Substate
    procedure ::  LETKF_Apply                =>  LEKF_Meas_Update_LETKF_Apply
    procedure ::  LETKF_Update_Ensemble      =>  LEKF_Meas_Update_LETKF_Update_Ensemble
    procedure ::  LETKF_Update_Cell          =>  LEKF_Meas_Update_LETKF_Update_Cell
  end type T_Meas_Update


  ! --- var -------------------------------
  
  ! timers:
  integer                                ::  itim_obs_collect
  integer                                ::  itim_obs_collect_grnd
  integer                                ::  itim_grnd_vn, itim_grnd_add
  integer                                ::  itim_obs_collect_sat
  !integer                                ::  itim_obs_bcast
  integer                                ::  itim_obs_exchange
  integer                                ::  itim_domains
  !integer                                ::  itim_exchange_modes
  integer                                ::  itim_solve_gain
  integer                                ::  itim_update_modes


contains


  ! ********************************************************************
  ! ***
  ! *** module
  ! ***
  ! ********************************************************************


  subroutine LEKF_Meas_Update_Timers( status )

    use GO, only : GO_Timer_Def

    ! --- in/out -------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Timers'

    ! --- local --------------------------
    
    ! --- begin ---------------------------

    ! define timers:
    call GO_Timer_Def( itim_obs_collect     , 'lekf update obs collect'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_obs_collect_grnd, 'lekf update obs collect grnd', status )
    IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_grnd_vn, 'lekf update obs grnd vn', status )
      IF_NOTOK_RETURN(status=1)
      call GO_Timer_Def( itim_grnd_add, 'lekf update obs grnd add', status )
      IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_obs_collect_sat , 'lekf update obs collect sat' , status )
    IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_obs_bcast     , 'lekf update obs bcast'     , status )
    !IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_obs_exchange  , 'lekf update obs bcast'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_domains       , 'lekf update domains'       , status )
    IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_exchange_modes, 'lekf update exchange modes', status )
    !IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_solve_gain    , 'lekf update solve gain'    , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_update_modes  , 'lekf update update modes'  , status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Timers


  ! ********************************************************************
  ! ***
  ! *** Meas_Update
  ! ***
  ! ********************************************************************


  subroutine LEKF_Meas_Update_Init( self, rcF, status )

    use GO             , only : TrcFile
    use Grid           , only : Init
    use MAORI          , only : MAORI_Data_Inquire, MAORI_Data_Set_Inquire
    use LE_Grid        , only : ugg
    use Indices        , only : nspec_all
#ifdef with_kf_meas_maori
    use LEKF_Data      , only : mad
#endif
    use LEKF_Meas_Tools, only : GetSpecApply
    use LEKF_Noise     , only : nnoise, GetNoiseApply
#ifdef with_kf_meas_sat
    use LEKF_Meas_Sat, only : meas_sat
    use LEKF_Meas_Sat, only : T_LEKF_Meas_Sat_Set
#endif

    ! --- in/out -------------------------

    class(T_Meas_Update), intent(out)   ::  self
    type(TrcFile), intent(in)           ::  rcF
    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Init'

    ! --- local --------------------------
    
!    integer       ::  ndom_x, ndom_y
!    integer       ::  idom_x, idom_y, idom
!    integer       ::  nx_loc, ny_loc
!    integer       ::  i0, j0
    
    integer                               ::  nset_maori
    integer                               ::  nset_sat
    integer                               ::  nset, iset
    character(len=64)                     ::  sname
    logical                               ::  analyse
    character(len=1024)                   ::  key
#ifdef with_kf_meas_sat
    integer                               ::  k
    type(T_LEKF_Meas_Sat_Set), pointer    ::  msat
#endif

    ! --- begin ---------------------------

    ! info ...
    write (gol,'(a," - init analysis method ...")') rname; call goPr
    
    ! form key:
    call rcF%Get( 'kf.meas.method', self%method, status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'(a," -   method: ",a)') rname, trim(self%method); call goPr

    ! use random errors?
    select case ( trim(self%method) )
      case ( 'EnKF' )
        self%with_v   = .true.
        self%with_rdc = .true.
      case ( 'LETKF' )
        self%with_v   = .false.
        self%with_rdc = .false.
      case default
        write (gol,'("unsupported method `",a,"`")') trim(self%method); call goErr
        TRACEBACK; status=1; return
    end select

    ! info ...
    write (gol,'(a," - init analysis domain ...")') rname; call goPr

    ! init storage for analysis observations for this domain
    call self%andom%Init( self%with_v, self%with_rdc, status )
    IF_NOTOK_RETURN(status=1)

    ! ***
    ! keep track of failed and succesful inversions?
    self%n_w_test_error=0
    self%n_w_test_pass =0
    ! info ...
    write (gol,'(a," - extract chemical/noise correlation info ...")') rname; call goPr
    
#ifdef with_kf_meas_maori
    ! get number of sets to put out:
    call MAORI_Data_Inquire( mad, status, nset=nset_maori )
    IF_NOTOK_RETURN(status=1)
#else
    nset_maori = 0
#endif

    ! SAT observations?
#ifdef with_kf_meas_sat
    nset_sat = meas_sat%nset
#else
    nset_sat = 0
#endif

    ! info ...
    write (gol,'(a," -   number of maori sets : ",i2)') rname, nset_maori; call goPr
    write (gol,'(a," -   number of sat   sets : ",i2)') rname, nset_sat  ; call goPr
    ! total:
    nset = nset_maori + nset_sat

    ! any?
    if ( nset > 0 ) then
      ! storage:
      allocate( self%xcorr(nset), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! init counter:
      iset = 0
      
#ifdef with_kf_meas_maori
      ! ~ maori

      ! loop over maori sets:
      do k = 1, nset_maori
        ! increase counter:
        iset = iset + 1

        ! set name:
        call MAORI_Data_Set_Inquire( mad, iset, status, &
                                       name=sname, assim_analyse=analyse )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        write (gol,'(a," -   meas init for maori set ",a,"; analyse ",l1)') &
                         rname, trim(sname), analyse; call goPr

        ! skip ?
        if ( .not. analyse ) cycle

        ! read noise correlation keys:
        call rcf%Get( 'maori.'//trim(sname)//'.assim.spec', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_spec(nspec_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetSpecApply( key, self%xcorr(iset)%analyse_spec, self%xcorr(iset)%analyse_sia, status )
        IF_NOTOK_RETURN(status=1)

        ! read noise correlation keys:
        call rcF%Get( 'maori.'//trim(sname)//'.assim.noise', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_noise(nnoise), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetNoiseApply( key, self%xcorr(iset)%analyse_noise, status )
        IF_NOTOK_RETURN(status=1)

      end do  ! maori sets
#endif

      ! ~ Sat
      
#ifdef with_kf_meas_sat
      ! loop over sat sets:
      do k = 1, meas_sat%nset
        ! increase counter:
        iset = iset + 1
        ! pointer:
        msat => meas_sat%set(k)
        ! analyse?
        if ( msat%analyse ) then

          ! read noise correlation keys:
          call rcF%Get( 'kf.meas.sat.'//trim(msat%name)//'.spec', key, status )
          IF_NOTOK_RETURN(status=1)
          ! storage for noise correlation weights:
          allocate( self%xcorr(iset)%analyse_spec(nspec_all), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! fill:
          call GetSpecApply( key, self%xcorr(iset)%analyse_spec, self%xcorr(iset)%analyse_sia, status )
          IF_NOTOK_RETURN(status=1)

          ! read noise correlation keys:
          call rcF%Get( 'kf.meas.sat.'//trim(msat%name)//'.noise', key, status )
          IF_NOTOK_RETURN(status=1)
          ! storage for noise correlation weights:
          allocate( self%xcorr(iset)%analyse_noise(nnoise), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! fill:
          call GetNoiseApply( key, self%xcorr(iset)%analyse_noise, status )
          IF_NOTOK_RETURN(status=1)

        end if ! analyse?
      end do ! sat sets
#endif

    end if ! any sets
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Init


  ! ***


  subroutine LEKF_Meas_Update_Done( self, status )
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(out)   ::  self
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Done'

    ! --- local ----------------------------------

    integer       ::  idom

    ! --- begin ----------------------------------
    write (gol,*) 'reached end of simulation',self%n_w_test_error,'faulty inversions'; call goPr
    write (gol,*) 'reached end of simulation',self%n_w_test_pass,'passed inversions'; call goPr
    ! done:
    call self%andom%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! any maori sets?
    if ( allocated(self%xcorr) ) then
      ! storage:
      deallocate( self%xcorr, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Done
  
  
  ! ***


  subroutine LEKF_Meas_Update_Apply( self, t1, t2, status )
  
    use GO, only : TDate
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(self%method) )
    
      ! original EnKF update:
      case ( 'EnKF' )
        ! analysis per subdomain
        call self%EnKF_Apply( t1, t2, status )
        IF_NOTOK_RETURN(status=1)
    
      ! Local Ensemble Transform KF
      case ( 'LETKF' )
        ! analysis per grid cell
        call self%LETKF_Apply( t1, t2, status )
        IF_NOTOK_RETURN(status=1)
      
      ! unknown ...
      case default
        write (gol,'("unsupported update form `",a,"`")') trim(self%method); call goErr
        TRACEBACK; status=1; return
        
    end select
    
    ! ok
    status = 0
    
  end subroutine LEKF_Meas_Update_Apply
  
  
  ! ***


  !
  ! Measurement update.
  !
  ! Collect observations y, repr. err. stdvs r, and simulations for local domain,
  ! and retrieve values from nearby observations in other domains:
  !   yd, Rd, Hd Xd(1:m)
  !
  ! For each analysis domain:
  !
  !   Compute mean and stdv:
  !     Sd, Hxd, HSd
  !
  !   Equation to solve (omitted some "d" for clearity):
  !     Kd =   PH' [ HPH' +R]^{-1}
  !     Kd = SS'H' [HSS'H'+R]^{-1}
  !     Kd = S(HS)' [(HS)(HS)'+R]^{-1}
  !     Kd [(HS)(HS)'+R] = S(HS)'
  !     [(HS)(HS)'+R] Kd' = (HS)S'   # solve for rows of K on local domain
  !
  !   Update local domain:
  !     Xd(1:m) = Xd(1:m) + Kd ( yd - Hxd + vd )
  !

  subroutine LEKF_Meas_Update_EnKF_Apply( self, t1, t2, status )
  
    use GO, only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO, only : TDate
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_EnKF_Apply'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a," - apply measurement update ...")') rname; call goPr
    
    ! info ...
    write (gol,'(a," -   collect observations, assgin to domains ...")') rname; call goPr

    ! start timing:
    call GO_Timer_Start( itim_obs_collect, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Retrieve all observations y, repr. err. stdvs r, 
    ! and simulations for local modes ;
    ! results are stored into analysis domain structures:
    call self%CollectObservations( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
    
    ! end timing:
    call GO_Timer_End( itim_obs_collect, status )
    IF_NOTOK_RETURN(status=1)

    ! start timing:
    call GO_Timer_Start( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! any obs?
    if ( self%andom%nobs > 0 ) then
      ! info ...
      write (gol,'(a," -     analyse domain ...")')  rname; call goPr
      ! apply:
      call self%EnKF_Apply_Domain( t2, status )
      IF_NOTOK_RETURN(status=1)
    else
      ! info ...
      write (gol,'(a," -     no observations in domain ...")')  rname; call goPr
    end if

    ! end timing:
    call GO_Timer_End( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_EnKF_Apply
  
  
  ! ***


  subroutine LEKF_Meas_Update_CollectObservations( self, t1, t2, status )

    use GO        , only : TDate, NewDate, operator(<), operator(<=), wrtgol
    use GO        , only : GO_Timer_Start, GO_Timer_End
    use GO        , only : goc
    !use Grid      , only : WithinDistance
    !use Grid      , only : In_Domain
    use LE_Grid   , only : ugg
    use MAORI     , only : MAORI_Data_Inquire
    use MAORI     , only : MAORI_Data_Set_Inquire
    use MAORI     , only : MAORI_Data_Loc_Inquire
    use MAORI     , only : MAORI_Data_Obs_Get
    use MAORI     , only : MAORI_Data_Obs_Put
    use MAORI     , only : MAORI_State_Obs_Get
    use MAORI     , only : MAORI_SAMPLE, MAORI_TYPE_NAME
    use MAORI     , only : MAORI_ASTAT_NODATA, MAORI_ASTAT_ANALYSED
    use LEKF_State, only : nmodes
    use LEKF_State, only : Ens
#ifdef with_kf_meas_maori
    use LEKF_Data , only : mad
#endif
#ifdef with_kf_meas_sat
    use LE_Output        , only : T_LE_Output_CSO_Data
    use LEKF_Data        , only : leo
    use LEKF_State       , only : substate_sat
    use LEKF_Meas_Sat    , only : meas_sat
    use LEKF_Meas_Sat    , only : T_LEKF_Meas_Sat_Set
#endif
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_CollectObservations'

    ! --- local ----------------------------------

    integer                     ::  nset_maori
    integer                     ::  iset
    character(len=64)           ::  sname
    integer                     ::  stype
    logical                     ::  analyse
    real                        ::  rho_km, rho_m
    integer                     ::  nloc, iloc
    integer                     ::  nloc_all
    integer                     ::  loc_id
    !integer                     ::  loc_id_range(2), loc_id_range_all(2)
    real, pointer                        ::  lon, lat
    integer                     ::  obs_nvar, obs_ivar
    real                        ::  y, r
    real                        ::  alfa
    integer                     ::  idom
    logical                     ::  nearby
    logical                     ::  inside
    integer                     ::  j
    ! real, allocatable           ::   vn_all(:,:)  ! (loc_id_range,nmodes)
    ! real, allocatable           ::   v_loc(:)  ! (nmodes)
    real, allocatable           ::  HX_loc(:)  ! (nmodes)

    integer                     ::  nmeas_tot,nmeas_tot_loc
    integer                     ::  nmeas_intimestep_tot
    integer                     ::  ipix
    type(TDate)                 ::  t
    integer                     ::  inds(0:7)
    
#ifdef with_kf_meas_sat
!    integer                               ::  displ
    ! real, allocatable                     ::  chi(:,:)    ! (nmodes,npix)
    ! real, allocatable                     ::  chi_tot(:,:)    ! (nmodes,nmeas_tot) ! only on root
    
    integer                               ::  k, i
    type(T_LEKF_Meas_Sat_Set), pointer    ::  msat
    type(T_LE_Output_CSO_Data), pointer   ::  satd
    real                                  ::  bbox(4)
    integer                               ::  glbid
	real,pointer                          ::  yloc(:), rloc(:,:)
	integer,pointer                          ::  nr
#endif

    !! testing ...
    !real, allocatable   ::  vv(:,:)  ! (nmodes_loc,npix)
    !integer             ::  ipix1, ipix2, ipixstep
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a," - start")') rname; call goPr

    ! reset to empty:
    call self%andom%Reset( status )
    IF_NOTOK_RETURN(status=1)
    
    !! storage for ensemble simulations,
    !! dummy in case no modes are allocated here (master pe with xb only):
    !allocate(  v_loc(max(1,nmodes)), stat=status )
    !IF_NOTOK_RETURN(status=1)
    !allocate( HX_loc(max(1,nmodes)), stat=status )
    !IF_NOTOK_RETURN(status=1)
    
    ! storage for ensemble simulations:
    ! if ( self%with_v ) then
      ! allocate(  v_loc(nmodes), stat=status )
      ! IF_NOTOK_RETURN(status=1)
    ! end if
    allocate( HX_loc(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! init counter:
    iset = 0

    ! *

    ! start timing:
    call GO_Timer_Start( itim_obs_collect_grnd, status )
    IF_NOTOK_RETURN(status=1)
    
#ifdef with_kf_meas_maori
    ! get number of maori sets:
    call MAORI_Data_Inquire( mad, status, nset=nset_maori )
    IF_NOTOK_RETURN(status=1)
#else
    nset_maori = 0
#endif
    
    ! info ..
    if ( nset_maori > 0 ) then
      write (gol,'(a," - loop over ",i0," maori sets ...")') rname, nset_maori; call goPr
    else
      write (gol,'(a," - no maori sets ...")') rname; call goPr
    end if

#ifdef with_kf_meas_maori
    ! loop over maori sets:
    do k = 1, nset_maori
      ! increase counter:
      iset = iset + 1

      ! get type of output, set name, and assimilation parameters:
      call MAORI_Data_Set_Inquire( mad, iset, status, &
                                       name=sname, type=stype, &
                                       assim_analyse=analyse, assim_rho=rho_km )
      IF_NOTOK_RETURN(status=1)
      ! convert ...
      rho_m = rho_km * 1.0e3
    
      ! info ..
      write (gol,'(a," -   set `",a,"` ...")') rname, trim(sname); call goPr
      
      ! skip if not to be analyzed ...
      if ( .not. analyse ) then
        ! info ...
        write (gol,'(a," -     not to be analyzed ...")') rname; call goPr
        ! next:
        cycle
      end if

      ! do something given the type:
      select case ( stype )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( MAORI_SAMPLE )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! get number of locations (on this domain),
          ! and min/max of location id's:
          call MAORI_Data_Set_Inquire( mad, iset, status, &
                                        nloc=nloc, obs_nvar=obs_nvar )
                                        !loc_id_range=loc_id_range, 
          IF_NOTOK_RETURN(status=1)
    
          ! info ..
          write (gol,'(a," -     number of locations: ",i0)') rname, nloc; call goPr
          
          ! global total on all domains:
          nloc_all = nloc
          call goc%AllReduce( 'sum', nloc_all, status )
          IF_NOTOK_RETURN(status=1)
          ! info ..
          write (gol,'(a," -     total number       : ",i0)') rname, nloc_all; call goPr
          
!          ! global location id range:
!          loc_id_range_all = loc_id_range
!          call goc%AllReduce( 'min', loc_id_range_all(1), status )
!          IF_NOTOK_RETURN(status=1)
!          call goc%AllReduce( 'max', loc_id_range_all(2), status )
!          IF_NOTOK_RETURN(status=1)
!          ! info ..
!          write (gol,'(a," -     location id range : ",i0,":",i0," of ",i0,":",i0)') &
!                            rname, loc_id_range, loc_id_range_all; call goPr
!          ! any values?
!          if ( nloc > 0 ) then
!            ! check ...
!            if ( loc_id_range_all(1) > loc_id_range_all(2) ) then
!              write (gol,'("strange loc_id range ",i0,":",i0)') loc_id_range_all; call goPr
!              TRACEBACK; status=1; return
!            end if
!            ! start timing:
!            call GO_Timer_Start( itim_grnd_vn, status )
!            IF_NOTOK_RETURN(status=1)
!            ! storage for random represenation error samples:
!            allocate(  vn_all(loc_id_range_all(1):loc_id_range_all(2),nmodes), stat=status )
!            IF_NOTOK_RETURN(status=1)
!            ! loop over modes:
!            do j = 1, nmodes
!              ! random error out N(0,1)
!              call Ens(j)%rnd%Get_Normal( vn_all(:,j), status )
!              IF_NOTOK_RETURN(status=1)
!            end do ! modes
!            ! end timing:
!            call GO_Timer_End( itim_grnd_vn, status )
!            IF_NOTOK_RETURN(status=1)
!          end if ! loc id range defined

          ! loop over locations:
          do iloc = 1, nloc

            ! get location:
            call MAORI_Data_Loc_Inquire( mad, iset, iloc, status, &
                                            loc_id=loc_id, lon=lon, lat=lat )
            IF_NOTOK_RETURN(status=1)
            
            !! check ...
            !if ( (loc_id < lbound(vn_all,1)) .or. (loc_id > ubound(vn_all,1)) ) then
            !  write (gol,'("location id ",i0," out of expected range ",i0,":",i0)') &
            !                  loc_id, lbound(vn_all,1), ubound(vn_all,2); call goErr
            !  TRACEBACK; status=1; return
            !end if

            ! loop over variables:
            do obs_ivar = 1, obs_nvar

              ! extract measurement:
              call MAORI_Data_Obs_Get( mad, iset, obs_ivar, status, &
                                        iloc=iloc, y=y, r=r, alfa=alfa  )
              IF_NOTOK_RETURN(status=1)
              
              ! nan ? then skip:
              if ( y <= 0.0 ) then
                ! set status flag:
                call MAORI_Data_Obs_Put( mad, iset, obs_ivar, status, &
                                   iloc=iloc, astat_ibset=MAORI_ASTAT_NODATA )
                IF_NOTOK_RETURN(status=1)
                ! skip:
                cycle
              end if
              
              ! no screening implemented anywore ..
              ! where has that gone?
              
              ! accepted for assimilation:
              call MAORI_Data_Obs_Put( mad, iset, obs_ivar, status, &
                                        iloc=iloc, astat_ibset=MAORI_ASTAT_ANALYSED )
              IF_NOTOK_RETURN(status=1)
              
              ! loop over local modes:
              do j = 1, nmodes
                ! extract simulation of measurements:
                call MAORI_State_Obs_Get( Ens(j)%mas, mad, iset, obs_ivar, status, &
                                            iloc=iloc, value=HX_loc(j) )
                IF_NOTOK_RETURN(status=1)
                ! random errors?
                ! if ( self%with_v ) then
                  ! ! random error out N(0,r)
                  ! !call Ens(j)%rnd%Get_Normal( v_loc(j), status, sigma=r )
                  ! !IF_NOTOK_RETURN(status=1)
                  ! ! use random number out of N(0,1) that is independent of (sub)domain:
                  ! v_loc(j) = vn_all(loc_id,j) * r
                ! end if
              end do
              
              !! testing ...
              ! write (gol,*) 'vvv loc_id=', loc_id, ' v=', vn_all(loc_id,:); call goPr

              !! loop over analysis domains:
              !do idom = 1, self%ndom
              
              !! within factor*rho [m] of current domain?
              !call ugg%WithinDistance( lon, lat, rho_frac_zero * rho_m, nearby, status )
              !IF_NOTOK_RETURN(status=1)
              !! nearby observation? then add to list:
              !if ( nearby ) then

                ! check if in domain (not in surrounding area):
                call ugg%InDomain( lon, lat, inside, status )
                IF_NOTOK_RETURN(status=1)
                ! expected to be inside ..
                if ( .not. inside ) then
                  write (gol,'("fould location outside domain")'); call goErr
                  TRACEBACK; status=1; return
                end if

                ! start timing:
                call GO_Timer_Start( itim_grnd_add, status )
                IF_NOTOK_RETURN(status=1)
                ! add to domain:
                call self%andom%Add_Observation( (/iset,iloc,obs_ivar/), inside, &
                                                         lon, lat, rho_m, &
                                                         y, r, HX_loc, &
                                                         status)!, v_loc=v_loc )
                IF_NOTOK_RETURN(status=1)

                !! info ...
                !write (gol,'(a," -   added obs (set ",i0,", loc ",i0,", var ",i0,") at (",f8.2,",",f8.2,") to domain ",i0)') &
                !                        rname, iset, iloc, obs_ivar, lon, lat, idom; call goPr

                ! end timing:
                call GO_Timer_End( itim_grnd_add, status )
                IF_NOTOK_RETURN(status=1)

              !end if ! nearby current domain
              
              !end do  ! analysis domains

            end do  ! observed variables

          end do  ! locations

          ! storage?
          ! if ( allocated(vn_all) ) then
           ! ! clear:
           ! deallocate(  vn_all, stat=status )
           ! IF_NOTOK_RETURN(status=1)
          ! end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("ERROR - unsupported maori type ",i6," (`",a,"`) for set ",i6)') &
                                   stype, trim(MAORI_TYPE_NAME(stype)), iset; call goPr
          TRACEBACK; status=1; return

      end select

    end do ! maori sets
#endif

    ! end start timing:
    call GO_Timer_End( itim_obs_collect_grnd, status )
    IF_NOTOK_RETURN(status=1)

    ! *
    
    ! begin start timing:
    call GO_Timer_Start( itim_obs_collect_sat, status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_kf_meas_sat
    ! info ..
    write (gol,'(a," - collect sat observations ...")') rname; call goPr
    ! loop over sat sets:
    do k = 1, meas_sat%nset
      ! set index:
      iset = iset + 1
      ! pointers:
      msat => meas_sat%set(k)
      satd => leo%cso_data(meas_sat%iout(k))
      ! info ..
      write (gol,'(a," -   sat ",i0," `",a,"`, iset=",i0," ...")') rname, k, trim(msat%name), iset; call goPr
      ! analyse?
      if ( len_trim(satd%orbit_filename) > 0 ) then

        ! target file:
       ! write (fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2,"_data.nc")') &
       ! trim(outputdir), trim(self%id_model), trim(self%id_expid), &
       ! trim(self%name), t%year, t%month, t%day, t%hour, t%min
       write (gol,*) 'new data for analysis, from', trim(satd%orbit_filename); call goPr
       ! info ..
       ! write (gol,'(a,":   create ",a," ...")') rname, trim(fname); call goPr

      ! write:

      if ( msat%analyse ) then

        ! info ...
        write (gol,'(a," -   to be analyzed ...")') rname; call goPr

        ! total number over all processors:
        nmeas_tot = satd%sdata%nglb
        call goc%ParInfo( satd%sdata%npix, status, ntot=nmeas_tot_loc )!, displ=displ )
        ! call goc%ParInfo( satd%npix, status, ntot=nmeas_tot )!, displ=displ )
        IF_NOTOK_RETURN(status=1)

        ! any measurements?
        if ( nmeas_tot_loc > 0 ) then

          ! info ...
          write (gol,'(a," -   number of pixels : ",i0," (local sum ",i0,"),(trop file ",i0,")")') rname, satd%sdata%npix,nmeas_tot_loc, nmeas_tot; call goPr
          
          ! random errors?
          ! if ( self%with_v ) then
            ! info ...
            ! generate random N(0,1) numbers per pixel and mode,
            ! get them now to have later v = chi * sigma
            ! storage:
            ! allocate( chi(nmodes,max(1,satd%npix)), stat=status )
            ! allocate( chi(nmodes,nmeas_tot), stat=status )
            ! IF_NOTOK_RETURN(status=1)
	  	  ! set to zero
	  	  ! chi(:,:) = 0.0

            ! generate random numbers on root only
            ! if (goc%root) then
              ! storage:
              ! allocate( chi_tot(nmodes,nmeas_tot), stat=status ) 
              ! IF_NOTOK_RETURN(status=1)
              ! write (gol,*) 'has shape',  shape(chi_tot(j,:)), 'chi is',nmodes,max(1,satd%npix); call goPr
              ! write (gol,*) 'has shape',  shape(chi(:,:)), 'chi is',nmodes,nmeas_tot; call goPr
              ! loop over modes:
              ! ! do j = 1, nmodes
                ! ! ! get random nunmbers:
                ! ! call satd%sdata%GetRandomAll( msat%rnd, chi(j,:), status )
                  ! ! ! write (gol,*) '(mode ( ', j, 'first 5',chi(j,1:5); call goPr
                  ! ! ! write (gol,*) '(mode ( ', j, 'mid somewhere 5',chi(j,1406:1406+5); call goPr
                  ! ! ! write (gol,*) '(mode ( ', j, 'last 5',chi(j,nmeas_tot-4:nmeas_tot); call goPr
                ! ! IF_NOTOK_RETURN(status=1)
              ! ! end do ! modes
            ! else
              ! ! dummy
                  ! allocate( chi_tot(1,1), stat=status ) 
                  ! IF_NOTOK_RETURN(status=1)
	  		! ! write (gol,*) 'test',  chi(1,1:10); call goPr
	  		! ! write (gol,*) 'test',  chi(1,nmeas_tot-10:nmeas_tot); call goPr
             ! end if ! root
          ! distribute:
          ! TBD: go_comm.F90, GO_Comm_ScatterV_r4_2d  net als GatherV_r4_2d
          !   denk dat we ook dims om moeten draaien: chi(nmodes,npix)
          !    want npix kan soms 0 zijn (nloc argument)
          ! call goc%ScatterV( chi_tot, chi, status, nloc=satd%npix )
          ! call goc%BCast( goc%root_id, chi, status)
                ! IF_NOTOK_RETURN(status=1)
			! write (gol,*) 'test return',  chi(1,1:10); call goPr
			! write (gol,*) 'test return',  chi(1,nmeas_tot-10:nmeas_tot); call goPr
          ! clear:
          ! deallocate( chi_tot, stat=status )
          ! IF_NOTOK_RETURN(status=1)
          
          ! ! end if ! with_v
          ! convert from km to m:
          rho_m = msat%rho * 1.0e3  ! m

          ! screening factor:
          alfa = msat%screening_factor

          ! loop over measurements assigned to this domain:
		  ! if ( self%with_v ) then
			! write (gol,*) 'vvv1 HXloc = ',' ipix, glbid, lon,lat,y,HX_loc(1),v_loc(1),r,chi(1,glbid)'; call goPr
		  ! else
		    ! write (gol,*) 'vvv1 HXloc = ',' ipix, glbid, lon,lat,y,HX_loc(1)'; call goPr
		  ! end if 
          do ipix = 1, satd%sdata%npix

            ! check if not accepted (validation, no data?):
            ! if ( satd%astatus%data(ipix) /= 0 ) then
               ! write (gol,'("sat observation status invalid: ",i0)') satd%astatus%data(ipix); call goErr
               ! TRACEBACK; status=1; return         
            ! end if
			! write (gol,*) 'reached ipix = ', ipix; call goPr
            ! extract measurement data:
            call satd%sdata%GetPixel( ipix, status, glbid=glbid, lon=lon, lat=lat, yr=yloc, R=rloc)
            IF_NOTOK_RETURN(status=1)
            ! write (gol,*) 'reached glbid = ', glbid; call goPr
            ! write (gol,*) 'reached lon = ', lon; call goPr
            ! write (gol,*) 'reached lat = ', lat; call goPr
            ! write (gol,*) 'reached yloc = ', yloc; call goPr
            ! write (gol,*) 'reached rloc = ', rloc; call goPr
            ! write (gol,*) 'reached satd%sdata%nr = ', satd%sdata%nr; call goPr
			! either set to first point or sum... in principle we only assimilate single values
			! write (gol,*) 'nr',satd%sdata%nr ; call goPr
			y = sum(yloc(1:satd%sdata%nr))
			! same for the covariance/std, keep option for future profiles
			r = 0.0d0
			do i = 1, satd%sdata%nr
				r = r + rloc(i,i)
			end do
            ! define how to extract value from state:
            inds(0:2) = (/ substate_sat, k, ipix /)

            ! loop over modes:
            do j = 1, nmodes
              ! extract simulation of measurements:
			  ! all profile pixel matching happens in le_output_sat.F90
              call Ens(j)%GetValue( inds, HX_loc(j), status )
              IF_NOTOK_RETURN(status=1)
              ! random erors?
              ! if ( self%with_v ) then
                ! ! random error out N(0,r):
                ! ! v_loc(j) = chi(j,ipix) * r
                ! v_loc(j) = chi(j,glbid) * r
              ! end if !with v
            end do
            
            !! testing ..
            ! write (gol,*) 'reached end of ensemble loop'; call goPr
            ! write (gol,*) 'vvv2 chi = ', ipix, glbid, chi(:,glbid); call goPr
			! if ( self%with_v ) then
				! write (gol,*) 'vvv1 HXloc = ', ipix, glbid, lon,lat,y,HX_loc(1),v_loc(1),r,chi(1,glbid); call goPr
            ! else
				! write (gol,*) 'vvv1 HXloc = ', ipix, glbid, lon,lat,y,HX_loc(1); call goPr
			! end if
            ! write (gol,*) 'vvv1 v = ', ipix, glbid, chi(:,ipix); call goPr

            ! check if in domain (not in surrounding area):
            call ugg%InDomain( lon, lat, inside, status )
            IF_NOTOK_RETURN(status=1)
            ! expected to be inside ..
            ! ! if ( .not. inside ) then
              ! ! call ugg%GetBoundingBox( bbox, status )
              ! ! IF_NOTOK_RETURN(status=1)
              ! ! write (gol,'("fould location outside domain:")'); call goErr
              ! ! write (gol,'("  iset    : ",i0)') iset; call goErr
              ! ! write (gol,'("  ipix    : ",i0)') ipix; call goErr
              ! ! write (gol,'("  glbid   : ",i0)') glbid; call goErr
              ! ! write (gol,'("  lon,lat : ",2f8.2)') lon,lat; call goErr
              ! ! write (gol,'("  bbox    : ",4f8.2)') bbox; call goErr
              ! ! TRACEBACK; status=1; return
            ! ! end if
			! choice to use the pixel in analysis:
			! if (satd%assim_flag%data(ipix) .ne. 2) then
				! satd%astatus%data(ipix) = 1
				! write (gol,*) 'astatus', satd%astatus%data(ipix),shape(satd%astatus%data(:)),satd%npix,nmeas_tot_loc; call goPr
				! write (gol,*) 'assim_flag', satd%assim_flag%data(ipix),shape(satd%astatus%data(:)),satd%npix,nmeas_tot_loc; call goPr

				
			! elseif (satd%assim_flag%data(ipix) .not. 2) then
				! satd%astatus%data(ipix) = 1
			! else
				! add to domain, identify by set number (used for chem/noise correlations)
				! and pixel number:
			call self%andom%Add_Observation( (/iset,ipix,glbid/), inside, &
													 lon, lat, rho_m, &
													 y, r, HX_loc, &
													 status)!, v_loc=v_loc )!, debug=.true. )
			IF_NOTOK_RETURN(status=1)
			! end if
          end do ! pixels
          ! testing
           ! call ugg%GetBoundingBox( bbox, status )
           ! IF_NOTOK_RETURN(status=1)
          ! write (gol,*)'(" end of domain "' ; call goPr
          ! write (gol,'("  iset    : ",i0)') iset; call goPr
          ! write (gol,'("  ipix    : ",i0)') ipix; call goPr
          ! write (gol,'("  glbid    : ",i0)') glbid; call goPr
          ! write (gol,'("  lon,lat : ",2f8.2)') lon,lat; call goPr
          ! write (gol,'("  bbox    : ",4f8.2)') bbox; call goPr
          ! clear:
          ! if ( self%with_v ) then
            ! deallocate( chi, stat=status )
            ! IF_NOTOK_RETURN(status=1)
          ! end if

        else

          ! info ...
          write (gol,'(a," -   no pixels for this timestep ...")') rname; call goPr

        end if  ! nmeas_tot > 0

      else

        ! info ...
        write (gol,'(a," -   not to be analyzed ...")') rname; call goPr

      end if  ! analyse?
    else
    
      ! info ..
      write (gol,'(a,":   no new data for this time, no analysis ...")') rname; call goPr
      
    end if ! filename check
    end do ! sat sets
#endif

    ! end start timing:
    call GO_Timer_End( itim_obs_collect_sat, status )
    IF_NOTOK_RETURN(status=1)

    ! *
        
    ! clear:
    deallocate( HX_loc, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! if ( self%with_v ) then
      ! deallocate( v_loc, stat=status )
      ! IF_NOTOK_RETURN(status=1)
    ! end if
    
    ! info ..
    write (gol,'(a," - exchange nearby observations between domains ...")') rname; call goPr
    
    ! start timing:
    call GO_Timer_Start( itim_obs_exchange, status )
    IF_NOTOK_RETURN(status=1)
    ! exchange:
    call self%andom%Exchange_Nearby( status )
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End( itim_obs_exchange, status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ..
    write (gol,'(a," - end")') rname; call goPr

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_CollectObservations
  
  
  ! ***


  !
  ! For each model domain:
  !
  !   Broadcast modes:
  !     Ed(1:m)
  !
  !   Compute mean and stdv:
  !     Sd, Hxd, HSd
  !
  !   Equation to solve (ommit "d" for clean notations):
  !     Kd =   PH' [ HPH' +R]^{-1}
  !     Kd = SS'H' [HSS'H'+R]^{-1}
  !     Kd = S(HS)' [(HS)(HS)'+R]^{-1}
  !     Kd [(HS)(HS)'+R] = S(HS)'
  !     [(HS)(HS)'+R] Kd' = (HS)S'   # solve for rows of K on analysis domain
  !
  !   Update local members for analysis domain:
  !     Xd(1:ml) = Xd(1:ml) + Kd ( yd - Hxd + v )
  !

  subroutine LEKF_Meas_Update_EnKF_Apply_Domain( self, t, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use Dims      , only : nspec
    use LEKF_State, only : substate_c, substate_cg, substate_aerh2o, substate_dc
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_EnKF_Apply_Domain'

    ! --- local ----------------------------------
    
    !type(T_AnDom), pointer    ::  dom
    integer                   ::  ispec

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a," - measurement update for domain")') rname; call goPr
    
    !! short:
    !dom => self%dom(idom)
    
    ! fill observation simulation entities:
    !    dom%HSd
    !    dom%HPHR
    !    dom%U
    ! allocated on return:
    call self%andom%Fill_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    !! testing ...
    !call self%andom%Dump( status, t=t )
    !IF_NOTOK_RETURN(status=1)
    
    ! call for specific elements;
    ! update concentration arrays per spec to limit memory usage:
    do ispec = 1, nspec
      call self%EnKF_Apply_Domain_SubState( substate_c     , ispec, t, status )
      IF_NOTOK_RETURN(status=1)
      call self%EnKF_Apply_Domain_SubState( substate_cg    , ispec, t, status )
      IF_NOTOK_RETURN(status=1)
    end do
    ! update full arrays, pass dummy ispec:
    call self%EnKF_Apply_Domain_SubState( substate_aerh2o, -999, t, status )
    IF_NOTOK_RETURN(status=1)
    call self%EnKF_Apply_Domain_SubState( substate_dc    , -999, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear current arrays
    call self%andom%Clear_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_EnKF_Apply_Domain
  
  
  ! *


  subroutine LEKF_Meas_Update_EnKF_Apply_Domain_SubState( self, substate, ispec, t, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use GO        , only : goc
    use Num       , only : LinAlg_Sym_FactorSolve
    use Dims      , only : nx, ny, nz, nspec
    use LEKF_noise, only : nnoise, nhist
    use LEKF_State, only : nmodes
    !use LEKF_State, only : imodes, imode_pe
    use LEKF_State, only : Ens
    use LEKF_State, only : substates, substate_c, substate_cg, substate_aerh2o, substate_dc
  
      use file_nc
      use le_grid, only : ledom => dom

    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    !integer, intent(in)                   ::  idom
    integer, intent(in)                   ::  substate
    integer, intent(in)                   ::  ispec
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_EnKF_Apply_Domain_SubState'

    ! --- local ----------------------------------
    
    !type(T_AnDom), pointer    ::  dom
!    integer                   ::  nxd, nyd
!    integer                   ::  i1, i2, j1, j2
    integer                   ::  j
!    integer                   ::  imode
    integer                   ::  n1
    integer                   ::  shp(4)
    real, allocatable         ::  HSSt(:,:)  ! (nobs,n1)
    real, allocatable         ::  Kdt(:,:)   ! (nobs,n1)
    real, allocatable         ::  Ed_c(:,:,:,:,:)   ! (nxd,nyd,nz,nspec,nmodes)
    real, allocatable         ::  Sd_c(:,:,:,:,:)   ! (nxd,nyd,nz,nspec,nmodes)
    real, allocatable         ::  xd_c(:,:,:,:)     ! (nxd,nyd,nz,nspec)
    real, allocatable         ::  Kd_c(:,:,:,:,:)   ! (nobs,nxd,nyd,nz,nspec)
    integer                   ::  iobs
    integer                   ::  iset
    real                      ::  res
    integer                   ::  iz
    integer                   ::  inoise, ihist
    
      !character(len=1024)   ::  dumpfile
      !integer               ::  iy
      !integer               ::  off(2)

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a," -   state element `",a,"` ...")') rname, trim(substates(substate)); call goPr
    write (gol,'(a," -     fill matrices ...")') rname; call goPr
    
    !! testing ...
    !call ledom%Get( status, off=off )
    !IF_NOTOK_RETURN(status=1)

    !! short:
    !dom => self%dom(idom)
    
    ! check ispec argument ...
    if ( (substate==substate_c) .or. (substate==substate_cg) ) then
      if ( (ispec < 1) .or. (ispec > nspec) ) then
        write (gol,'("ispec ",i0," out of range 1,..",i0," for substate ",i0)') &
                        ispec, nspec, substate; call goErr
        TRACEBACK; status=1; return
      end if
    else
      if ( ispec > 0 ) then
        write (gol,'("ispec ",i0," should be undefined (<0) for substate ",i0)') &
                        ispec, substate; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    !!...............................................
    !! start timing:
    !call GO_Timer_Start( itim_exchange_modes, status )
    !IF_NOTOK_RETURN(status=1)
    !!...............................................
    
    ! input shape:
    select case ( substate )
      case ( substate_c      ) ; shp = (/nx,ny,nz,1/)   ! single ispec
      case ( substate_cg     ) ; shp = (/nx,ny,1 ,1/)   ! single ispec
      case ( substate_aerh2o ) ; shp = (/nx,ny,nz,1    /)
      case ( substate_dc     ) ; shp = (/nx,ny,nnoise,nhist/)
      case default
        write (gol,'("unsupported substate ",i0)') substate; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! temporary storage:
    allocate( Ed_c(shp(1),shp(2),shp(3),shp(4),nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Sd_c(shp(1),shp(2),shp(3),shp(4),nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( xd_c(shp(1),shp(2),shp(3),shp(4)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Kd_c(self%andom%nobs,shp(1),shp(2),shp(3),shp(4)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of elements in 1d:
    n1 = size(xd_c)
    ! storage to solve gain:
    allocate( HSSt(self%andom%nobs,n1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Kdt(self%andom%nobs,n1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! init to zero, needed if 'AllReduce' is used to exchange the modes by a sum over the pe's:
    Ed_c = 0.0
    ! loop over modes:
    do j = 1, nmodes
      !! global mode index:
      !imode = imodes(j)
      ! copy local slab:
      select case ( substate )
        case ( substate_c      ) ; Ed_c(:,:,:,1,j) = Ens(j)%c     (:,:,:,ispec)
        case ( substate_cg     ) ; Ed_c(:,:,1,1,j) = Ens(j)%cg    (:,:  ,ispec)
        case ( substate_aerh2o ) ; Ed_c(:,:,:,1,j) = Ens(j)%aerh2o(:,:,:)
        case ( substate_dc     ) ; Ed_c(:,:,:,:,j) = Ens(j)%dc    (:,:,:,:)
        case default
          write (gol,'("unsupported substate ",i0)') substate; call goErr
          TRACEBACK; status=1; return
      end select
    end do ! local modes
    !
    !! sum all local arrays (with many zeros ...) using a 'reduce',
    !! and return the result to all processes (input is the same as output):
    !call goc%AllReduce( Ed_c, status )
    !IF_NOTOK_RETURN(status=1)
    
    ! ensemble mean:
    xd_c = sum( Ed_c, dim=5 )/nmodes
    ! ensemble covariance square root:
    do j = 1, nmodes
      Sd_c(:,:,:,:,j) = ( Ed_c(:,:,:,:,j) - xd_c )/sqrt(nmodes-1.0)
    end do

    !...............................................
    !! switch timing:
    !call GO_Timer_Switch( itim_exchange_modes, itim_solve_gain, status )
    !IF_NOTOK_RETURN(status=1)
    ! switch timing:
    call GO_Timer_Start( itim_solve_gain, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! solve gain from:
    !    [(HS)(HS)'+R] Kd' = (HS)S'
    
    ! fill right hand side:
    HSST = matmul( self%andom%HSd, transpose(reshape(Sd_c,(/n1,nmodes/))) )

    !! solve:
    !call LinAlg_SolveSym( dom%HPHR, HSSt, Kdt, status )
    !IF_NOTOK_RETURN(status=1)
    
    ! info ...
    write (gol,'(a," -     solve (",i0,",",i0,") system ...")') rname, shape(Kdt); call goPr

    ! solve using upper-triangular factorization:
    call LinAlg_Sym_FactorSolve( self%andom%U, HSSt, Kdt, status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ...
    write (gol,'(a," -     distribute result ...")') rname; call goPr
    
    ! reshape:
    Kd_c = reshape( Kdt, (/self%andom%nobs,shp(1),shp(2),shp(3),shp(4)/) )

    !...............................................
    ! switch timing:
    call GO_Timer_Switch( itim_solve_gain, itim_update_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
!    ! dump ...
!    write (dumpfile,'("HPHR__",i4.4,2i2.2,"_",2i2.2,"__dom_",i2.2,".nc")') &
!              t%year, t%month, t%day, t%hour, t%min, goc%id
!    call nc_dump( trim(dumpfile), self%andom%HPHR, 'HPHR', (/'nobs ','nobs2'/), status )
!    IF_NOTOK_RETURN(status=1)
!    ! dump ...
!    write (dumpfile,'("HSd__",i4.4,2i2.2,"_",2i2.2,"__dom_",i2.2,"__substate_",a,".nc")') &
!              t%year, t%month, t%day, t%hour, t%min, goc%id, trim(substates(substate))
!    call nc_dump( trim(dumpfile), self%andom%HSd, 'Sd', (/'nobs ','nmode'/), status )
!    IF_NOTOK_RETURN(status=1)
!    ! dump ...
!    write (dumpfile,'("Sd__",i4.4,2i2.2,"_",2i2.2,"__dom_",i2.2,"__substate_",a,".nc")') &
!              t%year, t%month, t%day, t%hour, t%min, goc%id, trim(substates(substate))
!    call nc_dump( trim(dumpfile), Sd_c, 'Sd', (/'x    ','y    ','z    ','s    ','nmode'/), status )
!    IF_NOTOK_RETURN(status=1)
!    ! dump ...
!    write (dumpfile,'("Kd__",i4.4,2i2.2,"_",2i2.2,"__dom_",i2.2,"__substate_",a,".nc")') &
!              t%year, t%month, t%day, t%hour, t%min, goc%id, trim(substates(substate))
!    call nc_dump( trim(dumpfile), Kd_c, 'Kd', (/'nobs','x   ','y   ','z   ','s   '/), status )
!    IF_NOTOK_RETURN(status=1)
    
    ! safety check ...
    if ( .not. self%with_rdc ) then
      write (gol,'("rho/dist/corr arrays not enabled")'); call goErr
      TRACEBACK; status=1; return
    end if
 
    ! loop over local modes:
    do j = 1, nmodes
      !! global mode index:
      !imode = imodes(j)
      ! update ensemble member:
      !    x := x + K (y-Hx+v)
      ! take care of spatial and chemical correlations:
      do iobs = 1, self%andom%nobs
        ! residue:
        !   d =      y           +          vj          -           HXj
        res = self%andom%y(iobs) + self%andom%v(iobs,j) - self%andom%HX(iobs,j)
        ! TESTING: no random errors ...
        !res = self%andom%y(iobs) + 0.0                  - self%andom%HX(iobs,j)
        !print *, 'xxx1 mode ', j,' iobs ', iobs
        !print *, '  x1 y = ', dom%y(iobs), 'hx=', dom%HX(iobs,j), 'v=', dom%v(iobs,j)
        !print *, '  x1 res = ', res
        ! current observation set:
        iset = self%andom%inds(iobs)%iset
        ! switch:
        select case ( substate )
          !~
          case ( substate_c )
            ! analyse?
            if ( self%xcorr(iset)%analyse_spec(ispec) ) then
              ! loop over levels:
              do iz = 1, nz
                ! add contribution, apply spatial localization:
                Ens(j)%c(:,:,iz,ispec) = Ens(j)%c(:,:,iz,ispec) + &
                             Kd_c(iobs,:,:,iz,1) * self%andom%an_corr(iobs,:,:) * res
              end do ! iz
            end if ! analyse spec?
          !~
          case ( substate_cg )
            ! analyse?
            if ( self%xcorr(iset)%analyse_spec(ispec) ) then
              ! add contribution, apply spatial localization:
              Ens(j)%cg(:,:,ispec) = Ens(j)%cg(:,:,ispec) + &
                           Kd_c(iobs,:,:,1,1) * self%andom%an_corr(iobs,:,:) * res
            end if ! analyse spec?
          !~
          case ( substate_aerh2o )
            ! analyse sia tracers? then also update aerh2o:
            if ( self%xcorr(iset)%analyse_sia ) then
              ! loop over levels:
              do iz = 1, nz
                ! add contribution, apply spatial localization:
                Ens(j)%aerh2o(:,:,iz) = Ens(j)%aerh2o(:,:,iz) + &
                               Kd_c(iobs,:,:,iz,1) * self%andom%an_corr(iobs,:,:) * res
              end do ! iz
            end if
          !~
          case ( substate_dc )
            ! loop over noise fields:
            do inoise = 1, nnoise
              !print *, 'xxx2 inoise ', inoise, self%xcorr(iset)%analyse_noise(inoise)
              !print *, 'xxx2   fc ', minval(Ens(j)%dc(:,:,inoise,:)), maxval(Ens(j)%dc(:,:,inoise,:))
              !print *, 'xxx2   nhist = ', nhist
              ! analyse?
              if ( self%xcorr(iset)%analyse_noise(inoise) ) then
                ! loop over history:
                do ihist = 1, nhist
                  !print *, 'xxx2   Kd   range ', minval(Kd_c(iobs,:,:,inoise,ihist)), maxval(Kd_c(iobs,:,:,inoise,ihist))
                  !print *, 'xxx2   corr range ', minval(dom%corr(iobs,:,:)), maxval(dom%corr(iobs,:,:))
                  ! add contribution, apply spatial localization:
                  Ens(j)%dc(:,:,inoise,ihist) = Ens(j)%dc(:,:,inoise,ihist) + &
                                 Kd_c(iobs,:,:,inoise,ihist) * self%andom%an_corr(iobs,:,:) * res
                  !
                  !! testing ...
                  !if ( ihist==1 ) then
                  !  do iy = 1, ny
                  !    !if ( self%andom%corr(iobs,1,iy) > 0.0 ) then
                  !    if ( off(2)+iy == 31 ) then
                  !      write (gol,*) 'aaa1 mode=',j,'; iobs=',iobs,'; giy=',off(2)+iy, &
                  !                       '; ivar=', self%andom%inds(iobs)%obs_ivar, &
                  !                       '; lon,lat=', self%andom%inds(iobs)%lon, self%andom%inds(iobs)%lat, &
                  !                       '; dc=',Ens(j)%dc(1,iy,inoise,ihist), &
                  !                       '; Kd=', Kd_c(iobs,1,iy,inoise,ihist), &
                  !                       '; corr=', self%andom%corr(iobs,1,iy), &
                  !                       '; res=', res, self%andom%y(iobs), self%andom%v(iobs,j), self%andom%HX(iobs,j); call goPr
                  !    end if
                  !  end do
                  !end if
                  !
                end do ! ihist
              end if  ! analyse
              !print *, 'xxx2   an ', minval(Ens(j)%dc(:,:,inoise,:)), maxval(Ens(j)%dc(:,:,inoise,:))
            end do ! inoise
          !~
          case default
            write (gol,'("unsupported substate ",i0)') substate; call goErr
            TRACEBACK; status=1; return
        end select
       end do ! obs
    end do ! modes

    !...............................................
    ! end timing:
    call GO_Timer_End( itim_update_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! ~
    
    ! clear:
    deallocate( HSST, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Kdt, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( Ed_c, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Sd_c, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( xd_c, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Kd_c, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_EnKF_Apply_Domain_SubState
  
  
  ! ***


  !
  ! Local Ensemble Transform KF
  !
  ! See doc/build/html/index.html for description.
  !

  subroutine LEKF_Meas_Update_LETKF_Apply( self, t1, t2, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use LE_Grid   , only : ugg
    use Dims      , only : nx, ny
    use Indices   , only : nspec, specname
    use LEKF_State, only : nmodes
    use LEKF_Noise, only : nnoise, noise_name
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_LETKF_Apply'

    ! --- local ----------------------------------
    
    integer     ::  i, j
    integer     ::  ispec
    integer     ::  inoise
    logical     ::  debug
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a," - apply measurement update ...")') rname; call goPr
    
    ! storage:
    allocate( self%xi(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%xi_a(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%xx(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ...
    write (gol,'(a," -   collect observations, assign to domains ...")') rname; call goPr

    ! start timing:
    call GO_Timer_Start( itim_obs_collect, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Retrieve all observations y, repr. err. stdvs r, 
    ! and simulations for local modes ;
    ! results are stored into analysis domain structures:
    call self%CollectObservations( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
    
    ! end timing:
    call GO_Timer_End( itim_obs_collect, status )
    IF_NOTOK_RETURN(status=1)

    ! start timing:
    call GO_Timer_Start( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! any obs in this domain?
    if ( self%andom%nobs > 0 ) then

      ! info ...
      write (gol,'(a," -     analyse concentrations ...")')  rname; call goPr
      
      ! loop over grid cells:
      do j = 1, ny
        do i = 1, nx
        
          ! loop over species:
          do ispec = 1, nspec
        
            ! collect observations within range of this cell,
            ! and that should be used to update these concentrations
            ! use these to fill andom%Y
            call self%andom%Fill_Y_etc( ugg%longitude(i,j), ugg%latitude(i,j), self%xcorr, status, ispec=ispec )
            IF_NOTOK_RETURN(status=1)

            ! any obs selected?
            if ( self%andom%nobsl > 0 ) then
              !! testing ..
              !write (gol,'(a," -       cell (",i3,",",i3,") using ",i4," obs for ",a)') &
              !                    rname, i, j, self%andom%nobsl, trim(specname(ispec)); call goPr
              ! update ensemble:
              call self%LETKF_Update_Ensemble( i, j, self%andom%wa, self%andom%WWa, status, ispec=ispec )
              IF_NOTOK_RETURN(status=1)
              
              !! testing ..
              !write (gol,'("break after first analysis")'); call goErr
              !TRACEBACK; status=1; return
            end if ! obs selected

          end do  ! ispec
  
        end do ! i
      end do ! j

      ! info ...
      write (gol,'(a," -     analyse aerosol water content ...")')  rname; call goPr
      
      ! loop over grid cells:
      do j = 1, ny
        do i = 1, nx
        
          ! collect observations within range of this cell,
          ! and that should be used to update these concentrations
          ! use these to fill andom%Y
          call self%andom%Fill_Y_etc( ugg%longitude(i,j), ugg%latitude(i,j), self%xcorr, status, ispec=ispec )
          IF_NOTOK_RETURN(status=1)

          ! any obs selected?
          if ( self%andom%nobsl > 0 ) then
            !! testing ..
            !write (gol,'(a," -       cell (",i3,",",i3,") using ",i4," obs")') &
            !                    rname, i, j, self%andom%nobsl; call goPr
            ! update ensemble:
            call self%LETKF_Update_Ensemble( i, j, self%andom%wa, self%andom%WWa, status, sia=.true. )
            IF_NOTOK_RETURN(status=1)
          end if ! obs selected

        end do ! i
      end do ! j

      ! info ...
      write (gol,'(a," -     analyse perturbation factors ...")')  rname; call goPr
      
      ! loop over grid cells:
      do j = 1, ny
        do i = 1, nx
        
          ! loop over noise elements:
          do inoise = 1, nnoise
            
            ! testing ...
            debug = inoise == 1
        
            ! collect observations within range of this cell,
            ! and that should be used to update these concentrations
            ! use these to fill andom%Y
            call self%andom%Fill_Y_etc( ugg%longitude(i,j), ugg%latitude(i,j), self%xcorr, status, &
                                           inoise=inoise )!, debug=debug )
            IF_NOTOK_RETURN(status=1)

            ! any obs selected?
            if ( self%andom%nobsl > 0 ) then
              ! testing ..
              if ( debug ) then
                write (gol,'(a," -       cell (",i3,",",i3,") using ",i4," obs for ",a)') &
                                     rname, i, j, self%andom%nobsl, trim(noise_name(inoise)); call goPr
              end if
              ! update ensemble:
              call self%LETKF_Update_Ensemble( i, j, self%andom%wa, self%andom%WWa, status, inoise=inoise )
              IF_NOTOK_RETURN(status=1)
            end if ! obs selected

          end do  ! ispec
  
        end do ! i
      end do ! j

      ! clear current arrays
      call self%andom%Clear_Y_etc( status )
      IF_NOTOK_RETURN(status=1)

    else
      ! info ...
      write (gol,'(a," -     no observations in domain ...")')  rname; call goPr
    end if

    ! clear:
    deallocate( self%xi, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%xi_a, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%xx, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! end timing:
    call GO_Timer_End( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_LETKF_Apply
  
  
  ! *


  !
  ! Update ensemble members using LETKF weights.
  !

  subroutine LEKF_Meas_Update_LETKF_Update_Ensemble( self, i, j, wa, WWa, status, &
                                                       ispec, sia, inoise )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use LEKF_State, only : nmodes
    use LEKF_State, only : Ens

    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    integer, intent(in)                   ::  i, j     ! grid cell indices
    real, intent(in)                      ::  wa(:)     ! (nmodes) analysis weights for mean
    real, intent(in)                      ::  WWa(:,:)  ! (nmodes) analysis weights for perturbations
    integer, intent(out)                  ::  status
    
    integer, intent(in), optional         ::  ispec   ! update state concentrations for species
    logical, intent(in), optional         ::  sia     ! update state aerosol water
    integer, intent(in), optional         ::  inoise  ! update state noise elements

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_LETKF_Update_Ensemble'

    ! --- local ----------------------------------
    
    integer                   ::  nk
    integer                   ::  k
    integer                   ::  imode

    ! --- begin ----------------------------------

    !...............................................
    ! start timing:
    call GO_Timer_Start( itim_update_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! check that [1;..;1] is an eigenvector of W:
    ! if ( any( sum(WWa,dim=2) - 1.0 > 0.001 ) ) then
    ! if ( any( ABS(sum(WWa,dim=2) - 1.0) > 0.02 ) ) then
    if ( any( ABS(sum(WWa,dim=2) - 1.0) > 0.1 ) ) then
      !write (gol,'("sum over columns in W should be 1:")'); call goPr
      self%n_w_test_error = self%n_w_test_error + 1
    else
      self%n_w_test_pass = self%n_w_test_pass + 1
    end if
    !  write (gol,'("sum over columns in W should be 1:")'); call goErr
    !  do imode = 1, nmodes
    !    write (gol,*) WWa(imode,:),' ; sum;', sum(WWa(imode,:)); call goErr
    !  end do
    !  write (gol,*), self%andom%lambda,';sum;','lambda'; call goErr
      !write (gol,*), self%andom%QQ,';','QQ'; call goErr
      !write (gol,*), self%andom%invPa,';','invPa'; call goErr
      
    !   write (gol,*), 'y0 ',maxval(self%andom%yo(1:self%andom%nobs)),';','max',minval(self%andom%yo(1:self%andom%nobs)),';','min'; call goErr
    !   write (gol,*), 'ybar ',maxval(self%andom%ybar(1:self%andom%nobs)),';','max',minval(self%andom%ybar(1:self%andom%nobs)),';','min'; call goErr
    !   write (gol,*), 'hx ',maxval(self%andom%HX(1:self%andom%nobs,1:nmodes)),';','max',minval(self%andom%HX(1:self%andom%nobs,1:nmodes)),';','min'; call goErr
    !   write (gol,*), 'YY ',maxval(self%andom%YY(1:self%andom%nobs,1:nmodes)),';','max',minval(self%andom%YY(1:self%andom%nobs,1:nmodes)),';','min'; call goErr
    !   write (gol,*), 'CC ',maxval(self%andom%CC(1:self%andom%nobs,1:nmodes)),';','max',minval(self%andom%CC(1:self%andom%nobs,1:nmodes)),';','min'; call goErr
       
    !  TRACEBACK; status=1; return
    !end if
    
    ! update tracers?
    if ( present(ispec) ) then
    
      ! ~ concentrations:
      !     c(nx,ny,nz,nspec)
      nk = size( Ens(1)%c, 3 )
      ! loop over dimensions:
      do k = 1, nk
        ! extract ensemble values:
        do imode = 1, nmodes
          self%xi(imode) = Ens(imode)%c(i,j,k,ispec)
        end do
        ! obtain analyzed values:
        call self%LETKF_Update_Cell( self%xi, wa, WWa, self%xi_a, self%xx, status )
        IF_NOTOK_RETURN(status=1)
        ! reset:
        do imode = 1, nmodes
          Ens(imode)%c(i,j,k,ispec) = self%xi_a(imode)
        end do
        !! testing ...
        !exit
      end do ! k

      ! ~ ground concentrations:
      !     cg(nx,ny,nspec)
      ! extract ensemble values:
      do imode = 1, nmodes
        self%xi(imode) = Ens(imode)%cg(i,j,ispec)
      end do
      ! obtain analyzed values:
      call self%LETKF_Update_Cell( self%xi, wa, WWa, self%xi_a, self%xx, status )
      IF_NOTOK_RETURN(status=1)
      ! reset:
      do imode = 1, nmodes
        Ens(imode)%cg(i,j,ispec) = self%xi_a(imode)
      end do
      
    end if ! ispec
    
    ! sia component?
    if ( present(sia) .and. sia ) then

      ! ~ aerosol water content:
      !      aerh2o(nx,ny,nz)
      nk = size( Ens(1)%aerh2o, 3 )
      ! loop over dimensions:
      do k = 1, nk
        ! extract ensemble values:
        do imode = 1, nmodes
          self%xi(imode) = Ens(imode)%aerh2o(i,j,k)
        end do
        ! obtain analyzed values:
        call self%LETKF_Update_Cell( self%xi, wa, WWa, self%xi_a, self%xx, status )
        IF_NOTOK_RETURN(status=1)
        ! reset:
        do imode = 1, nmodes
          Ens(imode)%aerh2o(i,j,k) = self%xi_a(imode)
        end do
      end do ! k
      
    end if ! sia
    
    ! noise components?
    if ( present(inoise) ) then

      ! ~ perturbation factors:
      !      dc(nx,ny,nnoise,nhist)
      nk = size( Ens(1)%dc, 4 )
      ! loop over dimensions:
      do k = 1, nk
       ! extract ensemble values:
        do imode = 1, nmodes
          self%xi(imode) = Ens(imode)%dc(i,j,inoise,k)
        end do
        ! obtain analyzed values:
        call self%LETKF_Update_Cell( self%xi, wa, WWa, self%xi_a, self%xx, status )
        IF_NOTOK_RETURN(status=1)
        ! reset:
        do imode = 1, nmodes
          Ens(imode)%dc(i,j,inoise,k) = self%xi_a(imode)
        end do
      end do ! k
      
    end if ! noise

    !...............................................
    ! end timing:
    call GO_Timer_End( itim_update_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_LETKF_Update_Ensemble
  
  
  ! *


  !
  ! Update ensemble elements using LETKF weights.
  !
  ! Mean and perturbations:
  !     x = 1/m sum_i x_i
  !     X = [ .., x_i - x, .. ]
  !
  ! The sum over the columns of X is zero:
  !     X [1;..;1] = sum_i (x_i - 1/m sum_i x_i) = 0
  !
  ! The transformation matrix:
  !     W = [ .., w_i, .. ]
  ! has property:
  !     W [1;..;1] = [1;..;1]
  ! Therefore:
  !     X W [1;..;1] = X [1;..;1] = 0
  ! thus:
  !     sum_i X w_i = 0
  !
  ! The weights 'wa' define the update of the mean:
  !     xa = x + X wa
  ! Update of ensemble:
  !     xa_i = xa + X wa_i
  !          = x + X ( wa + wa_i )
  !
  ! The new ensemble mean is indeed:
  !     1/m sum_i xa_i = xa + 1/m sum_i X wa_i = xa + 0
  !
  

  subroutine LEKF_Meas_Update_LETKF_Update_Cell( self, xi, wa, WWa, xi_a, xx, status )

    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    real, intent(in)                      ::  xi(:)       ! (nmodes) forecast ensemble values
    real, intent(in)                      ::  wa(:)       ! (nmodes) analysis weights for mean
    real, intent(in)                      ::  WWa(:,:)    ! (nmodes) analysis weights for perturbations
    real, intent(out)                     ::  xi_a(:)     ! (nmodes) analysis ensemble values
    real, intent(inout)                   ::  xx(:)       ! (nmodes) work array
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_LETKF_Update_Cell'

    ! --- local ----------------------------------
    
    real          ::  xf
    real          ::  xa

    ! --- begin ----------------------------------
    
    ! ensemble mean:
    xf = sum(xi) / size(xi)
    ! perturbations:
    xx = xi - xf
    ! analyzed mean:
    xa = xf + dot_product( xx, wa )
    ! analyzed ensemble members:
    xi_a = xa + matmul( xx, WWa )
    
    !! testing ..
    !write (gol,*) 'xi ; <xi> = '; call goPr
    !write (gol,*) xi, ' ; ', xf; call goPr
    !write (gol,*) 'xx ; sum(xx) = '; call goPr
    !write (gol,*) xx, ' ; ', sum(xx); call goPr
    !write (gol,*) 'xi_a ; <xi_a>, xa = '; call goPr
    !write (gol,*) xi_a, ' ; ', sum(xi_a)/size(xi_a), xa; call goPr
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_LETKF_Update_Cell



end module LEKF_Meas_Update

