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
#include "lekf.inc"
!
!###############################################################################


module LEKF_Meas_Update

  use GO, only : gol, goPr, goErr

  use Grid, only : TllGridInfo

  implicit none
  
  
  ! --- in/out -----------------------------

  private

  public    ::  T_Meas_Update
  public    ::  LEKF_Meas_Update_Timers
  

  ! --- const --------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_Update'
  
  ! truncate correlations after frac * rho:
  real, parameter               ::  rho_frac_zero = 3.5


  ! --- types --------------------------------
  
  ! observation inidces:
  type T_ObsInd
    ! maori obs identified by set, location, and observerd variable:
    integer                       ::  iset
    integer                       ::  iloc
    integer                       ::  obs_ivar
    ! inside domain?
    logical                       ::  inside
    ! representative location:
    real                          ::  lon, lat
  end type T_ObsInd
  
  ! observation selection:
  type T_AnDom
!    ! analysis domain:
!    type(TllGridInfo)             ::  lli
!    ! offset in global domain:
!    integer                       ::  i0, j0
    ! selected nearby observations:
    integer                       ::  nobs
    ! reference to observation data:
    type(T_ObsInd), pointer       ::  inds(:)  ! (>= n)
    ! observation and error stdv:
    real, pointer                 ::  y(:)
    real, pointer                 ::  r(:)
    real, pointer                 ::  v(:,:)    ! (nobs,nmodes)
    real, pointer                 ::  HX(:,:)   ! (nobs,nmodes)
    ! correlation length scale for locaization,
    ! only one value supported yet:
    real                          ::  rho_m   ! m
    ! spatial localization weight:
    real, pointer                 ::  corr(:,:,:)  ! (nobs,nlon,nlat)
    ! help array with distances:
    real, pointer                 ::  dist(:,:)  ! (nlon,nlat)
    ! analysis matrices:
    real, allocatable             ::  HSd(:,:)   ! (nobs,nmodes)
    real, allocatable             ::  HPHR(:,:)  ! (nobs,nobs)
    real, allocatable             ::  U(:,:)     ! (nobs,nobs)
  contains
    procedure ::  Init                  =>  AnDom_Init
    procedure ::  Done                  =>  AnDom_Done
    procedure ::  Reset                 =>  AnDom_Reset
    procedure ::  Add_Observation       =>  AnDom_Add_Observation
    procedure ::  Exchange_Nearby       =>  AnDom_Exchange_Nearby
    !procedure ::  BCast_Modes           =>  AnDom_BCast_Modes
    procedure ::  Fill_HPHR_etc         =>  AnDom_Fill_HPHR_etc
    procedure ::  Clear_HPHR_etc        =>  AnDom_Clear_HPHR_etc
  end type T_AnDom

  ! *
  
  ! chem and other correlation info:
  type T_XCorr
    ! decorrelation between obs and species:
    logical, allocatable              ::  analyse_spec(:)  ! (nspec_all)
    logical                           ::  analyse_sia
    ! decorrelation between obs and noise factors:
    logical, allocatable              ::  analyse_noise(:)  ! (nnoise)
  end type T_XCorr

  ! *
  
  type T_Meas_Update
    ! number of analysis domains:
    integer                                ::  ndom
    !! analysis domains:
    !type(T_AnDom), pointer                 ::  dom(:)  ! (ndom)
    ! analysis domain (subdomain and area around):
    type(T_AnDom)                          ::  andom
    ! maori info (chemical and noise correlations):
    type(T_XCorr), allocatable             ::  xcorr(:)  ! (nset)
  contains
    procedure ::  Init                  =>  LEKF_Meas_Update_Init
    procedure ::  Done                  =>  LEKF_Meas_Update_Done
    procedure ::  Apply                 =>  LEKF_Meas_Update_Apply
    procedure ::  CollectObservations   =>  LEKF_Meas_Update_CollectObservations
    procedure ::  Apply_Domain          =>  LEKF_Meas_Update_Apply_Domain
    procedure ::  Apply_Domain_Substate =>  LEKF_Meas_Update_Apply_Domain_Substate
  end type T_Meas_Update


  ! --- interfaces ------------------------
  
  interface HCorrGaussian
    module procedure HCorrGaussian_0d
    module procedure HCorrGaussian_2d
  end interface HCorrGaussian


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
  ! *** module
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
#ifdef with_kf_meas_omi_trc
    use LEKF_Meas_OMI_Trc, only : omi_trc
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
    
    integer               ::  nset_maori
    integer               ::  nset_omi
    integer               ::  nset, iset
    character(len=64)     ::  sname
    logical               ::  analyse
    character(len=1024)   ::  key

    ! --- begin ---------------------------

    ! info ...
    write (gol,'(a," - init analysis domain ...")') rname; call goPr

!    ! read distribution of sub-domains for analysis:
!    call rcF%Get( 'lekf.meas.andom.nx', ndom_x, status )
!    IF_NOTOK_RETURN(status=1)
!    call rcF%Get( 'lekf.meas.andom.ny', ndom_y, status )
!    IF_NOTOK_RETURN(status=1)
!    ! total number:
!    self%ndom = ndom_x * ndom_y
!      
!    ! allocate analysis domains:
!    allocate( self%dom(self%ndom), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! init index:
!    idom = 0
!    ! loop over domain grid:
!    do idom_y = 1, ndom_y
!      do idom_x = 1, ndom_x
!        ! increase counter:
!        idom = idom + 1
!        ! cells:
!        nx_loc = lli%nlon / ndom_x
!        ny_loc = lli%nlat / ndom_y
!        ! info ...
!        write (gol,'(a," - an domain ",i0," (",i0,",",i0,") shape (",i0,",",i0,")")') &
!                       rname, idom, idom_x, idom_y, nx_loc, ny_loc; call goPr
!        ! check ...
!        if ( modulo(lli%nlon,nx_loc) /= 0 ) then
!          write (gol,'(a," - number of x analysis domain ",i0," should match with grid size ",i0)') &
!                               rname, ndom_x, lli%nlon; call goErr
!          TRACEBACK; status=1; return
!        end if
!        ! check ...
!        if ( modulo(lli%nlat,ny_loc) /= 0 ) then
!          write (gol,'(a," - number of y analysis domain ",i0," should match with grid size ",i0)') &
!                               rname, ndom_y, lli%nlat; call goErr
!          TRACEBACK; status=1; return
!        end if
!        ! ofset:
!        i0 = (idom_x-1)*nx_loc
!        j0 = (idom_y-1)*ny_loc
!        ! init domain:
!        call self%dom(idom)%Init( lli%lon_deg(1)+i0*lli%dlon_deg, lli%dlon_deg, nx_loc, i0, &
!                                     lli%lat_deg(1)+j0*lli%dlat_deg, lli%dlat_deg, ny_loc, j0, status )
!        IF_NOTOK_RETURN(status=1)
!        ! info ..
!        write (gol,'(a," -   bounding box : [",f8.2,",",f8.2,"] x [",f8.2,",",f8.2,"]")') rname, &
!                        self%dom(idom)%lli%blon_deg(0), self%dom(idom)%lli%blon_deg(nx_loc), &
!                        self%dom(idom)%lli%blat_deg(0), self%dom(idom)%lli%blat_deg(ny_loc); call goPr
!      end do ! ix
!    end do ! iy

    ! init storage for analysis observations for this domain
    call self%andom%Init( status )
    IF_NOTOK_RETURN(status=1)

    ! ***
    
    ! info ...
    write (gol,'(a," - extract chemical/noise correlation info ...")') rname; call goPr
    
#ifdef with_kf_meas_maori
    ! get number of sets to put out:
    call MAORI_Data_Inquire( mad, status, nset=nset_maori )
    IF_NOTOK_RETURN(status=1)
#else
    nset_maori = 0
#endif

    ! OMI observations?
#ifdef with_kf_meas_omi_trc
    nset_omi = 1
#else
    nset_omi = 0
#endif

    ! info ...
    write (gol,'(a," -   number of maori sets : ",i2)') rname, nset_maori; call goPr
    write (gol,'(a," -   number of omi   sets : ",i2)') rname, nset_omi  ; call goPr
    ! total:
    nset = nset_maori + nset_omi

    ! any?
    if ( nset > 0 ) then
      ! storage:
      allocate( self%xcorr(nset), stat=status )
      IF_NOTOK_RETURN(status=1)
      
#ifdef with_kf_meas_maori
      ! ~ maori

      ! loop over maori sets:
      do iset = 1, nset_maori

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
        call rcF%Get( 'maori.'//trim(sname)//'.assim.spec', key, status )
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

      ! ~ OMI
      
#ifdef with_kf_meas_omi_trc
      ! OMI set:
      iset = nset_maori + 1
      ! analyse?
      if ( omi_trc%analyse ) then

        ! read noise correlation keys:
        call rcF%Get( 'kf.meas.omi_trc.spec', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_spec(nspec_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetSpecApply( key, self%xcorr(iset)%analyse_spec, self%xcorr(iset)%analyse_sia, status )
        IF_NOTOK_RETURN(status=1)

        ! read noise correlation keys:
        call rcF%Get( 'kf.meas.omi_trc.noise', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_noise(nnoise), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetNoiseApply( key, self%xcorr(iset)%analyse_noise, status )
        IF_NOTOK_RETURN(status=1)
      
      end if ! analyse?
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

!    ! loop:
!    do idom = 1, self%ndom
!      ! done:
!      call self%dom(idom)%Done( status )
!      IF_NOTOK_RETURN(status=1)
!    end do
!    ! clear:
!    deallocate( self%dom, stat=status )
!    IF_NOTOK_RETURN(status=1)
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


  !
  ! Measurement update.
  !
  ! Retrieve all observations y, repr. err. stdvs r, and simulations for local modes:
  !   yy, rr, H X(1:ml)
  ! Distribute over analysis domains for which an observation is relevant:
  !   yd, rd, Hd X(1:ml)
  ! Broadcast simulations to all pe:
  !   Hd X(1:m)
  !
  ! For each analysis domain:
  !
  !   Broadcast modes:
  !     Xd(1:m)
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

  subroutine LEKF_Meas_Update_Apply( self, t1, t2, status )
  
    use GO, only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO, only : TDate
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply'

    ! --- local ----------------------------------
    
!    integer     ::  idom

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

!    ! info ...
!    write (gol,'(a," -   broadcast ensemble simulations ...")') rname; call goPr
!
!    ! start timing:
!    call GO_Timer_Start( itim_obs_bcast, status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ensure that all nodes have the current simulations for all modes,
!    ! not influenced by upcoming measurement updates ...
!    ! loop over analysis domains:
!    do idom = 1, self%ndom
!      ! broadcast simulations for local ensemble members to other nodes:
!      call self%dom(idom)%BCast_Modes( status )
!      IF_NOTOK_RETURN(status=1)
!    end do ! domains
!
!    ! switch timing:
!    call GO_Timer_Switch( itim_obs_bcast, itim_domains, status )
!    IF_NOTOK_RETURN(status=1)
!   
!    ! info ...
!    write (gol,'(a," -   loop over domains ...")') rname; call goPr
!    ! measurement update per domain:
!    do idom = 1, self%ndom
!      ! any obs?
!      if ( self%dom(idom)%nobs > 0 ) then
!        ! info ...
!        write (gol,'(a," -     analyse domain ",i0," ...")')  rname, idom; call goPr
!        ! apply:
!        call self%Apply_Domain( idom, status )
!        IF_NOTOK_RETURN(status=1)
!      else
!        ! info ...
!        write (gol,'(a," -     no observations in domain ",i0," ...")')  rname, idom; call goPr
!      end if
!    end do

    ! start timing:
    call GO_Timer_Start( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! any obs?
    if ( self%andom%nobs > 0 ) then
      ! info ...
      write (gol,'(a," -     analyse domain ...")')  rname; call goPr
      ! apply:
      call self%Apply_Domain( status )
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

  end subroutine LEKF_Meas_Update_Apply
  
  
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
#ifdef with_kf_meas_omi_trc
    use LEKF_State, only : substate_omi_trc
    use LEKF_Meas_OMI_TRC, only : omi_trc
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
    real                        ::  lon, lat
    integer                     ::  obs_nvar, obs_ivar
    real                        ::  y, r
    real                        ::  alfa
    integer                     ::  idom
    logical                     ::  nearby
    logical                     ::  inside
    integer                     ::  j
    !real, allocatable           ::   vn_all(:,:)  ! (loc_id_range,nmodes)
    real, allocatable           ::   v_loc(:)  ! (nmodes)
    real, allocatable           ::  HX_loc(:)  ! (nmodes)

    integer                     ::  nmeas_tot
    integer                     ::  nmeas_intimestep_tot
    integer                     ::  ipix
    type(TDate)                 ::  t
    integer                     ::  inds(0:7)
    
    !! testing ...
    !real, allocatable   ::  vv(:,:)  ! (nmodes_loc,npix)
    !integer             ::  ipix1, ipix2, ipixstep
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a," - start")') rname; call goPr

!    ! loop over analysis domains:
!    do idom = 1, self%ndom
!      ! reset to empty:
!      call self%dom(idom)%Reset( status )
!      IF_NOTOK_RETURN(status=1)
!    end do ! idom
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
    allocate(  v_loc(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( HX_loc(nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)

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
    write (gol,'(a," - loop over ",i0," maori sets ...")') rname, nset_maori; call goPr

#ifdef with_kf_meas_maori
    ! loop over maori sets:
    do iset = 1, nset_maori

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
                ! random error out N(0,r)
                call Ens(j)%rnd%Get_Normal( v_loc(j), status, sigma=r )
                IF_NOTOK_RETURN(status=1)
                !! use random number out of N(0,1) that is independent of (sub)domain:
                !v_loc(j) = vn_all(loc_id,j) * r
              end do
              
              !! testing ...
              !write (gol,*) 'vvv loc_id=', loc_id, ' v=', vn_all(loc_id,:); call goPr

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
                                                         y, r, v_loc, HX_loc, &
                                                         status )
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

          !! storage?
          !if ( allocated(vn_all) ) then
          !  ! clear:
          !  deallocate(  vn_all, stat=status )
          !  IF_NOTOK_RETURN(status=1)
          !end if

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

#ifdef with_kf_meas_omi_trc
    ! info ..
    write (gol,'(a," - collect OMI observations ...")') rname; call goPr
    ! set index:
    iset = nset_maori + 1
    ! to be analysed?
    if ( omi_trc%analyse ) then
    
      ! info ...
      write (gol,'(a," -   OMI observations to be analyzed ...")') rname; call goPr

      ! total number over all processors:
      nmeas_tot = omi_trc%nmeas
      call goc%AllReduce( 'sum', nmeas_tot, status )
      IF_NOTOK_RETURN(status=1)
      ! total number of measurments in current timestep:
      nmeas_intimestep_tot = omi_trc%nmeas_intimestep
      call goc%AllReduce( 'sum',  nmeas_intimestep_tot, status )
      IF_NOTOK_RETURN(status=1)
 
      ! any measurements?
      if ( (nmeas_tot > 0) .and. (nmeas_intimestep_tot > 0) ) then
    
        ! info ...
        write (gol,'(a," -   number of pixels for this timestep : ",i0)') rname, omi_trc%nmeas; call goPr
          
!        ! local range of unique pixel id's:
!        loc_id_range(1) = minval(omi_trc%pixel_id(1:omi_trc%nmeas))
!        loc_id_range(2) = maxval(omi_trc%pixel_id(1:omi_trc%nmeas))
!        ! global location id range:
!        loc_id_range_all = loc_id_range
!        call goc%AllReduce( 'min', loc_id_range_all(1), status )
!        IF_NOTOK_RETURN(status=1)
!        call goc%AllReduce( 'max', loc_id_range_all(2), status )
!        IF_NOTOK_RETURN(status=1)
!        ! info ..
!        write (gol,'(a," -     location id range : ",i0,":",i0," of ",i0,":",i0)') &
!                          rname, loc_id_range, loc_id_range_all; call goPr
!        ! check ...
!        if ( loc_id_range_all(1) > loc_id_range_all(2) ) then
!          write (gol,'("strange loc_id range ",i0,":",i0)') loc_id_range_all; call goPr
!          TRACEBACK; status=1; return
!        end if
!        ! storage for random represenation error samples:
!        allocate(  vn_all(loc_id_range_all(1):loc_id_range_all(2),nmodes), stat=status )
!        IF_NOTOK_RETURN(status=1)
!        ! loop over modes:
!        do j = 1, nmodes
!          ! random error out N(0,1)
!          call Ens(j)%rnd%Get_Normal( vn_all(:,j), status )
!          IF_NOTOK_RETURN(status=1)
!        end do ! modes

        ! convert from km to m:
        rho_m = omi_trc%rho * 1.0e3  ! m
        
        
        !! ADHOC FOR TESTING REVERSE ORDER ...
        !allocate( vv(nmodes_loc,omi_trc%nmeas), stat=status )
        !IF_NOTOK_RETURN(status=1)
        !! loop over pixels:
        !do ipix = 1, omi_trc%nmeas
        !  ! first check if this obs needs processing
        !  if ( .not. omi_trc%intimestep(ipix) ) cycle
        !  ! extract repr.err. standard deviation:
        !  r    = omi_trc%sigma_vcd_trop(ipix)
        !  ! loop over modes:
        !  do j = 1, nmodes_loc
        !    ! generate random error with std.dev. for this pixel:
        !    call Ens(j)%rnd%Get_Normal( vv(j,ipix), status, sigma=r )
        !    IF_NOTOK_RETURN(status=1)
        !  end do ! local mode
        !end do ! ipix
    
        !! pixel order:
        !select case ( omi_trc%analyse_order )
        !  case ( 'normal' )
        !    ipix1    = 1
        !    ipix2    = omi_trc%nmeas
        !    ipixstep = 1
        !  case ( 'reverse' )
        !    ipix1    = omi_trc%nmeas
        !    ipix2    = 1
        !    ipixstep = -1
        !  case default
        !    write (gol,'("could not set pixel range from analyse_order `",a,"`")') trim(omi_trc%analyse_order); call goErr
        !    TRACEBACK; status=1; return
        !end select
        !! loop over measurements; testing different orders ...
        !do ipix = ipix1, ipix2, ipixstep

        ! loop over measurements assigned to this domain:
        do ipix = 1, omi_trc%nmeas

          ! first check if this obs needs processing
          if ( .not. omi_trc%intimestep(ipix) ) cycle

          ! Doube-check if outside (t1,t2] :
          t = NewDate( time6=omi_trc%date(:,ipix) )
          if ( (t <= t1) .or. (t2 < t) ) then
             write (gol,'("omi trc measurement not in time interval:")'); call goErr
             call wrtgol( '  model time interval  : ', t1, ' - ', t2 ); call goErr
             call wrtgol( '  pixel time           : ', t ); call goErr
             TRACEBACK; status=1; return         
          end if

          ! double-check if not accepted (validation, no data?):
          if ( omi_trc%status(ipix) /= 0 ) then
             write (gol,'("omi trc observation status invalid: ",i0)') omi_trc%status(ipix); call goErr
             TRACEBACK; status=1; return         
          end if

          ! check to see if obs has already been analysed ...
          if ( omi_trc%analysed(ipix) ) then
             write (gol,'("WARNING omi trc measurement seems to be analysed already:")'); call goErr
             write (gol,'("  measurement index    : ",i6)') ipix; call goErr
             call wrtgol( '  measurement time     : ', t ); call goErr
             call wrtgol( '  model time interval  : ', t1, ' - ', t2 ); call goErr
             TRACEBACK; status=1; return
          endif

          !! info ...
          !write (gol,'(a," -   accepted : ",i6,2l3)') ipix, meas_omi(i)%analyse, meas_omi(i)%accepted; call goPr
          
          ! set flags:
          omi_trc%screened(ipix) = .false.
          omi_trc%analysed(ipix) = .true.
          
          ! extract measurement:
          lon  = omi_trc%longitude(ipix)
          lat  = omi_trc%latitude(ipix)
          y    = omi_trc%vcd_trop(ipix)
          r    = omi_trc%sigma_vcd_trop(ipix)
          alfa = omi_trc%screening_factor

          ! define how to extract value from state:
          inds(0:1) = (/ substate_omi_trc, ipix /)
          
          ! location id:
          loc_id = omi_trc%pixel_id(ipix)

          ! loop over modes:
          do j = 1, nmodes
            ! extract simulation of measurements:
            call Ens(j)%GetValue( inds, HX_loc(j), status )
            IF_NOTOK_RETURN(status=1)
            ! random error out N(0,r)
            call Ens(j)%rnd%Get_Normal( v_loc(j), status, sigma=r )
            IF_NOTOK_RETURN(status=1)
            !! testing reversed order:
            !v_loc(j) = vv(j,ipix)
            !! use random number out of N(0,1) that is independent of (sub)domain:
            !v_loc(j) = vn_all(loc_id,j) * r
          end do
              
          !! testing ...
          !write (gol,*) 'vvv loc_id=', loc_id, ' v=', vn_all(loc_id,:); call goPr

          !! loop over analysis domains:
          !do idom = 1, self%ndom

          !! within factor*rho [m] of domain?
          !call ugg%WithinDistance( lon, lat, rho_frac_zero * rho_m, nearby, status )
          !IF_NOTOK_RETURN(status=1)

          !! nearby observation? then add to list:
          !if ( nearby ) then

            ! TESTING ...
            ! check if in domain (not in surrounding area):
            call ugg%InDomain( lon, lat, inside, status )
            IF_NOTOK_RETURN(status=1)
            ! expected to be inside ..
            if ( .not. inside ) then
              write (gol,'("fould location outside domain")'); call goErr
              TRACEBACK; status=1; return
            end if

            ! add to domain, identify by set number (used for chem/noise correlations)
            ! and pixel number:
            call self%andom%Add_Observation( (/iset,ipix,-999/), inside, &
                                                     lon, lat, rho_m, &
                                                     y, r, v_loc, HX_loc, &
                                                     status, debug=.true. )
            IF_NOTOK_RETURN(status=1)

            !! info ...
            !write (gol,'(a," -   added OMI pixel (index ",i0,") at (",f8.2,",",f8.2,") to domain ",i0)') &
            !                        rname, ipix, lon, lat, idom; call goPr

          !end if ! nearby domain

          !end do  ! analysis domains

        end do ! pixels
        
        !! clear:
        !deallocate( vv, stat=status )
        !IF_NOTOK_RETURN(status=1)

        !! clear:
        !deallocate(  vn_all, stat=status )
        !IF_NOTOK_RETURN(status=1)

      else
      
        ! info ...
        write (gol,'(a," -   no pixels for this timestep ...")') rname; call goPr

      end if

    else

      ! info ...
      write (gol,'(a," -   OMI observations not to be analyzed ...")') rname; call goPr

    end if
#endif

    ! end start timing:
    call GO_Timer_End( itim_obs_collect_sat, status )
    IF_NOTOK_RETURN(status=1)

    ! *
        
    ! clear:
    deallocate( HX_loc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( v_loc, stat=status )
    IF_NOTOK_RETURN(status=1)
    
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
  ! For each analysis domain:
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

  subroutine LEKF_Meas_Update_Apply_Domain( self, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use Dims      , only : nspec
    use LEKF_State, only : substate_c, substate_cg, substate_aerh2o, substate_dc
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    !integer, intent(in)                   ::  idom
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply_Domain'

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
    
    ! call for specific elements;
    ! update concentration arrays per spec to limit memory usage:
    do ispec = 1, nspec
      call self%Apply_Domain_SubState( substate_c     , ispec, status )
      IF_NOTOK_RETURN(status=1)
      call self%Apply_Domain_SubState( substate_cg    , ispec, status )
      IF_NOTOK_RETURN(status=1)
    end do
    ! update full arrays, pass dummy ispec:
    call self%Apply_Domain_SubState( substate_aerh2o, -999, status )
    IF_NOTOK_RETURN(status=1)
    call self%Apply_Domain_SubState( substate_dc    , -999, status )
    IF_NOTOK_RETURN(status=1)

    ! clear current arrays
    call self%andom%Clear_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Apply_Domain
  
  
  ! *


  subroutine LEKF_Meas_Update_Apply_Domain_SubState( self, substate, ispec, status )
  
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
  
      !use file_nc

    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    !integer, intent(in)                   ::  idom
    integer, intent(in)                   ::  substate
    integer, intent(in)                   ::  ispec
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply_Domain_SubState'

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

    ! --- begin ----------------------------------
    
    !! info ...
    !write (gol,'(a," -   state element `",a,"` ...")') rname, trim(substates(substate)); call goPr

    !! short:
    !dom => self%dom(idom)
    
!    ! domain size:
!    nxd = dom%lli%nlon
!    nyd = dom%lli%nlat
!    ! local range:
!    i1 = dom%i0+1 ; i2 = dom%i0+nxd
!    j1 = dom%j0+1 ; j2 = dom%j0+nyd
    
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

    ! solve using upper-triangular factorization:
    call LinAlg_Sym_FactorSolve( self%andom%U, HSSt, Kdt, status )
    IF_NOTOK_RETURN(status=1)
    
    ! reshape:
    Kd_c = reshape( Kdt, (/self%andom%nobs,shp(1),shp(2),shp(3),shp(4)/) )

    !...............................................
    ! switch timing:
    call GO_Timer_Switch( itim_solve_gain, itim_update_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    !! dump ...
    !write (dumpfile,'("HPHR__",i4.4,2i2.2,"_",2i2.2,"__idom_",i0,".nc")') &
    !          t%year, t%month, t%day, t%hour, t%min, idom
    !call nc_dump( trim(dumpfile), dom%HPHR, 'HPHR', (/'nobs','nobs2'/), status )
    !IF_NOTOK_RETURN(status=1)
    !! dump ...
    !write (dumpfile,'("Kd__",i4.4,2i2.2,"_",2i2.2,"__idom_",i0,"__substate_",a,".nc")') &
    !          t%year, t%month, t%day, t%hour, t%min, idom, trim(substates(substate))
    !call nc_dump( trim(dumpfile), Kd_c, 'Kd', (/'nobs','x','y','z','s'/), status )
    !IF_NOTOK_RETURN(status=1)
 
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
                             Kd_c(iobs,:,:,iz,1) * self%andom%corr(iobs,:,:) * res
              end do ! iz
            end if ! analyse spec?
          !~
          case ( substate_cg )
            ! analyse?
            if ( self%xcorr(iset)%analyse_spec(ispec) ) then
              ! add contribution, apply spatial localization:
              Ens(j)%cg(:,:,ispec) = Ens(j)%cg(:,:,ispec) + &
                           Kd_c(iobs,:,:,1,1) * self%andom%corr(iobs,:,:) * res
            end if ! analyse spec?
          !~
          case ( substate_aerh2o )
            ! analyse sia tracers? then also update aerh2o:
            if ( self%xcorr(iset)%analyse_sia ) then
              ! loop over levels:
              do iz = 1, nz
                ! add contribution, apply spatial localization:
                Ens(j)%aerh2o(:,:,iz) = Ens(j)%aerh2o(:,:,iz) + &
                               Kd_c(iobs,:,:,iz,1) * self%andom%corr(iobs,:,:) * res
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
                                 Kd_c(iobs,:,:,inoise,ihist) * self%andom%corr(iobs,:,:) * res
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

  end subroutine LEKF_Meas_Update_Apply_Domain_SubState


  ! ********************************************************************
  ! ***
  ! *** analysis subdomains
  ! ***
  ! ********************************************************************


  subroutine AnDom_Init( self, status )

    ! --- in/out -------------------------

    class(T_AnDom), intent(out)     ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Init'

    ! --- local --------------------------
    
    ! --- begin ---------------------------

    ! no observations yet:
    self%nobs = 0
    ! safety:
    nullify( self%inds )
    nullify( self%y    )
    nullify( self%r    )
    nullify( self%v    )
    nullify( self%HX   )
    nullify( self%corr )
    nullify( self%dist )
    
    ! ok
    status = 0

  end subroutine AnDom_Init


  ! ***


  subroutine AnDom_Done( self, status )
  
    use Grid, only : Done

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! observation selection allocated?
    if ( associated(self%inds) ) then
      ! clear:
      deallocate( self%inds, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%y   , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%r   , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%v   , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%HX  , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%corr, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%dist, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine AnDom_Done


  ! ***

  
  ! Reset storage for observations to empty.

  subroutine AnDom_Reset( self, status )
  
    use Grid, only : Done

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Reset'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! empty:
    self%nobs = 0

    ! ok
    status = 0

  end subroutine AnDom_Reset


  ! ***


  subroutine AnDom_Add_Observation( self, inds, inside, lon, lat, rho_m, &
                                       y, r, v_loc, Hx_loc, status, debug )
  
    use Grid      , only : RoundToResolution, DistanceGrid
    use LE_Grid   , only : ugg
    !use LEKF_State, only : nmodes_all, nmodes_loc, imodes
    use LEKF_State, only : nmodes

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(in)             ::  inds(:)  ! 3: iset, iloc, obs_ivar
    logical, intent(in)             ::  inside
    real, intent(in)                ::  lon, lat
    real, intent(in)                ::  rho_m ! m
    real, intent(in)                ::  y, r
    real, intent(in)                ::   v_loc(:) ! (nmodes)
    real, intent(in)                ::  Hx_loc(:) ! (nmodes)
    integer, intent(out)            ::  status
    
        logical, intent(in), optional  ::  debug

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Add_Observation'
    
    ! extra number of observations allocated:
    integer, parameter      ::  next = 100

    ! --- local ----------------------------------

    type(T_ObsInd), pointer   ::  new_inds(:)
    real, pointer             ::  new_y(:)
    real, pointer             ::  new_r(:)
    real, pointer             ::  new_v (:,:)
    real, pointer             ::  new_HX(:,:)
    real, pointer             ::  new_corr(:,:,:)

    integer                   ::  j
!    integer                   ::  imode
    real                      ::  lon0, lat0
    
        logical   ::  dbg

    ! --- begin ----------------------------------
    
    ! testing ...
    dbg = .false.
    if ( present(debug) ) dbg = debug
    
    ! new storage needed?
    if ( .not. associated(self%inds) ) then
      ! storage:
      allocate( self%inds(next), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%y(next)   , stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%r(next)   , stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%v (next,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%HX(next,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%corr(next,ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! help array
      allocate( self%dist(ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      ! more storage needed?
      if ( size(self%inds) == self%nobs ) then
        ! new storage:
        allocate( new_inds(self%nobs+next), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_y   (self%nobs+next), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_r   (self%nobs+next), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_v   (self%nobs+next,nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_HX  (self%nobs+next,nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_corr(self%nobs+next,ugg%nlon,ugg%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! only local modes are filled, ensure that the rest is zero:
        new_Hx = 0.0
        ! copy:
        new_inds(1:self%nobs)     = self%inds
        new_y   (1:self%nobs)     = self%y
        new_r   (1:self%nobs)     = self%r
        new_v   (1:self%nobs,:)   = self%v
        new_HX  (1:self%nobs,:)   = self%HX
        new_corr(1:self%nobs,:,:) = self%corr
        ! clear old:
        deallocate( self%inds, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%y   , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%r   , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%v   , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%HX  , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%corr, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! re-assign:
        self%inds => new_inds
        self%y    => new_y
        self%r    => new_r
        self%v    => new_v
        self%HX   => new_HX
        self%corr => new_corr
      end if
    end if
    
    ! increase counter:
    self%nobs = self%nobs + 1
    
    ! store location info:
    self%inds(self%nobs)%iset     = inds(1)
    self%inds(self%nobs)%iloc     = inds(2)
    self%inds(self%nobs)%obs_ivar = inds(3)
    self%inds(self%nobs)%inside   = inside
    self%inds(self%nobs)%lon      = lon
    self%inds(self%nobs)%lat      = lat
    ! observation:
    self%y   (self%nobs)          = y
    self%r   (self%nobs)          = r
    
    ! loop over modes:
    do j = 1, nmodes
      !! global mode index:
      !imode = imodes(j)
      ! copy elements:
      self%v (self%nobs,j) =  v_loc(j)
      self%HX(self%nobs,j) = HX_loc(j)
    end do ! j
    
    ! store length scale:
    if ( self%nobs == 1 ) then
      self%rho_m = rho_m
    else if ( rho_m /= self%rho_m ) then
      write (gol,'("only single localization length scale supported yet;")'); call goErr
      write (gol,'("  current value : ",f12.2," m")') self%rho_m; call goErr
      write (gol,'("  new value     : ",f12.2," m")') rho_m; call goErr
      TRACEBACK; status=1; return
    end if

    ! round to grid cell resolution, if inside domain this is a grid cell center:
    call ugg%RoundToResolution( lon, lat, lon0, lat0, status )!, debug=debug )
    IF_NOTOK_RETURN(status=1)
    ! fill localization field with distances to rounded location:
    call ugg%DistanceGrid( lon0, lat0, self%dist, status )!, debug=debug )
    IF_NOTOK_RETURN(status=1)

!    ! testing ...
!    if ( dbg ) then
!      print *, '  yyy1 added nobs ', self%nobs
!      print *, '    y1 inds ', self%inds(self%nobs)%iset, self%inds(self%nobs)%iloc, self%inds(self%nobs)%obs_ivar
!      print *, '    y1 inside ', self%inds(self%nobs)%inside
!      print *, '    y1 loc  ', self%inds(self%nobs)%lon, self%inds(self%nobs)%lat
!      print *, '    y1 loc0 ', lon0, lat0
!      print *, '    y1 dat  ', self%y(self%nobs), self%r(self%nobs)
!      print *, '    y1 v    ', self%v (self%nobs,:)
!      print *, '    y1 HX   ', self%HX(self%nobs,:)
!      print *, '    y1 dist ', self%corr(self%nobs,:,:)
!    end if

    ! convert to correlation:
    call HCorrGaussian( self%rho_m, self%dist, self%corr(self%nobs,:,:), status )
    IF_NOTOK_RETURN(status=1)

!    ! testing ...
!    if ( dbg ) then
!      print *, '    y1 corr ', self%corr(self%nobs,:,:)
!    end if

    ! ok
    status = 0

  end subroutine AnDom_Add_Observation


  ! ***
  
  
  !
  ! Exchange observations with other domains.
  ! - broadcast bounding box, receive all
  ! - loop over domains
  !   - select obs in current domain that are nearby
  ! - broadcast numbers of nearby obs, receive all
  ! - loop over source domains
  !   - loop over target domains
  !     - send/receive observations
  !     - add to current storage
  !

  subroutine AnDom_Exchange_Nearby( self, status )
  
    use GO        , only : goc
    use C3PO      , only : BBoxWithinDistance
    use LE_Grid   , only : ugg
    use LEKF_State, only : nmodes

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Exchange_Nearby'
    
    ! --- local ----------------------------------

    real                  ::  bbox(4)
    real, allocatable     ::  bboxs(:,:)  ! (4,npes)
    integer               ::  id
    integer               ::  id_send
    integer               ::  id_recv
    integer               ::  nobs0
    integer               ::  iobs
    logical, allocatable  ::  selected(:,:)  ! (nobs,npes)
    integer, allocatable  ::  nsend(:)       ! (npes)      # from this pe to others
    integer, allocatable  ::  nsends(:,:)    ! (npes,npes) # from each pe (2nd dim) to others (1st dim)
    integer               ::  nsnd
    
    ! offsets:
    integer               ::  i0set, i0loc, i0var
    integer               ::  i0lon, i0lat, i0rho, i0y, i0r, i0v, i0hx
    ! buffers:
    integer               ::  ni, nr
    integer, allocatable  ::  ibuffer(:,:)     ! (ni,nsend)
    real, allocatable     ::  rbuffer(:,:)     ! (nr,nsend)
    integer               ::  k
    integer               ::  itag, rtag
    
    ! --- begin ----------------------------------

    ! storage:
    allocate( bboxs(4,0:goc%npes-1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! bounding box of current grid:
    call ugg%GetBoundingBox( bbox, status )
    IF_NOTOK_RETURN(status=1)
    ! send to all, receive all in return:
    call goc%AllGather( bbox, bboxs, status )
    IF_NOTOK_RETURN(status=1)
    
    !! testing ...
    !write (gol,'("xxx bounding boxes:")'); call goPr
    !do id = 0, goc%npes-1
    !  write (gol,'("  x   domain ",i2," bbox [",f8.2,3(",",f8.2),"]")'), id, bboxs(:,id); call goPr
    !end do
    
    ! current number of observations (within this domain),
    ! total number will be increased with nearby stations:
    nobs0 = self%nobs
    
    ! storage:
    allocate( selected(max(1,nobs0),0:goc%npes-1), source=.false., stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( nsend(0:goc%npes-1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( nsends(0:goc%npes-1,0:goc%npes-1), stat=status )
    IF_NOTOK_RETURN(status=1)

    !! testing ..
    !write (gol,'(a," - select nearby observations ...")') rname; call goPr
    ! loop over processors (domains):
    do id = 0, goc%npes-1
    
      ! not this one?
      if ( id /= goc%id ) then
      
        !! testing ..
        !write (gol,'(a," -   domain ",i0)') rname, id; call goPr
        ! loop over observations:
        do iobs = 1, nobs0
        
          ! only if inside (for safety, all are inside before exchange)
          if ( .not. self%inds(iobs)%inside ) cycle
      
          ! select if current observation is nearby the other domain:
          call BboxWithinDistance( bboxs(:,id), self%inds(iobs)%lon, self%inds(iobs)%lat, &
                                     rho_frac_zero*self%rho_m, selected(iobs,id), status )
          IF_NOTOK_RETURN(status=1)
          
          !! testing ..
          !write (gol,'(a," -   observation ",i4," : ",l1)') rname, iobs, selected(iobs,id); call goPr

        end do  ! observations
        
      end if  ! other domain
      
    end do ! domains
    
    ! number of processes to be send per domain:
    nsend = count( selected, dim=1 )
    !! testing ...
    !write (gol,*) 'xx nsend = ', nsend; call goPr
    ! testing ...
    write (gol,'(a," - observations nearby to other domains:")') rname; call goPr
    do id = 0, goc%npes-1
      if ( id == goc%id ) cycle
      write (gol,'(a," -   found ",i4," observations nearby domain ",i2)') rname, nsend(id), id; call goPr
    end do
    
    ! broadcast and receive:
    call goc%AllGather( nsend, nsends, status )
    IF_NOTOK_RETURN(status=1)

    !! testing ...
    !write (gol,'("xxx nsends:")'); call goPr
    !do id = 0, goc%npes-1
    !  write (gol,*) '  x   ', id, ' : ', nsends(id,:); call goPr
    !end do
    
    ! number of integer values, offsets:
    ni = 0
    i0set = ni ; ni = ni + 1
    i0loc = ni ; ni = ni + 1
    i0var = ni ; ni = ni + 1
    ! number of real values, offsets:
    nr = 0
    i0lon = nr ; nr = nr + 1
    i0lat = nr ; nr = nr + 1
    i0rho = nr ; nr = nr + 1
    i0y   = nr ; nr = nr + 1
    i0r   = nr ; nr = nr + 1
    i0v   = nr ; nr = nr + nmodes
    i0hx  = nr ; nr = nr + nmodes

    ! testing ...
    write (gol,'(a," - send/recv observations ...")') rname; call goPr
    ! init tags:
    itag = 0
    rtag = 0
    ! loop over sending domains:
    do id_send = 0, goc%npes-1
      ! loop over receiving domains:
      do id_recv = 0, goc%npes-1

        ! short:
        nsnd = nsends(id_recv,id_send)
        ! exchange needed?
        if ( nsnd > 0 ) then

          !! testing ...
          !write (gol,'("  x  send ",i4," observations from ",i2," to ",i2)') &
          !                 nsnd, id_send, id_recv; call goPr

          ! storage:
          allocate( ibuffer(ni,nsnd), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( rbuffer(nr,nsnd), stat=status )
          IF_NOTOK_RETURN(status=1)

          ! sending pe?
          if ( goc%id == id_send ) then
            !! testing ...
            !write (gol,'("  x    collect observations ...")'); call goPr
            ! collect:
            k = 0
            do iobs = 1, nobs0
              ! selected to be send to receiving pe?
              if ( selected(iobs,id_recv) ) then
                ! increase counter:
                k = k + 1
                ! copy integer values:
                ibuffer(i0set+1,k) = self%inds(iobs)%iset
                ibuffer(i0loc+1,k) = self%inds(iobs)%iloc
                ibuffer(i0var+1,k) = self%inds(iobs)%obs_ivar
                ! copy real values:
                rbuffer(i0lon+1            ,k) = self%inds(iobs)%lon
                rbuffer(i0lat+1            ,k) = self%inds(iobs)%lat
                rbuffer(i0rho+1            ,k) = self%rho_m
                rbuffer(i0y  +1            ,k) = self%y(iobs)
                rbuffer(i0r  +1            ,k) = self%r(iobs)
                rbuffer(i0v  +1:i0v +nmodes,k) = self%v(iobs,:)
                rbuffer(i0hx +1:i0hx+nmodes,k) = self%v(iobs,:)
              end if ! selected
            end do ! observations
          end if ! sending pe

          !! testing ...
          !write (gol,'("  x    send/recv ...")'); call goPr
          ! increase tags:
          itag = itag + 2
          rtag = itag + 1
          ! send and receive:
          call goc%SendAndRecv( ibuffer, id_send, id_recv, itag, status )
          IF_NOTOK_RETURN(status=1)
          ! send and receive:
          call goc%SendAndRecv( rbuffer, id_send, id_recv, rtag, status )
          IF_NOTOK_RETURN(status=1)

          ! receiving?
          if ( goc%id == id_recv ) then
            !! testing ...
            !write (gol,'("  x    unpack ...")'); call goPr
            ! loop over received observations:
            do k = 1, nsnd
              ! add observation, not inside:
              call self%Add_Observation( ibuffer(:,k), .false., &
                                         rbuffer(i0lon+1            ,k), &
                                         rbuffer(i0lat+1            ,k), &
                                         rbuffer(i0rho+1            ,k), &
                                         rbuffer(i0y  +1            ,k), &
                                         rbuffer(i0r  +1            ,k), &
                                         rbuffer(i0v  +1:i0v +nmodes,k), &
                                         rbuffer(i0hx +1:i0hx+nmodes,k), &
                                         status )
              IF_NOTOK_RETURN(status=1)
            end do ! received obs
          end if

          ! clear:
          deallocate( ibuffer, stat=status )
          IF_NOTOK_RETURN(status=1)
          deallocate( rbuffer, stat=status )
          IF_NOTOK_RETURN(status=1)
          
          !! wait ...
          !call goc%Barrier( status )
          !IF_NOTOK_RETURN(status=1)

        end if ! exchange needed

      end do ! recv domain
    end do ! send domain
    
    !! testing ...
    !write (gol,'(a," - end send/recv")') rname; call goPr

    ! clear:
    deallocate( bboxs, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( selected, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( nsend, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( nsends, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    !! testing ..
    !write (gol,'("break")'); call goPr
    !TRACEBACK; status=1; return
    
    ! ok
    status = 0

  end subroutine AnDom_Exchange_Nearby

  ! ***


!  subroutine AnDom_BCast_Modes( self, status )
!  
!    use GO        , only : goc
!    use LEKF_State, only : nmodes_all, imode_pe
!
!    ! --- in/out ---------------------------------
!
!    class(T_AnDom), intent(inout)   ::  self
!    integer, intent(out)            ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/AnDom_BCast_Modes'
!    
!    ! --- local ----------------------------------
!
!    integer                   ::  j
!    integer                   ::  imode
!
!    ! --- begin ----------------------------------
!    
!    ! any observations?
!    if ( self%nobs > 0 ) then
!
!      ! loop over global modes:
!      do imode = 1, nmodes_all
!        ! broadcast from pe holding this mode,
!        ! will be received by all others:
!        call goc%BCast( imode_pe(imode), self%HX(:,imode), status )
!        IF_NOTOK_RETURN(status=1)
!      end do ! modes
!      
!    end if  ! any obs
!
!    ! ok
!    status = 0
!
!  end subroutine AnDom_BCast_Modes


  ! ***


  subroutine AnDom_Fill_HPHR_etc( self, status )
  
    use GO        , only : goc
    use Num       , only : LinAlg_Sym_Factorize
    use Grid      , only : ll_distance
    !use LEKF_State, only : nmodes_all, imode_pe
    use LEKF_State, only : nmodes

      use file_nc

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Fill_HPHR_etc'
    
    ! --- local ----------------------------------

    real, allocatable         ::  Hxd(:)     ! (nobs)
    integer                   ::  imode
    integer                   ::  iobs, iobs2
    real                      ::  d
    real                      ::  c
    character(len=1024)       ::  dumpfile

      !real                      ::  s2max    
      !real, allocatable         ::  A(:,:)

    ! --- begin ----------------------------------
    
    ! output storage:
    allocate( self%HSd(self%nobs,nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%HPHR(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%U(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! temporary storage:
    allocate( Hxd(self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ensemble mean:
    Hxd = sum( self%Hx(1:self%nobs,:), dim=2 )/nmodes
    ! ensemble covariance square root: Santiago EnKF-KA
    do imode = 1, nmodes
      self%HSd(:,imode) = ( self%Hx(1:self%nobs,imode) - Hxd )/sqrt(nmodes-1.0)
    end do

    ! clear:
    deallocate( Hxd, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill HPH:
    self%HPHR = matmul( self%HSd, transpose(self%HSd) )
    
    ! apply localization;
    ! loop over first observations:
    do iobs = 1, self%nobs
      ! loop over remaining observations:
      do iobs2 = iobs, self%nobs
        ! distance in m:
        d = ll_distance( self%inds(iobs )%lon, self%inds(iobs )%lat, &
                         self%inds(iobs2)%lon, self%inds(iobs2)%lat  )
        ! correlation:
        call HCorrGaussian( self%rho_m, d, c, status )
        IF_NOTOK_RETURN(status=1)
        ! apply as factor, use symetry:
        self%HPHR(iobs ,iobs2) = self%HPHR(iobs ,iobs2) * c
        self%HPHR(iobs2,iobs ) = self%HPHR(iobs2,iobs ) * c
      end do ! iobs2
    end do ! iobs
    
      !! testing ...
      !allocate( A(self%nobs,self%nobs), stat=status )
      !IF_NOTOK_RETURN(status=1)
      !! copy:
      !A = self%HPHR
    
    ! add R;
    ! loop over observations:
    do iobs = 1, self%nobs
      ! add variance to diagonal:
      self%HPHR(iobs,iobs) = self%HPHR(iobs,iobs) + self%r(iobs)**2
      !! TESTING: add minimum error ...
      !A(iobs,iobs) = A(iobs,iobs) + max(0.4,self%r(iobs))**2
    end do
    
    ! factorize:
    call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
    if ( status /= 0 ) then
      ! target file:
      write (dumpfile,'("HPHR__pe",i0,".nc")') goc%id
      ! info ...
      write (gol,'("Factorization failed, write matrix to:")'); call goErr
      write (gol,'("  ",a)') trim(dumpfile); call goErr
      ! dump ...
      call nc_dump( trim(dumpfile), self%HPHR, 'HPHR', (/'nobs','nobs2'/), status )
      IF_NOTOK_RETURN(status=1)

        !! info ...
        !write (gol,'("TESTING - increase diagonal ...")'); call goErr
        !! reset:
        !self%HPHR = A
        !! adhoc: increase diagonal to ensure positive definiteness;
        !! this is equivalent to assuming larger obs error ...
        !!~ maximum diagonal value:
        !s2max = 0.0
        !do iobs = 1, self%nobs
        !  s2max = max( s2max, self%HPHR(iobs,iobs) )
        !end do
        !!~ set minimum value of 1% of maximum for diagonal elements:
        !do iobs = 1, self%nobs
        !  self%HPHR(iobs,iobs) = max( 0.01*s2max, self%HPHR(iobs,iobs) )
        !end do
        !! dump:
        !write (dumpfile,'("HPHR__pe",i0,"_v2.nc")') goc%id
        !call nc_dump( trim(dumpfile), self%HPHR, 'HPHR', (/'nobs','nobs2'/), status )
        !IF_NOTOK_RETURN(status=1)
        !! factorize:
        !call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
        !! info ...
        !write (gol,'("TESTING - return status factorization: ",i0)') status; call goErr
    
      ! leave:
      TRACEBACK; status=1; return
    end if
    
!      ! testing ...
!      deallocate( A, stat=status )
!      IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine AnDom_Fill_HPHR_etc


  ! ***


  subroutine AnDom_Clear_HPHR_etc( self, status )
  
    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Clear_HPHR_etc'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! output storage:
    deallocate( self%HSd, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%HPHR, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%U, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AnDom_Clear_HPHR_etc
  
  
  ! ***


  subroutine HCorrGaussian_0d( rho, d, c, status )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)                ::  rho  ! lenght scale (m)
    real, intent(in)                ::  d    ! distance (m)
    real, intent(out)               ::  c    ! correlation [1]
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/HCorrGaussian_0d'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! zero length scale?
    if ( rho < 1.0e-2 ) then
      ! only correlated at zero distance:
      if ( d < 1.0e-2 ) then
        ! full correlation:
        c = 1.0
      else
        ! uncorrelated:
        c = 0.0
      end if
    else
      ! truncate:
      if ( d <= rho_frac_zero*rho ) then
        ! correlation as function of distance:
        c = exp( -0.5*(d/rho)**2 )
      else
        ! uncorrelated:
        c = 0.0
      end if
    end if

    ! ok
    status = 0

  end subroutine HCorrGaussian_0d
  
  
  ! ***


  subroutine HCorrGaussian_2d( rho, d, c, status )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)                ::  rho       ! lenght scale (m)
    real, intent(in)                ::  d(:,:)    ! distance (m)
    real, intent(out)               ::  c(:,:)    ! correlation [1]
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/HCorrGaussian_2d'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! zero length scale?
    if ( rho < 1.0e-2 ) then
      ! only correlated at zero distance:
      where ( d < 1.0e-2 )
        ! full correlation:
        c = 1.0
      elsewhere
        ! uncorrelated:
        c = 0.0
      end where
    else
      ! truncate:
      where ( d <= rho_frac_zero*rho )
        ! correlation as function of distance:
        c = exp( -0.5*(d/rho)**2 )
      elsewhere
        ! uncorrelated:
        c = 0.0
      end where
    end if

    ! ok
    status = 0

  end subroutine HCorrGaussian_2d


end module LEKF_Meas_Update

