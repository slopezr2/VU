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

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_StrError, NF90_NOERR

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
    ! maori obs identified by set, location, and observed variable:
    integer                       ::  iset
    integer                       ::  iloc
    integer                       ::  obs_ivar
    ! inside domain?
    logical                       ::  inside
    ! representative location:
    real                          ::  lon, lat
    real                          ::  lon0, lat0   ! rounded to grid cell center ...
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
    real, pointer                 ::  y(:)      ! (nobs)
    real, pointer                 ::  r(:)      ! (nobs)
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
    procedure ::  Dump                  =>  AnDom_Dump
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
      call self%Apply_Domain( t2, status )
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
#ifdef with_kf_meas_sat
    use LE_Output        , only : T_LE_Output_Sat_Data
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
    
#ifdef with_kf_meas_sat
!    integer                               ::  displ
    real, allocatable                     ::  chi(:,:)    ! (npix,nmodes)
    integer                               ::  k
    type(T_LEKF_Meas_Sat_Set), pointer    ::  msat
    type(T_LE_Output_Sat_Data), pointer   ::  satd
    real                                  ::  bbox(4)
    integer                               ::  glbid
#endif

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

#ifdef with_kf_meas_sat
    ! info ..
    write (gol,'(a," - collect sat observations ...")') rname; call goPr
    ! loop over sat sets:
    do k = 1, meas_sat%nset
      ! set index:
      iset = iset + 1
      ! pointers:
      msat => meas_sat%set(k)
      satd => leo%sat_data(meas_sat%iout(k))
      ! info ..
      write (gol,'(a," -   sat ",i0," `",a,"`, iset=",i0," ...")') rname, k, trim(msat%name), iset; call goPr
      ! analyse?
      if ( msat%analyse ) then

        ! info ...
        write (gol,'(a," -   to be analyzed ...")') rname; call goPr

        ! total number over all processors:
        call goc%ParInfo( satd%npix, status, ntot=nmeas_tot )!, displ=displ )
        IF_NOTOK_RETURN(status=1)

        ! any measurements?
        if ( nmeas_tot > 0 ) then

          ! info ...
          write (gol,'(a," -   number of pixels : ",i0," (local ",i0,")")') rname, nmeas_tot, satd%npix; call goPr
          
          ! generate random N(0,1) numbers per pixel and mode,
          ! get them now to have later v = chi * sigma
          ! storage:
          allocate( chi(max(1,satd%npix),nmodes), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! loop over modes:
          do j = 1, nmodes
            ! get random nunmbers:
            call satd%GetRandom( msat%rnd, chi(:,j), status )
            IF_NOTOK_RETURN(status=1)
          end do ! modes

          ! convert from km to m:
          rho_m = msat%rho * 1.0e3  ! m

          ! screening factor:
          alfa = msat%screening_factor

          ! loop over measurements assigned to this domain:
          do ipix = 1, satd%npix

            ! check if not accepted (validation, no data?):
            if ( satd%astatus%data(ipix) /= 0 ) then
               write (gol,'("sat observation status invalid: ",i0)') satd%astatus%data(ipix); call goErr
               TRACEBACK; status=1; return         
            end if

            ! extract measurement data:
            call satd%GetPixel( ipix, status, glbid=glbid, lon=lon, lat=lat, y=y, sigma=r )
            IF_NOTOK_RETURN(status=1)

            ! define how to extract value from state:
            inds(0:2) = (/ substate_sat, k, ipix /)

            ! loop over modes:
            do j = 1, nmodes
              ! extract simulation of measurements:
              call Ens(j)%GetValue( inds, HX_loc(j), status )
              IF_NOTOK_RETURN(status=1)
              ! random error out N(0,r):
              v_loc(j) = chi(ipix,j) * r
            end do
            
            !! testing ..
            !write (gol,*) 'vvv1 v = ', ipix, glbid, v_loc; call goPr

            ! check if in domain (not in surrounding area):
            call ugg%InDomain( lon, lat, inside, status )
            IF_NOTOK_RETURN(status=1)
            ! expected to be inside ..
            if ( .not. inside ) then
              call ugg%GetBoundingBox( bbox, status )
              IF_NOTOK_RETURN(status=1)
              write (gol,'("fould location outside domain:")'); call goErr
              write (gol,'("  iset    : ",i0)') iset; call goErr
              write (gol,'("  ipix    : ",i0)') ipix; call goErr
              write (gol,'("  lon,lat : ",2f8.2)') lon,lat; call goErr
              write (gol,'("  bbox    : ",4f8.2)') bbox; call goErr
              TRACEBACK; status=1; return
            end if

            ! add to domain, identify by set number (used for chem/noise correlations)
            ! and pixel number:
            call self%andom%Add_Observation( (/iset,ipix,glbid/), inside, &
                                                     lon, lat, rho_m, &
                                                     y, r, v_loc, HX_loc, &
                                                     status, debug=.true. )
            IF_NOTOK_RETURN(status=1)

          end do ! pixels

          ! clear:
          deallocate( chi, stat=status )
          IF_NOTOK_RETURN(status=1)

        else

          ! info ...
          write (gol,'(a," -   no pixels for this timestep ...")') rname; call goPr

        end if  ! nmeas_tot > 0

      else

        ! info ...
        write (gol,'(a," -   not to be analyzed ...")') rname; call goPr

      end if  ! analyse?
      
    end do ! sat sets
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

  subroutine LEKF_Meas_Update_Apply_Domain( self, t, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use Dims      , only : nspec
    use LEKF_State, only : substate_c, substate_cg, substate_aerh2o, substate_dc
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    type(TDate), intent(in)               ::  t
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
    
    !! testing ...
    !call self%andom%Dump( status, t=t )
    !IF_NOTOK_RETURN(status=1)
    
    ! call for specific elements;
    ! update concentration arrays per spec to limit memory usage:
    do ispec = 1, nspec
      call self%Apply_Domain_SubState( substate_c     , ispec, t, status )
      IF_NOTOK_RETURN(status=1)
      call self%Apply_Domain_SubState( substate_cg    , ispec, t, status )
      IF_NOTOK_RETURN(status=1)
    end do
    ! update full arrays, pass dummy ispec:
    call self%Apply_Domain_SubState( substate_aerh2o, -999, t, status )
    IF_NOTOK_RETURN(status=1)
    call self%Apply_Domain_SubState( substate_dc    , -999, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear current arrays
    call self%andom%Clear_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Apply_Domain
  
  
  ! *


  subroutine LEKF_Meas_Update_Apply_Domain_SubState( self, substate, ispec, t, status )
  
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
    ! store:
    self%inds(self%nobs)%lon0 = lon0
    self%inds(self%nobs)%lat0 = lat0
    
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
    real,allocatable          ::  inflation(:)
    !character(len=1024)       ::  dumpfile

      real                      ::  s2max,factor    
      real, allocatable         ::  A(:,:)

    ! --- begin ----------------------------------
    
    ! output storage:
    allocate( self%HSd(self%nobs,nmodes), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%HPHR(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%U(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate(inflation(self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    call RANDOM_NUMBER(inflation)
    ! temporary storage:
    allocate( Hxd(self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ensemble mean:
    Hxd = sum( self%Hx(1:self%nobs,:), dim=2 )/nmodes
    ! ensemble covariance square root Santiago EnKF-KA and inflation:
    do imode = 1, nmodes
      self%HSd(:,imode) = ( self%Hx(1:self%nobs,imode)+0*inflation*Hxd - Hxd )/sqrt(nmodes-1.0)
    end do

    ! clear:
    deallocate( Hxd, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( inflation, stat=status )
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
      allocate( A(self%nobs,self%nobs), stat=status )
      IF_NOTOK_RETURN(status=1)
      !! copy:
    ! add R;
    ! loop over observations:
    do iobs = 1, self%nobs
      ! add variance to diagonal:
      self%HPHR(iobs,iobs) = self%HPHR(iobs,iobs) + self%r(iobs)**2
      !! TESTING: add minimum error ...
      !A(iobs,iobs) = A(iobs,iobs) + max(0.4,self%r(iobs))**2
    end do
    A = self%HPHR
    ! factorize:
    call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
    if ( status /= 0 ) then
      write (gol,'("Factorization failed, dump matrices...")'); call goErr
      ! dump HPHR etc:
      call self%Dump( status )
      IF_NOTOK_RETURN(status=1)
      factor=0.01
      status=1
      do while (status/=0)
        write(*,*) 'Santiago Status', status
        write(*,*) 'Santiago aplying factor of ', factor 
   			self%HPHR = A
        s2max = 0.0
        do iobs = 1, self%nobs
          s2max = max( s2max, self%HPHR(iobs,iobs) )
        end do
        !!~ set minimum value of 1% of maximum for diagonal elements:
        do iobs = 1, self%nobs
          self%HPHR(iobs,iobs) = max( factor*s2max, self%HPHR(iobs,iobs) )
        end do
        call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
        factor=factor+0.01
			end do  
      ! leave:
!        if ( status /= 0 ) then
!        	self%HPHR = A
!          s2max = 0.0
!       	  do iobs = 1, self%nobs
!          	s2max = max( s2max, self%HPHR(iobs,iobs) )
!          end do
!          do iobs = 1, self%nobs
!          	self%HPHR(iobs,iobs) = max( 0.5*s2max, self%HPHR(iobs,iobs) )
!          end do
!          call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
!          if ( status /= 0 ) then
!         		self%HPHR = A
!          	s2max = 0.0
!       	  	do iobs = 1, self%nobs
!          		s2max = max( s2max, self%HPHR(iobs,iobs) )
!          	end do
!          	do iobs = 1, self%nobs
!          		self%HPHR(iobs,iobs) = max( 1*s2max, self%HPHR(iobs,iobs) )
!          	end do
!          	call LinAlg_Sym_Factorize( self%HPHR, self%U, status )	
!         		if ( status /= 0 ) then
!         			self%HPHR = A
!		        	s2max = 0.0
!		     	  	do iobs = 1, self%nobs
!		        		s2max = max( s2max, self%HPHR(iobs,iobs) )
!		        	end do
!		        	do iobs = 1, self%nobs
!		        		self%HPHR(iobs,iobs) = max( 1.5*s2max, self%HPHR(iobs,iobs) )
!		        	end do
!		        	call LinAlg_Sym_Factorize( self%HPHR, self%U, status )	
!		       		if ( status /= 0 ) then
!		       			self%HPHR = A
!				      	s2max = 0.0
!				   	  	do iobs = 1, self%nobs
!				      		s2max = max( s2max, self%HPHR(iobs,iobs) )
!				      	end do
!				      	do iobs = 1, self%nobs
!				      		self%HPHR(iobs,iobs) = max( 2.5*s2max, self%HPHR(iobs,iobs) )
!				      	end do
!				      	call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
!		       		
!		       		  if ( status /= 0 ) then
!		       		  	self%HPHR = A
!				      		s2max = 0.0
!				   	  		do iobs = 1, self%nobs
!				      			s2max = max( s2max, self%HPHR(iobs,iobs) )
!				      		end do
!				      		do iobs = 1, self%nobs
!				      			self%HPHR(iobs,iobs) = max( 10*s2max, self%HPHR(iobs,iobs) )
!				      		end do
!				      		call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
!		       		  	if ( status /= 0 ) then
!         		 				TRACEBACK; status=1; return
!         		 			endif
!         		 		endif	
!         		 	endif
!         		endif
!          endif		
!      	end if
    end if
    
!      ! testing ...
		deallocate( A, stat=status )
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


  subroutine AnDom_Dump( self, status, t )

    use NetCDF, only : NF90_Create, NF90_Close
    use NetCDF, only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF, only : NF90_Def_Dim
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_INT, NF90_FLOAT
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_EndDef
    use NetCDF, only : NF90_Put_Var
  
    use GO        , only : TDate
    use GO        , only : goc
    use LE_Config , only : outputdir
    use LEKF_State, only : nmodes
    use LE_Grid   , only : ugg

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)       ::  self
    integer, intent(out)                ::  status
    type(TDate), intent(in), optional   ::  t

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Dump'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  filename
    integer                   ::  cmode
    integer                   ::  ncid
    integer                   ::  dimid_obs
    integer                   ::  dimid_mode
    integer                   ::  dimid_lon, dimid_lat
    integer                   ::  varid_iset, varid_iloc, varid_ivar
    integer                   ::  varid_inside
    integer, allocatable      ::  inside(:)
    integer                   ::  varid_lon, varid_lat
    integer                   ::  varid_lon0, varid_lat0
    integer                   ::  varid_dist
    integer                   ::  varid_corr
    integer                   ::  varid_y, varid_r, varid_v, varid_HX
    integer                   ::  varid_HSd, varid_HPHR

    ! --- begin ----------------------------------
    
    ! any local obs?
    if ( self%nobs > 0 ) then

      ! target file:
      if ( present(t) ) then
        write (filename,'(a,"/AnDom__",i4.4,2i2.2,"_",2i2.2,"__dom_",i2.2,".nc")') &
                    trim(outputdir), t%year, t%month, t%day, t%hour, t%min, goc%id
      else
        write (filename,'(a,"/AnDom__dom_",i2.2,".nc")') trim(outputdir), goc%id
      end if

      ! set creation mode flag:
      cmode = NF90_CLOBBER       ! overwrite existing files

      ! create file:
      status = NF90_Create( filename, cmode, ncid )
      if ( status /= 0 ) then
         write (gol,'("creating file :")'); call goErr
         write (gol,'("  ",a)') trim(filename); call goErr
         TRACEBACK; status=1; return
      end if

      ! define dimensions:
      status = NF90_Def_Dim( ncid, 'obs', self%nobs, dimid_obs )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( ncid, 'mode', nmodes, dimid_mode )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( ncid, 'lon', ugg%nlon, dimid_lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( ncid, 'lat', ugg%nlat, dimid_lat )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variable:
      status = NF90_Def_Var( ncid, 'iset'     , NF90_INT  , (/dimid_obs/), varid_iset )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'iloc'     , NF90_INT  , (/dimid_obs/), varid_iloc )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'ivar'     , NF90_INT  , (/dimid_obs/), varid_ivar )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'inside'   , NF90_INT  , (/dimid_obs/), varid_inside )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'longitude', NF90_FLOAT, (/dimid_obs/), varid_lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'latitude' , NF90_FLOAT, (/dimid_obs/), varid_lat )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'lon0'     , NF90_FLOAT, (/dimid_obs/), varid_lon0 )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'lat0'     , NF90_FLOAT, (/dimid_obs/), varid_lat0 )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variable:
      status = NF90_Def_Var( ncid, 'y'     , NF90_FLOAT, (/dimid_obs/), varid_y )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'r'     , NF90_FLOAT, (/dimid_obs/), varid_r )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'v'     , NF90_FLOAT, (/dimid_obs,dimid_mode/), varid_v )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'HX'    , NF90_FLOAT, (/dimid_obs,dimid_mode/), varid_HX )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variable:
      status = NF90_Def_Var( ncid, 'dist'     , NF90_FLOAT, (/dimid_lon,dimid_lat/), varid_dist )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'corr'     , NF90_FLOAT, (/dimid_obs,dimid_lon,dimid_lat/), varid_corr )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variable:
      status = NF90_Def_Var( ncid, 'HSd'     , NF90_FLOAT, (/dimid_obs,dimid_mode/), varid_HSd )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( ncid, 'HPHR'    , NF90_FLOAT, (/dimid_obs,dimid_obs/), varid_HPHR )
      IF_NF90_NOTOK_RETURN(status=1)

      ! end defintion mode:
      status = NF90_EndDef( ncid )
      IF_NF90_NOTOK_RETURN(status=1)

      ! write variable:
      status = NF90_Put_Var( ncid, varid_iset, self%inds(1:self%nobs)%iset )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_iloc, self%inds(1:self%nobs)%iloc )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_ivar, self%inds(1:self%nobs)%obs_ivar )
      IF_NF90_NOTOK_RETURN(status=1)

      ! storage:
      allocate( inside(self%nobs), stat=status )
      IF_NF90_NOTOK_RETURN(status=1)
      ! fill:
      where ( self%inds(1:self%nobs)%inside )
        inside = 1
      elsewhere
        inside = 0
      end where
      ! write variable:
      status = NF90_Put_Var( ncid, varid_inside, inside )
      IF_NF90_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( inside, stat=status )
      IF_NF90_NOTOK_RETURN(status=1)

      ! write variable:
      status = NF90_Put_Var( ncid, varid_lon, self%inds(1:self%nobs)%lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_lat, self%inds(1:self%nobs)%lat )
      IF_NF90_NOTOK_RETURN(status=1)

      ! write variable:
      status = NF90_Put_Var( ncid, varid_lon0, self%inds(1:self%nobs)%lon0 )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_lat0, self%inds(1:self%nobs)%lat0 )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! write variable:
      status = NF90_Put_Var( ncid, varid_y, self%y(1:self%nobs) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_r, self%r(1:self%nobs) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_v, self%v(1:self%nobs,:) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_HX, self%HX(1:self%nobs,:) )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! write variable:
      status = NF90_Put_Var( ncid, varid_dist, self%dist )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_corr, self%corr(1:self%nobs,:,:) )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! write variable:
      status = NF90_Put_Var( ncid, varid_HSd, self%HSd(1:self%nobs,:) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( ncid, varid_HPHR, self%HPHR(1:self%nobs,1:self%nobs) )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! close:
      status = NF90_Close( ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
    end if

    ! ok
    status = 0

  end subroutine AnDom_Dump
  
  
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

    ! zero length scale? check on 100 m:
    if ( rho < 100.0 ) then
      ! only correlated at zero distance:
      if ( d < rho_frac_zero*100.0 ) then
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

    ! zero length scale? check on 100 m:
    if ( rho < 100.0 ) then
      ! only correlated at zero distance:
      where ( d < rho_frac_zero*100.0 )
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

