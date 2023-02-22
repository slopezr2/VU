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
    ! analysis domain:
    type(TllGridInfo)             ::  lli
    ! offset in global domain:
    integer                       ::  i0, j0
    ! selected nearby observations:
    integer                       ::  nobs
    ! reference to observation data:
    type(T_ObsInd), pointer       ::  inds(:)  ! (>= n)
    ! observation and error stdv:
    real, pointer                 ::  y(:)
    real, pointer                 ::  r(:)
    real, pointer                 ::  v(:,:)    ! (nobs,nmodes_all)
    real, pointer                 ::  HX(:,:)   ! (nobs,nmodes_all)
    ! correlation length scale for locaization,
    ! only one value supported yet:
    real                          ::  rho_m   ! m
    ! spatial localization weight:
    real, pointer                 ::  corr(:,:,:)  ! (nobs,nlon,nlat)
    ! help array with distances:
    real, pointer                 ::  dist(:,:)  ! (nlon,nlat)
    ! analysis matrices:
    real, allocatable             ::  HSd(:,:)   ! (nobs,nmodes_all)
    real, allocatable             ::  HPHR(:,:)  ! (nobs,nobs)
    real, allocatable             ::  U(:,:)     ! (nobs,nobs)
  contains
    procedure ::  Init                  =>  AnDom_Init
    procedure ::  Done                  =>  AnDom_Done
    procedure ::  Reset                 =>  AnDom_Reset
    procedure ::  Add_Observation       =>  AnDom_Add_Observation
    procedure ::  BCast_Modes           =>  AnDom_BCast_Modes
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
    ! analysis domains:
    type(T_AnDom), pointer                 ::  dom(:)  ! (ndom)
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
  integer                                ::  itim_obs_bcast
  integer                                ::  itim_domains
  integer                                ::  itim_exchange_modes
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
    call GO_Timer_Def( itim_obs_collect   , 'lekf update obs collect'   , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_obs_bcast     , 'lekf update obs bcast'     , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_domains       , 'lekf update domains'       , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_exchange_modes, 'lekf update exchange modes', status )
    IF_NOTOK_RETURN(status=1)
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

    use GO             , only : TrcFile, ReadRc
    use Grid           , only : Init
    use MAORI          , only : MAORI_Data_Inquire, MAORI_Data_Set_Inquire
    use LE_Grid        , only : lli
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
    
    integer       ::  ndom_x, ndom_y
    integer       ::  idom_x, idom_y, idom
    integer       ::  nx_loc, ny_loc
    integer       ::  i0, j0
    
    integer               ::  nset_maori
    integer               ::  nset_omi
    integer               ::  nset, iset
    character(len=64)     ::  sname
    logical               ::  analyse
    character(len=1024)   ::  key

    ! --- begin ---------------------------

    ! info ...
    write (gol,'(a," - init analysis domains ...")') rname; call goPr

    ! read distribution of sub-domains for analysis:
    call ReadRc( rcF, 'lekf.meas.andom.nx', ndom_x, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'lekf.meas.andom.ny', ndom_y, status )
    IF_NOTOK_RETURN(status=1)
    ! total number:
    self%ndom = ndom_x * ndom_y
      
    ! allocate analysis domains:
    allocate( self%dom(self%ndom), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! init index:
    idom = 0
    ! loop over domain grid:
    do idom_y = 1, ndom_y
      do idom_x = 1, ndom_x
        ! increase counter:
        idom = idom + 1
        ! cells:
        nx_loc = lli%nlon / ndom_x
        ny_loc = lli%nlat / ndom_y
        ! info ...
        write (gol,'(a," - an domain ",i0," (",i0,",",i0,") shape (",i0,",",i0,")")') &
                       rname, idom, idom_x, idom_y, nx_loc, ny_loc; call goPr
        ! check ...
        if ( modulo(lli%nlon,nx_loc) /= 0 ) then
          write (gol,'(a," - number of x analysis domain ",i0," should match with grid size ",i0)') &
                               rname, ndom_x, lli%nlon; call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( modulo(lli%nlat,ny_loc) /= 0 ) then
          write (gol,'(a," - number of y analysis domain ",i0," should match with grid size ",i0)') &
                               rname, ndom_y, lli%nlat; call goErr
          TRACEBACK; status=1; return
        end if
        ! ofset:
        i0 = (idom_x-1)*nx_loc
        j0 = (idom_y-1)*ny_loc
        ! init domain:
        call self%dom(idom)%Init( lli%lon_deg(1)+i0*lli%dlon_deg, lli%dlon_deg, nx_loc, i0, &
                                     lli%lat_deg(1)+j0*lli%dlat_deg, lli%dlat_deg, ny_loc, j0, status )
        IF_NOTOK_RETURN(status=1)
        ! info ..
        write (gol,'(a," -   bounding box : [",f8.2,",",f8.2,"] x [",f8.2,",",f8.2,"]")') rname, &
                        self%dom(idom)%lli%blon_deg(0), self%dom(idom)%lli%blon_deg(nx_loc), &
                        self%dom(idom)%lli%blat_deg(0), self%dom(idom)%lli%blat_deg(ny_loc); call goPr
      end do ! ix
    end do ! iy

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
        call ReadRc( rcF, 'maori.'//trim(sname)//'.assim.spec', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_spec(nspec_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetSpecApply( key, self%xcorr(iset)%analyse_spec, self%xcorr(iset)%analyse_sia, status )
        IF_NOTOK_RETURN(status=1)

        ! read noise correlation keys:
        call ReadRc( rcF, 'maori.'//trim(sname)//'.assim.noise', key, status )
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
        call ReadRc( rcF, 'kf.meas.omi_trc.spec', key, status )
        IF_NOTOK_RETURN(status=1)
        ! storage for noise correlation weights:
        allocate( self%xcorr(iset)%analyse_spec(nspec_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        call GetSpecApply( key, self%xcorr(iset)%analyse_spec, self%xcorr(iset)%analyse_sia, status )
        IF_NOTOK_RETURN(status=1)

        ! read noise correlation keys:
        call ReadRc( rcF, 'kf.meas.omi_trc.noise', key, status )
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

    ! loop:
    do idom = 1, self%ndom
      ! done:
      call self%dom(idom)%Done( status )
      IF_NOTOK_RETURN(status=1)
    end do
    ! clear:
    deallocate( self%dom, stat=status )
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
    
    integer     ::  idom

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

    ! info ...
    write (gol,'(a," -   broadcast ensemble simulations ...")') rname; call goPr

    ! start timing:
    call GO_Timer_Start( itim_obs_bcast, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ensure that all nodes have the current simulations for all modes,
    ! not influenced by upcoming measurement updates ...
    ! loop over analysis domains:
    do idom = 1, self%ndom
      ! broadcast simulations for local ensemble members to other nodes:
      call self%dom(idom)%BCast_Modes( status )
      IF_NOTOK_RETURN(status=1)
    end do ! domains

    ! switch timing:
    call GO_Timer_Switch( itim_obs_bcast, itim_domains, status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ...
    write (gol,'(a," -   loop over domains ...")') rname; call goPr
    ! measurement update per domain:
    do idom = 1, self%ndom
      ! any obs?
      if ( self%dom(idom)%nobs > 0 ) then
        ! info ...
        write (gol,'(a," -     analyse domain ",i0," ...")')  rname, idom; call goPr
        ! apply:
        call self%Apply_Domain( idom, status )
        IF_NOTOK_RETURN(status=1)
      else
        ! info ...
        write (gol,'(a," -     no observations in domain ",i0," ...")')  rname, idom; call goPr
      end if
    end do

    ! switch timing:
    call GO_Timer_End( itim_domains, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Apply
  
  
  ! ***


  subroutine LEKF_Meas_Update_CollectObservations( self, t1, t2, status )

    use GO        , only : TDate, NewDate, operator(<), operator(<=), wrtgol
    use Grid      , only : WithinDistance
    use Grid      , only : In_Domain
    use MAORI     , only : MAORI_Data_Inquire
    use MAORI     , only : MAORI_Data_Set_Inquire
    use MAORI     , only : MAORI_Data_Loc_Inquire
    use MAORI     , only : MAORI_Data_Obs_Get
    use MAORI     , only : MAORI_Data_Obs_Put
    use MAORI     , only : MAORI_State_Obs_Get
    use MAORI     , only : MAORI_SAMPLE, MAORI_TYPE_NAME
    use MAORI     , only : MAORI_ASTAT_NODATA
    use LEKF_State, only : nmodes_loc
    use LEKF_State, only : substate_omi_trc
    use LEKF_State, only : Ens
#ifdef with_kf_meas_maori
    use LEKF_Data , only : mad
#endif
#ifdef with_kf_meas_omi_trc
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
    real                        ::  lon, lat
    integer                     ::  obs_nvar, obs_ivar
    real                        ::  y, r
    real                        ::  alfa
    integer                     ::  idom
    logical                     ::  nearby
    logical                     ::  inside
    integer                     ::  j
    real, allocatable           ::   v_loc(:)  ! (nmodes_loc)
    real, allocatable           ::  HX_loc(:)  ! (nmodes_loc)

    integer                     ::  ipix
    type(TDate)                 ::  t
    integer                     ::  inds(0:7)
    
    !! testing ...
    !real, allocatable   ::  vv(:,:)  ! (nmodes_loc,npix)
    !integer             ::  ipix1, ipix2, ipixstep
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a," - start")') rname; call goPr

    ! loop over analysis domains:
    do idom = 1, self%ndom
      ! reset to empty:
      call self%dom(idom)%Reset( status )
      IF_NOTOK_RETURN(status=1)
    end do ! idom
    
    ! storage for ensemble simulations, 
    ! dummy in case no modes are allocated here (master pe with xb only):
    allocate(  v_loc(max(1,nmodes_loc)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( HX_loc(max(1,nmodes_loc)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! *

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

          ! get number of locations:
          call MAORI_Data_Set_Inquire( mad, iset, status, &
                                        nloc=nloc, obs_nvar=obs_nvar )
          IF_NOTOK_RETURN(status=1)
    
          ! info ..
          write (gol,'(a," -     number of locations: ",i0)') rname, nloc; call goPr

          ! loop over locations:
          do iloc = 1, nloc

            ! get location:
            call MAORI_Data_Loc_Inquire( mad, iset, iloc, status, &
                                             lon=lon, lat=lat )
            IF_NOTOK_RETURN(status=1)

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
              
              ! loop over local modes:
              do j = 1, nmodes_loc
                ! extract simulation of measurements:
                call MAORI_State_Obs_Get( Ens(j)%mas, mad, iset, obs_ivar, status, &
                                            iloc=iloc, value=HX_loc(j) )
                IF_NOTOK_RETURN(status=1)
                ! random error out N(0,r)
                call Ens(j)%rnd%Get_Normal( v_loc(j), status, sigma=r )
                IF_NOTOK_RETURN(status=1)
              end do

              ! loop over analysis domains:
              do idom = 1, self%ndom
              
                ! within factor*rho [m] of domain?
                call WithinDistance( self%dom(idom)%lli, lon, lat, rho_frac_zero * rho_m, nearby, status )
                IF_NOTOK_RETURN(status=1)

                ! nearby observation? then add to list:
                if ( nearby ) then

                  ! check if in domain (not in surrounding area):
                  call In_Domain( self%dom(idom)%lli, lon, lat, inside, status )
                  IF_NOTOK_RETURN(status=1)

                  ! add to domain:
                  call self%dom(idom)%Add_Observation( (/iset,iloc,obs_ivar/), inside, &
                                                           lon, lat, rho_m, &
                                                           y, r, v_loc, HX_loc, &
                                                           status )
                  IF_NOTOK_RETURN(status=1)

                  !! info ...
                  !write (gol,'(a," -   added obs (set ",i0,", loc ",i0,", var ",i0,") at (",f8.2,",",f8.2,") to domain ",i0)') &
                  !                        rname, iset, iloc, obs_ivar, lon, lat, idom; call goPr

                end if ! nearby domain

              end do  ! analysis domains

            end do  ! observed variables

          end do  ! locations

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("ERROR - unsupported maori type ",i6," (`",a,"`) for set ",i6)') &
                                   stype, trim(MAORI_TYPE_NAME(stype)), iset; call goPr
          TRACEBACK; status=1; return

      end select

    end do ! maori sets
#endif

    ! *
    
#ifdef with_kf_meas_omi_trc
    ! info ..
    write (gol,'(a," - collect OMI observations ...")') rname; call goPr
    ! set index:
    iset = nset_maori + 1
    ! to be analysed?
    if ( omi_trc%analyse ) then
    
      ! info ...
      write (gol,'(a," -   OMI observations to be analyzed ...")') rname; call goPr
 
      ! any measurements?
      if ( (omi_trc%nmeas > 0) .and. (omi_trc%nmeas_intimestep > 0) ) then
    
        ! info ...
        write (gol,'(a," -   number of pixels for this timestep : ",i0)') rname, omi_trc%nmeas; call goPr
     
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

        ! loop over measurements:
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

          ! extract measurement:
          lon  = omi_trc%longitude(ipix)
          lat  = omi_trc%latitude(ipix)
          y    = omi_trc%vcd_trop(ipix)
          r    = omi_trc%sigma_vcd_trop(ipix)
          alfa = omi_trc%screening_factor

          ! define how to extract value from state:
          inds(0:1) = (/ substate_omi_trc, ipix /)

          ! loop over local modes:
          do j = 1, nmodes_loc
            ! extract simulation of measurements:
            call Ens(j)%GetValue( inds, HX_loc(j), status )
            IF_NOTOK_RETURN(status=1)
            ! random error out N(0,r)
            call Ens(j)%rnd%Get_Normal( v_loc(j), status, sigma=r )
            IF_NOTOK_RETURN(status=1)
            !! testing reversed order:
            !v_loc(j) = vv(j,ipix)
          end do

          ! loop over analysis domains:
          do idom = 1, self%ndom

            ! within factor*rho [m] of domain?
            call WithinDistance( self%dom(idom)%lli, lon, lat, rho_frac_zero * rho_m, nearby, status )
            IF_NOTOK_RETURN(status=1)

            ! nearby observation? then add to list:
            if ( nearby ) then

              ! check if in domain (not in surrounding area):
              call In_Domain( self%dom(idom)%lli, lon, lat, inside, status )
              IF_NOTOK_RETURN(status=1)

              ! add to domain, identify by set number (used for chem/noise correlations)
              ! and pixel number:
              call self%dom(idom)%Add_Observation( (/iset,ipix,-999/), inside, &
                                                       lon, lat, rho_m, &
                                                       y, r, v_loc, HX_loc, &
                                                       status, debug=(idom==1) )
              IF_NOTOK_RETURN(status=1)

              !! info ...
              !write (gol,'(a," -   added OMI pixel (index ",i0,") at (",f8.2,",",f8.2,") to domain ",i0)') &
              !                        rname, ipix, lon, lat, idom; call goPr

            end if ! nearby domain

          end do  ! analysis domains

        end do ! pixels
        
        !! clear:
        !deallocate( vv, stat=status )
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

    ! *
        
    ! clear:
    deallocate( HX_loc, stat=status )
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

  subroutine LEKF_Meas_Update_Apply_Domain( self, idom, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use Dims      , only : nspec
    use LEKF_State, only : substate_c, substate_cg, substate_aerh2o, substate_dc
  
    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    integer, intent(in)                   ::  idom
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply_Domain'

    ! --- local ----------------------------------
    
    type(T_AnDom), pointer    ::  dom
    integer                   ::  ispec

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a," - measurement update for domain ",i0)') rname, idom; call goPr
    
    ! short:
    dom => self%dom(idom)
    
    ! fill observation simulation entities:
    !    dom%HSd
    !    dom%HPHR
    !    dom%U
    ! allocated on return:
    call dom%Fill_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! call for specific elements;
    ! update concentration arrays per spec to limit memory usage:
    do ispec = 1, nspec
      call self%Apply_Domain_SubState( idom, substate_c     , ispec, status )
      IF_NOTOK_RETURN(status=1)
      call self%Apply_Domain_SubState( idom, substate_cg    , ispec, status )
      IF_NOTOK_RETURN(status=1)
    end do
    ! update full arrays, pass dummy ispec:
    call self%Apply_Domain_SubState( idom, substate_aerh2o, -999, status )
    IF_NOTOK_RETURN(status=1)
    call self%Apply_Domain_SubState( idom, substate_dc    , -999, status )
    IF_NOTOK_RETURN(status=1)

    ! clear current arrays
    call dom%Clear_HPHR_etc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LEKF_Meas_Update_Apply_Domain
  
  
  ! *


  subroutine LEKF_Meas_Update_Apply_Domain_SubState( self, idom, substate, ispec, status )
  
    use GO        , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use GO        , only : TDate
    use GO        , only : goc
    use Num       , only : LinAlg_Sym_FactorSolve
    use Dims      , only : nx, ny, nz, nspec
    use LEKF_noise, only : nnoise, nhist
    use LEKF_State, only : nmodes_all, nmodes_loc
    use LEKF_State, only : imodes, imode_pe
    use LEKF_State, only : Ens
    use LEKF_State, only : substates, substate_c, substate_cg, substate_aerh2o, substate_dc
  
      !use file_nc

    ! --- in/out ---------------------------------

    class(T_Meas_Update), intent(inout)   ::  self
    integer, intent(in)                   ::  idom
    integer, intent(in)                   ::  substate
    integer, intent(in)                   ::  ispec
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Meas_Update_Apply_Domain_SubState'

    ! --- local ----------------------------------
    
    type(T_AnDom), pointer    ::  dom
    integer                   ::  nxd, nyd
    integer                   ::  i1, i2, j1, j2
    integer                   ::  j
    integer                   ::  imode
    integer                   ::  n1
    integer                   ::  shp(4)
    real, allocatable         ::  HSSt(:,:)  ! (nobs,n1)
    real, allocatable         ::  Kdt(:,:)   ! (nobs,n1)
    real, allocatable         ::  Ed_c(:,:,:,:,:)   ! (nxd,nyd,nz,nspec,nmodes_all)
    real, allocatable         ::  Sd_c(:,:,:,:,:)   ! (nxd,nyd,nz,nspec,nmodes_all)
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

    ! short:
    dom => self%dom(idom)
    
    ! domain size:
    nxd = dom%lli%nlon
    nyd = dom%lli%nlat
    ! local range:
    i1 = dom%i0+1 ; i2 = dom%i0+nxd
    j1 = dom%j0+1 ; j2 = dom%j0+nyd
    
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
    
    !...............................................
    ! start timing:
    call GO_Timer_Start( itim_exchange_modes, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! input shape:
    select case ( substate )
      case ( substate_c      ) ; shp = (/nxd,nyd,nz,1/)   ! single ispec
      case ( substate_cg     ) ; shp = (/nxd,nyd,1 ,1/)   ! single ispec
      case ( substate_aerh2o ) ; shp = (/nxd,nyd,nz,1    /)
      case ( substate_dc     ) ; shp = (/nxd,nyd,nnoise,nhist/)
      case default
        write (gol,'("unsupported substate ",i0)') substate; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! temporary storage:
    allocate( Ed_c(shp(1),shp(2),shp(3),shp(4),nmodes_all), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Sd_c(shp(1),shp(2),shp(3),shp(4),nmodes_all), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( xd_c(shp(1),shp(2),shp(3),shp(4)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Kd_c(dom%nobs,shp(1),shp(2),shp(3),shp(4)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of elements in 1d:
    n1 = size(xd_c)
    ! storage to solve gain:
    allocate( HSSt(dom%nobs,n1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Kdt(dom%nobs,n1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! init to zero, needed if 'AllReduce' is used to exchange the modes by a sum over the pe's:
    Ed_c = 0.0
    ! loop over local modes:
    do j = 1, nmodes_loc
      ! global mode index:
      imode = imodes(j)
      ! copy local slab:
      select case ( substate )
        case ( substate_c      ) ; Ed_c(:,:,:,1,imode) = Ens(j)%c     (i1:i2,j1:j2,:,ispec)
        case ( substate_cg     ) ; Ed_c(:,:,1,1,imode) = Ens(j)%cg    (i1:i2,j1:j2  ,ispec)
        case ( substate_aerh2o ) ; Ed_c(:,:,:,1,imode) = Ens(j)%aerh2o(i1:i2,j1:j2,:)
        case ( substate_dc     ) ; Ed_c(:,:,:,:,imode) = Ens(j)%dc    (i1:i2,j1:j2,:,:)
        case default
          write (gol,'("unsupported substate ",i0)') substate; call goErr
          TRACEBACK; status=1; return
      end select
    end do ! local modes
    !
    ! sum all local arrays (with many zeros ...) using a 'reduce',
    ! and return the result to all processes (input is the same as output):
    call goc%AllReduce( Ed_c, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ensemble mean:
    xd_c = sum( Ed_c, dim=5 )/nmodes_all
    ! ensemble covariance square root:
    do imode = 1, nmodes_all
      Sd_c(:,:,:,:,imode) = ( Ed_c(:,:,:,:,imode) - xd_c )/sqrt(nmodes_all-1.0)
    end do

    !...............................................
    ! switch timing:
    call GO_Timer_Switch( itim_exchange_modes, itim_solve_gain, status )
    IF_NOTOK_RETURN(status=1)
    !...............................................
    
    ! solve gain from:
    !    [(HS)(HS)'+R] Kd' = (HS)S'
    
    ! fill right hand side:
    HSST = matmul( dom%HSd, transpose(reshape(Sd_c,(/n1,nmodes_all/))) )

    !! solve:
    !call LinAlg_SolveSym( dom%HPHR, HSSt, Kdt, status )
    !IF_NOTOK_RETURN(status=1)

    ! solve using upper-triangular factorization:
    call LinAlg_Sym_FactorSolve( dom%U, HSSt, Kdt, status )
    IF_NOTOK_RETURN(status=1)
    
    ! reshape:
    Kd_c = reshape( Kdt, (/dom%nobs,shp(1),shp(2),shp(3),shp(4)/) )

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
    do j = 1, nmodes_loc
      ! global mode index:
      imode = imodes(j)
      ! update ensemble member:
      !    x := x + K (y-Hx+v)
      ! take care of spatial and chemical correlations:
      do iobs = 1, dom%nobs
        ! residue:
        !   d = y + vj - HXj
        res = dom%y(iobs) + dom%v(iobs,imode) - dom%HX(iobs,imode)
        !print *, 'xxx1 mode ', j,' iobs ', iobs
        !print *, '  x1 y = ', dom%y(iobs), 'hx=', dom%HX(iobs,imode), 'v=', dom%v(iobs,imode)
        !print *, '  x1 res = ', res
        ! current observation set:
        iset = dom%inds(iobs)%iset
        ! switch:
        select case ( substate )
          !~
          case ( substate_c )
            ! analyse?
            if ( self%xcorr(iset)%analyse_spec(ispec) ) then
              ! loop over levels:
              do iz = 1, nz
                ! add contribution, apply spatial localization:
                Ens(j)%c(i1:i2,j1:j2,iz,ispec) = Ens(j)%c(i1:i2,j1:j2,iz,ispec) + &
                             Kd_c(iobs,:,:,iz,1) * dom%corr(iobs,:,:) * res
              end do ! iz
            end if ! analyse spec?
          !~
          case ( substate_cg )
            ! analyse?
            if ( self%xcorr(iset)%analyse_spec(ispec) ) then
              ! add contribution, apply spatial localization:
              Ens(j)%cg(i1:i2,j1:j2,ispec) = Ens(j)%cg(i1:i2,j1:j2,ispec) + &
                           Kd_c(iobs,:,:,1,1) * dom%corr(iobs,:,:) * res
            end if ! analyse spec?
          !~
          case ( substate_aerh2o )
            ! analyse sia tracers? then also update aerh2o:
            if ( self%xcorr(iset)%analyse_sia ) then
              ! loop over levels:
              do iz = 1, nz
                ! add contribution, apply spatial localization:
                Ens(j)%aerh2o(i1:i2,j1:j2,iz) = Ens(j)%aerh2o(i1:i2,j1:j2,iz) + &
                               Kd_c(iobs,:,:,iz,1) * dom%corr(iobs,:,:) * res
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
                  Ens(j)%dc(i1:i2,j1:j2,inoise,ihist) = Ens(j)%dc(i1:i2,j1:j2,inoise,ihist) + &
                                 Kd_c(iobs,:,:,inoise,ihist) * dom%corr(iobs,:,:) * res
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


  subroutine AnDom_Init( self, west , dlon, nlon, i0_lon, &
                               south, dlat, nlat, i0_lat, status )

    use Grid   , only : Init

    ! --- in/out -------------------------

    class(T_AnDom), intent(out)     ::  self
    real, intent(in)                ::  west
    real, intent(in)                ::  dlon
    integer, intent(in)             ::  nlon
    integer, intent(in)             ::  i0_lon
    real, intent(in)                ::  south
    real, intent(in)                ::  dlat
    integer, intent(in)             ::  nlat
    integer, intent(in)             ::  i0_lat
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Init'

    ! --- local --------------------------
    
    ! --- begin ---------------------------

    ! offset:
    self%i0 = i0_lon
    self%j0 = i0_lat

    ! init sub domain:
    call Init( self%lli, west , dlon, nlon, &
                         south, dlat, nlat, status )
    IF_NOTOK_RETURN(status=1)
    
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

    ! done:
    call Done( self%lli, status )
    IF_NOTOK_RETURN(status=1)
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
    use LEKF_State, only : nmodes_all, nmodes_loc, imodes

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(in)             ::  inds(:)  ! 3: iset, iloc, obs_ivar
    logical, intent(in)             ::  inside
    real, intent(in)                ::  lon, lat
    real, intent(in)                ::  rho_m ! m
    real, intent(in)                ::  y, r
    real, intent(in)                ::   v_loc(:) ! (nmodes_loc)
    real, intent(in)                ::  Hx_loc(:) ! (nmodes_loc)
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
    integer                   ::  imode
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
      allocate( self%v (next,nmodes_all), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%HX(next,nmodes_all), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%corr(next,self%lli%nlon,self%lli%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! help array
      allocate( self%dist(self%lli%nlon,self%lli%nlat), stat=status )
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
        allocate( new_v   (self%nobs+next,nmodes_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_HX  (self%nobs+next,nmodes_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( new_corr(self%nobs+next,self%lli%nlon,self%lli%nlat), stat=status )
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
    
    ! loop over local modes:
    do j = 1, nmodes_loc
      ! global mode index:
      imode = imodes(j)
      ! copy elements:
      self%v (self%nobs,imode) =  v_loc(j)
      self%HX(self%nobs,imode) = HX_loc(j)
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
    call RoundToResolution( self%lli, lon, lat, lon0, lat0, status )!, debug=debug )
    IF_NOTOK_RETURN(status=1)
    ! fill localization field with distances to rounded location:
    call DistanceGrid( self%lli, lon0, lat0, self%dist, status )!, debug=debug )
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


  subroutine AnDom_BCast_Modes( self, status )
  
    use GO        , only : goc
    use LEKF_State, only : nmodes_all, imode_pe

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_BCast_Modes'
    
    ! --- local ----------------------------------

    integer                   ::  j
    integer                   ::  imode

    ! --- begin ----------------------------------
    
    ! any observations?
    if ( self%nobs > 0 ) then

      ! loop over global modes:
      do imode = 1, nmodes_all
        ! broadcast from pe holding this mode,
        ! will be received by all others:
        call goc%BCast( imode_pe(imode), self%HX(:,imode), status )
        IF_NOTOK_RETURN(status=1)
      end do ! modes
      
    end if  ! any obs

    ! ok
    status = 0

  end subroutine AnDom_BCast_Modes


  ! ***


  subroutine AnDom_Fill_HPHR_etc( self, status )
  
    use GO        , only : goc
    use Num       , only : LinAlg_Sym_Factorize
    use Grid      , only : ll_distance
    use LEKF_State, only : nmodes_all, imode_pe

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
    
      character(len=1024)   ::  dumpfile

    ! --- begin ----------------------------------
    
    ! output storage:
    allocate( self%HSd(self%nobs,nmodes_all), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%HPHR(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%U(self%nobs,self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! temporary storage:
    allocate( Hxd(self%nobs), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ensemble mean:
    Hxd = sum( self%Hx(1:self%nobs,:), dim=2 )/nmodes_all
    ! ensemble covariance square root:
    do imode = 1, nmodes_all
      self%HSd(:,imode) = ( self%Hx(1:self%nobs,imode) - Hxd )/sqrt(nmodes_all-1.0)
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
    
    ! add R:
    do iobs = 1, self%nobs
      self%HPHR(iobs,iobs) = self%HPHR(iobs,iobs) + self%r(iobs)**2
    end do
    
    ! factorize:
    call LinAlg_Sym_Factorize( self%HPHR, self%U, status )
    if ( status /= 0 ) then
      ! target file:
      write (dumpfile,'("HPHR__pe",i0,".nc")') goc%myid
      ! info ...
      write (gol,'("Factorization failed, write matrix to:")'); call goErr
      write (gol,'("  ",a)') trim(dumpfile); call goErr
      ! dump ...
      call nc_dump( trim(dumpfile), self%HPHR, 'HPHR', (/'nobs','nobs2'/), status )
      IF_NOTOK_RETURN(status=1)
      ! leave:
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine AnDom_Fill_HPHR_etc


  ! ***


  subroutine AnDom_Clear_HPHR_etc( self, status )
  
    use LEKF_State, only : nmodes_all

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

