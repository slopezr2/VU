!###############################################################################
!
! Measurement update per domain
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

module LEKF_Meas_AnDom

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_StrError, NF90_NOERR

  use Grid, only : TllGridInfo

  implicit none
  
  
  ! --- in/out -----------------------------

  private

  public    ::  rho_frac_zero
  public    ::  T_AnDom
  

  ! --- const --------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_AnDom'
  
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
    ! selected nearby observations:
    integer                       ::  nobs
    ! reference to observation data:
    type(T_ObsInd), pointer       ::  inds(:)  ! (>= n)
    ! observation and error stdv:
    real, pointer                 ::  y(:)      ! (nobs)
    real, pointer                 ::  r(:)      ! (nobs)
    real, pointer                 ::  HX(:,:)   ! (nobs,nmodes)
    ! correlation length scale:
    real, pointer                 ::  rho(:)    ! (nobs) [m]
    ! random errors?
    logical                       ::  with_v
    real, pointer                 ::  v(:,:)    ! (nobs,nmodes)
    !
    ! enable rho_m, corr, dist?
    logical                       ::  with_rdc
    ! correlation length scale for locaization,
    ! only one value supported yet:
    real                          ::  an_rho_m   ! m
    ! spatial localization weight:
    real, pointer                 ::  an_corr(:,:,:)  ! (nobs,nlon,nlat)
    ! help array with distances:
    real, pointer                 ::  an_dist(:,:)  ! (nlon,nlat)
    !
    ! EnKF update matrices:
    real, allocatable             ::  HSd(:,:)   ! (nobs,nmodes)
    real, allocatable             ::  HPHR(:,:)  ! (nobs,nobs)
    real, allocatable             ::  U(:,:)     ! (nobs,nobs)
    ! LETKF update entities:
    integer                       ::  nobsl
    real, allocatable             ::  yo(:)    ! (nobs)
    real, allocatable             ::  ybar(:)    ! (nobs)
    real, allocatable             ::  YY(:,:)    ! (nobs,nmodes)
    real, allocatable             ::  CC(:,:)    ! (nmodes,nobs)
    real, allocatable             ::  invPa(:,:) ! (nmodes,nmodes)
    real, allocatable             ::  Pa(:,:)    ! (nmodes,nmodes)
    real, allocatable             ::  lambda(:)   ! (nmodes)
    real, allocatable             ::  QQ(:,:)    ! (nmodes,nmodes)
    real, allocatable             ::  WWa(:,:)   ! (nmodes,nmodes)
    real, allocatable             ::  wa(:)      ! (nmodes)
  contains
    procedure ::  Init                  =>  AnDom_Init
    procedure ::  Done                  =>  AnDom_Done
    procedure ::  Reset                 =>  AnDom_Reset
    procedure ::  Add_Observation       =>  AnDom_Add_Observation
    procedure ::  Exchange_Nearby       =>  AnDom_Exchange_Nearby
    !procedure ::  BCast_Modes           =>  AnDom_BCast_Modes
    procedure ::  Fill_HPHR_etc         =>  AnDom_Fill_HPHR_etc
    procedure ::  Clear_HPHR_etc        =>  AnDom_Clear_HPHR_etc
    procedure ::  Fill_Y_etc            =>  AnDom_Fill_Y_etc
    procedure ::  Clear_Y_etc           =>  AnDom_Clear_Y_etc
    procedure ::  Dump                  =>  AnDom_Dump
  end type T_AnDom


  ! --- interfaces ------------------------
  
  interface HCorrGaussian
    module procedure HCorrGaussian_0d
    module procedure HCorrGaussian_2d
  end interface HCorrGaussian


contains


  ! ********************************************************************
  ! ***
  ! *** analysis subdomains
  ! ***
  ! ********************************************************************


  subroutine AnDom_Init( self, with_v, with_rdc, status )

    ! --- in/out -------------------------

    class(T_AnDom), intent(out)     ::  self
    logical                         ::  with_v
    logical                         ::  with_rdc
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
    nullify( self%HX   )
    nullify( self%rho  )
    
    ! store:
    self%with_rdc = with_rdc
    ! safety ...
    nullify( self%an_corr )
    nullify( self%an_dist )
    
    ! store:
    self%with_v = with_v
    ! safety ...
    nullify( self%v    )
    
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
      deallocate( self%rho , stat=status )
      IF_NOTOK_RETURN(status=1)
      ! random errors?
      if ( self%with_v ) then
        deallocate( self%v   , stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! simulations:
      deallocate( self%HX  , stat=status )
      IF_NOTOK_RETURN(status=1)
      ! correlations:
      if ( self%with_rdc ) then
        deallocate( self%an_corr, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%an_dist, stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
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
                                       y, r, Hx_loc, status, &
                                       v_loc )!, debug )
  
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
    real, intent(in)                ::  Hx_loc(:) ! (nmodes)
    integer, intent(out)            ::  status
    
    real, intent(in), optional      ::   v_loc(:) ! (nmodes)
    !logical, intent(in), optional   ::  debug

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Add_Observation'
    
    ! extra number of observations allocated:
    integer, parameter      ::  next = 1000

    ! --- local ----------------------------------

    !logical                   ::  dbg
    type(T_ObsInd), pointer   ::  new_inds(:)
    real, pointer             ::  new_y(:)
    real, pointer             ::  new_r(:)
    real, pointer             ::  new_rho(:)
    real, pointer             ::  new_v (:,:)
    real, pointer             ::  new_HX(:,:)
    real, pointer             ::  new_corr(:,:,:)

    integer                   ::  j
    real                      ::  lon0, lat0
    
    ! --- begin ----------------------------------
    
    !! testing ...
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! check ...
    if ( self%with_v .and. (.not. present(v_loc)) ) then
      write (gol,'("analysis requires random errors, but `v_loc` not supplied")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! new storage needed?
    if ( .not. associated(self%inds) ) then
      ! storage:
      allocate( self%inds(next), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%y(next)   , stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%r(next)   , stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%rho(next) , stat=status )
      IF_NOTOK_RETURN(status=1)
      ! random errors?
      if ( self%with_v ) then
        allocate( self%v (next,nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! simulations:
      allocate( self%HX(next,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! rho,dist,corr?
      if ( self%with_rdc ) then
        ! correlations:
        allocate( self%an_corr(next,ugg%nlon,ugg%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! help array
        allocate( self%an_dist(ugg%nlon,ugg%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
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
        allocate( new_rho (self%nobs+next), stat=status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_v ) then
          allocate( new_v   (self%nobs+next,nmodes), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        allocate( new_HX  (self%nobs+next,nmodes), stat=status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_rdc ) then
          allocate( new_corr(self%nobs+next,ugg%nlon,ugg%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! only local modes are filled, ensure that the rest is zero:
        new_Hx = 0.0
        ! copy:
        new_inds(1:self%nobs)       = self%inds
        new_y   (1:self%nobs)       = self%y
        new_r   (1:self%nobs)       = self%r
        new_rho (1:self%nobs)       = self%rho
        if ( self%with_v ) then
          new_v   (1:self%nobs,:)   = self%v
        end if
        new_HX  (1:self%nobs,:)     = self%HX
        if ( self%with_rdc ) then
          new_corr(1:self%nobs,:,:) = self%an_corr
        end if
        ! clear old:
        deallocate( self%inds, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%y   , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%r   , stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( self%rho , stat=status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_v ) then
          deallocate( self%v   , stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        deallocate( self%HX  , stat=status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_rdc ) then
          deallocate( self%an_corr, stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! re-assign:
        self%inds => new_inds
        self%y    => new_y
        self%r    => new_r
        self%rho  => new_rho
        if ( self%with_v ) then
          self%v    => new_v
        end if
        self%HX   => new_HX
        if ( self%with_rdc ) then
          self%an_corr => new_corr
        end if
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
    ! correlation length scale:
    self%rho (self%nobs)          = rho_m
    
    ! loop over modes:
    do j = 1, nmodes
      !! global mode index:
      !imode = imodes(j)
      ! copy elements:
      if ( self%with_v ) then
        self%v (self%nobs,j) =  v_loc(j)
      end if
      self%HX(self%nobs,j) = HX_loc(j)
    end do ! j
    
    ! round to grid cell resolution, if inside domain this is a grid cell center:
    call ugg%RoundToResolution( lon, lat, lon0, lat0, status )!, debug=debug )
    IF_NOTOK_RETURN(status=1)
    ! store:
    self%inds(self%nobs)%lon0 = lon0
    self%inds(self%nobs)%lat0 = lat0
    
    ! rho, dist, corr?
    if ( self%with_rdc ) then
    
      ! single rho for entire set ...
      if ( self%nobs == 1 ) then
        self%an_rho_m = rho_m
      else if ( rho_m /= self%an_rho_m ) then
        write (gol,'("only single localization length scale supported yet;")'); call goErr
        write (gol,'("  current value : ",f12.2," m")') self%an_rho_m; call goErr
        write (gol,'("  new value     : ",f12.2," m")') rho_m; call goErr
        TRACEBACK; status=1; return
      end if

      ! fill localization field with distances to rounded location:
      call ugg%DistanceGrid( lon0, lat0, self%an_dist, status )!, debug=debug )
      IF_NOTOK_RETURN(status=1)

      !! testing ...
      !if ( dbg ) then
      !  print *, '  yyy1 added nobs ', self%nobs
      !  print *, '    y1 inds ', self%inds(self%nobs)%iset, self%inds(self%nobs)%iloc, self%inds(self%nobs)%obs_ivar
      !  print *, '    y1 inside ', self%inds(self%nobs)%inside
      !  print *, '    y1 loc  ', self%inds(self%nobs)%lon, self%inds(self%nobs)%lat
      !  print *, '    y1 loc0 ', lon0, lat0
      !  print *, '    y1 dat  ', self%y(self%nobs), self%r(self%nobs)
      !  print *, '    y1 v    ', self%v (self%nobs,:)
      !  print *, '    y1 HX   ', self%HX(self%nobs,:)
      !  print *, '    y1 dist ', self%an_corr(self%nobs,:,:)
      !end if

      ! convert to correlation:
      call HCorrGaussian( self%an_rho_m, self%an_dist, self%an_corr(self%nobs,:,:), status )
      IF_NOTOK_RETURN(status=1)

      !! testing ...
      !if ( dbg ) then
      !  print *, '    y1 corr ', self%an_corr(self%nobs,:,:)
      !end if
      
    end if ! rdc

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
    real, allocatable     ::  v_loc(:)         ! (nmodes)
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
    
    ! storage:
    allocate( v_loc(nmodes), stat=status )
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
                                     rho_frac_zero*self%rho(iobs), selected(iobs,id), status )
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
    if ( self%with_v ) then
      i0v   = nr ; nr = nr + nmodes
    else
      i0v   = -999
    end if
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
                rbuffer(i0rho+1            ,k) = self%rho(iobs)
                rbuffer(i0y  +1            ,k) = self%y(iobs)
                rbuffer(i0r  +1            ,k) = self%r(iobs)
                if ( self%with_v ) then
                  rbuffer(i0v  +1:i0v +nmodes,k) = self%v(iobs,:)
                end if
                rbuffer(i0hx +1:i0hx+nmodes,k) = self%HX(iobs,:)
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
              ! random errors?
              if ( self%with_v ) then
                v_loc = rbuffer(i0v  +1:i0v +nmodes,k)
              end if
              ! add observation, not inside:
              call self%Add_Observation( ibuffer(:,k), .false., &
                                         rbuffer(i0lon+1            ,k), &
                                         rbuffer(i0lat+1            ,k), &
                                         rbuffer(i0rho+1            ,k), &
                                         rbuffer(i0y  +1            ,k), &
                                         rbuffer(i0r  +1            ,k), &
                                         rbuffer(i0hx +1:i0hx+nmodes,k), &
                                         status, &
                                         v_loc=v_loc )
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
    deallocate( v_loc, stat=status )
    IF_NOTOK_RETURN(status=1)

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
    !character(len=1024)       ::  dumpfile

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
    ! ensemble covariance square root:
    do imode = 1, nmodes
      self%HSd(:,imode) = ( self%Hx(1:self%nobs,imode) - Hxd )/sqrt(nmodes-1.0)
    end do

    ! clear:
    deallocate( Hxd, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill HPH:
    self%HPHR = matmul( self%HSd, transpose(self%HSd) )
    
    ! safety check ...
    if ( .not. self%with_rdc ) then
      write (gol,'("rho/dist/corr arrays not enabled")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! apply localization;
    ! loop over first observations:
    do iobs = 1, self%nobs
      ! loop over remaining observations:
      do iobs2 = iobs, self%nobs
        ! distance in m:
        d = ll_distance( self%inds(iobs )%lon, self%inds(iobs )%lat, &
                         self%inds(iobs2)%lon, self%inds(iobs2)%lat  )
        ! correlation:
        call HCorrGaussian( self%an_rho_m, d, c, status )
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
      ! info ...
      write (gol,'("Factorization failed, dump matrices ...")'); call goErr
      ! dump HPHR etc:
      call self%Dump( status )
      IF_NOTOK_RETURN(status=1)
      
      !! target file:
      !write (dumpfile,'("HPHR__pe",i0,".nc")') goc%id
      !! info ...
      !write (gol,'("Factorization failed, write matrix to:")'); call goErr
      !write (gol,'("  ",a)') trim(dumpfile); call goErr
      !! dump ...
      !call nc_dump( trim(dumpfile), self%HPHR, 'HPHR', (/'nobs ','nobs2'/), status )
      !IF_NOTOK_RETURN(status=1)

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
        !call nc_dump( trim(dumpfile), self%HPHR, 'HPHR', (/'nobs ','nobs2'/), status )
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


  subroutine AnDom_Fill_Y_etc( self, lon, lat, xcorr, status, &
                                  ispec, sia, inoise )!, debug )
  
    use Num       , only : LinAlg_Pod_Eigen
    use Grid      , only : ll_distance
    use LEKF_State     , only : nmodes
    use LEKF_Meas_Tools, only : T_XCorr

    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    real, intent(in)                ::  lon, lat     ! grid cell center (lon,lat)
    type(T_XCorr), intent(in)       ::  xcorr(:)     ! (nset) flags for correlation between observation set
    integer, intent(out)            ::  status
    
    integer, intent(in), optional   ::  ispec   ! only observations to update tracer specie
    logical, intent(in), optional   ::  sia     ! only observations to update sia component
    integer, intent(in), optional   ::  inoise  ! only observations to update noise element
    !logical, intent(in), optional   ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Fill_Y_etc'
    
    ! --- local ----------------------------------

    !logical                   ::  dbg
    integer                   ::  iobs
    integer                   ::  iset
    real                      ::  d
    real                      ::  cdist
    integer                   ::  imode
    
    !! testing ...
    !real, allocatable         ::  A(:,:)

    ! --- begin ----------------------------------
    
    !! testing ...
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! storage for maximum number of observations:
    if ( .not. allocated(self%YY) ) then
      allocate( self%yo(self%nobs), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%ybar(self%nobs), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%YY(self%nobs,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%CC(nmodes,self%nobs), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%invPa(nmodes,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%Pa(nmodes,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%lambda(nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%QQ(nmodes,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%WWa(nmodes,nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%wa(nmodes), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! no local observations yet:
    self%nobsl = 0
    
    ! loop:
    do iobs = 1, self%nobs
      ! current observation set:
      iset = self%inds(iobs)%iset
      ! skip if no correlation with state elements:
      if ( present(ispec) ) then
        if ( .not. xcorr(iset)%analyse_spec(ispec) ) cycle
      end if
      if ( present(sia) ) then
        if ( sia .and. (.not. xcorr(iset)%analyse_sia) ) cycle
      end if
      if ( present(inoise) ) then
        if ( .not. xcorr(iset)%analyse_noise(inoise) ) cycle
      end if
      ! distance in m:
      d = ll_distance( self%inds(iobs )%lon, self%inds(iobs )%lat, lon, lat )
      ! within 3.5*rho ?
      if ( d <= rho_frac_zero*self%rho(iobs) ) then
        ! increase counter:
        self%nobsl = self%nobsl + 1
        ! weight with distance:
        call HCorrGaussian( self%rho(iobs), d, cdist, status )
        IF_NOTOK_RETURN(status=1)
        ! observation:
        self%yo(self%nobsl) = self%y(iobs)
        ! average simulation over modes:
        self%ybar(self%nobsl) = sum( self%HX(iobs,:) )/nmodes
        ! fill simulation perturbation:
        !   Y = [ ..., Hxi - <Hxi>, ..]
        self%YY(self%nobsl,:) = self%HX(iobs,:) - self%ybar(self%nobsl)
        ! fill:
        !   C = Y^T R^{-1} o Lc
        self%CC(:,self%nobsl) = self%YY(self%nobsl,:) / self%r(iobs)**2 * cdist
        !! testing ...
        !if ( dbg ) then
        !  write (gol,*) '..y    select iobs ', iobs, '; yo=', self%yo(self%nobsl), &
        !                  '; ybar=', self%ybar(self%nobsl), '; Y=', self%YY(self%nobsl,:), &
        !                  '; cdist=', cdist; call goPr
        !end if
      end if ! within range
    end do ! iobs
    
    ! any local observations?
    if ( self%nobsl > 0 ) then
    
      ! fill inverse of Pa:
      !    Pa^{-1} = [ C Y + (m-1) I ]
      self%invPa = matmul( self%CC(:,1:self%nobsl), self%YY(1:self%nobsl,:) )
      ! add diagonal:
      do imode = 1, nmodes
        self%invPa(imode,imode) = self%invPa(imode,imode) + (nmodes-1)
      end do
    
      ! eigenvalue decomposition of inverse of analyzed covariance in ensemble space:
      !   Pa^{-1} = Q Lambda Q^T
      call LinAlg_Pod_Eigen( self%invPa, self%lambda, self%QQ, status )
      if ( status /= 0 ) then
        ! info ...
        write (gol,'("eigenvalue decomposition failed ...")'); call goErr
        TRACEBACK; status=1; return
      end if ! eigen
      
      ! eigenvalue decomposition of inverse of analyzed covariance in ensemble space,
      !   Pa^{-1} = Q Lambda^{-1} Q^T
      ! where the eigenvector matrix Q is orthonormal:
      !   Q^T Q = I, thus Q^{-1} = Q^T
      do imode = 1, nmodes
        self%Pa(:,imode) = self%QQ(:,imode) / self%lambda(imode)
      end do
      self%Pa = matmul( self%Pa, transpose(self%QQ) )
      
      ! Symetric square roots W of (m-1)Pa :
      !    (m-1) Pa = (m-1) Q Lambda^{-1} Q^T 
      !       W W^T = (m-1) Q Lambda^{-1/2} Q^T Q Lambda^{-1/2} Q^T
      !           W = sqrt(m-1) Q Lambda^{-1/2} Q^T
      do imode = 1, nmodes
        self%WWa(:,imode) =  sqrt(nmodes-1.0) * self%QQ(:,imode) / sqrt(self%lambda(imode))
      end do
      self%WWa = matmul( self%WWa, transpose(self%QQ) )
            
      ! analysis weights:
      !   w = Pa C ( y - ybar )
      self%wa = matmul( self%Pa, matmul( self%CC(:,1:self%nobsl), ( self%yo(1:self%nobsl) - self%ybar(1:self%nobsl) ) ) )
      !! testing ...
      !if ( dbg ) then
      !  write (gol,*) 'wa = ', self%wa ; call goPr
      !end if

      !! testing ...
      !write (gol,*) ''; call goPr
      !write (gol,*) 'Pa =   (eigenvalue ',1.0/(nmodes-1),')'; call goPr
      !do imode = 1, nmodes
      !  write (gol,*) self%Pa(imode,:), ' ; ', sum(self%Pa(imode,:)); call goPr
      !end do
      !write (gol,*) ''; call goPr
      !write (gol,*) 'lambda = '; call goPr
      !write (gol,*) self%lambda; call goPr
      !write (gol,*) ''; call goPr
      !write (gol,*) 'Q = '; call goPr
      !do imode = 1, nmodes
      !  write (gol,*) self%QQ(imode,:), ' ; ', sum(self%QQ(imode,:)); call goPr
      !end do
      !write (gol,*) ''; call goPr
      !write (gol,*) 'W =   (eigenvalue 1.0)'; call goPr
      !do imode = 1, nmodes
      !  write (gol,*) self%WWa(imode,:), ' ; ', sum(self%WWa(imode,:)); call goPr
      !end do
      !allocate( A(nmodes,nmodes), stat=status )
      !IF_NOTOK_RETURN(status=1)
      !write (gol,*) ''; call goPr
      !write (gol,*) 'abs( W W^T - (m-1)Pa ) = '; call goPr
      !A = matmul( self%WWa, transpose(self%WWa) )
      !do imode = 1, nmodes
      !  write (gol,*) abs( A(imode,:) - (nmodes-1)*self%Pa(imode,:) ); call goPr
      !end do
      !write (gol,*) ''; call goPr
      !deallocate( A, stat=status )
      !IF_NOTOK_RETURN(status=1)
      
    end if  ! nobsl > 0
    
    ! ok
    status = 0

  end subroutine AnDom_Fill_Y_etc


  ! ***


  subroutine AnDom_Clear_Y_etc( self, status )
  
    ! --- in/out ---------------------------------

    class(T_AnDom), intent(inout)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AnDom_Clear_Y_etc'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! output storage:
    deallocate( self%yo, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%ybar, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%YY, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%CC, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%invPa, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%Pa, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%lambda, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%QQ, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%WWa, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%wa, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AnDom_Clear_Y_etc


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
    integer                   ::  varid_rho
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
      status = NF90_Def_Var( ncid, 'rho'   , NF90_FLOAT, (/dimid_obs/), varid_rho )
      IF_NF90_NOTOK_RETURN(status=1)
      if ( self%with_v ) then
        status = NF90_Def_Var( ncid, 'v'     , NF90_FLOAT, (/dimid_obs,dimid_mode/), varid_v )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      status = NF90_Def_Var( ncid, 'HX'    , NF90_FLOAT, (/dimid_obs,dimid_mode/), varid_HX )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variable:
      if ( self%with_rdc ) then
        status = NF90_Def_Var( ncid, 'dist'     , NF90_FLOAT, (/dimid_lon,dimid_lat/), varid_dist )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Var( ncid, 'corr'     , NF90_FLOAT, (/dimid_obs,dimid_lon,dimid_lat/), varid_corr )
        IF_NF90_NOTOK_RETURN(status=1)
      end if

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
      status = NF90_Put_Var( ncid, varid_rho, self%rho(1:self%nobs) )
      IF_NF90_NOTOK_RETURN(status=1)
      if ( self%with_v ) then
        status = NF90_Put_Var( ncid, varid_v, self%v(1:self%nobs,:) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      status = NF90_Put_Var( ncid, varid_HX, self%HX(1:self%nobs,:) )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! write variable:
      if ( self%with_rdc ) then
        status = NF90_Put_Var( ncid, varid_dist, self%an_dist )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncid, varid_corr, self%an_corr(1:self%nobs,:,:) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      
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
  
  
  ! ********************************************************************
  ! ***
  ! *** correlation functions
  ! ***
  ! ********************************************************************

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


end module LEKF_Meas_AnDom

