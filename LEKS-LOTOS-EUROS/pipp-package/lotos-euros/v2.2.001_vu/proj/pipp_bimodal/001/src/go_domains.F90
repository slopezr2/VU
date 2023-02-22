!#######################################################################
!
! Data domain decomposition tools.
!
! Initialization and finalization
!
!   ! use communicator object:
!   use GO_Comm, only : goc
!
!    integer              ::  status
!    character(len=1024)  ::  msg
!    integer              ::  lenmsg
!
!    ! init communication:
!    call goc%Init( status )
!    if (status/=0) stop
!
!    ! init module, provide communicator:
!    call GO_Par_Init( MPI_COMM_WORLD, status )
!    if (status/=0) stop
!
!    ! ...
!
!    ! done with module:
!    call GO_Par_Done( status )
!    if (status/=0) stop
!
!  #ifdef _MPI
!    ! done with communication:
!    call MPI_Finalize( status )
!    if ( status /= MPI_SUCCESS ) then
!      call MPI_Error_String( status, msg, lenmsg, status )
!      write (*,'("Error from MPI_Done: ",a)') trim(msg)
!    end if
!  #endif
!
!
! Preprocessor macro's
!
!  _MPI   : should be defined to enable MPI code, e.g. 'f90 -D_MPI ...'
!
!
! Manuals
!
! - MPI 3.0 manual:
!     http://mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; errorcode=status; call MPI_Error_String(errorcode,gol,ngol,status); call goErr; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#######################################################################


module GO_Domains

  use GO_Print, only : ngol, gol, goPr, goErr
#ifdef _MPI
  use MPI_F08, only : MPI_SUCCESS, MPI_Error_String
  use MPI_F08, only : MPI_COMM
#endif
#ifdef with_netcdf
  use NetCDF, only : NF90_NOERR, nf90_strerror
#endif

  implicit none


  ! --- in/out ------------------------

  private

  public  ::  T_Domains


  ! --- const --------------------------

  character(len=*), parameter   ::  mname = 'GO_Domains'
  
  ! directions for data shifts:
  integer, parameter      ::  ndir = 2
  ! shifts per direction:
  integer, parameter      ::  dir_shift(ndir) = (/-1,+1/)


  ! --- local  --------------------------
  
!  ! mpi variables:
!  integer             ::  comm
!  integer             ::  nproc, me
!  ! flag:
!  logical             ::  module_setup


  ! --- types -----------------------------------
  
  ! domain ranges from all process to facilitate decomposition
  type T_Domains
    ! rank:
    integer                ::  ndim
    ! global bounds and shape:
    integer, allocatable   ::  glb_lbo(:)  ! (ndim)
    integer, allocatable   ::  glb_ubo(:)  ! (ndim)
    integer, allocatable   ::  glb_shp(:)  ! (ndim)
    ! global bounds:
    integer, allocatable   ::  glbo(:,:)  ! (ndim,0:nproc-1)
    integer, allocatable   ::  gubo(:,:)  ! (ndim,0:nproc-1)
    ! local bounds:
    integer, allocatable   ::  lbo(:,:)   ! (ndim,0:nproc-1)
    integer, allocatable   ::  ubo(:,:)   ! (ndim,0:nproc-1)
    ! offset:
    integer, allocatable   ::  off(:,:)   ! (ndim,0:nproc-1)
    ! shape:
    integer, allocatable   ::  shp(:,:)   ! (ndim,0:nproc-1)
    ! total number of local elements:
    integer, allocatable   ::  n(:)   ! (0:nproc-1)
    ! cartesian processor topology;
#ifdef _MPI
    !~ dedicated communicator:
    type(MPI_Comm)         ::  cart_comm
#endif
    !~ info on layout:
    integer, allocatable   ::  cart_dims(:)       ! (ndim)
    logical, allocatable   ::  cart_periods(:)    ! (ndim)
    integer, allocatable   ::  cart_coords(:)     ! (ndim)    
    integer, allocatable   ::  cart_src_id(:,:)   ! (ndim,ndir)
    integer, allocatable   ::  cart_dst_id(:,:)   ! (ndim,ndir)
    !
  contains
    procedure   ::  Init            => Domains_Init
    procedure   ::  Done            => Domains_Done
    procedure   ::  Get             => Domains_Get
    procedure   ::  GetDim          => Domains_GetDim
    procedure   ::  Inside          => Domains_Inside
    procedure   ::  Find            => Domains_Find
    procedure   ::  Intersection    => Domains_Intersection
    !
    procedure   ::  ExchangeHalo       => Domains_ExchangeHalo_3d_r
    !
    procedure   ::                        Domains_ExchangeHaloX_3d_l
    procedure   ::                        Domains_ExchangeHaloX_3d_r
    procedure   ::                        Domains_ExchangeHaloX_4d_r
    generic     ::  ExchangeHaloX      => Domains_ExchangeHaloX_3d_l, &
                                          Domains_ExchangeHaloX_3d_r, &
                                          Domains_ExchangeHaloX_4d_r
    !
    procedure   ::                        Domains_ExchangeHaloY_3d_l
    procedure   ::                        Domains_ExchangeHaloY_3d_r
    procedure   ::                        Domains_ExchangeHaloY_4d_r
    generic     ::  ExchangeHaloY      => Domains_ExchangeHaloY_3d_l, &
                                          Domains_ExchangeHaloY_3d_r, &
                                          Domains_ExchangeHaloY_4d_r
    !
    procedure   ::  TransferSlab       => Domains_TransferSlab_3d_r
    procedure   ::  TransferData       => Domains_TransferData_3d_r
    procedure   ::  TransferDataAll    => Domains_TransferDataAll_3d_r
    !
    procedure   ::                        ExchangeTest
    !
    procedure   ::                        Domains_Swap_2d_r8
    procedure   ::                        Domains_Swap_3d_r8
    procedure   ::                        Domains_Swap_4d_r8
    generic     ::  Swap               => Domains_Swap_2d_r8, &
                                          Domains_Swap_3d_r8, &
                                          Domains_Swap_4d_r8
    !
    procedure   ::                        Domains_GatherV_1d_i
    procedure   ::                        Domains_GatherV_1d_i_multiple
    procedure   ::                        Domains_GatherV_1d_r
    procedure   ::                        Domains_GatherV_1d_r_multiple
    procedure   ::                        Domains_GatherV_1d_r_multiple2d
    procedure   ::                        Domains_GatherV_2d_r
    procedure   ::                        Domains_GatherV_3d_r
    procedure   ::                        Domains_GatherV_4d_r
    generic     ::  GatherV            => Domains_GatherV_1d_i, &
                                          Domains_GatherV_1d_i_multiple, &
                                          Domains_GatherV_1d_r, &
                                          Domains_GatherV_1d_r_multiple, &
                                          Domains_GatherV_1d_r_multiple2d, &
                                          Domains_GatherV_2d_r, &
                                          Domains_GatherV_3d_r, &
                                          Domains_GatherV_4d_r
    !
    procedure   ::                        Domains_AllGather_1d_i4
    procedure   ::                        Domains_AllGather_1d_r8
    procedure   ::                        Domains_AllGather_3d_c8
    generic     ::  AllGather          => Domains_AllGather_1d_i4, &
                                          Domains_AllGather_1d_r8, &
                                          Domains_AllGather_3d_c8
    !
    procedure   ::                        Domains_ScatterV_3d_r
    procedure   ::                        Domains_ScatterV_2d_r
    generic   ::  ScatterV           => Domains_ScatterV_3d_r, &
                                                         Domains_ScatterV_2d_r      
    !
    procedure   ::  Extract            => Domains_Extract_1d_r8
    !
    procedure   ::                        Domains_Put_Var_2d
    procedure   ::                        Domains_Put_Var_3d
    procedure   ::                        Domains_Put_Var_4d
    generic     ::  Put_Var            => Domains_Put_Var_2d, &
                                          Domains_Put_Var_3d, &
                                          Domains_Put_Var_4d
    !
  end type T_Domains


  ! --- var -------------------------------------
  
  ! mpi error code; used as argument for 'MPI_Error_String' to avoid
  ! warnings about same argument 'status' being used for both 'errorcode' and 'ierror':
  integer                           ::  errorcode


  
contains


  ! ********************************************************************
  ! ***
  ! *** Domains
  ! ***
  ! ********************************************************************

    !
  ! initialize subdomains as ND array; provide:
  ! - processor topology 
  ! - global lower and upper bounds
  ! local grid sizes are set as equal as possible
  !

  subroutine Domains_Init( self, cart_dims, glb_lbo, glb_ubo, status )

    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_AllGather
    use MPI_F08, only : MPI_INTEGER
    use MPI_F08, only : MPI_Cart_Create, MPI_Cart_Get, MPI_Cart_Rank
#endif
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(out)     ::  self
    integer, intent(in)               ::  cart_dims(:)  ! (ndim)
    integer, intent(in)               ::  glb_lbo(:)  ! (ndim)
    integer, intent(in)               ::  glb_ubo(:)  ! (ndim)
    integer, intent(out)              ::  status

    ! --- local ----------------------------------
    
    integer                 ::  idim
    integer                 ::  idir
    integer                 ::  idom
    integer                 ::  gn
    integer, allocatable    ::  glbo(:)     ! (ndim)
    integer, allocatable    ::  gubo(:)     ! (ndim)
    integer, allocatable    ::  nn(:)
    integer                 ::  i, j
    integer                 ::  rnk
    integer, allocatable    ::  idmap(:,:)  ! (*cart_dims)
    character(len=16)       ::  fmt
    logical                 ::  has_west, has_east, has_south, has_north
    
    ! --- begin ----------------------------------
    
    ! problem size:
    self%ndim = size(cart_dims)

    ! storage for cartesian topology:
    allocate( self%cart_dims(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%cart_periods(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%cart_coords(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%cart_src_id(self%ndim,ndir), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%cart_dst_id(self%ndim,ndir), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy dims:
    self%cart_dims = cart_dims
    ! no periodicity:
    self%cart_periods = .false.

    ! check ...
    if ( product(cart_dims) /= goc%npes ) then
      write (gol,'("requested decomposition on ",i0," domains while number of processors is ",i0," ...")') &
                         product(cart_dims), goc%npes; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! define cartesian processor topology,
    ! allow re-ordering of the processors ;
    ! return dedicated communicator:
    call MPI_Cart_Create( goc%comm, self%ndim, self%cart_dims, self%cart_periods, .true., &
                           self%cart_comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
  
    ! get properties (dims, periods, and coords):
    call MPI_Cart_Get( self%cart_comm, self%ndim, &
                        self%cart_dims, self%cart_periods, self%cart_coords, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! dummy:
    self%cart_coords = 0
#endif

    ! info ...
    write (gol,'("cartesian topology:")'); call goPr
    write (gol,'("  ndim          = ",i0)') self%ndim; call goPr
    if ( self%ndim == 2 ) then
      write (gol,'("  dims          : ",i0,",",i0)') self%cart_dims; call goPr
      write (gol,'("  periodic      : ",l1,",",l1)') self%cart_periods; call goPr
      write (gol,'("  local coords  : ",i0,",",i0)') self%cart_coords; call goPr
      ! storage:
      allocate( idmap(0:self%cart_dims(1)-1,0:self%cart_dims(2)-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! list processor id's:
      write (gol,'("  processor ids:")'); call goPr
      ! loop over domains (zero based indices!):
      do j = 0, self%cart_dims(2)-1
        do i = 0, self%cart_dims(1)-1
#ifdef _MPI
          ! get rank within cartesian grid:
          call MPI_Cart_Rank( self%cart_comm, (/i,j/), rnk, status )
          IF_MPI_NOT_OK_RETURN(status=1)
#else
          ! serial, all on root:
          rnk = 0
#endif
          ! info ...
          write (gol,'("    domain (",i0,",",i0,") on pe ",i0)') i, j, rnk; call goPr
          ! store:
          idmap(i,j) = rnk
        end do ! i
      end do ! j
      ! show "map" of processor id's
      write (gol,'("  processor map:")'); call goPr
      write (fmt,'("(4x,",i0,"i3)")') self%cart_dims(1)
      do j = self%cart_dims(2)-1, 0, -1
        write (gol,fmt) idmap(:,j); call goPr
      end do ! j
      ! clear:
      deallocate( idmap, stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      write (gol,'("unsupported ndim ",i0)') self%ndim; call goPr
      TRACEBACK; status=1; return
    end if
  
    ! local domain at border?
    call self%Get( status, has_west=has_west, has_east=has_east, has_south=has_south, has_north=has_north )
    IF_NOT_OK_RETURN(status=1)
    ! info ...
    write (gol,'("  global domain borders:")'); call goPr
    write (gol,'("    west    : ",l1)') has_west; call goPr
    write (gol,'("    east    : ",l1)') has_east; call goPr
    write (gol,'("    south   : ",l1)') has_south; call goPr
    write (gol,'("    north   : ",l1)') has_north; call goPr

    ! info ..
    write (gol,'("communication with neighbours:")'); call goPr
    ! processor id's for communication over neighbours:
    do idim = 1, self%ndim
      ! info ..
      write (gol,'("dimension: ",i0)') idim; call goPr
      ! loop over directions:
      do idir = 1, ndir
#ifdef _MPI
        ! obtain souce and dest id for cummunication in this dimension;
        ! direction (dimension) is zero-based:
        call MPI_Cart_Shift( self%cart_comm, idim-1, dir_shift(idir), &
                              self%cart_src_id(idim,idir), self%cart_dst_id(idim,idir), status )
        IF_MPI_NOT_OK_RETURN(status=1)
#else
        ! no neighbours:
        self%cart_src_id(idim,idir) = -2
        self%cart_dst_id(idim,idir) = -2
#endif
        ! info ..
        write (gol,'("  direction: ",i0)') dir_shift(idir); call goPr
        write (gol,'("    receive from : ",i0)') self%cart_src_id(idim,idir); call goPr
        write (gol,'("    send to      : ",i0)') self%cart_dst_id(idim,idir); call goPr
      end do ! dirs
    end do ! dims

    ! *
    
    ! check ...
    if ( size(glb_lbo) /= self%ndim ) then
      write (gol,'("size of argument glb_lbo (",i0,") does not match with size of topology (",i0,")")') &
                      size(glb_lbo), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(glb_ubo) /= self%ndim ) then
      write (gol,'("size of argument glb_ubo (",i0,") does not match with size of topology (",i0,")")') &
                      size(glb_ubo), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage for global bounds and shape:
    allocate( self%glb_lbo(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%glb_ubo(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%glb_shp(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%glb_lbo = glb_lbo
    self%glb_ubo = glb_ubo
    ! fill shape:
    self%glb_shp = glb_ubo - glb_lbo + 1
    
    ! storage for global bounds on this pe:
    allocate( glbo(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( gubo(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, self%ndim
    
      ! storage for number of local cells per domain in this dimension:
      allocate( nn(0:self%cart_dims(idim)-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! global number of cells:
      gn = glb_ubo(idim) - glb_lbo(idim) + 1
      ! obtain local size in this dimension:
      call AssignLocalSizes( gn, cart_dims(idim), nn, status )
      IF_NOT_OK_RETURN(status=1)
      ! 0-based index within this dimension for local domain:
      idom = self%cart_coords(idim)
      ! global index space:
      gubo(idim) = glb_lbo(idim) - 1
      do i = 0, idom
        gubo(idim) = gubo(idim) + nn(i)
      end do
      glbo(idim) = gubo(idim) - nn(idom) + 1
      ! clear:
      deallocate( nn, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end do ! dims
    
    ! storage:
    allocate( self%glbo(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%gubo(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%off(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%lbo(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ubo(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%shp(self%ndim,0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%n(0:goc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
#ifdef _MPI
    ! exchange global bounds:
        
    call MPI_AllGather(      glbo, self%ndim, MPI_INTEGER, &
                        self%glbo, self%ndim, MPI_INTEGER, &
                        goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    call MPI_AllGather(      gubo, self%ndim, MPI_INTEGER, &
                        self%gubo, self%ndim, MPI_INTEGER, &
                        goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    self%glbo(:,0) = glbo
    self%gubo(:,0) = gubo
#endif
    
    ! local offset in global space:
    self%off = self%glbo - 1
    
    ! local shapes:
    self%shp = self%gubo - self%glbo + 1
    ! trap undefined per process:
    do i = 0, goc%npes-1
      ! any dimension undefined ? set all to zero:
      if ( any(self%shp(:,i) <= 0) ) self%shp(:,i) = 0
    end do
    
    ! total number of local elements,
    ! will be zero for processes with empty domain:
    self%n = product( self%shp, dim=1 )
    
    ! local bounds:
    self%lbo = 1
    self%ubo = self%shp
    
    ! clear:
    deallocate( glbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( gubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("domain definition:")'); call goPr
    if ( self%ndim == 1 ) then
      write (gol,'("  global: (",i0,":",i0,")")') self%glbo(1,goc%id), self%gubo(1,goc%id); call goPr
      write (gol,'("  local : (",i0,":",i0,")")') self%lbo(1,goc%id), self%ubo(1,goc%id); call goPr
      write (gol,'("  offset: ",i0)') self%off(:,goc%id); call goPr
      write (gol,'("  shape : ",i0)') self%shp(:,goc%id); call goPr
      write (gol,'("  size  : ",i0)') self%n(goc%id); call goPr
    else if ( self%ndim == 2 ) then
      write (gol,'("  global: (",i0,":",i0,",",i0,":",i0,")")') self%glbo(1,goc%id), self%gubo(1,goc%id), &
                                                                self%glbo(2,goc%id), self%gubo(2,goc%id); call goPr
      write (gol,'("  local : (",i0,":",i0,",",i0,":",i0,")")') self%lbo(1,goc%id), self%ubo(1,goc%id), &
                                                                self%lbo(2,goc%id), self%ubo(2,goc%id); call goPr
      write (gol,'("  offset: ",i0,",",i0)') self%off(:,goc%id); call goPr
      write (gol,'("  shape : ",i0,",",i0)') self%shp(:,goc%id); call goPr
      write (gol,'("  size  : ",i0)') self%n(goc%id); call goPr
    else
      write (gol,'("unsupported number of dimensions ",i0)') self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    !! test exchange of halo cells ...
    !call self%ExchangeTest( status )
    !IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Domains_Init
  

  ! *

  
  !
  ! Distribute 'gn' values over 'ndom' domains,
  ! return number of local cells 'nn(:)' for each domain
  !
  
  subroutine AssignLocalSizes( gn, ndom, nn, status )

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/AssignLocalSizes'
    
    ! --- in/out ---------------------------------
    
    integer, intent(in)               ::  gn
    integer, intent(in)               ::  ndom
    integer, intent(out)              ::  nn(ndom)
    integer, intent(out)              ::  status

    ! --- local ----------------------------------
    
    integer     ::  d
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! first estimate of local size: each domain same number;
    ! integer division will round to lower integer:
    nn = gn / ndom
    ! number of missing cells:
    d = gn - sum(nn)
    ! any missing?
    if ( d > 0 ) then
      ! add from top to bottom, first domain is usually more busy already
      do i = ndom, ndom-d+1, -1
        nn(i) = nn(i) + 1
      end do
    end if
    
    ! ok
    status = 0
    
  end subroutine AssignLocalSizes


  ! ***
  
  
  subroutine Domains_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(inout)      ::  self
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! storage for global bounds and shape:
    deallocate( self%glb_lbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%glb_ubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%glb_shp, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%glbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%gubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%lbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%ubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%off, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%shp, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%n, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
!    ! done with topology:
!    call self%topo%Done( status )
!    IF_NOT_OK_RETURN(status=1)

    ! clear cartesian topology:
    deallocate( self%cart_dims, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%cart_periods, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%cart_coords, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%cart_src_id, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%cart_dst_id, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Domains_Done


  ! ***
  
  
  subroutine Domains_Get( self, status, shp, off, glb_shp, glbo, gubo, &
                                 has_west, has_east, has_south, has_north )

    use GO_Comm, only : goc
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)          ::  self
    integer, intent(out)                  ::  status

    integer, intent(out), optional        ::  shp(:)
    integer, intent(out), optional        ::  off(:)
    integer, intent(out), optional        ::  glb_shp(:)
    integer, intent(out), optional        ::  glbo(:), gubo(:)
    logical, intent(out), optional        ::  has_west, has_east, has_south, has_north

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! return local shape?
    if ( present(shp) ) then
      ! check ..
      if ( size(shp) /= self%ndim ) then
        write (gol,'("argument shp has size ",i0," while dim is ",i0)') size(shp), self%ndim; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      shp = self%shp(:,goc%id)
    end if
    
    ! return local shape?
    if ( present(off) ) then
      ! check ..
      if ( size(off) /= self%ndim ) then
        write (gol,'("argument off has size ",i0," while dim is ",i0)') size(off), self%ndim; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      off = self%off(:,goc%id)
    end if
    
    ! return global lower bounds?
    if ( present(glbo) ) then
      ! check ..
      if ( size(glbo) /= self%ndim ) then
        write (gol,'("argument glbo has size ",i0," while dim is ",i0)') size(off), self%ndim; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      glbo = self%glbo(:,goc%id)
    end if
    
    ! return global upper bounds?
    if ( present(gubo) ) then
      ! check ..
      if ( size(gubo) /= self%ndim ) then
        write (gol,'("argument gubo has size ",i0," while dim is ",i0)') size(off), self%ndim; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      gubo = self%gubo(:,goc%id)
    end if
    
    ! return global shape?
    if ( present(glb_shp) ) then
      ! check ..
      if ( size(glb_shp) /= self%ndim ) then
        write (gol,'("argument glb_shp has size ",i0," while dim is ",i0)') size(glb_shp), self%ndim; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      glb_shp = self%glb_shp
    end if
    
    ! has west border ?
    if ( present(has_west) ) then
      has_west = self%cart_coords(1) == 0
    end if
    ! has east border ?
    if ( present(has_east) ) then
      has_east = self%cart_coords(1) == self%cart_dims(1)-1
    end if
    ! has south border ?
    if ( present(has_south) ) then
      if ( self%ndim >= 2 ) then
        has_south = self%cart_coords(2) == 0
      else
        has_south = .true.
      end if
    end if
    ! has north border ?
    if ( present(has_north) ) then
      if ( self%ndim >= 2 ) then
        has_north = self%cart_coords(2) == self%cart_dims(2)-1
      else
        has_north = .true.
      end if
    end if
    
    ! ok
    status = 0

  end subroutine Domains_Get


  ! ***
  
  
  subroutine Domains_GetDim( self, idim, status, cart_size, cart_coord )

    use GO_Comm, only : goc
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)          ::  self
    integer, intent(in)                   ::  idim
    integer, intent(out)                  ::  status

    integer, intent(out), optional        ::  cart_size
    integer, intent(out), optional        ::  cart_coord

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_GetDim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (idim < 1) .or. (idim > self%ndim) ) then
      write (gol,'("dimension ",i0," outside range 1,..,",i0)') idim, self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return size of cartesian domain layout in this dimension?
    if ( present(cart_size) ) cart_size = self%cart_dims(idim)
    ! return coordinate in cartesian domain layout?
    if ( present(cart_coord) ) cart_coord = self%cart_coords(idim)
    
    ! ok
    status = 0

  end subroutine Domains_GetDim
        

  ! ***
  
  !
  ! Compute local intersection bounds and shape.
  ! Return status:
  !   -1  : no intersection
  !    0  : ok
  !   else error
  !

  subroutine Domains_Intersection( self, glbo, gubo, lbo, ubo, shp, status )
  
    use GO_Comm, only : goc
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)         ::  self
    integer, intent(in)                  ::  glbo(:)
    integer, intent(in)                  ::  gubo(:)
    integer, intent(out)                 ::  lbo(:)
    integer, intent(out)                 ::  ubo(:)
    integer, intent(out)                 ::  shp(:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( any( (/ size(glbo), size(gubo), size(lbo), size(ubo), size(shp) /) /= self%ndim ) ) then
      write (gol,'("arguments should have size ndim = ",i0)') self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check on empty domains first:
    if ( any( self%shp(:,goc%id) <= 0 ) ) then

      ! empty local domain:
      status = -1

    else if ( any( gubo < glbo ) ) then

      ! empty inquired domain:
      status = -1

    else

      ! local bounds of intersection:
      lbo = max( self%glbo(:,goc%id), glbo ) - self%off(:,goc%id)
      ubo = min( self%gubo(:,goc%id), gubo ) - self%off(:,goc%id)
      ! local shape:
      shp = ubo - lbo + 1

      ! set return status:
      if ( any(shp <= 0) ) then
        status = -1
      else
        status = 0
      end if
      
    end if
    
  end subroutine Domains_Intersection


  ! ***
  
  
  ! return status:
  !   -1  : no
  !    0  : yes
  !   else error
  
  subroutine Domains_Inside( self, ind, status )
  
    use GO_Comm, only : goc
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)            ::  self
    integer, intent(in)                     ::  ind(:)  ! global indices
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Inside'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(ind) /= self%ndim ) then
      write (gol,'("argument ind has size ",i0," while ndim = ",i0)') size(ind), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! compare:
    if ( all(self%glbo(:,goc%id) <= ind) .and. all(ind <= self%gubo(:,goc%id)) ) then
      ! in domain:
      status = 0
    else
      ! outside ..
      status = -1
    end if
    
  end subroutine Domains_Inside
  
  
  ! index of domain holding cell with provided global indices
  
  subroutine Domains_Find( self, ind, iproc, status, locind )
  
    use GO_Comm, only : goc
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)            ::  self
    integer, intent(in)                     ::  ind(:)  ! global indices
    integer, intent(out)                    ::  iproc   ! 0:nproc-1
    integer, intent(out)                    ::  status

    integer, intent(out), optional          ::  locind(:)  ! local indices

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Find'
    
    ! --- local ----------------------------------
    
    integer     ::  k
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(ind) /= self%ndim ) then
      write (gol,'("argument ind has size ",i0," while ndim = ",i0)') size(ind), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! dummy result:
    iproc = -999
    ! loop:
    do k = 0, goc%npes-1
      ! compare:
      if ( all(self%glbo(:,k) <= ind) .and. all(ind <= self%gubo(:,k)) ) then
        ! found !
        iproc = k
        ! set output?
        if ( present(locind) ) then
          ! local indices are global minus offset:
          locind = ind - self%off(:,k)
        end if
        ! leave:
        exit
      end if  ! location on proc domain
    end do  ! procs
    ! check ...
    if ( iproc < 0 ) then
      write (gol,*) 'could not find domain holding global indices : ', ind; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine Domains_Find


  ! ***************************************************************************
  ! ***
  ! *** exchange halo cells
  ! ***
  ! ***************************************************************************
  
  
  !
  ! Exchange halo cells between domains.
  ! Input/output field should be passed as pointer to maintain the index space;
  ! index space includes halo:
  !   (1-nh:nx+nh,1-nh:ny+nh)
  ! The field is passed as pointer to maintain original index space.
  ! Shape (nx,ny) of core area (1:nx,1:ny) is taken from the index spaces 
  ! stored in the attributes.
  !
  
  subroutine Domains_ExchangeHalo_3d_r( self, nh, field, status )

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHalo_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    real, pointer                 ::  field(:,:,:)   ! (1-nh:nx+nh,1-nh:ny+nh,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! apply in two directions:
    call self%ExchangeHaloX( nh, field, status )
    IF_NOT_OK_RETURN(status=1)
    call self%ExchangeHaloY( nh, field, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHalo_3d_r
  
  ! *
  
  subroutine Domains_ExchangeHaloX_3d_r( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloX_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    real, pointer                 ::  field(:,:,:)   ! (1-nh:nx+nh,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  x1, x2
    integer                     ::  ny, nz
    real, allocatable           ::  bound_e(:,:,:)  ! (nh,:,:)
    real, allocatable           ::  bound_w(:,:,:)  ! (nh,:,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    x1 = self%lbo(1,goc%id)
    x2 = self%ubo(1,goc%id)
    ! other dimensions:
    ny = size(field,2)
    nz = size(field,3)

    ! check ...
    if ( (x1-nh < lbound(field,1)) .or. (x2+nh > ubound(field,1)) ) then
      write (gol,'("index space of field (",i0,":",i0,")")') lbound(field,1), ubound(field,1); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(1,goc%id), self%ubo(1,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(1,goc%id)-nh, self%ubo(1,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_w(nh,ny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_e(nh,ny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ west exchange:
    !                  x1   x1+1      cell
    !         + - - +-----+-----+---
    !                  v     v
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2   x2+1 x2+2       cell
    !
    idim = 1  ! west-east
    idir = 1  ! previous (=west)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_w = field(x1:x1+nh-1,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send west, receive east:
      call MPI_SendRecv( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send west:
      call MPI_Send( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! receive east:
      call MPI_Recv( bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x2+1:x2+nh,:,:) = bound_e
    end if
    !
    !~~ east exchange:
    !     x1-2  x1-1   x1   x1+1      cell
    !   + - - + - - +-----+-----+---
    !      ^     ^
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2                   cell
    !
    idim = 1  ! west-east
    idir = 2  ! next (=east)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_e = field(x2-nh+1:x2,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send east, receive west:
      call MPI_SendRecv( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send east:
      call MPI_Send( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ receive west:
      call MPI_Recv( bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x1-nh:x1-1,:,:) = bound_w
    end if

    ! clear:
    deallocate( bound_w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_e, stat=status )
    IF_NOT_OK_RETURN(status=1)

#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloX_3d_r
  
  ! *
    
  subroutine Domains_ExchangeHaloX_4d_r( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloX_4d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    real, pointer                 ::  field(:,:,:,:)   ! (1-nh:nx+nh,:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  x1, x2
    integer                     ::  ny, nz, nl
    real, allocatable           ::  bound_e(:,:,:,:)  ! (nh,:,:)
    real, allocatable           ::  bound_w(:,:,:,:)  ! (nh,:,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    x1 = self%lbo(1,goc%id)
    x2 = self%ubo(1,goc%id)
    ! other dimensions:
    ny = size(field,2)
    nz = size(field,3)
    nl = size(field,4)

    ! check ...
    if ( (x1-nh < lbound(field,1)) .or. (x2+nh > ubound(field,1)) ) then
      write (gol,'("index space of field (",i0,":",i0,")")') lbound(field,1), ubound(field,1); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(1,goc%id), self%ubo(1,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(1,goc%id)-nh, self%ubo(1,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_w(nh,ny,nz,nl), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_e(nh,ny,nz,nl), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ west exchange:
    !                  x1   x1+1      cell
    !         + - - +-----+-----+---
    !                  v     v
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2   x2+1 x2+2       cell
    !
    idim = 1  ! west-east
    idir = 1  ! previous (=west)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_w = field(x1:x1+nh-1,:,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send west, receive east:
      call MPI_SendRecv( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send west:
      call MPI_Send( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! receive east:
      call MPI_Recv( bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x2+1:x2+nh,:,:,:) = bound_e
    end if
    !
    !~~ east exchange:
    !     x1-2  x1-1   x1   x1+1      cell
    !   + - - + - - +-----+-----+---
    !      ^     ^
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2                   cell
    !
    idim = 1  ! west-east
    idir = 2  ! next (=east)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_e = field(x2-nh+1:x2,:,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send east, receive west:
      call MPI_SendRecv( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send east:
      call MPI_Send( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ receive west:
      call MPI_Recv( bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x1-nh:x1-1,:,:,:) = bound_w
    end if

    ! clear:
    deallocate( bound_w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_e, stat=status )
    IF_NOT_OK_RETURN(status=1)

#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloX_4d_r
  
  ! *

  subroutine Domains_ExchangeHaloY_3d_r( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloY_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    real, pointer                 ::  field(:,:,:)   ! (:,1-nh:ny+nh,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  y1, y2
    integer                     ::  nx, nz
    real, allocatable           ::  bound_n(:,:,:)  ! (1-nh:nx+nh,nh,:)
    real, allocatable           ::  bound_s(:,:,:)  ! (1-nh:nx+nh,nh,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    y1 = self%lbo(2,goc%id)
    y2 = self%ubo(2,goc%id)
    ! other dimensions:
    nx = size(field,1)
    nz = size(field,3)

    ! check ...
    if ( (y1-nh < lbound(field,2)) .or. (y2+nh > ubound(field,2))        ) then
      write (gol,'("y index space of field (",i0,":",i0,")")') lbound(field,2), ubound(field,2); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(2,goc%id)-nh, self%ubo(2,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_s(nx,nh,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_n(nx,nh,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ south exchange:
    !             
    !     cell  |
    !           +    +
    !     y1+1  | -> : y2+2
    !           +    + 
    !     y1    | -> : y2+1 
    !           +    + 
    !           :    | y2
    !           +    + 
    !                | 
    !                  cell
    !                  
    !
    idim = 2  ! south-north
    idir = 1  ! previous (=south)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_s = field(:,y1:y1+nh-1,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send south, receive north:
      call MPI_SendRecv( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send south:
      call MPI_Send( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send south, receive north:
      call MPI_Recv( bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y2+1:y2+nh,:) = bound_n
    end if
    !
    !~~ north exchange:
    !            
    !     cell |
    !          +    +
    !     y1+1 |    :
    !          +    + 
    !     y1   |    : 
    !          +    + 
    !          : <- | y2
    !          +    + 
    !          : <- | y2-1 
    !          +    +  
    !                 cell
    !
    idim = 2  ! south-north
    idir = 2  ! next (=north)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_n = field(:,y2-nh+1:y2,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send north, receive south:
      call MPI_SendRecv( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send north:
      call MPI_Send( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ receive south:
      call MPI_Recv( bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y1-nh:y1-1,:) = bound_s
    end if

    ! clear:
    deallocate( bound_s, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_n, stat=status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloY_3d_r
  
  ! *

  subroutine Domains_ExchangeHaloY_4d_r( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloY_4d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    real, pointer                 ::  field(:,:,:,:) ! (:,1-nh:ny+nh,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  y1, y2
    integer                     ::  nx, nz, nl
    real, allocatable           ::  bound_n(:,:,:,:)  ! (1-nh:nx+nh,nh,:)
    real, allocatable           ::  bound_s(:,:,:,:)  ! (1-nh:nx+nh,nh,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    y1 = self%lbo(2,goc%id)
    y2 = self%ubo(2,goc%id)
    ! other dimensions:
    nx = size(field,1)
    nz = size(field,3)
    nl = size(field,4)
    ! check ...
    if ( (y1-nh < lbound(field,2)) .or. (y2+nh > ubound(field,2))        ) then
      write (gol,'("y index space of field (",i0,":",i0,")")') lbound(field,2), ubound(field,2); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(2,goc%id)-nh, self%ubo(2,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_s(nx,nh,nz,nl), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_n(nx,nh,nz,nl), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ south exchange:
    !             
    !     cell  |
    !           +    +
    !     y1+1  | -> : y2+2
    !           +    + 
    !     y1    | -> : y2+1 
    !           +    + 
    !           :    | y2
    !           +    + 
    !                | 
    !                  cell
    !                  
    !
    idim = 2  ! south-north
    idir = 1  ! previous (=south)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_s = field(:,y1:y1+nh-1,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send south, receive north:
      call MPI_SendRecv( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send south:
      call MPI_Send( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send south, receive north:
      call MPI_Recv( bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y2+1:y2+nh,:,:) = bound_n
    end if
    !
    !~~ north exchange:
    !            
    !     cell |
    !          +    +
    !     y1+1 |    :
    !          +    + 
    !     y1   |    : 
    !          +    + 
    !          : <- | y2
    !          +    + 
    !          : <- | y2-1 
    !          +    +  
    !                 cell
    !
    idim = 2  ! south-north
    idir = 2  ! next (=north)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_n = field(:,y2-nh+1:y2,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send north, receive south:
      call MPI_SendRecv( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send north:
      call MPI_Send( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ receive south:
      call MPI_Recv( bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y1-nh:y1-1,:,:) = bound_s
    end if

    ! clear:
    deallocate( bound_s, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_n, stat=status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloY_4d_r
  
  ! *
  
  subroutine Domains_ExchangeHaloX_3d_l( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloX_3d_l'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    logical, pointer              ::  field(:,:,:)   ! (1-nh:nx+nh,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  x1, x2
    integer                     ::  ny, nz
    logical, allocatable        ::  bound_e(:,:,:)  ! (nh,:,:)
    logical, allocatable        ::  bound_w(:,:,:)  ! (nh,:,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    x1 = self%lbo(1,goc%id)
    x2 = self%ubo(1,goc%id)
    ! other dimensions:
    ny = size(field,2)
    nz = size(field,3)

    ! check ...
    if ( (x1-nh < lbound(field,1)) .or. (x2+nh > ubound(field,1)) ) then
      write (gol,'("index space of field (",i0,":",i0,")")') lbound(field,1), ubound(field,1); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(1,goc%id), self%ubo(1,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(1,goc%id)-nh, self%ubo(1,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_w(nh,ny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_e(nh,ny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'logical', 0, mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ west exchange:
    !                  x1   x1+1      cell
    !         + - - +-----+-----+---
    !                  v     v
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2   x2+1 x2+2       cell
    !
    idim = 1  ! west-east
    idir = 1  ! previous (=west)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_w = field(x1:x1+nh-1,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send west, receive east:
      call MPI_SendRecv( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send west:
      call MPI_Send( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! receive east:
      call MPI_Recv( bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x2+1:x2+nh,:,:) = bound_e
    end if
    !
    !~~ east exchange:
    !     x1-2  x1-1   x1   x1+1      cell
    !   + - - + - - +-----+-----+---
    !      ^     ^
    !  -+-----+-----+ - - + - - +
    !     x2-1   x2                   cell
    !
    idim = 1  ! west-east
    idir = 2  ! next (=east)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_e = field(x2-nh+1:x2,:,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send east, receive west:
      call MPI_SendRecv( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send east:
      call MPI_Send( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ receive west:
      call MPI_Recv( bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(x1-nh:x1-1,:,:) = bound_w
    end if

    ! clear:
    deallocate( bound_w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_e, stat=status )
    IF_NOT_OK_RETURN(status=1)

#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloX_3d_l
  
  ! *
  
  subroutine Domains_ExchangeHaloY_3d_l( self, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloY_3d_l'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  nh             ! number of halo cells
    logical, pointer              ::  field(:,:,:)   ! (:,1-nh:ny+nh,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  y1, y2
    integer                     ::  nx, nz
    logical, allocatable        ::  bound_n(:,:,:)  ! (1-nh:nx+nh,nh,:)
    logical, allocatable        ::  bound_s(:,:,:)  ! (1-nh:nx+nh,nh,:)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------

#ifdef _MPI
    
    ! short variables for local index spaces:
    y1 = self%lbo(2,goc%id)
    y2 = self%ubo(2,goc%id)
    ! other dimensions:
    nx = size(field,1)
    nz = size(field,3)

    ! check ...
    if ( (y1-nh < lbound(field,2)) .or. (y2+nh > ubound(field,2))        ) then
      write (gol,'("y index space of field (",i0,":",i0,")")') lbound(field,2), ubound(field,2); call goErr
      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(2,goc%id)-nh, self%ubo(2,goc%id)+nh; call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bound_s(nx,nh,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( bound_n(nx,nh,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! data type:
    call goc%GetDataType( 'logical', 0, mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! init tag number (not necessary?)
    tag = 0

    !~~ south exchange:
    !             
    !     cell  |
    !           +    +
    !     y1+1  | -> : y2+2
    !           +    + 
    !     y1    | -> : y2+1 
    !           +    + 
    !           :    | y2
    !           +    + 
    !                | 
    !                  cell
    !                  
    !
    idim = 2  ! south-north
    idir = 1  ! previous (=south)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_s = field(:,y1:y1+nh-1,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! send south, receive north:
      call MPI_SendRecv( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      ! send south:
      call MPI_Send( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      ! receive north:
      call MPI_Recv( bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y2+1:y2+nh,:) = bound_n
    end if
    !
    !~~ north exchange:
    !            
    !     cell |
    !          +    +
    !     y1+1 |    :
    !          +    + 
    !     y1   |    : 
    !          +    + 
    !          : <- | y2
    !          +    + 
    !          : <- | y2-1 
    !          +    +  
    !                 cell
    !
    idim = 2  ! south-north
    idir = 2  ! next (=north)
    ! increase tag:
    tag = tag + 1
    ! fill send buffer from local field:
    bound_n = field(:,y2-nh+1:y2,:)
    ! on some machines send/recv does not work for undefined dest/source ...
    if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send north, receive south:
      call MPI_SendRecv( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                         bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                         self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) >= 0) .and. (self%cart_src_id(idim,idir) < 0) ) then

      !~ send north:
      call MPI_Send( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
                     self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (self%cart_dst_id(idim,idir) < 0) .and. (self%cart_src_id(idim,idir) >= 0) ) then

      !~ send north:
      call MPI_Recv( bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
                     self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    end if
    ! any received?
    if ( self%cart_src_id(idim,idir) >= 0 ) then
      ! store in halo:
      field(:,y1-nh:y1-1,:) = bound_s
    end if

    ! clear:
    deallocate( bound_s, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( bound_n, stat=status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0
    
  end subroutine Domains_ExchangeHaloY_3d_l
  
  ! *
  
!  subroutine Domains_ExchangeHaloEdgeX_3d_r( self, nh, field, status )
!  
!    use GO_Comm, only : goc
!#ifdef _MPI
!    use MPI_F08, only : MPI_DataType
!    use MPI_F08, only : MPI_Status
!    use MPI_F08, only : MPI_SendRecv
!#endif
!
!    ! --- const ---------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloEdgeX_3d_r'
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_Domains), intent(in)  ::  self
!    integer, intent(in)           ::  nh             ! number of halo cells
!    real, pointer                 ::  field(:,:,:)   ! (0-nh:nx+nh,1:ny,:)
!    integer, intent(out)          ::  status
!    
!    ! --- local -----------------------------------
!  
!    integer                     ::  x1, x2, y1, y2
!    integer                     ::  nz
!    real, allocatable           ::  bound_e(:,:,:)  ! (nh,1:ny,:)
!    real, allocatable           ::  bound_w(:,:,:)  ! (nh,1:ny,:)
!    integer                     ::  idim
!    integer                     ::  idir
!    integer                     ::  tag
!#ifdef _MPI
!    type(MPI_DataType)          ::  mpi_dtype
!    type(MPI_Status)            ::  mpi_stat
!#endif
!    
!        integer     ::  i
!    
!    ! --- begin -----------------------------------
!
!#ifdef _MPI
!    
!    ! short variables for local index spaces:
!    x1 = self%lbo(1,goc%id)
!    x2 = self%ubo(1,goc%id)
!    y1 = self%lbo(2,goc%id)
!    y2 = self%ubo(2,goc%id)
!    ! other dimensions:
!    nz = size(field,3)
!
!    ! check ...
!    if ( (x1-1-nh < lbound(field,1)) .or. (x2+nh > ubound(field,1)) .or. &
!         (y1      < lbound(field,2)) .or. (y2    > ubound(field,2))        ) then
!      write (gol,'("index space of field (",i0,":",i0,",",i0,":",i0,")")') &
!              lbound(field,1), ubound(field,1), lbound(field,2), ubound(field,2); call goErr
!      write (gol,'("does not cover expected index space (",i0,":",i0,",",i0,":",i0,")")') &
!              self%lbo(1,goc%id), self%ubo(1,goc%id), self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
!      write (gol,'("with ",i0," halo (",i0,":",i0,",",i0,":",i0,")")') nh, &
!              self%lbo(1,goc%id)-1-nh, self%ubo(1,goc%id)+nh, self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! storage:
!    allocate( bound_w(nh,y1:y2,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( bound_e(nh,y1:y2,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! data type:
!    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! init tag number (not necessary?)
!    tag = 0
!
!    !~~ west exchange:
!    !                  x1    x1+1      cell
!    !       x1-2   x1-1   x1   x1+1     edge
!    !         + - - +-----+-----+---
!    !                     v
!    !  -+-----+-----+ - - +
!    !       x2-1   x2   x2+1           edge
!    !     x2-1   x2                   cell
!    !
!    idim = 1  ! west-east
!    idir = 1  ! previous (=west)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_w = field(x1:x1+nh-1,y1:y2,:)
!    ! send west, receive east:
!    call MPI_SendRecv( bound_w, size(bound_w), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_e, size(bound_e), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(x2+1:x2+nh,y1:y2,:) = bound_e
!    end if
!    !
!    !~~ east exchange:
!    !                  x1    x1+1      cell
!    !       x1-2   x1-1   x1   x1+1     edge
!    !         + - - +-----+-----+---
!    !         ^
!    !  -+-----+-----+ - - +
!    !       x2-1   x2   x2+1            edge
!    !     x2-1   x2                   cell
!    !
!    idim = 1  ! west-east
!    idir = 2  ! next (=east)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_e = field(x2-nh:x2-1,y1:y2,:)
!!    write (gol,*) 'rrr '; call goPr
!!    write (gol,*) 'rrr fill buffer ...'; call goPr
!!    do i = lbound(field,1), ubound(field,1)
!!      write (gol,*) 'rrr field(',self%off(1,goc%id)+i,') = ', field(i,10,1); call goPr
!!    end do
!!    write (gol,*) 'rrr bound_e = ', bound_e(:,10,1); call goPr
!!    write (gol,*) 'rrr send/recv ...'; call goPr
!    !~ send east, receive west:
!    call MPI_SendRecv( bound_e, size(bound_e), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_w, size(bound_w), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(x1-1-nh:x1-2,y1:y2,:) = bound_w
!!      write (gol,*) 'rrr bound_w = ', bound_w(:,10,1); call goPr
!!      do i = lbound(field,1), ubound(field,1)
!!        write (gol,*) 'rrr field(',self%off(1,goc%id)+i,') = ', field(i,10,1); call goPr
!!      end do
!    end if
!!    write (gol,*) 'rrr '; call goPr
!
!    ! storage:
!    deallocate( bound_w, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( bound_e, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!#endif
!
!    ! ok
!    status = 0
!    
!  end subroutine Domains_ExchangeHaloEdgeX_3d_r
!  
!  ! *
!  
!  subroutine Domains_ExchangeHaloEdgeY_3d_r( self, nh, field, status )
!  
!    use GO_Comm, only : goc
!#ifdef _MPI
!    use MPI_F08, only : MPI_DataType
!    use MPI_F08, only : MPI_Status
!    use MPI_F08, only : MPI_SendRecv
!#endif
!
!    ! --- const ---------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloEdgeY_3d_r'
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_Domains), intent(in)  ::  self
!    integer, intent(in)           ::  nh             ! number of halo cells
!    real, pointer                 ::  field(:,:,:)   ! (1:nx,0-nh:ny+nh,:)
!    integer, intent(out)          ::  status
!    
!    ! --- local -----------------------------------
!  
!    integer                     ::  x1, x2, y1, y2
!    integer                     ::  nz
!    real, allocatable           ::  bound_n(:,:,:)  ! (1:nx,nh,:)
!    real, allocatable           ::  bound_s(:,:,:)  ! (1:nx,nh,:)
!    integer                     ::  idim
!    integer                     ::  idir
!    integer                     ::  tag
!#ifdef _MPI
!    type(MPI_DataType)          ::  mpi_dtype
!    type(MPI_Status)            ::  mpi_stat
!#endif
!    
!    ! --- begin -----------------------------------
!
!#ifdef _MPI
!    
!    ! short variables for local index spaces:
!    x1 = self%lbo(1,goc%id)
!    x2 = self%ubo(1,goc%id)
!    y1 = self%lbo(2,goc%id)
!    y2 = self%ubo(2,goc%id)
!    ! other dimensions:
!    nz = size(field,3)
!
!    ! check ...
!    if ( (x1      < lbound(field,1)) .or. (x2    > ubound(field,1)) .or. &
!         (y1-1-nh < lbound(field,2)) .or. (y2+nh > ubound(field,2))        ) then
!      write (gol,'("index space of field (",i0,":",i0,",",i0,":",i0,")")') &
!              lbound(field,1), ubound(field,1), lbound(field,2), ubound(field,2); call goErr
!      write (gol,'("does not cover expected index space (",i0,":",i0,",",i0,":",i0,")")') &
!              self%lbo(1,goc%id), self%ubo(1,goc%id), self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
!      write (gol,'("with ",i0," halo (",i0,":",i0,",",i0,":",i0,")")') nh, &
!              self%lbo(1,goc%id), self%ubo(1,goc%id), self%lbo(2,goc%id)-1-nh, self%ubo(2,goc%id)+nh; call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! storage:
!    allocate( bound_s(x1:x2,nh,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( bound_n(x1:x2,nh,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! data type:
!    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! init tag number (not necessary?)
!    tag = 0
!
!    !~~ south exchange:
!    !          edge    
!    !     cell       |
!    !          y1+1  +
!    !     y1+1       |
!    !          y1    + -> + (y2+1)
!    !     y1         |    :      
!    !          y1-1  +    + y2
!    !                :    |      y2
!    !         (y1-2) +    + y2-1
!    !                     |      
!    !                            cell
!    !                       edge
!    !
!    idim = 2  ! south-north
!    idir = 1  ! previous (=south)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_s = field(x1:x2,y1:y1+nh-1,:)
!    ! send south, receive north:
!    call MPI_SendRecv( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(x1:x2,y2+1:y2+nh,:) = bound_n
!    end if
!    !
!    !~~ north exchange:
!    !          edge    
!    !     cell       |
!    !          y1+1  +
!    !     y1+1       |
!    !          y1    +    + (y2+1)
!    !     y1         |    :      
!    !          y1-1  +    + y2
!    !                :    |      y2
!    !         (y1-2) + <- + y2-1
!    !                     |      
!    !                            cell
!    !                       edge
!    !
!    idim = 2  ! south-north
!    idir = 2  ! next (=north)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_n = field(x1:x2,y2-nh:y2-1,:)
!    !~ send north, receive south:
!    call MPI_SendRecv( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(x1:x2,y1-1-nh:y1-2,:) = bound_s
!    end if
!
!    ! storage:
!    deallocate( bound_s, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( bound_n, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!
!#endif
!
!    ! ok
!    status = 0
!    
!  end subroutine Domains_ExchangeHaloEdgeY_3d_r
  
  ! *
  
  !
  ! For selected cartesian domain index 'coord',
  ! send slab of width 'nh' to, or receive a slab from
  ! a neighbour. What to do exacly is specified by the key:
  !   'from-west'   : get  slice from west
  !   'to-west'     : send slice to   west
  ! and similar for east/south/north.
  !
  
  subroutine Domains_TransferSlab_3d_r( self, coord, key, nh, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_SendRecv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_TransferSlab_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  coord          ! cartesian coordinate in dimension 1
    character(len=*), intent(in)  ::  key
    integer, intent(in)           ::  nh             ! number of halo cells
    real, intent(inout)           ::  field(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    logical                     ::  do_transfer
    integer                     ::  nx, ny, nz
    real, allocatable           ::  sendbuf(:,:,:)  ! (nh,ny,nz)
    real, allocatable           ::  recvbuf(:,:,:)  ! (nh,ny,nz)
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  source, dest
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif

        integer :: j
    
    ! --- begin -----------------------------------
    
#ifdef _MPI

    ! switch:
    select case ( key )
      case ( 'from-west', 'to-west', 'from-east', 'to-east' )
        idim = 1
      case ( 'from-south', 'to-south', 'from-north', 'to-north' )
        idim = 2
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! switch:
    select case ( key )
      case ( 'from-west', 'to-west', 'from-south', 'to-south' )
        do_transfer = (coord > 0) .and. &
                      ((self%cart_coords(idim) == coord-1) .or. (self%cart_coords(idim) == coord))
      case ( 'from-east', 'to-east', 'from-north', 'to-north' )
        do_transfer = (coord < self%cart_dims(idim)-1) .and. &
                      ((self%cart_coords(idim) == coord) .or. (self%cart_coords(idim) == coord+1))
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select
    
    !! testing ...
    !write (gol,*) '  www transfer ', trim(key), ' ', transfer; call goPr

    ! local domain involved in transfer?
    if ( do_transfer ) then

      ! dimensions:
      nx = size(field,1)
      ny = size(field,2)
      nz = size(field,3)

      ! storage:
      select case ( key )
        case ( 'from-west', 'to-west', 'from-east', 'to-east' )
          allocate( sendbuf(nh,ny,nz), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( recvbuf(nh,ny,nz), stat=status )
          IF_NOT_OK_RETURN(status=1)
        case ( 'from-south', 'to-south', 'from-north', 'to-north' )
          allocate( sendbuf(nx,nh,nz), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( recvbuf(nx,nh,nz), stat=status )
          IF_NOT_OK_RETURN(status=1)
        case default
          write (gol,'("unsupported key `",a,"`")'); call goErr
          TRACEBACK; status=1; return
      end select


      ! data type:
      call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
      IF_NOT_OK_RETURN(status=1)

      ! tag number:
      tag = 1
      
      ! send direction:
      select case ( key )
        !~  o  <--  coord   <--  o
        case ( 'to-west', 'from-east', 'to-south', 'from-north' )
          idir = 1
        !~  o  -->  coord   -->  o
        case ( 'from-west', 'to-east', 'from-south', 'to-north' )
          idir = 2
        case default
          write (gol,'("unsupported key `",a,"`")'); call goErr
          TRACEBACK; status=1; return
      end select

      ! fill send buffer:
      select case ( key )
        !~ 'coord' is the receiving domain:
        case ( 'from-west', 'from-east', 'from-south', 'from-north' )
          ! is this 'coord' ?
          if ( self%cart_coords(idim) == coord ) then
            ! this is the destination domain, define source pe:
            source  = self%cart_src_id(idim,idir)
            ! no data to be send:
            dest    = -2
          else
            ! define destination pe:
            dest    = self%cart_dst_id(idim,idir)
            ! no data is received:
            source  = -2
          end if
        !~ 'coord' is the sending domain:
        case ( 'to-west', 'to-east', 'to-south', 'to-north' )
          ! is this 'coord' ?
          if ( self%cart_coords(idim) == coord ) then
            ! define destination pe:
            dest    = self%cart_dst_id(idim,idir)
            ! no data is received:
            source  = -2
          else
            ! this is the destination domain, define source pe:
            source  = self%cart_src_id(idim,idir)
            ! no data to be send:
            dest    = -2
          end if
        !
        case default
          write (gol,'("unsupported key `",a,"`")'); call goErr
          TRACEBACK; status=1; return
      end select
      
      ! need to fill send buffer?
      if ( dest >= 0 ) then
        ! fill send buffer:
        select case ( key )
          case ( 'to-west', 'from-east' )
            sendbuf = field(1:nh,:,:)
          case ( 'from-west', 'to-east' )
            sendbuf = field(nx-nh+1:nx,:,:)
          case ( 'to-south', 'from-north' )
            sendbuf = field(:,1:nh,:)
          case ( 'from-south', 'to-north' )
            sendbuf = field(:,ny-nh+1:ny,:)
          case default
            write (gol,'("unsupported key `",a,"`")'); call goErr
            TRACEBACK; status=1; return
        end select
      end if
      
      !! testing ...
      !write (gol,*) '  www send to ', dest, ', receive from ', source; call goPr

      ! on some machines send/recv does not work for undefined dest/source ...
      if ( (dest >= 0) .and. (source >= 0) ) then

        ! send and receive:
        call MPI_SendRecv( sendbuf, size(sendbuf), mpi_dtype, dest  , tag, &
                           recvbuf, size(recvbuf), mpi_dtype, source, tag, &
                           self%cart_comm, mpi_stat, status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( (dest >= 0) .and. (source < 0) ) then

        ! send and receive:
        call MPI_Send( sendbuf, size(sendbuf), mpi_dtype, dest  , tag, &
                       self%cart_comm, status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( (dest < 0) .and. (source >= 0) ) then

        ! send and receive:
        call MPI_Recv( recvbuf, size(recvbuf), mpi_dtype, source, tag, &
                       self%cart_comm, mpi_stat, status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if
      ! anything received?
      if ( source >= 0 ) then
        ! fill result:
        select case ( key )
          case ( 'from-west', 'to-east' )
            field(1:nh,:,:) = recvbuf
          case ( 'from-east', 'to-west' )
            field(nx-nh+1:nx,:,:) = recvbuf
          case ( 'from-south', 'to-north' )
            field(:,1:nh,:) = recvbuf
          case ( 'from-north', 'to-south' )
            field(:,ny-nh+1:ny,:) = recvbuf
          case default
            write (gol,'("unsupported key `",a,"`")'); call goErr
            TRACEBACK; status=1; return
        end select
        !! testing ...
        !do j = 10, 1, -1
        !  write (gol,*) '  www recv ', j, field(1,j,1); call goPr
        !end do
      end if

      ! clear:
      deallocate( sendbuf, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( recvbuf, stat=status )
      IF_NOT_OK_RETURN(status=1)

    end if  ! send/recv pe
#endif

    ! ok
    status = 0
    
  end subroutine Domains_TransferSlab_3d_r
  
  ! *
  
  !
  ! Transfer field to or receive field from neighbour given key:
  !   'from-west' , 'to-west' ,
  !   'from-east' , 'to-east' ,
  !   'from-south', 'to-south',
  !   'from-north', 'to-north'
  ! Do this only for domains with the given coordinate id,
  ! where the dimension (west/east, or north/south) is given by the key.
  ! 
  
  subroutine Domains_TransferData_3d_r( self, coord, key, field, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send, MPI_Recv
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_TransferData_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  coord          ! cartesian coordinate in dimension 1
    character(len=*), intent(in)  ::  key
    real, intent(inout)           ::  field(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer                     ::  nx, ny, nz
    integer                     ::  idim
    integer                     ::  idir
    integer                     ::  source, dest
    integer                     ::  tag
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
    type(MPI_Status)            ::  mpi_stat
#endif
    
    ! --- begin -----------------------------------
    
#ifdef _MPI

    ! switch:
    select case ( key )
      case ( 'from-west', 'to-west', 'from-east', 'to-east' )
        idim = 1
      case ( 'from-south', 'to-south', 'from-north', 'to-north' )
        idim = 2
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! dimensions:
    nx = size(field,1)
    ny = size(field,2)
    nz = size(field,3)

    ! data type:
    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! tag number:
    tag = 1

    ! send direction:
    select case ( key )
      !~  o  <--  coord   <--  o
      case ( 'to-west', 'from-east', 'to-south', 'from-north' )
        idir = 1
      !~  o  -->  coord   -->  o
      case ( 'from-west', 'to-east', 'from-south', 'to-north' )
        idir = 2
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! fill send buffer:
    select case ( key )
      !~ 'coord' is the receiving domain:
      case ( 'from-west', 'from-south', 'from-east', 'from-north' )
        !~ is this 'coord' ?
        if ( self%cart_coords(idim) == coord ) then
          ! no data to be send:
          dest      = -2
          ! this is the destination domain, define source pe:
          source    = self%cart_src_id(idim,idir)
        !~ is this the source?
        else if ( self%cart_coords(idim) == coord-dir_shift(idir) ) then
          ! define destination pe:
          dest      = self%cart_dst_id(idim,idir)
          ! no data is received:
          source    = -2
        !~ other:
        else
          ! not involved in exchange:
          dest      = -2
          source    = -2
        end if
      !~ 'coord' is the sending domain:
      case ( 'to-west', 'to-south', 'to-east', 'to-north' )
        !~ is this 'coord' ?
        if ( self%cart_coords(idim) == coord ) then
          ! define destination pe:
          dest      = self%cart_dst_id(idim,idir)
          ! no data is received:
          source    = -2
        !~ is this the destination?
        else if ( self%cart_coords(idim) == coord-dir_shift(idir) ) then
          ! no data to be send:
          dest      = -2
          ! this is the destination domain, define source pe:
          source    = self%cart_src_id(idim,idir)
        !~ other:
        else
          ! not involved in exchange:
          dest      = -2
          source    = -2
        end if
      !~
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! send or receive:
    if ( (dest >= 0) .and. (source < 0) ) then
    
      ! send to dest:
      call MPI_Send( field, size(field), mpi_dtype, dest, tag, &
                       self%cart_comm, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (dest < 0) .and. (source >= 0) ) then

      ! recv from source:
      call MPI_Recv( field, size(field), mpi_dtype, source, tag, &
                       self%cart_comm, mpi_stat, status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( (dest >= 0) .and. (source >= 0) ) then
      write (gol,'("could not both send (to ",i0,") and recv (from ",i0,")")') dest, source; call goErr
      TRACEBACK; status=1; return
    end if

#endif

    ! ok
    status = 0
    
  end subroutine Domains_TransferData_3d_r
  
  ! *
  
  subroutine Domains_TransferDataAll_3d_r( self, key, values, status )
  
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_TransferDataAll_3d_r'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    character(len=*), intent(in)  ::  key
    real, intent(inout)           ::  values(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    integer               ::  idim
    integer               ::  cart_size, cart_coord
    integer               ::  coord

    ! --- begin -----------------------------------

    ! switch:
    select case ( key )
      case ( 'from-west', 'to-west', 'from-east', 'to-east' )
        idim = 1
      case ( 'from-south', 'to-south', 'from-north', 'to-north' )
        idim = 2
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! Info on domain layout in this dimension:
    call self%GetDim( idim, status, cart_size=cart_size, cart_coord=cart_coord )
    IF_NOT_OK_RETURN(status=1)

    ! switch:
    select case ( key )

      !~ low to high:
      case ( 'from-west', 'to-east', 'from-south', 'to-north' )
        ! loop over domains:
        do coord = 0, cart_size-1
          ! transfer values according to the description:
          call self%TransferData( coord, key, values, status )
          IF_NOT_OK_RETURN(status=1)
        end do

      !~ high to low:
      case ( 'to-west', 'from-east', 'to-south', 'from-north' )
        ! loop over domains:
        do coord = cart_size-1, 0, -1
          ! transfer values according to the description:
          call self%TransferData( coord, key, values, status )
          IF_NOT_OK_RETURN(status=1)
        end do

      !~
      case default
        write (gol,'("unsupported key `",a,"`")'); call goErr
        TRACEBACK; status=1; return
    end select

    
    ! ok
    status = 0
    
  end subroutine Domains_TransferDataAll_3d_r
  
  ! *
  
!  subroutine Domains_ExchangeHaloY_3d_r( self, nh, field, status )
!  
!    use GO_Comm, only : goc
!#ifdef _MPI
!    use MPI_F08, only : MPI_DataType
!    use MPI_F08, only : MPI_Status
!    use MPI_F08, only : MPI_SendRecv
!#endif
!
!    ! --- const ---------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Domains_ExchangeHaloY_3d_r'
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_Domains), intent(in)  ::  self
!    integer, intent(in)           ::  nh             ! number of halo cells
!    real, pointer                 ::  field(:,:,:)   ! (:,1-nh:ny+nh,:)
!    integer, intent(out)          ::  status
!    
!    ! --- local -----------------------------------
!  
!    integer                     ::  y1, y2
!    integer                     ::  nx, nz
!    real, allocatable           ::  bound_n(:,:,:)  ! (1-nh:nx+nh,nh,:)
!    real, allocatable           ::  bound_s(:,:,:)  ! (1-nh:nx+nh,nh,:)
!    integer                     ::  idim
!    integer                     ::  idir
!    integer                     ::  tag
!#ifdef _MPI
!    type(MPI_DataType)          ::  mpi_dtype
!    type(MPI_Status)            ::  mpi_stat
!#endif
!    
!    ! --- begin -----------------------------------
!
!#ifdef _MPI
!    
!    ! short variables for local index spaces:
!    y1 = self%lbo(2,goc%id)
!    y2 = self%ubo(2,goc%id)
!    ! other dimensions:
!    nx = size(field,1)
!    nz = size(field,3)
!
!    ! check ...
!    if ( (y1-nh < lbound(field,2)) .or. (y2+nh > ubound(field,2))        ) then
!      write (gol,'("y index space of field (",i0,":",i0,")")') lbound(field,2), ubound(field,2); call goErr
!      write (gol,'("does not cover expected index space (",i0,":",i0,")")') self%lbo(2,goc%id), self%ubo(2,goc%id); call goErr
!      write (gol,'("with ",i0," halo (",i0,":",i0,")")') nh, self%lbo(2,goc%id)-nh, self%ubo(2,goc%id)+nh; call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! storage:
!    allocate( bound_s(nx,nh,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( bound_n(nx,nh,nz), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! data type:
!    call goc%GetDataType( 'real', kind(field), mpi_dtype, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! init tag number (not necessary?)
!    tag = 0
!
!    !~~ south exchange:
!    !             
!    !     cell  |
!    !           +    +
!    !     y1+1  | -> : y2+2
!    !           +    + 
!    !     y1    | -> : y2+1 
!    !           +    + 
!    !           :    | y2
!    !           +    + 
!    !                | 
!    !                  cell
!    !                  
!    !
!    idim = 2  ! south-north
!    idir = 1  ! previous (=south)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_s = field(:,y1:y1+nh-1,:)
!    ! send south, receive north:
!    call MPI_SendRecv( bound_s, size(bound_s), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_n, size(bound_n), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(:,y2+1:y2+nh,:) = bound_n
!    end if
!    !
!    !~~ north exchange:
!    !            
!    !     cell |
!    !          +    +
!    !     y1+1 |    :
!    !          +    + 
!    !     y1   |    : 
!    !          +    + 
!    !          : <- | y2
!    !          +    + 
!    !          : <- | y2-1 
!    !          +    +  
!    !                 cell
!    !
!    idim = 2  ! south-north
!    idir = 2  ! next (=north)
!    ! increase tag:
!    tag = tag + 1
!    ! fill send buffer from local field:
!    bound_n = field(:,y2-nh+1:y2,:)
!    !~ send north, receive south:
!    call MPI_SendRecv( bound_n, size(bound_n), mpi_dtype, self%cart_dst_id(idim,idir), tag, &
!                       bound_s, size(bound_s), mpi_dtype, self%cart_src_id(idim,idir), tag, &
!                       self%cart_comm, mpi_stat, status )
!    IF_MPI_NOT_OK_RETURN(status=1)
!    ! any received?
!    if ( self%cart_src_id(idim,idir) >= 0 ) then
!      ! store in halo:
!      field(:,y1-nh:y1-1,:) = bound_s
!    end if
!
!#endif
!
!    ! ok
!    status = 0
!    
!  end subroutine Domains_ExchangeHaloY_3d_r
  
  
  ! *

  
  subroutine ExchangeTest( self, status )
  
    use GO_Comm, only : goc

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/ExchangeTest'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
  
    real, pointer        ::  field(:,:,:)
    integer              ::  nx, ny, nh
    character(len=32)    ::  fmt
    integer              ::  j
    
    ! --- begin -----------------------------------
    
    ! grid size:
    nx = self%shp(1,goc%id)
    ny = self%shp(2,goc%id)
    ! halo:
    nh = 2
    
    ! ~ cells
  
    ! test field including halo:
    allocate( field(1-nh:nx+nh,1-nh:ny+nh,1), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! fill with processor id in core:
    field = -9
    field(1:nx,1:ny,:) = goc%id

    ! exchange halo:
    call self%ExchangeHalo( nh, field, status )
    IF_NOT_OK_RETURN(status=1)

    ! show:
    write (gol,'(" ")'); call goPr
    write (gol,'("test exchange of size-",i0," halo ...")') nh; call goPr
    write (gol,'("-----------------------------------------------------------------------")'); call goPr
    write (fmt,'("(4x,",i0,"i3)")') nx+2*nh
    do j = ubound(field,2), lbound(field,2), -1
      write (gol,fmt) int(field(:,j,1)); call goPr
    end do
    write (gol,'("-----------------------------------------------------------------------")'); call goPr

    ! clear:
    deallocate( field, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    !! ~ edge x
    !
    !! test field including halo:
    !allocate( field(0-nh:nx+nh,1:ny,1), stat=status )
    !IF_NOT_OK_RETURN(status=1)
    !
    !! fill with processor id in core:
    !field = -9
    !field(0:nx,1:ny,:) = goc%id
    !
    !! exchange halo:
    !call self%ExchangeHaloEdgeX( nh, field, status )
    !IF_NOT_OK_RETURN(status=1)
    !
    !! show:
    !write (gol,'(" ")'); call goPr
    !write (gol,'("test exchange of size-",i0," halo at x-edge ...")') nh; call goPr
    !write (gol,'("-----------------------------------------------------------------------")'); call goPr
    !write (fmt,'("(4x,",i0,"i3)")') nx+1+2*nh
    !do j = ubound(field,2), lbound(field,2), -1
    !  write (gol,fmt) int(field(:,j,1)); call goPr
    !end do
    !write (gol,'("-----------------------------------------------------------------------")'); call goPr
    !
    !! clear:
    !deallocate( field, stat=status )
    !IF_NOT_OK_RETURN(status=1)
    !
    !! ~ edge y
    !
    !! test field including halo:
    !allocate( field(1:nx,0-nh:ny+nh,1), stat=status )
    !IF_NOT_OK_RETURN(status=1)
    !
    !! fill with processor id in core:
    !field = -9
    !field(1:nx,0:ny,:) = goc%id
    !
    !! exchange halo:
    !call self%ExchangeHaloEdgeY( nh, field, status )
    !IF_NOT_OK_RETURN(status=1)
    !
    !! show:
    !write (gol,'(" ")'); call goPr
    !write (gol,'("test exchange of size-",i0," halo at y-edge ...")') nh; call goPr
    !write (gol,'("-----------------------------------------------------------------------")'); call goPr
    !write (fmt,'("(4x,",i0,"i3)")') nx+1+2*nh
    !do j = ubound(field,2), lbound(field,2), -1
    !  write (gol,fmt) int(field(:,j,1)); call goPr
    !end do
    !write (gol,'("-----------------------------------------------------------------------")'); call goPr
    !write (gol,'(" ")'); call goPr
    !
    !! clear:
    !deallocate( field, stat=status )
    !IF_NOT_OK_RETURN(status=1)

    !! testing ...
    !call goc%Barrier( status )
    !IF_NOT_OK_RETURN(status=1)
    !write (gol,'("break")'); call goErr
    !TRACEBACK; status=1; return
    
    ! ok
    status = 0
    
  end subroutine ExchangeTest

  


  ! ***************************************************************************
  ! ***
  ! *** swap decomposition
  ! ***
  ! ***************************************************************************


  subroutine Domains_Swap_2d_r8( self, input, doms, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_2d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 2
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,goc%id) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(goc%id)) )
    allocate( sendcounts(0:goc%npes-1) )
    allocate( sdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,1(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(goc%id) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(goc%id)) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( rdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(goc%id) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
                        recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                        goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,1(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_2d_r8


  ! ***

  subroutine Domains_Swap_3d_r8( self, input, doms, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_3d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 3
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,goc%id) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(goc%id)) )
    allocate( sendcounts(0:goc%npes-1) )
    allocate( sdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2),&
                                                  ilbo(3):iubo(3)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(goc%id) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(goc%id)) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( rdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(goc%id) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
                        recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                        goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2),&
               ilbo(3):iubo(3)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_3d_r8


  ! ***


  subroutine Domains_Swap_4d_r8( self, input, doms, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_4d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 4
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:,:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,goc%id) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(goc%id)) )
    allocate( sendcounts(0:goc%npes-1) )
    allocate( sdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,3(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2),&
                                                  ilbo(3):iubo(3),&
                                                  ilbo(4):iubo(4)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(goc%id) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(goc%id)) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( rdispls   (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(goc%id) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(goc%id); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
                        recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                        goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,3(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2),&
               ilbo(3):iubo(3),&
               ilbo(4):iubo(4)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_4d_r8
  


  ! ***************************************************************************
  ! ***
  ! *** gather on root
  ! ***
  ! ***************************************************************************


  subroutine Domains_GatherV_1d_i( self, sendbuf, n, recvbuf, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_1d_i'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)      ::  self
    integer, intent(in)               ::  sendbuf(:)
    integer, intent(in)               ::  n
    integer, intent(out), allocatable ::  recvbuf(:)
    integer, intent(out)              ::  status
    
    ! --- local -----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
#endif
    integer                     ::  gn
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    !check
    if ( n < 0 .or. n > size(sendbuf) ) then
      write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf), n ;call goErr
      TRACEBACK;status=1;return
    end if 
    
    ! ~ prepare receive counts
    if ( goc%root ) then
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( recvcounts(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! total number over all processors, ans start index in output:
    call goc%Gather( n, recvcounts, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! free memory  for receive buffer  
    if ( allocated( recvbuf) ) deallocate(recvbuf)
    
    ! target?
    if ( goc%root ) then
      
      ! total output
      gn = sum(recvcounts)
      
      ! storage:s
      allocate( recvbuf(gn), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init displacement:
      displ = 0
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! displacement:
        displs(pid) = displ
        ! increase displacement:
        displ = displ + recvcounts(pid)
      end do ! iproc
    
    else
    
      ! dummy storage ...
      allocate( recvbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs (1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ~ 
    
#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_GatherV( sendbuf, n                 , mpi_dtype, &
                      recvbuf, recvcounts, displs, mpi_dtype, &
                      goc%root_id, goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)
        
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_1d_i

  
  ! * 
  
  
  subroutine Domains_GatherV_1d_i_multiple( self, sendbuf, n, recvbuf, ind , status)
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_1d_i_multiple'
      
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)       ::  self
    integer, intent(in)                ::  sendbuf(:,:)
    integer, intent(in)                ::  n
    integer, intent(out), allocatable  ::  recvbuf(:,:)
    integer, intent(in)                ::  ind
    integer, intent(out)               ::  status
    
    ! --- local -----------------------------------
    
    integer                            ::  i
    integer                            ::  nmulti, gn
    integer, allocatable               ::  recvbuf_tmp(:)
    integer, allocatable               ::  recvcounts(:)
    
        
    ! --- begin -----------------------------------
    
    ! check sizes
    if ( ind == 1 ) then
      ! check size of sendbuf
      if ( size(sendbuf(ind,:)) < 0 .or. size(sendbuf(ind,:)) < n )  then
        write( gol, '("Length of input array ",i0," not compliant with given input ",i0)') size(sendbuf(ind,:)), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti = size(sendbuf,1)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if

      ! total number over all processors, ans start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root ) then
        ! storage:
        allocate( recvbuf( nmulti,gn), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( nmulti,1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if

      
    else if ( ind == 2 ) then
      ! check size of sendbuf
      if ( size(sendbuf(:,ind)) < 0 .or. size(sendbuf(:,ind)) < n )  then
        write( gol, '("Length of input array ",i0," not compliant with given input ",i0)') size(sendbuf(:,ind)), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti = size(sendbuf,2)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! total number over all processors, ans start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root) then
        ! storage:
        allocate( recvbuf( gn,nmulti), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( 1, nmulti), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
    else 
      ! error
      write( gol, '("Higher multiple index given ",i0", regarding input array with dimension 2")' ) ind ; call goErr
      TRACEBACK;status=1;return
      
    end if 
    
    do i = 1, nmulti
      ! 
      if ( ind == 1 ) then

        ! gather over total dimension
        call self%GatherV( sendbuf(i,:), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(i,:) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else if ( ind == 2 ) then

        ! gather over total dimension
        call self%GatherV( sendbuf(:,i), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(:,i) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else 
        ! error
        write( gol, '("Higher multiple index given ",i0", regarding input array with dimension 2")' ) ind ; call goErr
        TRACEBACK;status=1;return

      end if   
    end do  ! multiple dimensions
    
    ! clear :
    deallocate( recvcounts, stat=status)
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_1d_i_multiple

  
  ! *
  
  
  subroutine Domains_GatherV_1d_r( self, sendbuf, n, recvbuf, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_1d_r'
      
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)    ::  self
    real, intent(in)                ::  sendbuf(:)
    integer, intent(in)             ::  n
    real, intent(out), allocatable  ::  recvbuf(:)
    integer, intent(out)            ::  status
    
    ! --- local -----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
#endif
    integer                     ::  gn
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    !check
    if ( n < 0 .or. n > size(sendbuf) ) then
      write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf), n ;call goErr
      TRACEBACK;status=1;return
    end if 
    
    ! ~ prepare receive counts
    if ( goc%root ) then
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( recvcounts(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! total number over all processors, ans start index in output:
    call goc%Gather( n, recvcounts, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! free memory  for receive buffer  
    if ( allocated( recvbuf) ) deallocate(recvbuf)
    
    ! target?
    if ( goc%root ) then
      
      ! total output
      gn = sum(recvcounts)
      
      ! storage:s
      allocate( recvbuf(gn), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init displacement:
      displ = 0
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! displacement:
        displs(pid) = displ
        ! increase displacement:
        displ = displ + recvcounts(pid)
      end do ! iproc
    
    else
    
      ! dummy storage ...
      allocate( recvbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs (1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ~ 
    
#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_GatherV( sendbuf, n                 , mpi_dtype, &
                      recvbuf, recvcounts, displs, mpi_dtype, &
                      goc%root_id, goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)
        
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_1d_r

    
  ! * 
  
  
  subroutine Domains_GatherV_1d_r_multiple( self, sendbuf, n, recvbuf, ind , status)
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_1d_r_multiple'
      
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)    ::  self
    real, intent(in)                ::  sendbuf(:,:)
    integer, intent(in)             ::  n
    real, intent(out), allocatable  ::  recvbuf(:,:)
    integer, intent(in)             ::  ind
    integer, intent(out)            ::  status
    
    ! --- local -----------------------------------
    
    integer                         ::  i
    integer                         ::  nmulti, gn
    real, allocatable               ::  recvbuf_tmp(:)
    integer, allocatable            ::  recvcounts(:)
        
    ! --- begin -----------------------------------
    
    ! check sizes
    if ( ind == 1 ) then
      ! check size of sendbuf
      if ( size(sendbuf,2) < 0 .or. size(sendbuf,2) < n )  then
        write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf,2), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti = size(sendbuf,1)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! total number over all processors, ans start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root ) then
        ! storage:
        allocate( recvbuf( nmulti,gn), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( nmulti,1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if     
      
    else if ( ind == 2 ) then
      ! check size of sendbuf
      if ( size(sendbuf,1) < 0 .or. size(sendbuf,1) > n )  then
        write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf,1), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti = size(sendbuf,2)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! total number over all processors, and start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root) then
        ! storage:
        allocate( recvbuf( gn,nmulti), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( 1, nmulti), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if      
      
    else 
      ! error
      write( gol, '("Higher multiple index given ",i0", regarding input array with dimension 2")' ) ind ; call goErr
      TRACEBACK;status=1;return
      
    end if 
    
    do i = 1, nmulti
      ! 
      if ( ind == 1 ) then
        call self%GatherV( sendbuf(i,:), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(i,:) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else if ( ind == 2 ) then

        ! gather over total dimension
        call self%GatherV( sendbuf(:,i), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(:,i) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else 
        ! error
        write( gol, '("Higher multiple index given ",i0", regarding input array with dimension 2")' ) ind ; call goErr
        TRACEBACK;status=1;return

      end if   
    end do  ! multiple dimensions
  
    ! clear :
    deallocate( recvcounts, stat=status)
    IF_NOT_OK_RETURN(status=1)  
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_1d_r_multiple

  
  ! * 
  
  
  subroutine Domains_GatherV_1d_r_multiple2d( self, sendbuf, n, recvbuf, inds , status)
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_1d_r_multiple2d'
      
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)    ::  self
    real, intent(in)                ::  sendbuf(:,:,:)
    integer, intent(in)             ::  n
    real, intent(out), allocatable  ::  recvbuf(:,:,:)
    integer, intent(in)             ::  inds(2)
    integer, intent(out)            ::  status
    
    ! --- local -----------------------------------
    
    integer                         ::  i,j
    integer                         ::  nmulti(2), gn
    real, allocatable               ::  recvbuf_tmp(:)
    integer, allocatable            ::  recvcounts(:)
        
    ! --- begin -----------------------------------
    
    ! check sizes
    if ( inds(1) == 1 .and. inds(2) == 2 ) then
      ! check size of sendbuf
      if ( size(sendbuf,3) < 0 .or. size(sendbuf,3) < n )  then
        write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf,3), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti(1) = size(sendbuf,1)
      nmulti(2) = size(sendbuf,2)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! total number over all processors, ans start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root ) then
        ! storage:
        allocate( recvbuf( nmulti(1),nmulti(2),gn), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( nmulti(1),nmulti(2),1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if     
      
    else if ( inds(1) == 1 .and. inds(2) == 3 ) then
      ! check size of sendbuf
      if ( size(sendbuf,2) < 0 .or. size(sendbuf,2) > n )  then
        write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf,2), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti(1) = size(sendbuf,1)
      nmulti(2) = size(sendbuf,3)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! total number over all processors, and start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root) then
        ! storage:
        allocate( recvbuf( nmulti(1),n,nmulti(2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( nmulti(1),1,nmulti(2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if      
      
    else if ( inds(1) == 2 .and. inds(2) == 3 ) then
      ! check size of sendbuf
      if ( size(sendbuf,1) < 0 .or. size(sendbuf,1) > n )  then
        write( gol, '("Length of input array ",i0,"not compliant with given input ",i0)') size(sendbuf,1), n ;call goErr
        TRACEBACK;status=1;return
      end if
      nmulti(1) = size(sendbuf,2)
      nmulti(2) = size(sendbuf,3)
      
      ! ~ prepare receive counts
      if ( goc%root ) then
        allocate( recvcounts(0:goc%npes-1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( recvcounts(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! total number over all processors, and start index in output:
      call goc%Gather( n, recvcounts, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! total output
      gn = sum(recvcounts)
      
      if ( goc%root) then
        ! storage:
        allocate( recvbuf( n,nmulti(1),nmulti(2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else 
        ! dummy storaga
        allocate( recvbuf( 1,nmulti(1),nmulti(2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if      
      
    else 
      ! error
      write( gol, '("Unknown combination of multiple indices given ",i0," ", i0, ", regarding input array with dimension 3")' ) inds(1),inds(2) ; call goErr
      TRACEBACK;status=1;return
      
    end if 
    
    do i = 1, nmulti(1)
    do j = 1, nmulti(2)
      ! 
      if ( inds(1) == 1 .and. inds(2) == 2 ) then
      
        ! gather over total dimension
        call self%GatherV( sendbuf(i,j,:), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(i,j,:) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else if ( inds(1) == 1 .and. inds(2) == 3 ) then
      
        ! gather over total dimension
        call self%GatherV( sendbuf(i,:,j), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(i,:,j) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else if ( inds(1) == 2 .and. inds(2) == 3 ) then
      
        ! gather over total dimension
        call self%GatherV( sendbuf(:,i,j), n, recvbuf_tmp, status )
        IF_NOT_OK_RETURN(status=1)
        recvbuf(:,i,j) = recvbuf_tmp
        
        ! clear:
        deallocate( recvbuf_tmp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
      else 
        ! error
        write( gol, '("Unknown combination of multiple indices given ",i0," ", i0, ", regarding input array with dimension 3")' ) inds(1),inds(2) ; call goErr
        TRACEBACK;status=1;return
      
      end if   
    end do  ! multiple dimensions (i)
    end do  ! multiple dimensions (j)
  
    ! clear :
    deallocate( recvcounts, stat=status)
    IF_NOT_OK_RETURN(status=1)  
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_1d_r_multiple2d


  ! *


  subroutine Domains_GatherV_2d_r( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_2d_r'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 2
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real, intent(in)              ::  input(:,:)
    real, intent(out)             ::  output(:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
#endif
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    real, allocatable           ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    real, allocatable           ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! destination?
    if ( goc%root ) then
      ! check ...
      if ( any( shape(output) /= gshp ) ) then
        write (gol,'("shape of output (",i0,",",i0,") while assumed global shape (",i0,",",i0,")")') &
                      shape(output), gshp; call goErr
        TRACEBACK; status=1; return
      end if      
      ! total:
      gn = product( gshp )
      ! offset, arguments are 1-based:
      goff = glbo - 1
    else
      ! dummy value to avoid warnings ...
      goff = 0
    end if
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id),&
                               self%lbo(2,goc%id):self%ubo(2,goc%id)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ~ prepare receive buffer
    
    ! target?
    if ( goc%root ) then
    
      ! storage:
      allocate( recvbuf(gn), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init displacement:
      displ = 0
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! source defined ?
        if ( self%n(pid) > 0 ) then
          ! all values:
          n = self%n(pid)
          !! info ...
          !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
        else
          ! empty:
          n = 0
        end if
        ! number of values:
        recvcounts(pid) = n
        ! displacement:
        displs(pid) = displ
        ! increase displacement:
        displ = displ + n
      end do ! iproc
      ! check ...
      if ( displ /= gn ) then
        write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
        TRACEBACK; status=1; return
      end if
    
    else
    
      ! dummy storage ...
      allocate( recvbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ~ 
    
#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_GatherV( sendbuf, sendcount         , mpi_dtype, &
                      recvbuf, recvcounts, displs, mpi_dtype, &
                      goc%root_id, goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! destination?
    if ( goc%root ) then
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! extract:
        n     = recvcounts(pid)
        displ = displs(pid)
        ! any received from this process ?
        if ( n > 0 ) then
          !! info ...
          !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
          !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
          !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
          ! copy shape, needed to avoid strange compilation error ... 
          lshp = self%shp(:,pid)
          ! extract output slab from buffer:
          output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1),&
                 self%glbo(2,pid)-goff(2):self%gubo(2,pid)-goff(2)) = &
              reshape( recvbuf(displ+1:displ+n), lshp )
        end if
      end do ! iproc
    end if  ! root

    ! clear:
    deallocate( recvbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_2d_r


  ! *
  
  
  subroutine Domains_GatherV_3d_r( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_3d_r'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 2
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real, intent(in)              ::  input(:,:,:)
    real, intent(out)             ::  output(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
#endif
    integer                     ::  nx, ny, nz
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    real, allocatable           ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    real, allocatable           ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! dims:
    nx = size(input,1)
    ny = size(input,2)
    nz = size(input,3)
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( (/nx,ny/) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,",",i0,",:) while expected (",i0,",",i0,",:)")') &
                      nx,ny, self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! destination?
    if ( goc%root ) then
      ! check ...
      if ( any( shape(output) /= (/gshp(1),gshp(2),nz/) ) ) then
        write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                      shape(output), gshp,nz; call goErr
        TRACEBACK; status=1; return
      end if      
      ! total:
      gn = product( gshp ) * nz
      ! offset, arguments are 1-based:
      goff = glbo - 1
    else
      ! dummy value to avoid warnings ...
      goff = 0
    end if
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id) * nz
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id),&
                               self%lbo(2,goc%id):self%ubo(2,goc%id),:), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ~ prepare receive buffer
    
    ! target?
    if ( goc%root ) then
    
      ! storage:
      allocate( recvbuf(gn), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init displacement:
      displ = 0
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! source defined ?
        if ( self%n(pid) > 0 ) then
          ! all values:
          n = self%n(pid) * nz
          !! info ...
          !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
        else
          ! empty:
          n = 0
        end if
        ! number of values:
        recvcounts(pid) = n
        ! displacement:
        displs(pid) = displ
        ! increase displacement:
        displ = displ + n
      end do ! iproc
      ! check ...
      if ( displ /= gn ) then
        write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
        TRACEBACK; status=1; return
      end if
    
    else
    
      ! dummy storage ...
      allocate( recvbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ~ 
    
#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_GatherV( sendbuf, sendcount         , mpi_dtype, &
                      recvbuf, recvcounts, displs, mpi_dtype, &
                      goc%root_id, goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! destination?
    if ( goc%root ) then
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! extract:
        n     = recvcounts(pid)
        displ = displs(pid)
        ! any received from this process ?
        if ( n > 0 ) then
          !! info ...
          !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
          !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
          !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
          ! copy shape, needed to avoid strange compilation error ... 
          lshp = self%shp(:,pid)
          ! extract output slab from buffer:
          output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1),&
                 self%glbo(2,pid)-goff(2):self%gubo(2,pid)-goff(2),:) = &
              reshape( recvbuf(displ+1:displ+n), (/lshp(1),lshp(2),nz/) )
        end if
      end do ! iproc
    end if  ! root

    ! clear:
    deallocate( recvbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_3d_r


  ! *
  
  
  subroutine Domains_GatherV_4d_r( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType, MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_GatherV_4d_r'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 2
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real, intent(in)              ::  input(:,:,:,:)
    real, intent(out)             ::  output(:,:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)          ::  mpi_dtype
#endif
    integer                     ::  nx, ny, nz, nt
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    real, allocatable           ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    real, allocatable           ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! dims:
    nx = size(input,1)
    ny = size(input,2)
    nz = size(input,3)
    nt = size(input,4)
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( (/nx,ny/) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,",",i0,",:) while expected (",i0,",",i0,",:)")') &
                      nx,ny, self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! destination?
    if ( goc%root ) then
      ! check ...
      if ( any( shape(output) /= (/gshp(1),gshp(2),nz,nt/) ) ) then
        write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,3(",",i0),")")') &
                      shape(output), gshp,nz,nt; call goErr
        TRACEBACK; status=1; return
      end if      
      ! total:
      gn = product( gshp ) * nz * nt
      ! offset, arguments are 1-based:
      goff = glbo - 1
    else
      ! dummy value to avoid warnings ...
      goff = 0
    end if
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id) * nz * nt
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id),&
                               self%lbo(2,goc%id):self%ubo(2,goc%id),:,:), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ~ prepare receive buffer
    
    ! target?
    if ( goc%root ) then
    
      ! storage:
      allocate( recvbuf(gn), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init displacement:
      displ = 0
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! source defined ?
        if ( self%n(pid) > 0 ) then
          ! all values:
          n = self%n(pid) * nz * nt
          !! info ...
          !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
        else
          ! empty:
          n = 0
        end if
        ! number of values:
        recvcounts(pid) = n
        ! displacement:
        displs(pid) = displ
        ! increase displacement:
        displ = displ + n
      end do ! iproc
      ! check ...
      if ( displ /= gn ) then
        write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
        TRACEBACK; status=1; return
      end if
    
    else
    
      ! dummy storage ...
      allocate( recvbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( recvcounts(0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( displs    (0:goc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ~ 
    
#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_GatherV( sendbuf, sendcount         , mpi_dtype, &
                      recvbuf, recvcounts, displs, mpi_dtype, &
                      goc%root_id, goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! destination?
    if ( goc%root ) then
      ! loop over source processes:
      do pid = 0, goc%npes-1
        ! extract:
        n     = recvcounts(pid)
        displ = displs(pid)
        ! any received from this process ?
        if ( n > 0 ) then
          !! info ...
          !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
          !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
          !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
          ! copy shape, needed to avoid strange compilation error ... 
          lshp = self%shp(:,pid)
          ! extract output slab from buffer:
          output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1),&
                 self%glbo(2,pid)-goff(2):self%gubo(2,pid)-goff(2),:,:) = &
              reshape( recvbuf(displ+1:displ+n), (/lshp(1),lshp(2),nz,nt/) )
        end if
      end do ! iproc
    end if  ! root

    ! clear:
    deallocate( recvbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_GatherV_4d_r
  


  ! ***************************************************************************
  ! ***
  ! *** gather and broadcast to all
  ! ***
  ! ***************************************************************************


  subroutine Domains_AllGather_1d_i4( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_INTEGER
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_1d_i4'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 4
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer(wpr), intent(in)      ::  input(:)
    integer(wpr), intent(out)     ::  output(:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    integer(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    integer(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( displs    (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, MPI_INTEGER, &
                         recvbuf, recvcounts, displs, MPI_INTEGER, &
                         goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_1d_i4
  
  
  ! ***


  subroutine Domains_AllGather_1d_r8( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_1d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:)
    real(wpr), intent(out)        ::  output(:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    real(wpr), allocatable      ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    real(wpr), allocatable      ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,") while expected (",i0,")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,") while assumed global shape (",i0,")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( displs    (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
                         recvbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                         goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_1d_r8


  ! ***
  

  subroutine Domains_AllGather_3d_c8( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DOUBLE_COMPLEX
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_3d_c8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 3
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    complex(wpr), intent(in)      ::  input(:,:,:)
    complex(wpr), intent(out)     ::  output(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    complex(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    complex(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! local input slab ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of input (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                      shape(input), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, goc%npes-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(goc%id)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),&
      !           self%lbo(3,goc%id),self%ubo(3,goc%id); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,goc%id):self%ubo(1,goc%id),&
                               self%lbo(2,goc%id):self%ubo(2,goc%id),&
                               self%lbo(3,goc%id):self%ubo(3,goc%id)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:goc%npes-1) )
    allocate( displs    (0:goc%npes-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                         goc%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, goc%npes-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1),&
               self%glbo(2,pid)-goff(2):self%gubo(2,pid)-goff(2),&
               self%glbo(3,pid)-goff(3):self%gubo(3,pid)-goff(3)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_3d_c8
  

!  ! ***************************************************************************
!  ! ***
!  ! *** scatter
!  ! ***
!  ! ***************************************************************************
!
  subroutine Domains_ScatterV_2d_r( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
   use MPI_F08, only : MPI_ScatterV
#endif

   ! --- const ---------------------------------
   
   character(len=*), parameter   ::  rname = mname//'/Domains_ScatterV_2d_r'
   
   ! rank and kind:
   integer, parameter  ::  ndim = 2
 
   ! --- in/out ---------------------------------
   
   class(T_Domains), intent(in)  ::  self
   real, intent(in)              ::  input(:,:)
   real, intent(out)             ::  output(:,:)
   integer, intent(out)          ::  status
   
   ! --- local -----------------------------------
   
#ifdef _MPI
   type(MPI_DataType)          ::  mpi_dtype
#endif
   integer                     ::  nx, ny, splitx,splity,i,j
   integer                     ::  gn
   integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
   integer                     ::  lshp(ndim)
   real, allocatable           ::  sendbuf(:)     ! (input n)
   integer, allocatable        ::  sendcounts(:)  ! (nproc)
   integer, allocatable        ::  displs(:)      ! (nproc)
   real, allocatable           ::  recvbuf(:)     ! (output n)
   integer                     ::  recvcount
   integer                     ::  pid
   integer                     ::  n
   integer                     ::  displ
   
   ! --- begin -----------------------------------
   
   ! dims:
   nx = size(output,1)
   ny = size(output,2)
   

   ! local output slab ?
   if ( all( self%shp(:,goc%id) > 0 ) ) then
     ! check ...
     if ( any( (/nx,ny/) /= self%shp(:,goc%id) ) ) then
       write (gol,'("shape of output (",i0,",",i0,",:) while expected (",i0,",",i0,",:)")') &
                     nx,ny, self%shp(:,goc%id); call goErr
       TRACEBACK; status=1; return
     end if
   !else
   !  ! info ...
   !  write (gol,'("empty output domain, careful ...")'); call goPr
   end if

   ! global bounds:
   glbo = self%glbo(:,0)
   gubo = self%gubo(:,0)
   do pid = 1, goc%npes-1
     if ( self%n(pid) > 0 ) then
       glbo = min( glbo, self%glbo(:,pid) )
       gubo = max( gubo, self%gubo(:,pid) )
     end if
   end do
   ! global shape:
   gshp = gubo - glbo + 1
   ! source?
   if ( goc%root ) then
     ! check ...
     if ( any( shape(input) /= (/gshp(1),gshp(2)/) ) ) then
       write (gol,'("shape of output (",i0,",",i0,") while assumed global shape (",i0,",",i0,")")') &
                     shape(input), gshp; call goErr
       TRACEBACK; status=1; return
     end if      
     ! total:
     gn = product( gshp ) 
     ! offset, arguments are 1-based:
     goff = glbo - 1
   else
     ! dummy value to avoid warnings ...
     goff = 0
   end if
   
      splitx=(gshp(1)/nx)
      splity=(gshp(2)/ny)

   ! ~ create receive buffer
   
   ! local size:
   n = self%n(goc%id) 
   ! store:
   recvcount = n
   ! defined ?
   if ( n > 0 ) then
     ! storage:
     allocate( recvbuf(n), stat=status )
     IF_NOT_OK_RETURN(status=1)
   else
     ! dummy:
     allocate( recvbuf(1), stat=status )
     IF_NOT_OK_RETURN(status=1)
   end if

   ! ~ create send buffer
   
   ! target?
   if ( goc%root ) then
   
     ! storage:
     allocate( sendbuf(gn), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( sendcounts(0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( displs    (0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     
     ! init displacement:
     displ = 0
     ! loop over target processes:
     i=1
     j=1    
     do pid = 0, goc%npes-1
       ! target has values ?
       if ( self%n(pid) > 0 ) then
         ! all values:
         n = self%n(pid) 
       else
         ! empty:
         n = 0
       end if
       ! number of values:
       sendcounts(pid) = n
       ! displacement:
       displs(pid) = displ
       
       ! any target data from this process ?
       if ( n > 0 ) then
         !! info ...
         !write (gol,'("store ",i0," elements in sendbuf(",i0,":",i0,") from input(",i0,":",i0,2(",",i0,":",i0),")")') &
         !        n, displ+1, displ+n, self%lbo(1,pid),self%ubo(1,pid),self%lbo(2,pid),self%ubo(2,pid),1,nz; call goPr
         ! copy shape, needed to avoid strange compilation error ... 
         lshp = self%shp(:,pid)
         ! copy input slab into buffer:
             
             
             
!          sendbuf(displ+1:displ+n) = reshape( input(self%lbo(1,pid):self%ubo(1,pid),self%lbo(2,pid):self%ubo(2,pid)), (/n/) )
         
         sendbuf(displ+1:displ+n) = reshape( input((j-1)*nx+1:j*nx,(i-1)*ny+1:i*ny), (/n/) )
         i=i+1
         if (i>splitx) then
            i=1
            j=j+1
         end if
         
!          write (gol,'("Santiago sendbuf ...")'); call goPr
!          write (gol,*) '  pid  : ', pid; call goPr
!          write (gol,*) '  self%lbo(1,pid)  : ', self%lbo(1,pid); call goPr
!          write (gol,*) '  self%ubo(1,pid)  : ', self%ubo(1,pid); call goPr
!          write (gol,*) '  self%lbo(2,pid)  : ', self%lbo(2,pid); call goPr
!          write (gol,*) '  self%ubo(2,pid)  : ', self%ubo(2,pid); call goPr
       end if
       
       ! increase displacement:
       displ = displ + n
     end do ! iproc
     ! check ...
     if ( displ /= gn ) then
       write (gol,'("prepared ",i0," values in sendbuf, but local output size is ",i0)') displ, gn; call goErr
       TRACEBACK; status=1; return
     end if
   
     ! fill:
     !sendbuf = reshape( input, (/gn/) )
    
   else
   
     ! dummy storage ...
     allocate( sendbuf(1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( sendcounts(0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( displs    (0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     
   end if

   ! ~ 
   
#ifdef _MPI
   ! data type:
   call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
   IF_NOT_OK_RETURN(status=1)

   ! info ...
   write (gol,'("send data ...")'); call goPr
   write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
   write (gol,*) '   sendcounts  : ', sendcounts; call goPr
   write (gol,*) '   displs      : ', displs; call goPr
   write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
   write (gol,*) '   recvcount   : ', recvcount; call goPr
   ! transfer data:
   call MPI_ScatterV( sendbuf, sendcounts, displs, mpi_dtype, &
                      recvbuf, recvcount         , mpi_dtype, &
                      goc%root_id, goc%comm, status )
   IF_MPI_NOT_OK_RETURN(status=1)
#else
   !! info ...
   !write (gol,'("copy data ...")'); call goPr
   ! copy:
   recvbuf = sendbuf
#endif

   ! ~ unpack

   ! defined ?
   if ( n > 0 ) then
     !! testing ...
     !write (gol,'("fill ",i0," elements in output(",i0,":",i0,2(",",i0,":",i0),") from recvbuf")') &
     !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),1,nz; call goPr
     ! copy buffer into target array:
     output(self%lbo(1,goc%id):self%ubo(1,goc%id),&
            self%lbo(2,goc%id):self%ubo(2,goc%id)) = &
               reshape(recvbuf,(/self%shp(1,goc%id),self%shp(2,goc%id)/))
   end if ! n > 0

   ! clear:
   deallocate( sendbuf, stat=status )
   IF_NOT_OK_RETURN(status=1)
   deallocate( sendcounts, stat=status )
   IF_NOT_OK_RETURN(status=1)
   deallocate( displs, stat=status )
   IF_NOT_OK_RETURN(status=1)
   
   ! clear:
   deallocate( recvbuf )
   IF_NOT_OK_RETURN(status=1)
   
   ! ok
   status = 0
   
 end subroutine Domains_ScatterV_2d_r  



  subroutine Domains_ScatterV_3d_r( self, input, output, status )
  
    use GO_Comm, only : goc
#ifdef _MPI
    use MPI_F08, only : MPI_DataType
   use MPI_F08, only : MPI_ScatterV
#endif

   ! --- const ---------------------------------
   
   character(len=*), parameter   ::  rname = mname//'/Domains_ScatterV_3d_r'
   
   ! rank and kind:
   integer, parameter  ::  ndim = 2
 
   ! --- in/out ---------------------------------
   
   class(T_Domains), intent(in)  ::  self
   real, intent(in)              ::  input(:,:,:)
   real, intent(out)             ::  output(:,:,:)
   integer, intent(out)          ::  status
   
   ! --- local -----------------------------------
   
#ifdef _MPI
   type(MPI_DataType)          ::  mpi_dtype
#endif
   integer                     ::  nx, ny, nz
   integer                     ::  gn
   integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
   integer                     ::  lshp(ndim)
   real, allocatable           ::  sendbuf(:)     ! (input n)
   integer, allocatable        ::  sendcounts(:)  ! (nproc)
   integer, allocatable        ::  displs(:)      ! (nproc)
   real, allocatable           ::  recvbuf(:)     ! (output n)
   integer                     ::  recvcount
   integer                     ::  pid
   integer                     ::  n
   integer                     ::  displ
   
   ! --- begin -----------------------------------
   
   ! dims:
   nx = size(output,1)
   ny = size(output,2)
   nz = size(output,3)
   
   ! local output slab ?
   if ( all( self%shp(:,goc%id) > 0 ) ) then
     ! check ...
     if ( any( (/nx,ny/) /= self%shp(:,goc%id) ) ) then
       write (gol,'("shape of output (",i0,",",i0,",:) while expected (",i0,",",i0,",:)")') &
                     nx,ny, self%shp(:,goc%id); call goErr
       TRACEBACK; status=1; return
     end if
   !else
   !  ! info ...
   !  write (gol,'("empty output domain, careful ...")'); call goPr
   end if

   ! global bounds:
   glbo = self%glbo(:,0)
   gubo = self%gubo(:,0)
   do pid = 1, goc%npes-1
     if ( self%n(pid) > 0 ) then
       glbo = min( glbo, self%glbo(:,pid) )
       gubo = max( gubo, self%gubo(:,pid) )
     end if
   end do
   ! global shape:
   gshp = gubo - glbo + 1
   ! source?
   if ( goc%root ) then
     ! check ...
     if ( any( shape(input) /= (/gshp(1),gshp(2),nz/) ) ) then
       write (gol,'("shape of input (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                     shape(input), gshp,nz; call goErr
       TRACEBACK; status=1; return
     end if      
     ! total:
     gn = product( gshp ) * nz
     ! offset, arguments are 1-based:
     goff = glbo - 1
   else
     ! dummy value to avoid warnings ...
     goff = 0
   end if
   
   ! ~ create receive buffer
   
   ! local size:
   n = self%n(goc%id) * nz
   ! store:
   recvcount = n
   ! defined ?
   if ( n > 0 ) then
     ! storage:
     allocate( recvbuf(n), stat=status )
     IF_NOT_OK_RETURN(status=1)
   else
     ! dummy:
     allocate( recvbuf(1), stat=status )
     IF_NOT_OK_RETURN(status=1)
   end if

   ! ~ create send buffer
   
   ! target?
   if ( goc%root ) then
   
     ! storage:
     allocate( sendbuf(gn), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( sendcounts(0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( displs    (0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     
     ! init displacement:
     displ = 0
     ! loop over target processes:
     do pid = 0, goc%npes-1
       ! target has values ?
       if ( self%n(pid) > 0 ) then
         ! all values:
         n = self%n(pid) * nz
       else
         ! empty:
         n = 0
       end if
       ! number of values:
       sendcounts(pid) = n
       ! displacement:
       displs(pid) = displ
       
       ! any target data from this process ?
       if ( n > 0 ) then
         !! info ...
         !write (gol,'("store ",i0," elements in sendbuf(",i0,":",i0,") from input(",i0,":",i0,2(",",i0,":",i0),")")') &
         !        n, displ+1, displ+n, self%lbo(1,pid),self%ubo(1,pid),self%lbo(2,pid),self%ubo(2,pid),1,nz; call goPr
         ! copy shape, needed to avoid strange compilation error ... 
         lshp = self%shp(:,pid)
         ! copy input slab into buffer:
         sendbuf(displ+1:displ+n) = reshape( input(self%lbo(1,pid):self%ubo(1,pid),&
                                                   self%lbo(2,pid):self%ubo(2,pid),1:nz), &
                                             (/n/) )
       end if
       
       ! increase displacement:
       displ = displ + n
     end do ! iproc
     ! check ...
     if ( displ /= gn ) then
       write (gol,'("prepared ",i0," values in sendbuf, but local output size is ",i0)') displ, gn; call goErr
       TRACEBACK; status=1; return
     end if
   
     ! fill:
     sendbuf = reshape( input, (/gn/) )

   else
   
     ! dummy storage ...
     allocate( sendbuf(1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( sendcounts(0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     allocate( displs    (0:goc%npes-1), stat=status )
     IF_NOT_OK_RETURN(status=1)
     
   end if

   ! ~ 
   
#ifdef _MPI
   ! data type:
   call goc%GetDataType( 'real', kind(sendbuf), mpi_dtype, status )
   IF_NOT_OK_RETURN(status=1)

   ! info ...
   write (gol,'("send data ...")'); call goPr
   write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
   write (gol,*) '   sendcounts  : ', sendcounts; call goPr
   write (gol,*) '   displs      : ', displs; call goPr
   write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
   write (gol,*) '   recvcount   : ', recvcount; call goPr
   ! transfer data:
   call MPI_ScatterV( sendbuf, sendcounts, displs, mpi_dtype, &
                      recvbuf, recvcount         , mpi_dtype, &
                      goc%root_id, goc%comm, status )
   IF_MPI_NOT_OK_RETURN(status=1)
#else
   !! info ...
   !write (gol,'("copy data ...")'); call goPr
   ! copy:
   recvbuf = sendbuf
#endif

   ! ~ unpack

   ! defined ?
   if ( n > 0 ) then
     !! testing ...
     !write (gol,'("fill ",i0," elements in output(",i0,":",i0,2(",",i0,":",i0),") from recvbuf")') &
     !        n, self%lbo(1,goc%id),self%ubo(1,goc%id),self%lbo(2,goc%id),self%ubo(2,goc%id),1,nz; call goPr
     ! copy buffer into target array:
     output(self%lbo(1,goc%id):self%ubo(1,goc%id),&
            self%lbo(2,goc%id):self%ubo(2,goc%id),1:nz) = &
               reshape(recvbuf,(/self%shp(1,goc%id),self%shp(2,goc%id),nz/))
   end if ! n > 0

   ! clear:
   deallocate( sendbuf, stat=status )
   IF_NOT_OK_RETURN(status=1)
   deallocate( sendcounts, stat=status )
   IF_NOT_OK_RETURN(status=1)
   deallocate( displs, stat=status )
   IF_NOT_OK_RETURN(status=1)
   
   ! clear:
   deallocate( recvbuf )
   IF_NOT_OK_RETURN(status=1)
   
   ! ok
   status = 0
   
 end subroutine Domains_ScatterV_3d_r  






  ! ***************************************************************************
  ! ***
  ! *** extract local slab
  ! ***
  ! ***************************************************************************


  subroutine Domains_Extract_1d_r8( self, input, output, status )
  
    use GO_Comm, only : goc

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Extract_1d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:)  ! global array
    real(wpr), intent(out)        ::  output(:)  ! local array
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------

    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  pid
    
    ! --- begin -----------------------------------
    
    ! local slab has non-zero shape ?
    if ( all( self%shp(:,goc%id) > 0 ) ) then

      ! global bounds:
      glbo = self%glbo(:,0)
      gubo = self%gubo(:,0)
      do pid = 1, goc%npes-1
        if ( self%n(pid) > 0 ) then
          glbo = min( glbo, self%glbo(:,pid) )
          gubo = max( gubo, self%gubo(:,pid) )
        end if
      end do
      ! global shape:
      gshp = gubo - glbo + 1
      ! check ...
      if ( any( shape(input) /= gshp ) ) then
        write (gol,'("shape of input (",i0,") while expected global shape (",i0,")")') &
                      shape(input), gshp; call goErr
        TRACEBACK; status=1; return
      end if
      ! offset, arguments are 1-based:
      goff = glbo - 1

      ! check ...
      if ( any( shape(output) /= self%shp(:,goc%id) ) ) then
        write (gol,'("shape of output (",i0,") while expected local shape (",i0,")")') &
                      shape(output), self%shp(:,goc%id); call goErr
        TRACEBACK; status=1; return
      end if

      !! testing ..
      !print *, me, ': extract ', self%glbo(1,goc%id)-goff(1), self%gubo(1,goc%id)-goff(1)
      
      ! copy local slab:
      output = input(self%glbo(1,goc%id)-goff(1):self%gubo(1,goc%id)-goff(1))

    end if

    ! ok
    status = 0
    
  end subroutine Domains_Extract_1d_r8
  

  ! ***************************************************************************
  ! ***
  ! *** write field to nc file
  ! ***
  ! ***************************************************************************


  !
  ! Collect local 2D domains on root, put into netcdf file.
  ! Use irec<0 for variable without time dimension.
  !

  subroutine Domains_Put_Var_2d( self, ncid, varid, irec, values, status )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_Put_Var
#endif
    use GO_Comm, only : goc

    ! --- in/out ---------------------------------

    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  ncid
    integer, intent(in)           ::  varid
    integer, intent(in)           ::  irec
    real, intent(in)              ::  values(:,:)
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Domains_Put_Var_2d'

    ! --- local ----------------------------------
    
    integer               ::  shp(2)
    real, allocatable     ::  glb_values(:,:)

    ! --- begin ----------------------------------
    
    ! global shape:
    call self%Get( status, glb_shp=shp )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    if ( goc%root ) then
      allocate( glb_values(shp(1),shp(2)), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( glb_values(1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! collect:
    call self%GatherV( values, glb_values, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! root?
    if ( goc%root ) then
    
      ! without time dimension?
      if ( irec < 0 ) then

        ! write field:
#ifdef with_netcdf
        status = NF90_Put_Var( ncid, varid, glb_values )
        IF_NF90_NOT_OK_RETURN(status=1)
#endif

      else

        ! write slab:
#ifdef with_netcdf
        status = NF90_Put_Var( ncid, varid, glb_values, &
                               start=(/1,1,irec/), count=(/shp(1),shp(2),1/) )
        IF_NF90_NOT_OK_RETURN(status=1)
#endif

      end if ! time dim?
      
    end if  ! root

    ! storage:
    deallocate( glb_values, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_Put_Var_2d


  ! *
  

  !
  ! Collect local 2D domains on root, put into netcdf file.
  ! Use irec<0 for variable without time dimension.
  !

  subroutine Domains_Put_Var_3d( self, ncid, varid, irec, values, status )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_Put_Var
#endif
    use GO_Comm, only : goc

    ! --- in/out ---------------------------------

    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  ncid
    integer, intent(in)           ::  varid
    integer, intent(in)           ::  irec
    real, intent(in)              ::  values(:,:,:)
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Domains_Put_Var_3d'

    ! --- local ----------------------------------
    
    integer               ::  nz
    integer               ::  shp(2)
    real, allocatable     ::  glb_values(:,:,:)

    ! --- begin ----------------------------------
    
    ! global shape:
    call self%Get( status, glb_shp=shp )
    IF_NOT_OK_RETURN(status=1)

    ! other dims:
    nz = size(values,3)
    
    ! storage:
    if ( goc%root ) then
      allocate( glb_values(shp(1),shp(2),nz), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( glb_values(1,1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! collect:
    call self%GatherV( values, glb_values, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! root?
    if ( goc%root ) then
    
      ! without time dimension?
      if ( irec < 0 ) then

        ! write field:
#ifdef with_netcdf
        status = NF90_Put_Var( ncid, varid, glb_values )
        IF_NF90_NOT_OK_RETURN(status=1)
#endif

      else

        ! write slab:
#ifdef with_netcdf
        status = NF90_Put_Var( ncid, varid, glb_values, &
                               start=(/1,1,1,irec/), count=(/shp(1),shp(2),nz,1/) )
        IF_NF90_NOT_OK_RETURN(status=1)
#endif

      end if ! time dim?
      
    end if  ! root

    ! storage:
    deallocate( glb_values, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_Put_Var_3d


  !
  ! Collect local 2D domains on root, put into netcdf file.
  ! Use irec<0 for variable without time dimension.
  !

  subroutine Domains_Put_Var_4d( self, ncid, varid, irec, values, status )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_Put_Var
#endif
    use GO_Comm, only : goc

    ! --- in/out ---------------------------------

    class(T_Domains), intent(in)  ::  self
    integer, intent(in)           ::  ncid
    integer, intent(in)           ::  varid
    integer, intent(in)           ::  irec
    real, intent(in)              ::  values(:,:,:,:)
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Domains_Put_Var_4d'

    ! --- local ----------------------------------
    
    integer               ::  nz, ns
    integer               ::  is
    integer               ::  shp(2)
    real, allocatable     ::  glb_values(:,:,:)

    ! --- begin ----------------------------------
    
    ! global shape:
    call self%Get( status, glb_shp=shp )
    IF_NOT_OK_RETURN(status=1)

    ! other dims:
    nz = size(values,3)
    ns = size(values,4)
    
    ! storage:
    if ( goc%root ) then
      allocate( glb_values(shp(1),shp(2),nz), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( glb_values(1,1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! loop over fourth dimension
    do is = 1, ns
    
      ! collect:
      call self%GatherV( values(:,:,:,is), glb_values, status )
      IF_NOT_OK_RETURN(status=1)

      ! root?
      if ( goc%root ) then

        ! without time dimension?
        if ( irec < 0 ) then

          ! write field:
#ifdef with_netcdf
          status = NF90_Put_Var( ncid, varid, glb_values , &
                                 start=(/1,1,1,is/), count=(/shp(1),shp(2),nz,1/))
          IF_NF90_NOT_OK_RETURN(status=1)
#endif

        else

          ! write slab:
#ifdef with_netcdf
          status = NF90_Put_Var( ncid, varid, glb_values, &
                                 start=(/1,1,1,is,irec/), count=(/shp(1),shp(2),nz,1,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
#endif

        end if ! time dim?

      end if  ! root
    
    end do  ! fourth dimension (nspec?)
    
    ! storage:
    deallocate( glb_values, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Domains_Put_Var_4d


end module GO_Domains
