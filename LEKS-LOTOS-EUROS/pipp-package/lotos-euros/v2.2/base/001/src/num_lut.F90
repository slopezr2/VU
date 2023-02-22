!#################################################################
!
! NAME
!   Num_LUT - Look-Up-Table
!
! DESCRIPTION
!
!   Changes compared to GO_LUT :
!    - axes are not expected to be linear increasing 
!      but assumed to be irregular spaced;
!    - selection is based on linear interpolation
!    - interpolation indices and weights are stored internally,
!      this allows multiple applications to different tables
!      with the same shape
!
!   
! USAGE
!
!   ! object to interpolate in table:
!   type(T_LookUp)              ::  lookup
!    
!   ! init 3D lookup:
!   call lookup%Init( 3, status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! define axis by passing 3 1D arrays:
!   call lookup%SetAx( 1, (/1.0,2.0,3.0,../), status )  
!   IF_NOTOK_RETURN(status=1)
!   call lookup%SetAx( 2, (/23.0,38.0,42.0,.../), status )  
!   IF_NOTOK_RETURN(status=1)
!   call lookup%SetAx( 3, (/0.0,100.0,350.0,.../), status )  
!   IF_NOTOK_RETURN(status=1)
!
!   ! set internally stored indices and weights for linear interpolation
!   ! (3D in table) ; constant values if outside bounds:
!   call lookup%InterpolSet( (/2.2,20.0,75.0/), status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! apply to one or more arrays:
!   call lookup%InterpolApply( data1(:,:,:), value1, status )
!   IF_NOTOK_RETURN(status=1)
!   call lookup%InterpolApply( data2(:,:,:), value2, status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! done with lookup:
!   call lookup%Done( status )
!   IF_NOTOK_RETURN(status=1)
!
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#################################################################


module Num_LUT

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public    ::  T_LookUp_1D, T_LookUp
  

  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  mname = 'Num_LUT'
  
    
  ! --- types ------------------------------------
  
  ! NOTE: parameterized types are in F2003 standard,
  ! but not supported by compilers yet ...  
  
  type :: T_LookUp_1D!( knd )
    !! parameters:
    !integer, kind                 ::  knd
    !! axis values:
    !real(kind=knd), allocatable   ::  x(:)
    ! axis values:
    integer                       ::  n
    real, allocatable             ::  x(:)
    ! speedup search:
    integer                       ::  ilast
  contains
    procedure   ::  Init            => LookUp_1D_Init
    procedure   ::  Done            => LookUp_1D_Done
    procedure   ::  InterpolWeights => LookUp_1D_InterpolWeights
  end type T_LookUp_1D
    
  ! *
  
  type T_LookUp!( ndim, knd )
    ! parameters:
    !integer, len                     ::  ndim
    !integer, kind                    ::  knd
    !integer                          ::  ndim
    ! axis:
    !type(T_LookUp_1D(knd))           ::  ax(ndim)
    !     logical                     ::  defined(ndim)
    ! parameters:
    integer                          ::  ndim
    ! axis:
    type(T_LookUp_1D), allocatable   ::  ax(:)  ! (ndim)
    logical, allocatable             ::  defined(:)  ! (ndim)
    integer, allocatable             ::  shp(:)
    ! interpolation:
    integer                          ::  np
    integer, allocatable             ::  ii(:,:)    ! (ndim,np)
    real, allocatable                ::  ww(:)      ! (np)
    logical                          ::  setup
  contains
    procedure   ::  Init            => LookUp_Init
    procedure   ::  Done            => LookUp_Done
    procedure   ::  SetAx           => LookUp_SetAx
    procedure   ::  InterpolSet     => LookUp_InterpolSet_3d
    procedure   ::  InterpolApply   => LookUp_InterpolApply_3d
  end type T_LookUp


contains


  ! ====================================================================
  ! ===
  ! === looking up in table
  ! ===
  ! ====================================================================


  subroutine LookUp_1D_Init( self, x, status )     

    ! --- in/out ---------------------------------
    
    class(T_LookUp_1D), intent(out)       ::  self
    real, intent(in)                      ::  x(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_1D_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! size:
    self%n = size(x)
    ! storage:
    allocate( self%x(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! copy:
    self%x = x
    
    ! init search index:
    self%ilast = 1

    ! ok
    status = 0

  end subroutine LookUp_1D_Init
  
  
  ! ***


  subroutine LookUp_1D_Done( self, status )     

    ! --- in/out ---------------------------------
    
    class(T_LookUp_1D), intent(inout)       ::  self
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_1D_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%x, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LookUp_1D_Done
  
  
  ! ***


  subroutine LookUp_1D_InterpolWeights( self, x0, ii, ww, status )
  
    use Num_Tools, only : Interval    

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_1D_InterpolWeights'
    
    ! --- in/out ---------------------------------
    
    class(T_LookUp_1D), intent(inout)       ::  self
    real, intent(in)                        ::  x0
    integer, intent(out)                    ::  ii(2)
    real, intent(out)                       ::  ww(2)
    integer, intent(out)                    ::  status

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! interval index:
    call Interval( self%x, x0, self%ilast, status )
    if ( status < 0 ) then
      ! left from first element:
      ii(1) = 1 ; ww(1) = 1.0
      ii(2) = 1 ; ww(2) = 0.0
    else if ( status > 0 ) then
      ! right from last element:
      ii(1) = self%n ; ww(1) = 1.0
      ii(2) = self%n ; ww(2) = 0.0
    else
      ! linear weights:
      ii(1) = self%ilast   ; ww(1) = (self%x(self%ilast+1) - x0) / (self%x(self%ilast+1) - self%x(self%ilast))
      ii(2) = self%ilast+1 ; ww(2) = 1.0 - ww(1)
    end if
    
    ! check ..
    if ( any(ww < 0.0) .or. any(ww > 1.0) ) then
      write (gol,'("wrong interpol weights: ",2f12.4)'); call goErr
      write (gol,'("  x0 : ",es12.4)') x0; call goErr
      write (gol,'("  x  : ")'); call goErr
      write (gol,*) self%x; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine LookUp_1D_InterpolWeights


  ! ====================================================================
  ! ===
  ! === looking up in nD table
  ! ===
  ! ====================================================================

  subroutine LookUp_Init( self, ndim, status )     

    ! --- in/out ---------------------------------
    
    class(T_LookUp), intent(out)          ::  self
    integer, intent(in)                   ::  ndim
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! store:
    self%ndim = ndim
    
    ! storage:
    allocate( self%ax(self%ndim), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%defined(self%ndim), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%shp(self%ndim), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! init flags:
    self%defined = .false.
    
    ! number of points in ND array involveld with linear interpolation:
    self%np = 2**self%ndim
    ! storage:
    allocate( self%ii(self%ndim,self%np), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%ww(self%np), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! not filled yet:
    self%setup = .false.

    ! ok
    status = 0

  end subroutine LookUp_Init
  
  
  ! ***


  subroutine LookUp_Done( self, status )     

    ! --- in/out ---------------------------------
    
    class(T_LookUp), intent(inout)          ::  self
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_Done'

    ! --- local ----------------------------------
    
    integer   ::  idim

    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%ii, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%ww, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over axes:
    do idim = 1, self%ndim
      ! enabled ?
      if ( self%defined(idim) ) then
        ! clear:
        call self%ax(idim)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if ! enabled
    end do  ! dims

    ! clear:
    deallocate( self%ax, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%defined, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%shp, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LookUp_Done
  
  
  ! ***


  subroutine LookUp_SetAx( self, idim, x, status )     

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_SetAx'
    
    ! --- in/out ---------------------------------
    
    class(T_LookUp), intent(inout)          ::  self
    integer, intent(in)                     ::  idim
    real, intent(in)                        ::  x(:)
    integer, intent(out)                    ::  status

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( (idim < 1) .or. (self%ndim < idim) ) then
      write (gol,'("requested dim ",i0," outside range 1:",i0)') idim, self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( self%defined(idim) ) then
      write (gol,'("ax ",i0," already defined")') idim; call goErr
      TRACEBACK; status=1; return
    end if
    ! init:
    call self%ax(idim)%Init( x, status )
    IF_NOTOK_RETURN(status=1)
    ! set flag:
    self%defined(idim) = .true.
    ! store size:
    self%shp(idim) = size(x)
    
    ! ok
    status = 0
    
  end subroutine LookUp_SetAx
  
  
  ! ***


  subroutine LookUp_InterpolSet_3d( self, x0, status )
  
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_InterpolSet_3d'
    
    integer, parameter  ::  nd  = 3

    ! --- in/out ---------------------------------
    
    class(T_LookUp), intent(inout)          ::  self
    real, intent(in)                        ::  x0(nd)
    integer, intent(out)                    ::  status

    ! --- local ----------------------------------

    integer                  ::  idim
    integer                  ::  ip
    integer                  ::  ii(2,nd)
    real                     ::  ww(2,nd)
    integer                  ::  i, j, k
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. all(self%defined) ) then
      write (gol,'("not all dimensions defined : ",7l2)') self%defined; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( nd /= self%ndim ) then
      write (gol,'("input vector x0 has size ",i0," while ndim is ",i0)') nd, self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop over dims:
    do idim = 1, nd
      ! get indices and weights for linear interpolation within this dimension:
      call self%ax(idim)%InterpolWeights( x0(idim), ii(:,idim), ww(:,idim), status )
      IF_NOTOK_RETURN(status=1)
    end do  ! dims
    
    ! init counter:
    ip = 0
    ! loop:
    do i = 1, 2
      do j = 1, 2
        do k = 1, 2
          ! next:
          ip = ip + 1
          ! fill indices in array:
          self%ii(1,ip) = ii(i,1)
          self%ii(2,ip) = ii(j,2)
          self%ii(3,ip) = ii(k,3)
          ! weights:
          self%ww(ip) = ww(i,1) * ww(j,2) * ww(k,3)
        end do
      end do
    end do
    ! (re)set flag:
    self%setup = .true.
    
    ! ok
    status = 0
    
  end subroutine LookUp_InterpolSet_3d
  
  
  ! ***


  subroutine LookUp_InterpolApply_3d( self, data, value, status )
  
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LookUp_InterpolApply_3d'
    
    ! --- in/out ---------------------------------
    
    class(T_LookUp), intent(in)         ::  self
    real, intent(in)                    ::  data(:,:,:)
    real, intent(out)                   ::  value
    integer, intent(out)                ::  status

    ! --- local ----------------------------------

    integer     ::  ip
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%setup ) then
      write (gol,'("interpolation weights not setup yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( any( self%shp /= shape(data) ) ) then
      write (gol,'("shape of data table (",i0,2(",",i0),") does not match with definition (",i0,2(",",i0),")")') &
                      self%shp, shape(data); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init weighted sum:
    value = 0.0
    ! loop over interpolation points:
    do ip = 1, self%np
      ! add contribution for this point:
      value = value + data(self%ii(1,ip),self%ii(2,ip),self%ii(3,ip)) * self%ww(ip)
    end do
    
    ! ok
    status = 0
    
  end subroutine LookUp_InterpolApply_3d


end module Num_LUT

