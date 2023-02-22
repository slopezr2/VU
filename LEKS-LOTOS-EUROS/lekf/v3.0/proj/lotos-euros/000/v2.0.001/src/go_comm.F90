!#######################################################################
!
! GO_Comm - communication through MPI.
!
!
! USAGE
!
!   Exchange of data between two executables, both using this module.
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    model 1 "LE"                                                model 2 "OPS"
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     use GO_Comm, only : T_GO_Comm                               use GO_Comm, only : T_GO_Comm
!
!     use LE_Dims, only : nx, ny                                  use OPS_Dims, only : nx, ny
!     use LE_Data, only : uwind   ! (nx,ny)                       use OPS_Data, only : conc   ! (nx,ny)
!
!     integer            ::  status                               integer            ::  status
!     type(T_GO_Comm)    ::  goc                                  type(T_GO_Comm)    ::  goc
!     integer            ::  tag                                  integer            ::  tag
!     integer            ::  ops_nx, ops_ny                       integer            ::  le_nx, le_ny
!     real, allocatable  ::  ops_conc(:,:)                        real, allocatable  ::  le_uwind(:,:)
!
!     call GO_Comm_Init( goc, 'LE', status )                      call GO_Comm_Init( goc, 'OPS', status )
!
!     ! send own dimensions:                                      ! receive dimensions of other application:
!     tag = 101                                                   tag = 101
!     call GO_Comm_Send( goc, 'OPS', tag, nx, status )            call GO_Comm_Recv( goc, 'LE', tag, le_nx, status )
!     call GO_Comm_Send( goc, 'OPS', tag, ny, status )            call GO_Comm_Recv( goc, 'LE', tag, le_ny, status )
!
!     ! recive dimensions of other application:                   ! send own dimensions:
!     tag = 101                                                   tag = 101
!     call GO_Comm_Recv( goc, 'OPS', tag, ops_nx, status )        call GO_Comm_Send( goc, 'LE', tag, nx, status )
!     call GO_Comm_Recv( goc, 'OPS', tag, ops_ny, status )        call GO_Comm_Send( goc, 'LE', tag, ny, status )
!
!     ! storage:                                                  ! storage:
!     allocate( ops_conc(ops_nx,ops_ny) )                         allocate( le_uwind(le_nx,le_ny) )
!
!     ! send field:                                               ! receive field:
!     tag = 123                                                   tag = 123
!     call GO_Comm_Send( goc, 'OPS', tag, uwind, status )         call GO_Comm_Recv( goc, 'LE', tag, le_uwind, status )
!
!     ! receive field:                                            ! send field:
!     tag = 456                                                   tag = 456
!     call GO_Comm_Recv( goc, 'OPS', tag, ops_conc, status)       call GO_Comm_Send( goc, 'OPS', tag, conc, status )
!
!     call GO_Comm_Done( goc, status )                            call GO_Comm_Done( goc, status )
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! MACRO'S
!
!   _MPI       : should be defined to enable MPI library
!
!
! SEE ALSO
!   MPI Tutorial:
!     https://computing.llnl.gov/tutorials/mpi/
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!!
#define IF_MPI_NOT_OK_ABORT(comm,exitcode) if (status/=MPI_SUCCESS) then; TRACEBACK; call MPI_Abort(comm,exitcode,status); end if
!
#include "go.inc"
!
!#######################################################################

module GO_Comm

  use GO_Print, only : gol, goPr, goErr
#ifdef _MPI
  use MPI, only : MPI_SUCCESS, MPI_Abort
#endif

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_GO_Comm

  public  ::  goc
  
  public  ::  GO_Comm_Init, GO_Comm_Done
  public  ::  GO_Comm_GetID
  public  ::  GO_Comm_Send, GO_Comm_Recv



  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  mname = 'GO_Comm'
  
  ! character lengths:
  integer, parameter  ::  LEN_NAME = 32


  ! --- types-------------------------------------
  
  ! application info:
  type T_Appl
    ! application name:
    character(len=LEN_NAME)         ::  name
  end type T_Appl
  
  ! comunicator:
  type T_GO_Comm
    ! own name:
    character(len=LEN_NAME)         ::  myname
    ! MPI communicator:
    integer                         ::  comm
    ! size and rank:
    integer                         ::  npes
    integer                         ::  myid
    ! info for each of the applications, including own:
    type(T_Appl), allocatable       ::  appl(:)   ! (0:npes-1)
    ! flags:
    logical                         ::  root
    !
  contains
    procedure ::  Init              =>  GO_Comm_Init
    procedure ::  Done              =>  GO_Comm_Done
    !
    procedure ::                        GO_Comm_AllReduce_r4_1d
    procedure ::                        GO_Comm_AllReduce_r4_2d
    procedure ::                        GO_Comm_AllReduce_r4_3d
    procedure ::                        GO_Comm_AllReduce_r4_4d
    procedure ::                        GO_Comm_AllReduce_r4_5d
    procedure ::                        GO_Comm_AllReduce_InPlace_r4_5d
    generic   ::  AllReduce         =>  GO_Comm_AllReduce_r4_1d, &
                                        GO_Comm_AllReduce_r4_2d, &
                                        GO_Comm_AllReduce_r4_3d, &
                                        GO_Comm_AllReduce_r4_4d, &
                                        GO_Comm_AllReduce_r4_5d, &
                                        GO_Comm_AllReduce_InPlace_r4_5d
    !
    procedure ::                        GO_Comm_BCast_r4_1d
    procedure ::                        GO_Comm_BCast_r4_2d
    procedure ::                        GO_Comm_BCast_r4_3d
    procedure ::                        GO_Comm_BCast_r4_4d
    generic   :: BCast              =>  GO_Comm_BCast_r4_1d, &
                                        GO_Comm_BCast_r4_2d, &
                                        GO_Comm_BCast_r4_3d, &
                                        GO_Comm_BCast_r4_4d
    !
  end type T_GO_Comm


  ! --- interfaces -------------------------------
  
  interface GO_Comm_Recv
    module procedure GO_Comm_Recv_i4_0d
    module procedure GO_Comm_Recv_r4_0d
    module procedure GO_Comm_Recv_r4_2d
    module procedure GO_Comm_Recv_r4_3d
  end interface
  
  interface GO_Comm_Send
    module procedure GO_Comm_Send_i4_0d
    module procedure GO_Comm_Send_r4_0d
    module procedure GO_Comm_Send_r4_2d
    module procedure GO_Comm_Send_r4_3d
  end interface


  ! --- variables -------------------------------
  
  type(T_GO_Comm)   ::  goc


contains


  ! ====================================================================


  subroutine GO_Comm_Init( self, name, status )

#ifdef _MPI
    use MPI, only : MPI_COMM_WORLD
    use MPI, only : MPI_CHARACTER
    use MPI, only : MPI_Init
    use MPI, only : MPI_Comm_Size, MPI_Comm_Rank
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(out)     ::  self
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Init'
    
    ! --- local ----------------------------------
    
    integer           ::  id
    integer           ::  l
    character(len=2)  ::  cc
    
    ! --- begin ----------------------------------

#ifdef _MPI
    ! init MPI interface:
    call MPI_Init( status )
    if ( status /= 0 ) then
      write (gol,'("could not initialize MPI environment; exit code : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! use global communicator:
    self%comm = MPI_COMM_WORLD

    ! size and rank:
    call MPI_Comm_Size( self%comm, self%npes, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    call MPI_Comm_Rank( self%comm, self%myid, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)

#else
    ! serial dummy values:
    self%comm = 0
    self%npes = 1
    self%myid = 0
#endif
    
    ! set flag:
    self%root = self%myid == 0
    
    ! store name:
    self%myname = trim(name)
    ! eventually replace '##' by id number:
    l = index( self%myname, '##' )
    if ( l > 0 ) then
      write (cc,'(i2.2)') self%myid
      self%myname(l:l+2) = cc
    end if
    
    ! storage for info on applications:
    allocate( self%appl(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! store own:
    self%appl(self%myid)%name = trim(self%myname)
#ifdef _MPI
    ! exchange names; loop over all pes:
    do id = 0, self%npes-1
      ! send from id all other:
      call MPI_Bcast( self%appl(id)%name, len(self%appl(id)%name), MPI_CHARACTER, &
                          id, self%comm, status )
      IF_MPI_NOT_OK_ABORT(self%comm,1)
    end do
#endif
    
    ! info ...
    write (gol,'(a," - connected to the following pes:")') trim(self%myname); call goPr
    do id = 0, self%npes-1
      if ( id == self%myid ) cycle
      write (gol,'(a," -   pe ",i3," with application `",a,"`")') &
               trim(self%myname), id, trim(self%appl(id)%name); call goPr
    end do
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Init
  
  
  ! ***


  subroutine GO_Comm_Done( self, status )

    !use MPI, only : MPI_Finalize    ! not availble from module yet
  
    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

#ifdef _MPI  
    ! done with MPI interface:
    call MPI_Finalize( status )
    IF_NOT_OK_RETURN(status=1)
#endif
    
    ! clear:
    deallocate( self%appl, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! reset name:
    self%myname = ''
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Done
  
  
  ! ***
  
  
  ! Return process ID of requested application.
  
  subroutine GO_Comm_GetID( self, name, id, status )

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  name
    integer, intent(out)                ::  id
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GetID'
    
    ! --- local ----------------------------------
    
    integer           ::  k
    
    ! --- begin ----------------------------------

    ! init as dummy ..
    id = -999
    ! loop over applications:
    do k = 0, self%npes-1
      ! compare:
      if ( trim(name) == trim(self%appl(k)%name) ) then
        ! found:
        id = k
        ! leave:
        exit
      end if
    end do   ! pe's
    ! not found ?
    if ( id < 0 ) then
      write (gol,'("could not find application with name `",a,"`")') trim(name); call goErr
      write (gol,'("available applications:")'); call goErr
      do k = 0, self%npes-1
        write (gol,'("  - ",a)') trim(self%appl(k)%name); call goErr
      end do
      TRACEBACK; status=1; return
    end if    

    ! reset name:
    self%myname = ''
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_GetID


  ! ***************************************************************************
  ! ***
  ! *** tools
  ! ***
  ! ***************************************************************************
  

  subroutine GO_Comm_Barrier( self, status )
  
#ifdef _MPI
    use MPI, only : MPI_Barrier
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Barrier'
    
    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    integer, intent(out)                ::  status
    
    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
#ifdef _MPI
    ! start:
    call MPI_Barrier( self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif

    ! ok
    status = 0
    
  end subroutine GO_Comm_Barrier
  
  
  ! ***
  
  
  subroutine GO_Comm_Abort( self, ierr, status )

#ifdef _MPI
    use MPI      , only : MPI_Abort
#endif
    use GO_System, only : goExit
  
    ! --- in/out ----------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    integer, intent(in)                 ::  ierr
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Comm_Abort'          
    
    ! --- begin -----------------------------
    
#ifdef _MPI
    ! start:
    call MPI_Abort( self%comm, ierr, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#else
    ! system exit (non-standard routine!)
    call goExit( 1 )
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Abort


  ! ********************************************************************
  ! ***
  ! *** send
  ! ***
  ! ********************************************************************
  
    
  ! Send field to other application ;
  ! at the other side, a call to 'GO_Comm_Recv' should be present.
  
  subroutine GO_Comm_Send_i4_0d( self, appl, tag, value, status )

#ifdef _MPI
    use MPI, only : MPI_INTEGER4
    !use MPI, only : MPI_Send
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    integer(4), intent(in)              ::  value
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Send_i4_0d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
    
    ! --- begin ----------------------------------

    ! target id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! send ...
    call MPI_Send( value, 1, MPI_INTEGER4, id, tag, self%comm, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Send_i4_0d
     
  ! ***
    
  subroutine GO_Comm_Send_r4_0d( self, appl, tag, value, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_Send
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(in)                 ::  value
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Send_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
    
    ! --- begin ----------------------------------

    ! target id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! send ...
    call MPI_Send( value, 1, MPI_REAL, id, tag, self%comm, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Send_r4_0d  
  
! ***

  
  subroutine GO_Comm_Send_r4_2d( self, appl, tag, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_Send
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(in)                 ::  values(:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Send_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
    
    ! --- begin ----------------------------------

    ! target id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! send ...
    call MPI_Send( values, size(values), MPI_REAL, id, tag, self%comm, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Send_r4_2d

! ***
  
  subroutine GO_Comm_Send_r4_3d( self, appl, tag, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_Send
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(in)                 ::  values(:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Send_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
    
    ! --- begin ----------------------------------

    ! target id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! send ...
    call MPI_Send( values, size(values), MPI_REAL, id, tag, self%comm, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Send_r4_3d


  ! ********************************************************************
  ! ***
  ! *** broadcast
  ! ***
  ! ********************************************************************

  
  subroutine GO_Comm_BCast_r4_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_1d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! send values from root to all other pe's:    
    call MPI_BCast( values, size(values), MPI_REAL, rootid, goc%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_BCast_r4_1d
  
  
  ! ***

  
  subroutine GO_Comm_BCast_r4_2d( self, rootid, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! send values from root to all other pe's:    
    call MPI_BCast( values, size(values), MPI_REAL, rootid, goc%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_BCast_r4_2d
  
  
  ! ***

  
  subroutine GO_Comm_BCast_r4_3d( self, rootid, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! send values from root to all other pe's:    
    call MPI_BCast( values, size(values), MPI_REAL, rootid, goc%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_BCast_r4_3d
  
  
  ! ***

  
  subroutine GO_Comm_BCast_r4_4d( self, rootid, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    !use MPI, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! send values from root to all other pe's:    
    call MPI_BCast( values, size(values), MPI_REAL, rootid, goc%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_BCast_r4_4d



  ! ********************************************************************
  ! ***
  ! *** receive
  ! ***
  ! ********************************************************************
  

  ! Receive field from other application ;
  ! at the other side, a call to 'GO_Comm_Send' should be present.
  
  subroutine GO_Comm_Recv_i4_0d( self, appl, tag, value, status )

#ifdef _MPI
    use MPI, only : MPI_INTEGER4
    use MPI, only : MPI_STATUS_SIZE
    !use MPI, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    integer(4), intent(out)             ::  value
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Recv_i4_0d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
#ifdef _MPI
    integer           ::  mpi_status(MPI_STATUS_SIZE)
#endif
    
    ! --- begin ----------------------------------

    ! source id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! Recv ...
    call MPI_Recv( value, 1, MPI_INTEGER4, id, tag, self%comm, mpi_status, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#else
    ! dummy ...
    value = 0
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Recv_i4_0d
  
  ! ***
      
  subroutine GO_Comm_Recv_r4_0d( self, appl, tag, value, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_STATUS_SIZE
    !use MPI, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(out)                ::  value
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Recv_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
#ifdef _MPI
    integer           ::  mpi_status(MPI_STATUS_SIZE)
#endif
    
    ! --- begin ----------------------------------

    ! source id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! Recv ...
    call MPI_Recv( value, 1, MPI_REAL, id, tag, self%comm, mpi_status, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#else
    ! dummy ...
    value = 0.0
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Recv_r4_0d

  ! ***
    
  subroutine GO_Comm_Recv_r4_2d( self, appl, tag, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_STATUS_SIZE
    !use MPI, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(out)                ::  values(:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Recv_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
#ifdef _MPI
    integer           ::  mpi_status(MPI_STATUS_SIZE)
#endif
    
    ! --- begin ----------------------------------

    ! source id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! Recv ...
    call MPI_Recv( values, size(values), MPI_REAL, id, tag, self%comm, mpi_status, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#else
    ! dummy ...
    values = 0.0
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Recv_r4_2d

  ! ***
    
  subroutine GO_Comm_Recv_r4_3d( self, appl, tag, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_STATUS_SIZE
    !use MPI, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------
    
    type(T_GO_Comm), intent(inout)      ::  self
    character(len=*), intent(in)        ::  appl
    integer, intent(in)                 ::  tag
    real(4), intent(out)                ::  values(:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Recv_r4_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  id
#ifdef _MPI
    integer           ::  mpi_status(MPI_STATUS_SIZE)
#endif
    
    ! --- begin ----------------------------------

    ! source id:
    call GO_Comm_GetID( self, appl, id, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef _MPI
    ! Recv ...
    call MPI_Recv( values, size(values), MPI_REAL, id, tag, self%comm, mpi_status, status ) 
    IF_MPI_NOT_OK_ABORT(self%comm,1)
#else
    ! dummy ...
    values = 0.0
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_Recv_r4_3d


  ! ********************************************************************
  ! ***
  ! *** allreduce
  ! ***
  ! ********************************************************************

  
  subroutine GO_Comm_AllReduce_r4_1d( self, input, output, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(in)                 ::  input(:)
    real(4), intent(out)                ::  output(:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r4_1d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum inputs from all pe, store in output on each pe:    
    call MPI_AllReduce( input, output, size(input), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#else

    ! single pe, just copy:
    output = input

#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_r4_1d
  

  ! ***

  
  subroutine GO_Comm_AllReduce_r4_2d( self, input, output, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(in)                 ::  input(:,:)
    real(4), intent(out)                ::  output(:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r4_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum inputs from all pe, store in output on each pe:    
    call MPI_AllReduce( input, output, size(input), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#else

    ! single pe, just copy:
    output = input

#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_r4_2d
  

  ! ***

  
  subroutine GO_Comm_AllReduce_r4_3d( self, input, output, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(in)                 ::  input(:,:,:)
    real(4), intent(out)                ::  output(:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r4_3d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum inputs from all pe, store in output on each pe:    
    call MPI_AllReduce( input, output, size(input), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#else

    ! single pe, just copy:
    output = input

#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_r4_3d
  

  ! ***
  

  subroutine GO_Comm_AllReduce_r4_4d( self, input, output, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(in)                 ::  input(:,:,:,:)
    real(4), intent(out)                ::  output(:,:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r4_4d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum inputs from all pe, store in output on each pe:    
    call MPI_AllReduce( input, output, size(input), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#else

    ! single pe, just copy:
    output = input

#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_r4_4d
  

  ! ***
  

  subroutine GO_Comm_AllReduce_r4_5d( self, input, output, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(in)                 ::  input(:,:,:,:,:)
    real(4), intent(out)                ::  output(:,:,:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r4_5d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum inputs from all pe, store in output on each pe:    
    call MPI_AllReduce( input, output, size(input), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#else

    ! single pe, just copy:
    output = input

#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_r4_5d
  

  ! ***
  
  
  !
  ! Reduce "in place" (input and output buffer are the same)
  !

  subroutine GO_Comm_AllReduce_InPlace_r4_5d( self, values, status )

#ifdef _MPI
    use MPI, only : MPI_REAL
    use MPI, only : MPI_SUM
    use MPI, only : MPI_IN_PLACE
    !use MPI, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------
    
    class(T_GO_Comm), intent(inout)     ::  self
    real(4), intent(inout)              ::  values(:,:,:,:,:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_InPlace_r4_5d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! sum output from all pe, return in output on each pe;
    ! use special parameter 'MPI_IN_PLACE' for the send buffer,
    ! receive buffer is the input/output array:    
    call MPI_AllReduce( MPI_IN_PLACE, values, size(values), MPI_REAL, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_ABORT(self%comm,1)
    
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Comm_AllReduce_InPlace_r4_5d


!  ! ********************************************************************
!  ! ***
!  ! *** AllToAllV
!  ! ***
!  ! ********************************************************************
!  
!
!  ! ***
!  
!  
!  !
!  ! All-to-All with variable amounts of data.
!  ! On entry, the last dimension of values is assumed to be scattered,
!  ! thus filled with useful values on different process ;
!  ! the 'start' and 'count' specify for the last dimension
!  ! the first value that is filled and the number of values.
!  ! The shape of values should be the same on every processor.
!  ! On leaving, all data is exchanged with all other processes.
!  ! This is an "in place" version where input and output buffer are the same array.
!  !
!
!  subroutine GO_Comm_AllToAllV_r4_5d( self, values, start, count, status )
!
!#ifdef _MPI
!    use MPI, only : MPI_INT, MPI_REAL
!    use MPI, only : MPI_IN_PLACE
!    !use MPI, only : MPI_AllToAll
!    !use MPI, only : MPI_AllToAllV
!#endif
!
!    ! --- in/out ---------------------------------
!    
!    class(T_GO_Comm), intent(inout)     ::  self
!    real(4), intent(inout)              ::  values(:,:,:,:,:)
!    integer, intent(in)                 ::  start
!    integer, intent(in)                 ::  count
!    integer, intent(out)                ::  status
!    
!    ! --- const ----------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllToAllV_r4_5d'
!    
!    integer, parameter  ::  rank = 5
!    
!    ! --- local ----------------------------------
!    
!    integer                   ::  shp(rank)
!    integer                   ::  n
!    integer, allocatable      ::  sendcounts(:)
!    integer, allocatable      ::  sdispls(:)
!    integer, allocatable      ::  recvcounts(:)
!    integer, allocatable      ::  rdispls(:)
!    integer                   ::  i, j
!
!    ! --- begin ----------------------------------
!    
!#ifdef _MPI
!
!    ! shape of data:
!    shp = shape(values)
!    ! number of data values in first dimensions:
!    n = product(shp(1:rank-1))
!    
!    ! storage:
!    allocate( sendcounts(self%npes), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( sdispls(self%npes), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( recvcounts(self%npes), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( rdispls(self%npes), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    ! dummy because in-place call is used:
!    sendcounts = -999
!    sdispls    = -999
!    ! init arrays, will be exchanged:
!    recvcounts = 0
!    rdispls    = 0
!    ! local values:
!    recvcounts(self%myid) = n * count
!    rdispls   (self%myid) = n * (start-1)
!    ! exchange:
!    call MPI_AllToAll( MPI_IN_PLACE, 1, MPI_INT,
!                       recvcounts  , 1, MPI_INT, &
!                       self%comm, status )
!    IF_MPI_NOT_OK_ABORT(self%comm,1)
!    call MPI_AllToAll( MPI_IN_PLACE, 1, MPI_INT,
!                       rdispls     , 1, MPI_INT, &
!                       self%comm, status )
!    IF_MPI_NOT_OK_ABORT(self%comm,1)
!    ! check ...
!    if ( sum(recvcounts) /= size(values) ) then
!      write (gol,'("amount of received data does not match with array shape:")'); call goErr
!      write (gol,'("  shape              : ",i0,4(",",i0)') shape(values); call goErr
!      write (gol,'("  size of first dims : ",i0)') n; call goErr
!      write (gol,'("  received (in last dim):")'); call goErr
!      do i = 1, self%npes
!        write (gol,'("    pe ",i4," : ",i0," (",i0,")")') i-1, recvcounts(i), recvcounts(i)/n; call goErr
!      end do
!      TRACEBACK; status=1; return
!    end if
!    ! check ...
!    do i = 2, self%npes
!      if ( rdispls(i) /= rdispls(i-1)+recvcounts(i) ) then
!        write (gol,'("displacements and received counts do not match for pe ",i0,":")') i-1; call goErr
!        do j = 1, self%npes
!          write (gol,'("  pe ",i4," displ and count : ",i0," ",i0)') j-1, rdispls(j), recvcounts(j); call goErr
!        end do
!        TRACEBACK; status=1; return
!      end if
!    end do
!    
!    ! send local data, amount and location definded by recvcounts(myid) and rdispls(myid) ;
!    ! receive data from all other processes:
!    call MPI_AllToAllV( MPI_IN_PLACE, sendcounts, sdispls, MPI_REAL,
!                        values      , recvcounts, rdispls, MPI_REAL, &
!                        self%comm, status )
!    IF_MPI_NOT_OK_ABORT(self%comm,1)
!    
!    ! clear:
!    deallocate( sendcounts, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( sdispls, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( recvcounts, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( rdispls, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    
!#endif
!    
!    ! ok
!    status = 0
!    
!  end subroutine GO_Comm_AllToAllV_r4_5d


end module GO_Comm

