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
!     use GO_Comm, only : TgoComm                               use GO_Comm, only : TgoComm
!
!     use LE_Dims, only : nx, ny                                  use OPS_Dims, only : nx, ny
!     use LE_Data, only : uwind   ! (nx,ny)                       use OPS_Data, only : conc   ! (nx,ny)
!
!     integer            ::  status                               integer            ::  status
!     class(TgoComm)    ::  goc                                  class(TgoComm)    ::  goc
!     integer            ::  tag                                  integer            ::  tag
!     integer            ::  ops_nx, ops_ny                       integer            ::  le_nx, le_ny
!     real, allocatable  ::  ops_conc(:,:)                        real, allocatable  ::  le_uwind(:,:)
!
!     call GO_Comm_Init( goc, 'LE', status )                      call GO_Comm_Init( goc, 'OPS', status )
!
!     ! send own dimensions:                                      ! receive dimensions of other application:
!     tag = 101                                                   tag = 101
!     call GO_Comm_Send( goc, 'OPS', nx, status )                 call GO_Comm_Recv( goc, 'LE', le_nx, status )
!     call GO_Comm_Send( goc, 'OPS', ny, status )                 call GO_Comm_Recv( goc, 'LE', le_ny, status )
!
!     ! recive dimensions of other application:                   ! send own dimensions:
!     tag = 101                                                   tag = 101
!     call GO_Comm_Recv( goc, 'OPS', ops_nx, status )             call GO_Comm_Send( goc, 'LE', nx, status )
!     call GO_Comm_Recv( goc, 'OPS', ops_ny, status )             call GO_Comm_Send( goc, 'LE', ny, status )
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
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; errorcode=status; call MPI_Error_String(errorcode,gol,ngol,status); call goErr; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#######################################################################

module GO_Comm

  use GO_Print, only : gol, ngol, goPr, goErr
#ifdef _MPI
  use MPI_F08, only : MPI_SUCCESS, MPI_Error_String
  use MPI_F08, only : MPI_Comm
  use MPI_F08, only : MPI_Info
#endif

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  goc


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'GO_Comm'

  ! character lengths:
  integer, parameter  ::  LEN_NAME = 32


  ! --- types-------------------------------------

  !! application info:
  !type T_Appl
  !  ! application name:
  !  character(len=LEN_NAME)         ::  name
  !end type T_Appl

  ! comunicator:
  type TgoComm
    !! own name:
    !character(len=LEN_NAME)         ::  name
#ifdef _MPI
    ! MPI communicator;
    ! use 'comm%mpi_val' for the integer:
    type(MPI_Comm)                  ::  comm
    ! MPI info structure;
    ! use 'info%mpi_val' for the integer:
    type(MPI_Info)                  ::  info
#endif
    ! size and rank:
    integer                         ::  npes
    integer                         ::  id
    ! root?
    integer                         ::  root_id
    logical                         ::  root
    !! info for each of the applications, including own:
    !type(T_Appl), allocatable       ::  appl(:)   ! (0:npes-1)
    !
  contains
    procedure   ::  Init        =>  GO_Comm_Init
    procedure   ::  Done        =>  GO_Comm_Done
    !procedure   ::  SetName     =>  GO_Comm_SetName
    !procedure   ::  GetID       =>  GO_Comm_GetID
    procedure   ::  Abort       =>  GO_Comm_Abort
    procedure   ::  Barrier     =>  GO_Comm_Barrier
    !
#ifdef _MPI
    procedure   ::  GetDataType =>  GO_Comm_GetDataType
    procedure   ::  GetOper     =>  GO_Comm_GetOper
#endif
    !
    procedure                       GO_Comm_SendAndRecv_i4_2d
    procedure                       GO_Comm_SendAndRecv_r4_2d
    generic     ::  SendAndRecv =>  GO_Comm_SendAndRecv_i4_2d, &
                                    GO_Comm_SendAndRecv_r4_2d
    !
    procedure                       GO_Comm_Reduce_i
    generic     ::  Reduce      =>  GO_Comm_Reduce_i
    !
    procedure                       GO_Comm_AllReduce_i
    procedure                       GO_Comm_AllReduce_r
    procedure                       GO_Comm_AllReduce_l
    procedure                       GO_Comm_AllReduce_r4_1d
    procedure                       GO_Comm_AllReduce_r4_2d
    procedure                       GO_Comm_AllReduce_r4_3d
    procedure                       GO_Comm_AllReduce_r4_4d
    procedure                       GO_Comm_AllReduce_r4_5d
    procedure                       GO_Comm_AllReduce_InPlace_r4_5d
    generic     ::  AllReduce   =>  GO_Comm_AllReduce_i, &
                                    GO_Comm_AllReduce_r, &
                                    GO_Comm_AllReduce_l, &
                                    GO_Comm_AllReduce_r4_1d, &
                                    GO_Comm_AllReduce_r4_2d, &
                                    GO_Comm_AllReduce_r4_3d, &
                                    GO_Comm_AllReduce_r4_4d, &
                                    GO_Comm_AllReduce_r4_5d, &
                                    GO_Comm_AllReduce_InPlace_r4_5d
    !
    procedure ::                    GO_Comm_BCast_r4_1d
    procedure ::                    GO_Comm_BCast_r4_2d
    procedure ::                    GO_Comm_BCast_r4_3d
    procedure ::                    GO_Comm_BCast_r4_4d
    generic   :: BCast          =>  GO_Comm_BCast_r4_1d, &
                                    GO_Comm_BCast_r4_2d, &
                                    GO_Comm_BCast_r4_3d, &
                                    GO_Comm_BCast_r4_4d
    procedure ::  Gather        =>  GO_Comm_Gather_i
    !
    procedure ::                    GO_Comm_AllGather_i
    procedure ::                    GO_Comm_AllGather_r
    procedure ::                    GO_Comm_AllGather_i1
    procedure ::                    GO_Comm_AllGather_r1
    generic   ::  AllGather     =>  GO_Comm_AllGather_i, &
                                    GO_Comm_AllGather_r, &
                                    GO_Comm_AllGather_i1, &
                                    GO_Comm_AllGather_r1
    !
    procedure                       GO_Comm_GatherV_i_1d
    procedure                       GO_Comm_GatherV_r4_1d
    procedure                       GO_Comm_GatherV_r4_2d
    procedure                       GO_Comm_GatherV_r4_3d
    generic     ::  GatherV     =>  GO_Comm_GatherV_i_1d, &
                                    GO_Comm_GatherV_r4_1d, &
                                    GO_Comm_GatherV_r4_2d, &
                                    GO_Comm_GatherV_r4_3d
    !
    procedure ::  ParInfo       =>  GO_Comm_ParInfo
    !
  end type TgoComm


  ! --- var -------------------------------------

  ! mpi error code; used as argument for 'MPI_Error_String' to avoid
  ! warnings about same argument 'status' being used for both 'errorcode' and 'ierror':
  integer                           ::  errorcode

  ! communicator:
  type(TgoComm)       ::  goc



contains


  ! ====================================================================


  subroutine GO_Comm_Init( self, status )

#ifdef _MPI
    use MPI_F08, only : MPI_COMM_WORLD
    use MPI_F08, only : MPI_INFO_NULL
    use MPI_F08, only : MPI_Init
    use MPI_F08, only : MPI_Comm_Size, MPI_Comm_Rank
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(out)       ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Init'

    ! --- local ----------------------------------

    integer         ::  id

    ! --- begin ----------------------------------

#ifdef _MPI
    ! init MPI interface:
    call MPI_Init( ierror=status )
    if ( status /= MPI_SUCCESS ) then
      write (gol,'("could not initialize MPI environment; exit code : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! use global communicator:
    self%comm = MPI_COMM_WORLD

    ! dummy info:
    self%info = MPI_INFO_NULL

    ! size and rank:
    call MPI_Comm_Size( self%comm, self%npes, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    call MPI_Comm_Rank( self%comm, self%id, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe:
    goc%npes = 1
    goc%id = 0

#endif

    !! storage for info on applications:
    !allocate( self%appl(0:self%npes-1), stat=status )
    !IF_NOT_OK_RETURN(status=1)

    !! empty name:
    !self%name = ''

    ! set root id:
    self%root_id = 0
    ! root?
    self%root = self%id == self%root_id

    ! ok
    status = 0

  end subroutine GO_Comm_Init


  ! ***


  subroutine GO_Comm_Done( self, status )

#ifdef _MPI
    use MPI_F08, only : MPI_Finalize
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! done with MPI interface:
    call MPI_Finalize( ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    !! clear:
    !deallocate( self%appl, stat=status )
    !IF_NOT_OK_RETURN(status=1)

    !! reset name:
    !self%name = ''

    ! ok
    status = 0

  end subroutine GO_Comm_Done


!  ! ***
!
!
!  subroutine GO_Comm_SetName( self, name, status )
!
!#ifdef _MPI
!    use MPI_F08, only : MPI_COMM_WORLD
!    use MPI_F08, only : MPI_CHARACTER
!#endif
!
!    ! --- in/out ---------------------------------
!
!    class(TgoComm), intent(out)       ::  self
!    character(len=*), intent(in)      ::  name
!    integer, intent(out)              ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/GO_Comm_SetName'
!
!    ! --- local ----------------------------------
!
!    integer         ::  id
!
!    ! --- begin ----------------------------------
!
!    ! store name:
!    self%name = trim(name)
!
!    ! store own:
!    self%appl(self%id)%name = trim(name)
!#ifdef _MPI
!    ! exchange names; loop over all pes:
!    do id = 0, self%npes-1
!      ! send from this id to all other:
!      call MPI_BCast( self%appl(id)%name, len(self%appl(id)%name), MPI_CHARACTER, &
!                          id, self%comm, status )
!      IF_MPI_NOT_OK_RETURN(status=1)
!    end do
!#endif
!
!    ! info ...
!    write (gol,'(a," - connected to the following pes:")') trim(self%name); call goPr
!    do id = 0, self%npes-1
!      if ( id == self%id ) cycle
!      write (gol,'(a," -   pe ",i3," with application `",a,"`")') &
!               trim(self%name), id, trim(self%appl(id)%name); call goPr
!    end do
!
!    ! ok
!    status = 0
!
!  end subroutine GO_Comm_SetName
!
!
!  ! ***
!
!
!  ! Return process ID of requested application.
!
!  subroutine GO_Comm_GetID( self, name, id, status )
!
!    ! --- in/out ---------------------------------
!
!    class(TgoComm), intent(inout)       ::  self
!    character(len=*), intent(in)        ::  name
!    integer, intent(out)                ::  id
!    integer, intent(out)                ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GetID'
!
!    ! --- local ----------------------------------
!
!    integer           ::  k
!
!    ! --- begin ----------------------------------
!
!    ! init as dummy ..
!    id = -999
!    ! loop over applications:
!    do k = 0, self%npes-1
!      ! compare:
!      if ( trim(name) == trim(self%appl(k)%name) ) then
!        ! found:
!        id = k
!        ! leave:
!        exit
!      end if
!    end do   ! pe's
!    ! not found ?
!    if ( id < 0 ) then
!      write (gol,'("could not find application with name `",a,"`")') trim(name); call goErr
!      write (gol,'("available applications:")'); call goErr
!      do k = 0, self%npes-1
!        write (gol,'("  - ",a)') trim(self%appl(k)%name); call goErr
!      end do
!      TRACEBACK; status=1; return
!    end if
!
!    ! ok
!    status = 0
!
!  end subroutine GO_Comm_GetID


  ! ***


  subroutine GO_Comm_Barrier( self, status )

#ifdef _MPI
    use MPI_F08, only : MPI_Barrier
#endif

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Barrier'

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    integer, intent(out)                ::  status

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

#ifdef _MPI
    ! start:
    call MPI_Barrier( self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_Barrier


  ! ***


  subroutine GO_Comm_Abort( self, errcode )

#ifdef _MPI
    use MPI_F08, only : MPI_Abort
#endif

    ! --- in/out ----------------------------

    class(TgoComm), intent(inout)       ::  self
    integer, intent(in)                 ::  errcode

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/GO_Comm_Abort'

    ! --- local -----------------------------

    integer        ::  status

    ! --- begin -----------------------------

#ifdef _MPI
    ! abort, get status back since even this might fail ...
    call MPI_Abort( goc%comm, errcode, ierror=status )
    if ( status /= MPI_SUCCESS ) then
      errorcode=status
      call MPI_Error_String( errorcode, gol, ngol, status ); call goErr
      TRACEBACK
      ! alternative abort via system exit (non-standard routine!);
      ! this might not kill all other processes:
      call Exit( status )
    end if
#else
    ! system exit (non-standard routine!)
    call Exit( errcode )
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_Abort


  ! ====================================================================


#ifdef _MPI

  !
  ! Return 'MPI_DataType' for specified type ('real', 'integer') and kind.
  ! Note that the result does not use the communicatator,
  ! but since this object is needed for each MPI call anayway it is an
  ! easy way to ship it ...
  !

  subroutine GO_Comm_GetDataType( self, typ, knd, dtype, status )

    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_LOGICAL
    use MPI_F08, only : MPI_INTEGER
    use MPI_F08, only : MPI_REAL, MPI_DOUBLE_PRECISION

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    character(len=*), intent(in)        ::  typ
    integer, intent(in)                 ::  knd
    type(MPI_DataType), intent(out)     ::  dtype
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GetDataType'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! switch:
    select case ( trim(typ) )

      !~ logical:
      case ( 'logical' )
        ! set result:
        dtype = MPI_LOGICAL

      !~ integer:
      case ( 'integer' )
        ! switch:
        select case ( knd )
          case ( 4 )
            dtype = MPI_INTEGER
          case default
            write (gol,'("could not set MPI data type for ",a," kind ",i0)') trim(typ), knd; call goErr
            TRACEBACK; status=1; return
        end select

      !~ reals:
      case ( 'real' )
        ! switch:
        select case ( knd )
          case ( 4 )
            dtype = MPI_REAL
          case ( 8 )
            dtype = MPI_DOUBLE_PRECISION
          case default
            write (gol,'("could not set MPI data type for ",a," kind ",i0)') trim(typ), knd; call goErr
            TRACEBACK; status=1; return
        end select

      !~
      case default
        write (gol,'("could not set MPI data type for ",a," variables")') trim(typ); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine GO_Comm_GetDataType


  ! *


  !
  ! Translate collective operator description ('min','max')
  ! to an MPI_Op type.
  !

  subroutine GO_Comm_GetOper( self, oper, op, status )

    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_MIN, MPI_MAX, MPI_SUM
    use MPI_F08, only : MPI_LAND, MPI_LOR

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    character(len=*), intent(in)        ::  oper
    type(MPI_Op), intent(out)           ::  op
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GetOper'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! operator:
    select case ( oper )
      case ( 'min' ) ; op = MPI_MIN
      case ( 'max' ) ; op = MPI_MAX
      case ( 'sum' ) ; op = MPI_SUM
      case ( 'and' ) ; op = MPI_LAND
      case ( 'or'  ) ; op = MPI_LOR
      case default
        write (gol,'("unsupported oper `",a,"`")') trim(oper); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine GO_Comm_GetOper

#endif


  ! ********************************************************************
  ! ***
  ! *** send and receive
  ! ***
  ! ********************************************************************


  subroutine GO_Comm_SendAndRecv_i4_2d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_INTEGER
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    integer, intent(in)                 ::  values(:,:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_SendAndRecv_i4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_Status)  ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! sending or receiving pe?
    if ( goc%id == source ) then

      !! testing ...
      !write (gol,'("  x    send ",i0," integers (tag ",i0,")")') size(values), tag; call goPr

      ! send from this (source) pe to destination pe:
      call MPI_Send( values, size(values), MPI_DTYPE, &
                       dest, tag, self%comm, &
                       ierror=status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( goc%id == dest ) then

      !! testing ...
      !write (gol,'("  x    recv ",i0," integers (tag ",i0,")")') size(values), tag; call goPr

      ! receive on this (destination) pe from source pe:
      call MPI_Recv( values, size(values), MPI_DTYPE, &
                       source, tag, self%comm, &
                       mpi_stat, ierror=status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    end if
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_SendAndRecv_i4_2d
  
  
  ! *


  subroutine GO_Comm_SendAndRecv_r4_2d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    real(4), intent(in)                 ::  values(:,:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_SendAndRecv_r4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_Status)  ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! sending or receiving pe?
    if ( goc%id == source ) then

      !! testing ...
      !write (gol,'("  x    send ",i0," reals (tag ",i0,")")') size(values), tag; call goPr

      ! send from this (source) pe to destination pe:
      call MPI_Send( values, size(values), MPI_DTYPE, &
                       dest, tag, self%comm, &
                       status )
      IF_MPI_NOT_OK_RETURN(status=1)

    else if ( goc%id == dest ) then

      !! testing ...
      !write (gol,'("  x    recv ",i0," reals (tag ",i0,")")') size(values), tag; call goPr

      ! receive on this (destination) pe from source pe:
      call MPI_Recv( values, size(values), MPI_DTYPE, &
                       source, tag, self%comm, &
                       mpi_stat, ierror=status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
    end if
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_SendAndRecv_r4_2d


  ! ********************************************************************
  ! ***
  ! *** reduce
  ! ***
  ! ********************************************************************


  !
  ! reduce values to root;
  ! 'oper' is one of: 'sum' | 'min' | 'max'
  !
  
  subroutine GO_Comm_Reduce_i( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_Reduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    character(len=*), intent(in)        ::  oper
    integer, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Reduce_i'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    integer               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call goc%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_Reduce( value_loc, value, 1, dtype, op, self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_Reduce_i


  ! ********************************************************************
  ! ***
  ! *** reduce, and send result to all
  ! ***
  ! ********************************************************************


  subroutine GO_Comm_AllReduce_i( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    character(len=*), intent(in)        ::  oper
    integer, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_i'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    integer               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call goc%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllReduce_i

  ! *

  subroutine GO_Comm_AllReduce_r( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    character(len=*), intent(in)        ::  oper
    real, intent(inout)                 ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_r'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    real                  ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call goc%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllReduce_r

  ! *

  subroutine GO_Comm_AllReduce_l( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
    character(len=*), intent(in)        ::  oper
    logical, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllReduce_l'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    logical               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call goc%GetDataType( 'logical', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call goc%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllReduce_l


  ! .......................................................................
  ! ...
  ! ... allreduce (from lekf v3.0.000, should become similar as above ...
  ! ...
  ! .......................................................................


  subroutine GO_Comm_AllReduce_r4_1d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

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
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

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
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

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
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

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
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

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
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_SUM
    use MPI_F08, only : MPI_IN_PLACE
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(inout)       ::  self
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
    IF_MPI_NOT_OK_RETURN(status=1)

#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllReduce_InPlace_r4_5d


  ! ********************************************************************
  ! ***
  ! *** broadcast
  ! ***
  ! ********************************************************************


  subroutine GO_Comm_BCast_r4_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_1d'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), MPI_REAL, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_BCast_r4_1d


  ! ***


  subroutine GO_Comm_BCast_r4_2d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), MPI_REAL, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_BCast_r4_2d


  ! ***


  subroutine GO_Comm_BCast_r4_3d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), MPI_REAL, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_BCast_r4_3d


  ! ***


  subroutine GO_Comm_BCast_r4_4d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_REAL
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), MPI_REAL, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_BCast_r4_4d


  ! ********************************************************************
  ! ***
  ! *** gather
  ! ***
  ! ********************************************************************

  !
  ! Fill on root an array with length npes
  ! with values collected from everywhere.
  !

  subroutine GO_Comm_Gather_i( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_INTEGER
    use MPI_F08, only : MPI_Gather
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    integer, intent(in)               ::  value
    integer, intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_Gather_i'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! collect from all pe and broadcast result:
    call MPI_Gather( value, 1, MPI_DTYPE, &
                     values, 1, MPI_DTYPE, &
                     self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_Gather_i

  ! *


  ! ********************************************************************
  ! ***
  ! *** gather and broadcast
  ! ***
  ! ********************************************************************

  !
  ! Fill on each processor an array with length npes
  ! with values broadcasted from everywhere.
  !

  subroutine GO_Comm_AllGather_i( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_INTEGER
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    integer, intent(in)               ::  value
    integer, intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllGather_i'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( size(values) /= self%npes ) then
      write (gol,'("output array has size ",i0," while npes=",i0)') size(values), self%npes; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, 1, MPI_DTYPE, &
                        values, 1, MPI_DTYPE, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllGather_i

  ! *

  subroutine GO_Comm_AllGather_r( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    real(4), intent(in)               ::  value
    real(4), intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllGather_r'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( size(values) /= self%npes ) then
      write (gol,'("output array has size ",i0," while npes=",i0)') size(values), self%npes; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, 1, MPI_DTYPE, &
                        values, 1, MPI_DTYPE, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllGather_r

  ! *

  subroutine GO_Comm_AllGather_i1( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_INTEGER
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    integer, intent(in)               ::  value(:)     ! (n)
    integer, intent(out)              ::  values(:,:)  ! (n,1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllGather_i1'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( (size(values,1) /= size(value)) .or. (size(values,2) /= self%npes) ) then
      write (gol,'("output array has shape (",i0,",",i0,") while input has shape (",i0,") and npes=",i0)') &
          size(value), shape(values), self%npes; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, size(value), MPI_DTYPE, &
                        values, size(value), MPI_DTYPE, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values(:,1) = value(:)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllGather_i1

  ! *

  subroutine GO_Comm_AllGather_r1( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    real(4), intent(in)               ::  value(:)     ! (n)
    real(4), intent(out)              ::  values(:,:)  ! (n,1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_AllGather_r1'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! check ...
    if ( (size(values,1) /= size(value)) .or. (size(values,2) /= self%npes) ) then
      write (gol,'("output array has shape (",i0,",",i0,") while input has shape (",i0,") and npes=",i0)') &
          size(value), shape(values), self%npes; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, size(value), MPI_DTYPE, &
                        values, size(value), MPI_DTYPE, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values(:,1) = value(:)
#endif

    ! ok
    status = 0

  end subroutine GO_Comm_AllGather_r1



  ! ********************************************************************
  ! ***
  ! *** GatherV
  ! ***
  ! ********************************************************************
  
  !
  ! Collect send buffers in 1D array on root.
  ! If send is supposed to be empty, use optional nloc=0 to specify this.
  !

  subroutine GO_Comm_GatherV_i_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_INTEGER
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    integer, intent(in)                 ::  send(:)   ! (max(1,nloc))
    integer, intent(out)                ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GatherV_i_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( goc%root ) then
      if ( size(recv) /= ntot ) then
        write (gol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, MPI_DTYPE, &
                      recv, recvcounts, displs, MPI_DTYPE, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine GO_Comm_GatherV_i_1d
  
  ! *

  subroutine GO_Comm_GatherV_r4_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    real(4), intent(in)                 ::  send(:)   ! (max(1,nloc))
    real(4), intent(out)                ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GatherV_r4_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( goc%root ) then
      if ( size(recv) /= ntot ) then
        write (gol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, MPI_DTYPE, &
                      recv, recvcounts, displs, MPI_DTYPE, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine GO_Comm_GatherV_r4_1d
  
  ! *

  subroutine GO_Comm_GatherV_r4_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    real(4), intent(in)                 ::  send(:,:)   ! (m,max(1,nloc))
    real(4), intent(out)                ::  recv(:,:)   ! (m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GatherV_r4_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m
    integer                 ::  n
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dim:
    m = size(send,1)
    
#ifdef _MPI

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( goc%root ) then
      ! check ...
      if ( any( shape(recv) /= (/m,ntot/) ) ) then
        write (gol,'("receive buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(recv), m,ntot; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, m*n, MPI_DTYPE, &
                      recv, m*recvcounts, m*displs, MPI_DTYPE, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine GO_Comm_GatherV_r4_2d
  
  ! *

  subroutine GO_Comm_GatherV_r4_3d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DTYPE => MPI_REAL
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)          ::  self
    real(4), intent(in)                 ::  send(:,:,:)   ! (p,m,max(1,nloc))
    real(4), intent(out)                ::  recv(:,:,:)   ! (p,m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_GatherV_r4_2d'

    ! --- local ----------------------------------
    
    integer                 ::  p, m
    integer                 ::  n
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dims:
    p = size(send,1)
    m = size(send,2)
    
#ifdef _MPI

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( goc%root ) then
      ! check ...
      if ( any( shape(recv) /= (/p,m,ntot/) ) ) then
        write (gol,'("receive buffer has shape (",i0,2(",",i0),") while (m,ntot) is (",i0,2(",",i0),")")') &
                        shape(recv), p,m,ntot; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, p*m*n, MPI_DTYPE, &
                      recv, p*m*recvcounts, p*m*displs, MPI_DTYPE, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine GO_Comm_GatherV_r4_3d



  ! ********************************************************************
  ! ***
  ! *** i/o tools
  ! ***
  ! ********************************************************************

  !
  ! Data with different size n per processor.
  ! Returns:
  ! - total size (output dimension) 
  ! - start index in output array (netcdf put_var argument)
  !

  subroutine GO_Comm_ParInfo( self, n, status, &
                                ntot, istart, recvcounts, displs, displ )

    ! --- in/out ---------------------------------

    class(TgoComm), intent(in)        ::  self
    integer, intent(in)               ::  n
    integer, intent(out)              ::  status
    
    integer, intent(out), optional    ::  ntot
    integer, intent(out), optional    ::  istart
    integer, intent(out), optional    ::  recvcounts(:)  ! (npes)
    integer, intent(out), optional    ::  displs(:)  ! (npes)
    integer, intent(out), optional    ::  displ

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GO_Comm_ParInfo'

    ! --- local ----------------------------------
    
    integer, allocatable    ::  ns(:)
    integer                 ::  i

    ! --- begin ----------------------------------
    
    ! storage:
    allocate( ns(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of data values for each processor, broadcast to all:
    call self%AllGather( n, ns, status )
    IF_NOT_OK_RETURN(status=1)

    ! total number of data values:
    if ( present(ntot) ) ntot = sum(ns)

    ! start in output array:
    if ( present(istart) ) then
      if ( self%id == 0 ) then
        istart = 1
      else
        istart = sum(ns(0:self%id-1)) + 1
      end if      
    end if
    
    ! receive counts:
    if ( present(recvcounts) ) recvcounts = ns

    ! displacements:
    if ( present(displs) ) then
      do i = 1, self%npes
        if ( i == 1 ) then
          displs(i) = 0
        else
          displs(i) = sum(ns(0:i-2))
        end if
      end do
    end if
    ! current displacement:
    if ( present(displ) ) then
      if ( self%id == 0 ) then
        displ = 0
      else
        displ = sum(ns(0:self%id-1))
      end if
    end if

    ! clear:
    deallocate( ns, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine GO_Comm_ParInfo



end module GO_Comm

