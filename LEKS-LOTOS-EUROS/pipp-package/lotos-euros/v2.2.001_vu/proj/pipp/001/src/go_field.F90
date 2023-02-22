!#######################################################################
!
! GO_Fields - storage and operation on 3D fields
!
!
! DESCRIPTION
!
!   Provides classes with data fields and operators to store
!   and setup 3D data fields.
!
!
! INHERITENCE
!
!   T_FieldBase
!
!     T_Field
!       T_Constant_Field
!         T_Constant_Field_Series
!         T_Accumulated_Field_Series
!       T_Instant_Field
!         T_Instant_Field_Series
!
!     T_VectorField
!       T_Instant_VectorField
!         T_Instant_VectorField_Series
!
!
! USAGE FOR CONSTANT FIELDS
!
!   ! A constant field remains valid for a certain time range,
!   ! which might be from now to infinity. If the field is not
!   ! valid anymore, a new value should be read.
!
!   use GO_Field, only : T_Constant_Field_Series
!   use GO_Date , only : TDate
!
!   ! dimensions:
!   integer, parameter     ::  nx = 4, ny = 3, nz = 1
!
!   ! data:
!   class(T_Constant_Field_Series)  ::  F
!   type(TDate)                     ::  t
!   type(TDate)                     ::  tt(2)
!   real                            ::  data(nx,ny,nz)
!   integer                         ::  status
!
!   ! Initialize time series of 2D field.
!   ! Define the content by a name and units.
!   call F%Init( 'orography', 'm', status )
!   if (status/=0) stop
!
!   ! current time, F%data should became valid for this:
!   t = NewDate( 2011, 03, 31, 18, 58 )
!
!   ! setup for time t; return status -1 if new data should be put:
!   call F%Setup( t, status )
!   if ( status > 0 ) then
!     stop 'error'
!   else if ( status < 0 ) then
!     ! obtain data from somewhere, also interval for which it is valid:
!     call read_input_data( ..., t, ..., data, tt, status )
!     if (status/=0) stop
!     ! info ...
!     call wrtgol( '    put field valid for ', tt(1), ' to ', tt(2) ); call goPr
!     ! store; optionally define the lower and upper bounds.
!     call F%Setup_Put( data, tt, status, &
!                          lbo=(/1,1,1/), ubo=(/nx,ny,nz/) )
!     if (status/=0) stop
!   end if
!
!   ! clear:
!   call F%Done( status )
!   if (status/=0) stop
!
!
!
! USAGE FOR TIMESERIES OF INSTANT FIELDS
!
!   !
!   ! Instant fields are valid for a single moment.
!   !
!   ! A series of fields could be defined which are valid
!   ! for any requested time, constructed from interpolation
!   ! in time between time steps for which samples are available.
!   ! The temporal interpolation should be defined first,
!   ! and based on this the setup methods will request the
!   ! appropriate samples to perform the interpolation.
!   ! If possible, samples that have been read before are re-used.
!   !
!
!   use GO_Field
!
!   ! dimensions:
!   integer, parameter     ::  nx = 4, ny = 3, nz = 2
!
!   ! data:
!   class(T_Instant_Field_Series)   ::  F
!   type(TDate)                     ::  t
!   integer                         ::  nreceive
!   integer                         ::  ireceive
!   real                            ::  receive_t
!   real                            ::  data(nx,ny,nz)
!   integer                         ::  status
!
!   ! Initialize time series of 3D field.
!   ! Define the content by a name and units.
!   call F%Init( 'temperature', 'K', status )
!   if (status/=0) stop
!
!   ! Define the temporal interpolation.
!   ! Currently only linear interpolation is supported:
!   !
!   !   interpolation=linear
!   !     Linear interpolation requires a step and units ;
!   !     start depends on units but is probably midnight:
!   !       step=3;units=hour
!   ! 
!   call F%Set_Interp( '3 hour', status )
!   if (status/=0) stop
!
!   ! current time, F%data should became valid for this:
!   t = NewDate( 2011, 03, 31, 18, 58 )
!
!   ! setup for time t; return status -1 if new data should be put:
!   call F%Setup( t, status )
!   if ( status > 0 ) then
!     stop 'error'
!   else if ( status < 0 ) then
!     ! prepare to receive new data:
!     call F%Setup_Prepare( nreceive, status )
!     if (status/=0) stop
!     ! loop over fields to be received:
!     do ireceive = 1, nreceive
!       ! get time value:
!       call F%Setup_InqTime( ireceive, receive_t, status )
!       if (status/=0) stop
!       ! obtain data from somewhere:
!       call read_input_data( ..., receive_t, ..., data, status )
!       if (status/=0) stop
!       ! info ...
!       call wrtgol( '    put field valid for ', receive_t ); call goPr
!       ! store; optionally define the lower and upper bounds.
!       call F%Setup_Put( ireceive, data, receive_t, status, &
!                              lbo=(/1,1,1/), ubo=(/nx,ny,nz/) )
!       if (status/=0) stop
!     end do
!   end if
!
!   ! clear:
!   call F%Done( status )
!   if (status/=0) stop
!
!
! TOOLS
!
!   For each of the types the following tools are available.
!
!   ! allocate external work array with same index range
!   ! as the data in the ojbect:
!   ...
!   real, allocatable   ::  tmp(:,:,:)
!   ...
!   call F%AllocExt( tmp, status )
!   if (status/=0) stop
!   ...
!   deallocate( tmp )
!
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#######################################################################


module GO_Field

  use GO_Print, only : gol, goPr, goErr
  use GO_Date , only : TDate

  implicit none
  
  
  ! ---- in/out ----------------------------------
  
  private
  
  public  ::  T_FieldBase

  public  ::  T_Field
  public  ::  T_Constant_Field
  public  ::  T_Constant_Field_Series
  public  ::  T_Accumulated_Field_Series
  
  public  ::  T_Instant_Field
  public  ::  T_Instant_Field_Series

  public  ::  T_VectorField
  public  ::  T_Instant_VectorField
  public  ::  T_Instant_VectorField_Series

  
  ! ---- const -----------------------------------
  
  ! module name:
  character(len=*), parameter        ::  mname = 'GO_Fields'
  
  ! character lengths:
  integer, parameter      ::  LEN_NAME  = 512
  integer, parameter      ::  LEN_UNITS = 64
  integer, parameter      ::  LEN_INPUT = 1024
  
  ! key values:
  integer, parameter      ::  INSTF_DATA1 = 1
  integer, parameter      ::  INSTF_DATA2 = 2
  integer, parameter      ::  ACCUF_CURR  = 3
  integer, parameter      ::  ACCUF_PREV  = 4
  
  
  ! --- types ------------------------------------

  type T_FieldBase
    ! description:
    character(len=LEN_NAME)   ::  name
    character(len=LEN_UNITS)  ::  units
    ! input description:
    character(len=LEN_INPUT)  ::  input
    ! data shape and bounds:
    integer                   ::  shp(3)
    integer                   ::  lbo(3), ubo(3)
    ! in use in this run ?
    logical                   ::  enabled
    ! filled ?
    logical                   ::  with_data
  contains
    procedure ::                  FieldBase_Init
    procedure ::                  FieldBase_Done
    procedure :: Check         => FieldBase_Check
    procedure :: Set           => FieldBase_Set
    procedure ::                  FieldBase_Get
    procedure ::                  FieldBase_SetupShape
  end type T_FieldBase
  
  ! *

  type, extends(T_FieldBase) :: T_Field
    ! data, pointers will refer to this array:
    real, pointer             ::  data(:,:,:)
  contains
    procedure :: Init          => Field_Init
    procedure :: Done          => Field_Done
    procedure :: Get           => Field_Get
    procedure :: Alloc         => Field_Alloc
    procedure :: AllocExt      => Field_AllocExt
    procedure :: Put           => Field_Put
  end type T_Field
  
  ! *
  
  type, extends(T_Field) :: T_Constant_Field
    ! timerange for which field is valid:
    type(TDate)               ::  tt(2)
  contains
    procedure :: Init           => Constant_Field_Init
    procedure :: Done           => Constant_Field_Done
    procedure :: CheckSample    => Constant_Field_CheckSample
    procedure :: PutSample      => Constant_Field_PutSample
    procedure :: GetTargetPtr   => Constant_Field_GetTargetPtr
  end type T_Constant_Field
  
  ! *
  
  type, extends(T_Constant_Field) :: T_Constant_Field_Series
    ! target time for next setup:
    type(TDate)               ::  setup_tt(2)
  contains
    procedure :: Setup          => Constant_Field_Series_Setup
    procedure :: Setup_Put      => Constant_Field_Series_Setup_Put
  end type T_Constant_Field_Series
  
  ! *
  
  type, extends(T_Constant_Field) :: T_Accumulated_Field_Series
    ! temporal interpolation:
    integer                   ::  interpolation_step
    character(len=32)         ::  interpolation_step_units
    character(len=64)         ::  interpolation
    real                      ::  interpolation_dt_sec
    ! previous data, interval for which it is valid:
    real, allocatable         ::  data_prev(:,:,:)   ! units/s
    type(TDate)               ::  tt_prev(2)
    logical                   ::  with_prev
    ! flags to know if and which data is to be receied:
    type(TDate)               ::  setup_tt(2)
    integer                   ::  nreceive
    type(TDate)               ::  receive_time(2)
    integer                   ::  receive_targ(2)
  contains
    procedure :: Init           => Accumulated_Field_Series_Init
    procedure :: Done           => Accumulated_Field_Series_Done
    procedure :: Set_Interp     => Accumulated_Field_Series_Set_Interp
    procedure :: Setup          => Accumulated_Field_Series_Setup
    procedure :: Setup_Prepare  => Accumulated_Field_Series_Setup_Prepare
    procedure :: Setup_InqTime  => Accumulated_Field_Series_Setup_InqTime
    procedure :: Setup_Put      => Accumulated_Field_Series_Setup_Put
    procedure :: Get_Tendency   => Accumulated_Field_Series_Get_Tendency
  end type T_Accumulated_Field_Series
  
  ! *
  
  type, extends(T_Field) :: T_Instant_Field
    ! time for which field is valid:
    type(TDate)               ::  t
  contains
    procedure :: Init           => Instant_Field_Init
    procedure :: Done           => Instant_Field_Done
    procedure :: CheckSample    => Instant_Field_CheckSample
    procedure :: PutSample      => Instant_Field_PutSample
    procedure :: GetTargetPtr   => Instant_Field_GetTargetPtr
  end type T_Instant_Field
  
  ! *
  
  type, extends(T_Instant_Field) :: T_Instant_Field_Series
    ! temporal interpolation:
    integer                   ::  interpolation_step
    character(len=32)         ::  interpolation_step_units
    character(len=64)         ::  interpolation
    !! tendency, interval for which it is valid:
    !real, pointer             ::  ddata_dt(:,:,:)   ! units/s
    ! data between which is interpolated:
    real, allocatable         ::  data1(:,:,:)
    real, allocatable         ::  data2(:,:,:)
    ! tendendy, allocated if necessary,
    ! pointer might refer to this array:
    real, pointer             ::  tend(:,:,:)
    ! interval for which it is valid:
    type(TDate)               ::  tt(2)
    !logical                   ::  with_tend
    logical                   ::  with_data1
    logical                   ::  with_data2
    ! flags to know if and which data is to be receied:
    type(TDate)               ::  setup_t
    integer                   ::  nreceive
    type(TDate)               ::  receive_time(2)
    integer                   ::  receive_targ(2)
  contains
    procedure :: Init           => Instant_Field_Series_Init
    procedure :: Done           => Instant_Field_Series_Done
    procedure :: Set_Interp     => Instant_Field_Series_Set_Interp
    procedure :: Interpol       => Instant_Field_Series_Interpol
    procedure :: Setup          => Instant_Field_Series_Setup
    procedure :: Setup_Prepare  => Instant_Field_Series_Setup_Prepare
    procedure :: Setup_InqTime  => Instant_Field_Series_Setup_InqTime
    procedure :: Setup_Put      => Instant_Field_Series_Setup_Put
    procedure :: Get_Tendency   => Instant_Field_Series_Get_Tendency
  end type T_Instant_Field_Series

  ! ***

  type, extends(T_FieldBase) :: T_VectorField
    ! data, pointers will refer to this array:
    real, pointer             ::  udata(:,:,:)
    real, pointer             ::  vdata(:,:,:)
  contains
    procedure :: Init          => VectorField_Init
    procedure :: Done          => VectorField_Done
    procedure :: Get           => VectorField_Get
    procedure :: Alloc         => VectorField_Alloc
    procedure :: AllocExt      => VectorField_AllocExt
    procedure :: Put           => VectorField_Put
  end type T_VectorField
  
  ! *
  
  type, extends(T_VectorField) :: T_Instant_VectorField
    ! time for which field is valid:
    type(TDate)               ::  t
  contains
    procedure :: Init           => Instant_VectorField_Init
    procedure :: Done           => Instant_VectorField_Done
    procedure :: CheckSample    => Instant_VectorField_CheckSample
    procedure :: PutSample      => Instant_VectorField_PutSample
  end type T_Instant_VectorField
  
  ! *
  
  type, extends(T_Instant_VectorField) :: T_Instant_VectorField_Series
    ! temporal interpolation:
    integer                   ::  interpolation_step
    character(len=32)         ::  interpolation_step_units
    character(len=64)         ::  interpolation
    ! tendency, interval for which it is valid;
    ! for vector interpolation, also length arrays are needed:
    type(TDate)               ::  t0
    !real, allocatable         ::  udata0(:,:,:)       ! units
    !real, allocatable         ::  vdata0(:,:,:)       ! units
    !real, allocatable         ::  dudata0_dt(:,:,:)   ! units
    !real, allocatable         ::  dvdata0_dt(:,:,:)   ! units
    !real, allocatable         ::  ldata0(:,:,:)       ! units
    !real, allocatable         ::  dldata0_dt(:,:,:)   ! units
    real, allocatable         ::  udata1(:,:,:)       ! units
    real, allocatable         ::  vdata1(:,:,:)       ! units
    real, allocatable         ::  ldata1(:,:,:)       ! units
    real, allocatable         ::  udata2(:,:,:)       ! units
    real, allocatable         ::  vdata2(:,:,:)       ! units
    real, allocatable         ::  ldata2(:,:,:)       ! units
    real, allocatable         ::  udata_(:,:,:)       ! units
    real, allocatable         ::  vdata_(:,:,:)       ! units
    real, allocatable         ::  ldata_(:,:,:)       ! units
    type(TDate)               ::  tt(2)
    !logical                   ::  with_tend
    logical                   ::  with_data1
    logical                   ::  with_data2
    ! flags to know if and which data is to be receied:
    type(TDate)               ::  setup_t
    integer                   ::  nreceive
    type(TDate)               ::  receive_time(2)
    integer                   ::  receive_targ(2)
  contains
    procedure :: Init           => Instant_VectorField_Series_Init
    procedure :: Done           => Instant_VectorField_Series_Done
    procedure :: Set_Interp     => Instant_VectorField_Series_Set_Interp
    procedure :: Interpol       => Instant_VectorField_Series_Interpol
    procedure :: Setup          => Instant_VectorField_Series_Setup
    procedure :: Setup_Prepare  => Instant_VectorField_Series_Setup_Prepare
    procedure :: Setup_InqTime  => Instant_VectorField_Series_Setup_InqTime
    procedure :: Setup_Put      => Instant_VectorField_Series_Setup_Put
  end type T_Instant_VectorField_Series
  

contains


  ! ====================================================================
  ! ===
  ! === FieldBase
  ! ===
  ! ====================================================================


  subroutine FieldBase_Init( self, name, units, status )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(out)               ::  self
    character(len=*), intent(in)                  ::  name, units
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Init'
      
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    self%units = trim(units)
    
    ! dummy value:
    self%input = ''
    
    ! no data yet:
    self%shp = -999
    self%lbo =    1
    self%ubo = -999

    ! no content yet:
    self%enabled   = .false.
    self%with_data = .false.

    ! ok
    status = 0
    
  end subroutine FieldBase_Init


  ! ***


  subroutine FieldBase_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(inout)             ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Done'
      
    ! --- begin ----------------------------------
    
    ! reset data:
    call FieldBase_Reset( self, status )
    IF_NOTOK_RETURN(status=1)

    ! empty:
    self%name = ''
    self%units = ''
    
    ! no data:
    self%shp = -999
    self%lbo =    1
    self%ubo = -999
    
    ! ok
    status = 0
    
  end subroutine FieldBase_Done


  ! ***


  subroutine FieldBase_Reset( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(inout)             ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Reset'
      
    ! --- begin ----------------------------------
    
    ! no content anymore:
    self%with_data = .false.

    ! ok
    status = 0
    
  end subroutine FieldBase_Reset


  ! ***


  !
  ! Check data in field:
  !   enabled   : compare enabled flag 
  ! For messages:
  !   silent    : do not raise error messages
  ! Return status:
  !   <1 if check failed
  !

  subroutine FieldBase_Check( self, status, enabled, filled, silent )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(inout)             ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
    logical, intent(in), optional                 ::  filled
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Check'
      
    ! --- local ----------------------------------
    
    logical       ::  shout
      
    ! --- begin ----------------------------------
    
    ! shout?
    shout = .true.
    if ( present(silent) ) shout = .not. silent
    
    ! check enabled flag ?
    if ( present(enabled) ) then
      ! compare:
      if ( self%enabled .neqv. enabled ) then
        ! not as expected ...
        if ( shout ) then
          if ( enabled ) then
            write (gol,'("field `",a,"` is not enabled as requested.")') trim(self%name); call goErr
          else
            write (gol,'("field `",a,"` is enabled but should be disabled.")') trim(self%name); call goErr
          end if
          TRACEBACK
        end if
        ! return with warning status:
        status=-1; return
      end if
    end if
    
    ! check filled flag ?
    if ( present(filled) ) then
      ! compare:
      if ( self%with_data .neqv. filled ) then
        ! not as expected ...
        if ( shout ) then
          if ( filled ) then
            write (gol,'("field `",a,"` is not filled as requested.")') trim(self%name); call goErr
          else
            write (gol,'("field `",a,"` is filled but should be not.")') trim(self%name); call goErr
          end if
          TRACEBACK
        end if
        ! return with warning status:
        status=-1; return
      end if
    end if

    ! ok
    status = 0
    
  end subroutine FieldBase_Check


  ! ***


  subroutine FieldBase_Set( self, status, enabled, input, lbo, ubo, nlev )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(inout)             ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
    character(len=*), intent(in), optional        ::  input
    integer, intent(in), optional                 ::  lbo(3), ubo(3)
    integer, intent(in), optional                 ::  nlev
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Set'
      
    ! --- begin ----------------------------------
    
    ! reset:
    if ( present(enabled) ) self%enabled = enabled
    if ( present(input  ) ) self%input   = trim(input)
    
    ! bounds?
    if ( present(lbo) .or. present(ubo) ) then
      ! check ...
      if ( (.not. present(lbo)) .or. (.not. present(ubo)) ) then
        write (gol,'("both lbo and ubo should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      self%lbo = lbo
      self%ubo = ubo
      ! set shape:
      self%shp = ubo - lbo + 1
    end if
    
    ! number of levels provided?
    if ( present(nlev) ) then
      ! setup data for levels 1:nlev :
      self%lbo(3) = 1
      self%lbo(3) = nlev
      ! set shape:
      self%shp(3) = nlev
    end if

    ! ok
    status = 0
    
  end subroutine FieldBase_Set


  ! ***


  subroutine FieldBase_Get( self, status, &
                         enabled, filled, &
                         name, units, input, &
                         lbo, ubo, shp )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(in)                ::  self
    integer, intent(out)                          ::  status
    logical, intent(out), optional                ::  enabled
    logical, intent(out), optional                ::  filled
    character(len=*), intent(out), optional       ::  name, units, input
    integer, intent(out), optional                ::  lbo(3), ubo(3), shp(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_Get'
      
    ! --- begin ----------------------------------
    
    ! description:
    if ( present(name   ) ) name    = self%name
    if ( present(units  ) ) units   = self%units
    
    ! flags:
    if ( present(enabled) ) enabled = self%enabled
    if ( present(filled ) ) filled  = self%with_data

    ! check ..
    if ( any( (/ present(input), &
                 present(lbo), present(ubo), present(shp) /) ) &
         .and. (.not. self%enabled) ) then
      write (gol,'("requested properties of `",a,"` field while not enabled")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    else
      ! input source:
      if ( present(input  ) ) input   = self%input
      ! data bounds:
      if ( present(lbo    ) ) lbo     = self%lbo
      if ( present(ubo    ) ) ubo     = self%ubo
      if ( present(shp    ) ) shp     = self%shp
    end if

    ! ok
    status = 0
    
  end subroutine FieldBase_Get


  ! ***
  

  subroutine FieldBase_SetupShape( self, status, lbo, ubo, change_levs )
    
    ! --- in/out ---------------------------------
    
    class(T_FieldBase), intent(inout)   ::  self
    integer, intent(out)                ::  status
    integer, intent(in), optional       ::  lbo(3), ubo(3)
    logical, intent(in), optional       ::  change_levs
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/FieldBase_SetupShape'
      
    ! --- local ----------------------------------
    
    integer         ::  i
    logical         ::  change_levels
      
    ! --- begin ----------------------------------
    
    ! change levels?
    change_levels = .false.
    if ( present(change_levs) ) change_levels = change_levs
    
    ! enabled ?
    if ( .not. self%enabled ) then
      write (gol,'("field `",a,"` not enabled")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! index space not defined yet?
    if ( present(lbo) .or. present(ubo) ) then
      ! check ...
      if ( (.not. present(lbo)) .or. (.not. present(ubo)) ) then
        write (gol,'("both lbo and ubo should be provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! loop over dimensions:
      do i = 1, 3
        ! not defined yet ?
        if ( self%shp(i) < 0 .or. (i == 3 .and. change_levels) ) then
          ! set bounds:
          self%lbo(i) = lbo(i)
          self%ubo(i) = ubo(i)
          ! reset shape:
          self%shp(i) = ubo(i) - lbo(i) + 1
        end if
      end do  ! dims
      
    else
      if ( any(self%shp < 0) ) then
        ! provided ?
        write (gol,'("index space not defined yet:")'); call goErr
        write (gol,'("  lbo  : ",3i6)') self%lbo; call goErr
        write (gol,'("  ubo  : ",3i6)') self%ubo; call goErr
        write (gol,'("  shp  : ",3i6)') self%shp; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine FieldBase_SetupShape


  ! ====================================================================
  ! ===
  ! === Field
  ! ===
  ! ====================================================================


  subroutine Field_Init( self, name, units, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(out)                   ::  self
    character(len=*), intent(in)                  ::  name, units
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Init'
      
    ! --- begin ----------------------------------
    
    ! init base:
    call FieldBase_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)
    
    ! no data yet:
    nullify( self%data )

    ! ok
    status = 0
    
  end subroutine Field_Init


  ! ***


  subroutine Field_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)                 ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Done'
      
    ! --- begin ----------------------------------
    
    ! reset data:
    call Field_Reset( self, status )
    IF_NOTOK_RETURN(status=1)

    ! done with base:
    call FieldBase_Done( self, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Field_Done


  ! ***


  subroutine Field_Reset( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)                 ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Reset'
      
    ! --- begin ----------------------------------
    
    ! no data:
    if ( associated(self%data) ) then
      deallocate( self%data, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! reset base:
    call FieldBase_Reset( self, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Field_Reset


  ! ***


  subroutine Field_Get( self, status, &
                         enabled, filled, &
                         name, units, input, &
                         lbo, ubo, shp, &
                         pdata )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(in)                    ::  self
    integer, intent(out)                          ::  status
    logical, intent(out), optional                ::  enabled
    logical, intent(out), optional                ::  filled
    character(len=*), intent(out), optional       ::  name, units, input
    integer, intent(out), optional                ::  lbo(3), ubo(3), shp(3)
    real, pointer, optional                       ::  pdata(:,:,:)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Get'
      
    ! --- begin ----------------------------------
    
    ! get from base:
    call FieldBase_Get( self, status, &
                         enabled=enabled, filled=filled, &
                         name=name, units=units, input=input, &
                         lbo=lbo, ubo=ubo, shp=shp )
    IF_NOTOK_RETURN(status=1)
    
    ! check ..
    if ( present(pdata) ) then
      ! check ..
      if ( .not. self%enabled ) then
        write (gol,'("requested properties of `",a,"` field while not enabled")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! set pointer:
      pdata   => self%data
    end if

    ! ok
    status = 0
    
  end subroutine Field_Get


  ! ***
  

  subroutine Field_Alloc( self, status, lbo, ubo, change_levs )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)       ::  self
    integer, intent(out)                ::  status
    integer, intent(in), optional       ::  lbo(3), ubo(3)
    logical, intent(in), optional       ::  change_levs
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_Alloc'
      
    ! --- local ----------------------------------
    
    integer         ::  i
    logical         ::  change_levels
      
    ! --- begin ----------------------------------
    
    ! change levels?
    change_levels = .false.
    if ( present(change_levs) ) change_levels = change_levs
    
    ! setup shp attributes:
    call FieldBase_SetupShape( self, status, lbo=lbo, ubo=ubo, change_levs=change_levels )
    IF_NOTOK_RETURN(status=1)
    
    ! need to allocate ?
    if ( .not. associated(self%data) ) then
      ! storage:
      allocate( self%data(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
    else if ( change_levels .and. (any(lbound(self%data) /= self%lbo) .or. any(ubound(self%data) /= self%ubo) ) ) then
      ! allocate with new shape
      deallocate( self%data, stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%data(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine Field_Alloc
  

  ! ***
  
  
  ! allocate external array with same index range
  
  subroutine Field_AllocExt( self, data, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)                 ::  self
    real, intent(inout), allocatable              ::  data(:,:,:)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Field_AllocExt'
      
    ! --- begin ----------------------------------
    
    ! enabled ?
    if ( .not. self%enabled ) then
      write (gol,'("field `",a,"` not enabled")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( data(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Field_AllocExt
  

  ! ***


  subroutine Field_Put( self, data, status, lbo, ubo, change_levs )
  
    ! --- in/out ---------------------------------
    
    class(T_Field), intent(inout)       ::  self
    real, intent(in)                    ::  data(:,:,:)
    integer, intent(out)                ::  status
    integer, intent(in), optional       ::  lbo(3), ubo(3)
    logical, intent(in), optional       ::  change_levs

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Field_Put'

    ! --- local ----------------------------------
    
    !integer     ::  i, j, k
    
    ! --- begin ----------------------------------
    
    ! allocate data if necessary, setup index space if not done yet:
    call self%Alloc( status, lbo=lbo, ubo=ubo, change_levs=change_levs )
    IF_NOTOK_RETURN(status=1)
    
    ! check ...
    if ( any(shape(data) /= self%shp) ) then
      write (gol,'("shape of data does not match with definition:")'); call goErr
      write (gol,'("  data     : ",3i6)') shape(data); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=1; return
    end if

    ! copy data:
    self%data = data
    ! set flag:
    self%with_data = .true.

    !! testing ...
    !write (gol,*) 'fff ', trim(self%name),'(5,28,1) = ', self%data(5,28,1); call goPr

    !! check ..
    !if ( any(isnan(self%data)) ) then
    !  write (gol,'("found NaN values in `",a,"`")') trim(self%name); call goErr
    !  do k = self%lbo(3), self%ubo(3)
    !    do j = self%lbo(2), self%ubo(2)
    !      do i = self%lbo(1), self%ubo(1)
    !        if ( isnan(self%data(i,j,k)) ) then
    !          write (gol,'("  cell ",3i6)') i, j, k; call goErr
    !        end if
    !      end do ! i
    !    end do ! j
    !  end do ! k
    !  TRACEBACK; status=1; return
    !end if

    ! ok
    status = 0
    
  end subroutine Field_Put
  


  

  ! ====================================================================
  ! ===
  ! === Constant_Field
  ! ===
  ! ====================================================================


  subroutine Constant_Field_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field), intent(out)   ::  self
    character(len=*), intent(in)          ::  name, units
    integer, intent(out)                  ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Constant_Field_Init'

    ! --- begin ----------------------------------
    
    ! initialize the Constant field:
    call Field_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no time range assigned yet:
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()

    ! ok
    status = 0
    
  end subroutine Constant_Field_Init
  

  ! ***


  subroutine Constant_Field_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field), intent(inout)    ::  self
    integer, intent(out)                     ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Constant_Field_Done'
      
    ! --- begin ----------------------------------
    
    ! no current data:
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    
    ! done with field:
    call Field_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Constant_Field_Done


  !
  ! Check data in field:
  !   enabled   : compare enabled flag 
  !   tt(2)     : required time interval, should be within valid interval
  ! For messages:
  !   silent    : do not raise error messages
  ! Return status:
  !   <1 if check failed

  subroutine Constant_Field_CheckSample( self, status, enabled, tt, silent )
  
    use GO_Date, only : TDate, operator(<), wrtgol
    
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field), intent(inout)        ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
    type(TDate), intent(in), optional             ::  tt(2)
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Constant_Field_CheckSample'
      
    ! --- local ----------------------------------
    
    logical       ::  shout
      
    ! --- begin ----------------------------------
    
    ! original:
    call FieldBase_Check( self, status, enabled=enabled, silent=silent )
    ! leave with warning status ?
    if ( status < 0 ) return
    ! failed:
    IF_NOTOK_RETURN(status=1)
    
    ! shout?
    shout = .true.
    if ( present(silent) ) shout = .not. silent
    
    ! check time ?
    if ( present(tt) ) then
      ! filled at all?
      if ( .not. self%with_data ) then
        ! not as expected ...
        if ( shout ) then
          write (gol,'("no data filled yet")'); call goErr
          TRACEBACK
        end if
        ! return with warning status:
        status=-1; return
      else
        ! compare:
        if ( (tt(1) < self%tt(1)) .or. (self%tt(2) < tt(2)) ) then
          ! not as expected ...
          if ( shout ) then
            call wrtgol( 'expected timerange ', tt ); call goErr
            call wrtgol( '  but sample valid for ', self%tt ); call goErr
            TRACEBACK
          end if
          ! return with warning status:
          status=-1; return
        end if
      end if
    end if

    ! ok
    status = 0
    
  end subroutine Constant_Field_CheckSample
  

  ! ***


  subroutine Constant_Field_PutSample( self, data, tt, status, lbo, ubo )
  
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field), intent(inout)  ::  self
    real, intent(in)                        ::  data(:,:,:)
    type(TDate), intent(in)                 ::  tt(2)
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Constant_Field_PutSample'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store data:
    call self%Put( data, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! store time range:
    self%tt = tt

    ! ok
    status = 0
    
  end subroutine Constant_Field_PutSample
  

  ! ***


  subroutine Constant_Field_GetTargetPtr( self, pdata, tt, status, lbo, ubo )
  
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field), intent(inout)  ::  self
    real, pointer                           ::  pdata(:,:,:)
    type(TDate), intent(in)                 ::  tt(2)
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Constant_Field_GetTargetPtr'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! allocate data if not done yet, eventually using provided bounds:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer:
    call self%Get( status, pdata=pdata )
    IF_NOTOK_RETURN(status=1)

    ! store time value:
    self%tt = tt
    ! set flag already now, user is responsible to fill the array ..
    self%with_data = .true.

    ! ok
    status = 0
    
  end subroutine Constant_Field_GetTargetPtr
  

  ! ====================================================================
  ! ===
  ! === Constant_Field_Series
  ! ===
  ! ====================================================================


  subroutine Constant_Field_Series_Setup( self, tt, status )
  
    use GO_Date, only : operator(<)
    use GO_Date, only : wrtgol, Pretty
 
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field_Series), intent(inout)   ::  self
    type(TDate), intent(in)                         ::  tt(2)
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Constant_Field_Series_Setup'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target time:
    self%setup_tt = tt
    
    ! no data present yet ?
    if ( .not. self%with_data ) then
      !! info ...
      !write (gol,'(a,": no data filled yet")') rname; call goPr
      ! new data required ...
      status = -1; return
    else
      ! data defined; not valid for this time ?
      if ( (tt(1) < self%tt(1)) .or. (self%tt(2) < tt(2)) ) then
        !! info ...
        !write (gol,'(a,": target interval [",a,",",a,"] outside current interval [",a,",",a,"]; need to receive new data")') &
        !          rname, trim(Pretty(tt(1))), trim(Pretty(tt(2))), &
        !          trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
        ! new data required:
        status = -1; return
      else
        !! info ...
        !write (gol,'(a,": target interval [",a,",",a,"] inside current interval [",a,",",a,"]; keep current data")') &
        !          rname, trim(Pretty(tt(1))), trim(Pretty(tt(2))), &
        !          trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine Constant_Field_Series_Setup
  

  ! ***


  subroutine Constant_Field_Series_Setup_Put( self, data, tt, status, lbo, ubo )
  
    use GO_Date, only : operator(<), wrtgol
    
    ! --- in/out ---------------------------------
    
    class(T_Constant_Field_Series), intent(inout)   ::  self
    real, intent(in)                                ::  data(:,:,:)
    type(TDate), intent(in)                         ::  tt(2)
    integer, intent(out)                            ::  status
    integer, intent(in), optional                   ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Constant_Field_Series_Setup_Put'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (self%setup_tt(1) < tt(1)) .or. (tt(2) < self%setup_tt(2)) ) then
      write (gol,'("supplied time range of constant data does not include target time:")'); call goErr
      call wrtgol( '  time range   : ', tt(1), ' to ', tt(2) ); call goErr
      call wrtgol( '  target time  : ', self%setup_tt(1), ' to ', self%setup_tt(2) ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! store data, use index space if necessary:
    call Constant_Field_PutSample( self, data, tt, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Constant_Field_Series_Setup_Put
  


  ! ====================================================================
  ! ===
  ! === Accumulated_Field_Series
  ! ===
  ! ====================================================================


  subroutine Accumulated_Field_Series_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(out)  ::  self
    character(len=*), intent(in)                    ::  name, units
    integer, intent(out)                            ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Init'

    ! --- begin ----------------------------------
    
    ! initialize the constant field:
    call Constant_Field_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no previous data yet:
    self%with_prev = .false.

    ! no time stepping set yet:
    self%interpolation_step = -999
    self%interpolation_step_units = 'None'
    self%interpolation = 'None'

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Init
  

  ! ***


  subroutine Accumulated_Field_Series_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(inout)    ::  self
    integer, intent(out)                                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Done'
      
    ! --- begin ----------------------------------
    
    ! no tendency data anymore:
    if ( self%with_prev ) then
      ! clear:
      deallocate( self%data_prev, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! reset:
      self%with_prev = .false.
    end if
    ! dummy:
    self%tt_prev(1) = AnyDate()
    self%tt_prev(2) = AnyDate()
    
    ! done with field:
    call Constant_Field_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Done
  

  ! ***


  subroutine Accumulated_Field_Series_Set_Interp( self, interp, status )
  
    use GO_String, only : goReadFromLine
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(inout)  ::  self
    character(len=*), intent(in)                      ::  interp
    integer, intent(out)                              ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Set_Interp'

    ! --- local ----------------------------------
    
    character(len=64)        ::    line

    ! --- begin ----------------------------------
    
    ! save:
    self%interpolation = trim(interp)
    
    ! copy:
    line = trim(interp)
    
    ! format:  3 hour
    ! read time step:
    call goReadFromLine( line, self%interpolation_step, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! rest is step units:
    self%interpolation_step_units = trim(line)
    
    ! time step in second, used to determine if previous accumulated field is needed too:
    select case ( trim(self%interpolation_step_units) )
      case ( 'hour' )
        self%interpolation_dt_sec = self%interpolation_step * 3600.0
      case default
        write (gol,'("interpolation step units `",a,"` not supported yet")') self%interpolation_step; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! no previous data yet:
    self%tt_prev(1) = AnyDate()
    self%tt_prev(2) = AnyDate()
    self%with_prev = .false.

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Set_Interp
  

  ! ***


  !
  ! Check if data has been setup for requested time interval.
  ! Return status:
  !   <0  :  need to setup new data
  !    0  :  data already valid
  !   >0  :  error
  !

  subroutine Accumulated_Field_Series_Setup( self, tt, status )
  
    use GO_Date, only : operator(/=), operator(<)
    use GO_Date, only : wrtgol, Pretty
 
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(inout)::  self
    type(TDate), intent(in)                         ::  tt(2)
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Setup'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target time:
    self%setup_tt = tt
    
    ! no data present yet ?
    if ( .not. self%with_data ) then
      !! info ...
      !write (gol,'(a,": no data filled yet")') rname; call goPr
      ! new data required ...
      status = -1; return
    else
      ! data defined; also previous data defined?
      if ( self%with_prev ) then
        ! requested interval should be between end of previous and end of current;
        ! check if it is partly outside:
        if ( (tt(1) < self%tt_prev(2)) .or. (self%tt(2) < tt(2)) ) then
          ! info ...
          write (gol,'(a,": target interval [",a,",",a,"] not inside difference between")') &
                    rname, trim(Pretty(tt(1))), trim(Pretty(tt(2))); call goPr
          write (gol,'(a,": previous inteval [",a,",",a,"] and current [",a,",",a,"]; need to receive new data")') &
                    rname, trim(Pretty(self%tt_prev(1))), trim(Pretty(self%tt_prev(2))), &
                    trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
          ! new data required:
          status = -1; return
        else
          ! info ...
          write (gol,'(a,": target interval [",a,",",a,"] inside difference between")') &
                    rname, trim(Pretty(tt(1))), trim(Pretty(tt(2))); call goPr
          write (gol,'(a,": previous inteval [",a,",",a,"] and current [",a,",",a,"]")') &
                    rname, trim(Pretty(self%tt_prev(1))), trim(Pretty(self%tt_prev(2))), &
                    trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
        end if
      else
        ! only a current interval; requested interval should be inside:
        if ( (tt(1) < self%tt(1)) .or. (self%tt(2) < tt(2)) ) then
          ! info ...
          write (gol,'(a,": target interval [",a,",",a,"] not inside interval [",a,",",a,"]; need to receive new data")') &
                    trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
          ! new data required:
          status = -1; return
        else
          ! info ...
          write (gol,'(a,": target interval [",a,",",a,"] inside interval [",a,",",a,"]")') &
                    trim(Pretty(self%tt(1))), trim(Pretty(self%tt(2))); call goPr
        end if
      end if  ! prev defined
    end if

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Setup
  

  ! ***


  subroutine Accumulated_Field_Series_Setup_Prepare( self, tref, nreceive, status )
  
    use GO_Date, only : TDate
    use GO_Date, only : operator(+), operator(-), operator(*)
    use GO_Date, only : operator(/=)!, operator(==), operator(<), operator(<=)
    use GO_Date, only : wrtgol
    use GO_Date, only : Get_Surrounding_Interval
 
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(inout)::  self
    type(TDate), intent(in)                         ::  tref
    integer, intent(out)                            ::  nreceive
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Setup_Prepare'
    
    ! --- local ----------------------------------
    
    type(TDate)   ::  tmid
    type(TDate)   ::  xtt(2)

    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': prepare setup for ', self%setup_tt ); call goPr
    
    ! some fields need to be changed; get mid time:
    tmid = self%setup_tt(1) + 0.5*(self%setup_tt(2) - self%setup_tt(1))
    ! get interpolation points given temporal resolution:
    call Get_Surrounding_Interval( tmid, &
                                   tref, self%interpolation, &
                                   xtt, status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    call wrtgol( rname//': surrounding interval ', xtt ); call goPr
    
    ! current data present?
    if ( self%with_data ) then
      ! current data is expected to end at start of surrounding interval; check:
      if ( self%tt(2) /= xtt(1) ) then
        call wrtgol( 'current interval ', self%tt ); call goErr
        call wrtgol( 'is expected to end at new surrounding interval ', xtt ); call goErr
        TRACEBACK; status=1; return
      end if
      ! no previous data allocated yet?
      if ( .not. self%with_prev ) then
        ! storage:
        allocate( self%data_prev(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! reset flag:
        self%with_prev = .true.
      end if
      ! copy:
      self%data_prev = self%data
      self%tt_prev   = self%tt
      ! info ...
      call wrtgol( rname//':   copied current data to previous, valid for ', self%tt_prev ); call goPr
      write (gol,'(a,":   retrieve new current data ...")') rname; call goPr
      ! need to obtain new current data that ends at surrounding interval:
      self%nreceive = 1
      self%receive_time(1) = xtt(2)
      self%receive_targ(1) = ACCUF_CURR
    else
      ! info ...
      write (gol,'(a,":   retrieve new previous and current data ...")') rname; call goPr
      ! retrieve both current and previous ;
      ! first retrieve current, if that is already valid for the requested time resolution
      ! then the previous will be skipped (by setting AnyDate in the time)
      self%nreceive = 2
      self%receive_time(1) = xtt(2)
      self%receive_targ(1) = ACCUF_CURR
      self%receive_time(2) = xtt(1)
      self%receive_targ(2) = ACCUF_PREV
    end if
    
    ! return value:
    nreceive = self%nreceive

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Setup_Prepare
  

  ! ***


  subroutine Accumulated_Field_Series_Setup_InqTime( self, ip, t, status )
 
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(in)   ::  self
    integer, intent(in)                             ::  ip
    type(TDate), intent(out)                        ::  t
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Setup_InqTime'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract:
    t = self%receive_time(ip)
    
    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Setup_InqTime
  

  ! ***


  subroutine Accumulated_Field_Series_Setup_Put( self, ip, data, tt, status, lbo, ubo )
  
    use GO_Date, only : operator(/=), operator(==), operator(-), rTotal, AnyDate
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(inout)::  self
    integer, intent(in)                             ::  ip
    real, intent(in)                                ::  data(:,:,:)
    type(TDate), intent(in)                         ::  tt(2)
    integer, intent(out)                            ::  status
    integer, intent(in), optional                   ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Setup_Put'

    ! --- local ----------------------------------
    
    real        ::  dt_sec

    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') ip, self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! allocate data if not done yet, setup index space if necessary:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
  
    ! check ...
    if ( any(shape(data) /= self%shp) ) then
      write (gol,'("shape of data does not match with definition:")'); call goErr
      write (gol,'("  data     : ",3i6)') shape(data); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=1; return
    end if

    ! check ...
    if ( tt(2) /= self%receive_time(ip) ) then
      write (gol,'("end time of received field ",i6," does not match with setup:")') ip; call goErr
      call wrtgol( '  setup time     : ', self%receive_time(ip) ); call goErr
      call wrtgol( '  received times : ', tt ); call goErr
      TRACEBACK; status=1; return
    end if

    ! what to do with this data ?
    select case ( self%receive_targ(ip) )
      !~~ use for previous data?
      case ( ACCUF_PREV )
        ! setup previous data if necessary:
        if ( .not. allocated(self%data_prev) ) then
          ! storage:
          allocate( self%data_prev(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! store:
        self%data_prev = data
        self%tt_prev   = tt
        self%with_prev = .true.
        ! info ...
        call wrtgol( rname//': stored previous data for ', self%tt_prev ); call goPr
        !write (gol,'(a,": average value : ",f12.4)') rname, sum(self%data)/size(self%data); call goPr
        !
      !~~ store as current
      case ( ACCUF_CURR )
        ! store data:
        self%data = data
        ! store time:
        self%tt = tt
        ! set flag:
        self%with_data = .true.
        ! info ...
        call wrtgol( rname//': stored current data for ', self%tt ); call goPr
        !write (gol,'(a,": average value : ",f12.4)') rname, sum(self%data)/size(self%data); call goPr
        !
        ! should match with previous if available:
        if ( self%with_prev ) then
          ! start of current interval is either same as start of previous (same accumulation),
          ! or same as end time of previous (new accumulation); otherwise error:
          if ( self%tt(1) /= self%tt_prev(1) ) then
            ! new accumulutation?
            if ( self%tt(1) == self%tt_prev(2) ) then
              ! reset previous to zero:
              self%data_prev = 0.0
              self%tt_prev   = (/self%tt(1),self%tt(1)/)
              ! info ...
              write (gol,'(a,": new accumulation,  reset previous data to zero ..")') rname; call goPr
            else
              ! strange ...
              write (gol,'("start time of current field should be start or end of previous:")'); call goErr
              call wrtgol( '  time range of previous data : ', self%tt_prev ); call goErr
              call wrtgol( '  time range of current data  : ', self%tt ); call goErr
              TRACEBACK; status=1; return
            end if
          end if ! different start times
        end if ! prev data present
        !
        ! no need for previous data in ip=2 if time step is already the interpolation step;
        ! check if second step is needed:
        if ( (ip == 1) .and. (self%nreceive > 1) ) then
          ! current time step:
          dt_sec = rTotal( tt(2)-tt(1), 'sec' )
          ! match?
          if ( abs( dt_sec - self%interpolation_dt_sec ) < 0.01 ) then
            ! info ...
            write (gol,'(a,": accumulation already valid for interpolation interval, no need for previous field ...")') rname; call goPr
            ! reset request time for second field:
            self%receive_time(2) = AnyDate()
          end if
        end if
        !
      !~~ strange ...
      case default
        ! something went wrong; bug ?
        write (gol,'("unsupported destination code ",i6," for interpolation point ",i6)') self%receive_targ(ip), ip; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! reset destination to dummy value:
    self%receive_targ(ip) = -1

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Setup_Put


  ! ***


  subroutine Accumulated_Field_Series_Get_Tendency( self, tr, ddata_dt, status )

    use GO_Date, only : TDate, operator(<=), operator(<), operator(/=), operator(-), rTotal, wrtgol

    ! --- in/out ---------------------------------
    
    class(T_Accumulated_Field_Series), intent(in)   ::  self
    type(TDate), intent(in)                         ::  tr(2)
    real, intent(out)                               ::  ddata_dt(:,:,:)
    integer, intent(out)                            ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Accumulated_Field_Series_Get_Tendency'
      
    ! --- local ----------------------------------
    
    real        ::  dt_sec
    
    ! --- begin ----------------------------------

    ! check output shape:
    if ( any(shape(ddata_dt) /= self%shp) ) then
      write (gol,'("shape of data does not match with definition:")'); call goErr
      write (gol,'("  data     : ",3i6)') shape(ddata_dt); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=1; return
    end if

    ! also previous field present?
    if ( self%with_prev ) then
    
      ! check start times:
      if ( self%tt(1) /= self%tt_prev(1) ) then
        write (gol,'("start times of accumulated fields should match:")'); call goErr
        call wrtgol( '  time range of previous data : ', self%tt_prev ); call goErr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! check end times:
      if ( self%tt(2) <= self%tt_prev(1) ) then
        write (gol,'("end time of current fields should be after end of previous:")'); call goErr
        call wrtgol( '  time range of previous data : ', self%tt_prev ); call goErr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! check interval:
      if ( (tr(1) < self%tt_prev(2)) .or. (self%tt(2) < tr(2)) ) then
        write (gol,'("requested interval should be within end times of previous and current data:")'); call goErr
        call wrtgol( '  time range of previous data : ', self%tt_prev ); call goErr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        call wrtgol( '  requested time range        : ', tr ); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! length of interval:
      dt_sec = rTotal( self%tt(2)-self%tt_prev(2), 'sec' )
      
      ! check ...
      if ( abs(dt_sec - self%interpolation_dt_sec) > 0.01 ) then
        write (gol,'("interval between previous and current accumulation does not match with requested temporal resolution:")'); call goPr
        call wrtgol( '  time range of previous data : ', self%tt_prev ); call goErr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        write (gol,'("  assumed temporal resolution : ",f10.2," sec")') self%interpolation_dt_sec; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! compute result:
      ddata_dt = ( self%data - self%data_prev ) / dt_sec
      
    else
    
      ! only current data; check interval:
      if ( (tr(1) < self%tt(1)) .or. (self%tt(2) < tr(2)) ) then
        write (gol,'("requested interval should be within current data:")'); call goErr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        call wrtgol( '  requested time range        : ', tr ); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! length of interval:
      dt_sec = rTotal( self%tt(2)-self%tt(1), 'sec' )
      
      ! check ...
      if ( abs(dt_sec - self%interpolation_dt_sec) > 0.01 ) then
        write (gol,'("interval of current accumulation does not match with requested temporal resolution:")'); call goPr
        call wrtgol( '  time range of current data  : ', self%tt ); call goErr
        write (gol,'("  assumed temporal resolution : ",f10.2," sec")') self%interpolation_dt_sec; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! compute result:
      ddata_dt = self%data / dt_sec

    end if

    ! ok
    status = 0
    
  end subroutine Accumulated_Field_Series_Get_Tendency
  


  

  ! ====================================================================
  ! ===
  ! === Instant_Field
  ! ===
  ! ====================================================================


  subroutine Instant_Field_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(out)   ::  self
    character(len=*), intent(in)          ::  name, units
    integer, intent(out)                  ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Init'

    ! --- begin ----------------------------------
    
    ! initialize the instant field:
    call Field_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no time assigned yet:
    self%t = AnyDate()

    ! ok
    status = 0
    
  end subroutine Instant_Field_Init
  

  ! ***


  subroutine Instant_Field_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(inout)    ::  self
    integer, intent(out)                     ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Done'
      
    ! --- begin ----------------------------------
    
    ! reset data:
    call Instant_Field_Reset( self, status )
    IF_NOTOK_RETURN(status=1)

    ! done with field:
    call Field_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Done
  

  ! ***


  subroutine Instant_Field_Reset( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(inout)    ::  self
    integer, intent(out)                     ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Reset'
      
    ! --- begin ----------------------------------
    
    ! no current data:
    self%t = AnyDate()
    
    ! reset field:
    call Field_Reset( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Reset


  !
  ! Check data in field:
  !   enabled   : compare enabled flag 
  !   t         : time
  ! For messages:
  !   silent    : do not raise error messages
  ! Return status:
  !   <1 if check failed

  subroutine Instant_Field_CheckSample( self, status, enabled, t, silent )
  
    use GO_Date, only : TDate, operator(/=), wrtgol
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(inout)         ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
    type(TDate), intent(in), optional             ::  t
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_CheckSample'
      
    ! --- local ----------------------------------
    
    logical       ::  shout
      
    ! --- begin ----------------------------------
    
    ! original:
    call FieldBase_Check( self, status, enabled=enabled, silent=silent )
    ! leave with warning status ?
    if ( status < 0 ) return
    ! failed:
    IF_NOTOK_RETURN(status=1)
    
    ! shout?
    shout = .true.
    if ( present(silent) ) shout = .not. silent
    
    ! check time ?
    if ( present(t) ) then
      ! filled at all?
      if ( .not. self%with_data ) then
        ! not as expected ...
        if ( shout ) then
          write (gol,'("no data filled yet")'); call goErr
          TRACEBACK
        end if
        ! return with warning status:
        status=-1; return
      else
        ! compare:
        if ( self%t /= t ) then
          ! not as expected ...
          if ( shout ) then
            call wrtgol( 'expected time ', t, ' but sample valid for ', self%t ); call goErr
            TRACEBACK
          end if
          ! return with warning status:
          status=-1; return
        end if ! compare
      end if ! filled?
    end if ! present t ?

    ! ok
    status = 0
    
  end subroutine Instant_Field_CheckSample
  

  ! ***


  subroutine Instant_Field_PutSample( self, data, t, status, lbo, ubo, change_levs )
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(inout)   ::  self
    real, intent(in)                        ::  data(:,:,:)
    type(TDate), intent(in)                 ::  t
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  lbo(3), ubo(3)
    logical, intent(in), optional           ::  change_levs

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_PutSample'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store data, use index space if not defined yet:
    call self%Put( data, status, lbo=lbo, ubo=ubo, change_levs=change_levs )
    IF_NOTOK_RETURN(status=1)

    ! store time value:
    self%t = t

    ! ok
    status = 0
    
  end subroutine Instant_Field_PutSample
  

  ! ***


  subroutine Instant_Field_GetTargetPtr( self, pdata, t, status, lbo, ubo )
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field), intent(inout)   ::  self
    real, pointer                           ::  pdata(:,:,:)
    type(TDate), intent(in)                 ::  t
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_GetTargetPtr'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! allocate data if not done yet, eventually using provided bounds:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer:
    call self%Get( status, pdata=pdata )
    IF_NOTOK_RETURN(status=1)

    ! store time value:
    self%t = t
    ! set flag already now, user is responsible to fill the array ..
    self%with_data = .true.

    ! ok
    status = 0
    
  end subroutine Instant_Field_GetTargetPtr
  


  ! ====================================================================
  ! ===
  ! === Instant_Field_Series
  ! ===
  ! ====================================================================


  subroutine Instant_Field_Series_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(out)  ::  self
    character(len=*), intent(in)                ::  name, units
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Init'

    ! --- begin ----------------------------------
    
    ! initialize the instant field:
    call Instant_Field_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no interpolation points yet:
    self%with_data1 = .false.
    self%with_data2 = .false.
    ! no tendencay array yet:
    nullify( self%tend )

    ! no interpolation set yet:
    self%interpolation_step = -999
    self%interpolation_step_units = 'None'
    self%interpolation = 'None'

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Init
  

  ! ***


  subroutine Instant_Field_Series_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    integer, intent(out)                            ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Done'
      
    ! --- begin ----------------------------------
    
    ! reset data:
    call Instant_Field_Series_Reset( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with field:
    call Instant_Field_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Done
  

  ! ***


  subroutine Instant_Field_Series_Reset( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    integer, intent(out)                            ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Reset'
      
    ! --- begin ----------------------------------
    
    ! no interpolation data anymore:
    if ( self%with_data1 ) then
      deallocate( self%data1, stat=status )
      IF_NOTOK_RETURN(status=1)
      self%with_data1 = .false.
    end if
    if ( self%with_data2 ) then
      deallocate( self%data2, stat=status )
      IF_NOTOK_RETURN(status=1)
      self%with_data2 = .false.
    end if
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    
    ! extra ..
    if ( associated(self%tend) ) then
      deallocate( self%tend, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! reset field:
    call Instant_Field_Reset( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Reset
  

  ! ***


  subroutine Instant_Field_Series_Set_Interp( self, interp, status )
  
    use GO_String, only : goReadFromLine
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)  ::  self
    character(len=*), intent(in)                  ::  interp
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Set_Interp'

    ! --- local ----------------------------------
    
    character(len=64)        ::    line

    ! --- begin ----------------------------------
    
    ! save:
    self%interpolation = trim(interp)
    
    ! copy:
    line = trim(interp)
    
    ! format:  3 hour
    ! read time step:
    call goReadFromLine( line, self%interpolation_step, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! rest is step units:
    self%interpolation_step_units = trim(line)
    
    ! no tendency data yet:
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    !self%with_tend = .false.
    self%with_data1 = .false.
    self%with_data2 = .false.

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Set_Interp
  

  ! ***


  ! Interpolate to t given available data and tendency.
  ! Tool for use within this module only.
  
  subroutine Instant_Field_Series_Interpol( self, t, status )
  
    use GO_Date, only : operator(/=), operator(-), rTotal
    use GO_Date, only : wrtgol, operator(==), operator(<), wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Interpol'

    ! --- local ----------------------------------
    
    real      ::  dt_sec
    real      ::  alfa
    
    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': interpolate between [', self%tt(1), ',', self%tt(2), '] to ', t ); call goPr
    
    ! interval bounds?
    if ( self%with_data1 .and. (t == self%tt(1)) ) then
    
      ! copy:
      self%data = self%data1
    
    else
    
      ! check ...
      if ( (.not. self%with_data1) .or. (.not. self%with_data2) ) then
        write (gol,'("not all data present: data1 (",l1,"), data2 (",l1,")")') &
                self%with_data1, self%with_data2; call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( (t < self%tt(1)) .or. (self%tt(2) < t) ) then
        call wrtgol( 'target time ', t ); call goErr
        call wrtgol( 'outside data interval ', self%tt ); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! data defined by start and end;
      ! interpolate:  data(t) = data(t1) * (t2-t)/(t2-t1) + ddata(self%t) * (t-t1)/(t2-t1)
      !                                        alfa                            (1-alfa)
      ! total timestep in seconds:
      dt_sec = rTotal( self%tt(2) - self%tt(1), 'sec' )
      ! interpolation weight:
      alfa = min( max( 0.0, rTotal( self%tt(2) - t, 'sec' ) / dt_sec ), 1.0 )
      ! update data:
      self%data = self%data1 * alfa + self%data2 * (1.0-alfa)

    end if ! copy or interpolate
    
    ! set flag:
    self%with_data = .true.
    ! update time:
    self%t = t

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Interpol
  

  ! ***


  subroutine Instant_Field_Series_Setup( self, t, status )
  
    use GO_Date, only : operator(/=), operator(<=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target time:
    self%setup_t = t
    
    ! no data present yet ?
    if ( .not. self%with_data ) then
      !! info ...
      !write (gol,'(a,": no data filled yet")') rname; call goPr
      ! new data required ...
      status = -1; return
    else
      ! data defined; not valid for this time ?
      if ( self%t /= t ) then
        !! info ...
        !write (gol,'(a,": data not valid for target time")') rname; call goPr
        ! tendency not defined yet ?
        if ( .not. self%with_data2 ) then
          !! info ...
          !write (gol,'(a,": no tendency data yet; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        else if ( (self%tt(1) <= t) .and. (t <= self%tt(2)) ) then
          !! info ...
          !call wrtgol( rname//': target time ', t, ' in interpolation interval ', self%tt ); call goPr
          ! interpolate:
          call self%Interpol( t, status )
          IF_NOTOK_RETURN(status=1)
        else 
          !! info ...
          !write (gol,'(a,": target time outside interpolation interval; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        end if
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup
  

  ! ***


  subroutine Instant_Field_Series_Setup_Prepare( self, tref, nreceive, status, reset_levels )
  
    use GO_Date, only : TDate
    use GO_Date, only : operator(==), operator(/=), operator(<), operator(<=)
    use GO_Date, only : wrtgol
    use GO_Date, only : Get_Surrounding_Interval
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    type(TDate), intent(in)                         ::  tref
    integer, intent(out)                            ::  nreceive
    integer, intent(out)                            ::  status
    logical, intent(in), optional                   ::  reset_levels

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_Prepare'
    
    ! --- local ----------------------------------
    
    type(TDate)   ::  tt(2)
    logical       ::  do_reset_levels

    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': prepare setup for ', self%setup_t ); call goPr
    
    ! prepare reset?
    do_reset_levels = .false.
    if ( present(reset_levels) ) do_reset_levels = reset_levels
    
    ! reset 3d data dimension?
    if ( do_reset_levels ) then
      ! clear data, reset time values:
      call Instant_Field_Series_Reset( self, status )
      IF_NOTOK_RETURN(status=1)
      ! reset 3rd dimension:
      self%ubo(3) = -999
      self%shp(3) = self%ubo(3) - self%lbo(3) + 1
    end if
    
    ! some fields need to be changed; 
    ! get interpolation points given temporal resolution:
    call Get_Surrounding_Interval( self%setup_t, &
                                   tref, self%interpolation, &
                                   tt, status )
    IF_NOTOK_RETURN(status=1)
    !! info ...
    !call wrtgol( rname//': surrounding interval ', tt ); call goPr

    ! single interpolation point (at target time) or interval ?
    !~ single point
    if ( tt(1) == tt(2) ) then
      ! already data present ?
      if ( self%with_data ) then
        ! try if current data is already valid for target time;
        ! this routine should not have been called ...
        if ( self%t == self%setup_t ) then
          call wrtgol( 'no need to prepare setup, current data already valid for ', self%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! idem for interpolation:
        if ( self%with_data2 .and. (self%tt(1) <= self%setup_t) .and. (self%setup_t <= self%tt(2)) ) then
          call wrtgol( 'no need to prepare setup, interpolation interval includes ', self%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! receive a single field valid for the target time, use it for the base:
      self%nreceive = 1
      self%receive_time(self%nreceive) = tt(1)
      self%receive_targ(self%nreceive) = INSTF_DATA1
    !~ interval, interpolate to target time
    else

      !~ no first data present yet ?
      if ( .not. self%with_data1 ) then
        ! retreive first field:
        self%nreceive = 1
        self%receive_time(self%nreceive) = tt(1)
        self%receive_targ(self%nreceive) = INSTF_DATA1
      !~ data1 present; valid for target time already?
      else if ( self%tt(1) == tt(1) ) then
        ! no field to be received for this:
        self%nreceive = 0
      !~ second field available?
      else if ( self%with_data2 ) then
        ! expected to be valid ...
        if ( tt(1) /= self%tt(2) ) then
          call wrtgol( 'start of target interval ', tt ); call goPr
          call wrtgol( 'expected to be end of current interval ', self%tt ); call goPr
          TRACEBACK; status=1;
        end if
        ! copy:
        self%data1 = self%data2
        ! no field to be received for this:
        self%nreceive = 0
        self%tt(1) = self%tt(2)
      !~ otherwise receive new field:
      else
        ! retreive first field:
        self%nreceive = 1
        self%receive_time(self%nreceive) = tt(1)
        self%receive_targ(self%nreceive) = INSTF_DATA1
      end if

      !~ no second data present yet ?
      if ( .not. self%with_data2 ) then
        ! retreive field for data2:
        self%nreceive = self%nreceive + 1
        self%receive_time(self%nreceive) = tt(2)
        self%receive_targ(self%nreceive) = INSTF_DATA2
      !~ data2 present; valid for target time already?
      else if ( self%tt(2) == tt(2) ) then
        ! nothing to be done ...
      !~ otherwise receive new field:
      else
        ! retreive field for data2:
        self%nreceive = self%nreceive + 1
        self%receive_time(self%nreceive) = tt(2)
        self%receive_targ(self%nreceive) = INSTF_DATA2
      end if

    end if  ! single field or interval
    
    ! return value:
    nreceive = self%nreceive

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_Prepare
  

  ! ***


  subroutine Instant_Field_Series_Setup_InqTime( self, ip, t, status )
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(in)       ::  self
    integer, intent(in)                             ::  ip
    type(TDate), intent(out)                        ::  t
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_InqTime'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract:
    t = self%receive_time(ip)
    
    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_InqTime
  

  ! ***

  !
  ! Return status:
  !  -1 : shape of data does not match field description ;
  !       this might be used to reset the field with 'reset_levels'
  !       argument to 'Setup_Prepare'

  subroutine Instant_Field_Series_Setup_Put( self, ip, data, t, status, lbo, ubo )
  
    use GO_Date, only : rTotal, operator(-), operator(==), operator(/=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    integer, intent(in)                             ::  ip
    real, intent(in)                                ::  data(:,:,:)
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status
    integer, intent(in), optional                   ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Setup_Put'

    ! --- local ----------------------------------

    real        ::  dt_sec
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') ip, self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! allocate data if not done yet, setup index space if necessary:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
  
    ! check ...
    if ( any(shape(data) /= self%shp) ) then
      write (gol,'("shape of data does not match with definition:")'); call goErr
      write (gol,'("  data     : ",3i6)') shape(data); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=-1; return
    end if

    ! check ...
    if ( t /= self%receive_time(ip) ) then
      write (gol,'("time of received field ",i6," does not match with setup:")') ip; call goErr
      call wrtgol( '  setup time    : ', self%receive_time(ip) ); call goErr
      call wrtgol( '  received time : ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! what to do with this data ?
    select case ( self%receive_targ(ip) )

      !~~ store as first:
      case ( INSTF_DATA1 )
        ! setup tendency if necessary:
        if ( .not. self%with_data1 ) then
          ! storage:
          allocate( self%data1(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! store:
        self%data1 = data
        ! store time:
        self%tt(1) = t
        ! set flag:
        self%with_data1 = .true.
        !! info ...
        !call wrtgol( rname//': stored data1 valid for ', self%tt(1) ); call goPr

      !~~ store as second:
      case ( INSTF_DATA2 )
        ! setup tendency if necessary:
        if ( .not. self%with_data2 ) then
          ! storage:
          allocate( self%data2(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! store:
        self%data2 = data
        ! store time:
        self%tt(2) = t
        ! set flag:
        self%with_data2 = .true.
        !! info ...
        !call wrtgol( rname//': stored data2 valid for ', self%tt(2) ); call goPr

      !~~ strange ...
      case default
        ! something went wrong; bug ?
        write (gol,'("unsupported destination code ",i6," for interpolation point ",i6)') self%receive_targ(ip), ip; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! reset destination to dummy value:
    self%receive_targ(ip) = -1

    ! last one received?
    if ( ip == self%nreceive ) then
      ! interpolate:
      call self%Interpol( self%setup_t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Setup_Put


  ! ***

  subroutine Instant_Field_Series_Get_Tendency( self, tr, pdata, status )

    use GO_Date, only : TDate, operator(<), wrtgol
    use GO_Date, only : rTotal, operator(-)

    ! --- in/out ---------------------------------
    
    class(T_Instant_Field_Series), intent(inout)    ::  self
    type(TDate), intent(in)                         ::  tr(2)
    real, pointer                                   ::  pdata(:,:,:)
    integer, intent(out)                            ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_Field_Series_Get_Tendency'
    
    ! --- local ----------------------------------

    real        ::  dt_sec
      
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%with_data2 ) then
      write (gol,'("no tendency information stored yet (first step of time series?)")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check time ..
    if ( (tr(1) < self%tt(1)) .or. (self%tt(2) < tr(2)) ) then
      write (gol,'("requested time interval is not within current interpolation step:")'); call goErr
      call wrtgol( '  interpolation step : ', self%tt ); call goErr
      call wrtgol( '  requested tendency : ', tr ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! setup tendency if necessary:
    if ( .not. associated(self%tend) ) then
      ! storage:
      allocate( self%tend(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! timestep in seconds:
    dt_sec = rTotal( self%tt(2) - self%tt(1), 'sec' )
    ! fill:
    self%tend = ( self%data2 - self%data1 )/dt_sec
    ! set pointer:
    pdata => self%tend

    ! ok
    status = 0
    
  end subroutine Instant_Field_Series_Get_Tendency
  

  ! ====================================================================
  ! ===
  ! === VectorField
  ! ===
  ! ====================================================================


  subroutine VectorField_Init( self, name, units, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(out)             ::  self
    character(len=*), intent(in)                  ::  name, units
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_Init'
      
    ! --- begin ----------------------------------
    
    ! init base:
    call FieldBase_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)
    
    ! no data yet:
    nullify( self%udata )
    nullify( self%vdata )

    ! ok
    status = 0
    
  end subroutine VectorField_Init


  ! ***


  subroutine VectorField_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(inout)           ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_Done'
      
    ! --- begin ----------------------------------
    
    ! no data:
    if ( associated(self%udata) ) then
      deallocate( self%udata, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    if ( associated(self%vdata) ) then
      deallocate( self%vdata, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! done with base:
    call FieldBase_Done( self, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine VectorField_Done


  ! ***


  subroutine VectorField_Reset( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(inout)           ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_Reset'
      
    ! --- begin ----------------------------------
    
    ! no data:
    if ( associated(self%udata) ) then
      deallocate( self%udata, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    if ( associated(self%vdata) ) then
      deallocate( self%vdata, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! done with base:
    call FieldBase_Reset( self, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine VectorField_Reset


  ! ***


  subroutine VectorField_Get( self, status, &
                         enabled, filled, &
                         name, units, input, &
                         lbo, ubo, shp, &
                         pudata, pvdata )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(in)              ::  self
    integer, intent(out)                          ::  status
    logical, intent(out), optional                ::  enabled
    logical, intent(out), optional                ::  filled
    character(len=*), intent(out), optional       ::  name, units, input
    integer, intent(out), optional                ::  lbo(3), ubo(3), shp(3)
    real, pointer, optional                       ::  pudata(:,:,:)
    real, pointer, optional                       ::  pvdata(:,:,:)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_Get'
      
    ! --- begin ----------------------------------
    
    ! get from base:
    call FieldBase_Get( self, status, &
                         enabled=enabled, filled=filled, &
                         name=name, units=units, input=input, &
                         lbo=lbo, ubo=ubo, shp=shp )
    IF_NOTOK_RETURN(status=1)
    
    ! check ..
    if ( any( (/ present(pudata), present(pvdata) /) ) ) then
      ! check ..
      if ( .not. self%enabled ) then
        write (gol,'("requested properties of `",a,"` VectorField while not enabled")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! set pointers:
      if ( present(pudata) ) pudata => self%udata
      if ( present(pvdata) ) pvdata => self%vdata
    end if

    ! ok
    status = 0
    
  end subroutine VectorField_Get


  ! ***
  

  subroutine VectorField_Alloc( self, status, lbo, ubo )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(inout) ::  self
    integer, intent(out)                ::  status
    integer, intent(in), optional       ::  lbo(3), ubo(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_Alloc'
      
    ! --- local ----------------------------------
    
    integer         ::  i
      
    ! --- begin ----------------------------------
    
    ! setup shp attributes:
    call FieldBase_SetupShape( self, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! need to allocate ?
    if ( .not. associated(self%udata) ) then
      ! storage:
      allocate( self%udata(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! need to allocate ?
    if ( .not. associated(self%vdata) ) then
      ! storage:
      allocate( self%vdata(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine VectorField_Alloc
  

  ! ***
  
  
  ! allocate external arrays with same index range
  
  subroutine VectorField_AllocExt( self, udata, vdata, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(inout)           ::  self
    real, intent(inout), allocatable              ::  udata(:,:,:)
    real, intent(inout), allocatable              ::  vdata(:,:,:)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VectorField_AllocExt'
      
    ! --- begin ----------------------------------
    
    ! enabled ?
    if ( .not. self%enabled ) then
      write (gol,'("VectorField `",a,"` not enabled")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( udata(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( vdata(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine VectorField_AllocExt
  

  ! ***


  subroutine VectorField_Put( self, udata, vdata, status, lbo, ubo )
  
    ! --- in/out ---------------------------------
    
    class(T_VectorField), intent(inout) ::  self
    real, intent(in)                    ::  udata(:,:,:)
    real, intent(in)                    ::  vdata(:,:,:)
    integer, intent(out)                ::  status
    integer, intent(in), optional       ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/VectorField_Put'

    ! --- local ----------------------------------
    
    !integer     ::  i, j, k
    
    ! --- begin ----------------------------------
    
    ! allocate data if necessary, setup index space if not done yet:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! check ...
    if ( any(shape(udata) /= self%shp) ) then
      write (gol,'("shape of udata does not match with definition:")'); call goErr
      write (gol,'("  udata    : ",3i6)') shape(udata); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=1; return
    end if
    if ( any(shape(vdata) /= self%shp) ) then
      write (gol,'("shape of vdata does not match with definition:")'); call goErr
      write (gol,'("  vdata    : ",3i6)') shape(vdata); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=1; return
    end if

    ! store data:
    self%udata = udata
    self%vdata = vdata
    ! set flag:
    self%with_data = .true.

    !! testing ...
    !write (gol,*) 'fff ', trim(self%name),'(5,28,1) = ', self%udata(5,28,1), self%vdata(5,28,1); call goPr

    !! check ..
    !if ( any(isnan(self%udata)) .or. any(isnan(self%vdata)) ) then
    !  write (gol,'("found NaN values in `",a,"`")') trim(self%name); call goErr
    !  do k = self%lbo(3), self%ubo(3)
    !    do j = self%lbo(2), self%ubo(2)
    !      do i = self%lbo(1), self%ubo(1)
    !        if ( isnan(self%udata(i,j,k)) .or. isnan(self%vdata(i,j,k)) ) then
    !          write (gol,'("  cell ",3i6)') i, j, k; call goErr
    !        end if
    !      end do ! i
    !    end do ! j
    !  end do ! k
    !  TRACEBACK; status=1; return
    !end if

    ! ok
    status = 0
    
  end subroutine VectorField_Put
  


  

  ! ====================================================================
  ! ===
  ! === Instant_VectorField
  ! ===
  ! ====================================================================


  subroutine Instant_VectorField_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField), intent(out)   ::  self
    character(len=*), intent(in)                ::  name, units
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Init'

    ! --- begin ----------------------------------
    
    ! initialize the instant field:
    call VectorField_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no time assigned yet:
    self%t = AnyDate()

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Init
  

  ! ***


  subroutine Instant_VectorField_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField), intent(inout)    ::  self
    integer, intent(out)                           ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Done'
      
    ! --- begin ----------------------------------
    
    ! no current data:
    self%t = AnyDate()
    
    ! done with field:
    call VectorField_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Done


  ! ***


  subroutine Instant_VectorField_Reset( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField), intent(inout)    ::  self
    integer, intent(out)                           ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Reset'
      
    ! --- begin ----------------------------------
    
    ! no current data:
    self%t = AnyDate()
    
    ! done with field:
    call VectorField_Reset( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Reset

  !
  ! Check data in field:
  !   enabled   : compare enabled flag 
  !   t         : time
  ! For messages:
  !   silent    : do not raise error messages
  ! Return status:
  !   <1 if check failed

  subroutine Instant_VectorField_CheckSample( self, status, enabled, t, silent )
  
    use GO_Date, only : TDate, operator(/=), wrtgol
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField), intent(inout)   ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
    type(TDate), intent(in), optional             ::  t
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_CheckSample'
      
    ! --- local ----------------------------------
    
    logical       ::  shout
      
    ! --- begin ----------------------------------
    
    ! original:
    call FieldBase_Check( self, status, enabled=enabled, silent=silent )
    ! leave with warning status ?
    if ( status < 0 ) return
    ! failed:
    IF_NOTOK_RETURN(status=1)
    
    ! shout?
    shout = .true.
    if ( present(silent) ) shout = .not. silent
    
    ! check time ?
    if ( present(t) ) then
      ! filled at all?
      if ( .not. self%with_data ) then
        ! not as expected ...
        if ( shout ) then
          write (gol,'("no data filled yet")'); call goErr
          TRACEBACK
        end if
        ! return with warning status:
        status=-1; return
      else
        ! compare:
        if ( self%t /= t ) then
          ! not as expected ...
          if ( shout ) then
            call wrtgol( 'expected time ', t, ' but sample valid for ', self%t ); call goErr
            TRACEBACK
          end if
          ! return with warning status:
          status=-1; return
        end if ! compare
      end if ! filled?
    end if ! present t ?

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_CheckSample
  

  ! ***


  subroutine Instant_VectorField_PutSample( self, udata, vdata, t, status, lbo, ubo )
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField), intent(inout)   ::  self
    real, intent(in)                              ::  udata(:,:,:)
    real, intent(in)                              ::  vdata(:,:,:)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
    integer, intent(in), optional                 ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_PutSample'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store data, use index space if not defined yet:
    call self%Put( udata, vdata, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! store time value:
    self%t = t

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_PutSample
  


  ! ====================================================================
  ! ===
  ! === Instant_VectorField_Series
  ! ===
  ! ====================================================================


  subroutine Instant_VectorField_Series_Init( self, name, units, status )
  
    use GO_String, only : goVarValue
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(out)  ::  self
    character(len=*), intent(in)                      ::  name, units
    integer, intent(out)                              ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Init'

    ! --- begin ----------------------------------
    
    ! initialize the instant field:
    call Instant_VectorField_Init( self, name, units, status )
    IF_NOTOK_RETURN(status=1)

    ! no interpolation points yet:
    self%with_data1 = .false.
    self%with_data2 = .false.
    
    ! no interpolation set yet:
    self%interpolation_step = -999
    self%interpolation_step_units = 'None'
    self%interpolation = 'None'

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Init
  

  ! ***


  subroutine Instant_VectorField_Series_Done( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    integer, intent(out)                                  ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Done'
      
    ! --- begin ----------------------------------
    
    ! no interpolation data anymore:
    if ( self%with_data1 ) then
      deallocate( self%udata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%udata_, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata_, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata_, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    if ( self%with_data2 ) then
      deallocate( self%udata2, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata2, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata2, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    
    ! done with field:
    call Instant_VectorField_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Done
  

  ! ***


  subroutine Instant_VectorField_Series_Reset( self, status )
  
    use GO_Date, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    integer, intent(out)                                  ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Reset'
      
    ! --- begin ----------------------------------
    
    ! no interpolation data anymore:
    if ( self%with_data1 ) then
      deallocate( self%udata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata1, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%udata_, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata_, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata_, stat=status )
      IF_NOTOK_RETURN(status=1)
      self%with_data1 = .false.
    end if
    if ( self%with_data2 ) then
      deallocate( self%udata2, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%vdata2, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%ldata2, stat=status )
      IF_NOTOK_RETURN(status=1)
      self%with_data2 = .false.
    end if
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    
    
    ! reset field:
    call Instant_VectorField_Reset( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Reset
  

  ! ***
  ! ***


  subroutine Instant_VectorField_Series_Set_Interp( self, interp, status )
  
    use GO_String, only : goReadFromLine
    use GO_Date  , only : AnyDate
    
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)  ::  self
    character(len=*), intent(in)                        ::  interp
    integer, intent(out)                                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Set_Interp'

    ! --- local ----------------------------------
    
    character(len=64)        ::    line

    ! --- begin ----------------------------------
    
    ! save:
    self%interpolation = trim(interp)
    
    ! copy:
    line = trim(interp)
    
    ! format:  3 hour
    ! read time step:
    call goReadFromLine( line, self%interpolation_step, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! rest is step units:
    self%interpolation_step_units = trim(line)
    
    ! no tendency data yet:
    self%tt(1) = AnyDate()
    self%tt(2) = AnyDate()
    !self%with_tend = .false.
    self%with_data1 = .false.
    self%with_data2 = .false.

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Set_Interp
  

  ! ***


  ! Interpolate to t given available data and tendency.
  ! Tool for use within this module only.
  
  subroutine Instant_VectorField_Series_Interpol( self, t, status )
  
    use GO_Date, only : operator(==), operator(/=), operator(<)
    use GO_Date, only : operator(-)
    use GO_Date, only : rTotal
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    type(TDate), intent(in)                               ::  t
    integer, intent(out)                                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Interpol'

    ! --- local ----------------------------------
    
    real                  ::  dt_sec
    real                  ::  alfa
    real, allocatable     ::  len_(:,:,:)

    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': interpolate between [', self%tt(1), ',', self%tt(2), '] to ', t ); call goPr
    
    ! interval bounds?
    if ( self%with_data1 .and. (t == self%tt(1)) ) then
    
      ! copy data:
      self%udata_ = self%udata1
      self%vdata_ = self%vdata1
      ! copy length:
      self%ldata_ = self%ldata1
    
    else
    
      ! check ...
      if ( (.not. self%with_data1) .or. (.not. self%with_data2) ) then
        write (gol,'("not all data present: data1 (",l1,"), data2 (",l1,")")') &
                self%with_data1, self%with_data2; call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( (t < self%tt(1)) .or. (self%tt(2) < t) ) then
        call wrtgol( 'target time ', t ); call goErr
        call wrtgol( 'outside data interval ', self%tt ); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! data defined by start and end;
      ! interpolate:  data(t) = data(t1) * (t2-t)/(t2-t1) + ddata(self%t) * (t-t1)/(t2-t1)
      !                                        alfa                            (1-alfa)
      ! total timestep in seconds:
      dt_sec = rTotal( self%tt(2) - self%tt(1), 'sec' )
      ! interpolation weight:
      alfa = min( max( 0.0, rTotal( self%tt(2) - t, 'sec' ) / dt_sec ), 1.0 )
      ! interpolate data:
      self%udata_ = self%udata1 * alfa + self%udata2 * (1.0-alfa)
      self%vdata_ = self%vdata1 * alfa + self%vdata2 * (1.0-alfa)
      ! interpolate length:
      self%ldata_ = self%ldata1 * alfa + self%ldata2 * (1.0-alfa)

    end if ! copy or interpolate

    ! storage for current length:
    allocate( len_(self%shp(1),self%shp(2),self%shp(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! compute current length:
    len_ = sqrt( self%udata_**2 + self%vdata_**2 )
    ! scale vectors to interpolated length:
    where ( len_ > 0.0 )
      self%udata = self%udata_ * self%ldata_/len_
      self%vdata = self%vdata_ * self%ldata_/len_
    end where
    ! set flag:
    self%with_data = .true.
    ! update time:
    self%t = t

    ! clear:
    deallocate( len_, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Interpol
  

  ! ***


  subroutine Instant_VectorField_Series_Setup( self, t, status )
  
    use GO_Date, only : operator(/=), operator(<=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    type(TDate), intent(in)                               ::  t
    integer, intent(out)                                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Setup'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target time:
    self%setup_t = t
    
    ! no data present yet ?
    if ( .not. self%with_data ) then
      !! info ...
      !write (gol,'(a,": no data filled yet")') rname; call goPr
      ! new data required ...
      status = -1; return
    else
      ! data defined; not valid for this time ?
      if ( self%t /= t ) then
        !! info ...
        !write (gol,'(a,": data not valid for target time")') rname; call goPr
        ! tendency not defined yet ?
        if ( .not. self%with_data2 ) then
          !! info ...
          !write (gol,'(a,": no tendency data yet; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        else if ( (self%tt(1) <= t) .and. (t <= self%tt(2)) ) then
          !! info ...
          !call wrtgol( rname//': target time ', t, ' in interpolation interval ', self%tt ); call goPr
          ! interpolate:
          call Instant_VectorField_Series_Interpol( self, t, status )
          IF_NOTOK_RETURN(status=1)
        else 
          !! info ...
          !write (gol,'(a,": target time outside interpolation interval; need to receive new data")') rname; call goPr
          ! new data required:
          status = -1; return
        end if
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Setup
  

  ! ***


  subroutine Instant_VectorField_Series_Setup_Prepare( self, tref, nreceive, status, reset_levels )
  
    use GO_Date, only : TDate
    use GO_Date, only : operator(==), operator(/=), operator(<), operator(<=)
    use GO_Date, only : wrtgol
    use GO_Date, only : Get_Surrounding_Interval
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    type(TDate), intent(in)                               ::  tref
    integer, intent(out)                                  ::  nreceive
    integer, intent(out)                                  ::  status
    logical, intent(in), optional                         ::  reset_levels

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Setup_Prepare'
    
    ! --- local ----------------------------------
    
    type(TDate)   ::  tt(2)
    logical       ::  do_reset_levels

    ! --- begin ----------------------------------
    
    !! info ...
    !call wrtgol( rname//': prepare setup for ', self%setup_t ); call goPr
    
    ! prepare reset?
    do_reset_levels = .false.
    if ( present(reset_levels) ) do_reset_levels = reset_levels
    
    ! reset 3d data dimension?
    if ( do_reset_levels ) then
      ! clear data, reset time values:
      call Instant_VectorField_Series_Reset( self, status )
      IF_NOTOK_RETURN(status=1)
      ! reset 3rd dimension:
      self%ubo(3) = -999
      self%shp(3) = self%ubo(3) - self%lbo(3) + 1
    end if
    
    ! some fields need to be changed; 
    ! get interpolation points given temporal resolution:
    call Get_Surrounding_Interval( self%setup_t, &
                                   tref, self%interpolation, &
                                   tt, status )
    IF_NOTOK_RETURN(status=1)
    !! info ...
    !call wrtgol( rname//': surrounding interval ', tt ); call goPr

    ! single interpolation point (at target time) or interval ?
    !~ single point
    if ( tt(1) == tt(2) ) then
      ! already data present ?
      if ( self%with_data ) then
        ! try if current data is already valid for target time;
        ! this routine should not have been called ...
        if ( self%t == self%setup_t ) then
          call wrtgol( 'no need to prepare setup, current data already valid for ', self%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! idem for interpolation:
        if ( self%with_data2 .and. (self%tt(1) <= self%setup_t) .and. (self%setup_t <= self%tt(2)) ) then
          call wrtgol( 'no need to prepare setup, interpolation interval includes ', self%setup_t ); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! receive a single field valid for the target time, use it for the base:
      self%nreceive = 1
      self%receive_time(self%nreceive) = tt(1)
      !self%receive_targ(self%nreceive) = INSTF_BASE
      self%receive_targ(self%nreceive) = INSTF_DATA1
    !~ interval, interpolate to target time
    else
      !~ no first data present yet ?
      if ( .not. self%with_data1 ) then
        ! retreive first field:
        self%nreceive = 1
        self%receive_time(self%nreceive) = tt(1)
        self%receive_targ(self%nreceive) = INSTF_DATA1
      !~ data1 present; valid for target time already?
      else if ( self%tt(1) == tt(1) ) then
        ! no field to be received for this:
        self%nreceive = 0
      !~ second field available?
      else if ( self%with_data2 ) then
        ! expected to be valid ...
        if ( tt(1) /= self%tt(2) ) then
          call wrtgol( 'start of target interval ', tt ); call goPr
          call wrtgol( 'expected to be end of current interval ', self%tt ); call goPr
          TRACEBACK; status=1;
        end if
        ! copy:
        self%udata1 = self%udata2
        self%vdata1 = self%vdata2
        self%ldata1 = self%ldata2
        self%tt(1) = self%tt(2)
        ! no field to be received for this:
        self%nreceive = 0
      !~ otherwise receive new field:
      else
        ! retreive first field:
        self%nreceive = 1
        self%receive_time(self%nreceive) = tt(1)
        self%receive_targ(self%nreceive) = INSTF_DATA1
      end if

      !~ no second data present yet ?
      if ( .not. self%with_data2 ) then
        ! retreive field for data2:
        self%nreceive = self%nreceive + 1
        self%receive_time(self%nreceive) = tt(2)
        self%receive_targ(self%nreceive) = INSTF_DATA2
      !~ data2 present; valid for target time already?
      else if ( self%tt(2) == tt(2) ) then
        ! nothing to be done ...
      !~ otherwise receive new field:
      else
        ! retreive field for data2:
        self%nreceive = self%nreceive + 1
        self%receive_time(self%nreceive) = tt(2)
        self%receive_targ(self%nreceive) = INSTF_DATA2
      end if

    end if  ! single field or interval
    
    ! return value:
    nreceive = self%nreceive

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Setup_Prepare
  

  ! ***


  subroutine Instant_VectorField_Series_Setup_InqTime( self, ip, t, status )
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(in)   ::  self
    integer, intent(in)                               ::  ip
    type(TDate), intent(out)                          ::  t
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Setup_InqTime'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract:
    t = self%receive_time(ip)
    
    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Setup_InqTime
  

  ! ***


  subroutine Instant_VectorField_Series_Setup_Put( self, ip, udata, vdata, t, status, lbo, ubo )
  
    use GO_Date, only : rTotal, operator(-), operator(==), operator(/=)
    use GO_Date, only : wrtgol
 
    ! --- in/out ---------------------------------
    
    class(T_Instant_VectorField_Series), intent(inout)    ::  self
    integer, intent(in)                                   ::  ip
    real, intent(in)                                      ::  udata(:,:,:)
    real, intent(in)                                      ::  vdata(:,:,:)
    type(TDate), intent(in)                               ::  t
    integer, intent(out)                                  ::  status
    integer, intent(in), optional                         ::  lbo(3), ubo(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Instant_VectorField_Series_Setup_Put'

    ! --- local ----------------------------------

    real        ::  dt_sec
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ip < 1) .or. (ip > self%nreceive) ) then
      write (gol,'("index ",i6," of interpolation point outside bounds 1 .. ",i6)') ip, self%nreceive; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! allocate data if not done yet, setup index space if necessary:
    call self%Alloc( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
  
    ! check ...
    if ( any(shape(udata) /= self%shp) ) then
      write (gol,'("shape of udata does not match with definition:")'); call goErr
      write (gol,'("  udata    : ",3i6)') shape(udata); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=-1; return
    end if
    if ( any(shape(vdata) /= self%shp) ) then
      write (gol,'("shape of vdata does not match with definition:")'); call goErr
      write (gol,'("  vdata    : ",3i6)') shape(vdata); call goErr
      write (gol,'("  defined  : ",3i6)') self%shp; call goErr
      TRACEBACK; status=-1; return
    end if

    ! check ...
    if ( t /= self%receive_time(ip) ) then
      write (gol,'("time of received field ",i6," does not match with setup:")') ip; call goErr
      call wrtgol( '  setup time    : ', self%receive_time(ip) ); call goErr
      call wrtgol( '  received time : ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! what to do with this data ?
    select case ( self%receive_targ(ip) )

      !~~ store as first:
      case ( INSTF_DATA1 )
        ! setup tendency if necessary:
        if ( .not. self%with_data1 ) then
          ! storage for first data:
          allocate( self%udata1(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%vdata1(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%ldata1(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! storage for interpolation:
          allocate( self%udata_(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%vdata_(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%ldata_(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! store:
        self%udata1 = udata
        self%vdata1 = vdata
        ! compute length:
        self%ldata1 = sqrt( self%udata1**2 + self%vdata1**2 )
        ! store time:
        self%tt(1) = t
        ! set flag:
        self%with_data1 = .true.
        !! info ...
        !call wrtgol( rname//': stored data1 valid for ', self%tt(1) ); call goPr

      !~~ store as second:
      case ( INSTF_DATA2 )
        ! setup tendency if necessary:
        if ( .not. self%with_data2 ) then
          ! storage for second data:
          allocate( self%udata2(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%vdata2(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
          allocate( self%ldata2(self%lbo(1):self%ubo(1),self%lbo(2):self%ubo(2),self%lbo(3):self%ubo(3)), stat=status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! store:
        self%udata2 = udata
        self%vdata2 = vdata
        ! compute length:
        self%ldata2 = sqrt( self%udata2**2 + self%vdata2**2 )
        ! store time:
        self%tt(2) = t
        ! set flag:
        self%with_data2 = .true.
        !! info ...
        !call wrtgol( rname//': stored data2 valid for ', self%tt(2) ); call goPr

      !~~ strange ...
      case default
        ! something went wrong; bug ?
        write (gol,'("unsupported destination code ",i6," for interpolation point ",i6)') self%receive_targ(ip), ip; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! reset destination to dummy value:
    self%receive_targ(ip) = -1

    ! last one received?
    if ( ip == self%nreceive ) then
      ! interpolate:
      call self%Interpol( self%setup_t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine Instant_VectorField_Series_Setup_Put
  

end module GO_Field



!!#######################################################################
!!###
!!### Test
!!###
!!#######################################################################
!!
!#define if (status/=0) stop if (status/=0) then; TRACEBACK; stop; end if
!#define IF_ERROR_STOP if (status >0) then; TRACEBACK; stop; end if
!!
!! f90 -o test.x go_fu.F90 go_print.F90 go_string.F90 go_date.F90 go_field.F90
!
!program test
!
!  use GO_Print , only : gol, goPr, goErr
!  use GO_Date  , only : TDate, NewDate, wrtgol, IsAnyDate
!  use GO_Fields
!  
!  implicit none
!  
!  ! dimensions:
!  integer, parameter     ::  nx = 4, ny = 3, nz = 2
!
!  ! data:
!  class(T_Instant_Field_Series)  ::  self
!  integer                        ::  hour
!  type(TDate)                    ::  t
!  integer                        ::  nreceive
!  type(TDate)                    ::  receive_time
!  integer                        ::  ireceive
!  real                           ::  data_in(nx,ny,nz)
!  integer                        ::  status
!  
!  ! info ...
!  write (gol,'("test:")'); call goPr
!  write (gol,'("test: testing GO_Fields module")'); call goPr
!  write (gol,'("test:")'); call goPr
!
!  ! Initialize 3D field; provide name and units,
!  ! and description of temporal interpolation.
!  ! optionally provide shape or lower/upper bounds:
!  call Instant_Field_Series_Init( self, 'temperature', 'K', shape(data_in), &
!                 'interpolation=linear;step=3;units=hour', status )
!  if (status/=0) stop
!  
!  ! loop over hours:
!  do hour = 0, 6
!  
!    write (gol,'("test:")'); call goPr
!    write (gol,'("test: >>> hour ",i2)') hour; call goPr
!    write (gol,'("test:")'); call goPr
!
!    ! current time:
!    t = NewDate( 2011, 03, 31, hour, 0 )
!
!    ! info ...
!    call wrtgol( 'test: target time : ', t ); call goPr
!
!    ! info ...
!    write (gol,'("test: setup ...")'); call goPr
!    ! setup for time t, or return status -1 if new data should be received:
!    call Instant_Field_Series_Setup( self, t, status )
!    IF_ERROR_STOP
!    if ( status < 0 ) then
!      ! info ...
!      write (gol,'("test: something to be done ...")'); call goPr
!      ! interpolate, or prepare to receive new data:
!      call Instant_Field_Series_Setup_Prepare( self, nreceive, status )
!      if (status/=0) stop
!      ! info ...
!      write (gol,'("test: number of fields to receive : ",i6)') nreceive; call goPr
!      ! might need to receive new data ...
!      do ireceive = 1, nreceive
!        ! get time value:
!        call Instant_Field_Series_Setup_InqTime( self, ireceive, receive_time, status )
!        if (status/=0) stop
!        ! obtain data from somewhere; here fill with hour plus decimal minutes ...
!        data_in = receive_time%hour + receive_time%min/60.0/10.0
!        ! info ...
!        call wrtgol( 'test: put field valid for ', receive_time ); call goPr
!        write (gol,'("test: average value : ",f12.4)') sum(data_in)/size(data_in); call goPr
!        ! store:
!        call Instant_Field_Series_Setup_Put( self, ireceive, data_in, receive_time, status )
!        if (status/=0) stop
!      end do
!    end if
!    ! info ...
!    call wrtgol( 'test: data valid for : ', self%t ); call goPr
!    write (gol,'("test: average value : ",f12.4)') sum(self%data)/size(self%data); call goPr
!    
!  end do
!
!  ! clear:
!  call Instant_Field_Series_Done( self, status )
!  if (status/=0) stop
!  
!  ! info ...
!  write (gol,'("test:")'); call goPr
!  write (gol,'("test: end.")'); call goPr
!  write (gol,'("test:")'); call goPr
!  
!end program test


