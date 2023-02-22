!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_Variable

  use GO, only : gol, goPr, goErr, goLabel
  use GO, only : TrcFile
  use GO, only : T_FieldBase
  use GO, only : T_Field
  use GO, only : T_Constant_Field, T_Constant_Field_Series, T_Accumulated_Field_Series
  use GO, only : T_Instant_Field, T_Instant_Field_Series
  use GO, only : T_VectorField, T_Instant_VectorField, T_Instant_VectorField_Series

  ! --- in/out -----------------------------------
  
  private

  public  ::  T_Variable

  public  ::  LEN_VARNAME
  public  ::  LEN_LEVTYPE
  public  ::  LEN_UNITS

  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_Variable'

  ! keyword lengths:
  integer, parameter          ::  LEN_VARNAME   = 32
  integer, parameter          ::  LEN_LONG_NAME = 256
  integer, parameter          ::  LEN_UNITS     = 64
  integer, parameter          ::  LEN_GRIDTYPE  = 16
  integer, parameter          ::  LEN_LEVTYPE   = 16
  integer, parameter          ::  LEN_DATATYPE  = 32
  integer, parameter          ::  LEN_TIME      = 1
  integer, parameter          ::  LEN_INPUT     = 1024
  integer, parameter          ::  LEN_ROUTINE   = 64
  
  ! maximum number of dependencies:
  integer, parameter          ::  maxdep = 10


  ! --- types ------------------------------------
  
  !
  ! Variable with gridded field.
  ! Contains type derived from 'T_Field' ;
  ! in future, shared entities might go into this type.
  !

  type T_Variable
    ! id:
    character(len=LEN_VARNAME)                ::  name
    character(len=LEN_LONG_NAME)              ::  long_name
    character(len=LEN_UNITS)                  ::  units
    ! truncation:
    logical                                   ::  with_minimum, with_maximum
    real                                      ::  minimum, maximum
    ! flag:
    logical                                   ::  enabled
    ! grid description:
    character(len=LEN_GRIDTYPE)               ::  gridtype
    character(len=LEN_GRIDTYPE)               ::  gridmapping
    character(len=LEN_LEVTYPE)                ::  levtype
    ! flags used for output:
    integer                                   ::  rank
    logical                                   ::  const
    ! data type for actual storage:
    character(len=LEN_DATATYPE)                    ::  datatype
    type(T_Field), pointer                         ::  F
    type(T_Constant_Field), pointer                ::  cF
    type(T_Constant_Field_Series), pointer         ::  cFS
    type(T_Accumulated_Field_Series), pointer      ::  aFS
    type(T_Instant_Field), pointer                 ::  iF
    type(T_Instant_Field_Series), pointer          ::  iFS
    type(T_VectorField), pointer                   ::  VF
    type(T_Instant_VectorField), pointer           ::  iVF
    type(T_Instant_VectorField_Series), pointer    ::  iVFS
    ! times:
    character(len=LEN_TIME)                   ::  time
    ! settings used while running:
    type(TrcFile), pointer                    ::  rcF
    ! processing:
    character(len=LEN_INPUT)                  ::  input
    character(len=LEN_ROUTINE)                ::  routine
    ! dependencies:
    integer                                   ::  ndep
    character(len=LEN_VARNAME)                ::  deps(maxdep)
    ! output variables:
    integer                                   ::  nout
    character(len=LEN_VARNAME)                ::  outs(maxdep)
  contains
    procedure :: Init          => Variable_Init
    procedure :: Done          => Variable_Done
    procedure :: Set           => Variable_Set
    procedure :: Get           => Variable_Get
    procedure :: Valid         => Variable_Valid
    procedure :: Setup         => Variable_Setup
    procedure :: Truncate      => Variable_Truncate
  end type T_Variable



contains


  ! ====================================================================
  ! ===
  ! === Variable
  ! ===
  ! ====================================================================


  subroutine Variable_Init( self, rcF, name, status )
  
    use GO            , only : TrcFile, ReadRc
    use GO            , only : goSplitString, goReadFromLine, goTranslate
    use GO            , only : T_FieldBase
    use LE_Grid       , only : ugg
    use LE_Data_Common, only : nlev, nlev_top
    
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(out)                ::  self
    type(TrcFile), pointer                        ::  rcF
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Init'
      
    ! --- local ----------------------------------
    
    character(len=64)                 ::  rcbase
    character(len=1024)               ::  line
    character(len=32)                 ::  part
    character(len=1024)               ::  names
    class(T_FieldBase), pointer       ::  FB
    integer                           ::  lbo(3), ubo(3)
    integer                           ::  k
      
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a,": init variable `",a,"`")') rname, trim(name); call goPr
    
    ! store:
    self%name = trim(name)
    
    ! base of rc keys:
    write (rcbase,'("data.",a)') trim(self%name)
    
    ! read setting:
    call ReadRc( rcF, trim(rcbase)//'.long_name', self%long_name, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read setting:
    call ReadRc( rcF, trim(rcbase)//'.units', self%units, status )
    IF_NOTOK_RETURN(status=1)
    
    ! range defined? 
    !   0.0 100.0
    !   0.0 Inf 
    !   -Inf 0.0
    call ReadRc( rcF, trim(rcbase)//'.range', line, status )
    IF_NOTOK_RETURN(status=1)
    ! first part:
    call goReadFromLine( line, part, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! special?
    if ( trim(part) == '-Inf' ) then
      ! no truncation:
      self%with_minimum = .false.
    else
      ! enable truncation:
      self%with_minimum = .true.
      ! extract threshold:
      call goReadFromLine( part, self%minimum, status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! second part:
    call goReadFromLine( line, part, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! special?
    if ( trim(part) == 'Inf' ) then
      ! no truncation:
      self%with_maximum = .false.
    else
      ! enable truncation:
      self%with_maximum = .true.
      ! extract threshold:
      call goReadFromLine( part, self%maximum, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! read setting:
    call ReadRc( rcF, trim(rcbase)//'.gridtype', self%gridtype, status )
    IF_NOTOK_RETURN(status=1)
    
    ! optional mapping:
    call ReadRc( rcF, trim(rcbase)//'.gridmapping', self%gridmapping, status, default='None' )
    IF_ERROR_RETURN(status=1)
    
    ! read setting:
    call ReadRc( rcF, trim(rcbase)//'.levtype', self%levtype, status )
    IF_NOTOK_RETURN(status=1)
    
    ! translate to rank:
    select case ( trim(self%levtype) )
      case ( 'sfc' )
        self%rank = 2
      case ( 'input_levels', 'input_halflevels', &
             'levels', 'halflevels', &
             'levels_top', 'halflevels_top' )
        self%rank = 3
      case default
        write (gol,'("could not set rank for levtype `",a,"`")') trim(self%levtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! read setting useful for output, by default false (output will have time dim):
    call ReadRc( rcF, trim(rcbase)//'.const', self%const, status, default=.false. )
    IF_ERROR_RETURN(status=1)
    
    ! read setting:
    call ReadRc( rcF, trim(rcbase)//'.datatype', self%datatype, status )
    IF_NOTOK_RETURN(status=1)
    ! init pointers:
    nullify( self%F    )
    nullify( self%cF   )
    nullify( self%cFS  )
    nullify( self%iF   )
    nullify( self%iFS  )
    nullify( self%iVF  )
    nullify( self%iVFS )
    ! init pointers:
    nullify( FB )
    ! switch:
    select case ( trim(self%datatype) )
      !
      case ( 'field' )
        ! storage:
        allocate( self%F, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%F%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB => self%F
      !
      case ( 'constant_field' )
        ! storage:
        allocate( self%cF, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%cF%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB => self%cF
      !
      case ( 'constant_field_series' )
        ! storage:
        allocate( self%cFS, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%cFS%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB => self%cFS
      !
      case ( 'accumulated_field_series' )
        ! storage:
        allocate( self%aFS, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%aFS%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB => self%aFS
      !
      case ( 'instant_field' )
        ! storage:
        allocate( self%iF, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%iF%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB  => self%iF
      !
      case ( 'instant_field_series' )
        ! storage:
        allocate( self%iFS, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%iFS%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB  => self%iFS
      !
      case ( 'vectorfield' )
        ! storage:
        allocate( self%VF, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%VF%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB  => self%VF
      !
      case ( 'instant_vectorfield' )
        ! storage:
        allocate( self%iVF, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%iVF%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB  => self%iVF
      !
      case ( 'instant_vectorfield_series' )
        ! storage:
        allocate( self%iVFS, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init:
        call self%iVFS%Init( trim(self%name), trim(self%units), status )
        IF_NOTOK_RETURN(status=1)
        ! pointer to ancestor:
        FB  => self%iVFS
      !
      case default
        ! unkown ...
        write (gol,'("unsupported datatype `",a,"`")') trim(self%datatype); call goErr
        TRACEBACK; status=1; return
      !
    end select
    
    ! read setting, defaults to current time:
    call ReadRc( rcF, trim(rcbase)//'.time', self%time, status, default='c' )
    IF_ERROR_RETURN(status=1)
    
    ! fill horizontal index space:
    select case ( trim(self%gridtype) )
      case ( 'cell', 'center' )
        lbo = (/        1,        1, 1 /)
        ubo = (/ ugg%nlon, ugg%nlat, 1 /)
      case ( 'corner' )
        lbo = (/        0,        0, 1 /)
        ubo = (/ ugg%nlon, ugg%nlat, 1 /)
      case ( 'u-edge' )
        lbo = (/        0,        1, 1 /)
        ubo = (/ ugg%nlon, ugg%nlat, 1 /)
      case ( 'v-edge' )
        lbo = (/        1,        0, 1 /)
        ubo = (/ ugg%nlon, ugg%nlat, 1 /)
      case ( 'cell_bnds' )
        lbo = (/          0,          0, 1 /)
        ubo = (/ ugg%nlon+1, ugg%nlat+1, 1 /)
      case ( 'corner_bnds' )
        lbo = (/         -1,         -1, 1 /)
        ubo = (/ ugg%nlon+1, ugg%nlat+1, 1 /)
      case default
        write (gol,'("unsupported gridtype `",a,"`")') trim(self%gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    ! set vertical index space:
    select case ( trim(self%levtype) )
      case ( 'sfc' )
        lbo(3) = 1
        ubo(3) = 1
      case ( 'levels' )
        lbo(3) = 1
        ubo(3) = nlev
      case ( 'halflevels' )
        lbo(3) = 0
        ubo(3) = nlev
      case ( 'levels_top' )
        lbo(3) = 1
        ubo(3) = nlev_top
      case ( 'halflevels_top' )
        lbo(3) = 0
        ubo(3) = nlev_top
      case ( 'input_levels' )
        lbo(3) = 1
        ubo(3) = -999
      case ( 'input_halflevels' )
        lbo(3) = 0
        ubo(3) = -999
      case default
        write (gol,'("could not set rank for levtype `",a,"`")') trim(self%levtype); call goErr
        TRACEBACK; status=1; return
    end select
    ! store:
    call FB%Set( status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! point to settings file, used at run time too:
    self%rcF => rcF
    
    ! read rckey with input definition, defaults to empty
    call ReadRc( rcF, trim(rcbase)//'.input', self%input, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! defined?
    if ( len_trim(self%input) > 0 ) then
      ! input description:
      call ReadRc( rcf, trim(self%input)//'.input', line, status )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call FB%Set( status, input=trim(line) )
      IF_NOTOK_RETURN(status=1)
      ! with temporal interpolation?
      if ( associated(self%aFS) .or. associated(self%iFS) .or. associated(self%iVFS) ) then
        ! read temporal interpolation description:
        call ReadRc( rcf, trim(self%input)//'.tinterp', line, status )
        IF_NOTOK_RETURN(status=1)
        ! store:
        if ( associated(self%aFS) ) then
          call self%aFS%Set_Interp( trim(line), status )
          IF_NOTOK_RETURN(status=1)
        else if ( associated(self%iFS) ) then
          call self%iFS%Set_Interp( trim(line), status )
          IF_NOTOK_RETURN(status=1)
        else if ( associated(self%iVFS) ) then
          call self%iVFS%Set_Interp( trim(line), status )
          IF_NOTOK_RETURN(status=1)
        else
          write (gol,'("no pointer associated, could net set temporal interpolation")'); call goErr
          TRACEBACK; status=1; return
        end if
      end if ! with tinterp
    end if
    
    ! read routine call, default empty:
    !    MyRoutine( invar1, ... [; outvar1, ..] )
    call ReadRc( rcF, trim(rcbase)//'.call', line, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! defined ?
    if ( len_trim(line) > 0 ) then
      ! check on presence of '(..)'
      k = len_trim(line)
      if ( (index(line(1:k),'(') < 1) .or. (line(k:k) /= ')') ) then
        write (gol,'("call description should have form: function(var,..)")'); call goErr
        write (gol,'("but did not find `(..)` in: ",a)') trim(line); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract routine name:
      call goReadFromLine( line, self%routine, status, sep='(' )
      IF_NOTOK_RETURN(status=1)
      ! remove closing parenthese:
      call goTranslate( line, ')', ' ', status )
      IF_NOTOK_RETURN(status=1)
      ! extract input variables:
      call goReadFromLine( line, names, status, sep=';' )
      IF_NOTOK_RETURN(status=1)
      ! split into input variable names:
      call goSplitString( names, self%ndep, self%deps, status, sep=',' )
      IF_NOTOK_RETURN(status=1)
      ! remainder ?
      if ( len_trim(line) > 0 ) then
        ! split into output variable names:
        call goSplitString( line, self%nout, self%outs, status, sep=',' )
        IF_NOTOK_RETURN(status=1)
      else
        ! no output names:
        self%nout = 0
      end if
    else
      ! none:
      self%routine = ''
      self%ndep = 0
    end if

    ! ok
    status = 0
    
  end subroutine Variable_Init


  ! ***


  subroutine Variable_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(inout)              ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Done'
      
    ! --- begin ----------------------------------

    ! switch:
    select case ( trim(self%datatype) )
      !
      case ( 'field' )
        ! done:
        call self%F%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%F, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'constant_field' )
        ! done:
        call self%cF%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%cF, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'constant_field_series' )
        ! done:
        call self%cFS%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%cFS, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'accumulated_field_series' )
        ! done:
        call self%aFS%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%aFS, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'instant_field' )
        ! done:
        call self%iF%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%iF, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'instant_field_series' )
        ! done:
        call self%iFS%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%iFS, stat=status )
      !
      case ( 'vectorfield' )
        ! done:
        call self%VF%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%VF, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'instant_vectorfield' )
        ! done:
        call self%iVF%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%iVF, stat=status )
        IF_NOTOK_RETURN(status=1)
      !
      case ( 'instant_vectorfield_series' )
        ! done:
        call self%iVFS%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( self%iVFS, stat=status )
      !
      case default
        write (gol,'("unsupported datatype `",a,"`")') trim(self%datatype); call goErr
        TRACEBACK; status=1; return
      !
    end select
    
    ! ok
    status = 0
    
  end subroutine Variable_Done


  ! ***


  !
  ! Get elements of variable.
  !
  
  subroutine Variable_Get( self, status, &
                                FieldBase, &
                                Field, Constant_Field, Accumulated_Field_Series, &
                                Instant_Field, Instant_Field_Series, &
                                VectorField, Instant_VectorField, &
                                enabled, levtype, gridtype, lbo, ubo, units, &
                                check_units, check_lbo, check_ubo, &
                                pdata, pudata, pvdata )

    use GO            , only : T_FieldBase
    use LE_Data_Common, only : LE_Data_CompareUnits
  
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(in)                         ::  self
    integer, intent(out)                                  ::  status
    class(T_FieldBase), pointer, optional                 ::  FieldBase
    class(T_Field), pointer, optional                     ::  Field
    class(T_Instant_Field), pointer, optional             ::  Instant_Field
    class(T_Instant_Field_Series), pointer, optional      ::  Instant_Field_Series
    class(T_Constant_Field), pointer, optional            ::  Constant_Field
    class(T_Accumulated_Field_Series), pointer, optional  ::  Accumulated_Field_Series
    class(T_VectorField), pointer, optional               ::  VectorField
    class(T_Instant_VectorField), pointer, optional       ::  Instant_VectorField
    logical, intent(out), optional                        ::  enabled
    character(len=*), intent(out), optional               ::  levtype
    character(len=*), intent(out), optional               ::  gridtype
    integer, intent(out), optional                        ::  lbo(3), ubo(3)
    character(len=*), intent(out), optional               ::  units
    character(len=*), intent(in), optional                ::  check_units
    integer, intent(in), optional                         ::  check_lbo(3), check_ubo(3)
    real, pointer, optional                               ::  pdata(:,:,:)
    real, pointer, optional                               ::  pudata(:,:,:)
    real, pointer, optional                               ::  pvdata(:,:,:)
      
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Get'
      
    ! --- local ----------------------------------
    
    class(T_FieldBase), pointer         ::  FB
    class(T_Field), pointer             ::  F
    class(T_VectorField), pointer       ::  VF
    integer                             ::  F_lbo(3), F_ubo(3)
    
    ! --- begin ----------------------------------
    
    ! init as undefined:
    nullify( FB )
    nullify( F )
    nullify( VF )
    ! switch:
    select case ( trim(self%datatype) )
      case ( 'field'                      ) ; FB => self%F    ; F  => self%F
      case ( 'constant_field'             ) ; FB => self%cF   ; F  => self%cF
      case ( 'constant_field_series'      ) ; FB => self%cFS  ; F  => self%cFS
      case ( 'accumulated_field_series'   ) ; FB => self%aFS  ; F  => self%aFS
      case ( 'instant_field'              ) ; FB => self%iF   ; F  => self%iF
      case ( 'instant_field_series'       ) ; FB => self%iFS  ; F  => self%iFS
      case ( 'vectorfield'                ) ; FB => self%VF   ; VF => self%VF
      case ( 'instant_vectorfield'        ) ; FB => self%iVF  ; VF => self%iVF
      case ( 'instant_vectorfield_series' ) ; FB => self%iVFS ; VF => self%iVFS
      case default
        write (gol,'("could not associate FieldBase pointer for variable `",a,"`")') trim(self%name); call goErr
        TRACEBACK; status=1; return
    end select

    ! return pointer to fieldbase ?
    if ( present(FieldBase) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'field'                      ) ; FieldBase => self%F
        case ( 'constant_field'             ) ; FieldBase => self%cF
        case ( 'constant_field_series'      ) ; FieldBase => self%cFS
        case ( 'accumulated_field_series'   ) ; FieldBase => self%aFS
        case ( 'instant_field'              ) ; FieldBase => self%iF
        case ( 'instant_field_series'       ) ; FieldBase => self%iFS
        case ( 'vectorfield'                ) ; FieldBase => self%VF
        case ( 'instant_vectorfield'        ) ; FieldBase => self%iVF
        case ( 'instant_vectorfield_series' ) ; FieldBase => self%iVFS
        case default
          write (gol,'("could not associate FieldBase pointer for variable `",a,"`")') trim(self%name); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to field?
    if ( present(Field) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'field'                    ) ; Field => self%F
        case ( 'constant_field'           ) ; Field => self%cF
        case ( 'constant_field_series'    ) ; Field => self%cFS
        case ( 'accumulated_field_series' ) ; Field => self%aFS
        case ( 'instant_field'            ) ; Field => self%iF
        case ( 'instant_field_series'     ) ; Field => self%iFS
        case default
          write (gol,'("could not associate Field pointer for variable `",a,"`")') trim(self%name); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to field?
    if ( present(Instant_Field) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'instant_field'         ) ; Instant_Field => self%iF
        case ( 'instant_field_series'  ) ; Instant_Field => self%iFS
        case default
          write (gol,'("could not associate Instant_Field pointer for variable `",a,"`")') trim(self%name); call goErr
          write (gol,'("defined datatype `",a,"`")') trim(self%datatype); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to field?
    if ( present(Instant_Field_Series) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'instant_field_series'  ) ; Instant_Field_Series => self%iFS
        case default
          write (gol,'("could not associate Instant_Field_Series pointer for variable `",a,"`")') trim(self%name); call goErr
          write (gol,'("defined datatype `",a,"`")') trim(self%datatype); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to field?
    if ( present(Constant_Field) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'constant_field'           ) ; Constant_Field => self%cF
        case ( 'constant_field_series'    ) ; Constant_Field => self%cFS
        case ( 'accumulated_field_series' ) ; Constant_Field => self%aFS
        case default
          write (gol,'("could not associate Constant_Field pointer for variable `",a,"`")') trim(self%name); call goErr
          write (gol,'("defined datatype `",a,"`")') trim(self%datatype); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to field?
    if ( present(Accumulated_Field_Series) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'accumulated_field_series'  ) ; Accumulated_Field_Series => self%aFS
        case default
          write (gol,'("could not associate Accumulated_Field_Series pointer for variable `",a,"`")') trim(self%name); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to vector field?
    if ( present(VectorField) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'vectorfield'                 ) ; VectorField => self%VF
        case ( 'instant_vectorfield'         ) ; VectorField => self%iVF
        case ( 'instant_vectorfield_series'  ) ; VectorField => self%iVFS
        case default
          write (gol,'("could not associate VectorField pointer for variable `",a,"`")') trim(self%name); call goErr
          write (gol,'("defined datatype `",a,"`")') trim(self%datatype); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return pointer to vector field?
    if ( present(Instant_VectorField) ) then
      ! switch:
      select case ( trim(self%datatype) )
        case ( 'instant_vectorfield'         ) ; Instant_VectorField => self%iVF
        case ( 'instant_vectorfield_series'  ) ; Instant_VectorField => self%iVFS
        case default
          write (gol,'("could not associate Instant_VectorField pointer for variable `",a,"`")') trim(self%name); call goErr
          write (gol,'("defined datatype `",a,"`")') trim(self%datatype); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! Field

    ! return element of field?
    if ( any((/present(enabled),present(lbo),present(ubo)/)) ) then
      ! get flag:
      call FB%FieldBase_Get( status, enabled=enabled, lbo=lbo, ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    end if ! present enabled
    
    ! attributes:
    if ( present(units   ) ) units    = trim(self%units)
    if ( present(levtype ) ) levtype  = trim(self%levtype)
    if ( present(gridtype) ) gridtype = trim(self%gridtype)
    
    ! check units?
    if ( present(check_units) ) then
      ! compare, status<0 if not equivalent:
      call LE_Data_CompareUnits( trim(self%units), trim(check_units), status )
      IF_ERROR_RETURN(status=1)
      if ( status < 0 ) then
        write (gol,'("variable `",a,"` has units `",a,"` while expected `",a,"`")') &
                        trim(self%name), trim(self%units), trim(check_units); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! check index space?
    if ( present(check_lbo) ) then
      ! get variable bounds:
      call FB%FieldBase_Get( status, lbo=F_lbo, ubo=F_ubo )
      IF_NOTOK_RETURN(status=1)
      ! index space should include requested lower bounds:
      if ( any( check_lbo < F_lbo ) ) then
        write (gol,'("requested lower bounds (",3i6,") below index space of `",a,"` :")') &
                 check_lbo, trim(self%name); call goErr
        write (gol,'("  lower bounds : (",3i6,")")') F_lbo; call goErr
        write (gol,'("  upper bounds : (",3i6,")")') F_ubo; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check index space?
    if ( present(check_ubo) ) then
      ! get variable bounds:
      call FB%FieldBase_Get( status, lbo=F_lbo, ubo=F_ubo )
      IF_NOTOK_RETURN(status=1)
      ! index space should include requested lower bounds:
      if ( any( check_ubo > F_ubo ) ) then
        write (gol,'("requested upper bounds (",3i6,") above index space of `",a,"` :")') &
                 check_ubo, trim(self%name); call goErr
        write (gol,'("  lower bounds : (",3i6,")")') F_lbo; call goErr
        write (gol,'("  upper bounds : (",3i6,")")') F_ubo; call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! return field data?
    if ( present(pdata) ) then
      ! check ...
      if ( .not. associated(F) ) then
        write (gol,'("could not set pointer pdata since F is not associated")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! get flag:
      call F%Get( status, pdata=pdata )
      IF_NOTOK_RETURN(status=1)
    end if ! present enabled

    ! return vector data?
    if ( any( (/present(pudata),present(pvdata)/) ) ) then
      ! check ...
      if ( .not. associated(VF) ) then
        write (gol,'("could not set pudata/pvdata since VF is not associated")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! get flag:
      call VF%Get( status, pudata=pudata, pvdata=pvdata )
      IF_NOTOK_RETURN(status=1)
    end if ! present enabled
    
    ! ok
    status = 0
    
  end subroutine Variable_Get


  ! ***


  subroutine Variable_Set( self, status, enabled )
  
    use GO, only : T_FieldBase
  
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(inout)              ::  self
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Set'
      
    ! --- local ----------------------------------
    
    class(T_FieldBase), pointer             ::  FB
    
    ! --- begin ----------------------------------

    ! return enabled flag?
    if ( present(enabled) ) then
      ! set flag:
      self%enabled = enabled
      ! pointer to FieldBase:
      call self%Get( status, FieldBase=FB )
      IF_NOTOK_RETURN(status=1)
      ! set flag:
      call FB%Set( status, enabled=enabled )
      IF_NOTOK_RETURN(status=1)
    end if ! present enabled

    ! ok
    status = 0
    
  end subroutine Variable_Set


  ! ***


  !
  ! Check if variable is valid for requested interval;
  ! return status:
  !   -1    : field not valid 
  !    0    : field valid
  !    1    : error (not enabled?)
  !
  
  subroutine Variable_Valid( self, t1, t2, status, silent )
  
    use GO, only : TDate, operator(+), operator(-), operator(*)

    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(inout)              ::  self
    type(TDate), intent(in)                       ::  t1, t2
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Valid'
      
    ! --- local ----------------------------------
    
    class(T_FieldBase), pointer               ::  FB
    class(T_Constant_Field), pointer          ::  cF
    class(T_Instant_Field), pointer           ::  iF
    class(T_Instant_VectorField), pointer     ::  iVF
    logical                                   ::  enabled
    logical                                   ::  filled
    type(TDate)                               ::  tt(2)
    type(TDate)                               ::  tmid
    
    ! --- begin ----------------------------------

    ! switch:
    select case ( trim(self%datatype) )

      !~ field without time info:
      case ( 'field', 'vectorfield' )
        ! get pointer:
        call self%Get( status, FieldBase=FB )
        IF_NOTOK_RETURN(status=1)
        ! should be enabled:
        call FB%Check( status, enabled=.true. )
        IF_NOTOK_RETURN(status=1)
        ! filled ? do not shout if not filled;
        ! check on error status only:
        call FB%Check( status, filled=.true., silent=silent )
        IF_ERROR_RETURN(status=1)

      !~ field constant over interval
      case ( 'constant_field', 'constant_field_series', 'accumulated_field_series' )
        ! get pointer:
        call self%Get( status, Constant_Field=cF )
        IF_NOTOK_RETURN(status=1)
        ! should be enabled:
        call cF%CheckSample( status, enabled=silent )
        IF_NOTOK_RETURN(status=1)
        ! fill interval:
        tt(1) = t1
        tt(2) = t2
        ! valid for this timestep? do not shout errors,
        ! check on errors only:
        call cF%CheckSample( status, tt=tt, silent=silent )
        IF_ERROR_RETURN(status=1)
      
      !~ field valid for instant time:
      case ( 'instant_field', 'instant_field_series'  )
        ! get pointer:
        call self%Get( status, Instant_Field=iF )
        IF_NOTOK_RETURN(status=1)
        ! should be enabled:
        call iF%CheckSample( status, enabled=.true. )
        IF_NOTOK_RETURN(status=1)
        ! mid time:
        tmid = t1 + (t2-t1)*0.5
        ! valid for this time? do not shout errors,
        ! check on errors only:
        call iF%CheckSample( status, t=tmid, silent=silent )
        IF_ERROR_RETURN(status=1)
      
      !~ field valid for instant time:
      case ( 'instant_vectorfield', 'instant_vectorfield_series'  )
        ! get pointer:
        call self%Get( status, Instant_VectorField=iVF )
        IF_NOTOK_RETURN(status=1)
        ! should be enabled:
        call iVF%CheckSample( status, enabled=.true. )
        IF_NOTOK_RETURN(status=1)
        ! mid time:
        tmid = t1 + (t2-t1)*0.5
        ! valid for this time? do not shout errors,
        ! check on errors only:
        call iVF%CheckSample( status, t=tmid, silent=silent )
        IF_ERROR_RETURN(status=1)

      case default
        write (gol,'("unsupported data type `",a,"` for variable `",a,"`")') &
                trim(self%datatype), trim(self%name); call goErr
        TRACEBACK; status=1; return
    end select

    ! do not set status to zero, keep latest return status ...
    
  end subroutine Variable_Valid


  ! ***


  subroutine Variable_Setup( self, tref, t1, t2, status )
  
    use GO  , only : TDate, operator(+), operator(-), operator(*)
  
    use LE_Data_Meteo, only : LE_Data_Setup_cFS_sfc
    use LE_Data_Meteo, only : LE_Data_Setup_cFS_3D
    use LE_Data_Meteo, only : LE_Data_Setup_aFS_sfc
    use LE_Data_Meteo, only : LE_Data_Setup_iFS_sfc
    use LE_Data_Meteo, only : LE_Data_Setup_iFS_3D
    use LE_Data_Meteo, only : LE_Data_Setup_iVFS_sfc
    use LE_Data_Meteo, only : LE_Data_Setup_iVFS_3D
    
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(inout)              ::  self
    type(TDate), intent(in)                       ::  tref
    type(TDate), intent(in)                       ::  t1, t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Setup'
      
    ! --- local ----------------------------------
    
    logical                ::  halflevels
    logical                ::  enabled
    type(TDate)            ::  tmid
      
    ! --- begin ----------------------------------
    
    ! mid time:
    tmid = t1 + (t2-t1)*0.5

    ! switch:
    select case ( trim(self%datatype) )

      !case ( 'field' )
      
      !case ( 'constant_field' )
      
      !~
      case ( 'constant_field_series' )
      
        ! switch:
        select case ( trim(self%levtype) )
          !~ 2D fields:
          case ( 'sfc' )
            ! setup, just mid time to find record:
            call LE_Data_Setup_cFS_sfc( self%cFS, tref, (/t1,t2/), self%rcF, &
                                          self%gridtype, self%gridmapping, status )
            IF_NOTOK_RETURN(status=1)
          !~ 3D fields at corners:
          case ( 'input_levels', 'input_halflevels' )
            ! set flag:
            halflevels = trim(self%levtype) == 'input_halflevels'
            ! setup, just mid time to find record:
            call LE_Data_Setup_cFS_3D( self%cFS, tref, (/t1,t2/), self%rcF, &
                                         self%gridtype, self%gridmapping, status, &
                                         halflevels=halflevels )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("unsupported levtype `",a,"`")') trim(self%levtype); call goErr
            TRACEBACK; status=1; return
        end select
      
      !~
      case ( 'accumulated_field_series' )
      
        ! switch:
        select case ( trim(self%levtype) )
          !~ 2D fields:
          case ( 'sfc' )
            ! setup, just mid time to find record:
            call LE_Data_Setup_aFS_sfc( self%aFS, tref, (/t1,t2/), self%rcF, &
                                          self%gridtype, self%gridmapping, status )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("unsupported levtype `",a,"`")') trim(self%levtype); call goErr
            TRACEBACK; status=1; return
        end select

      !~
      case ( 'instant_field_series' )
      
        ! switch:
        select case ( trim(self%levtype) )
          !~ 2D fields:
          case ( 'sfc' )
            ! setup, just mid time to find record:
            call LE_Data_Setup_iFS_sfc( self%iFS, tref, tmid, self%rcF, &
                                          self%gridtype, self%gridmapping, status )
            IF_NOTOK_RETURN(status=1)
          !~ 3D fields at corners:
          case ( 'input_levels', 'input_halflevels' )
            ! set flag:
            halflevels = trim(self%levtype) == 'input_halflevels'
            ! setup, just mid time to find record:
            call LE_Data_Setup_iFS_3D( self%iFS, tref, tmid, self%rcF, &
                                         self%gridtype, self%gridmapping, status, &
                                         halflevels=halflevels )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("unsupported levtype `",a,"`")') trim(self%levtype); call goErr
            TRACEBACK; status=1; return
        end select

      !~
      case ( 'instant_vectorfield_series' )
      
        ! switch:
        select case ( trim(self%levtype) )
          !~ 2D fields:
          case ( 'sfc' )
            ! setup, just mid time to find record:
            call LE_Data_Setup_iVFS_sfc( self%iVFS, tref, tmid, self%rcF, &
                                          self%gridtype, self%gridmapping, status )
            IF_NOTOK_RETURN(status=1)
          !~ 3D fields at corners:
          case ( 'input_levels', 'input_halflevels' )
            ! set flag:
            halflevels = trim(self%levtype) == 'input_halflevels'
            ! setup, just mid time to find record:
            call LE_Data_Setup_iVFS_3D( self%iVFS, tref, tmid, self%rcF, &
                                         self%gridtype, self%gridmapping, status, &
                                         halflevels=halflevels )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("unsupported levtype `",a,"`")') trim(self%levtype); call goErr
            TRACEBACK; status=1; return
        end select

      case default
        write (gol,'("unsupported datatype `",a,"`")') trim(self%datatype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine Variable_Setup


  ! ***


  subroutine Variable_Truncate( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Variable), intent(inout)              ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variable_Truncate'
      
    ! --- local ----------------------------------
    
    real, pointer     ::  pdata(:,:,:)
    real, pointer     ::  pudata(:,:,:), pvdata(:,:,:)
      
    ! --- begin ----------------------------------

    ! switch:
    select case ( trim(self%datatype) )

      !~
      case ( 'field', 'constant_field', 'constant_field_series', 'accumulated_field_series', &
                      'instant_field', 'instant_field_series')
    
        ! pointer to data:
        call self%Get( status, pdata=pdata )
        IF_NOTOK_RETURN(status=1)

        ! truncate below threshold?
        if ( self%with_minimum ) then
          ! reset data:
          pdata = max( self%minimum, pdata )
        end if

        ! truncate above threshold?
        if ( self%with_maximum ) then
          ! reset data:
          pdata = min( pdata, self%maximum )
        end if

      !~
      case ( 'vectorfield', 'instant_vectorfield', 'instant_vectorfield_series')
    
        ! pointer to data:
        call self%Get( status, pudata=pudata, pvdata=pvdata )
        IF_NOTOK_RETURN(status=1)

        ! truncate below threshold?
        if ( self%with_minimum ) then
          ! reset data:
          pudata = max( self%minimum, pudata )
          pvdata = max( self%minimum, pvdata )
        end if

        ! truncate above threshold?
        if ( self%with_maximum ) then
          ! reset data:
          pudata = min( pudata, self%maximum )
          pvdata = min( pvdata, self%maximum )
        end if

      case default
        write (gol,'("unsupported datatype `",a,"`")') trim(self%datatype); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine Variable_Truncate
    

end module LE_Data_Variable
