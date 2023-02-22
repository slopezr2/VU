!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status>0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_File_Nc

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  maxdim
  public  ::  show_filename_on_open
  public  ::  T_File_Nc
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Nc'

  ! max number of dimensions:
  integer, parameter  ::  maxdim = 7
  
  ! enable to write message when file is opened:
  logical             ::  show_filename_on_open = .false.


  ! --- types ----------------------------------------
  
  type T_File_Nc
    character(len=1024)     ::  filename
    integer                 ::  ncid
    integer                 ::  formatNum
  contains
    !
    procedure   ::  Open                 => File_NC_Open
    procedure   ::  Create               => File_NC_Create
    procedure   ::  EndDef               => File_NC_EndDef
    procedure   ::  Close                => File_NC_Close
    procedure   ::  Inquire              => File_NC_Inquire
    !
    procedure   ::                          File_NC_Inquire_Dimension_ind
    procedure   ::                          File_NC_Inquire_Dimension_name
    generic     ::  Inquire_Dimension    => File_NC_Inquire_Dimension_ind, &
                                            File_NC_Inquire_Dimension_name
    !
    procedure   ::  Copy_Dimension       => File_NC_Copy_Dimension
    !
    procedure   ::  Inq_VarID            => File_NC_Inq_VarID
    procedure   ::  Inq_Var              => File_NC_Inq_Var
    procedure   ::  Inq_VarPacking       => File_NC_Inq_VarPacking
    procedure   ::  Inq_VarUnits         => File_NC_Inq_VarUnits
    procedure   ::  Inq_VarMissing       => File_NC_Inq_VarMissing
    procedure   ::  Inq_VarCoordinates   => File_NC_Inq_VarCoordinates
    procedure   ::  Define_Variable_Copy => File_NC_Define_Variable_Copy
    procedure   ::  Copy_Atts            => File_NC_Copy_Atts
    procedure   ::  Extend_History       => File_NC_Extend_History
    !
    procedure   ::                          File_NC_Get_Var_Attr_s
    procedure   ::                          File_NC_Get_Var_Attr_t
    generic     ::  Get_Var_Attr         => File_NC_Get_Var_Attr_s, &
                                            File_NC_Get_Var_Attr_t
    !
    procedure   ::                          File_NC_Get_Var_i_1d
    procedure   ::                          File_NC_Get_Var_i_2d
    procedure   ::                          File_NC_Get_Var_i_3d
    procedure   ::                          File_NC_Get_Var_r_1d
    procedure   ::                          File_NC_Get_Var_r_2d
    procedure   ::                          File_NC_Get_Var_r_3d
    procedure   ::                          File_NC_Get_Var_t_2d
    procedure   ::                          File_NC_Get_Var_t_1d
    generic     ::  Get_Var              => File_NC_Get_Var_i_1d, &
                                            File_NC_Get_Var_i_2d, &
                                            File_NC_Get_Var_i_3d, &
                                            File_NC_Get_Var_r_1d, &
                                            File_NC_Get_Var_r_2d, &
                                            File_NC_Get_Var_r_3d, &
                                            File_NC_Get_Var_t_1d, &
                                            File_NC_Get_Var_t_2d
    !
    procedure   ::                          File_NC_Inq_TimeRecord
    procedure   ::                          File_NC_Inq_TimeBoundsRecord
    generic     ::  Inq_TimeRecord       => File_NC_Inq_TimeRecord, &
                                            File_NC_Inq_TimeBoundsRecord
    !
    procedure   ::                          File_Nc_Get_Levs_bnds
    procedure   ::                          File_Nc_Get_Levs_h
    generic     ::  Get_Levs             => File_Nc_Get_Levs_bnds, &
                                            File_Nc_Get_Levs_h
  end type T_File_Nc

  

contains


  ! ********************************************************************
  ! ***
  ! *** ncfile
  ! ***
  ! ********************************************************************


  subroutine File_NC_Open( self, filename, status )
  
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_NOWRITE
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(out)   ::  self
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Open'
    
    ! --- local ----------------------------------
    
    logical      ::  exist
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! check ..
    inquire( file=trim(self%filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("WARNING - file not found : ",a)') trim(self%filename); call goPr
      status=-1; return
    end if
    
    ! Export filename to log?
    if ( show_filename_on_open ) then
      write (gol,'("Attempting to read from file: ", a)') trim(filename); call goPr
    end if
    
    ! open file for reading:
    status = NF90_Open( trim(filename), NF90_NOWRITE, self%ncid )
    if ( status /= NF90_NOERR ) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("opening filename: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    status = NF90_Inquire( self%ncid, formatNum=self%formatNum )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine File_NC_Open


  ! ***


  subroutine File_NC_Create( self, filename, status )
  
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_CLOBBER
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(out)    ::  self
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Create'
    
    ! --- local ----------------------------------
    
    integer     ::  cmode
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! creation mode:
    cmode = 0
    cmode = cmode + NF90_CLOBBER     ! overwrite existing
    cmode = cmode + self%formatNum   ! NetCDF flavour
    
    ! open file for reading:
    status = NF90_Create( trim(filename), cmode, self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine File_NC_Create


  ! ***
  
  
  subroutine File_NC_EndDef( self, status )

    use NetCDF, only : NF90_EndDef
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(inout)  ::  self
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_EndDef( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_NC_EndDef


  ! ***
  
  
  subroutine File_NC_Close( self, status )

    use NetCDF, only : NF90_Close
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(inout)  ::  self
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Close'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_Close( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_NC_Close


  ! ***
  
  
  subroutine File_NC_Inquire( self, status, ndim, nvar )

    use NetCDF, only : NF90_Inquire
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)     ::  self
    integer, intent(out)            ::  status
    integer, intent(out), optional  ::  ndim
    integer, intent(out), optional  ::  nvar

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inquire'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! obtain info:
    status = NF90_Inquire( self%ncid, nDimensions=ndim, nVariables=nvar )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inquire


  ! ***
  
  
  subroutine File_NC_Inquire_Dimension_ind( self, idim, status, name, length )

    use NetCDF, only : NF90_Inquire_Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)              ::  self
    integer, intent(in)                       ::  idim
    integer, intent(out)                      ::  status
    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  length

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inquire_Dimension_ind'
    
    ! --- local ----------------------------------
    
    integer     ::  dimid
    
    ! --- begin ----------------------------------
    
    ! variable number is variable id ?
    dimid = idim
    ! obtain info:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=name, len=length )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inquire_Dimension_ind


  ! ***
  
  
  subroutine File_NC_Inquire_Dimension_name( self, name, status, length )

    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)              ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status
    integer, intent(out), optional            ::  length

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inquire_Dimension_name'
    
    ! --- local ----------------------------------
    
    integer     ::  dimid
    
    ! --- begin ----------------------------------

    ! obtain dimension id for requested name:
    status = NF90_Inq_DimID( self%ncid, name, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain info:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=length )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inquire_Dimension_name


  ! ***
  
  
  subroutine File_NC_Copy_Dimension( self, dimname, outfile, status )

    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Def_Dim
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)              ::  self
    character(len=*), intent(in)              ::  dimname
    class(T_File_Nc), intent(inout)           ::  outfile
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Copy_Dimension'
    
    ! --- local ----------------------------------
    
    integer     ::  dimid
    integer     ::  length
    
    ! --- begin ----------------------------------
    
    ! obtain dimension id for requested name:
    status = NF90_Inq_DimID( self%ncid, dimname, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain info:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=length )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! create same dimension in output file,
    ! returned dimension id is not used:
    status = NF90_Def_Dim( outfile%ncid, dimname, length, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Copy_Dimension


  ! ***

  
  !
  ! Return id of variable identified by it's name or an
  ! attribute value. The description is a ';' seperated list 
  ! of 'key=value' pairs, where 'var_name' should be used
  ! for the variable name and all other names are assumed
  ! to be attributes:
  !
  !    standard_name=pressure;long_name=air pressure;var_name=p
  !
  ! Duplicates are allowed, for example:
  !
  !    long_name=air pressure;long_name=Pressure
  !
  ! If a match is found the variable id is returned.
  !
  ! Return status:
  !   < 0  : no variable found, and quiet=.true.
  !     0  : variable found, varid filled
  !   > 0  : variable not found and quiet=.false. (default), or other error
  !
  
  subroutine File_NC_Inq_VarID( self, description, varid, status, quiet )

    use NetCDF, only : NF90_INQ_VarID
    use NetCDF, only : NF90_Get_Att
    
    use GO, only : goReadFromLine, goNtsTrim
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  varid
    integer, intent(out)              ::  status
    logical, intent(in), optional     ::  quiet

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_VarID'
    
    ! --- local ----------------------------------
    
    logical             ::  verb
    logical             ::  found
    integer             ::  nvar, ivar
    character(len=1024) ::  line
    character(len=64)   ::  name
    character(len=256)  ::  value
    character(len=256)  ::  attr_value
    
    ! --- begin ----------------------------------
    
    ! verbose?
    verb = .true.
    if ( present(quiet) ) verb = .not. quiet
    
    ! set flag:
    found = .false.
    
    ! copy:
    line = description
    ! loop over parts:
    do
      ! empty?
      if ( len_trim(line) == 0 ) exit

      ! split part:
      call goReadFromLine( line, value, status, sep=';' )
      IF_NOT_OK_RETURN(status=1)
      ! split-off first part:
      call goReadFromLine( value, name, status, sep='=' )
      IF_NOT_OK_RETURN(status=1)

      ! switch:
      select case ( name )
        !~
        case ( 'var_name' )
          ! if variable is present with name equal to the value,
          ! then the variable id is returned and status is not error;
          ! otherwise (no variable with this name) an error status is returned:
          status = NF90_INQ_VarID( self%ncid, trim(value), varid )
          found = status == NF90_NOERR
        !~
        case default
          ! number of variables:
          call self%Inquire( status, nvar=nvar )
          IF_NOT_OK_RETURN(status=1)
          ! loop:
          do ivar = 1, nvar
            ! variable id is number:
            varid = ivar
            ! get name if present, error status if not exists:
            status = NF90_Get_Att( self%ncid, varid, trim(name), attr_value )
            if ( status /= NF90_NOERR ) cycle
            ! fix nul-terminated strings:
            attr_value = goNtsTrim( attr_value )
            ! match?
            found = trim(attr_value) == trim(value)
            ! leave ?
            if ( found ) exit
          end do
        !~
      end select
      
      ! leave?
      if ( found ) exit
      
    end do  ! key=value pairs
      
    ! check ...
    if ( .not. found ) then
      if ( verb ) then
        write (gol,'("no variable found matching description: ",a)') trim(description); call goErr
        write (gol,'("  file: ",a)') trim(self%filename); call goErr
        TRACEBACK; status=1; return
      else
        status = -1; return
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_VarID
  
  
  ! ***

  
  ! 
  ! Inquire variable (by some name):
  !   variable id is numbered 1,..,nvar
  !   time_dim_name   : name of time dimension, or 'None' if not present
  !   level_dim_name  : name of level dimension, or 'None' if not present
  !
  
  subroutine File_NC_Inq_Var( self, varid, status, &
                                name, ndim, time_dim_name, &
                                level_dim_name, &
                                dim, dim_length )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    character(len=*), intent(out), optional     ::  name
    integer, intent(out), optional              ::  ndim
    integer, intent(in), optional               ::  dim
    integer, intent(out), optional              ::  dim_length
    character(len=*), intent(out), optional     ::  time_dim_name
    character(len=*), intent(out), optional     ::  level_dim_name

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_Var'
    
    ! --- local ----------------------------------
    
    integer                 ::  ndims
    integer, allocatable    ::  dimids(:)
    integer                 ::  idim
    character(len=32)       ::  dim_name
    integer                 ::  n

    ! --- begin ----------------------------------
    
    ! return name?
    if ( present(name) ) then
      ! obtain info:
      status = NF90_Inquire_Variable( self%ncid, varid, name=name )
      IF_NF90_NOT_OK_RETURN(status=1)   
    end if
    
    ! return number of dimensions?
    if ( present(ndim) ) then
      ! get  number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndim )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! return dimension length?
    if ( present(dim_length) ) then
      ! check ...
      if ( .not. present(dim) ) then
        write (gol,'("dim_length requested without dim (1..ndim) provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! get number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndims )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! storage:
      allocate( dimids(ndims), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get length of dimension:
      status = NF90_Inquire_Dimension( self%ncid, dimids(dim), len=dim_length )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! dim_name
    
   ! return name of time or level dimension:
    if ( present(time_dim_name) .or. present(level_dim_name) ) then
      ! default:
      if ( present(time_dim_name ) ) time_dim_name  = 'None'
      if ( present(level_dim_name) ) level_dim_name = 'None'
      ! get  number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndims )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! storage:
      allocate( dimids(ndims), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! loop:
      do idim = 1, ndims
        ! get name of dimension:
        status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dim_name )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! assign to requested dimension names for level or time if possible;
        ! in future, also check on 'standard_name' attribute of coordinate (if present):
        select case ( trim(dim_name) )
          !~ spatial:
          case ( 'longitude', 'latitude', 'lon', 'lat', 'slon', 'slat', &
                 'rgrid', 'rlon', 'rlat', 'srlon', 'srlat', &
                 'x', 'y', &
                 'west_east', 'west_east_stag', &
                 'south_north', 'south_north_stag' )
            ! no time or level ...
          !~ level:
          case ( 'level', 'level1', 'lev', 'depth', 'soil', 'soil1', &
                 'bottom_top', 'bottom_top_stag', 'soil_layers_stag', 'lev_2', 'lev_3', &
                 'hlevel', 'slevel', 'height', 'height_5' )
            ! found level dim:
            if ( present(level_dim_name) ) level_dim_name = trim(dim_name)
          !~ known time dimensions:
          case ( 'time', 'Time', 'XTIME' )
            ! found time dim:
            if ( present(time_dim_name ) ) time_dim_name  = trim(dim_name)
            ! leave:
            exit
          !~ unknown ..
          case default
            write (gol,'("could not decide if dimension ",i0," named `",a,"`")') idim, trim(dim_name); call goErr
            write (gol,'("is a time, a level, or something else")'); call goErr
            TRACEBACK; status=1; return
        end select
      end do  ! dim names
      ! clear:
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! dim_name
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_Var
  
  ! *
  
  ! ***
  
  
  subroutine File_NC_Define_Variable_Copy( self, varname, outfile, varid, status )

    use NetCDF, only : NF90_Inquire_Dimension, NF90_Inq_DimID
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Def_Var

    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)              ::  self
    character(len=*), intent(in)              ::  varname
    class(T_File_Nc), intent(inout)           ::  outfile
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Define_Variable_Copy'
    
    ! --- local ----------------------------------
    
    integer             ::  varid_in
    integer             ::  xtype
    integer             ::  ndim, idim
    integer             ::  dimids_in(maxdim)
    integer             ::  dimids(maxdim)
    character(len=64)   ::  dimname
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    status = NF90_INQ_VarID( self%ncid, varname, varid_in )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain info:
    status = NF90_Inquire_Variable( self%ncid, varid_in, xtype=xtype, &
                                           ndims=ndim, dimids=dimids_in )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! get dimension id's in output file:
    do idim = 1, ndim
      ! get name:
      status = NF90_Inquire_Dimension( self%ncid, dimids_in(idim), name=dimname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get id in output file:
      status = NF90_Inq_DimID( outfile%ncid, dimname, dimids(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! define new variable:
    status = NF90_Def_Var( outfile%ncid, varname, xtype, dimids(1:ndim), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! copy attributes:
    call self%Copy_Atts( varname, outfile, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Define_Variable_Copy


  ! ***
  
  
  !
  ! Copy attributes from one file to another ;
  ! 'varname' should be the variable name,
  ! or empty or '/' for global attributes.
  !
  
  subroutine File_NC_Copy_Atts( self, varname, outfile, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName, NF90_Copy_Att
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)     ::  self
    character(len=*), intent(in)    ::  varname
    class(T_File_Nc), intent(inout)  ::  outfile
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Copy_Atts'
    
    ! --- local ----------------------------------
    
    integer             ::  varid, varid_out
    integer             ::  natt, iatt
    character(len=64)   ::  attname
    
    ! --- begin ----------------------------------
    
    ! get variable id's:
    if ( (len_trim(varname) == 0) .or. (trim(varname) == '/') ) then
      ! set variable id's:
      varid = NF90_GLOBAL
      varid_out = NF90_GLOBAL
      ! count number of attributes:
      status = NF90_Inquire( self%ncid, nAttributes=natt )
      IF_NF90_NOT_OK_RETURN(status=1)
    else
      ! get variable id's:
      status = NF90_INQ_VarID( self%ncid, varname, varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_INQ_VarID( outfile%ncid, varname, varid_out )
      IF_NF90_NOT_OK_RETURN(status=1) 
      ! count number of attributes:
      status = NF90_Inquire_Variable( self%ncid, varid, nAtts=natt )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! loop over attributes:
    do iatt = 1, natt
      ! name:
      status = NF90_Inq_AttName( self%ncid, varid, iatt, attname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! copy:
      status = NF90_Copy_Att( self%ncid, varid, trim(attname), &
                               outfile%ncid, varid_out )
      IF_NF90_NOT_OK_RETURN(status=1) 
    end do ! attributes
    
    ! ok
    status = 0
    
  end subroutine File_NC_Copy_Atts
  
  
  ! ***
  
  
  !
  ! Prepend line to history:
  !   Wed Apr 23 16:32:03 2014: Message\n
  !
  ! Only allowed in definition mode.
  !
  
  subroutine File_NC_Extend_History( self, message, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Get_Att, NF90_Put_Att
    
    use GO, only : TDate, SystemDate
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(inout)  ::  self
    character(len=*), intent(in)    ::  message
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Extend_History'
    
    ! month names:
    character(len=3),parameter :: monthname(12) = &
                         (/'Jan','Feb','Mrc','Apr','May','Jun',&
                           'Jul','Aug','Sep','Oct','Nov','Dec'/)
   
    ! newline character:
    character(len=1), parameter  ::  newline = char(10)
     
    ! --- local ----------------------------------
    
    character(len=64)     ::  attname
    character(len=64)     ::  tstamp
    character(len=1024)   ::  history
    type(TDate)           ::  t
    
    ! --- begin ----------------------------------
    
    ! target attribute:
    attname = 'history'
    
    ! current:
    status = NF90_Get_Att( self%ncid, NF90_GLOBAL, attname, history )
    if ( status /= NF90_NOERR ) history = ''
    
    ! current time:
    t = SystemDate()
    ! write time stamp:
    write (tstamp,'("Day ",a," ",i2," ",i2,2(":",i2.2)," ",i4)') &
             monthname(t%month), t%day, t%hour, t%min, t%sec, t%year
    ! fill history:
    if ( len_trim(history)== 0 ) then
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message)
    else
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message), &
                                               newline, trim(history)
    end if
    
    ! write:
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, attname, history )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Extend_History
  
  
  ! *
  
  ! Return packing parameters 'add_offset' and 'scale_factor'.
  ! Sometimes files are found that have only one of these,
  ! the missing parameter is then set to 0.0 or 1.0 respectively.
  ! Return status -1 if none of the parameters is found.
  !

  subroutine File_NC_Inq_VarPacking( self, varid, add_offset, scale_factor, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  add_offset, scale_factor
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_VarPacking'
    
    ! --- local ----------------------------------

    logical     ::  any_present
    
    ! --- begin ----------------------------------
    
    ! by default nothing found:
    any_present = .false.
    
    ! try to get packing variable, ok if not present:
    status = NF90_Get_Att( self%ncid, varid, 'add_offset', add_offset )
    if ( status == NF90_NOERR ) then
      ! reset flag:
      any_present = .true.
    else if ( status == NF90_ENOTATT ) then
      ! default value:
      add_offset = 0.0
    else
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! try to get packing variable, ok if not present:
    status = NF90_Get_Att( self%ncid, varid, 'scale_factor', scale_factor )
    if ( status == NF90_NOERR ) then
      ! reset flag:
      any_present = .true.
    else if ( status == NF90_ENOTATT ) then
      ! default value:
      scale_factor = 1.0
    else
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! any found?
    if ( any_present ) then
      status = 0
    else
      status = -1
    end if
    
  end subroutine File_NC_Inq_VarPacking
  
  
  ! *
  
  ! Return units, eventually use from description

  subroutine File_NC_Inq_VarUnits( self, varid, description, units, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
    
    use GO, only : goVarValue, goNtsTrim
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    integer, intent(in)               ::  varid
    character(len=*), intent(in)      ::  description
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_VarUnits'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get units from attribute:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    !~ attribute found and read:
    if ( status == NF90_NOERR ) then
      ! remove nul-characters:
      units = goNtsTrim( units )
      ! empty?
      if ( len_trim(units) == 0 ) then
        ! try to read from description, warning status<0 if not defined:
        call goVarValue( description, ';', 'units', '=', units, status )
        IF_ERROR_RETURN(status=1)
        ! not found?
        if ( status < 0 ) then
          write (gol,'("empty units attribute in variable, and no explicit specification in description either:")'); call goErr
          write (gol,'("  filename             : ",a)') trim(self%filename); call goErr
          write (gol,'("  variable description : ",a)') trim(description); call goErr
          TRACEBACK; status=1; return
        end if
      end if
    !~ attribute not found ...
    else if ( status == NF90_ENOTATT ) then
      ! try to read from description, warning status<0 if not defined:
      call goVarValue( description, ';', 'units', '=', units, status )
      IF_ERROR_RETURN(status=1)
      ! not found?
      if ( status < 0 ) then
        write (gol,'("no units attribute in variable, and no explicit specification in description either:")'); call goErr
        write (gol,'("  filename             : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description : ",a)') trim(description); call goErr
        TRACEBACK; status=1; return
      end if
    !~ other error ...
    else
      ! some error ...
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_VarUnits
  
  
  ! *
  
  ! Return string attribute value ;
  ! use  description='' for global attribute:

  subroutine File_NC_Get_Var_Attr_s( self, description, aname, avalue, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
    use NetCDF, only : NF90_GLOBAL
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    character(len=*), intent(in)      ::  aname
    character(len=*), intent(out)     ::  avalue
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_Attr_s'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    if ( len_trim(description) > 0 ) then
      call File_NC_Inq_VarID( self, description, varid, status )
      IF_NOT_OK_RETURN(status=1)
    else
      varid = NF90_GLOBAL
    end if
    
    ! get attribute value:
    status = NF90_Get_Att( self%ncid, varid, trim(aname), avalue )
    !~ attribute not found ?
    if ( status == NF90_ENOTATT ) then
      write (gol,'("attribute not found:")'); call goErr
      write (gol,'("  filename             : ",a)') trim(self%filename); call goErr
      write (gol,'("  variable description : ",a)') trim(description); call goErr
      write (gol,'("  attribute name       : ",a)') trim(aname); call goErr
      TRACEBACK; status=1; return
    !~ error?
    else if ( status /= NF90_NOERR ) then
      ! some error ...
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_Attr_s
  
  ! *
  
  ! return time value from attribute ;
  ! use  description='' for global attribute:

  subroutine File_NC_Get_Var_Attr_t( self, description, aname, t, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
    use NetCDF, only : NF90_GLOBAL
    use GO    , only : TDate, goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    character(len=*), intent(in)      ::  aname
    type(TDate), intent(out)          ::  t
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_Attr_t'
    
    ! --- local ----------------------------------
    
    character(len=64)   ::  line
    integer             ::  varid
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    if ( len_trim(description) > 0 ) then
      call File_NC_Inq_VarID( self, description, varid, status )
      IF_NOT_OK_RETURN(status=1)
    else
      varid = NF90_GLOBAL
    end if
    
    ! get attribute value:
    status = NF90_Get_Att( self%ncid, varid, trim(aname), line )
    !~ attribute not found ?
    if ( status == NF90_ENOTATT ) then
      write (gol,'("attribute not found:")'); call goErr
      write (gol,'("  filename             : ",a)') trim(self%filename); call goErr
      write (gol,'("  variable description : ",a)') trim(description); call goErr
      write (gol,'("  attribute name       : ",a)') trim(aname); call goErr
      TRACEBACK; status=1; return
    !~ error?
    else if ( status /= NF90_NOERR ) then
      ! some error ...
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract time:
    call goReadFromLine( line, t, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_Attr_t


  ! *
  
  
  ! Return missing value, status -1 if not found.

  subroutine File_NC_Inq_VarMissing( self, varid, missing_value, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  missing_value
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_VarMissing'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! try to get first packing variable:
    status = NF90_Get_Att( self%ncid, varid, 'missing_value', missing_value )
    ! no error, missing value available
    if ( status == NF90_NOERR ) then
      ! missing value available,
      ! remain in output argument
      
    ! attribute not found ?
    else if ( status == NF90_ENOTATT ) then
      ! no missing value, set default values:
      missing_value = -999.9
      ! warning status:
      status = -1; return
    !
    else
      ! some error ...
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_VarMissing
  
  
  ! *
  
  
  ! Return coordinates attribute, status -1 if not found:
  !
  !   float HGT(Time, south_north, west_east) ;
  !     HGT:coordinates = "XLONG XLAT" ;
  !

  subroutine File_NC_Inq_VarCoordinates( self, varid, coordinates, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    integer, intent(in)               ::  varid
    character(len=*), intent(out)     ::  coordinates
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_VarCoordinates'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! try to get attribute:
    status = NF90_Get_Att( self%ncid, varid, 'coordinates', coordinates )
    ! no error, missing value available
    if ( status == NF90_NOERR ) then
      ! missing value available,
      ! remain in output argument
      
    ! attribute not found ?
    else if ( status == NF90_ENOTATT ) then
      ! no missing value, set default values:
      coordinates = ''
      ! warning status:
      status = -1; return
    !
    else
      ! some error ...
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_VarCoordinates
 

  !  *


  subroutine File_NC_Get_Var_i_1d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_i_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_i_1d
  
  ! *
  
  subroutine File_NC_Get_Var_i_2d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_i_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_i_2d
  
  
  ! *
  
  
  subroutine File_NC_Get_Var_i_3d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:,:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_i_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_i_3d


  ! *
  

  subroutine File_NC_Get_Var_r_1d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_r_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_r_1d
  
  ! *
  
  subroutine File_NC_Get_Var_r_2d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_r_2d'
    
    ! --- local ----------------------------------
    
    logical                 ::  combine
    integer, allocatable    ::  xstart(:)
    integer, allocatable    ::  xcount(:)
    integer                 ::  x1, nx
    integer                 ::  varid
    integer, allocatable    ::  dimids(:)
    real                    ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check start:
    if ( any((/present(start),present(count)/)) .and. &
         (.not. all((/present(start),present(count)/))) ) then
      write (gol,'("specify both start and count")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! combine slabs?
    combine = .false.
    if ( present(start) ) combine = start(1) < 1
    ! switch:
    if ( combine ) then

      ! storage:
      allocate( xstart(size(start)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xcount(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dimids(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      xstart = start
      xcount = count

      ! start index:
      x1 = xstart(1)
    
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get input dimension length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(1), len=nx )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! set input start and count for first slab:
      xstart(1) = nx + x1
      xcount(1) = nx - xstart(1) + 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(1:xcount(1),:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        write (gol,*) ' start                  : ', xstart; call goErr
        write (gol,*) ' count                  : ', xcount; call goErr
        TRACEBACK; status=1; return
      end if

      ! set input start and count for second slab:
      xstart(1) = 1
      xcount(1) = size(values,1) + x1 - 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(2-x1:size(values,1),:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        write (gol,*) ' start                  : ', xstart; call goErr
        write (gol,*) ' count                  : ', xcount; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! clear:
      deallocate( xstart, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xcount, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)

    else

      ! read:
      status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        if ( present(start) ) then
          write (gol,*) ' start                  : ', start; call goErr
        end if
        if ( present(count) ) then
          write (gol,*) ' count                  : ', count; call goErr
        end if
        TRACEBACK; status=1; return
      end if
      
    end if
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_r_2d
  
  ! *
  
  subroutine File_NC_Get_Var_r_3d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:,:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_r_3d'
    
    ! --- local ----------------------------------
    
    logical                 ::  combine
    integer, allocatable    ::  xstart(:)
    integer, allocatable    ::  xcount(:)
    integer                 ::  x1, nx
    integer                 ::  varid
    integer, allocatable    ::  dimids(:)
    real                    ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call File_NC_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check start:
    if ( any((/present(start),present(count)/)) .and. &
         (.not. all((/present(start),present(count)/))) ) then
      write (gol,'("specify both start and count")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! combine slabs?
    combine = .false.
    if ( present(start) ) combine = start(1) < 1
    ! switch:
    if ( combine ) then

      ! storage:
      allocate( xstart(size(start)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xcount(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dimids(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      xstart = start
      xcount = count

      ! start index:
      x1 = xstart(1)
    
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get input dimension length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(1), len=nx )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! set input start and count for first slab:
      xstart(1) = nx + x1
      xcount(1) = nx - xstart(1) + 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(1:xcount(1),:,:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        write (gol,*) ' start                  : ', xstart; call goErr
        write (gol,*) ' count                  : ', xcount; call goErr
        TRACEBACK; status=1; return
      end if

      ! set input start and count for second slab:
      xstart(1) = 1
      xcount(1) = size(values,1) + x1 - 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(2-x1:size(values,1),:,:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        write (gol,*) ' start                  : ', xstart; call goErr
        write (gol,*) ' count                  : ', xcount; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! clear:
      deallocate( xstart, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xcount, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)

    else
    
      ! read:
      status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
      if ( status /= NF90_NOERR ) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("while reading:")'); call goErr
        write (gol,'("  file name              : ",a)') trim(self%filename); call goErr
        write (gol,'("  variable description   : ",a)') trim(description); call goErr
        if ( present(start) ) then
          write (gol,*) ' start                  : ', start; call goErr
        end if
        if ( present(count) ) then
          write (gol,*) ' count                  : ', count; call goErr
        end if
        TRACEBACK; status=1; return
      end if
      
    end if

    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_r_3d


  ! *
  

  subroutine File_NC_Get_Var_t_1d( self, description, tt, status, &
                                     start, count )

    use GO, only : TDate, ExpandTime
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    type(TDate), intent(out)          ::  tt(:)
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_t_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  nt
    real, allocatable   ::  values(:)
    character(len=64)   ::  units
    integer             ::  it
    
    ! --- begin ----------------------------------
    
    ! size:
    nt = size(tt)
    ! storage:
    allocate( values(nt), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    call self%Get_Var( description, values, units, status, &
                           start=start, count=count )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over elements:
    do it = 1, nt
      ! exapand:
      call ExpandTIme( values(it), units, 'standard', tt(it), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_t_1d


  ! *
  

  subroutine File_NC_Get_Var_t_2d( self, description, tt, status, &
                                     start, count )

    use GO, only : TDate, ExpandTime
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  description
    type(TDate), intent(out)          ::  tt(:,:)
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_t_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  nt, nv
    real, allocatable   ::  values(:,:)
    character(len=64)   ::  units
    integer             ::  it, iv
    
    ! --- begin ----------------------------------
    
    ! size:
    nt = size(tt,1)
    nv = size(tt,2)
    ! storage:
    allocate( values(nt,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    call self%Get_Var( description, values, units, status, &
                           start=start, count=count )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over elements:
    do it = 1, nt
      do iv = 1, nv
        ! exapand:
        call ExpandTIme( values(it,iv), units, 'standard', tt(it,iv), status )
        IF_NOT_OK_RETURN(status=1)
      end do ! iv
    end do ! it
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_t_2d
  
  
  ! ***
  
  
  !
  ! Return record index in time axis.
  !
  ! Optional climatology mask,
  ! specify which year/month/day/hour/minute/second always match,
  ! for example to only check month and hour:
  !    climat=(/ .true., .false., .true., .false., .true., .true. /)
  !
  ! Optional keys:
  !    time_name=XTIME
  !
  
  subroutine File_NC_Inq_TimeRecord( self, name, t, irec, status, &
                                      climat, keys )

    use GO, only : goVarValue
    use GO, only : TDate, ExpandTime, Get, operator(==),operator(<=), Set
    use GO, only : wrtgol
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(inout)   ::  self
    character(len=*), intent(in)      ::  name
    type(TDate), intent(in)           ::  t
    integer, intent(out)              ::  irec
    integer, intent(out)              ::  status

    logical, intent(in), optional           ::  climat(6)
    character(len=*), intent(in), optional  ::  keys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_TimeRecord'
    
    ! --- local ----------------------------------
    
    character(len=64)       ::  var_name
    integer                 ::  varid
    integer                 ::  ntime
    integer                 ::  itime
    real, allocatable       ::  time(:)
    character(len=64)       ::  units
    type(TDate)             ::  t_rec
    character(len=10)       ::  label
    logical                 ::  cmask(6)
    integer                 ::  time6(6), time6_rec(6)
    
    ! --- begin ----------------------------------
    
    ! climat mask:
    cmask = .false.
    if ( present(climat) ) cmask = climat
    
    ! time variable name:
    var_name = trim(name)
    ! optional reset:
    if ( present(keys) ) then
      call goVarValue( keys, ';', 'time_name', '=', var_name, status )
      IF_ERROR_RETURN(status=1)
    end if

    ! number of time records:
    call self%Inquire_Dimension( name, status, length=ntime )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( time(ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read values:
    call self%Get_Var( 'var_name='//trim(var_name), time, units, status )
    IF_NOT_OK_RETURN(status=1)

    ! target values:
    call Get( t, time6=time6 )
          
    ! init index:
    irec = -999
    ! loop over time records in file:
    do itime = 1, ntime
      ! exapand:
      call ExpandTIme( time(itime), units, 'standard', t_rec, status )
      IF_NOT_OK_RETURN(status=1)
      ! extract values:
      call Get( t_rec, time6=time6_rec  )
      ! compare for complete match,
      ! but allow values to be different when masked:
      if ( all( (time6_rec == time6) .or. cmask ) ) then
        ! store record index:
        irec = itime
        ! leave loop over records:
        exit
      end if ! match?
    end do ! time records in file
    
    ! check ...
    if ( irec < 0 ) then
      write (gol,'("time record not found:")'); call goErr
      write (gol,'("  filename    : ",a)') trim(self%filename); call goErr
      call wrtgol( '  requested   : ', t ); call goErr
      if ( any(cmask) ) then
        write (gol,'("  climat mask : ",6l2)') cmask; call goErr
      end if
      write (gol,'("  variable    : ",a)') trim(name); call goErr
      write (gol,'("    units     : ",a)') trim(units); call goErr
      write (gol,'("  records     : ",i0)') ntime; call goErr
      do itime = 1, ntime
        call ExpandTIme( time(itime), units, 'standard', t_rec, status )
        IF_NOT_OK_RETURN(status=1)
        ! label "  nnnn  : " ; when changing, also change length!
        write (label,'("  ",i4,"  : ")') itime
        call wrtgol( label, t_rec ); call goErr
      end do
      TRACEBACK; status = -1; return
    end if
    
    ! clear:
    deallocate( time )
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_TimeRecord
  
  
  ! ***
  
  
  !
  ! Return record index in time axis
  !
  ! Optional climatology mask,
  ! specify which year/month/day/hour/minute/second always match,
  ! for example to only check month and hour:
  !    climat=(/ .true., .false., .true., .false., .true., .true. /)
  !
  ! Optional keys:
  !    time_name=XTIME
  !
  
  subroutine File_NC_Inq_TimeBoundsRecord( self, name, tt, irec, tt_rec, status, &
                                             climat, keys )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att

    use GO, only : goVarValue
    use GO, only : TDate, ExpandTime, operator(<=), operator(.in.), Set
    use GO, only : wrtgol
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(inout)   ::  self
    character(len=*), intent(in)      ::  name
    type(TDate), intent(in)           ::  tt(2)
    integer, intent(out)              ::  irec
    type(TDate), intent(out)          ::  tt_rec(2)
    integer, intent(out)              ::  status

    logical, intent(in), optional           ::  climat(6)
    character(len=*), intent(in), optional  ::  keys

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inq_TimeBoundsRecord'
    
    ! --- local ----------------------------------
    
    character(len=64)       ::  var_name
    integer                 ::  varid
    integer                 ::  ntime
    integer                 ::  itime
    character(len=64)       ::  var_name_bounds
    real, allocatable       ::  time(:)  ! (ntime)
    real, allocatable       ::  time_bounds(:,:)  ! (2,ntime)
    character(len=64)       ::  units
    character(len=8)        ::  label
    logical                 ::  cmask(6)
    integer                 ::  dyear
    
    ! --- begin ----------------------------------
    
    ! climat mask:
    cmask = .false.
    if ( present(climat) ) cmask = climat
    ! not yet ...
    if ( any(cmask) ) then
      write (gol,'("climatology not supported yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! time variable name:
    var_name = trim(name)
    ! optional reset:
    if ( present(keys) ) then
      call goVarValue( keys, ';', 'time_name', '=', var_name, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get variable id of time variable:
    call self%Inq_VarID( 'var_name='//trim(var_name), varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! time variable should have 'bounds' attribute
    ! with name of time bounds variable:
    status = NF90_Get_Att( self%ncid, varid, 'bounds', var_name_bounds )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! number of time records:
    call self%Inquire_Dimension( name, status, length=ntime )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( time(ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( time_bounds(2,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read values:
    call self%Get_Var( 'var_name='//trim(name), time, units, status )
    IF_NOT_OK_RETURN(status=1)
    call self%Get_Var( 'var_name='//trim(var_name_bounds), time_bounds, units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( any( time_bounds(2,:) <= time_bounds(1,:) ) &
         .or. any( time < time_bounds(1,:) ) &
         .or. any( time > time_bounds(2,:) )   ) then
      write (gol,'("time bounds should specify non-zero intervals")'); call goErr
      write (gol,'("including the time value:")'); call goPr
      do itime = 1, ntime
        write (gol,'("  ",i4,f16.4," [",f16.4,",",f16.4,"]")') itime, time(itime), time_bounds(:,itime); call goPr
      end do
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( ntime > 1 ) then
      if ( any( time_bounds(1,2:ntime) /= time_bounds(2,1:ntime-1) ) ) then
        write (gol,'("time bounds should specify a sequence of intervals:")'); call goErr
        do itime = 1, ntime
          write (gol,'("  ",i4,f16.4," [",f16.4,",",f16.4,"]")') itime, time(itime), time_bounds(:,itime); call goPr
        end do
        TRACEBACK; status=1; return
      end if
    end if

    ! init index:
    irec = -999
    ! loop:
    do itime = 1, ntime
      ! expand:
      call ExpandTime( time_bounds(1,itime), units, 'standard', tt_rec(1), status )
      IF_NOT_OK_RETURN(status=1)
      call ExpandTime( time_bounds(2,itime), units, 'standard', tt_rec(2), status )
      IF_NOT_OK_RETURN(status=1)
      ! some masks are supported ...
      if ( any(cmask) ) then
        ! do not check on year?
        if ( all( cmask .eqv. (/.true.,.false.,.false.,.false.,.false.,.false./) ) ) then
          ! reset year in source record interval:
          dyear = tt_rec(2)%year - tt_rec(1)%year
          call Set( tt_rec(1), year=tt(1)%year )
          call Set( tt_rec(2), year=tt(1)%year+dyear )
          ! compare:
          if ( tt .in. tt_rec ) then
            irec = itime
            exit
          end if
        else
          write (gol,'("climatology not supported yet: ",6l2)') cmask; call goErr
          TRACEBACK; status=1; return
        end if
      else
        ! compare:
        if ( tt .in. tt_rec ) then
          irec = itime
          exit
        end if
      end if
    end do
    
    ! check ...
    if ( irec < 0 ) then
      write (gol,'("time record not found:")'); call goErr
      write (gol,'("  filename    : ",a)') trim(self%filename); call goErr
      call wrtgol( '  requested   : ', tt(1), ' - ', tt(2) ); call goErr
      if ( any(cmask) ) then
        write (gol,'("  climat mask : ",6l2)') cmask; call goErr
      end if
      write (gol,'("  variable    : ",a)') trim(name); call goErr
      write (gol,'("    units     : ",a)') trim(units); call goErr
      write (gol,'("   records    : ")'); call goErr
      do itime = 1, ntime
        call ExpandTime( time_bounds(1,itime), units, 'standard', tt_rec(1), status )
        IF_NOT_OK_RETURN(status=1)
        call ExpandTime( time_bounds(2,itime), units, 'standard', tt_rec(2), status )
        IF_NOT_OK_RETURN(status=1)
        write (label,'("  ",i2.2,"  : ")') itime
        call wrtgol( label, tt_rec ); call goErr
      end do
      TRACEBACK; status = -1; return
    end if
    
    ! clear:
    deallocate( time, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( time_bounds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_TimeBoundsRecord
  
  
  ! ***


  subroutine File_Nc_Get_Levs_bnds( self, levs, status )
  
    use GO, only : goReadFromLine

    use NetCDF, only : NF90_INQ_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)        ::  self
    class(T_Levs_Hyb), intent(out)      ::  levs
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Get_Levs_bnds'
    
    ! names:
    character(len=*), parameter  ::  dimname_lev  = 'lev'
    character(len=*), parameter  ::  dimname_bnds = 'bnds'
    character(len=*), parameter  ::  varext_bnds  = '_'//dimname_bnds
    
    ! --- local ----------------------------------

    integer                 ::  dimid, varid
    integer                 ::  nlev, nv
    real, allocatable       ::  ap(:), b(:)
    real, allocatable       ::  ap_bnds(:,:), b_bnds(:,:)
    real                    ::  p0
    character(len=256)      ::  formula_terms, fterms
    character(len=32)       ::  cname, cvalue
    character(len=32)       ::  varname_ap, varname_a, varname_b, varname_ps, varname_p0
    
    ! --- begin ----------------------------------
    
    ! dimension with number of levels:
    status = NF90_INQ_DimID( self%ncid, dimname_lev, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! count:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=nlev )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! dimension with number of boundaries:
    status = NF90_INQ_DimID( self%ncid, dimname_bnds, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! count:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=nv )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( ap(nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( b (nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( ap_bnds(nv,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( b_bnds (nv,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! variable with level indices:
    status = NF90_INQ_VarID( self%ncid, dimname_lev, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read formula description:
    status = NF90_Get_Att( self%ncid, varid, 'formula_terms', formula_terms )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! init values:
    varname_a  = 'None'
    varname_ap = 'None'
    varname_b  = 'None'
    varname_ps = 'None'
    varname_p0 = 'None'
    ! copy:
    fterms = formula_terms
    ! loop:
    do
      ! all processed ?
      if ( len_trim(fterms) == 0 ) exit
      ! extract first name and value:
      !    ap: ap b: b ps: ps
      call goReadFromLine( fterms, cname, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      call goReadFromLine( fterms, cvalue, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      ! switch:
      select case ( trim(cname) )
        case ( 'ap:' ) ; varname_ap = trim(cvalue)
        case ( 'a:'  ) ; varname_a  = trim(cvalue)
        case ( 'p0:' ) ; varname_p0 = trim(cvalue)
        case ( 'b:'  ) ; varname_b  = trim(cvalue)
        case ( 'ps:' ) ; varname_ps = trim(cvalue)
        case default
          write (gol,'("unsupported name `",a,"` in formula_terms `",a,"`")') &
                          trim(cname), trim(formula_terms); call goErr
          TRACEBACK; status=1; return
      end select
    end do
    
    ! ap or a*p0 ?
    if ( trim(varname_ap) /= 'None' ) then

      ! variable with coeffs:
      status = NF90_INQ_VarID( self%ncid, trim(varname_ap), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read:
      status = NF90_Get_Var( self%ncid, varid, ap )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! variable with coeff bounds:
      status = NF90_INQ_VarID( self%ncid, trim(varname_ap)//varext_bnds, varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read:
      status = NF90_Get_Var( self%ncid, varid, ap_bnds )
      IF_NF90_NOT_OK_RETURN(status=1)

    else if ( (trim(varname_a) /= 'None') .and. (trim(varname_p0) /= 'None') ) then

      ! variable with p0:
      status = NF90_INQ_VarID( self%ncid, trim(varname_p0), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read:
      status = NF90_Get_Var( self%ncid, varid, p0 )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! variable with coeffs:
      status = NF90_INQ_VarID( self%ncid, trim(varname_a), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read 'a', fill into 'ap':
      status = NF90_Get_Var( self%ncid, varid, ap )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! multipy:
      ap = ap * p0

      ! variable with coeff bounds:
      status = NF90_INQ_VarID( self%ncid, trim(varname_a)//varext_bnds, varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read 'a', fill into 'ap':
      status = NF90_Get_Var( self%ncid, varid, ap_bnds )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! multipy:
      ap_bnds = ap_bnds * p0

    else
      write (gol,'("either ap or a and p0 should be defined in formula terms: ",a)') &
                      trim(formula_terms); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! variable with coeffs:
    status = NF90_INQ_VarID( self%ncid, trim(varname_b), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, b )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! variable with coeff bounds:
    status = NF90_INQ_VarID( self%ncid, trim(varname_b)//varext_bnds, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, b_bnds )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! init level definition:
    call levs%Init( ap, b, ap_bnds, b_bnds, status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( ap, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( b, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( ap_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( b_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Nc_Get_Levs_bnds

  ! *
  

  !
  ! Optional 'keys' line with 'key=value' pairs seperated by ';'
  ! to define names of dimensions or variables that are not
  ! included in the file because cf-conventions are not followed exactly,
  ! or because cf-conventions do not prescribe the information needed:
  !
  ! ~ Define one of these keys to use the (2D) 'hya_bnds' and 'hyb_bnds' arrays;
  !   the name itself is ignored, only tested if it is defined:
  !
  !     hyb_ap_bnds_name=hya_bnds
  !     hyb_b_bnds_name=hyb_bnds
  !
  ! ~ Otherwise, coeffs are expected to be in 1D arrays only.
  !   First check if names for the interface arrays are provided:
  !
  !     hyb_api_name=hyai  : hybride ap levels (Pa) at interfaces
  !     hyb_bi_name=hybi   : hybride b levels (1) at interfaces
  !
  !   or otherwise the name of a half-level coordinate in which
  !   formula_terms are expected to be found:
  !
  !     hlevel_name=levi   : half-level dimension name
  !
  !   Then also try if variable names for the mid-level coeff are defined;
  !   if not, then the mid levels are computed as average of the interfaces:
  !
  !     hyb_ap_name=hyam   : hybride ap levels (Pa) at mid levels
  !     hyb_b_name=hybm    : hybride b levels (1) at mid levels
  !
  
  subroutine File_Nc_Get_Levs_h( self, vid, keys, levs, status )
  
    use GO, only : goVarValue, goReadFromLine

    use NetCDF, only : NF90_INQ_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    class(T_File_Nc), intent(in)            ::  self
    integer, intent(in)                     ::  vid
    character(len=*)                        ::  keys
    class(T_Levs_Hyb), intent(out)          ::  levs
    integer, intent(out)                    ::  status
 
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Get_Levs_h'
      
    ! --- local ----------------------------------

    character(len=32)       ::  varname_ap_bnds, varname_b_bnds
    integer                 ::  ndims
    integer, allocatable    ::  dimids(:)
    integer                 ::  dimid
    integer                 ::  varid
    character(len=32)       ::  dimname
    integer                 ::  nlev
    real, allocatable       ::  ap(:), b(:)
    real, allocatable       ::  ap_i(:), b_i(:)
    character(len=256)      ::  formula_terms
    character(len=32)       ::  cname, cvalue
    character(len=32)       ::  varname_ap, varname_b, varname_ps
    integer                 ::  ncoef
    real                    ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! quick test ...
    varname_ap_bnds = 'None'
      call goVarValue( keys, ';', 'hyb_ap_bnds_name', '=', varname_ap_bnds, status )
      IF_ERROR_RETURN(status=1)
    varname_b_bnds = 'None'
      call goVarValue( keys, ';', 'hyb_b_bnds_name', '=', varname_b_bnds, status )
      IF_ERROR_RETURN(status=1)
      
    ! switch:
    if ( (trim(varname_ap_bnds) /= 'None') .or. (trim(varname_b_bnds) /= 'None') ) then
    
      ! ** coeff in '_bnds' arrays
      
      ! read using pre-defined names
      call self%Get_Levs( levs, status )
      IF_ERROR_RETURN(status=1)
      
    else
    
      ! ** coeff in 'interface' arrays
      
      ! ~ mid levels

      ! number of levels explicitly provided?
      nlev = -999
        call goVarValue( keys, ';', 'nlev', '=', nlev, status )
        IF_ERROR_RETURN(status=1)

      ! nlev to be obtained from dimensions?
      if ( nlev < 0 ) then

        ! level dimension name explicitly provided?
        dimname = 'None'
          call goVarValue( keys, ';', 'level_name', '=', dimname, status )
          IF_ERROR_RETURN(status=1)

        ! get level dimension name from variable dimensions?
        if ( trim(dimname) == 'None' ) then
        
          ! no idea ...
          write (gol,'("could not guess level dimension, provide `level_name` keyword in variable description")'); call goErr
          TRACEBACK; status=1; return 

          !... This did not not work, no level dimension in (x,y,t) surface pressure
          !! name explicitly provided?
          !! get  number of dimensions:
          !status = NF90_Inquire_Variable( self%ncid, vid, ndims=ndims )
          !IF_NF90_NOT_OK_RETURN(status=1)
          !! storage:
          !allocate( dimids(ndims), stat=status )
          !IF_NOT_OK_RETURN(status=1)
          !! get dimension id's:
          !status = NF90_Inquire_Variable( self%ncid, vid, dimids=dimids )
          !IF_NF90_NOT_OK_RETURN(status=1)

          !! check ..
          !if ( (idim < 0) .or. (idim > ndims) ) then
          !  write (gol,'("level dimension number shoulb be in 1,..,",i0," but requested ",i0)') ndims, idim; call goErr
          !  TRACEBACK; status=1; return
          !end if

          !! name and length of level dimension:
          !status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dimname, len=nlev )
          !IF_NF90_NOT_OK_RETURN(status=1)

          !! clear:
          !deallocate( dimids, stat=status )
          !IF_NOT_OK_RETURN(status=1)

        else

          ! get dimension id:
          status = NF90_INQ_DimID( self%ncid, trim(dimname), dimid )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! length of level dimension:
          status = NF90_Inquire_Dimension( self%ncid, dimid, len=nlev )
          IF_NF90_NOT_OK_RETURN(status=1)

        end if  ! dimname from variable dimensions

      end if  ! nlev from dimension

      ! ~ interface coeff

      ! storage:
      allocate( ap_i(0:nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( b_i (0:nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! names of interface coeffs taken from formula terms,
      ! unless specified explicitly by arguments:
      varname_ap = 'None'
        call goVarValue( keys, ';', 'hyb_api_name', '=', varname_ap, status )
        IF_ERROR_RETURN(status=1)
      varname_b = 'None'
        call goVarValue( keys, ';', 'hyb_bi_name', '=', varname_b, status )
        IF_ERROR_RETURN(status=1)

      ! not defined yet?
      if ( (varname_ap == 'None') .or. (varname_b == 'None') ) then

        ! halflevel name is not provided in cf-conventions yet,
        ! so should be provided as key; default guess is name used by leip:
        dimname = 'hlevel'
          call goVarValue( keys, ';', 'hlevel_name', '=', dimname, status )
          IF_ERROR_RETURN(status=1)

        ! variable with level indices:
        status = NF90_INQ_VarID( self%ncid, trim(dimname), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read formula description:
        status = NF90_Get_Att( self%ncid, varid, 'formula_terms', formula_terms )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! init values:
        varname_ap = 'None'
        varname_b  = 'None'
        varname_ps = 'None'
        ! loop:
        do
          ! all processed ?
          if ( len_trim(formula_terms) == 0 ) exit
          ! extract first name and value:
          !    ap: ap b: b ps: ps
          call goReadFromLine( formula_terms, cname, status, sep=' ' )
          IF_NOT_OK_RETURN(status=1)
          call goReadFromLine( formula_terms, cvalue, status, sep=' ' )
          IF_NOT_OK_RETURN(status=1)
          ! switch:
          select case ( trim(cname) )
            case ( 'ap:' ) ; varname_ap = trim(cvalue)
            case ( 'b:'  ) ; varname_b  = trim(cvalue)
            case ( 'ps:' ) ; varname_ps = trim(cvalue)
            case default
              write (gol,'("unsupported name `",a,"` in formula_terms")') trim(cname); call goErr
              TRACEBACK; status=1; return
          end select
        end do

      end if  ! hyb_ap/hyb_b or formula terms

      ! variable with coeffs:
      status = NF90_INQ_VarID( self%ncid, trim(varname_ap), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! length, in some files more coeff are stored than needed:
      call self%Inq_Var( varid, status, dim=1, dim_length=ncoef )
      IF_NOT_OK_RETURN(status=1)
      ! read coeffs near surface:
      status = NF90_Get_Var( self%ncid, varid, ap_i, start=(/ncoef-nlev/), count=(/nlev+1/) )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! packed?
      call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
      IF_ERROR_RETURN(status=1)
      if ( status == 0 ) then
        ! unpack:
        ap_i = add_offset + scale_factor * ap_i
      end if

      ! variable with coeffs:
      status = NF90_INQ_VarID( self%ncid, trim(varname_b), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! length, in some files more coeff are stored than needed:
      call self%Inq_Var( varid, status, dim=1, dim_length=ncoef )
      IF_NOT_OK_RETURN(status=1)
      ! read coeffs near surface:
      status = NF90_Get_Var( self%ncid, varid, b_i, start=(/ncoef-nlev/), count=(/nlev+1/) )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! packed?
      call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
      IF_ERROR_RETURN(status=1)
      if ( status == 0 ) then
        ! unpack:
        b_i = add_offset + scale_factor * b_i
      end if

      ! ~ mid level coeff

      ! storage:
      allocate( ap(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( b (nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! names of interface coeffs taken from formula terms,
      ! unless specified explicitly by arguments:
      varname_ap = 'None'
        call goVarValue( keys, ';', 'hyb_ap_name', '=', varname_ap, status )
        IF_ERROR_RETURN(status=1)
      varname_b = 'None'
        call goVarValue( keys, ';', 'hyb_b_name', '=', varname_b, status )
        IF_ERROR_RETURN(status=1)

      !! not defined yet?
      !if ( (varname_ap == 'None') .or. (varname_b == 'None') ) then
      !
      !  ! variable with level indices, first try if name is provided:
      !  status = NF90_INQ_VarID( self%ncid, trim(dimname), varid )
      !  IF_NF90_NOT_OK_RETURN(status=1)
      !  ! read formula description:
      !  status = NF90_Get_Att( self%ncid, varid, 'formula_terms', formula_terms )
      !  IF_NF90_NOT_OK_RETURN(status=1)
      !  ! init values:
      !  varname_ap = 'None'
      !  varname_b  = 'None'
      !  varname_ps = 'None'
      !  ! loop:
      !  do
      !    ! all processed ?
      !    if ( len_trim(formula_terms) == 0 ) exit
      !    ! extract first name and value:
      !    !    ap: ap b: b ps: ps
      !    call goReadFromLine( formula_terms, cname, status, sep=' ' )
      !    IF_NOT_OK_RETURN(status=1)
      !    call goReadFromLine( formula_terms, cvalue, status, sep=' ' )
      !    IF_NOT_OK_RETURN(status=1)
      !    ! switch:
      !    select case ( trim(cname) )
      !      case ( 'ap:' ) ; varname_ap = trim(cvalue)
      !      case ( 'b:'  ) ; varname_b  = trim(cvalue)
      !      case ( 'ps:' ) ; varname_ps = trim(cvalue)
      !      case default
      !        write (gol,'("unsupported name `",a,"` in formula_terms")') trim(cname); call goErr
      !        TRACEBACK; status=1; return
      !    end select
      !  end do
      !  
      !end if  ! hyb_ap/hyb_b or formula terms

      ! compute?
      if ( varname_ap == 'None' ) then
        ! average value for full levels:
        ap = 0.5 * ( ap_i(0:nlev-1) + ap_i(1:nlev) )
      else
        ! variable with coeffs:
        status = NF90_INQ_VarID( self%ncid, trim(varname_ap), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! length, in some files more coeff are stored than needed:
        call self%Inq_Var( varid, status, dim=1, dim_length=ncoef )
        IF_NOT_OK_RETURN(status=1)
        ! read coeffs near surface:
        status = NF90_Get_Var( self%ncid, varid, ap, start=(/ncoef+1-nlev/), count=(/nlev/) )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! packed?
        call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
        IF_ERROR_RETURN(status=1)
        if ( status == 0 ) then
          ! unpack:
          ap = add_offset + scale_factor * ap
        end if

      end if

      ! compute?
      if ( varname_b == 'None' ) then
        ! average value for full levels:
        b = 0.5 * ( b_i(0:nlev-1) + b_i(1:nlev) )
      else
        ! variable with coeffs:
        status = NF90_INQ_VarID( self%ncid, trim(varname_b), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! length, in some files more coeff are stored than needed:
        call self%Inq_Var( varid, status, dim=1, dim_length=ncoef )
        IF_NOT_OK_RETURN(status=1)
        ! read coeffs near surface:
        status = NF90_Get_Var( self%ncid, varid, b, start=(/ncoef+1-nlev/), count=(/nlev/) )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! packed?
        call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
        IF_ERROR_RETURN(status=1)
        if ( status == 0 ) then
          ! unpack:
          b = add_offset + scale_factor * b
        end if
      end if

      !~ definition

      ! init level definition:
      call levs%Init( ap, b, ap_i, b_i, status )
      IF_NOT_OK_RETURN(status=1)

      !~ done

      ! clear:
      deallocate( ap, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( b, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( ap_i, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( b_i, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if   ! 'bnds' or 'interface' arrays
    
    ! ok
    status = 0
        
  end subroutine File_Nc_Get_Levs_h
  
  

end module C3PO_File_Nc
