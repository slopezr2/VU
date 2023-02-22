!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
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
  public  ::  T_File_Nc
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Nc'

  ! max number of dimensions:
  integer, parameter  ::  maxdim = 7


  ! --- types ----------------------------------------
  
  type T_File_Nc
    character(len=1024)     ::  filename
    integer                 ::  ncid
    integer                 ::  unlimitedDimId
    integer                 ::  formatNum
    ! level i/o
    integer                 ::  dimid_lev
    integer                 ::  varid_lev
    integer                 ::  varid_hya
    integer                 ::  varid_hyb
    integer                 ::  dimid_levi
    integer                 ::  varid_levi
    integer                 ::  varid_hyai
    integer                 ::  varid_hybi
  contains
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
    procedure   ::  Inquire_Variable     => File_NC_Inquire_Variable
    procedure   ::  Define_Variable_Copy => File_NC_Define_Variable_Copy
    procedure   ::  Copy_Atts            => File_NC_Copy_Atts
    procedure   ::  Extend_History       => File_NC_Extend_History
    procedure   ::  Get_VarID            => File_NC_Get_VarID
    !
    procedure   ::                          File_NC_Get_Var_1d
    procedure   ::                          File_NC_Get_Var_2d
    procedure   ::                          File_NC_Get_Var_3d
    procedure   ::                          File_NC_Get_Var_4d
    generic     ::  Get_Var              => File_NC_Get_Var_1d, &
                                            File_NC_Get_Var_2d, &
                                            File_NC_Get_Var_3d, &
                                            File_NC_Get_Var_4d
    !
    procedure   ::                          File_NC_Put_Var_3d
    procedure   ::                          File_NC_Put_Var_4d
    generic     ::  Put_Var              => File_NC_Put_Var_3d, &
                                            File_NC_Put_Var_4d
    !
    procedure   ::  Inq_TimeRecord       => File_NC_Inq_TimeRecord
    !
    procedure   ::  Get_Levs             => File_Nc_Get_Levs
    procedure   ::  Def_Levs             => File_Nc_Def_Levs
    procedure   ::  Put_Levs             => File_Nc_Put_Levs
    procedure   ::  Copy_Variable_CoarsenLevs => File_Nc_Copy_Variable_CoarsenLevs
  end type T_File_Nc

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif
  

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
    
    XTYPE(T_File_Nc), intent(out)   ::  self
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
      write (gol,'("file not found : ",a)') trim(self%filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! open file for reading:
    status = NF90_Open( trim(filename), NF90_NOWRITE, self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! info ...
    status = NF90_Inquire( self%ncid, unlimitedDimId=self%unlimitedDimId, &
                             formatNum=self%formatNum )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine File_NC_Open


  ! ***


  subroutine File_NC_Create( self, filename, status, formatNum )
  
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_CLOBBER
    use NetCDF, only : NF90_CLASSIC_MODEL
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(out)    ::  self
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status
    integer, intent(in), optional   ::  formatNum

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Create'
    
    ! --- local ----------------------------------
    
    integer     ::  cmode
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! optional:
    if ( present(formatNum) ) then
      self%formatNum = formatNum
    else
      self%formatNum = NF90_CLASSIC_MODEL
    end if

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
    
    XTYPE(T_File_Nc), intent(inout)  ::  self
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
    
    XTYPE(T_File_Nc), intent(inout)  ::  self
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
    
    XTYPE(T_File_Nc), intent(in)     ::  self
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
    
    XTYPE(T_File_Nc), intent(in)              ::  self
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
    
    XTYPE(T_File_Nc), intent(in)              ::  self
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
    
    XTYPE(T_File_Nc), intent(in)              ::  self
    character(len=*), intent(in)              ::  dimname
    XTYPE(T_File_Nc), intent(inout)           ::  outfile
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
    ! is this the unlimitted dimension?
    if ( dimid == self%unlimitedDimId ) then
      ! dummy length to indicate unlimited dimension:
      length = 0
    else
      ! obtain actual length:
      status = NF90_Inquire_Dimension( self%ncid, dimid, len=length )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! create same dimension in output file,
    ! returned dimension id is not used:
    status = NF90_Def_Dim( outfile%ncid, dimname, length, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Copy_Dimension


  ! ***
  
  
  !
  ! INquire variable properties.
  ! Variable id is the same as sequential number 1,..,nvar .
  !
  
  subroutine File_NC_Inquire_Variable( self, varid, status, &
                                         name, xtype, dtype )

    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName, NF90_Inquire_Attribute
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)               ::  self
    integer, intent(in)                       ::  varid
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  xtype   ! storage type
    integer, intent(out), optional            ::  dtype   ! unpacked datatype

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Inquire_Variable'
    
    ! --- local ----------------------------------

    integer             ::  nattr, iattr
    character(len=64)   ::  attrname
    
    ! --- begin ----------------------------------
    
    ! values returned by standard function:
    if ( present(name) ) then
      ! obtain info:
      status = NF90_Inquire_Variable( self%ncid, varid, name=name, xtype=xtype )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! return unpacked type?
    if ( present(dtype) ) then
      ! default is storage type, return number of attributes too:
      status = NF90_Inquire_Variable( self%ncid, varid, xtype=dtype, nAtts=nattr )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! search for packing attributes:
      do iattr = 1, nattr
        ! get name:
        status = NF90_Inq_AttName( self%ncid, varid, iattr, name=attrname )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! match?
        if ( (trim(attrname) == 'add_offset') .or. (trim(attrname) == 'scale_factor') ) then
          ! get data type of packing attributes:
          status = NF90_Inquire_Attribute( self%ncid, varid, trim(attrname), xtype=dtype )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! leave:
          exit
        end if  ! packing attribute
      end do ! attributes
    end if  ! return dtype
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inquire_Variable


  ! ***
  
  
  subroutine File_NC_Define_Variable_Copy( self, varname, outfile, varid, status )

    use NetCDF, only : NF90_Inquire_Dimension, NF90_Inq_DimID
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Def_Var

    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)              ::  self
    character(len=*), intent(in)              ::  varname
    XTYPE(T_File_Nc), intent(inout)           ::  outfile
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
    
    XTYPE(T_File_Nc), intent(in)     ::  self
    character(len=*), intent(in)    ::  varname
    XTYPE(T_File_Nc), intent(inout)  ::  outfile
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
    
    XTYPE(T_File_Nc), intent(inout)  ::  self
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
    character(len=4000)   ::  history
    character(len=4000)   ::  history_new
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
      write (history_new,'(a,": ",a,a,a)') trim(tstamp), trim(message)
    else
      write (history_new,'(a,": ",a,a,a)') trim(tstamp), trim(message), &
                                               newline, trim(history)
    end if
    
    ! write:
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, trim(attname), trim(history_new) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_NC_Extend_History
  
  
  ! ***
  
  
  !
  ! Return id of variable identified by 'standard_name' or 'var_name'
  !
  
  subroutine File_NC_Get_VarID( self, nametype, name, varid, status )

    use NetCDF, only : NF90_INQ_VarID
    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)   ::  self
    character(len=*), intent(in)      ::  nametype
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  varid
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_VarID'
    
    ! --- local ----------------------------------
    
    logical             ::  found
    integer             ::  nvar, ivar
    character(len=64)   ::  standard_name
    
    ! --- begin ----------------------------------
    
    ! set flag:
    found = .false.

    ! switch:
    select case ( nametype )
      !~
      case ( 'var_name' )
        ! id of input variable:
        status = NF90_INQ_VarID( self%ncid, name, varid )
        found = status == NF90_NOERR
      !~
      case ( 'standard_name', 'long_name' )
        ! number of variables:
        call self%Inquire( status, nvar=nvar )
        IF_NOT_OK_RETURN(status=1)
        ! loop:
        do ivar = 1, nvar
          ! variable id is number:
          varid = ivar
          ! get standard name if possible:
          status = NF90_Get_Att( self%ncid, varid, trim(nametype), standard_name )
          if ( status /= NF90_NOERR ) cycle
          ! match?
          found = trim(standard_name) == trim(name)
          ! leave ?
          if ( found ) exit
        end do
      !~
      case default
        write (gol,'("unsupported name type `",a,"`")') trim(nametype); call goErr
        TRACEBACK; status=1; return
    end select
    ! check ...
    if ( .not. found ) then
      write (gol,'("variable with `",a,"` equal to `",a,"` not found")') trim(nametype), trim(name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_VarID
  
  ! *
  
  subroutine File_NC_Get_Var_1d( self, varid, values, status, &
                                   start, count, units )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)                ::  self
    integer, intent(in)                         ::  varid
    real, intent(out)                           ::  values(:)
    integer, intent(out)                        ::  status

    integer, intent(in), optional               ::  start(:), count(:)
    character(len=*), intent(out), optional     ::  units

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_1d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! return units?
    if ( present(units) ) then
      ! get units from attribute:
      status = NF90_Get_Att( self%ncid, varid, 'units', units )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_1d
  
  ! *
  
  subroutine File_NC_Get_Var_2d( self, varid, values, status, &
                                   start, count, unpack, units )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)                ::  self
    integer, intent(in)                         ::  varid
    real, intent(out)                           ::  values(:,:)
    integer, intent(out)                        ::  status

    integer, intent(in), optional               ::  start(:), count(:)
    logical, intent(in), optional               ::  unpack
    character(len=*), intent(out), optional     ::  units

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_2d'
    
    ! --- local ----------------------------------
    
    real             ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! return units?
    if ( present(units) ) then
      ! get units from attribute:
      status = NF90_Get_Att( self%ncid, varid, 'units', units )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! unpack?
    if ( present(unpack) .and. unpack ) then
      ! packing attributes:
      status = NF90_Get_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Get_Att( self%ncid, varid, 'scale_factor', scale_factor  )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! unpack:
      values = add_offset + scale_factor * values
    end if ! unpack
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_2d
  
  ! *
  
  subroutine File_NC_Get_Var_3d( self, varid, values, status, &
                                   start, count, unpack, units )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)                ::  self
    integer, intent(in)                         ::  varid
    real, intent(out)                           ::  values(:,:,:)
    integer, intent(out)                        ::  status

    integer, intent(in), optional               ::  start(:), count(:)
    logical, intent(in), optional               ::  unpack
    character(len=*), intent(out), optional     ::  units

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_3d'
    
    ! --- local ----------------------------------
    
    real             ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! return units?
    if ( present(units) ) then
      ! get units from attribute:
      status = NF90_Get_Att( self%ncid, varid, 'units', units )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! unpack?
    if ( present(unpack) .and. unpack ) then
      ! packing attributes:
      status = NF90_Get_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Get_Att( self%ncid, varid, 'scale_factor', scale_factor  )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! unpack:
      values = add_offset + scale_factor * values
    end if ! unpack
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_3d
  
  ! *
  
  subroutine File_NC_Get_Var_4d( self, varid, values, status, &
                                   start, count, unpack, units )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)                ::  self
    integer, intent(in)                         ::  varid
    real, intent(out)                           ::  values(:,:,:,:)
    integer, intent(out)                        ::  status

    integer, intent(in), optional               ::  start(:), count(:)
    logical, intent(in), optional               ::  unpack
    character(len=*), intent(out), optional     ::  units

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_Var_4d'
    
    ! --- local ----------------------------------
    
    real             ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! return units?
    if ( present(units) ) then
      ! get units from attribute:
      status = NF90_Get_Att( self%ncid, varid, 'units', units )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! unpack?
    if ( present(unpack) .and. unpack ) then
      ! packing attributes:
      status = NF90_Get_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Get_Att( self%ncid, varid, 'scale_factor', scale_factor  )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! unpack:
      values = add_offset + scale_factor * values
    end if ! unpack
    
    ! ok
    status = 0
    
  end subroutine File_NC_Get_Var_4d
  
  
  ! ***
  
  ! *
  
  subroutine File_NC_Put_Var_3d( self, varid, values, status, &
                                   start, count, pack )

    use NetCDF, only : NF90_SHORT
    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)   ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  values(:,:,:)
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    logical, intent(in), optional     ::  pack

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Put_Var_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  xtype
    real(4)             ::  packmax
    real(4)             ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! pack ?
    if ( present(pack) .and. pack ) then
      ! cannot pack single slab, needs all ..
      if ( present(start) .or. present(count) ) then
        write (gol,'("could not pack in presence of start or count")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! data type of packed variable:
      status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! maxium value:
      select case ( xtype )
        case ( NF90_SHORT )
          packmax = huge(int(1,kind=2))
        case default
          write (gol,'("unsupported xtype ",i6," for packing data type")') xtype; call goErr
          TRACEBACK; status=1; return
      end select
      ! pacing parameters:
      add_offset   = minval(values)
      scale_factor = ( maxval(values) - add_offset ) / packmax
      ! pack values:
      values = ( values - add_offset )/scale_factor
      ! convert and write:
      select case ( xtype )
        case ( NF90_SHORT )
          status = NF90_Put_Var( self%ncid, varid, int(values,kind=2), start=start, count=count )
          IF_NF90_NOT_OK_RETURN(status=1)
        case default
          write (gol,'("unsupported xtype ",i6," for packing data type")') xtype; call goErr
          TRACEBACK; status=1; return
      end select
      ! write attributes:
      status = NF90_Put_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( self%ncid, varid, 'scale_factor', scale_factor  )
      IF_NF90_NOT_OK_RETURN(status=1)

    else

      ! write:
      status = NF90_Put_Var( self%ncid, varid, values, start=start, count=count )
      IF_NF90_NOT_OK_RETURN(status=1)

    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Put_Var_3d
  
  ! *
  
  subroutine File_NC_Put_Var_4d( self, varid, values, status, &
                                   start, count, pack )

    use NetCDF, only : NF90_SHORT
    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)   ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  values(:,:,:,:)
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    logical, intent(in), optional     ::  pack

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Put_Var_4d'
    
    ! --- local ----------------------------------
    
    integer             ::  xtype
    real(4)             ::  packmax
    real(4)             ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! pack ?
    if ( present(pack) .and. pack ) then

      ! cannot pack single slab, needs all ..
      if ( present(start) .or. present(count) ) then
        write (gol,'("could not pack in presence of start or count")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! data type of packed variable:
      status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! maxium value:
      select case ( xtype )
        case ( NF90_SHORT )
          packmax = huge(int(1,kind=2))
        case default
          write (gol,'("unsupported xtype ",i6," for packing data type")') xtype; call goErr
          TRACEBACK; status=1; return
      end select
      ! pacing parameters:
      add_offset   = minval(values)
      scale_factor = ( maxval(values) - add_offset ) / packmax
      ! pack values:
      values = ( values - add_offset )/scale_factor
      ! convert and write:
      select case ( xtype )
        case ( NF90_SHORT )
          status = NF90_Put_Var( self%ncid, varid, int(values,kind=2), start=start, count=count )
          IF_NF90_NOT_OK_RETURN(status=1)
        case default
          write (gol,'("unsupported xtype ",i6," for packing data type")') xtype; call goErr
          TRACEBACK; status=1; return
      end select
      ! write attributes:
      status = NF90_Put_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( self%ncid, varid, 'scale_factor', scale_factor  )
      IF_NF90_NOT_OK_RETURN(status=1)

    else

      ! write:
      status = NF90_Put_Var( self%ncid, varid, values, start=start, count=count )
      IF_NF90_NOT_OK_RETURN(status=1)

    end if
    
    ! ok
    status = 0
    
  end subroutine File_NC_Put_Var_4d
  
  
  ! ***
  
  
  !
  ! Return record index in time axis
  !
  
  subroutine File_NC_Inq_TimeRecord( self, name, t, irec, status )

    use GO, only : TDate, TIncrDate, operator(+), operator(*), operator(==)
    use GO, only : Extract_Ref_and_Step
    use GO, only : wrtgol
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)   ::  self
    character(len=*), intent(in)      ::  name
    type(TDate), intent(in)           ::  t
    integer, intent(out)              ::  irec
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_NC_Get_VarID'
    
    ! --- local ----------------------------------
    
    integer                 ::  ntime
    integer                 ::  varid
    integer                 ::  itime
    real, allocatable       ::  time(:)
    character(len=64)       ::  units
    type(TDate)             ::  t_ref
    type(TIncrDate)         ::  t_step
    type(TDate)             ::  t_rec
    character(len=8)        ::  label
    
    ! --- begin ----------------------------------
    
    ! number of time records:
    call self%Inquire_Dimension( name, status, length=ntime )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( time(ntime) )
    ! get variable id:
    call self%Get_VarID( 'var_name', name, varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! read values:
    call self%Get_Var( varid, time, status, units=units )
    IF_NOT_OK_RETURN(status=1)
    ! extract ref time and step:
    call Extract_Ref_and_Step( units, 'standard', t_ref, t_step, status )
    IF_NOT_OK_RETURN(status=1)
    ! init index:
    irec = -999
    ! loop:
    do itime = 1, ntime
      ! convert:
      t_rec = t_ref + t_step * time(itime)
      ! compare:
      if ( t_rec == t ) then
        irec = itime
        exit
      end if
    end do
    
    ! check ...
    if ( irec < 0 ) then
      write (gol,'("time record not found:")'); call goErr
      write (gol,'("  filename    : ",a)') trim(self%filename); call goErr
      call wrtgol( '  requested   : ', t ); call goErr
      write (gol,'("  variable    : ",a)') trim(name); call goErr
      write (gol,'("    units     : ",a)') trim(units); call goErr
      call wrtgol( '    ref       : ', t_ref ); call goErr
      call wrtgol( '    step      : ', t_step ); call goErr
      write (gol,'("   records    : ")'); call goErr
      do itime = 1, ntime
        t_rec = t_ref + t_step * time(itime)
        write (label,'("  ",i2.2,"  : ")') itime
        call wrtgol( label, t_rec ); call goErr
      end do
      status = -1; return
    end if
    
    ! clear:
    deallocate( time )
    
    ! ok
    status = 0
    
  end subroutine File_NC_Inq_TimeRecord
  
  
  ! ***


  subroutine File_Nc_Get_Levs( self, levelname, hlevelname, levs, status )
  
    use GO, only : goReadFromLine

    use NetCDF, only : NF90_INQ_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)        ::  self
    character(len=*), intent(in)        ::  levelname
    character(len=*), intent(in)        ::  hlevelname
    XTYPE(T_Levs_Hyb), intent(out)      ::  levs
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Get_Levs'
    
    ! --- local ----------------------------------

    integer                 ::  dimid, varid
    integer                 ::  nlev, nv
    real, allocatable       ::  ap(:), b(:)
    real, allocatable       ::  ap_i(:), b_i(:)
    character(len=256)      ::  formula_terms
    character(len=32)       ::  cname, cvalue
    character(len=32)       ::  varname_ap, varname_b, varname_ps
    
    ! --- begin ----------------------------------
    
    ! dimension with number of levels:
    status = NF90_INQ_DimID( self%ncid, levelname, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! count:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=nlev )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ~ mid levels
    
    ! storage:
    allocate( ap(nlev) )
    allocate( b (nlev) )
    
    ! variable with level indices:
    status = NF90_INQ_VarID( self%ncid, levelname, varid )
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
    
    ! variable with coeffs:
    status = NF90_INQ_VarID( self%ncid, trim(varname_ap), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, ap )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! variable with coeffs:
    status = NF90_INQ_VarID( self%ncid, trim(varname_b), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, b )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ~ interfaces
    
    ! storage:
    allocate( ap_i(0:nlev) )
    allocate( b_i (0:nlev) )
    
    ! variable with level indices:
    status = NF90_INQ_VarID( self%ncid, hlevelname, varid )
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
    
    ! variable with coeffs:
    status = NF90_INQ_VarID( self%ncid, trim(varname_ap), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, ap_i )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! variable with coeffs:
    status = NF90_INQ_VarID( self%ncid, trim(varname_b), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, b_i )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    !~ definition
        
    ! init level definition:
    call levs%Init( ap, b, ap_i, b_i, status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( ap )
    deallocate( b )
    deallocate( ap_i )
    deallocate( b_i )
    
    ! ok
    status = 0
    
  end subroutine File_Nc_Get_Levs


  ! ***
  
  
  subroutine File_Nc_Def_Levs( self, hyb, levelname, hlevelname, status )

    use NetCDF, only : NF90_SHORT, NF90_INT, NF90_FLOAT
    use NetCDF, only : NF90_Def_Dim
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att

    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)     ::  self
    XTYPE(T_Levs_Hyb), intent(in)       ::  hyb
    character(len=*), intent(in)        ::  levelname
    character(len=*), intent(in)        ::  hlevelname
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Def_Levs'
    
    ! --- local ----------------------------------
    
    integer               ::  varid

    ! --- begin ----------------------------------

    ! ~ midpoints

    ! define dimension:
    status = NF90_Def_Dim( self%ncid, levelname, hyb%nlev, self%dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! index variable:
    status = NF90_Def_Var( self%ncid, levelname, NF90_INT, &
                              (/self%dimid_lev/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid level at layer midpoints' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'positive', 'down' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'formula_terms', 'ap: hyam b: hybm ps: sp' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_lev = varid
        
    ! ap values:
    status = NF90_Def_Var( self%ncid, 'hyam', NF90_FLOAT, &
                              (/self%dimid_lev/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid A coefficient at layer midpoints' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_hya = varid
    
    ! b values:
    status = NF90_Def_Var( self%ncid, 'hybm', NF90_FLOAT, &
                              (/self%dimid_lev/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid B coefficient at layer midpoints' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_hyb = varid

    ! ~ interfaces

    ! define dimension:
    status = NF90_Def_Dim( self%ncid, hlevelname, hyb%nlev+1, self%dimid_levi )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! index variable:
    status = NF90_Def_Var( self%ncid, hlevelname, NF90_INT, &
                              (/self%dimid_levi/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid level at layer interfaces' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'positive', 'down' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'formula_terms', 'ap: hyai b: hybi ps: sp' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_levi = varid
        
    ! ap values:
    status = NF90_Def_Var( self%ncid, 'hyai', NF90_FLOAT, &
                              (/self%dimid_levi/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid A coefficient at layer interfaces' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', 'Pa' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_hyai = varid
    
    ! b values:
    status = NF90_Def_Var( self%ncid, 'hybi', NF90_FLOAT, &
                              (/self%dimid_levi/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'hybrid B coefficient at layer interfaces' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !
    self%varid_hybi = varid
        
    
    ! ok
    status = 0
    
  end subroutine File_Nc_Def_Levs


  ! ***
  
  
  subroutine File_Nc_Put_Levs( self, hyb, status )

    use NetCDF, only : NF90_Put_Var

    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(inout)     ::  self
    XTYPE(T_Levs_Hyb), intent(in)       ::  hyb
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Put_Levs'
    
    ! --- local ----------------------------------
    
    integer, allocatable   ::  lev(:)
    integer, allocatable   ::  levi(:)
    integer                ::  k
    
    ! --- begin ----------------------------------
    
    ! storage:
    allocate( lev(hyb%nlev) )
    ! fill:
    do k = 1, hyb%nlev
      lev(k) = k
    end do
    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_lev, lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( lev )

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_hya, hyb%ap )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_hyb, hyb%b )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( levi(0:hyb%nlev) )
    ! fill:
    do k = 0, hyb%nlev
      levi(k) = k
    end do
    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_levi, levi )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( levi )

    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_hyai, hyb%ap_i )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! write:
    status = NF90_Put_Var( self%ncid, self%varid_hybi, hyb%b_i )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Nc_Put_Levs


  ! ***
  
  
  subroutine File_Nc_Copy_Variable_CoarsenLevs( self, varname, outfile, &
                                      levelname, hlevelname, &
                                      levs_in, levs_out, mapping, status )

    use NetCDF, only : NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
    use NetCDF, only : NF90_Inq_AttName, NF90_Inquire_Attribute, NF90_Get_Att

    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Nc), intent(in)      ::  self
    character(len=*), intent(in)      ::  varname
    XTYPE(T_File_Nc), intent(inout)   ::  outfile
    character(len=*), intent(in)      ::  levelname
    character(len=*), intent(in)      ::  hlevelname
    XTYPE(T_Levs_Hyb), intent(in)     ::  levs_in
    XTYPE(T_Levs_Hyb), intent(in)     ::  levs_out
    integer, intent(in)               ::  mapping(:,:)  ! (nlev,2) per level range
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Nc_Copy_Variable_CoarsenLevs'
    
    ! --- local ----------------------------------
    
    integer                ::  varid, varid_out, varid_sp
    integer                ::  xtype
    integer                ::  ndim, idim
    integer                ::  dimids(maxdim)
    integer                ::  shp(maxdim)
    logical                ::  levdim(maxdim)
    character(len=64)      ::  dimname
    integer                ::  ilev 
    integer                ::  k
    integer                ::  irec
    
    logical                ::  packed
    integer                ::  xtype_packed

    real, allocatable      ::  data_1d(:)

    real, allocatable      ::  data_2d(:,:)

    real, allocatable      ::  data_3d(:,:,:), data_3d_out(:,:,:)
    real, allocatable      ::    sp_3d1(:,:)
    real, allocatable      ::     w_3d1(:,:), ws_3d1(:,:)
    real, allocatable      ::    dp_3d(:,:,:)

    real, allocatable      ::  data_4d_out(:,:,:,:)
    
    ! --- begin ----------------------------------
    
    ! id of input variable:
    status = NF90_INQ_VarID( self%ncid, varname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! obtain info on dimensions:
    status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndim, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain shape:
    do idim = 1, ndim
      ! get length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! init flags for level dimension:
    levdim = .false.
    do idim = 1, ndim
      ! get name:
      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dimname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! compare:
      levdim(idim) = trim(dimname) == trim(levelname)
    end do
    ! check ...
    if ( count(levdim) > 1 ) then
      write (gol,'("could not coarsen for more ",i2," level dimensions")') count(levdim); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! obtain info on data type, might be packed:
    call self%Inquire_Variable( varid, status, xtype=xtype_packed, dtype=xtype )
    IF_NOT_OK_RETURN(status=1)
    ! packed?
    packed = xtype_packed /= xtype
    
    ! id of output variable:
    status = NF90_INQ_VarID( outfile%ncid, varname, varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! switch:
    !............................................  
    if ( any(levdim) ) then
    !............................................  
    
      ! switch on rank:
      select case ( ndim )

        !~ 3D
        case ( 3 )

          ! storage for result:
          allocate( data_3d_out(shp(1),levs_out%nlev,shp(3)), stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! which ?
          if ( levdim(2) ) then

            ! storage:
            allocate( data_3d(shp(1),shp(2),shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! storage for surface pressure and delta-pressure:
            allocate( sp_3d1(shp(1),       shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( dp_3d (shp(1),shp(2),shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! storage for weights:
            allocate(  w_3d1(shp(1),       shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( ws_3d1(shp(1),       shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)

            ! read, unpack immediatelly if necessary:
            call self%Get_Var( varid, data_3d, status, unpack=packed )
            IF_NOT_OK_RETURN(status=1)            

            ! id of surface pressure:
            status = NF90_INQ_VarID( self%ncid, 'sp', varid_sp )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! read surface pressure
            call self%Get_Var( varid_sp, sp_3d1, status, unpack=packed )
            IF_NOT_OK_RETURN(status=1)
            ! fill delta-pressure:
            call levs_in%Ps_to_dP( sp_3d1, dp_3d, status )
            IF_NOT_OK_RETURN(status=1)
            ! loop over target levels:
            do ilev = 1, levs_out%nlev
              ! init sum:
              data_3d_out(:,ilev,:) = 0.0
              ! init weight sum:
              ws_3d1 = 0.0
              ! loop over source levels:
              do k = mapping(ilev,1), mapping(ilev,2)
                ! weight:
                w_3d1 = dp_3d(:,k,:)
                ! add contribution:
                data_3d_out(:,ilev,:) = data_3d_out(:,ilev,:) + data_3d(:,k,:) * w_3d1
                ! update sum:
                ws_3d1 = ws_3d1 + w_3d1
              end do ! source levels
              ! average:
              data_3d_out(:,ilev,:) = data_3d_out(:,ilev,:) / ws_3d1
            end do  ! target levels

            ! clear:
            deallocate( data_3d, stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( sp_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)
            deallocate( dp_3d , stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate(  w_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)
            deallocate( ws_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)

          else
            ! not yet ...
            write (gol,'("unsupported level dimensions : ",3l2)') levdim(1:ndim); call goErr
            TRACEBACK; status=1; return
          end if

          ! write, eventually packed:
          call outfile%Put_Var( varid_out, data_3d_out, status, pack=packed )
          IF_NF90_NOT_OK_RETURN(status=1)
          
          ! clear:
          deallocate( data_3d_out, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~ 4D
        case ( 4 )

          ! storage for result:
          allocate( data_4d_out(shp(1),shp(2),levs_out%nlev,shp(4)), stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! which ?
          if ( levdim(3) ) then
          
            ! storage:
            allocate( data_3d(shp(1),shp(2),shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)

            ! storage for surface pressure and delta-pressure:
            allocate( sp_3d1(shp(1),shp(2)       ), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( dp_3d (shp(1),shp(2),shp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)

            ! storage for weights:
            allocate(  w_3d1(shp(1),shp(2)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( ws_3d1(shp(1),shp(2)), stat=status )
            IF_NOT_OK_RETURN(status=1)

            ! loop over records:
            do irec = 1, shp(4)

              ! read, unpack immediatelly if necessary:
              call self%Get_Var( varid, data_3d, status, unpack=packed, &
                                    start=(/     1,     1,     1,irec/), &
                                    count=(/shp(1),shp(2),shp(3),   1/) )
              IF_NOT_OK_RETURN(status=1)            

              ! id of surface pressure:
              status = NF90_INQ_VarID( self%ncid, 'sp', varid_sp )
              IF_NF90_NOT_OK_RETURN(status=1)
              ! read surface pressure
              call self%Get_Var( varid_sp, sp_3d1, status, unpack=packed, &
                                    start=(/     1,     1,irec/), &
                                    count=(/shp(1),shp(2),   1/)    )
              IF_NOT_OK_RETURN(status=1)
              ! fill delta-pressure:
              call levs_in%Ps_to_dP( sp_3d1, dp_3d, status )
              IF_NOT_OK_RETURN(status=1)
              ! loop over target levels:
              do ilev = 1, levs_out%nlev
                ! init sums:
                data_4d_out(:,:,ilev,irec) = 0.0
                ! init weight sum:
                ws_3d1 = 0.0
                ! loop over source levels:
                do k = mapping(ilev,1), mapping(ilev,2)
                  ! weight:
                  w_3d1 = dp_3d(:,:,k)
                  ! add contribution:
                  data_4d_out(:,:,ilev,irec) = data_4d_out(:,:,ilev,irec) + data_3d(:,:,k) * w_3d1
                  ! update sum:
                  ws_3d1 = ws_3d1 + w_3d1 
                end do ! source levels
                ! average:
                data_4d_out(:,:,ilev,irec) = data_4d_out(:,:,ilev,irec) / ws_3d1
              end do  ! target levels

            end do ! records

            ! clear weights:
            deallocate(  w_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)
            deallocate( ws_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)

            ! clear pressures:
            deallocate( sp_3d1, stat=status )
            IF_NOT_OK_RETURN(status=1)
            deallocate( dp_3d , stat=status )
            IF_NOT_OK_RETURN(status=1)
            
            ! clear:
            deallocate( data_3d, stat=status )
            IF_NOT_OK_RETURN(status=1)

          else
            ! not yet ...
            write (gol,'("unsupported level dimensions : ",4l2)') levdim(1:ndim); call goErr
            TRACEBACK; status=1; return
          end if
          
          ! write all together, eventually packed:
          call outfile%Put_Var( varid_out, data_4d_out, status, pack=packed )
          IF_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data_4d_out, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~
        case default
          write (gol,'("unsupported ndim ",i6)') ndim; call goErr
          TRACEBACK; status=1; return
      end select

    !............................................  
    else  ! no level dimensions ...
    !............................................  

      ! switch on rank:
      select case ( ndim )

        !~ 1D
        case ( 1 )

          ! storage:
          allocate( data_1d(shp(1)), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! read:
          status = NF90_Get_Var( self%ncid, varid, data_1d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! write:
          status = NF90_Put_Var( outfile%ncid, varid_out, data_1d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! clear:
          deallocate( data_1d, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~ 2D
        case ( 2 )

          ! storage:
          allocate( data_2d(shp(1),shp(2)), stat=status )
          ! read:
          IF_NOT_OK_RETURN(status=1)
          status = NF90_Get_Var( self%ncid, varid, data_2d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! write:
          status = NF90_Put_Var( outfile%ncid, varid_out, data_2d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! clear:
          deallocate( data_2d, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~ 3D
        case ( 3 )

          ! storage:
          allocate( data_3d(shp(1),shp(2),shp(3)), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! read:
          status = NF90_Get_Var( self%ncid, varid, data_3d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! write:
          status = NF90_Put_Var( outfile%ncid, varid_out, data_3d )
          IF_NF90_NOT_OK_RETURN(status=1)
          ! clear:
          deallocate( data_3d, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~ 4D
        case ( 4 )

          ! storage:
          allocate( data_3d(shp(1),shp(2),shp(3)), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop:
          do irec = 1, shp(4)
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_3d, &
                                    start=(/1,1,1,irec/), &
                                    count=(/shp(1),shp(2),shp(3),1/) )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_3d, &
                                    start=(/1,1,1,irec/), &
                                    count=(/shp(1),shp(2),shp(3),1/) )
            IF_NF90_NOT_OK_RETURN(status=1)
          end do ! records
          ! clear:
          deallocate( data_3d, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~
        case default
          write (gol,'("unsupported ndim ",i6)') ndim; call goErr
          TRACEBACK; status=1; return
      end select

    !............................................  
    end if  ! level dims?  
    !............................................    

    ! ok
    status = 0
    
  end subroutine File_Nc_Copy_Variable_CoarsenLevs
  
  

end module C3PO_File_Nc
