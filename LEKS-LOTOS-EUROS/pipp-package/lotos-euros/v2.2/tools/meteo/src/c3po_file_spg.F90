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
!######################################################################

module C3PO_File_Spg

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError
  
  use C3PO_File_Nc, only : T_File_Nc, maxdim

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  dimname_c2, dimname_sp
  
  public  ::  T_File_Spg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Spg'

  ! sh parameters:    
  character(len=*), parameter  ::  dimname_c2   = 'nc2'
  character(len=*), parameter  ::  dimname_sp   = 'nsp'
!  character(len=*), parameter  ::  dimname_rgrid = 'rgrid'
!  character(len=*), parameter  ::  dimname_nv    = 'nv'
!  character(len=*), parameter  ::  varname_lon   = 'longitude'
!  character(len=*), parameter  ::  varname_lat   = 'latitude'
!  character(len=*), parameter  ::  varext_bnds   = '_bnds'
  integer, parameter           ::  nc2 = 2
  

  ! --- types ----------------------------------------
  
#ifdef without_f2003
  type T_File_Spg
    character(len=1024)     ::  filename
    integer                 ::  ncid
    integer                 ::  formatNum
#else
  type, extends(T_File_Nc) :: T_File_Spg
#endif
    ! spg dims:
    integer                 ::  dimid_c2
    integer                 ::  dimid_sp
  contains
#ifdef without_f2003
    procedure   ::  Open                    => File_Spg_Open
    procedure   ::  Create                  => File_Spg_Create
    procedure   ::  EndDef                  => File_Spg_EndDef
    procedure   ::  Close                   => File_Spg_Close
    procedure   ::  Inquire                 => File_Spg_Inquire
    procedure   ::  Inquire_Dimension       => File_Spg_Inquire_Dimension
    procedure   ::  Copy_Dimension          => File_Spg_Copy_Dimension
    procedure   ::  Inquire_Variable        => File_Spg_Inquire_Variable
    procedure   ::  Define_Variable_Copy    => File_Spg_Define_Variable_Copy
    procedure   ::  Copy_Atts               => File_Spg_Copy_Atts
    procedure   ::  Extend_History          => File_Spg_Extend_History
#endif
    procedure   ::  Extend_Variable_Att       => File_Spg_Extend_Variable_Att
    procedure   ::  Inquire_Variable_Spectral => File_Spg_Inquire_Variable_Spectral
    procedure   ::  Rename_Variable           => File_Spg_Rename_Variable
    procedure   ::  Get_Grid                  => File_Spg_Get_Grid
    procedure   ::  Def_Grid                  => File_Spg_Def_Grid
    procedure   ::  Put_Grid                  => File_Spg_Put_Grid
    procedure   ::  Copy_Variable             => File_Spg_Copy_Variable
    procedure   ::  Put_Variable_Gradients    => File_Spg_Put_Variable_Gradients
    !procedure   ::  Define_Variable_Copy_Rgg  => File_Spg_Define_Variable_Copy_Rgg
  end type T_File_Spg

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif
    
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** rgg file
  ! ***
  ! ********************************************************************


#ifdef without_f2003

  subroutine File_Spg_Open( self, filename, status )
  
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_NOWRITE
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(out)   ::  self
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Open'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! open file for reading:
    status = NF90_Open( trim(filename), NF90_NOWRITE, self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! info ...
    status = NF90_Inquire( self%ncid, formatNum=self%formatNum )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine File_Spg_Open


  ! ***


  subroutine File_Spg_Create( self, filename, status, formatNum )
  
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_CLOBBER
    use NetCDF, only : NF90_CLASSIC_MODEL
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(out)   ::  self
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status
    integer, intent(in), optional   ::  formatNum

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Create'
    
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
  
  end subroutine File_Spg_Create


  ! ***
  
  
  subroutine File_Spg_EndDef( self, status )

    use NetCDF, only : NF90_EndDef
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout) ::  self
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_EndDef( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_Spg_EndDef


  ! ***
  
  
  subroutine File_Spg_Close( self, status )

    use NetCDF, only : NF90_Close
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout)   ::  self
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Close'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_Close( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Close


  ! ***
  
  
  subroutine File_Spg_Inquire( self, status, ndim, nvar )

    use NetCDF, only : NF90_Inquire
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)    ::  self
    integer, intent(out)            ::  status
    integer, intent(out), optional  ::  ndim
    integer, intent(out), optional  ::  nvar

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Inquire'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! obtain info:
    status = NF90_Inquire( self%ncid, nDimensions=ndim, nVariables=nvar )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Inquire


  ! ***
  
  
  subroutine File_Spg_Inquire_Dimension( self, idim, status, name )

    use NetCDF, only : NF90_Inquire_Dimension
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)              ::  self
    integer, intent(in)                       ::  idim
    integer, intent(out)                      ::  status
    character(len=*), intent(out), optional   ::  name

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Inquire_Dimension'
    
    ! --- local ----------------------------------
    
    integer     ::  dimid
    
    ! --- begin ----------------------------------
    
    ! variable number is variable id ?
    dimid = idim
    ! obtain info:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=name )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Inquire_Dimension


  ! ***
  
  
  subroutine File_Spg_Copy_Dimension( self, dimname, outfile, status )

    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Def_Dim
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)              ::  self
    character(len=*), intent(in)              ::  dimname
    XTYPE(T_File_Spg), intent(out)             ::  outfile
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Copy_Dimension'
    
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
    
  end subroutine File_Spg_Copy_Dimension


  ! ***
  
  
  subroutine File_Spg_Inquire_Variable( self, ivar, status, name )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    integer, intent(in)                       ::  ivar
    integer, intent(out)                      ::  status
    character(len=*), intent(out), optional   ::  name

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Inquire_Variable'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    
    ! --- begin ----------------------------------
    
    ! variable number is variable id:
    varid = ivar
    
    ! obtain info:
    status = NF90_Inquire_Variable( self%ncid, varid, name=name )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Inquire_Variable


  ! ***
  

  subroutine File_Spg_Define_Variable_Copy( self, varname, outfile, varid, status )

    use NetCDF, only : NF90_Inquire_Dimension, NF90_Inq_DimID
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Def_Var

    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    character(len=*), intent(in)              ::  varname
    XTYPE(T_File_Spg), intent(out)            ::  outfile
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Define_Variable_Copy'
    
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
    
  end subroutine File_Spg_Define_Variable_Copy


  ! ***
  
  
  !
  ! Copy attributes from one file to another ;
  ! 'varname' should be the variable name,
  ! or empty or '/' for global attributes.
  !
  
  subroutine File_Spg_Copy_Atts( self, varname, outfile, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName, NF90_Copy_Att
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)     ::  self
    character(len=*), intent(in)      ::  varname
    XTYPE(T_File_Spg), intent(out)    ::  outfile
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Copy_Atts'
    
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
    
  end subroutine File_Spg_Copy_Atts
  
  
  ! ***
  
  
  !
  ! Prepend line to history:
  !   Wed Apr 23 16:32:03 2014: Message\n
  !
  ! Only allowed in definition mode.
  !
  
  subroutine File_Spg_Extend_History( self, message, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Get_Att, NF90_Put_Att
    
    use GO, only : TDate, SystemDate
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout)   ::  self
    character(len=*), intent(in)      ::  message
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Extend_History'
    
    ! month names:
    character(len=3),parameter :: monthname(12) = &
                         (/'Jan','Feb','Mrc','Apr','May','Jun',&
                           'Jul','Aug','Sep','Oct','Nov','Dec'/)
   
    ! newline character:
    character(len=1), parameter  ::  newline = char(10)
     
    ! --- local ----------------------------------
    
    character(len=64)     ::  attname
    character(len=64)     ::  tstamp
    character(len=1024)   ::  history, history_curr
    type(TDate)           ::  t
    
    ! --- begin ----------------------------------
    
    ! target attribute:
    attname = 'history'
    
    ! current:
    status = NF90_Get_Att( self%ncid, NF90_GLOBAL, attname, history_curr )
    if ( status /= NF90_NOERR ) history_curr = ''
    
    ! current time:
    t = SystemDate()
    ! write time stamp:
    write (tstamp,'("Day ",a," ",i2," ",i2,2(":",i2.2)," ",i4)') &
             monthname(t%month), t%day, t%hour, t%min, t%sec, t%year
    ! fill history:
    if ( len_trim(history_curr)== 0 ) then
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message)
    else
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message), &
                                               newline, trim(history_curr)
    end if
    
    ! write:
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, attname, trim(history) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Extend_History
  
  
#endif

  ! ***
  
  
  !
  ! Add text to existing variable attribute.
  !
  
  subroutine File_Spg_Extend_Variable_Att( self, varname, prefix, attname, postfix, status )

    use NetCDF, only : NF90_INQ_VarID
    use NetCDF, only : NF90_Get_Att, NF90_Put_Att
    
    use GO, only : TDate, SystemDate
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout)    ::  self
    character(len=*), intent(in)        ::  varname
    character(len=*), intent(in)        ::  prefix
    character(len=*), intent(in)        ::  attname
    character(len=*), intent(in)        ::  postfix
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Extend_Variable_Att'

    ! --- local ----------------------------------
    
    integer                 ::  varid
    character(len=1024)     ::  value
    character(len=1024)     ::  newvalue
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    status = NF90_INQ_VarID( self%ncid, trim(varname), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! current:
    status = NF90_Get_Att( self%ncid, varid, trim(attname), value )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! extend:
    write (newvalue,'(3a)') prefix, trim(value), postfix
    
    ! write:
    status = NF90_Put_Att( self%ncid, varid, attname, trim(newvalue) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Extend_Variable_Att
  
  
  ! ***
  
  
  subroutine File_Spg_Inquire_Variable_Spectral( self, ivar, status, spectral )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    integer, intent(in)                       ::  ivar
    integer, intent(out)                      ::  status
    logical, intent(out), optional            ::  spectral

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Inquire_Variable_Spectral'
    
    ! --- local ----------------------------------
    
    integer             ::  ndim, idim
    integer             ::  dimids(maxdim)
    character(len=64)   ::  dimname
    integer             ::  varid
    logical             ::  has_c2, has_sp
    
    ! --- begin ----------------------------------
    
    ! variable number is variable id:
    varid = ivar
    
    ! spectral field ?
    if ( present(spectral) ) then
      ! init flags:
      has_c2 = .false.
      has_sp = .false.
      ! number of dims:
      status = NF90_Inquire_Variable( self%ncid, varid, &
                                       ndims=ndim, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! loop:
      do idim = 1, ndim
        ! get name:
        status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dimname )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! compare:
        if ( trim(dimname) == trim(dimname_c2) ) has_c2 = .true.
        if ( trim(dimname) == trim(dimname_sp) ) has_sp = .true.
      end do 
      ! combine:
      spectral = has_c2 .and. has_sp
    end if
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Inquire_Variable_Spectral


  ! ***

  
  
  subroutine File_Spg_Rename_Variable( self, varname, newname, status )

    use NetCDF, only : NF90_INQ_VarID, NF90_Rename_Var
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  newname
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Rename_Variable'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    
    ! --- begin ----------------------------------
    
    ! id of input variable:
    status = NF90_INQ_VarID( self%ncid, varname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! rename:
    status = NF90_Rename_Var( self%ncid, varid, newname )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Rename_Variable


  ! ***
  
  
  subroutine File_Spg_Get_Grid( self, spg, status )

    use NetCDF, only : NF90_INQ_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Get_Var
  
    use C3PO_Grid, only : T_Grid_Spg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)      ::  self
    type(T_Grid_Spg), intent(out)      ::  spg
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Get_Grid'
    
    ! --- local ----------------------------------

    integer                 ::  dimid!, varid
    integer                 ::  n
!    real, allocatable       ::  band_lats(:)
!    real, allocatable       ::  band_lats_bnds(:,:)
!    integer, allocatable    ::  band_nlon(:)
!    real, allocatable       ::  lons(:)
!    real, allocatable       ::  lons_bnds(:,:)
    
    ! --- begin ----------------------------------
    
    ! dimension with real/imag :
    status = NF90_INQ_DimID( self%ncid, dimname_c2, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! count:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=n )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! check ...
    if ( n /= nc2 ) then
      write (gol,'("found ",a," dimension with length ",i6,", found ",i6)') trim(dimname_c2), nc2, n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! dimension with spectral coeff:
    status = NF90_INQ_DimID( self%ncid, dimname_sp, dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! count:
    status = NF90_Inquire_Dimension( self%ncid, dimid, len=n )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! init grid defintion:
    call spg%Init_np( n, status )

    ! ok
    status = 0
    
  end subroutine File_Spg_Get_Grid


  ! ***
  
  
  subroutine File_Spg_Def_Grid( self, spg, status )

    use NetCDF, only : NF90_Def_Dim

    use C3PO_Grid, only : T_Grid_Spg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout)    ::  self
    XTYPE(T_Grid_Spg), intent(in)       ::  spg
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Def_Grid'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! define dimension:
    status = NF90_Def_Dim( self%ncid, dimname_c2, nc2, self%dimid_c2 )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! define dimension:
    status = NF90_Def_Dim( self%ncid, dimname_sp, spg%np, self%dimid_sp )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine File_Spg_Def_Grid


  ! ***
  
  
  subroutine File_Spg_Put_Grid( self, spg, status )

    use C3PO_Grid, only : T_Grid_Spg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(inout)    ::  self
    type(T_Grid_Spg), intent(in)        ::  spg
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Put_Grid'
    
    ! --- local ----------------------------------
    
    integer   ::  id
    
    ! --- begin ----------------------------------
    
    ! no coordinate variables yet;
    ! do something to avoid warnings about unused arguments:
    id = self%ncid
    id = spg%T

    ! ok
    status = 0
    
  end subroutine File_Spg_Put_Grid


  ! ***
  
  
  subroutine File_Spg_Copy_Variable( self, varname, outfile, status )

    use NetCDF, only : NF90_INT, NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)     ::  self
    character(len=*), intent(in)      ::  varname
    XTYPE(T_File_Spg), intent(inout)  ::  outfile
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Copy_Variable'
    
    ! --- local ----------------------------------
    
    integer             ::  varid, varid_out
    integer             ::  xtype
    integer             ::  ndim, idim
    integer             ::  dimids(maxdim)
    integer             ::  shp(maxdim)
    
    integer(4), allocatable   ::  data_i4_1d(:)
    real(4), allocatable      ::  data_r4_1d(:)
    real(8), allocatable      ::  data_r8_1d(:)
    
    integer(4), allocatable   ::  data_i4_2d(:,:)
    real(4), allocatable      ::  data_r4_2d(:,:)
    real(8), allocatable      ::  data_r8_2d(:,:)
    
    integer(4), allocatable   ::  data_i4_3d(:,:,:)
    real(4), allocatable      ::  data_r4_3d(:,:,:)
    real(8), allocatable      ::  data_r8_3d(:,:,:)

    ! --- begin ----------------------------------
    
    ! id of input variable:
    status = NF90_INQ_VarID( self%ncid, varname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! obtain info:
    status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype, &
                                       ndims=ndim, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain shape:
    do idim = 1, ndim
      ! get length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! id of output variable:
    status = NF90_INQ_VarID( outfile%ncid, varname, varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! switch on rank:
    select case ( ndim )

      !~ 1D
      case ( 1 )
        ! create storage and read:
        select case ( xtype )
        
          !~ i4
          case ( NF90_INT )
            ! storage:
            allocate( data_i4_1d(shp(1)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_i4_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_i4_1d )
        
          !~ r4
          case ( NF90_FLOAT )
            ! storage:
            allocate( data_r4_1d(shp(1)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r4_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r4_1d )
        
          !~ r8
          case ( NF90_DOUBLE )
            ! storage:
            allocate( data_r8_1d(shp(1)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r8_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_1d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r8_1d )

          !~ unkown
          case default
            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
            TRACEBACK; status=1; return
        end select

      !~ 2D
      case ( 2 )
        ! create storage and read:
        select case ( xtype )
        
          !~ i4
          case ( NF90_INT )
            ! storage:
            allocate( data_i4_2d(shp(1),shp(2)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_i4_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_i4_2d )
        
          !~ r4
          case ( NF90_FLOAT )
            ! storage:
            allocate( data_r4_2d(shp(1),shp(2)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r4_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r4_2d )
        
          !~ r8
          case ( NF90_DOUBLE )
            ! storage:
            allocate( data_r8_2d(shp(1),shp(2)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r8_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_2d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r8_2d )

          !~ unkown
          case default
            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
            TRACEBACK; status=1; return
        end select

      !~ 3d
      case ( 3 )
        ! create storage and read:
        select case ( xtype )
        
          !~ i4
          case ( NF90_INT )
            ! storage:
            allocate( data_i4_3d(shp(1),shp(2),shp(3)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_i4_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_i4_3d )
        
          !~ r4
          case ( NF90_FLOAT )
            ! storage:
            allocate( data_r4_3d(shp(1),shp(2),shp(3)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r4_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r4_3d )
        
          !~ r8
          case ( NF90_DOUBLE )
            ! storage:
            allocate( data_r8_3d(shp(1),shp(2),shp(3)) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r8_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! write:
            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r8_3d )

          !~ unkown
          case default
            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
            TRACEBACK; status=1; return
        end select
        
      !~
      case default
        write (gol,'("unsupported ndim ",i6)') ndim; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Copy_Variable


  ! ***
  
  
  subroutine File_Spg_Put_Variable_Gradients( self, spg, varname, &
                                               outfile, dvarnames, status )

    use NetCDF, only : NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
    
    use C3PO_Grid, only : T_Grid_Spg
  
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Spg_Put_Variable_Gradients'
    
    ! number of directions
    integer, parameter      ::  ndir = 2
    
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)       ::  self
    XTYPE(T_Grid_Spg), intent(in)       ::  spg
    character(len=*), intent(in)        ::  varname
    XTYPE(T_File_Spg), intent(inout)    ::  outfile
    character(len=*), intent(in)        ::  dvarnames(ndir)
    integer, intent(out)                ::  status

    ! --- local ----------------------------------
    
    integer             ::  varid, dvarids(ndir)
    integer             ::  idir
    integer             ::  xtype
    integer             ::  ndim, idim
    integer             ::  dimids(maxdim)
    integer             ::  shp(maxdim)
    character(len=64)   ::  dimnames(maxdim)
    integer             ::  k

    ! real data in default working precision:
    real, allocatable   ::  data_r_3d(:,:,:), ddata_r_3d(:,:,:,:)
    
    ! --- begin ----------------------------------
    
    ! id of input variable:
    status = NF90_INQ_VarID( self%ncid, varname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! obtain info:
    status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype, &
                                       ndims=ndim, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain shape:
    do idim = 1, ndim
      ! get length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), &
                                         len=shp(idim), name=dimnames(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! id of output variable:
    do idir = 1, ndir
      status = NF90_INQ_VarID( outfile%ncid, dvarnames(idir), dvarids(idir) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! switch on rank:
    select case ( ndim )

      !~ 3D
      case ( 3 )
        ! create storage and read:
        select case ( xtype )
        
          !~ r4
          case ( NF90_FLOAT, NF90_DOUBLE )
            ! check ...
            if ( (trim(dimnames(1)) /= dimname_c2) .or. (trim(dimnames(2)) /= dimname_sp) ) then
              write (gol,'("unexpected dim names for spectral field ",a)') trim(varname); call goErr
              do idim = 1, ndim
                write (gol,'("  ",i4," ",a)') idim, trim(dimnames(idim)); call goPr
              end do
              TRACEBACK; status=1; return
            end if
            ! storage:
            allocate(  data_r_3d(shp(1),shp(2),shp(3)) )
            allocate( ddata_r_3d(shp(1),shp(2),shp(3),ndir) )
            ! read:
            status = NF90_Get_Var( self%ncid, varid, data_r_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop:
            do k = 1, shp(3)
              ! compute spectral coeff of gradients:
              call spg%Nabla( data_r_3d(:,:,k), ddata_r_3d(:,:,k,:), status )
              IF_NOT_OK_RETURN(status=1)
            end do
            ! loop over directions:
            do idir = 1, ndir
              ! write:
              status = NF90_Put_Var( outfile%ncid, dvarids(idir), ddata_r_3d(:,:,:,idir) )
              IF_NF90_NOT_OK_RETURN(status=1)
            end do
            ! clear:
            deallocate( data_r_3d )
            deallocate( ddata_r_3d )

          !~ unkown
          case default
            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
            TRACEBACK; status=1; return
        end select
        
      !~
      case default
        write (gol,'("unsupported ndim ",i6)') ndim; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine File_Spg_Put_Variable_Gradients


  ! ***
  
  
end module C3PO_File_Spg
