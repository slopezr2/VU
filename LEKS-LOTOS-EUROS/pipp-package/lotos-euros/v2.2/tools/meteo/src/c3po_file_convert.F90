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

module C3PO_File_Convert

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  Spg_Copy_Dimension_Rgg
  public  ::  Spg_Define_Variable_Copy_Rgg, Spg_Define_Variable_RegriddedCopy_Rgg
  public  ::  Spg_Copy_Variable_Rgg, Spg_Regrid_Variable_Rgg
  public  ::  Spg_Copy_Atts_Rgg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Convert'

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif

  
contains


  ! ********************************************************************
  ! ***
  ! *** spectral field evaluated to reduced-gaussian-grid
  ! ***
  ! ********************************************************************


  subroutine Spg_Copy_Dimension_Rgg( self, dimname, outfile, status )

    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Def_Dim
    
    use C3PO_File_Spg, only : T_File_Spg
    use C3PO_File_Rgg, only : T_File_Rgg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    character(len=*), intent(in)              ::  dimname
    XTYPE(T_File_Rgg), intent(inout)          ::  outfile
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Copy_Dimension_Rgg'
    
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
    
  end subroutine Spg_Copy_Dimension_Rgg
  
  
  !
  ! Copy attributes from one file to another ;
  ! 'varname' should be the variable name,
  ! or empty or '/' for global attributes.
  !
  
  subroutine Spg_Copy_Atts_Rgg( self, varname, outfile, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName, NF90_Copy_Att
    
    use C3PO_File_Spg, only : T_File_Spg
    use C3PO_File_Rgg, only : T_File_Rgg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)     ::  self
    character(len=*), intent(in)      ::  varname
    XTYPE(T_File_Rgg), intent(inout)  ::  outfile
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Copy_Atts_Rgg'
    
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
    
  end subroutine Spg_Copy_Atts_Rgg


  ! ***
  

  subroutine Spg_Define_Variable_Copy_Rgg( self, varname, outfile, varid, status )

    use NetCDF, only : NF90_Inquire_Dimension, NF90_Inq_DimID
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Def_Var
    
    use C3PO_File_Nc , only : maxdim
    use C3PO_File_Rgg, only : T_File_Rgg
    use C3PO_File_Spg, only : T_File_Spg

    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    character(len=*), intent(in)              ::  varname
    XTYPE(T_File_Rgg), intent(out)            ::  outfile
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Define_Variable_Copy_Rgg'
    
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
    call Spg_Copy_Atts_Rgg( self, varname, outfile, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Spg_Define_Variable_Copy_Rgg


  ! ***
  

  subroutine Spg_Define_Variable_RegriddedCopy_Rgg( self, varname, rgridname, outfile, &
                                                  varid, status )

    use NetCDF, only : NF90_Inquire_Dimension, NF90_Inq_DimID
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Def_Var
    use NetCDF, only : NF90_Inq_AttName, NF90_Put_Att, NF90_Copy_Att
    
    use C3PO_File_Nc , only : maxdim
    use C3PO_File_Rgg, only : T_File_Rgg
    use C3PO_File_Spg, only : T_File_Spg, dimname_c2, dimname_sp

    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)             ::  self
    character(len=*), intent(in)              ::  varname
    character(len=*), intent(in)              ::  rgridname
    XTYPE(T_File_Rgg), intent(inout)          ::  outfile
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Define_Variable_RegriddedCopy_Rgg'
    
    ! --- local ----------------------------------
    
    integer             ::  varid_in
    integer             ::  xtype
    integer             ::  ndim, idim
    integer             ::  dimids_in(maxdim)
    integer             ::  dimids_out(maxdim)
    character(len=64)   ::  dimname
    integer             ::  natt, iatt
    character(len=64)   ::  attname
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    status = NF90_INQ_VarID( self%ncid, varname, varid_in )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain info:
    status = NF90_Inquire_Variable( self%ncid, varid_in, xtype=xtype, &
                                           ndims=ndim, dimids=dimids_in )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! get equivalent dimension id's in output file:
    dimids_out = -999
    do idim = 1, ndim
      ! get name:
      status = NF90_Inquire_Dimension( self%ncid, dimids_in(idim), name=dimname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! the first two are the skip the spectral dimensions,
      ! replace this by the single rgrid dimension;
      ! skip the  first one (real/imag dim):
      if ( trim(dimname) == dimname_c2 ) cycle
      ! for the second one, fill the rgrid dim as the first:
      if ( trim(dimname) == dimname_sp ) dimname = rgridname
      ! get id in output file, fill one place before:
      status = NF90_Inq_DimID( outfile%ncid, dimname, dimids_out(idim-1) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    ! decrease:
    ndim = ndim - 1
    
    ! define new variable:
    status = NF90_Def_Var( outfile%ncid, varname, xtype, dimids_out(1:ndim), varid )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! count number of attributes:
    status = NF90_Inquire_Variable( self%ncid, varid_in, nAtts=natt )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over attributes:
    do iatt = 1, natt
      ! name:
      status = NF90_Inq_AttName( self%ncid, varid_in, iatt, attname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! switch:
      select case ( trim(attname) )
        case ( 'grid_type' )
          status = NF90_Put_Att( outfile%ncid, varid, trim(attname), 'gaussian reduced' )
          IF_NF90_NOT_OK_RETURN(status=1)
        case ( 'truncation', 'axis' )
          ! skip ...
        case default
          status = NF90_Copy_Att( self%ncid, varid_in, trim(attname), &
                                   outfile%ncid, varid )
          IF_NF90_NOT_OK_RETURN(status=1) 
      end select
    end do ! attributes
    
    ! ok
    status = 0
    
  end subroutine Spg_Define_Variable_RegriddedCopy_Rgg


  ! ***
  
  
  subroutine Spg_Copy_Variable_Rgg( self, varname, outfile, status )

    use NetCDF, only : NF90_INT, NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
    
    use C3PO_File_Nc , only : maxdim
    use C3PO_File_Rgg, only : T_File_Rgg
    use C3PO_File_Spg, only : T_File_Spg
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)     ::  self
    character(len=*), intent(in)      ::  varname
    XTYPE(T_File_Rgg), intent(inout)  ::  outfile
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Copy_Variable_Rgg'
    
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
    
  end subroutine Spg_Copy_Variable_Rgg
  
  
  ! ***


  subroutine Spg_Regrid_Variable_Rgg( spgfile, spg, varname, &
                                     rggfile, rgg, status )

    use Binas, only : deg2rad

    use NetCDF, only : NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
    
    use C3PO_Grid, only : T_Grid_Spg
    use C3PO_Grid, only : T_Grid_Rgg
    
    use C3PO_File_Nc , only : maxdim
    use C3PO_File_Rgg, only : T_File_Rgg
    use C3PO_File_Spg, only : T_File_Spg, dimname_c2, dimname_sp
  
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spg_Regrid_Variable_Rgg'
    
    ! --- in/out ---------------------------------
    
    XTYPE(T_File_Spg), intent(in)       ::  spgfile
    XTYPE(T_Grid_Spg), intent(in)       ::  spg
    character(len=*), intent(in)        ::  varname
    XTYPE(T_File_Rgg), intent(inout)    ::  rggfile
    XTYPE(T_Grid_Rgg), intent(in)       ::  rgg
    integer, intent(out)                ::  status

    ! --- local ----------------------------------
    
    integer             ::  varid, varid_out
    integer             ::  xtype
    integer             ::  ndim, idim
    integer             ::  dimids(maxdim)
    integer             ::  shp(maxdim)
    character(len=64)   ::  dimnames(maxdim)
    integer             ::  j
    integer             ::  k
    real                ::  lat_rad
    real                ::  start_lon_rad
    integer             ::  i0, nlon, nlon_full

    ! real data in default working precission:
    real, allocatable   ::  data_r_3d(:,:,:), data_r_3d_out(:,:)
    
    ! --- begin ----------------------------------
    
    ! id of input variable:
    status = NF90_INQ_VarID( spgfile%ncid, varname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! obtain info:
    status = NF90_Inquire_Variable( spgfile%ncid, varid, xtype=xtype, &
                                       ndims=ndim, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! obtain shape:
    do idim = 1, ndim
      ! get length:
      status = NF90_Inquire_Dimension( spgfile%ncid, dimids(idim), &
                                         len=shp(idim), name=dimnames(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! id of output variable:
    status = NF90_INQ_VarID( rggfile%ncid, varname, varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! switch on rank:
    select case ( ndim )

      !~ 3D
      case ( 3 )
        ! create storage and read:
        select case ( xtype )
        
          !~ r4, r8
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
            allocate( data_r_3d(shp(1),shp(2),shp(3)) )
            allocate( data_r_3d_out(rgg%npoint,shp(3)) )
            ! read:
            status = NF90_Get_Var( spgfile%ncid, varid, data_r_3d )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop over latitude bands:
            do j = 1, rgg%nlat
              ! latitude in radians:
              lat_rad = rgg%band_lats(j) * deg2rad
              ! zero based index of first point of band in rgrid:
              i0 = rgg%band_i0(j)
              ! number of lon values:
              nlon = rgg%band_nlon(j)
              ! number of lon values on full circle,
              ! compute using distance between first and last lon:
              nlon_full = nint( 360.0/(rgg%longitude(i0+nlon)-rgg%longitude(i0+1)) * (nlon-1) )
              ! first lon on band:
              start_lon_rad = rgg%longitude(i0+1) * deg2rad
              ! loop over other dims ;
              ! in future, more efficient evaluation over multiple dims:
              do k = 1, shp(3)
                ! evaluate longitude points:
                call spg%Eval_Lons_Ulat( data_r_3d(:,:,k), lat_rad, nlon_full, start_lon_rad, &
                                          nlon, data_r_3d_out(i0+1:i0+nlon,k), status )
                IF_NOT_OK_RETURN(status=1)
              end do ! dims
            end do  ! lat bands
            ! write:
            status = NF90_Put_Var( rggfile%ncid, varid_out, data_r_3d_out )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! clear:
            deallocate( data_r_3d )
            deallocate( data_r_3d_out )

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
    
  end subroutine Spg_Regrid_Variable_Rgg


end module C3PO_File_Convert

