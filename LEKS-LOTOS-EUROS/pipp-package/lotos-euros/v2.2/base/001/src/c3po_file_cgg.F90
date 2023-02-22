!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!   C3PO_File_Cgg  : access files with data defined on carthesian-gaussian-grid
!
! EXAMPLE
!
!    use GO     , only : TDate, NewDate
!    use C3PO   , only : T_File_Cgg
!    use C3PO   , only : T_Grid_Cgg
!    use C3PO   , only : T_Levs_Hyb
!
!    ! variables:
!    type(T_File_Cgg)      ::  infile
!    type(T_Grid_Cgg)      ::  cgg
!    type(T_Levs_Hyb)      ::  levs
!    type(TDate)           ::  t
!    integer               ::  irec
!    real, allocatable     ::  values(:,:)  ! (nlon,nlat,nlev)
!    integer               ::  status
!
!    ! open file:
!    call infile%Open( 'temperature.nc', status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! init grid definition from file:
!    call infile%Get_Grid( cgg, status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! get level definition
!    call infile%Get_Levs( levs, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! target time:
!    t = NewDate( 2000, 01, 02, 00, 00 )
!    
!    ! index of time record:
!    call infile%Inq_TimeRecord( 'time', t, irec, status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! storage:
!    allocate( values_in(cgg%nlon,cgg%nlat,levs%nlev) )
!
!    ! read data slab:
!    call infile%Get_Var( trim(nametype), trim(name), values_in, units_in, status, &
!                              start=(/1,1,irec/), count=(/cgg%nlon,cgg%nlat,levs%nlev,1/) )
!    IF_NOT_OK_RETURN(status=1)
!      
!    ! clear input field:
!    deallocate( values_in )
!    
!    ! clear input grid definition:
!    call cgg%Done( status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! clear input level definition:
!    call levs%Done( status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! close input file:
!    call infile%Close( status )
!    IF_NOT_OK_RETURN(status=1)
!
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

module C3PO_File_Cgg

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError
  
  use C3PO_File_Nc, only : T_File_Nc, maxdim

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_File_Cgg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Cgg'

  ! cgg parameters:
  character(len=*), parameter  ::  dimname_nv   = 'nv'
  character(len=*), parameter  ::  varext_bnds  = '_bnds'
  integer, parameter           ::  nv = 2
  

  ! --- types ----------------------------------------
  
  type, extends(T_File_Nc) :: T_File_Cgg
    ! cgg dims:
    integer                 ::  dimid_lat
    integer                 ::  dimid_lon
    integer                 ::  dimid_nv
    ! cgg variables:
    integer                 ::  varid_longitude
    integer                 ::  varid_latitude
    integer                 ::  varid_longitude_bnds
    integer                 ::  varid_latitude_bnds
  contains
    procedure   ::  Get_Grid                => File_Cgg_Get_Grid
!    procedure   ::  Def_Grid                => File_Cgg_Def_Grid
!    procedure   ::  Put_Grid                => File_Cgg_Put_Grid
!    procedure   ::  Copy_Variable_Selection => File_Cgg_Copy_Variable_Selection
  end type T_File_Cgg

  
  

contains


  ! ********************************************************************
  ! ***
  ! *** cgg file
  ! ***
  ! ********************************************************************


  
  subroutine File_Cgg_Get_Grid( self, vid, cgg, status )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Get_Var

    use Num          , only : Interp_Lin, EXTRAPOL_LINEAR
    use C3PO_Grid_CGG, only : T_Grid_CGG
  
    ! --- in/out ---------------------------------
    
    class(T_File_Cgg), intent(in)       ::  self
    integer, intent(in)                 ::  vid
    class(T_Grid_CGG), intent(out)      ::  cgg
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Cgg_Get_Grid'
    
    ! --- local ----------------------------------

    integer                 ::  ndims
    integer, allocatable    ::  dimids(:)
    integer                 ::  idim
    integer                 ::  dimid
    character(len=32)       ::  dimname
    integer                 ::  varid
    integer                 ::  nlon, nlat, nv
    real, allocatable       ::  lats(:)
    real, allocatable       ::  lats_bnds(:,:)
    real, allocatable       ::  lons(:)
    real, allocatable       ::  lons_bnds(:,:)
    real, allocatable       ::  ii(:)
    integer                 ::  i
    
    ! --- begin ----------------------------------
    
    ! number of boundary values:
    nv = 2
    
    ! get  number of dimensions:
    status = NF90_Inquire_Variable( self%ncid, vid, ndims=ndims )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( dimids(ndims), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! get dimension id's:
    status = NF90_Inquire_Variable( self%ncid, vid, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)

    !~ longitudes
    
    ! select first dimension:
    dimid = dimids(1)
    ! get dimension name and length:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=dimname, len=nlon )
    IF_NF90_NOT_OK_RETURN(status=1)
        
    ! storage:
    allocate( lons(nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! variable with longitudes:
    status = NF90_INQ_VarID( self%ncid, dimname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, lons )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! some files have lons (...,358,359,0,1,2,...)
    if ( lons(1) > lons(nlon) ) then
      ! shift:
      do i = 1, nlon
        if ( lons(i) > lons(nlon) ) then
          lons(i) = lons(i) - 360.0
        end if
      end do
    end if

    ! range:
    allocate( ii(nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do i = 1, nlon
      ii(i) = i
    end do
    ! storage:
    allocate( lons_bnds(nv,nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! interpolate boundaries:
    call Interp_Lin( ii, lons, ii-0.5, lons_bnds(1,:), status, extrapol=EXTRAPOL_LINEAR )
    IF_NOT_OK_RETURN(status=1)
    call Interp_Lin( ii, lons, ii+0.5, lons_bnds(2,:), status, extrapol=EXTRAPOL_LINEAR )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ latitudes
    
    ! select second dimension:
    dimid = dimids(2)
    ! get dimension name and length:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=dimname, len=nlat )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( lats(nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! variable with latgitudes:
    status = NF90_INQ_VarID( self%ncid, dimname, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, lats )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! range:
    allocate( ii(nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do i = 1, nlat
      ii(i) = i
    end do
    ! storage:
    allocate( lats_bnds(nv,nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! interpolate boundaries:
    call Interp_Lin( ii, lats, ii-0.5, lats_bnds(1,:), status, extrapol=EXTRAPOL_LINEAR )
    IF_NOT_OK_RETURN(status=1)
    call Interp_Lin( ii, lats, ii+0.5, lats_bnds(2,:), status, extrapol=EXTRAPOL_LINEAR )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~
    
    ! define grid:
    call cgg%Init( lons, lons_bnds, lats, lats_bnds, status )
    IF_NOT_OK_RETURN(status=1)   
    
    ! clear:
    deallocate( lons, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lons_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( dimids, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine File_Cgg_Get_Grid


  ! ***
  
  
!  subroutine File_Cgg_Def_Grid( self, cgg, status )
!
!    use NetCDF, only : NF90_SHORT, NF90_INT, NF90_FLOAT
!    use NetCDF, only : NF90_Def_Dim
!    use NetCDF, only : NF90_Def_Var
!    use NetCDF, only : NF90_Put_Att
!
!    use C3PO_Grid_RGG, only : T_Grid_CGG
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Cgg), intent(inout)    ::  self
!    class(T_Grid_CGG), intent(in)       ::  cgg
!    integer, intent(out)                ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Cgg_Def_Grid'
!    
!    ! --- local ----------------------------------
!    
!    integer               ::  varid
!    character(len=1024)   ::  line
!
!    ! --- begin ----------------------------------
!
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_lat, cgg%nlat, self%dimid_ulat )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_lon, cgg%nulon, self%dimid_ulon )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_rgrid, cgg%npoint, self%dimid_rgrid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_nv, cgg%nv, self%dimid_nv )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! ~ index
!    
!    ! add index variable:
!    status = NF90_Def_Var( self%ncid, dimname_rgrid, NF90_INT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'reduced grid index' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    write (line,'(a," ",a)') dimname_lat, dimname_lon
!    status = NF90_Put_Att( self%ncid, varid, 'compress', trim(line) )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    write (line,'("zero based index that specifies location in regular (",a,",",a,") grid : ", &
!                & " j = rgrid(n) / ",a," ; i = rgrid(n) - ",a," * j")') &
!             dimname_lat, dimname_lon, dimname_lon, dimname_lon
!    IF_NF90_NOT_OK_RETURN(status=1)
!    status = NF90_Put_Att( self%ncid, varid, 'description', trim(line) )
!    !
!    self%varid_irgrid = varid
!        
!    ! ~ centers
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lon, NF90_FLOAT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'longitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'longitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_east' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', varname_lon//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_longitude = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lat, NF90_FLOAT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', varname_lat//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_latitude = varid
!
!    ! ~ boundary values
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lon//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'longitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_east' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_longitude_bnds = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lat//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_latitude_bnds = varid
!
!    ! ~ latitude bands (extra info)
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, dimname_lat, NF90_FLOAT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', dimname_lat//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_lats = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, dimname_lat//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_lats_bnds = varid
!    
!    ! number of longitudes per row:
!    status = NF90_Def_Var( self%ncid, dimname_lat//'_nlon', NF90_SHORT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'number of longitudes on latitude band' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_nlon = varid
!
!    ! index in reduced grid:
!    status = NF90_Def_Var( self%ncid, dimname_lat//'_i0', NF90_INT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'zero-based index in reduced grid of first point on latitude band' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_i0 = varid
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Cgg_Def_Grid
!
!
!  ! ***
!  
!  
!  subroutine File_Cgg_Put_Grid( self, cgg, status )
!
!    use NetCDF, only : NF90_Put_Var
!
!    use C3PO_Grid_RGG, only : T_Grid_CGG
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Cgg), intent(inout)    ::  self
!    class(T_Grid_CGG), intent(in)       ::  cgg
!    integer, intent(out)                ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Cgg_Put_Grid'
!    
!    ! --- local ----------------------------------
!    
!    ! --- begin ----------------------------------
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_irgrid, cgg%irgrid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_longitude, cgg%longitude )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_latitude, cgg%latitude )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_longitude_bnds, cgg%longitude_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_latitude_bnds, cgg%latitude_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_lats, cgg%band_lats )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_lats_bnds, cgg%band_lats_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_nlon, cgg%band_nlon )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_i0, cgg%band_i0 )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Cgg_Put_Grid
!
!
!  ! ***
!  
!  
!  subroutine File_Cgg_Copy_Variable_Selection( self, varname, outfile, selection, status )
!
!    use NetCDF, only : NF90_INT, NF90_FLOAT, NF90_DOUBLE
!    use NetCDF, only : NF90_Inquire_Dimension
!    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
!    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Cgg), intent(in)     ::  self
!    character(len=*), intent(in)      ::  varname
!    class(T_File_Cgg), intent(inout)  ::  outfile
!    integer, intent(in)               ::  selection(:)  ! 1-based indices in rgrid
!    integer, intent(out)              ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Cgg_Copy_Variable_Selection'
!    
!    ! --- local ----------------------------------
!    
!    integer             ::  varid, varid_out
!    integer             ::  xtype
!    integer             ::  ndim, idim
!    integer             ::  dimids(maxdim)
!    integer             ::  shp(maxdim)
!    logical             ::  rgridim(maxdim)
!    character(len=64)   ::  dimname
!    integer             ::  nselected
!    integer             ::  i
!    
!    integer(4), allocatable   ::  data_i4_1d(:), data_i4_1d_out(:)
!    real(4), allocatable      ::  data_r4_1d(:), data_r4_1d_out(:)
!    real(8), allocatable      ::  data_r8_1d(:), data_r8_1d_out(:)
!    
!    integer(4), allocatable   ::  data_i4_2d(:,:), data_i4_2d_out(:,:)
!    real(4), allocatable      ::  data_r4_2d(:,:), data_r4_2d_out(:,:)
!    real(8), allocatable      ::  data_r8_2d(:,:), data_r8_2d_out(:,:)
!    
!    integer(4), allocatable   ::  data_i4_3d(:,:,:), data_i4_3d_out(:,:,:)
!    real(4), allocatable      ::  data_r4_3d(:,:,:), data_r4_3d_out(:,:,:)
!    real(8), allocatable      ::  data_r8_3d(:,:,:), data_r8_3d_out(:,:,:)
!    
!    ! --- begin ----------------------------------
!    
!    ! id of input variable:
!    status = NF90_INQ_VarID( self%ncid, varname, varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! obtain info:
!    status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype, &
!                                       ndims=ndim, dimids=dimids )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! obtain shape:
!    do idim = 1, ndim
!      ! get length:
!      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), len=shp(idim) )
!      IF_NF90_NOT_OK_RETURN(status=1)
!    end do
!    
!    ! init flags for rgrid dim:
!    rgridim = .false.
!    do idim = 1, ndim
!      ! get name:
!      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dimname )
!      IF_NF90_NOT_OK_RETURN(status=1)
!      ! compare:
!      rgridim(idim) = trim(dimname) == trim(dimname_rgrid)
!    end do
!    ! check ...
!    if ( count(rgridim) > 1 ) then
!      write (gol,'("could not select points for ",i2," rgrid dimensions")') count(rgridim); call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! id of output variable:
!    status = NF90_INQ_VarID( outfile%ncid, varname, varid_out )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! number of selected grid points:
!    nselected = size( selection )
!    
!    ! switch on rank:
!    select case ( ndim )
!
!      !~ 1D
!      case ( 1 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_i4_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_i4_1d_out(i) = data_i4_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_i4_1d_out(shp(1)) )
!              ! copy:
!              data_i4_1d_out = data_i4_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_1d )
!            deallocate( data_i4_1d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_r4_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_r4_1d_out(i) = data_r4_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_r4_1d_out(shp(1)) )
!              ! copy:
!              data_r4_1d_out = data_r4_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_1d )
!            deallocate( data_r4_1d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_r8_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_r8_1d_out(i) = data_r8_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_r8_1d_out(shp(1)) )
!              ! copy:
!              data_r8_1d_out = data_r8_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_1d )
!            deallocate( data_r8_1d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!
!      !~ 2D
!      case ( 2 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_i4_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_i4_2d_out(i,:) = data_i4_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_i4_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_i4_2d_out = data_i4_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_2d )
!            deallocate( data_i4_2d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r4_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r4_2d_out(i,:) = data_r4_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r4_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_r4_2d_out = data_r4_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_2d )
!            deallocate( data_r4_2d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r8_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r8_2d_out(i,:) = data_r8_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r8_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_r8_2d_out = data_r8_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_2d )
!            deallocate( data_r8_2d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!
!      !~ 3D
!      case ( 3 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_i4_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_i4_3d_out(i,:,:) = data_i4_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_i4_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_i4_3d_out = data_i4_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_3d )
!            deallocate( data_i4_3d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r4_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r4_3d_out(i,:,:) = data_r4_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r4_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_r4_3d_out = data_r4_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_3d )
!            deallocate( data_r4_3d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r8_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r8_3d_out(i,:,:) = data_r8_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r8_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_r8_3d_out = data_r8_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_3d )
!            deallocate( data_r8_3d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!        
!      !~
!      case default
!        write (gol,'("unsupported ndim ",i6)') ndim; call goErr
!        TRACEBACK; status=1; return
!    end select
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Cgg_Copy_Variable_Selection


end module C3PO_File_Cgg
