!###############################################################################
!
! Landuse_Traffic  -   traffic intensity map
!
! Usage:
!
!    ! read traffic maps:
!    call Landuse_Traffic_Init( rcF, status )
!    if (status/=0) stop
!
!    ! clear:
!    call Landuse_Traffic_Done( status )
!    if (status/=0) stop
!
! Rcfile:
!
!  ! list with traffic intensity files:
!  landuse.traffic.files     :  highway.nc  nonhighway.nc
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
#define IF_NF90_NOTOK_STOP if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; stop; end if
!
#include "le.inc"
!
!###############################################################################

module LE_LandUse_Traffic

  use GO, only : gol, goPr, goErr

  use NetCDF, only : NF90_NOERR, nf90_strerror

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private

  public  ::  T_Traffic_Map
  public  ::  n_traffic_map, traffic_map
  
  public  ::  ROADTYPE_HIGHWAYPLUS
  public  ::  ROADTYPE_NONHIGHWAY
  public  ::  ROADTYPE_MAX
  
  public  ::  VEHICLETYPE_TRUCK
  public  ::  VEHICLETYPE_CAR
  public  ::  VEHICLETYPE_MOTORCYCLE
  public  ::  VEHICLETYPE_MAX
  
  public  ::  Landuse_Traffic_Init, Landuse_Traffic_Done


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Landuse_Traffic'
  
  ! road types:
  integer, parameter  ::  ROADTYPE_HIGHWAYPLUS = 1
  integer, parameter  ::  ROADTYPE_NONHIGHWAY  = 2
  integer, parameter  ::  ROADTYPE_MAX = ROADTYPE_NONHIGHWAY
  
  ! vehicle types:
  integer, parameter  ::  VEHICLETYPE_TRUCK      = 1
  integer, parameter  ::  VEHICLETYPE_CAR        = 2
  integer, parameter  ::  VEHICLETYPE_MOTORCYCLE = 3
  integer, parameter  ::  VEHICLETYPE_MAX = VEHICLETYPE_MOTORCYCLE

  ! maximum number of maps:
  integer, parameter    ::  max_traffic_map = 6
  

  ! --- types -----------------------------
  
  type T_Traffic_Map
    ! road type:
    integer               ::  roadtype
    ! vehicle type:
    integer               ::  vehicletype
    ! intensity:
    real, pointer         ::  km(:,:)   ! (nx,ny)
    character(len=32)     ::  units
    character(len=512)    ::  long_name
  end type T_Traffic_Map


  ! --- var -----------------------------
  
  ! maps:
  integer               ::  n_traffic_map
  type(T_Traffic_Map)   ::  traffic_map(max_traffic_map)



contains


  ! ====================================================================
  
  
  subroutine Landuse_Traffic_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine
    
    !use file_nc, only : nc_dump
    
    ! --- in/out ------------------------
    
    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_LandUse_Traffic_Init'

    ! --- local --------------------------

    character(len=1025)   ::  line
    character(len=512)    ::  fname

    !integer               ::  k
    !real, allocatable     ::  kml(:,:,:)

    ! --- begin -------------------------
    
    ! input files:
    call ReadRc( rcF, 'landuse.traffic.files', line, status )
    IF_NOTOK_RETURN(status=1)
    
    ! init map counter:
    n_traffic_map = 0
    ! loop over files:
    do
      ! nothing left ?
      if ( len_trim(line) == 0 ) exit
      ! extract:
      call goReadFromLine( line, fname, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! read from file, add result to maps array:
      call Landuse_Traffic_Add( fname, status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    !! dump ..
    !print *, 'dumping kml ...'
    !allocate( kml(nx,ny,n_traffic_map) )
    !do k = 1, n_traffic_map
    !  kml(:,:,k) = traffic_map(k)%km
    !  write (*,'("  ",a)') trim(traffic_map(k)%long_name)
    !end do
    !call nc_dump( 'km2.nc', kml, 'kml', (/'x','y','map'/), status )
    !IF_NOTOK_RETURN(status=1)
    !deallocate( kml )
    !stop 'break after dump'
    
    ! ok
    status = 0
    
  end subroutine Landuse_Traffic_Init
  
  
  ! ***
  
  
  subroutine Landuse_Traffic_Done( status )

    ! --- in/out ------------------------
    
    integer, intent(out)            ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_LandUse_Traffic_Done'

    ! --- local --------------------------

    integer               ::  imap

    ! --- begin -------------------------
    
    ! loop over maps in use:
    do imap = 1, n_traffic_map
      ! clear:
      deallocate( traffic_map(imap)%km )
    end do
    ! nothing loaded:
    n_traffic_map = 0
    
    ! ok
    status = 0
    
  end subroutine Landuse_Traffic_Done
  
  
  ! ***
  
  
  subroutine Landuse_Traffic_Add( fname, status )

    use NetCDF, only : NF90_NoWrite
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att

    use GO     , only : goCharToString

    use LE_Grid, only : ugg
    use C3PO   , only : T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors

    
    ! --- in/out ------------------------
    
    character(len=*), intent(in)    ::  fname
    integer, intent(out)            ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_LandUse_Traffic_Add'

    ! --- local --------------------------

    logical                         ::  exist
    integer                         ::  ncid
    integer                         ::  dimid
    integer                         ::  varid
    
    integer                         ::  nv
    integer                         ::  nlon
    real, allocatable               ::  lons(:)         ! (nlon)
    real, allocatable               ::  lons_bnds(:,:)  ! (nv, nlon)
    integer                         ::  ilon
    integer                         ::  nlat
    real, allocatable               ::  lats(:)         ! (nlat)
    real, allocatable               ::  lats_bnds(:,:)  ! (nv,nlat)
    integer                         ::  ilat

    !integer                         ::  ncnt
    !character(len=3), allocatable   ::  cnts(:)
    !integer                         ::  icnt

    integer                         ::  npoint
    integer, allocatable            ::  points(:)
    integer                         ::  point
    integer                         ::  ipoint

    integer                         ::  ncls, len_cls
    character(len=32), allocatable  ::  cls(:)
    character(len=1), allocatable   ::  cls_in(:,:)
    integer                         ::  icls

    real, allocatable               ::  kmg(:)   ! (npoint)
    character(len=32)               ::  units_in
    character(len=512)              ::  long_name
    integer                         ::  roadtype
    integer                         ::  vehicletype
    
    real, allocatable               ::  km_in(:,:) ! (nlon,nlat)
    type(T_Grid_Ugg)                ::  grid_in
    
    ! --- begin -------------------------

    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! open netcdf file:
    status = NF90_Open( trim(fname), NF90_NoWrite, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! nv-boundary dimension
    status = NF90_Inq_DimID( ncid, 'nv', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nv )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! lon dimension:
    status = NF90_Inq_DimID( ncid, 'lon', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlon )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lon axis:
    allocate( lons(nlon) )
    allocate( lons_bnds( nv, nlon) )
    status = NF90_Inq_VarID( ncid, 'lon', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, lons )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lon_bounds', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, lons_bnds )
    IF_NF90_NOTOK_RETURN(status=1)

    ! lat dimension:
    status = NF90_Inq_DimID( ncid, 'lat', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lat axis:
    allocate( lats(nlat) )
    allocate( lats_bnds(nv,nlat) )
    status = NF90_Inq_VarID( ncid, 'lat', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, lats )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lat_bounds', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, lats_bnds )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! Make grid-type definition
    call grid_in%Init( lons, lons_bnds, lats, lats_bnds, status )
    IF_NOTOK_RETURN(status=1)

    ! gathtered points:
    status = NF90_Inq_DimID( ncid, 'point', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=npoint )
    IF_NF90_NOTOK_RETURN(status=1)
    ! 1d indices:
    allocate( points(npoint) )
    status = NF90_Inq_VarID( ncid, 'point', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, points )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! number of road/vehicle classes:
    status = NF90_Inq_DimID( ncid, 'class', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=ncls )
    IF_NF90_NOTOK_RETURN(status=1)
    ! label length:
    status = NF90_Inq_DimID( ncid, 'len_class', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=len_cls )
    IF_NF90_NOTOK_RETURN(status=1)
    ! read as zero-terminated-string:
    allocate( cls_in(len_cls,ncls) )
    status = NF90_Inq_VarID( ncid, 'class', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, cls_in, &
                             start=(/1,1/), count=(/len_cls,ncls/) )
    IF_NF90_NOTOK_RETURN(status=1)
    ! class labels indices:
    allocate( cls(ncls) )
    call goCharToString( cls_in, cls, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( cls_in )
    
    ! storage:
    allocate( kmg(npoint) )
    allocate( km_in(grid_in%nlon,grid_in%nlat) )
    
    ! loop over classes:
    do icls = 1, ncls
      
      ! read traffic data:
      status = NF90_Inq_VarID( ncid, trim(cls(icls)), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, kmg )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Att( ncid, varid, 'long_name', long_name )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Att( ncid, varid, 'units', units_in )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! fill road and vehicle type:
      select case ( trim(long_name) )
        ! * v1
        case ( 'total distance driven on Highway roads by heavy_duty vehicles' )
           roadtype    = ROADTYPE_HIGHWAYPLUS
           vehicletype = VEHICLETYPE_TRUCK
        case ( 'total distance driven on Rural roads by heavy_duty vehicles', &
               'total distance driven on Urban roads by heavy_duty vehicles' )
           roadtype    = ROADTYPE_NONHIGHWAY
           vehicletype = VEHICLETYPE_TRUCK
        case ( 'total distance driven on Highway roads by light_duty vehicles' )
           roadtype    = ROADTYPE_HIGHWAYPLUS
           vehicletype = VEHICLETYPE_CAR
        case ( 'total distance driven on Rural roads by light_duty vehicles', &
               'total distance driven on Urban roads by light_duty vehicles' )
           roadtype    = ROADTYPE_NONHIGHWAY
           vehicletype = VEHICLETYPE_CAR
        ! * v2
        case ( 'total distance driven on highways (plus country specific part of inner roads) roads by truck vehicles' )
           roadtype    = ROADTYPE_HIGHWAYPLUS
           vehicletype = VEHICLETYPE_TRUCK
        case ( 'total distance driven on highways (plus country specific part of inner roads) roads by car vehicles' )
           roadtype    = ROADTYPE_HIGHWAYPLUS
           vehicletype = VEHICLETYPE_CAR
        case ( 'total distance driven on highways (plus country specific part of inner roads) roads by motorcycle vehicles' )
           roadtype    = ROADTYPE_HIGHWAYPLUS
           vehicletype = VEHICLETYPE_MOTORCYCLE
        case ( 'total distance driven on non-highways roads by truck vehicles' )
           roadtype    = ROADTYPE_NONHIGHWAY
           vehicletype = VEHICLETYPE_TRUCK
        case ( 'total distance driven on non-highways roads by car vehicles' )
           roadtype    = ROADTYPE_NONHIGHWAY
           vehicletype = VEHICLETYPE_CAR
        case ( 'total distance driven on non-highways roads by motorcycle vehicles' )
           roadtype    = ROADTYPE_NONHIGHWAY
           vehicletype = VEHICLETYPE_MOTORCYCLE
        ! *
        case default
          write (gol,'("unsuported long_name `",a,"`")') trim(long_name); call goErr
          TRACEBACK; status=1; return
      end select
      
      ! init 2d field to zero:
      km_in = 0.0
      ! loop over points:
      do ipoint = 1, npoint
        ! 0-based 1D index:
        point = points(ipoint)
        ! extract 1-based indices:
        !icnt =        point                          /(nlat*nlon)    + 1
        ilat = modulo(point,grid_in%nlat*grid_in%nlon)/ grid_in%nlon  + 1
        ilon = modulo(point,             grid_in%nlon)                + 1
        !! info ..
        !if ( ipoint <= 10 ) print *, 'ppp ', ipoint, ';', icnt, ilat, ilon
        ! add contribution from this country:
        km_in(ilon,ilat) = km_in(ilon,ilat) + kmg(ipoint)
      end do
      
      ! convert if necessary:
      select case ( trim(units_in) )
        case ( 'Gm' )
          ! mlj km, expected this ...
        case default
          write (gol,'("unsupported traffic units `",a,"`")') trim(units_in); call goErr
          TRACEBACK; status=1; return
      end select
      
      ! convert input to aera average:
      call grid_in%AreaOper( km_in, '/', 'm2', status )
      IF_NOTOK_RETURN(status=1)
      
      ! increase counter:
      n_traffic_map = n_traffic_map + 1
      ! check ...
      if ( n_traffic_map > max_traffic_map ) then
        write (gol,'("number of maps exceeds maxium of ",i6)') size(traffic_map); call goErr
        write (gol,'("increase parameter `max_traffic_map` in module `",a,"`")') mname; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! fill:
      traffic_map(n_traffic_map)%roadtype    = roadtype
      traffic_map(n_traffic_map)%vehicletype = vehicletype
      traffic_map(n_traffic_map)%long_name   = trim(long_name)
      traffic_map(n_traffic_map)%units = trim(units_in)

      ! storage for target grid:
      allocate( traffic_map(n_traffic_map)%km(ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! regrid, maintain average:
      call Grid_Convertors%Ugg_AreaAver( grid_in, km_in, ugg, traffic_map(n_traffic_map)%km, status )
      IF_NOTOK_RETURN(status=1)

      ! convert to km/cell :
      call ugg%AreaOper( traffic_map(n_traffic_map)%km, '*', 'm2', status )
      IF_NOTOK_RETURN(status=1)
            
    end do

    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    !deallocate( cnts )
    deallocate( points )
    deallocate( cls )
    deallocate( kmg )
    deallocate( km_in )
    
    ! clear
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0

  end subroutine Landuse_Traffic_Add


  ! ***
  

end module LE_LandUse_Traffic
